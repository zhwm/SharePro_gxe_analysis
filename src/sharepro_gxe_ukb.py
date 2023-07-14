import pandas as pd
import argparse
import os
import time
import numpy as np
from scipy.special import softmax, expit
from scipy.stats import chi2, entropy
import scipy.sparse as sparse

np.set_printoptions(precision=4, linewidth=200)


def title():
    print('**********************************************************************')
    print('* SharePro for genome-wide GxE analysis                              *')
    print('* Version 1.0.0                                                      *')
    print('* (C) Wenmin Zhang (wenmin.zhang@mail.mcgill.ca)                     *')
    print('**********************************************************************')
    print()


# obtained from https://storage.googleapis.com/broad-alkesgroup-public/UKBB_LD/readme_ld.txt
def load_ld_npz(ld_prefix):
    # load the SNPs metadata
    gz_file = '%s.gz' % (ld_prefix)
    df_ld_snps = pd.read_table(gz_file, sep='\s+')
    df_ld_snps.rename(columns={'rsid': 'SNP', 'chromosome': 'CHR', 'position': 'BP', 'allele1': 'A1', 'allele2': 'A2'},
                      inplace=True, errors='ignore')
    assert 'SNP' in df_ld_snps.columns
    assert 'CHR' in df_ld_snps.columns
    assert 'BP' in df_ld_snps.columns
    assert 'A1' in df_ld_snps.columns
    assert 'A2' in df_ld_snps.columns
    df_ld_snps.index = df_ld_snps['CHR'].astype(str) + '.' + df_ld_snps['BP'].astype(str) + '.' + df_ld_snps[
        'A1'] + '.' + df_ld_snps['A2']
    # load the LD matrix
    npz_file = '%s.npz' % (ld_prefix)
    try:
        R = sparse.load_npz(npz_file).toarray()
        R += R.T
    except ValueError:
        raise IOError('Corrupt file: %s' % (npz_file))
    df_R = pd.DataFrame(R, index=df_ld_snps.index, columns=df_ld_snps.index)
    return df_R, df_ld_snps


def get_HESS_h2_z(LD, Z, N, ptLD=0.2, ptp=1e-5):
    """calculate local heritabilities"""
    zsquare = Z ** 2
    pvalmin = chi2.sf(max(zsquare), 1)
    idx_retain = []
    idx_exclude = [i for i in range(len(Z))]
    while len(idx_exclude) > 0:  # P + T
        maxid = idx_exclude[np.argmax(zsquare[idx_exclude])]  # find the idx with smallest p-value
        pval = chi2.sf(zsquare[maxid], 1)
        idx_retain.append(maxid)
        idx_exclude = [i for i in idx_exclude if i not in np.where(np.abs(LD[maxid, :]) > ptLD)[0]]
        if pval > ptp:
            break
    Indidx = np.sort(idx_retain)  # obtain independent signals
    P = len(Indidx)
    LD_id = LD[np.ix_(Indidx, Indidx)]
    R_inv = np.linalg.inv(LD_id)
    vec_id = Z[Indidx]
    h2_hess = (np.dot(np.dot(vec_id.transpose(), R_inv), vec_id) - P) / (N - P)
    var_b = np.max(zsquare[Indidx]) / N
    if h2_hess < 0.0001:
        h2_hess = 0.0001
    if h2_hess > 0.9:
        h2_hess = 0.9
    return h2_hess, var_b, pvalmin


def get_ld_ukb(ldfile, z, start, end, args):
    """Get matched GWAS summary statistics for each block"""
    if args.LDdir.count(args.LDdir[0]) == len(args.LDdir):
        df_R, df_ld_snps = zip(*[load_ld_npz(os.path.join(args.LDdir[0], ldfile))] * len(args.LDdir))  # same ld
    else:
        df_R, df_ld_snps = zip(*[load_ld_npz(os.path.join(args.LDdir[i], ldfile)) for i in range(len(args.zdir))])
    idx = sorted(set.intersection(*[set(i.index) for i in df_R]).intersection(set(z.index)))
    if len(idx) < 10:
        print("Less than 10 variants matched in the range of {} to {}, skipping".format(start, end))
        return [], [], [], [], [], [], [], []
    pos = [int(i.split('.')[1]) for i in idx]
    effidx = [i for i in range(len(idx)) if ((pos[i] >= start) & (pos[i] < end))]
    effnum = len(effidx)
    if effnum <= 10:
        print('Less than 10 effective variants in the range of {} to {}, skipping'.format(start, end))
        return [], [], [], [], [], [], [], []
    print('{} variants in the range of {} to {}'.format(effnum, start, end))
    Z = z.loc[idx].values
    XX = np.ones(Z.shape) * args.N
    ytX = Z * np.sqrt(args.N)
    XtX = [i.loc[idx, idx].values * j for i, j in zip(df_R, args.N)]
    hess, varb, pchi = zip(*[get_HESS_h2_z(df_R[i].loc[idx, idx].values, Z[:, i], args.N[i])
                             for i in range(len(args.zdir))])
    return idx, effidx, XX, ytX, XtX, hess, varb, pchi


class SparseReg(object):
    def __init__(self, P, K, XX, h2, var_b, sigma):
        """initialize and set hyperparameters"""
        self.p = P
        self.k = K
        self.beta_mu = np.zeros((self.p, self.k))
        self.beta_prior_tau = np.tile((1.0 / var_b * np.array([k + 1 for k in range(self.k)])), (self.p, 1))
        self.y_tau = 1.0 / (1 - h2)
        self.sigma = sigma
        self.beta_post_tau = np.tile(XX.reshape(-1, 1), (1, self.k)) * self.y_tau + self.beta_prior_tau
        self.delta = np.zeros((self.p, self.k))
        self.u1 = 0.5 * np.log(self.beta_post_tau / self.beta_prior_tau) + np.log(sigma / (1 - sigma))

    def infer_q_beta(self, ytX, XtX, GAM, k):
        """perform variational updates for the k-th effect"""
        idxall = [x for x in range(self.k)]
        idxall.remove(k)
        beta_all_k = (GAM[:, idxall] * self.beta_mu[:, idxall] * self.delta[:, idxall]).sum(axis=1)
        self.beta_mu[:, k] = (ytX - np.dot(beta_all_k, XtX)) / self.beta_post_tau[:, k] * self.y_tau
        u = self.u1[:, k] + 0.5 * self.beta_mu[:, k] ** 2 * self.beta_post_tau[:, k]
        self.delta[:, k] = expit(u)
        return u

    def get_elbo(self, XX, ytX, XtX, GAM):
        """get elbo"""
        beta_all = (GAM * self.delta * self.beta_mu).sum(axis=1)
        ll1 = self.y_tau * np.dot(beta_all, ytX)
        ll2 = - 0.5 * self.y_tau * (((GAM * self.delta * (self.beta_mu ** 2 +
                                                          1 / self.beta_post_tau)).sum(axis=1) * XX).sum())
        W = GAM * self.delta * self.beta_mu
        WtRW = np.dot(np.dot(W.transpose(), XtX), W)
        ll3 = - 0.5 * self.y_tau * (WtRW.sum() - np.diag(WtRW).sum())
        ll = ll1 + ll2 + ll3
        mklbeta = - 0.5 * (GAM * self.delta * (self.beta_prior_tau * (self.beta_mu ** 2 + 1 / self.beta_post_tau) +
                                               np.log(self.beta_post_tau / self.beta_prior_tau) - 1)).sum()
        nozerodelta = self.delta[self.delta != 0]
        nozeroGAM = GAM[self.delta != 0]
        noonedelta = self.delta[self.delta != 1]
        nooneGAM = GAM[self.delta != 1]
        mkldelta = (nozeroGAM * nozerodelta * np.log(self.sigma / nozerodelta)).sum() - \
                   (nooneGAM * (1 - noonedelta) * np.log((1 - self.sigma) / (1 - noonedelta))).sum()
        return ll, mklbeta, mkldelta


class SharePro(object):
    def __init__(self, P, K, XX, h2, varb, sigma):
        """initialize and set hyperparameters"""
        self.num = XX.shape[1]
        self.SR = [SparseReg(P, K, XX[:, i], h2, varb, sigma) for i in range(self.num)]
        self.p = P
        self.k = K
        self.gamma = np.zeros((self.p, self.k))
        self.sigma = sigma
        self.prior_pi = np.ones((self.p,)) * (1 / self.p)
        self.usum1 = self.num * np.log(1 - self.sigma) + np.log(self.prior_pi.transpose())

    def infer_q_s(self, ytX, XtX):
        """perform variational updates"""
        for k in range(self.k):
            idxall = [x for x in range(self.k)]
            idxall.remove(k)
            u12 = np.array([self.SR[i].infer_q_beta(ytX[:, i], XtX[i], self.gamma, k) for i in range(len(self.SR))])
            usum = np.log(1 + np.exp(-u12))
            uall = u12.sum(axis=0) + usum.sum(axis=0) + self.usum1
            self.gamma[:, k] = softmax(uall)

    def get_elbo(self, XX, ytX, XtX):
        llsum, mklbetasum, mkldeltasum = zip(*[self.SR[i].get_elbo(XX[:, i], ytX[:, i], XtX[i], self.gamma)
                                               for i in range(len(self.SR))])
        gammaterm1 = (self.gamma * np.tile(self.prior_pi.reshape(-1, 1), (1, self.k))).sum()
        gammaterm2 = (self.gamma[self.gamma != 0] * np.log(self.gamma[self.gamma != 0])).sum()
        mklgamma = gammaterm1 - gammaterm2
        ll = sum(llsum)
        mklbeta = sum(mklbetasum)
        mkldelta = sum(mkldeltasum)
        elbo = ll + mklbeta + mkldelta + mklgamma
        return ll, mklbeta, mkldelta, mklgamma, elbo

    def train(self, XX, ytX, XtX, maxite=50, eps=0.1, verbose=True, loss=0.0):
        for ite in range(maxite):
            self.infer_q_s(ytX, XtX)
            ll, mklbeta, mkldelta, mklgamma, elbo = self.get_elbo(XX, ytX, XtX)
            if verbose:
                print('*' * 70)
                print('Iteration-->{} . Likelihood: {:.1f} . KL_b: {:.1f} . KL_c: {:.1f} . KL_s: {:.1f} . ELBO: {:.1f}'
                      .format(ite, ll, mklbeta, mkldelta, mklgamma, elbo))
            if abs(elbo - loss) < eps:
                break
            if ite == (maxite - 1):
                print("Algorithm not converged. Please make sure matched summary statistics and LD were provided!")
            loss = elbo

    def multiply_specific(self, i):
        """calculate c1(1-c2)*s"""
        allidx = [x for x in range(self.num)]
        allidx.remove(i)
        matspecific = self.SR[i].delta * (np.prod(np.array([(1 - self.SR[x].delta) for x in allidx]), axis=0))
        return matspecific

    def multiply_delta(self):
        """calculate c1*c2*s"""
        return np.prod(np.array([i.delta for i in self.SR]), axis=0)

    def get_summary(self, cthres=0.95, ethres=50):
        """get variant and effect level summary"""
        matidx = np.argsort(-self.gamma, axis=0)
        variantidx = np.argmax(self.gamma, axis=1).tolist()
        vgamma = np.max(self.gamma, axis=1)
        mat_eff = np.zeros((self.p, self.k))  # effective gamma
        mat_eff[range(self.p), variantidx] = self.gamma[range(self.p), variantidx]
        # matdelta = self.multiply_delta()
        # mat_specific = [self.multiply_specific(i) for i in range(self.num)]
        csum = mat_eff.sum(axis=0).round(2)
        print("Attainable coverage for effect groups: {}".format(csum))  # statistical evidence
        eff = {}
        eff_gamma = {}
        eff_mu = {}
        # eff_share = {}
        # eff_specific = {}
        # eff_pdiff = {}
        # eff_tau = {}
        # eff_c = {}
        for k in range(self.k):
            if csum[k] > cthres:
                if entropy(mat_eff[:, k]) < np.log(ethres):
                    for p in range(self.p):
                        if np.sum(mat_eff[matidx[0:p, k], k]) > cthres * csum[k] or mat_eff[matidx[p, k], k] < 0.01:
                            eff[k] = matidx[0:p, k].tolist()
                            eff_gamma[k] = mat_eff[eff[k], k].round(4).tolist()
                            effmuk = [i.beta_mu[eff[k][0], k].round(4) for i in self.SR]
                            eff_mu[k] = effmuk
                            # effgamma_n = eff_gamma[k] / csum[k]
                            # eff_share[k] = sum(np.multiply(matdelta[eff[k], k], effgamma_n)).round(4)
                            # eff_specific[k] = [sum(np.multiply(i[eff[k], k], effgamma_n)).round(4) for i in mat_specific]
                            # eff_c[k] = [sum(np.multiply(i.delta[eff[k], k], effgamma_n)).round(4) for i in self.SR]
                            # efftauk = [i.beta_post_tau[eff[k][0], k].round(4) for i in self.SR]
                            # eff_pdiff[k] = chi2.sf((effmuk[0] - effmuk[-1]) ** 2 / (1 / efftauk[0] + 1 / efftauk[-1]), 1)
                            # eff_tau[k] = efftauk
                            break
        return variantidx, vgamma, eff, eff_gamma, eff_mu


def check_effect(eff, XtX):
    cidx = [v[0] for k, v in eff.items()]
    sLD = XtX[cidx][:, cidx]
    maxld = np.max(np.abs(np.tril(sLD, k=-1)))/np.diag(sLD)[0]
    return maxld


def ukb(args):
    print("Using genome-wide mode with --ukb")
    z = pd.concat([pd.read_csv(i, sep="\s+", header=None, index_col=0) for i in args.zdir], axis=1, join='inner')
    print("summary statistics {} loaded at {}".format(z.shape, time.strftime("%Y-%m-%d %H:%M")))
    ldlists = pd.read_csv(args.ukb, sep='\s+', dtype={'ld': str, 'start': int, 'end': int})
    print("LD list with {} LD blocks loaded\n".format(len(ldlists)))
    allv, allgamma = [], []
    cs, cs_prob, cs_z, cs_eff, cs_pdiff, cs_tau = [], [], [], [], [], []
    for ite in range(len(ldlists)):
        ld, start, end = ldlists.iloc[ite, 0:3]
        ldfile = ld.replace('.npz', '')
        idx, effidx, XX, ytX, XtX, hess, varb, pvalvec = get_ld_ukb(ldfile, z, start, end, args)
        ldlists.at[ite, 'h2'] = ','.join(['{:.2e}'.format(i) for i in hess])
        ldlists.at[ite, 'varb'] = ','.join(['{:.2e}'.format(i) for i in varb])
        if len(effidx) == 0:
            continue
        if args.hess is not None:
            h2 = args.hess
        else:
            h2 = min(hess)
        if args.varb is not None:
            b2 = args.varb
        else:
            b2 = max(varb)
        #model = SharePro(len(idx), args.K, XX, h2, b2, sigma=args.sigma)
        #model.train(XX, ytX, XtX, verbose=args.verbose)
        #variantidx, vgamma, eff, eff_gamma, eff_share, eff_specific, eff_c, eff_mu, eff_tau, eff_pdiff = \
        #    model.get_summary(cthres=args.cthres, ethres=args.ethres)
        maxld = 1.0
        ldthres = 0.5
        while maxld > ldthres:
            model = SharePro(len(idx), args.K, XX, h2, b2, sigma=args.sigma)
            model.train(XX, ytX, XtX, verbose=args.verbose)
            variantidx, vgamma, eff, eff_gamma, eff_mu = model.get_summary(cthres=args.cthres, ethres=args.ethres)
            if len(eff) > 0.0:
                maxld = np.max([check_effect(eff, i) for i in XtX])
            else:
                maxld = 0.0
            b2 = b2 * 0.1
        eff_pdiff = {key: chi2.sf((val[0] - val[-1]) ** 2 / (1 / args.N[0] + 1 / args.N[-1]), 1)
                     for key, val in eff_mu.items()}
        allv.extend([idx[i] for i in effidx])
        allgamma.extend(['{:.2e}'.format(vgamma[i]) for i in effidx])
        for e in eff:
            if eff[e][0] in effidx:
                mcs_idx = [idx[j] for j in eff[e]]
                print('The {}-th effect contains effective variants:'.format(e))
                print('causal variants: {}'.format(mcs_idx))
                print('variant probabilities for this effect: {}'.format(eff_gamma[e]))
                #print('shared probability for this effect: {}'.format(eff_share[e]))
                #print('specific probability for this effect: {}'.format(eff_specific[e]))
                #print('probability of effect for traits: {}'.format(eff_c[e]))
                print('causal effect size for traits: {}'.format(eff_mu[e]))
                print('p value for effect size difference: {}'.format(eff_pdiff[e]))
                print()
                cs.append(mcs_idx)
                cs_prob.append(eff_gamma[e])
                cs_z.append(z.loc[idx[eff[e][0]]].tolist())
                #cs_share.append(eff_share[e])
                #cs_specific.append(eff_specific[e])
                #cs_c.append(eff_c[e])
                cs_eff.append(eff_mu[e])
                cs_pdiff.append(eff_pdiff[e])
                #cs_tau.append(eff_tau[e])
    df_z = z.loc[allv].copy()
    df_z.columns = args.zdir
    df_z['vProb'] = allgamma
    df_z.index.name = 'SNP'
    df_z.to_csv(os.path.join(args.save, "{}.snp".format(args.prefix)), sep='\t', header=True, index=True)
    allcs = pd.DataFrame({"cs": ['/'.join(i) for i in cs],
                          #"zscore": [','.join([str(j) for j in i]) for i in cs_z],
                          "p_diff": cs_pdiff,
                          "beta": [','.join([str(j) for j in i]) for i in cs_eff],
                          #"share": cs_share,
                          #"tau": [','.join([str(j) for j in i]) for i in cs_tau],
                          #"specific": [','.join([str(j) for j in i]) for i in cs_specific],
                          #"causalProb": [','.join([str(j) for j in i]) for i in cs_c],
                          "variantProb": ['/'.join([str(j) for j in i]) for i in cs_prob]})
    allcs.to_csv(os.path.join(args.save, "{}.cs".format(args.prefix)), sep='\t', header=True, index=False)
    ldlists.to_csv(os.path.join(args.save, "{}.h2".format(args.prefix)), sep='\t', header=True, index=False)


parser = argparse.ArgumentParser(description='SharePro Commands:')
parser.add_argument('--ukb', type=str, default=None, help='genome-wide mode: path to LD lists', required=True)
parser.add_argument('--zdir', type=str, default=None, nargs='+', help='path to zscores files', required=True)
parser.add_argument('--LDdir', type=str, default=None, nargs='+', help='path to ld files', required=True)
parser.add_argument('--N', type=int, default=None, nargs='+', help='sample size', required=True)
parser.add_argument('--save', type=str, default=None, help='path to save result', required=True)
parser.add_argument('--prefix', type=str, default=None, help='prefix for result files', required=True)
parser.add_argument('--verbose', action="store_true", help='options for displaying more information')
parser.add_argument('--K', type=int, default=None, help='largest number of effect', required=True)
parser.add_argument('--sigma', type=float, default=1e-2, help='prior probabilities for shared effect')
parser.add_argument('--hess', type=float, default=None, help='HESS estimator will be used as default')
parser.add_argument('--varb', type=float, default=None,
                    help='prior effect size variance derived from HESS will be used')
parser.add_argument('--ptLD', type=float, default=0.2, help='P+T LD cutoff')
parser.add_argument('--ptp', type=float, default=1e-5, help='P+T p value cutoff')
parser.add_argument('--cthres', type=float, default=0.95, help='coverage level for effect group')
parser.add_argument('--ethres', type=float, default=50.0, help='entropy level for effect group')


args = parser.parse_args()
title()
if not os.path.exists(args.save):
    os.makedirs(args.save)
ukb(args)
