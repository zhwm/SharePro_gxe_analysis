import sys
import argparse
import logging
import pandas as pd
import numpy as np
from scipy.stats import chi2
sys.path.append('/home/wmzh22/scratch/SharePro_gxe/src/')
from SharePro.sharepro_gxe import SharePro,SparseReg,get_bhat,adaptive_train

def parse_args():
    parser = argparse.ArgumentParser(description='SharePro Commands:')
    parser.add_argument('--z', type=str, default=None, help='path to matched summary statistics', required=True)
    parser.add_argument('--ld', type=str, default=None, help='path to matched ld matrix', required=True)
    parser.add_argument('--save', type=str, default=None, help='path to save results', required=True)
    parser.add_argument('--sigma', type=float, default=0.99, help='prior shared prob')
    parser.add_argument('--K', type=int, default=10, help='largest number of causal signals')
    parser.add_argument('--maxite', type=int, default=100, help='max number of iterations')
    parser.add_argument('--eps', type=float, default=1e-2, help='convergence criterion')
    parser.add_argument('--ubound', type=int, default=1e10, help='upper bound for inconvergence')
    parser.add_argument('--cthres', type=float, default=0.95, help='attainable coverage threshold for effect groups')
    parser.add_argument('--pthres', type=float, default=0.8, help='purity threshold for effect groups')
    args = parser.parse_args()
    return args

def print_args(args):
    for arg in vars(args):
        logging.info(f"{arg}: {getattr(args, arg)}")

def main(args):
    zfile = pd.read_csv(args.z, sep='\s+')
    ld = pd.read_csv(args.ld, sep='\s+', header=None).values
    bhatlst = [get_bhat(zfile['BETA1'].values, zfile['SE1'].values, zfile['NMISS1'].values),
        get_bhat(zfile['BETA2'].values, zfile['SE2'].values, zfile['NMISS2'].values)]
    Nlst = [zfile['NMISS1'].values, zfile['NMISS2'].values]
    eff, eff_gamma, eff_purity, PIP_overall = adaptive_train(bhatlst, ld, Nlst, args.sigma, args.K, args.maxite, args.eps, args.ubound, args.cthres, args.pthres)
    df_res = pd.DataFrame({'SNP': zfile['SNP']})
    df_res['PIP'] = PIP_overall
    df_res['cs'] = 0
    df_res['cs_variants'] = 'NA'
    df_res['cs_purity'] = 0
    df_res['nlog10pdiff'] = 0
    cstop = []
    pdiff = np.nan_to_num(chi2.logsf((zfile['BETA1'] - zfile['BETA2'])**2 / (zfile['SE1']**2 + zfile['SE2']**2), 1) / -np.log(10))
    df_res['nlog10pGxE'] = ['{:.4f}'.format(val) for val in pdiff]
    for e,val in eff.items():
        mcs_idx = [zfile['SNP'][j] for j in val]
        logging.info(f'The {e}-th effect group contains effective variants:')
        logging.info(f'causal variants: {mcs_idx}')
        logging.info(f'variant probabilities for this effect group: {eff_gamma[e]}')
        logging.info(f'variant purity for this effect group: {eff_purity[e]}')
        #logging.info('negative log GxE p-value for this effect group: {}\n'.format(df_res.loc[val[0], "nlog10pGxE"].values))
        df_res.loc[val[0], 'cs'] = e+1
        cstop.append(val[0])
        df_res.loc[val[0], 'cs_variants'] = '/'.join(mcs_idx)
        df_res.loc[val[0], 'cs_purity'] = eff_purity[e]
        df_res.loc[val[0], 'nlog10pdiff'] = df_res.loc[val[0],'nlog10pGxE']
    df_res.to_csv('{}.sharepro.gxe.txt'.format(args.save), sep='\t', header=True, index=False)

if __name__ == '__main__':
    args = parse_args()
    logging.basicConfig(filename='{}.sharepro.gxe.log'.format(args.save), level=logging.INFO, filemode='w', format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    print_args(args)
    main(args)








