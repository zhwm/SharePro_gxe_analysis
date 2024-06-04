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
    bhatlst = [get_bhat(zfile['Beta_Marginal'].values, zfile['SE_Beta_Marginal'].values, zfile['N_Samples'].values)]
    Nlst = [zfile['N_Samples'].values]
    eff, eff_gamma, eff_purity, PIP_overall = adaptive_train(bhatlst, ld, Nlst, args.sigma, args.K, args.maxite, args.eps, args.ubound, args.cthres, args.pthres)
    df_res = pd.DataFrame({'SNP': zfile['SNPID']})
    df_res['PIP'] = PIP_overall
    df_res['cs'] = 0
    for e,val in eff.items():
        mcs_idx = [df_res['SNP'][j] for j in val]
        logging.info(f'The {e}-th effect group contains effective variants:')
        logging.info(f'causal variants: {mcs_idx}')
        logging.info(f'variant probabilities for this effect group: {eff_gamma[e]}')
        logging.info(f'variant purity for this effect group: {eff_purity[e]}\n')
        df_res.loc[val[0], 'cs'] = e+1
        df_res.loc[val[0], 'cs_variants'] = '/'.join(mcs_idx)
        df_res.loc[val[0], 'cs_purity'] = eff_purity[e]
    df_res.to_csv('{}.sharepro.combine.txt'.format(args.save), sep='\t', header=True, index=False)

if __name__ == '__main__':
    args = parse_args()
    logging.basicConfig(level=logging.INFO)
    print_args(args)
    main(args)