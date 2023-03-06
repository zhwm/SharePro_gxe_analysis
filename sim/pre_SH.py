import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='Prepare for SharePro simulation')
parser.add_argument('--dir', type=str, default=None, help='path to GWAS result')
parser.add_argument('--loci', type=str, default=None, help='name of loci')
args = parser.parse_args()

for i in range(1, 51):
    CLss = pd.read_csv(os.path.join(args.dir, 'CL{}.qassoc.gxe'.format(i)), sep='\s+')
    CMss = pd.read_csv(os.path.join(args.dir, 'CM{}.qassoc.gxe'.format(i)), sep='\s+')
    CSss = pd.read_csv(os.path.join(args.dir, 'CS{}.qassoc.gxe'.format(i)), sep='\s+')
    ACLss = pd.read_csv(os.path.join(args.dir, 'CL{}'.format(i)), sep='\s+')
    ACMss = pd.read_csv(os.path.join(args.dir, 'CM{}'.format(i)), sep='\s+')
    ACSss = pd.read_csv(os.path.join(args.dir, 'CS{}'.format(i)), sep='\s+')
    CLss['CZ'] = (CLss['BETA1'] / CLss['SE1']).round(4)
    CLss['LZ'] = (CLss['BETA2'] / CLss['SE2']).round(4)
    CMss['MZ'] = (CMss['BETA2'] / CMss['SE2']).round(4)
    CSss['SZ'] = (CSss['BETA2'] / CSss['SE2']).round(4)
    ACLss['Z'] = (ACLss['Beta_Marginal'] / ACLss['SE_Beta_Marginal']).round(4)
    ACMss['Z'] = (ACMss['Beta_Marginal'] / ACMss['SE_Beta_Marginal']).round(4)
    ACSss['Z'] = (ACSss['Beta_Marginal'] / ACSss['SE_Beta_Marginal']).round(4)
    CLss[['SNP', 'CZ']].to_csv(os.path.join(args.dir, 'C{}.z'.format(i)), sep='\t', index=False, header=False)
    CLss[['SNP', 'LZ']].to_csv(os.path.join(args.dir, 'L{}.z'.format(i)), sep='\t', index=False, header=False)
    CMss[['SNP', 'MZ']].to_csv(os.path.join(args.dir, 'M{}.z'.format(i)), sep='\t', index=False, header=False)
    CSss[['SNP', 'SZ']].to_csv(os.path.join(args.dir, 'S{}.z'.format(i)), sep='\t', index=False, header=False)
    ACLss[['SNPID', 'Z']].to_csv(os.path.join(args.dir, 'ACL{}.z'.format(i)), sep='\t', index=False, header=False)
    ACMss[['SNPID', 'Z']].to_csv(os.path.join(args.dir, 'ACM{}.z'.format(i)), sep='\t', index=False, header=False)
    ACSss[['SNPID', 'Z']].to_csv(os.path.join(args.dir, 'ACS{}.z'.format(i)), sep='\t', index=False, header=False)

SHzld = pd.DataFrame({'z': ['C{}.z,L{}.z'.format(i, i) for i in range(1, 51)]})  # zld file for sharepro
SHzld['ld'] = '../{}.ld,../{}.ld'.format(args.loci, args.loci)
SHzld.to_csv(os.path.join(args.dir, 'CL.zld'), sep='\t', index=False, header=True)

SHzld = pd.DataFrame({'z': ['C{}.z,M{}.z'.format(i, i) for i in range(1, 51)]})  # zld file for sharepro
SHzld['ld'] = "../{}.ld,../{}.ld".format(args.loci, args.loci)
SHzld.to_csv(os.path.join(args.dir, 'CM.zld'), sep='\t', index=False, header=True)

SHzld = pd.DataFrame({'z': ['C{}.z,S{}.z'.format(i, i) for i in range(1, 51)]})  # zld file for sharepro
SHzld['ld'] = '../{}.ld,../{}.ld'.format(args.loci, args.loci)
SHzld.to_csv(os.path.join(args.dir, 'CS.zld'), sep='\t', index=False, header=True)

SPzld = pd.DataFrame({'z': ['C{}.z'.format(i) for i in range(1, 51)]})
SPzld['ld'] = '../{}.ld'.format(args.loci)
SPzld.to_csv(os.path.join(args.dir, 'SPC.zld'), sep='\t', index=False, header=True)

SPzld = pd.DataFrame({'z': ['ACL{}.z'.format(i) for i in range(1, 51)]})
SPzld['ld'] = '../{}.ld'.format(args.loci)
SPzld.to_csv(os.path.join(args.dir, 'ACL.zld'), sep='\t', index=False, header=True)

SPzld = pd.DataFrame({'z': ['ACM{}.z'.format(i) for i in range(1, 51)]})
SPzld['ld'] = '../{}.ld'.format(args.loci)
SPzld.to_csv(os.path.join(args.dir, 'ACM.zld'), sep='\t', index=False, header=True)

SPzld = pd.DataFrame({'z': ['ACS{}.z'.format(i) for i in range(1, 51)]})
SPzld['ld'] = '../{}.ld'.format(args.loci)
SPzld.to_csv(os.path.join(args.dir, 'ACS.zld'), sep='\t', index=False, header=True)
