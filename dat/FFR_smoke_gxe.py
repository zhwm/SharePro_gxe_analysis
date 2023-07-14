import pandas as pd
import numpy as np
from scipy.stats import chi2

c = pd.read_csv('FFR_CS_INT.txt.fastGWA', sep='\t')
e = pd.read_csv('FFR_ES_INT.txt.fastGWA', sep='\t')
n = pd.read_csv('FFR_NS_INT.txt.fastGWA', sep='\t')
c.index = c['CHR'].astype(str) + '.' + c['POS'].astype(str) + '.' + c['A2'] + '.' + c['A1']
e.index = e['CHR'].astype(str) + '.' + e['POS'].astype(str) + '.' + e['A2'] + '.' + e['A1']
n.index = n['CHR'].astype(str) + '.' + n['POS'].astype(str) + '.' + n['A2'] + '.' + n['A1']
fl = pd.read_csv('FFR.lst', sep='\t', header=None)
fl.columns = ['SNP']
fl.index = fl['SNP']
fl[['CHR', 'POS', 'A1', 'A2']] = c.loc[fl['SNP'], ['CHR', 'POS', 'A1', 'A2']].values
fl['betaCE'] = c.loc[fl['SNP'], 'BETA'] * (np.sqrt(2 * (1 - c.loc[fl['SNP'], 'AF1']) * c.loc[fl['SNP'], 'AF1'])) - \
               e.loc[fl['SNP'], 'BETA'] * (np.sqrt(2 * (1 - e.loc[fl['SNP'], 'AF1']) * e.loc[fl['SNP'], 'AF1']))
fl['betaCN'] = c.loc[fl['SNP'], 'BETA'] * (np.sqrt(2 * (1 - c.loc[fl['SNP'], 'AF1']) * c.loc[fl['SNP'], 'AF1'])) - \
               n.loc[fl['SNP'], 'BETA'] * (np.sqrt(2 * (1 - n.loc[fl['SNP'], 'AF1']) * n.loc[fl['SNP'], 'AF1']))
fl['betaEN'] = e.loc[fl['SNP'], 'BETA'] * (np.sqrt(2 * (1 - e.loc[fl['SNP'], 'AF1']) * e.loc[fl['SNP'], 'AF1'])) - \
               n.loc[fl['SNP'], 'BETA'] * (np.sqrt(2 * (1 - n.loc[fl['SNP'], 'AF1']) * n.loc[fl['SNP'], 'AF1']))
fl['pCE'] = chi2.sf((c.loc[fl['SNP'], 'BETA'] - e.loc[fl['SNP'], 'BETA']) ** 2 / (
            c.loc[fl['SNP'], 'SE'] ** 2 + e.loc[fl['SNP'], 'SE'] ** 2), 1)
fl['pCN'] = chi2.sf((c.loc[fl['SNP'], 'BETA'] - n.loc[fl['SNP'], 'BETA']) ** 2 / (
            c.loc[fl['SNP'], 'SE'] ** 2 + n.loc[fl['SNP'], 'SE'] ** 2), 1)
fl['pEN'] = chi2.sf((e.loc[fl['SNP'], 'BETA'] - n.loc[fl['SNP'], 'BETA']) ** 2 / (
            e.loc[fl['SNP'], 'SE'] ** 2 + n.loc[fl['SNP'], 'SE'] ** 2), 1)
fl.to_csv('FFR_smoke_gxesnp.txt', sep='\t', header=True, index=False)
