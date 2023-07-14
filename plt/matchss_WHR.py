import pandas as pd
import numpy as np
from scipy.stats import chi2

cs = pd.read_csv('../doc/WHR_cgene.txt', sep='\t')


def get_AF_diff(path1, path2, cs):
    female = pd.read_csv(path1, sep='\t')
    female.index = female['CHR'].astype(str) + '.' + female['POS'].astype(str) + '.' + female['A2'] + '.' + female['A1']
    male = pd.read_csv(path2, sep='\t')
    male.index = male['CHR'].astype(str) + '.' + male['POS'].astype(str) + '.' + male['A2'] + '.' + male['A1']
    FE = female.loc[cs['Top_variant'].values,]
    MA = male.loc[cs['Top_variant'].values,]
    cs['AFdiff'] = np.abs(FE['AF1']-MA['AF1']).values
    mmaf = ((FE['AF1'] * FE['N'] + MA['AF1'] * MA['N'])/(FE['N'] + MA['N'])).values
    cs['pAFdiff'] = -np.log10(chi2.sf(cs['AFdiff'] ** 2 / mmaf * (1 - mmaf) / (1/FE['N'].values + 1/MA['N'].values),1))
    return cs


cs = get_AF_diff('../dat/GxE/WHR_FE_BMI_INT.txt.fastGWA', '../dat/GxE/WHR_MA_BMI_INT.txt.fastGWA', cs)
cs['betaF'] = [float(i.split(',')[0]) for i in cs['beta']]
cs['betaM'] = [float(i.split(',')[1]) for i in cs['beta']]


def match_hormone_ss(path, sex, trait, cs):
    rsl = [i.split('/')[0] for i in cs['rsid']]
    M_SHBG = pd.read_csv(path, sep='\t')
    M_SHBG.index = M_SHBG['variant_id']
    cs['b{}_{}'.format(sex, trait)] = M_SHBG.loc[rsl, 'beta'].values
    cs['se{}_{}'.format(sex, trait)] = M_SHBG.loc[rsl, 'standard_error'].values
    cs['mf{}_{}'.format(sex, trait)] = M_SHBG.loc[rsl, 'effect_allele_frequency'].values
    return cs


cs = match_hormone_ss('../dat/GxE/GCST90012108_buildGRCh37.tsv.gz', 'M', 'SHBG', cs)
cs = match_hormone_ss('../dat/GxE/GCST90012106_buildGRCh37.tsv.gz', 'F', 'SHBG', cs)
cs = match_hormone_ss('../dat/GxE/GCST90012103_buildGRCh37.tsv.gz', 'M', 'Testo', cs)
cs = match_hormone_ss('../dat/GxE/GCST90012102_buildGRCh37.tsv.gz', 'F', 'Testo', cs)


def match_fastgwa_ss(path, sex, trait, cs):
    rsl = [i.split('/')[0] for i in cs['rsid']]
    F_TG = pd.read_csv(path, sep='\t')
    F_TG.index = F_TG['SNP']
    cs['b{}_{}'.format(sex, trait)] = F_TG.loc[rsl, 'BETA'].values
    cs['se{}_{}'.format(sex, trait)] = F_TG.loc[rsl, 'SE'].values
    cs['mf{}_{}'.format(sex, trait)] = F_TG.loc[rsl, 'AF1'].values
    return cs


cs = match_fastgwa_ss('../dat/GxE/TG_FE_INT.txt.fastGWA', 'F', 'TG', cs)
cs = match_fastgwa_ss('../dat/GxE/TG_MA_INT.txt.fastGWA', 'M', 'TG', cs)
cs = match_fastgwa_ss('../dat/GxE/LDL_FE_INT.txt.fastGWA', 'F', 'LDL', cs)
cs = match_fastgwa_ss('../dat/GxE/LDL_MA_INT.txt.fastGWA', 'M', 'LDL', cs)
cs = match_fastgwa_ss('../dat/GxE/HDL_MA_INT.txt.fastGWA', 'M', 'HDL', cs)
cs = match_fastgwa_ss('../dat/GxE/HDL_FE_INT.txt.fastGWA', 'F', 'HDL', cs)

cs.to_csv('../doc/WHR_sex_hormone_LDL_TG_HDL.txt', sep='\t', header=True, index=False)

