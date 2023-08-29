import math, random, warnings, copy, time
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from scipy.special import erf
import sklearn.linear_model as lm
import sklearn.model_selection as mod_sel
from sklearn.metrics import roc_auc_score, roc_curve, auc
from sklearn.decomposition import PCA
from sklearn import cluster, mixture
import matplotlib.pyplot as plt
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neighbors import kneighbors_graph
from sklearn.neighbors import NearestNeighbors
from sklearn.manifold import TSNE
from sklearn.exceptions import ConvergenceWarning
import umap
from marker_count import MkrCnt_Ref
#from MarkerCount.marker_count import MkrCnt_Ref

CLUSTERING_AGO = 'lv'
try:
    from sknetwork.clustering import Louvain
except ImportError:
    print('WARNING: sknetwork not installed. Will used GMM for clustering.')
    CLUSTERING_AGO = 'gmm'

SEABORN = True
try:
    import seaborn as sns
except ImportError:
    SEABORN = False
    print('WARNING: seaborn not installed. Install seaborn to check out the summary.')    
    
MHC1_prefix_lst = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G']
MHC2_prefix_lst = ['HLA-DM', 'HLA-DO', 'HLA-DP', 'HLA-DQ']
HLA_DR_prefix_lst = ['HLA-DR']
PNSH12 = '101000'
LOGREG_TOL = 1e-3
LOGREG_MAX_ITER = 1000
LOGREG_NCV = 5
LOGREG_PARAM_GRID = { 'C': [0.5,1,2,4,8,16], \
                   'l1_ratio': [0.2, 0.5, 0.9, 1] } 
CLUSTER_BASIS_CORECTION_DIST_SCALE = [2]
CLUSTER_BASIS_CORECTION_N_NEIGHBORS = 3
CLUSTER_BASIS_CORECTION_MIN_PCT_TO_INVALIDATE = 0.00
SEPARABILITY_THRESHOLD = 0.95
SEPARABILITY_MIN_NUM_CELLS = 100
SEPARABILITY_AUC_INIT_VALUE = 2
# PTH_CUTOFF_MARGIN_MULTIPLIER = 3
PCT_THRESHOLD_MAX = 0.6
PCT_THRESHOLD_MIN = 0.2


def check_start_with_a_key(s, key):
    if s[:len(key)] == key:
        return True
    else:
        return False

def check_start_with(s, keys):
    for key in keys:
        r = check_start_with_a_key(s, key)
        if r:
            break
    return r
 
def get_markers_from_df(df, target_cell_lst, pnsh12 = PNSH12, verbose = True):
    
    pos = bool(int(pnsh12[0]))
    neg = bool(int(pnsh12[1]))
    sec = bool(int(pnsh12[2]))
    hla_dr = bool(int(pnsh12[3]))
    mhc1 = bool(int(pnsh12[4]))
    mhc2 = bool(int(pnsh12[5]))

    mkr_lst_dict = {}
    for target_cell in target_cell_lst:
        dfs = df.loc[df['cell_type_minor'] == target_cell, :]

        cell_type_minor_lst =  list(set(dfs['cell_type_subset'].unique()))
        cell_type_minor_lst.sort()
        cell_type_minor_lst

        exp_lst = []
        if pos: exp_lst.append('pos')
        if neg: exp_lst.append('neg')
        if sec: exp_lst.append('sec')

        for c in cell_type_minor_lst:
            mkr_lst = []
            for e in exp_lst:
                b = (dfs['cell_type_subset'] == c) & (dfs['exp'] == e)
                if np.sum(b) == 1:
                    idx = dfs.index.values[b][0]
                    mkrs = []
                    items = dfs.loc[idx,'markers'].split(',')
                    for item in items:
                        mkrs.append(item.strip())            
                    mkr_lst = mkr_lst + mkrs
                elif np.sum(b) > 1:
                    print('ERROR: ----')

            if not mhc1:
                mkr_lst2 = []
                for mkr in mkr_lst:
                    if not check_start_with(mkr, MHC1_prefix_lst):
                        mkr_lst2.append(mkr)
                mkr_lst = copy.deepcopy(mkr_lst2)
                
            if not hla_dr:
                mkr_lst2 = []
                for mkr in mkr_lst:
                    if not check_start_with(mkr, HLA_DR_prefix_lst):
                        mkr_lst2.append(mkr)
                mkr_lst = copy.deepcopy(mkr_lst2)

            if not mhc2:
                mkr_lst2 = []
                for mkr in mkr_lst:
                    if not check_start_with(mkr, MHC2_prefix_lst):
                        mkr_lst2.append(mkr)
                mkr_lst = copy.deepcopy(mkr_lst2)

            mkr_lst_dict[c] = list(set(mkr_lst))

        # if verbose: print_mkrs(mkr_lst_dict) 
        
    return mkr_lst_dict


def remove_common( mkr_dict, prn = True ):

    cts = list(mkr_dict.keys())
    mkrs_all = []
    for c in cts:
        mkrs_all = mkrs_all + mkr_dict[c]
    mkrs_all = list(set(mkrs_all))
    df = pd.DataFrame(index = mkrs_all, columns = cts)
    df.loc[:,:] = 0

    for c in cts:
        df.loc[mkr_dict[c], c] = 1
    Sum = df.sum(axis = 1)
    
    to_del = []
    s = ''
    for c in cts:
        b = (df[c] > 0) & (Sum == 1)
        mkrs1 = list(df.index.values[b])
        if prn & (len(mkr_dict[c]) != len(mkrs1)):
            s = s + '%s: %i > %i, ' % (c, len(mkr_dict[c]), len(mkrs1))
        
        if len(mkrs1) == 0:
            to_del.append(c)
        else:
            mkr_dict[c] = mkrs1

    if prn & len(s) > 0:
        print(s[:-2])

    if len(to_del) > 0:
        for c in cts:
            if c in to_del:
                del mkr_dict[c]
                
    return mkr_dict


def remove_common3( mkr_dict, prn = True ):

    cts = list(mkr_dict.keys())
    to_del = []
    s = ''
    for c in cts:
        mkrs1 = mkr_dict[c]
        mkrs2 = []
        for c2 in cts:
            if c2 != c:
                mkrs2 = mkrs2 + mkr_dict[c2]
        mkrs2 = list(set(mkrs2))
        mkrs1 = list(set(mkrs1) - set(mkrs2))

        if prn & (len(mkr_dict[c]) != len(mkrs1)):
            s = s + '%s: %i > %i, ' % (c, len(mkr_dict[c]), len(mkrs1))
        
        if len(mkrs1) == 0:
            to_del.append(c)
        else:
            mkr_dict[c] = mkrs1

    if prn & len(s) > 0:
        print(s[:-2])

    if len(to_del) > 0:
        for c in cts:
            if c in to_del:
                del mkr_dict[c]
                
    return mkr_dict

        
def remove_common2( mkr_dict, prn = True ):

    cts = list(mkr_dict.keys())
    to_del = []
    s = ''

    mkr_common = mkr_dict[cts[0]]
    for c in cts:
        mkr_common = list(set( mkr_dict[c] ).intersection(mkr_common) )
    
    for c in cts:
        mkrs1 = mkr_dict[c]
        mkrs1 = list(set(mkrs1) - set(mkr_common))

        if prn & (len(mkr_dict[c]) != len(mkrs1)):
            s = s + '%s: %i > %i, ' % (c, len(mkr_dict[c]), len(mkrs1))
        
        if len(mkrs1) == 0:
            to_del.append(c)
        else:
            mkr_dict[c] = mkrs1

    if prn & len(s) > 0:
        print(s[:-2])

    if len(to_del) > 0:
        for c in cts:
            if c in to_del:
                del mkr_dict[c]
                
    return mkr_dict

def load_marker_file( file ):

    df = pd.read_csv(file, sep = '\t')

    b = ~df['markers'].isnull()
    mkr_lst = list(df.loc[b, 'markers'])
    mkr_lst_new = []
    for s in mkr_lst:
        mkr_lst_new.append(s.upper())
    df.loc[b, 'markers'] = mkr_lst_new

    return df
        
def get_markers_major_type(file, target_cells = [], pnsh12 = PNSH12,
                                                rem_common = True, verbose = False):
    
    if verbose: print('Load markers .. ', end = '', flush = True)
    df = load_marker_file( file )

    if len(target_cells) == 0:
        target_cells = list(df['cell_type_major'].unique())
    
    major_type_lst = list(df['cell_type_major'].unique())
        
    mkr_lst = {}
    mkr_lst_neg = {}
    for c in target_cells:
        if c in major_type_lst:
            b = df['cell_type_major'] == c
            cell_type_lst = list(df.loc[b, 'cell_type_minor'].unique())
            pnsh12_t = '%s0%s' % (pnsh12[0], pnsh12[2:])
            mkrs = get_markers_from_df(df, cell_type_lst, pnsh12 = pnsh12_t, verbose = verbose)
            mkr_c = []
            for key in list(mkrs.keys()):
                mkr_c = mkr_c + mkrs[key]
            mkr_lst[c] = list(set(mkr_c))
        
            mkr_c_neg = []
            if pnsh12[1] == '1':
                pnsh12_n = '010%s' % (pnsh12[3:])
                mkrs_neg = get_markers_from_df(df, cell_type_lst, pnsh12 = pnsh12_n, verbose = verbose)
                for key in list(mkrs_neg.keys()):
                    mkr_c_neg = list(set(mkr_c_neg).intersection(mkrs_neg[key]))
            mkr_lst_neg[c] = list(set(mkr_c_neg))
            
    if rem_common:
        mkr_lst = remove_common( mkr_lst, prn = verbose )
        mkr_lst_neg = remove_common( mkr_lst_neg, prn = verbose )
        
    sm = ''
    cell_types = list(mkr_lst.keys())
    cell_types.sort()
    for key in cell_types:
        sm = sm + '%s,' % key

    if verbose: print(' %i types. \n%s' % (len(mkr_lst.keys()), sm[:-1]))
        
    return mkr_lst, mkr_lst_neg

def get_markers_cell_type(file, target_cells = [], pnsh12 = PNSH12,
                          rem_common = True, verbose = False):
    
    if verbose: print('Load markers .. ', end = '', flush = True)
    df = load_marker_file( file )

    if len(target_cells) == 0:
        target_cells = list(df['cell_type_major'].unique())
    
    major_type_lst = list(df['cell_type_major'].unique())
        
    mkr_lst = {}
    mkr_lst_neg = {}
    pnsh12_t = '%s0%s' % (pnsh12[0], pnsh12[2:])
    for c in target_cells:
        if c in major_type_lst:
            b = df['cell_type_major'] == c
            cell_type_lst = list(df.loc[b, 'cell_type_minor'].unique())
            for c2 in cell_type_lst:

                mkrs = get_markers_from_df(df, [c2], pnsh12 = pnsh12_t, verbose = verbose)
                if 'common' in list(mkrs.keys()):
                    mkr_lst[c2] = mkrs['common']
                else:
                    mkr_c = []
                    for key in list(mkrs.keys()):
                        mkr_c = mkr_c + mkrs[key]
                    mkr_lst[c2] = list(set(mkr_c))
    
                mkr_c_neg = []
                if pnsh12[1] == '1':
                    pnsh12_n = '010%s' % (pnsh12[3:])
                    mkrs_neg = get_markers_from_df(df, [c2], pnsh12 = pnsh12_n, verbose = verbose)
                    if 'common' in list(mkrs_neg.keys()):
                        mkr_c_neg = mkrs_neg['common']
                    else:
                        cnt = 0
                        for key in list(mkrs_neg.keys()):
                            if cnt == 0:
                                mkr_c_neg = copy.deepcopy(mkrs_neg[key])
                            else:
                                mkr_c_neg = list(set(mkr_c_neg).intersection(mkrs_neg[key]))
                            cnt += 1
                mkr_lst_neg[c2] = list(set(mkr_c_neg))

    if rem_common:
        mkr_lst = remove_common( mkr_lst, prn = verbose )
        mkr_lst_neg = remove_common( mkr_lst_neg, prn = verbose )
        
    sm = ''
    cell_types = list(mkr_lst.keys())
    cell_types.sort()
    for key in cell_types:
        sm = sm + '%s,' % key
        
    if verbose: print(' %i types. \n%s' % (len(mkr_lst.keys()), sm[:-1]))
        
    return mkr_lst, mkr_lst_neg


def get_markers_minor_type(file, target_cells = [], pnsh12 = PNSH12, verbose = False):
    
    if verbose: print('Load markers .. ', end = '', flush = True)
    df = load_marker_file( file )
    
    if len(target_cells) == 0:
        target_cells = list(df['cell_type_major'].unique())
    
    major_type_lst = list(df['cell_type_major'].unique())
    
    mkr_lst = {}
    for c in target_cells:
        if c in major_type_lst:
            b = df['cell_type_major'] == c
            cell_type_lst = list(df.loc[b, 'cell_type'].unique())
            mkrs = get_markers_from_df(df, cell_type_lst, pnsh12 = pnsh12, verbose = verbose)
            mkr_lst.update(mkrs)
       
    Keys = mkr_lst.keys()
    if ('NCR+ ILC3' in Keys) & ('NCR- ILC3' in Keys):
        mkr_lst['ILC3'] = mkr_lst['NCR+ ILC3']
        mkr_lst['ILC3'] = mkr_lst['ILC3'] + mkr_lst['NCR- ILC3']
        mkr_lst['ILC3'] = list(set(mkr_lst['ILC3']))
        del mkr_lst['NCR+ ILC3']
        del mkr_lst['NCR- ILC3']
        
    if ('MDSC (Granulocytic)' in Keys) & ('MDSC (Monocytic)' in Keys):
        mkr_lst['MDSC'] = mkr_lst['MDSC (Monocytic)']
        mkr_lst['MDSC'] = mkr_lst['MDSC'] + mkr_lst['MDSC (Granulocytic)']
        mkr_lst['MDSC'] = list(set(mkr_lst['MDSC']))
        del mkr_lst['MDSC (Granulocytic)']
        del mkr_lst['MDSC (Monocytic)']
        
    if ('Macrophage (M2A)' in Keys) & ('Macrophage (M2B)' in Keys) & \
        ('Macrophage (M2C)' in Keys) & ('Macrophage (M2D)' in Keys):
        mkr_lst['Macrophage (M2)'] = mkr_lst['Macrophage (M2A)']
        mkr_lst['Macrophage (M2)'] = mkr_lst['Macrophage (M2)'] + mkr_lst['Macrophage (M2B)']
        mkr_lst['Macrophage (M2)'] = mkr_lst['Macrophage (M2)'] + mkr_lst['Macrophage (M2C)']
        mkr_lst['Macrophage (M2)'] = mkr_lst['Macrophage (M2)'] + mkr_lst['Macrophage (M2D)']
        mkr_lst['Macrophage (M2)'] = list(set(mkr_lst['Macrophage (M2)']))
        del mkr_lst['Macrophage (M2A)']
        del mkr_lst['Macrophage (M2B)']
        del mkr_lst['Macrophage (M2C)']
        del mkr_lst['Macrophage (M2D)']
        
    if ('Macrophage' in Keys) & ('Macrophage (M1)' in Keys) & ('Macrophage (M2)' in Keys):
        mkr_lst['Macrophage (M2)'] = mkr_lst['Macrophage (M2)'] + mkr_lst['Macrophage']
        mkr_lst['Macrophage (M2)'] = list(set(mkr_lst['Macrophage (M2)']))
        mkr_lst['Macrophage (M1)'] = mkr_lst['Macrophage (M1)'] + mkr_lst['Macrophage']
        mkr_lst['Macrophage (M1)'] = list(set(mkr_lst['Macrophage (M1)']))
        del mkr_lst['Macrophage']

    sm = ''
    cell_types = list(mkr_lst.keys())
    cell_types.sort()
    for key in cell_types:
        sm = sm + '%s,' % key
        
    if verbose: print(' %i types. \n%s' % (len(mkr_lst.keys()), sm[:-1]))
        
    return mkr_lst


def get_cell_type_dict(file):
    
    df = load_marker_file( file )
    
    target_cells = list(df['cell_type_major'].unique())    
        
    cell_type_dict = {}
    for c in target_cells:
        b = df['cell_type_major'] == c
        cell_type_lst = list(df.loc[b, 'cell_type_minor'].unique())
        for c2 in cell_type_lst:
            cell_type_dict[c2] = c
        
    return cell_type_dict

def get_minor_type_dict(file):
    
    df = load_marker_file( file )
    
    target_cells = list(df['cell_type_minor'].unique())    
        
    cell_type_dict = {}
    for c in target_cells:
        b = df['cell_type_minor'] == c
        cell_type_lst = list(df.loc[b, 'cell_type_subset'].unique())
        for c2 in cell_type_lst:
            cell_type_dict[c2] = c
        
    Keys = cell_type_dict.keys()
    if ('ILC3 (NCR+)' in Keys) & ('ILC3 (NCR-)' in Keys):
        cell_type_dict['ILC3'] = 'ILC'
        
    if ('MDSC (Granulocytic)' in Keys) & ('MDSC (Monocytic)' in Keys):
        cell_type_dict['MDSC'] = 'MDSC'
        
    if ('Macrophage (M2A)' in Keys) & ('Macrophage (M2B)' in Keys) & \
        ('Macrophage (M2C)' in Keys) & ('Macrophage (M2D)' in Keys):
        cell_type_dict['Macrophage (M2)'] = 'Macrophage'
    
    return cell_type_dict


def comb_markers(mkr_lst_in, maj_dict = None, min_dict = None):
    
    mkr_lst = copy.deepcopy(mkr_lst_in)
    Keys = mkr_lst.keys()
    if ('ILC3 (NCR+)' in Keys) & ('ILC3 (NCR-)' in Keys):
        mkr_lst['ILC3'] = mkr_lst['ILC3 (NCR+)']
        mkr_lst['ILC3'] = mkr_lst['ILC3'] + mkr_lst['ILC3 (NCR-)']
        mkr_lst['ILC3'] = list(set(mkr_lst['ILC3']))
        del mkr_lst['ILC3 (NCR+)']
        del mkr_lst['ILC3 (NCR-)']
        
        if maj_dict is not None:
            maj_dict['ILC3'] = maj_dict['ILC3 (NCR+)']
            del maj_dict['ILC3 (NCR+)']
            del maj_dict['ILC3 (NCR-)']
        
        if min_dict is not None:
            min_dict['ILC3'] = min_dict['ILC3 (NCR+)']
            del min_dict['ILC3 (NCR+)']
            del min_dict['ILC3 (NCR-)']
        
    '''
    if ('MDSC (Granulocytic)' in Keys) & ('MDSC (Monocytic)' in Keys):
        mkr_lst['MDSC'] = mkr_lst['MDSC (Monocytic)']
        mkr_lst['MDSC'] = mkr_lst['MDSC'] + mkr_lst['MDSC (Granulocytic)']
        mkr_lst['MDSC'] = list(set(mkr_lst['MDSC']))
        del mkr_lst['MDSC (Granulocytic)']
        del mkr_lst['MDSC (Monocytic)']
        
        if maj_dict is not None:
            maj_dict['MDSC'] = maj_dict['MDSC (Monocytic)']
            del maj_dict['MDSC (Granulocytic)']
            del maj_dict['MDSC (Monocytic)']
        
        if min_dict is not None:
            min_dict['MDSC'] = min_dict['MDSC (Monocytic)']
            del min_dict['MDSC (Granulocytic)']
            del min_dict['MDSC (Monocytic)']
    '''
    
    if ('Macrophage (M2A)' in Keys) & ('Macrophage (M2B)' in Keys) & \
        ('Macrophage (M2C)' in Keys) & ('Macrophage (M2D)' in Keys):
        mkr_lst['Macrophage (M2)'] = mkr_lst['Macrophage (M2A)']
        mkr_lst['Macrophage (M2)'] = mkr_lst['Macrophage (M2)'] + mkr_lst['Macrophage (M2B)']
        mkr_lst['Macrophage (M2)'] = mkr_lst['Macrophage (M2)'] + mkr_lst['Macrophage (M2C)']
        mkr_lst['Macrophage (M2)'] = mkr_lst['Macrophage (M2)'] + mkr_lst['Macrophage (M2D)']
        mkr_lst['Macrophage (M2)'] = list(set(mkr_lst['Macrophage (M2)']))
        del mkr_lst['Macrophage (M2A)']
        del mkr_lst['Macrophage (M2B)']
        del mkr_lst['Macrophage (M2C)']
        del mkr_lst['Macrophage (M2D)']
        
        if maj_dict is not None:
            maj_dict['Macrophage (M2)'] = maj_dict['Macrophage (M2A)']
            del maj_dict['Macrophage (M2A)']
            del maj_dict['Macrophage (M2B)']
            del maj_dict['Macrophage (M2C)']
            del maj_dict['Macrophage (M2D)']
        
        if min_dict is not None:
            min_dict['Macrophage (M2)'] = min_dict['Macrophage (M2A)']
            del min_dict['Macrophage (M2A)']
            del min_dict['Macrophage (M2B)']
            del min_dict['Macrophage (M2C)']
            del min_dict['Macrophage (M2D)']
        
    if ('Macrophage' in Keys) & ('Macrophage (M1)' in Keys) & ('Macrophage (M2)' in Keys):
        mkr_lst['Macrophage (M2)'] = mkr_lst['Macrophage (M2)'] + mkr_lst['Macrophage']
        mkr_lst['Macrophage (M2)'] = list(set(mkr_lst['Macrophage (M2)']))
        mkr_lst['Macrophage (M1)'] = mkr_lst['Macrophage (M1)'] + mkr_lst['Macrophage']
        mkr_lst['Macrophage (M1)'] = list(set(mkr_lst['Macrophage (M1)']))
        del mkr_lst['Macrophage']

        if maj_dict is not None:
            del maj_dict['Macrophage']
        
        if min_dict is not None:
            del min_dict['Macrophage']
            
    elif ('common' in Keys) & ('Macrophage (M1)' in Keys) & ('Macrophage (M2)' in Keys):       
        mkr_lst['Macrophage (M2)'] = mkr_lst['Macrophage (M2)'] + mkr_lst['common']
        mkr_lst['Macrophage (M2)'] = list(set(mkr_lst['Macrophage (M2)']))
        mkr_lst['Macrophage (M1)'] = mkr_lst['Macrophage (M1)'] + mkr_lst['common']
        mkr_lst['Macrophage (M1)'] = list(set(mkr_lst['Macrophage (M1)']))
        del mkr_lst['common']

        if maj_dict is not None:
            del maj_dict['common']
        
        if min_dict is not None:
            del min_dict['common']
        
    if (maj_dict is not None) & (min_dict is not None):
        return mkr_lst, maj_dict, min_dict
    else:
        return mkr_lst


def get_markers_minor_type2(file, target_cells = [], pnsh12 = PNSH12,
                            rem_common = True, verbose = False):
    
    if verbose: print('Load markers .. ', end = '', flush = True)
    df = load_marker_file( file )
    
    if len(target_cells) == 0:
        target_cells = list(df['cell_type_minor'].unique())
    
    major_type_lst = list(df['cell_type_minor'].unique())
    
    mkr_lst = {}
    mkr_lst_neg = {}
    pnsh12_t = '%s0%s' % (pnsh12[0], pnsh12[2:])
    for c in target_cells:
        if c in major_type_lst:
            b = df['cell_type_minor'] == c
            cell_type_lst = list(df.loc[b, 'cell_type_minor'].unique())
            mkrs = get_markers_from_df(df, cell_type_lst, pnsh12 = pnsh12_t, verbose = verbose)
            mkr_lst.update(mkrs)
        
            if pnsh12[1] == '1':
                pnsh12_n = '010%s' % (pnsh12[3:])
                mkrs_neg = get_markers_from_df(df, cell_type_lst, pnsh12 = pnsh12_n, verbose = verbose)
                mkr_lst_neg.update(mkrs_neg)
        
    if len(mkr_lst.keys()) == 0:
        return mkr_lst, mkr_lst_neg
       
    mkr_lst = comb_markers(mkr_lst)
    mkr_lst_neg = comb_markers(mkr_lst_neg)
    
    if rem_common:
        mkr_lst = remove_common( mkr_lst, prn = verbose )
        mkr_lst_neg = remove_common( mkr_lst_neg, prn = verbose )
        
    sm = ''
    cell_types = list(mkr_lst.keys())
    cell_types.sort()
    for key in cell_types:
        sm = sm + '%s,' % key
        
    if verbose: print(' %i types. \n%s' % (len(mkr_lst.keys()), sm[:-1]))
        
    return mkr_lst, mkr_lst_neg


def print_mkrs(mkr_lst_dict):
    
    for key in list(mkr_lst_dict.keys()):
        if len(mkr_lst_dict[key]) > 0:
            s = mkr_lst_dict[key][0]
            for mkr in mkr_lst_dict[key][1:]:
                s = s + ',%s' % mkr
            print('%s (%i): %s' % (key, len(mkr_lst_dict[key]), s))
        else:
            print('%s (%i): ' % (key, len(mkr_lst_dict[key])))

            
def load_markers_all(file, target_cells = [], pnsh12 = PNSH12, verbose = False):
    
    if verbose: print('Load markers .. ', end = '', flush = True)
    df = load_marker_file( file )
    
    if len(target_cells) == 0:
        target_cells = list(df['cell_type_major'].unique())
    
    major_type_lst = list(df['cell_type_major'].unique())
    
    mkr_lst = {}
    mkr_lst_neg = {}
    major_dict = {}
    minor_dict = {}
    pnsh12_t = '%s0%s' % (pnsh12[0], pnsh12[2:])
    for c in target_cells:
        if c in major_type_lst:
            b = df['cell_type_major'] == c
            cell_type_lst = list(df.loc[b, 'cell_type_minor'].unique())
            for c2 in cell_type_lst:
                mkrs = get_markers_from_df(df, [c2], pnsh12 = pnsh12_t, verbose = verbose)
                mkr_lst.update(mkrs)

                if pnsh12[1] == '1':
                    pnsh12_n = '010%s' % (pnsh12[3:])
                    mkrs_neg = get_markers_from_df(df, [c2], pnsh12 =pnsh12_n, verbose = verbose)
                    mkr_lst_neg.update(mkrs_neg)
                
                for key in mkrs.keys():
                    minor_dict[key] = c2
                    major_dict[key] = c
            
    if len(mkr_lst.keys()) == 0:
        return mkr_lst, mkr_lst_neg
       
    mkr_lst, major_dict, minor_dict = comb_markers(mkr_lst, major_dict, minor_dict)
    mkr_lst_neg = comb_markers(mkr_lst_neg)

    sm = ''
    cell_types = list(mkr_lst.keys())
    cell_types.sort()
    for key in cell_types:
        sm = sm + '%s,' % key
        
    if verbose: print(' %i types. \n%s' % (len(mkr_lst.keys()), sm[:-1]))
        
    return mkr_lst, mkr_lst_neg, major_dict, minor_dict


def save_markers_all(mkr_lst, mkr_lst_neg, major_dict, minor_dict, file):
    
    cols = ['cell_type_major', 'cell_type_minor', 'cell_type_subset', 'exp', 'markers'] 
    df_pos = pd.DataFrame(columns = cols)
    df_neg = pd.DataFrame(columns = cols)

    mkr_lst_tmp = mkr_lst
    ct_lst = list(mkr_lst_tmp.keys())
    mkrs = []
    ct_lst_new = []
    for c in ct_lst:
        if len(mkr_lst_tmp[c]) > 0:
            mkr_lst_tmp[c].sort()
            s = mkr_lst_tmp[c][0]
            if len(mkr_lst_tmp[c]) > 1:
                for mkr in mkr_lst_tmp[c][1:]:
                    s = s + ',%s' % mkr
            mkrs.append(s)
            ct_lst_new.append(c)
    
    df_pos['cell_type_subset'] = ct_lst_new
    df_pos['markers'] = mkrs
    df_pos['exp'] = ['pos']*len(mkrs)
    
    mkr_lst_tmp = mkr_lst_neg
    ct_lst = list(mkr_lst_tmp.keys())
    mkrs = []
    ct_lst_new = []
    for c in ct_lst:
        if len(mkr_lst_tmp[c]) > 0:
            mkr_lst_tmp[c].sort()
            s = mkr_lst_tmp[c][0]
            if len(mkr_lst_tmp[c]) > 1:
                for mkr in mkr_lst_tmp[c][1:]:
                    s = s + ',%s' % mkr
            mkrs.append(s)
            ct_lst_new.append(c)
            
    
    df_neg['cell_type_subset'] = ct_lst_new
    df_neg['markers'] = mkrs
    df_neg['exp'] = ['neg']*len(mkrs)
    
    df = pd.concat([df_pos, df_neg], axis = 0)
    
    df['cell_type_major'] = df['cell_type_subset'].copy(deep = True)
    df['cell_type_minor'] = df['cell_type_subset'].copy(deep = True)
    df['cell_type_major'].replace(major_dict, inplace = True)
    df['cell_type_minor'].replace(minor_dict, inplace = True)
    df.sort_values(by = ['cell_type_major', 'cell_type_major', 'cell_type_subset'], inplace = True)

    if file is None:
        return df
    else:
        df.to_csv(file, sep = '\t', index = False)
        return df
        

'''
pos = True
neg = False
sec = True
hla_dr = True
mhc1 = True
mhc2 = True
'''    

def Tiden_check_key_genes( gene_lst, key_genes ):
    
    b = True
    for g in list(key_genes):
        if g not in gene_lst:
            b = False
            break
    return b

def Tiden_print_error(error_code = 0):    
    if error_code == 1:
        print('ERROR: One or more of key genes for T cell subtyping (CD4, CD8A, CD8B) are not in the gene list.')
    else:
        print('ERROR: X_cell_by_gene must be a DataFrame with its columns being gene names having CD4, CD8A, CD8B')
        print('   Or, gene_names must be provided with its length equal to the column size of X_cell_by_gene, containing CD4, CD8A, CD8B.')

        
def get_stat2(df_score):

    df = df_score
    maxv = list(df.max(axis = 1))
    subtype = list(df.idxmax(axis = 1))
    #tc_subtype = [trans_dict[k] for k in tc_subtype]

    maxv2 = []
    idx2 = []
    subtype_lst = list(df.columns.values)

    for i in range(df.shape[0]):
        x = np.array(df.iloc[i])
        odr = (-x).argsort()
        if len(x) > 1:
            maxv2.append(x[odr[1]])
            idx2.append(subtype_lst[odr[1]])
        else:
            maxv2.append(0)
            idx2.append(None)

    # df_res = pd.DataFrame({'CD4+ T cell subtype': tc_subtype, 'NegLogPval': neg_log_pval}, index = df.index.values)
    df_res = pd.DataFrame({'cell_type': subtype, 'cell_type(rev)': subtype, 
                           'cell_type(1st)': subtype, 'cell_type(2nd)': idx2, 
                           'Clarity': ['-']*df.shape[0], 'Score': maxv, 'Score(2nd)': maxv2}, 
                          index = df.index.values)
    df_res['dScore'] = df_res['Score'] - df_res['Score(2nd)']
    
    return df_res


def plot_GSA_score_hist(df_score, title = 'GSA', histtype = 'bar'):
    
    df_sum = get_stat2(df_score)
    score_name = 'Score'
    target_types = list(df_sum['cell_type(1st)'].unique())
    thresholds = {}
    m1m2_ratio = {}
    plot_hist = True
    
    for t in target_types:
        target = t

        b1 = (df_sum['cell_type(1st)'] == target) #& (df_sum['-logP'] > -np.log10(pval_th))
        v1 = df_sum.loc[b1, score_name]
        m1 = df_sum.loc[b1, score_name].max() #.mean()
        n1 = np.sum(b1)

        b2 = (df_sum['cell_type(2nd)'] == target) # & (df_sum['-logP(2nd)'] > -np.log10(pval_th))
        if np.sum(b2) > 0:
            v2 = df_sum.loc[b2, '%s(2nd)' % score_name]
            m2 = df_sum.loc[b2, '%s(2nd)' % score_name].max() #.mean()
            n2 = np.sum(b2)
            
        if plot_hist:
            x1 = df_sum.loc[b1, score_name]
            mnv = x1.min()
            mxv = x1.max()
            X = x1
            if np.sum(b2) > 0:
                x2 = df_sum.loc[b2, '%s(2nd)' % score_name]
                mnv = min(mnv, x2.min())
                mxv = max(mxv, x2.min())
                X = np.array([x1,x2])
                
            if mnv < mxv:
                bins = np.arange(mnv, mxv, (mxv-mnv)/50)
                plt.figure(figsize = (5,3))
                # df_sum.loc[b1, score_name].hist(bins = 50, log = True, alpha = 0.7)
                plt.hist(X, bins = bins, log = True, alpha = 0.8, histtype = histtype)
                plt.title('Histogram for %s %s score' % (t, title))
                plt.xlabel('score')
                plt.ylabel('Number of cells')
                if np.sum(b2) > 0:
                    # df_sum.loc[b2, '%s(2nd)' % score_name].hist(bins = 50, log = True, alpha = 0.5)
                    # plt.hist(x2, bins = bins, log = True, alpha = 0.5, density = True)
                    plt.legend(['Primary', 'Secondary'])
                else:
                    plt.legend(['Primary'])
            else:
                print('Histogram for %s %s not available.' % (t, title))
    return 


def plot_GSA_score_violin(df_score, title = 'GSA', split = True, 
                          scale = 'width', inner = 'quartile', log = True, width = 1):
    
    df_sum = get_stat2(df_score)
    if log:
        df_score = np.log2(df_score + 1)
        
    score_name = 'Score'
    target_types = list(df_score.columns.values) # list(df_sum['cell_type(1st)'].unique())
    thresholds = {}
    m1m2_ratio = {}
    plot_hist = True
    
    cnt = 0
    for t in target_types:
        target = t

        b1 = (df_sum['cell_type(1st)'] == target) 
        b2 = (df_sum['cell_type(1st)'] != target) 
            
        X = None
        if plot_hist:
            x1 = df_score.loc[b1, t]
            if len(x1) > 10:
                if np.sum(b2) > 10:
                    X = pd.DataFrame( {'Score': x1} )
                    X['cell type'] = 'Target'
                    X['Cell type'] = t
                    
                    x2 = df_score.loc[b2, t] # df_sum.loc[b2, '%s(2nd)' % score_name]
                    # X = np.array([x1,x2])
                    X2 = pd.DataFrame( {'Score': x2} )
                    X2['cell type'] = 'Non-target'
                    X2['Cell type'] = t
                    X = pd.concat([X,X2], axis = 0)
            
                    if X is not None:
                        if cnt == 0:
                            df = X
                        else:
                            df = pd.concat([df,X], axis = 0)
                        cnt += 1
            
    if (cnt > 0) & SEABORN:
        nn = len(list(df['Cell type'].unique()))
        plt.figure(figsize = (1.3*nn, 4))
        sns.violinplot(x="Cell type", y="Score", hue="cell type", inner = inner,
                    data=df, palette="muted", split=split, scale = scale, width = width)
        plt.xticks(rotation = 20, ha = 'center', fontsize = 12)   
        plt.yticks(fontsize = 12)
        plt.title(title, fontsize = 14)
        # plt.xlabel('Cell type', fontsize = 12)
        plt.xlabel(None)
        if log:
            plt.ylabel('Log2(1+Score)', fontsize = 12)
        else:
            plt.ylabel('Score', fontsize = 12)
        plt.show()
    return 


def plot_roc_result(df_score, y_true, method = 'gmm'):
    
    label = list(df_score.columns.values)
    
    bs = y_true == label[0] 
    for l in label[1:]:
        bs = bs | (y_true == l)
    
    plt.figure(figsize=(5, 5), dpi=100)
    clst = ['darkorange', 'red', 'gold', 'yellow', 'firebrick', 'orange', \
            'magenta', 'crimson', 'violet', 'mediumorchid', 'purple', \
            'blueviolet', 'lime', 'turquoise']*4

    all_cells = ''
    fprs = []
    tprs = []
    aucs = []
    cells= []
    for k, l in enumerate(label):
        if (k%5 == 0) & (k>0):
            all_cells = all_cells + '\n'
        all_cells = all_cells + '%s, ' % l
        
        y_conf_1 = df_score.loc[bs, l]
        
        label_tmp = copy.deepcopy(label)
        label_tmp.remove(l)
        y_conf_0 = df_score.loc[bs, label_tmp].max(axis = 1)

        y_odd = y_conf_1 - y_conf_0

        bn = ~np.isnan(y_odd)
        y_odd = y_odd[bn]
        
        y = y_true[bs&bn]

        target = l
        if (y == target).sum() > 1:
            try:
                fpr, tpr, _ = roc_curve(y.ravel(), y_odd.ravel(), pos_label = target)
                roc_auc = auc(fpr, tpr)
                fprs.append(fpr)
                tprs.append(tpr)
                aucs.append(roc_auc)
                cells.append(target)
            except:
                pass
                # print(y)
                # print(y_odd)

    odr = np.array(aucs).argsort()
    for k, o in enumerate(reversed(odr)):
        title = "%s (AUC = %0.3f)" % (label[o], aucs[o])
        plt.plot(fprs[o], tprs[o], color=clst[k], label=title)
            
    plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate", fontsize = 12) #, fontsize = 11)
    plt.ylabel("True Positive Rate", fontsize = 12) #, fontsize = 11)

    if method == 'gmm':
        Method = 'Gaussian Mixture Model'
    elif method == 'logreg':
        Method = 'Logistic Regression Model'
    else:
        Method = method

    plt.title("ROC for identifying \n[%s] \nusing %s" % (all_cells[:-2], Method)) #, fontsize = 11 )
    plt.legend(loc="lower right") #, fontsize = 10)
    # plt.legend(loc="upper left", bbox_to_anchor=(1.03, 1)) # , fontsize = 12)
    plt.show()
    
    return cells, aucs


def show_summary( df_pred, summary, pth_mult = 3, level = 3, 
                    split = False, scale = 'area', inner = 'quartile', 
                    vlog = False, vwidth = 1 ):
    
    smry_res = summary['GSA_summary']
    pval_th, pth_mult, pth_min = summary['parameters']
    
    ## Cluster filter
    y_clust = df_pred['cluster'].astype(int) - 1
    score = smry_res['Major_Type']['-logP']
    pth = pval_th # smry_res['Pval_threshold']
    #'''
    if level > 0:
        find_pct_cutoff(score, y_clust, pth = -np.log10(pth), 
                        pct_margin_multi = pth_mult, pct_min_rf = pth_min, 
                        print_report = True)
    #'''

    ## Separability report

    ## GSA ROC
    smry_score = summary['GSA_scores']
    keys = list(smry_score.keys())

    key = 'Major_Type'
    df_score = smry_score[key]
    ys = smry_res[key]['cell_type']
    if level > 1:
        cells, aucs = plot_roc_result(df_score, ys, 'GSA')    

    ## ID model ROC
    smry_score = summary['Identification_model_scores']
    keys = list(smry_score.keys())

    key = 'Major_Type'
    df_score = smry_score[key]
    ys = smry_res[key]['cell_type']
    if level > 1:
        cells, aucs = plot_roc_result(df_score, ys, 'gmm')    

    ## GSA score histogram
    if level > 2:
        smry_score = summary['GSA_scores']
        for key in smry_score.keys():
            df_score = smry_score[key]
            if isinstance(df_score, pd.DataFrame):
                plot_GSA_score_violin(df_score, title = key, split = split, 
                                      scale = scale, inner = inner, 
                                      log = vlog, width = vwidth)

        '''
        ## ID model score histogram
        smry_score = summary['Identification_model_scores']
        for key in smry_score.keys():
            df_score = smry_score[key]
            if isinstance(df_score, pd.DataFrame):
                plot_GSA_score_violin(df_score, title = 'identifier', histtype = 'bar')
        '''
    return

    
def X_normalize(X):    
    return X.div(X.sum(axis=1)*0.0001 + 0.0001, axis = 0)


def X_scale(X, max_val = 10):    
    m = X.mean(axis = 0)
    s = X.std(axis = 0)
    
    Xs = X.sub(m).mul((s > 0)/(s+ 0.0001))
    Xs.clip(upper = max_val, lower = -max_val, inplace = True)
    
    return Xs


def select_variable_genes(X, log_transformed = False):
    
    if log_transformed:
        Xs = 10**X
        Xs = Xs - Xs.min().min()
    else:
        Xs = X
        
    m = Xs.mean(axis = 0)
    s = Xs.var(axis = 0)

    lm = np.log10(m+0.000001)
    ls = np.log10(s+0.000001)

    z = np.polyfit(lm,ls, 2)
    p = np.poly1d(z)
    s_fit = 10**(p(lm))

    Xt = Xs.sub(m).mul(1/np.sqrt(s_fit))
    Xt.clip(upper = np.sqrt(Xt.shape[0]), inplace = True)

    sr = Xt.std(axis = 0)
    odr = np.array(sr).argsort()
    s_th = sr[odr[-2000]]
    genes = list(Xt.columns.values[sr >= s_th])
    
    return genes


def clustering_alg(X_pca, clust_algo = 'lv', N_clusters = 25, resolution = 1, N_neighbors = 10):
    
    if clust_algo[:2] == 'gm':
        gmm = mixture.GaussianMixture(n_components = int(N_clusters), random_state = 0)
        cluster_label = gmm.fit_predict(X_pca)
        return cluster_label, gmm
    elif clust_algo[:2] == 'km':
        km = cluster.KMeans(n_clusters = int(N_clusters), random_state = 0)
        km.fit(X_pca)
        cluster_label = km.labels_
        return cluster_label, km
    else:
        adj = kneighbors_graph(X_pca, int(N_neighbors), mode='connectivity', include_self=True)
        louvain = Louvain(resolution = resolution)
        cluster_label = louvain.fit_transform(adj)        
        return cluster_label, louvain


def GSA_cell_subtyping( X_cell_by_gene, X_pca, mkrs_pos, mkrs_neg, verbose = False ):
    
    X = X_cell_by_gene    
    genes = list(X.columns.values)

    mkr_lst_dict = copy.deepcopy(mkrs_pos)
    mkr_lst_dict_neg = copy.deepcopy(mkrs_neg)
    for key in list(mkr_lst_dict.keys()):
        mkrs = mkr_lst_dict[key]
        mkrs2 = list(set(mkrs).intersection(genes))    
        if len(mkrs2) < len(mkrs):
            mkr_lst_dict[key] = mkrs2
            s = ''
            cnt = 0
            for mkr in mkrs:
                if mkr not in mkrs2:
                    s = s + '%s,' % mkr
                    cnt += 1
            if len(s) > 1:
                s = s[:-1]
                if verbose:
                    print('   WARNING: %15s pos_mkrs in db: %3i, where %2i missing (%s)' % (key, len(mkrs), cnt, s))
            '''
            else:
                if verbose:
                    print('%15s pos_mkrs in db = %2i' % (key, len(mkrs)))
            #'''
        if key in mkr_lst_dict_neg.keys():
            mkrs3 = mkr_lst_dict_neg[key]
            if len(mkrs3) > 0:
                mkrs4 = list(set(mkrs3).intersection(genes))    
                if len(mkrs4) < len(mkrs3):
                    mkr_lst_dict_neg[key] = mkrs2
                    s = ''
                    cnt = 0
                    for mkr in mkrs3:
                        if mkr not in mkrs4:
                            s = s + '%s,' % mkr
                            cnt += 1
                    if len(s) > 1:
                        s = s[:-1]
                        if verbose:
                            print('   WARNING: %15s neg_mkrs in db: %3i, where %2i missing (%s)' % (key, len(mkrs), cnt, s))

            
    if verbose: print('GSA .. ', end = '')
    start_time = time.time()
    
    ## Get stats for CD4+ T cells only

    Xb = X > 0
    N = Xb.shape[1]
    k = Xb.sum(axis = 1)

    mkr_stat = {}
    df = pd.DataFrame(index = Xb.index.values)
    dfn = pd.DataFrame(index = Xb.index.values)

    n_mkr = {}
    for key in list(mkr_lst_dict.keys()):
        n_mkr[key] = len(mkr_lst_dict[key])
        if key in mkr_lst_dict_neg.keys():
            n_mkr[key] = n_mkr[key] + len(mkr_lst_dict_neg[key])
            
    for key in list(mkr_lst_dict.keys()):
        mkrs = mkr_lst_dict[key]
        b = Xb[mkrs]
        n = b.sum(axis = 1)
        M = len(mkrs)
        
        if key in mkr_lst_dict_neg.keys():
            mkrs_neg = mkr_lst_dict_neg[key]
            b_neg = ~Xb[mkrs_neg]
            n = n + b_neg.sum(axis = 1)
            M = M + len(mkrs_neg)

        neg_log_pval = -hypergeom.logsf(n-1, N, M, k)*np.log10(np.exp(1))
        df[key] = neg_log_pval
        dfn[key] = n

    ## KNN reg.
    '''
    neigh = KNeighborsRegressor(n_neighbors=10)
    neigh.fit(X_pca, df)
    df2 = neigh.predict(X_pca)
    df.loc[:,:] = df2
    #'''
    '''
    neigh = KNeighborsRegressor(n_neighbors=12)
    neigh.fit(X_pca, df)
    dist, neighbors = neigh.kneighbors(X_pca)

    idx = df.index.values
    dfr = df.copy(deep=True)
    for i in range(neighbors.shape[0]):
        neighbor_lst = list(neighbors[i,:])
        dfr.iloc[i] = df.loc[idx[neighbor_lst]].median(axis = 0) 
    df.loc[:,:] = dfr
    #'''

    ## Get final results
    num = list(np.arange(df.shape[1]))
    name = list(df.columns.values)
    trans_dict = dict(zip(num,name))
    
    if len(mkr_lst_dict.keys()) == 1:
        c = name[0]
        neg_log_pval = df[c]
        subtype = [c]*df.shape[0]
        
        # df_res = pd.DataFrame({'CD4+ T cell subtype': tc_subtype, 'NegLogPval': neg_log_pval}, index = df.index.values)
        df_res = pd.DataFrame({'cell_type': subtype, 'cell_type(rev)': subtype, 
                               'cell_type(1st)': subtype, 'Overlap': dfn[c]},
                               index = df.index.values)
#                                'cell_type(2nd)': idx2, 'Overlap(2nd)': dfn[c],
#                                'Clarity': ['-']*df.shape[0], '-logP': maxv, '-logP(2nd)': maxv2}, 
        df_res['cell_type(2nd)'] = None
        df_res['Overlap(2nd)'] = None
        df_res['Clarity'] = ['-']*df.shape[0]
        df_res['-logP'] = df[c]
        df_res['-logP(2nd)'] = 0
                               
        df_res['-logP-logP(2nd)'] = df[c]
    else:
        neg_log_pval = df.max(axis = 1)
        subtype = list(df.idxmax(axis = 1))
        #tc_subtype = [trans_dict[k] for k in tc_subtype]

        maxv = []
        maxv2 = []
        idx2 = []
        Nn = []
        N1 = []
        N2 = []
        subtype_lst = list(df.columns.values)

        for i in range(df.shape[0]):
            x = np.array(df.iloc[i])
            n = np.array(dfn.iloc[i])
            odr = (-x).argsort()
            maxv.append(x[odr[0]])
            maxv2.append(x[odr[1]])
            idx2.append(subtype_lst[odr[1]])

            k1 = n_mkr[subtype_lst[odr[0]]]
            k2 = n_mkr[subtype_lst[odr[1]]]
            n1 = n[odr[0]]
            n2 = n[odr[1]]
            s1 = '%i/%i' % (n1, k1)
            s2 = '%i/%i' % (n2, k2)
            N1.append(s1)
            N2.append(s2)
            Nn.append(n1)

        # df_res = pd.DataFrame({'CD4+ T cell subtype': tc_subtype, 'NegLogPval': neg_log_pval}, index = df.index.values)
        df_res = pd.DataFrame({'cell_type': subtype, 'cell_type(rev)': subtype, 
                               'cell_type(1st)': subtype, 'Overlap': N1,
                               'cell_type(2nd)': idx2, 'Overlap(2nd)': N2,
                               'Clarity': ['-']*df.shape[0], '-logP': maxv, '-logP(2nd)': maxv2}, 
                              index = df.index.values)
        df_res['-logP-logP(2nd)'] = df_res['-logP'] - df_res['-logP(2nd)']

        b = np.array(Nn) == 0
        df_res.loc[b,'cell_type'] = 'unassigned'
        df_res.loc[b,'cell_type(rev)'] = 'unassigned'
        df_res.loc[b,'cell_type(1st)'] = 'unassigned'
        df_res.loc[b,'cell_type(2nd)'] = 'unassigned'
       
    # if verbose: print('done. (%i s)' % round(time.time() - start_time))
    df_score = df
    return df_res, df_score


GMM_MIN_SCORE = -100

def C_major_gmm(X_pca, ys, class_names, method = 'gmm', N_components = 8, verbose = False ):
    
    if verbose: print('Fitting GMM .. ', end = '')
    start_time = time.time()

    bs = (ys == class_names[0])
    for cname in list(class_names[1:]):
        bs = bs | (ys == cname)
        
    X = np.array(X_pca.loc[bs,:])
    y = ys[bs]

    ## Select training cells
    df_score = pd.DataFrame(index = X_pca.index.values)
    '''
    sz = []
    for cname in list(class_names):
        bs = y == cname
        sz.append(np.sum(bs))
    szmx = np.array(sz).max()
    if verbose: print('Nc: %i, Max: %i, Med: %i, Min: %i' % (len(sz), szmx, np.median(sz), np.min(sz)))
    '''
    cnt = 0
    for cname in list(class_names):
    
        X_tmp = X[y == cname,:]
        y_tmp = y[y == cname]
        
        # N_comp = min(max(int(N_components * len(y_tmp)/szmx),1), len(y_tmp))
        # if len(y_tmp) < 50: N_comp = 1
        #'''
        sqrtN = max(int(np.sqrt(X_pca.shape[0])/2), 1)
        if len(y_tmp) <= N_components:
            N_comp = 1
        else:
            N_comp = min(N_components, sqrtN)
        #'''

        if len(y_tmp) > N_comp:
            if method == 'gmm':
                gmm = mixture.GaussianMixture(n_components = int(N_comp), random_state = 0)
            else:
                gmm = mixture.BayesianGaussianMixture(n_components = int(N_comp), 
                                                      max_iter = 1000, random_state = 0)

            gmm.fit(X_tmp)
            y_conf = gmm.score_samples(X_pca)
            df_score['%s' % cname] = y_conf
            cnt += 1
    
    if cnt > 1:
        y_pred = df_score.idxmax(axis = 1)
    else:
        y_pred = pd.Series(index = X_pca.index.values, dtype = float)
        if cnt == 1:
            y_pred[:] = df_score.columns.values[0]
        else:
            y_pred[:] = 'unassigned'
        
    
    # if verbose: print('done. (%i s)' % round(time.time() - start_time))
    return y_pred, df_score


def C_major_logreg(X_pca, ys, class_names, verbose = False ):
    
    if verbose: print('Fitting Logistic Regression model .. ', end = '')
    start_time = time.time()

    bs = (ys == class_names[0])
    for cname in list(class_names[1:]):
        bs = bs | (ys == cname)
        
    X = np.array(X_pca.loc[bs,:])
    y = ys[bs]

    NCV = LOGREG_NCV
    MAX_ITER = LOGREG_MAX_ITER
    param_grid = LOGREG_PARAM_GRID

    n_samples, n_features = X.shape
    cv = mod_sel.StratifiedKFold(n_splits=NCV)
    classifier = lm.LogisticRegression(penalty = 'elasticnet', 
                                       max_iter = MAX_ITER, solver = 'saga', 
                                       class_weight = 'balanced',tol = LOGREG_TOL)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # warnings.filterwarnings("ignore", category=ConvergenceWarning)
        gs = mod_sel.GridSearchCV(classifier, param_grid, cv=cv, scoring='balanced_accuracy', refit = True, n_jobs = NCV, verbose = 0)

        gs.fit(X,y)
        print(gs.best_params_, gs.best_score_)

    classifier = gs.best_estimator_
    y_pred = classifier.predict(X_pca)
    y_prob = classifier.predict_proba(X_pca)

    cols = list(classifier.classes_)
    df_score = pd.DataFrame(np.log10(y_prob + 1e-100), index = X_pca.index.values, columns = cols)
            
    y_pred = df_score.idxmax(axis = 1)
    
    # if verbose: print('done. (%i s)' % round(time.time() - start_time))
    return y_pred, df_score
      
    
def get_threshold(values, target_FPR = 0.05, upper = True):
    
    z = np.array(values)
    odr = z.argsort()
    n = int(round(len(odr)*target_FPR))
    if n >= len(odr):
        n = int(len(odr)-1)
    if upper:
        th = z[odr[-n]]
    else:
        th = z[odr[n]]
    return th

def get_threshold_N(values, NN = 1, upper = True):
    
    z = np.array(values)
    odr = z.argsort()
    n = int(NN)
    if n >= len(odr):
        n = int(len(odr)-1)
    if upper:
        th = z[odr[-n]]
    else:
        th = z[odr[n]]
    return th

#'''

def get_normal_pdf( x, mu, var, nbins):
    
    y = np.array(x)
    mn_x = y.min()
    mx_x = y.max()
    L = 100
    # dx = len(y)*(mx_x-mn_x)/L
    dx = (mx_x-mn_x)/nbins
    xs = np.arange(mn_x,mx_x, dx )
    pdf = (dx*len(y))*np.exp(-((xs-mu)**2)/(2*var+1e-10))/(np.sqrt(2*math.pi*var)+1e-10) + 1e-10
    return pdf, xs


def get_threshold_using_param( param, Target_FPR = 0.1 ):
    
    w0, m0, v0, w1, m1, v1 = param 
    mxv = m1+np.sqrt(v1)
    mnv = m0
    z = np.arange(mnv,mxv, (mxv-mnv)/1000)

    e0 = 1 - 0.5*(1 + erf((z - m0)/np.sqrt(v0)))
    e1 = 1 - 0.5*(1 + erf((z - m1)/np.sqrt(v1)))
    fpr = e0/e1
    
    if fpr.min() < Target_FPR:
        i = np.abs(fpr-Target_FPR).argmin()
        threshold = z[i]
    else: 
        i = fpr.argmin()
        threshold = z[i]
            
    return threshold

def bimodal_fit( x_score ):
    
    x = x_score

    gmm = mixture.GaussianMixture(n_components = 2, random_state = 0)
    y = gmm.fit_predict(np.array(x).reshape(-1, 1))

    mns = [m[0] for m in gmm.means_]
    cvs = [cv[0,0] for cv in gmm.covariances_]

    wgs = gmm.weights_           
    if mns[0] < mns[1]:
        w0, w1 = wgs[0], wgs[1]
        m0, m1 = mns[0], mns[1]
        v0, v1 = cvs[0], cvs[1]
    else:
        w0, w1 = wgs[1], wgs[0]
        m0, m1 = mns[1], mns[0]
        v0, v1 = cvs[1], cvs[0]

    # param = [w0, m0, v0, w1, m1, v1]
            
    return w0, m0, v0, w1, m1, v1
    
    
def get_threshold_from_GSA_result(df_sum, pval_th = 0.05, target_FPR = 0.05, 
                                   verbose = False, plot_hist = False):
    
    target_types = list(df_sum['cell_type(1st)'].unique())
    thresholds = {}
    m1m2_ratio = {}
    
    if verbose:
        print('Thresholds:')
        
    if len(target_types) == 1:
        
        t = target_types[0]
        x_score = df_sum['-logP'] # np.log10(df_sum['-logP'] + 1e-3)
        param = bimodal_fit( x_score )
        thresh = get_threshold_using_param( param, Target_FPR = target_FPR )
        w0, m0, v0, w1, m1, v1 = param 
        # m0 = 10**m0
        # m1 = 10**m1
        # thresh = 10**thresh
        thresholds[t] = thresh
        m1m2_ratio[t] = m1/m0 # (10**m1)/(10**m0)

        if verbose:
            print('   %16s: %5.3f, %5.3f > %5.2f, %5.2f > %i/%i' % \
              (t, m1, m0, thresh, m1/m0, 
               np.sum(x_score >= thresh), len(x_score)))
    else:
        for t in target_types:
            if t != 'unassigned':
                target = t

                b1 = (df_sum['cell_type(1st)'] == target) & \
                     (df_sum['-logP'] >= -np.log10(pval_th))
                v1 = df_sum.loc[b1, '-logP']
                m1 = df_sum.loc[b1, '-logP'].mean()
                n1 = np.sum(b1)

                b2 = (df_sum['cell_type(2nd)'] == target) & \
                     (df_sum['-logP(2nd)'] > -np.log10(pval_th))
                if np.sum(b2) > 0:
                    v2 = df_sum.loc[b2, '-logP(2nd)']
                    m2 = df_sum.loc[b2, '-logP(2nd)'].mean()
                    n2 = np.sum(b2)

                    NN = int(n1*target_FPR)
                    thresh = get_threshold_N(v2, NN = NN, upper = True)
                    '''
                    th_lst = list((500 - np.arange(500))*v1.max()/500)
                    thresh = -np.log10(pval_th)
                    for th in th_lst:
                        z1 = np.sum(v1 > th)/len(v1)
                        z2 = np.sum(v2 > th)/len(v2)
                        if z2/(z1+z2) > target_FPR:
                            thresh = th
                    '''
                    thresh = max(thresh, -np.log10(pval_th))
                    # th_min = get_threshold(v1, target_FPR = 0.1, upper = False)
                    # if thresh < th_min: 
                    #     thresh = th_min
                    if thresh > m1: thresh = m1
                    elif (thresh < m2) & (m1 > m2): thresh = m2
                    
                else:
                    if len(v1) > 1:
                        thresh = get_threshold(v1, target_FPR = 0.1, upper = False)
                        if thresh > m1: thresh = m1
                        thresh = max(thresh, -np.log10(pval_th))
                        m2 = 0.01
                    else:
                        thresh = 0
                        m2 = 0.01

                thresholds[t] = thresh
                m1m2_ratio[t] = m1/m2

                if verbose:
                    print('   %16s: %5.3f, %5.3f > %5.2f, %5.2f > %i/%i' % \
                      (t, m1, m2, thresh, m1/m2, np.sum(v1 >= thresh), len(v1) ) )
            
    return thresholds, m1m2_ratio
#'''

def get_stat(df_score):

    df = df_score
    maxv = list(df.max(axis = 1))
    subtype = list(df.idxmax(axis = 1))
    #tc_subtype = [trans_dict[k] for k in tc_subtype]

    maxv2 = []
    idx2 = []
    subtype_lst = list(df.columns.values)

    for i in range(df.shape[0]):
        x = np.array(df.iloc[i])
        odr = (-x).argsort()
        if len(x) > 1:
            maxv2.append(x[odr[1]])
            idx2.append(subtype_lst[odr[1]])
        else:
            maxv2.append(0)
            idx2.append(None)

    # df_res = pd.DataFrame({'CD4+ T cell subtype': tc_subtype, 'NegLogPval': neg_log_pval}, index = df.index.values)
    df_res = pd.DataFrame({'cell_type': subtype, 'cell_type(rev)': subtype, 
                           'cell_type(1st)': subtype, 'cell_type(2nd)': idx2, 
                           'Clarity': ['-']*df.shape[0], 'Score': maxv, 'Score(2nd)': maxv2}, 
                          index = df.index.values)
    df_res['dScore'] = df_res['Score'] - df_res['Score(2nd)']
    
    return df_res


def get_threshold_from_GMM_result(df_sum, target_FPR = 0.05, 
                                   verbose = False, plot_hist = False):
    
    target_types = list(df_sum['cell_type(1st)'].unique())
    thresholds = {}
    m1m2_ratio = {}
    
    if verbose:
        print('Thresholds:')
        
    for t in target_types:
        target = t

        b1 = (df_sum['cell_type(rev)'] == target) #& (df_sum['-logP'] > -np.log10(pval_th))
        v1 = df_sum.loc[b1, 'Score']
        m1 = df_sum.loc[b1, 'Score'].max() #.mean()
        n1 = np.sum(b1)

        ##############################
        b2 = (df_sum['cell_type(2nd)'] == target) & \
             (df_sum['cell_type(rev)'] != 'unassigned') 
        if np.sum(b2) > 0:
            v2 = df_sum.loc[b2, 'Score(2nd)']
            m2 = df_sum.loc[b2, 'Score(2nd)'].max() #.mean()
            n2 = np.sum(b2)

            '''
            nx = int(n1*target_FPR)
            odr = np.array(v2).argsort()
            if nx == 0: nx = 1
            elif nx >= len(v2): nx = 0
            thresh = v2[odr[-nx]] 
            #'''
            NN = int(n1*target_FPR)
            thresh = get_threshold_N(v2, NN = NN, upper = True)
            #thresh = np.max(v2)
            # get_threshold(v2, target_FPR = target_FPR, upper = True)
        else:
            m2 = df_sum.loc[b1, 'Score'].mean()
            v2 = df_sum.loc[b1, 'Score']
            thresh = m2 #get_threshold(v1, target_FPR = target_FPR, upper = False)
            print('WARNING: thresholds for %s was set to %6.2f - 10 = %6.2f' % (t, m1, thresh))

        thresholds[t] = thresh
        m1m2_ratio[t] = m1-m2

        if verbose:
            print('   %16s: %6.2f, %6.2f > Threshold: %6.2f > %i/%i,  %i/%i' % \
              (t, m1, m2, thresh, np.sum(v1 >= thresh), len(v1), np.sum(v2 >= thresh), len(v2)))
            
    return thresholds, m1m2_ratio


def find_idx(bvec):
    bv = list(bvec)
    wh = []
    for k, b in enumerate(bv):
        if b: wh.append(k)
    return wh

def find_means(X, y_clust):
    
    clust_lst = list(set(y_clust))
    clust_lst.sort()
    n_clust = len(clust_lst)
    mean_mat = np.zeros([n_clust, X.shape[1]])
    csize_vec = np.zeros(n_clust)
    # df = pd.DataFrame(index = clust_lst)
    
    for i in clust_lst:
        b = y_clust == i
        mn = X.loc[b,:].mean(axis = 0)
        mean_mat[i,:] = mn
        csize_vec[i] = np.sum(b)
    
    csize_vec = csize_vec/np.sum(csize_vec)
    return mean_mat, csize_vec


def find_neighbor_list(mns, c_mult = 1):
    
    dist_mat = np.ones([mns.shape[0], mns.shape[0]])*1e6
    dist_vec = np.zeros([mns.shape[0]])
    for i in range(mns.shape[0]):
        for j in range(mns.shape[0]):
            if i > j:
                d = np.sum((mns[i,:] - mns[j,:])**2)
                dist_mat[i,j] = d
                dist_mat[j,i] = d

    #'''
    mxd = dist_mat.max().max()
    if mxd > 1e6:
        for i in range(mns.shape[0]):
            dist_mat[i,i] = mxd

    dist_vec = dist_mat.min(axis = 1)
    # dist_vec.sort()

    dth = np.median(dist_vec)*c_mult
    neighbor_mat = (dist_mat < dth)
    #'''
    '''
    for i in range(mns.shape[0]):
        for j in range(mns.shape[0]):
            if i != j:
                if (dist_mat[i,j] < std[i]*N_mult) | (dist_mat[i,j] < std[j]*N_mult):
                    neighbor_mat[i,j] = True
    '''

    n_lst_dict = {}
    for i in range(mns.shape[0]):
        # b = neighbor_mat[i,:]
        # n_lst_dict[i] = find_idx(b)
        odr = dist_mat[i,:].argsort()
        n_lst_dict[i] = []
        for o in odr:
            if dist_mat[i,o] > dth:
                break
            else:
                n_lst_dict[i].append(o)
        
    return n_lst_dict, dist_vec


def cluster_basis_correction(X_pca, y, cobj, y_clust, pmaj = 0.7, 
                             verbose = False, 
                             Nnei = CLUSTER_BASIS_CORECTION_N_NEIGHBORS):
    
    ys = pd.Series(y).copy(deep = True)
    
    cs = CLUSTER_BASIS_CORECTION_DIST_SCALE
    mns, svec = find_means(X_pca, y_clust)
    
    for loop, cmult in enumerate(cs):
        n_lst_dict, min_dist_vec = find_neighbor_list(mns, cmult)
        
        # print('Loop = %i' % loop)
        cluster_lst = list(set(y_clust))
        for clst in cluster_lst:
            b = y_clust == clst
            #'''
            neis = ''
            if len(n_lst_dict[clst]) > 0:
                for c in n_lst_dict[clst][:int(Nnei)]:
                    neis = neis + '%i,' % c
                    # b = b | (y_clust == c)
            #'''

            if np.sum(b) > 1:
                ctbl = pd.Series(y[b]).value_counts()
                idx = ctbl.index.values
                N = ctbl.sum()

                if N < X_pca.shape[0]*CLUSTER_BASIS_CORECTION_MIN_PCT_TO_INVALIDATE:
                    ys[b] = 'unassigned'
                    if verbose:
                        print('   Small:  %s(%s,%i) -> unassigned' % 
                              (clst, idx[0], N))
                else: 
                    majs = []
                    sz = []
                    T = ''
                    NT = ''
                    bn = np.full(len(b), False)
                    nx = min(len(n_lst_dict[clst]), Nnei)
                    for c in n_lst_dict[clst][:int(Nnei)]:
                        bt = (y_clust == c)
                        bn = y_clust == clst # bn | bt
                        ctbl2 = pd.Series(y[bt]).value_counts()
                        idx2 = ctbl2.index.values
                        majs.append(idx2[0])
                        sz.append(np.sum(bt))
                        T = T + '%s,' % idx2[0]
                        NT = NT + '%i,' % np.sum(bt)

                    if len(majs) == 0:
                        ys[b] = idx[0]
                        if verbose: print('   No nei: %s(%s,%i/%i), %s(%s/%s) -> %s' % 
                              (clst, idx[0], ctbl[idx[0]], N, neis[:-1], T[:-1], NT[:-1], idx[0]))

                    elif len(set(majs)) == 1:
                        Nxx = np.sum(sz)
                        if majs[0] == idx[0]:
                            ys[b] = idx[0]                                            
                            # if verbose: print('  Asis0: %s(%s,%i/%i), %s(%s/%s) -> %s' % 
                            #       (clst, idx[0], ctbl[idx[0]], N, neis[:-1], T[:-1], NT[:-1], majs[0]))
                        else:
                            Nxx = np.min(sz)
                            Nss = np.sum(sz)
                            if (Nxx > N) | (Nss > (len(sz)*N)):
                                ys[b] = majs[0]
                                if verbose: print('   Change: %s(%s,%i/%i), %s(%s/%s) -> %s' % 
                                      (clst, idx[0], ctbl[idx[0]], N, neis[:-1], T[:-1], NT[:-1], majs[0]))
                            else:
                                if verbose: print('   Skip:   %s(%s,%i/%i), %s(%s/%s) ' % 
                                      (clst, idx[0], ctbl[idx[0]], N, neis[:-1], T[:-1], NT[:-1]))

                    elif idx[0] in majs:
                        bxx = np.array(majs) == idx[0]
                        # if np.sum(bxx) > 1:
                        Nxx = np.max(np.array(sz)[bxx])
                        if Nxx > N:
                            ys[b] = idx[0]                                            
                            # if verbose: print('  Asis1a: %s(%s,%i/%i), %s(%s/%s) -> %s' % 
                            #       (clst, idx[0], ctbl[idx[0]], N, neis[:-1], T[:-1], NT[:-1], majs[0]))
                        elif (ctbl[idx[0]] > N*pmaj):
                            ys[b] = idx[0]
                            # if verbose: print('  Asis1: %s(%s,%i/%i), %s(%s/%s) -> %s' % 
                            #       (clst, idx[0], ctbl[idx[0]], N, neis[:-1], T[:-1], NT[:-1], majs[0]))
                        else:
                            if verbose: print('   Skip:   %s(%s,%i/%i), %s(%s/%s) ' % 
                                  (clst, idx[0], ctbl[idx[0]], N, neis, T, NT))

                    else: 
                        if (ctbl[idx[0]] >= N*pmaj):
                            ys[b] = idx[0]
                            # if verbose: print('  Asis2: %s(%s,%i/%i), %s(%s/%s) -> %s' % 
                            #       (clst, idx[0], ctbl[idx[0]], N, neis[:-1], T[:-1], NT[:-1], majs[0]))
                        elif len(idx) > 1:
                            if (idx[1] == 'unassigned') & \
                               ((ctbl[idx[0]]+ctbl[idx[1]]) >= N*pmaj):
                                ys[b] = idx[0]
                                # if verbose: print('  Asis3: %s(%s,%i/%i), %s(%s/%s) -> %s' % 
                                #       (clst, idx[0], ctbl[idx[0]], N, neis[:-1], T[:-1], NT[:-1], majs[0]))
                            else:
                                if verbose: print('   Skip:   %s(%s,%i/%i), %s(%s/%s) ' % 
                                      (clst, idx[0], ctbl[idx[0]], N, neis[:-1], T[:-1], NT[:-1]))

        y = ys.copy(deep = True)
        
    return ys


def find_pct_cutoff(score, y_clust, pth = 1.3, pct_margin_multi = 3, 
                    pct_min_rf = 1, print_report = False):
    
    clust_lst = list(set(y_clust))
    # clust_lst.sort()

    nn = []
    for c in clust_lst:
        b = y_clust == c
        nn.append((score[b] >= pth).sum()/np.sum(b))
    nn.sort()

    y = np.array(nn)
    x = np.arange(len(nn))
    
    L = np.sum(y < PCT_THRESHOLD_MAX)
    if len(y) - L < 10:
        L = max(len(y) - 10, 0)
        
    if L >= 0:
        z = np.polyfit(x[L:],y[L:], 1, w = y[L:])
        p = np.poly1d(z)
        s = p(x)
        margin = np.abs(y[L:]*(y[L:] - s[L:])).max()*pct_margin_multi

        for i in reversed(range(len(nn))):
            # if np.abs(s[i] - y[i]) > margin:
            if (s[i] - y[i]) > margin:
                break

        pct_thresh = max(s[i] - margin, pct_min_rf)
        pct_thresh = min(pct_thresh, PCT_THRESHOLD_MAX)
        
        if print_report:
            plt.figure(figsize = (5,5))
            plt.plot(x, y)
            plt.plot(x, s)
            a = pct_thresh
            plt.plot([0,len(nn)-1], [a,a])
            a = pct_thresh + margin
            plt.plot([0,len(nn)-1], [a,a], '--')
            plt.xlabel('Order', fontsize = 12)
            plt.ylabel('P[Pval>=Th]', fontsize = 12)
            plt.title('Cluster filter stats.', fontsize = 14)
            # plt.grid()
            pass            
    else:
        pct_thresh = PCT_THRESHOLD_MAX
        margin = 0
    
    return pct_thresh, margin


def run_gsa_and_clf(X_pca, Xs, cobj, y_clust, mkr_lst, mkr_lst_neg, method = 'gmm', 
                     N_components = 8, pval_th = 0.05, pct_margin_multi =3, pct_min_rf = 1,
                     Target_FPR = 0.05, pmaj = 0, minor_id_sep = True, min_logP_diff = 3,
                     thresholding = False, pct_cutoff = False, print_report = False):
    
    df_res, df_GSA_score = GSA_cell_subtyping( Xs, X_pca, mkr_lst, mkr_lst_neg, 
                                               verbose = print_report )

    #######################################################
    ## Compute GSA thresholds Even if thresholding is False
    
    log_pv_th = -np.log10(pval_th)

    if pct_cutoff:
        score = df_res['-logP'].copy(deep=True)
        pct_thresh, margin = find_pct_cutoff(score, y_clust, log_pv_th,
                                             pct_margin_multi = pct_margin_multi, 
                                             pct_min_rf = pct_min_rf, print_report = False)
        # pct_thresh = min(max(pct_thresh, 1-pmaj), PCT_THRESHOLD_MAX)
    else:
        pct_thresh = 0

    ## Find the GLOBAL threshold for each cell type
    th_dict, m1m2_ratio = get_threshold_from_GSA_result(df_res, 
                                            pval_th = pval_th, 
                                            target_FPR = Target_FPR, 
                                            verbose = print_report)
        
    pct_thresh2 = min(pct_min_rf, pct_thresh)
    clst_to_exclude = []
        
    ## Apply GSA thresholds If thresholding is True
    if thresholding:
        
        cluster_lst = list(set(y_clust))
        cluster_lst.sort()

        if pct_cutoff & print_report: 
            print('P[Pv<Th] cutoff = %4.2f(%4.2f)' % (pct_thresh, margin), end = '')
        
        cnt = 0
        for clst in cluster_lst:
        ## For each cluster,
            ## select cells in the cluster
            b = y_clust == clst
            pct = (df_res.loc[b, '-logP'] >= log_pv_th).sum()/np.sum(b)
            if pct < pct_thresh:
                df_res.loc[b, 'cell_type(rev)'] = 'unassigned'
                df_res.loc[b, 'Clarity'] = 'Unclear'
                
                clst_to_exclude.append(clst)
                cnt += 1
                if print_report:
                    print(', %s(%4.2f)' % (str(clst), pct), end = '')
                    # print('   Cluster %2s : %4.2f' % (str(clst), pct))
            else:
                ## Find majority cell
                cnt_tbl = df_res.loc[b, 'cell_type(1st)'].value_counts()
                # idx = cnt_tbl.index.values
                majority = cnt_tbl.index.values[0]
                if cnt_tbl[majority] >= cnt_tbl.sum()*pmaj:
                    df_res.loc[b, 'cell_type(rev)'] = majority
                    df_res.loc[b, 'Clarity'] = '-'
                else:
                    if minor_id_sep:
                        df_res.loc[b, 'cell_type(rev)'] = majority
                        df_res.loc[b, 'Clarity'] = '-'
                    else:
                        pass
                
        if pct_cutoff & print_report: 
            print(' --> %i cluster(s) among %i excluded. ' % (cnt, len(cluster_lst)))
                
            '''
            ## Find majority cell
            cnt_tbl = df_res.loc[b, 'cell_type(1st)'].value_counts()
            idx = cnt_tbl.index.values
            majority = idx[0]
            
            if m1m2_ratio[majority] < 1:
                df_res.loc[b, 'cell_type(rev)'] = 'unassigned'
                df_res.loc[b, 'Clarity'] = 'Unclear'
            else:
                # bt1 = df_res['-logP'] < max(thresh, -np.log10(pval_th))
                bt = df_res['-logP'] < max(th_dict[majority], -np.log10(pval_th))
                df_res.loc[b&bt, 'cell_type(rev)'] = 'unassigned'
                df_res.loc[b&bt, 'Clarity'] = 'Unclear'                
            #'''
    # else:
    celltypes = list(th_dict.keys())
    df_res['>=th'] = True
    for ct in th_dict.keys():
        bt = df_res['cell_type(rev)'] == ct
        bh = df_res['-logP'] < th_dict[ct]
        if np.sum(bt&bh) >= 10:
            df_res.loc[bt&bh, '>=th'] = False
        elif np.sum(bt) < 10:
            df_res.loc[bt, '>=th'] = False
    pass
    
    ys = df_res['cell_type(rev)'].copy(deep=True)
    
    ## Only for major type identification
    if False: # (pmaj > 0):
        ys = cluster_basis_correction(X_pca, ys, 
                       cobj, y_clust, pmaj = pmaj, 
                       verbose = True)
        
    df_res['cell_type(rev2)'] = ys
    
    class_names = list(df_GSA_score.columns.values)     
    if method == 'logreg':
        y_pred, df_score = C_major_logreg(X_pca, ys, class_names, 
                                          verbose = print_report )
    elif method == 'gmm':
        y_pred, df_score = C_major_gmm(X_pca, ys, class_names, 
                                   method = method, 
                                   N_components = N_components, 
                                   verbose = print_report )
            
        #'''
        ## reject if GMM score > 0
        # MinV = df_score.min().min()
        for i in range(df_score.shape[0]):
            score = df_score.iloc[i]
            b = score > 0
            score[b] = 0 # MinV
            df_score.iloc[i] = score
        #'''
    else: 
        y_pred = ys
        df_score = None
                    
    for clst in clst_to_exclude:
        b = y_clust == clst
        y_pred[b] = 'unassigned'
            
    # df_score = df_score + df_GSA_score.clip(upper = 5)
    
    df_res['cell_type(rev3)'] = y_pred
    
    if print_report:
        b_cur = y_pred == 'unassigned'
        print('Num of unassigned cells: %i among %i' % (np.sum(b_cur), len(b_cur)))
        
    if method is not None: 
        if (thresholding): # (df_score.shape[1] > 1)
            if (df_score.shape[1] == 1):
                y_pred = ys
            else:
                df_summary = get_stat(df_score)                
                th_dict, diff_dict = get_threshold_from_GMM_result(df_summary, 
                                                    target_FPR = Target_FPR, 
                                                    verbose = print_report, 
                                                    plot_hist = False ) # print_report) 

                label_lst = list(y_pred.unique())
                for label in label_lst:
                    if (label != 'unassigned') & (label in th_dict.keys()):
                        b1 = y_pred == label
                        b2 = df_summary['Score'] <  th_dict[label]
                        b3 = df_res['-logP-logP(2nd)'] < min_logP_diff
                        y_pred[b1&b2&b3] = 'unassigned'

                bx = ys == 'unassigned' ##
                bx = y_pred == 'unassigned' # not working
                for clst in cluster_lst:
                    b = y_clust == clst
                    ## If the majority of a cluster is unassigned,
                    ## set all the cells in the cluster unassigned
                    # if False: # np.sum(b&bx) > np.sum(b)*pmaj:
                    #     y_pred[b] = 'unassigned'
                    if np.sum(b&(~bx)) > 0:
                        cnt_tbl = y_pred[b&(~bx)].value_counts()
                        # idx = cnt_tbl.index.values
                        majority = cnt_tbl.index.values[0]
                        if cnt_tbl[majority] >= cnt_tbl.sum()*pmaj:
                            y_pred[b] = majority
                        else:
                            pass
                        
    
    df_res['cell_type(rev4)'] = y_pred
    
    if (pmaj > 0):
        y_pred = cluster_basis_correction(X_pca, y_pred, 
                       cobj, y_clust, pmaj = pmaj, 
                       verbose = False)
        
    df_res['cell_type'] = y_pred
    return df_res, df_score, df_GSA_score


def check_if_separable_pairwise(df_score, y_true):
    
    label = list(df_score.columns.values)        
    aucs = pd.DataFrame( np.ones([len(label),len(label)])*SEPARABILITY_AUC_INIT_VALUE, 
                         index = label, columns = label )
    
    for k, lr in enumerate(label):
        for m, lc in enumerate(label):
            if lr != lc:
                br = y_true == lr
                bc = y_true == lc
                if (np.sum(br) < SEPARABILITY_MIN_NUM_CELLS) | (np.sum(bc) < SEPARABILITY_MIN_NUM_CELLS):
                    aucs.loc[lr,lc] = 1
                    pass
                else:
                    bs = br | bc
                    y_conf_1 = df_score.loc[bs, lr]
                    y_conf_0 = df_score.loc[bs, lc]
                    y_odd = y_conf_1 - y_conf_0
                    y = y_true[bs]

                    bn = (~np.isnan(y_odd))
                    y_odd = y_odd[bn]
                    y = y[bn]

                    target = lr
                    try:
                        fpr, tpr, _ = roc_curve(y.ravel(), y_odd.ravel(), pos_label = target)
                        roc_auc = auc(fpr, tpr)
                        aucs.loc[lr,lc] = roc_auc
                    except:
                        print('WARNING: cannot determine the separability for %s' % target)        
    return aucs

def separability_check_pairwise(df_score, ys, mkr_lst, 
                                dict_celltype_comb, verbose = False):
    
    aucs = check_if_separable_pairwise(df_score, ys)
    idxc_vec = aucs.idxmin(axis = 1)
    mina_vec = aucs.min(axis = 1)
    idxr = mina_vec.idxmin()
    idxc = idxc_vec[idxr]
    
    min_auc = aucs.loc[idxr, idxc]
    b_ok = True
    if min_auc < SEPARABILITY_THRESHOLD:
        b_ok = False
    
    if b_ok:
        if verbose: 
            print('Separability check passed: Min.AUC = %6.4f btn %s and %s' % (min_auc, idxr, idxc))
    else:
        if verbose: 
            print('Separability check failed: Min.AUC = %6.4f btn %s and %s' % (min_auc, idxr, idxc))
        names = '%s and %s' % (idxr, idxc)
        new_name = '%s_or_%s' % (idxr, idxc)
        genes = list(set(mkr_lst[idxr] + mkr_lst[idxc]))
        to_comb_name = [idxr, idxc]
        
        dict_celltype_comb[new_name] = to_comb_name
        mkr_lst[new_name] = genes
        
        del mkr_lst[idxr]                        
        del mkr_lst[idxc]                        
        if idxr in dict_celltype_comb.keys():
            del dict_celltype_comb[idxr]
        if idxc in dict_celltype_comb.keys():
            del dict_celltype_comb[idxc]
                
        if verbose: 
            print('Separability check WARNING: %s are not clearly separable.' % names)
            # print('%s are not clearly separable -> combined into one major type.' % names)
    
    return b_ok, mkr_lst, dict_celltype_comb


def rem_minor( mkr_dict, to_rem_lst ):

    cts = list(mkr_dict.keys())
    if len(to_rem_lst) > 0:
        for c in cts:
            if c in to_rem_lst:
                del mkr_dict[c]
                
    return mkr_dict


def get_maj(ivec, cto, p_cells_dict, p_min = 0.1):

    items = list(set(ivec))
    if len(items) == 1:
        return cto
    
    Num = np.zeros(len(items))
    Score = np.zeros(len(items))
    for k, item in enumerate(items):
        b = ivec == item
        Num[k] = np.sum(b)
        '''
        if item in list(p_cells_dict.keys()):
            Score[k] = np.sum(b)*p_cells_dict[item]
        else: 
            Score[k] = 0
        #'''
    k = np.argmax(Num)

    b = False
    if items[k] == STR_UNASSIGNED:
        odr = (-Num).argsort()        
        if len(odr) > 1:
            if Num[odr[1]] >= round(np.sum(Num)*(p_min)):
                k = odr[1]
            # elif Num[k] <= round(np.sum(Num)*(1-p_min)):
            #     return cto
       
    return  items[k]


def apply_knn(nei_indices, cell_type, score, pmaj = 0.5, nlog_pval_th = 2):

    cell_type_new = cell_type.copy(deep = True)
    idx_ary = cell_type.index.values
    for k in range( len(cell_type) ):
        
        cts = cell_type[ idx_ary[ list(nei_indices.iloc[k]) ] ]
        if len(set(list(cts))) <= 1:
            pass
        else:
            cnt_tbl = cts.value_counts()
            majority = cnt_tbl.index.values[0]

            if majority == 'unassigned':
                cell_type_new[idx_ary[k]] = cnt_tbl.index.values[1]
            else: 
                if (cell_type[k] == 'unassigned') | (score[k] < nlog_pval_th):
                    Pmaj = 0
                else:
                    Pmaj = pmaj

                if cnt_tbl[majority] > cnt_tbl.sum()*Pmaj:
                    cell_type_new[idx_ary[k]] = majority
            
    return cell_type_new
    

def apply_knn2(nei_indices, cell_type, score, pmaj = 0.5, nlog_pval_th = 2):

    cell_type_new = cell_type.copy(deep = True)
    idx_ary = cell_type.index.values
    for k in range( len(cell_type) ):
        if (cell_type[k] == 'unassigned') | (score[k] < nlog_pval_th):            
            cts = cell_type[ idx_ary[ list(nei_indices.iloc[k]) ] ]
            cnt_tbl = cts.value_counts()
            majority = cnt_tbl.index.values[0]
            if majority == 'unassigned':
                if len(cnt_tbl.index.values) > 1:
                    cell_type_new[idx_ary[k]] = cnt_tbl.index.values[1]
            else: 
                cell_type_new[idx_ary[k]] = majority
            
    return cell_type_new
    

def MarkerCount_GSA( X_cell_by_gene, marker_file, log_transformed = True,
                   Clustering_algo = 'lv', Clustering_resolution = 1, 
                   Clustering_base = 'pca', N_pca_components = 10, 
                   model = 'gmm', N_gmm_components = 10, 
                   ident_level = 3, minor_id_sep = True, mkr_selector = PNSH12,
                   Target_FPR = 0.05, pval_th = 0.05, pmaj = 0.8, 
                   pth_margin_mult = 3, pth_min = 0.25, min_logP_diff = 3, 
                   rem_cmn_mkr = False, target_cell_types = [],
                   minor_types_to_exclude = [], N_neighbors_minor = 31, N_neighbors_subset = 11, 
                   pval_th2 = 0.01, cycling_cell = False, 
                   ref_correct = True, verbose = True ):
    
    gene_names = None
    method = model
    N_components = N_gmm_components
    start_time = time.time()
    dict_summary_res = {}
    dict_summary_score = {}
    dict_summary_score2 = {}
    dict_summary_score3 = {}
    
    if Target_FPR >= 1:
        thresholding = False
    else:
        thresholding = True
        
    if verbose: print('MarkerCount_GSA running ..')
    else: print('MarkerCount_GSA running ..', end = '', flush = True)
        
    if isinstance(X_cell_by_gene, pd.DataFrame):
        Xs = X_cell_by_gene.copy(deep=True)
        genes_lst_org = list(X_cell_by_gene.columns.values)

        genes_lst = []
        for g in genes_lst_org:
            genes_lst.append(g.upper())
        rend = dict(zip(genes_lst_org, genes_lst))
        Xs.rename(columns = rend, inplace = True)
            
    else:
        Tiden_print_error()
        return -1

    ### Correct Major type
    if verbose: print('PCA .. ', end = '', flush = True)
    N_components_pca = N_pca_components
    pca = PCA(n_components = int(N_components_pca), copy = True, random_state = 0)
    if Xs.shape[1] < 3000:
        if not log_transformed:
            Xx = np.log2(1 + X_normalize(Xs))
        else:
            Xx = Xs
        # Xx = (X_scale(Xx, max_val = 10))
        # b = Xx.sum(axis = 0) == 0
        # Xx = Xx.loc[:,b]
        res = pca.fit(Xx)
    else:
        gene_sel = select_variable_genes(Xs, log_transformed = log_transformed)
        if not log_transformed:
            Xx = np.log2(1 + X_normalize(Xs.loc[:,gene_sel]))
        else:
            Xx = Xs.loc[:,gene_sel]
        # Xx = (X_scale(Xx, max_val = 10))
        # b = Xx.sum(axis = 0) == 0
        # Xx = Xx.loc[:,b]
        res = pca.fit(Xx)
        
    ## Get 2D UMAP projection of the selected GEP
    X_pca = Xx.dot(res.components_.T) 

    N_clusters = int(25*np.sqrt(Clustering_resolution))
    MaxNcomp = N_clusters
    sqrtN = int(np.sqrt(X_pca.shape[0])/2)
    N_comp = min(MaxNcomp, sqrtN)
    '''
    gmm = mixture.GaussianMixture(n_components = N_comp, covariance_type = 'spherical', random_state = 0)
    y_clust = gmm.fit_predict(X_pca)
    '''
    if verbose: print('Clustering .. ', end = '', flush = True)
    
    if Clustering_base == 'pca':
        X_umap = X_pca
    else:
        umap_obj = umap.UMAP(random_state=2)
        X_umap =  umap_obj.fit_transform(X_pca)
        X_umap = pd.DataFrame(X_umap, index = X_pca.index.values)
        # X_umap =  = TSNE(n_components = 3, earning_rate=300, init='pca').fit_transform(X_pca.T)
    
    y_clust, cobj = clustering_alg(X_umap, clust_algo = Clustering_algo, N_clusters = N_comp, 
                             resolution = Clustering_resolution)
    
    rend = {}
    for i in list(X_pca.columns.values):
        rend[i] = str(i)
    X_pca.rename(columns = rend, inplace = True)

    # dict_summary_res['Pval_threshold'] = pval_th
    # dict_summary_res['CObj'] = cobj    
    dict_summary_res['Xpca'] = X_pca
    
    if verbose: print('Nc = %i, ' % len(set(y_clust)), end = '', flush = True)
       
    ##### Major type identification for target cell selection
    dict_celltype_comb = {}

    if True:
        target_cell = target_cell_types
        mkr_lst, mkr_lst_neg = get_markers_major_type(marker_file, target_cell, 
                                        pnsh12 = mkr_selector, rem_common = rem_cmn_mkr, verbose = verbose)
        if len(target_cell_types) == 0:
            target_cell_types = list(mkr_lst.keys())
        else:
            sm = ''
            for key in list(set(target_cell_types).difference(list(mkr_lst.keys()))):
                sm = sm + '%s,' % key
            if verbose & (len(sm) > 0): 
                print('INFO: The marker db does not contain %s' % sm[:-1])
            target_cell_types = list(set(target_cell_types).intersection(mkr_lst.keys()))
            
        mkrs_all = []
        for key in mkr_lst.keys():
            mkrs_all = mkrs_all + mkr_lst[key]
        mkrs_all = list(set(mkrs_all))
        mkrs_exist = list(set(mkrs_all).intersection(genes_lst))
        
        pct_r = len(mkrs_exist)/len(mkrs_all)
        if pct_r == 0:
            print('ERROR: No marker genes found in the data.')
            return None, None
        
        pval_th = min( pval_th/pct_r, 0.1 )
        pct_th_min = pth_min*pct_r
        
        if verbose & (pct_r < 0.9): 
            # print('PCT reduction factor = %4.2f' % pct_r)
            print('WARNING: Too many markers not present. pv_th, min_pth -> %5.3f, %5.3f' % \
                  (pval_th, pct_r*PCT_THRESHOLD_MIN))
            # print('INFO: Too many markers not present. min_pth -> %5.3f' % (pct_r*PCT_THRESHOLD_MIN))

        for loop in range(len(target_cell_types)):
        # loop = 0
        # if True:
            df_res_major, df_score_major, df_GSA_score_major = \
                run_gsa_and_clf(X_pca, Xs, cobj, y_clust, mkr_lst, mkr_lst_neg, 
                                 method = 'gmm', N_components = N_components, pval_th = pval_th, 
                                 pct_margin_multi = pth_margin_mult, pct_min_rf = pct_th_min,
                                 Target_FPR = Target_FPR, pmaj = pmaj, minor_id_sep = minor_id_sep,
                                 min_logP_diff = min_logP_diff, thresholding = thresholding,
                                 pct_cutoff = True, print_report = verbose)

            y_pred_major = df_res_major['cell_type']
            ys_major = df_res_major['cell_type(rev)']

            if (len(target_cell_types) > 1):
                b_ok, mkr_lst, dict_celltype_comb = \
                    separability_check_pairwise( df_score_major, ys_major, 
                                        mkr_lst, dict_celltype_comb, verbose)
                if b_ok: break
            else:
                break
                
        y_pred_correct = y_pred_major
        b_cur = y_pred_correct == 'unassigned'
        if verbose:
            print('Num of unassigned cells: %i among %i' % (np.sum(b_cur), len(b_cur)))

        if verbose:
            print('Major cell type identification done. (%i)' % \
                   round(time.time() - start_time))
            # plot_roc_result(df_score_major, ys_major, method)

        df_pred = pd.DataFrame({'cell_type_major': y_pred_correct},
                                index = df_res_major.index.values)
        df_pred['cell_type_minor'] = 'unassigned'
        df_pred['cell_type_subset'] = 'unassigned'
        df_pred['cluster'] = y_clust + 1
        df_pred['cluster'] = df_pred['cluster'].astype(str)

        df_pred['cell_type(1st)'] = df_res_major['cell_type(1st)']
        df_pred['Confidence(1st)'] = df_res_major['-logP']
        df_pred['cell_type(2nd)'] = df_res_major['cell_type(2nd)']
        df_pred['Confidence(2nd)'] = df_res_major['-logP(2nd)']
        
        dict_summary_res['Major_Type'] = df_res_major
        dict_summary_score['Major_Type'] = df_score_major
        dict_summary_score2['Major_Type'] = df_GSA_score_major
        
    #'''
    ##### Minor type identification for selected cells
    if ident_level > 1:
                    
        df_pred['cell_type_minor(1st)'] = df_pred['cell_type(1st)']
        df_pred['Confidence_minor(1st)'] = df_pred['Confidence(1st)']
        df_pred['cell_type_minor(2nd)'] = df_pred['cell_type(1st)']
        df_pred['Confidence_minor(2nd)'] = df_pred['Confidence(1st)']
        
        ## Cluster basis target cells Selection 
        b_sel = y_pred_correct != 'unassigned'
        b_sel[:] = False
        cluster_lst = list(set(y_clust))
        for clst in cluster_lst:
            bc = y_clust == clst
            ba = y_pred_correct == 'unassigned'
            if np.sum(bc&ba) <= np.sum(bc)*pmaj:
                b_sel[bc] = True

        ## Run GSA and perform identification
        X_sel = Xs.loc[b_sel,:]
        X_pca_sel = X_pca.loc[b_sel,:]
        y_clust_sel = y_clust[b_sel]
        y_pred_sel = y_pred_correct[b_sel]

        ## Get cell-type lists
        target_cell_lst = {} # target_cell_types
        for c in mkr_lst.keys():
            if c in dict_celltype_comb.keys():
                target_cell_lst[c] = dict_celltype_comb[c]
            elif c in target_cell_types:
                target_cell_lst[c] = [c]

        mkr_lst = {}
        mkr_lst_neg = {}
        for c in target_cell_lst.keys():
            target_cell = target_cell_lst[c]
            mkr_lst_tmp, mkr_lst_neg_tmp = get_markers_cell_type(marker_file, target_cell, 
                                            pnsh12 = mkr_selector, rem_common = rem_cmn_mkr,
                                                                 verbose = False)
            mkr_lst.update(mkr_lst_tmp)
            mkr_lst_neg.update(mkr_lst_neg_tmp)

        if len(minor_types_to_exclude) > 0:
            mkr_lst = rem_minor( mkr_lst, minor_types_to_exclude )
            mkr_lst_neg = rem_minor( mkr_lst_neg, minor_types_to_exclude )
        
        ### To use in the subset identification
        mkr_lst_sav = copy.deepcopy(mkr_lst)

        if not minor_id_sep:
            #'''
            sm = ''
            for key in mkr_lst.keys():
                sm = sm + '%s,' % key
            if verbose: print('Minor type identification: %s' % (sm[:-1]))

            df_res_sel, df_score_sel, df_GSA_score_sel = \
                run_gsa_and_clf(X_pca_sel, X_sel, cobj, y_clust_sel, mkr_lst, mkr_lst_neg, 
                                 method = None, N_components = N_components, pval_th = pval_th, 
                                 pct_margin_multi = pth_margin_mult, pct_min_rf = pct_th_min,
                                 Target_FPR = Target_FPR, pmaj = 0, minor_id_sep = minor_id_sep,
                                 min_logP_diff = min_logP_diff, thresholding = False,
                                 print_report = verbose)

            y_pred_sel = df_res_sel['cell_type']
            ys_sel = df_res_sel['cell_type(rev)']

            if (N_neighbors_minor > 0) & (len(y_pred_sel) >= N_neighbors_minor):
                if verbose:
                    print('Applying KNN rule to correct minor type .. ', end = '')

                ## KNN correction
                nbrs = NearestNeighbors(n_neighbors= int(N_neighbors_minor), algorithm='ball_tree').fit(X_pca_sel)
                distances, nei_indices = nbrs.kneighbors(X_pca_sel)
                nei_indices = pd.DataFrame( nei_indices, index = X_pca_sel.index.values )

                gsa_score = df_res_sel['-logP']
                y_pred_sel = apply_knn(nei_indices, y_pred_sel, gsa_score, pmaj = 0.5, nlog_pval_th = -np.log10(pval_th))        
                if verbose: print('done. (%i) ' % (round(time.time() - start_time)))
                
            df_pred.loc[b_sel, 'cell_type_minor'] = y_pred_sel
            
            df_pred.loc[b_sel, 'cell_type_minor(1st)'] = df_res_sel['cell_type(1st)']
            df_pred.loc[b_sel, 'Confidence_minor(1st)'] = df_res_sel['-logP']
            df_pred.loc[b_sel, 'cell_type_minor(2nd)'] = df_res_sel['cell_type(2nd)']
            df_pred.loc[b_sel, 'Confidence_minor(2nd)'] = df_res_sel['-logP(2nd)']
        
            dict_summary_res['Minor_Type'] = df_res_sel
            dict_summary_score['Minor_Type'] = df_score_sel
            dict_summary_score2['Minor_Type'] = df_GSE_score_sel

            if verbose:
                print('Minor type identification done. (%i)' % round(time.time() - start_time))
                # plot_roc_result(df_score_sel, ys_sel, method)
            #'''
        else:
            idx_sel = X_pca_sel.index.values

            cnt = 0
            for tc in target_cell_lst.keys():
                target_cell = target_cell_lst[tc]
                mkr_lst, mkr_lst_neg = get_markers_cell_type(marker_file, target_cell, 
                                                pnsh12 = mkr_selector, rem_common = rem_cmn_mkr,
                                                             verbose = False)

                if len(minor_types_to_exclude) > 0:
                    mkr_lst = rem_minor( mkr_lst, minor_types_to_exclude )
                    mkr_lst_neg = rem_minor( mkr_lst_neg, minor_types_to_exclude )
        
                sm = ''
                for key in mkr_lst.keys():
                    sm = sm + '%s,' % key
                if verbose: print('%s minor type identification: %s' % (tc, sm[:-1]))

                b_cur = y_pred_sel == 'unassigned'
                b_cur[:] = False
                cluster_lst = list(set(y_clust_sel))

                for clst in cluster_lst:
                    bc = y_clust_sel == clst
                    ba = y_pred_sel == tc
                    if np.sum(bc&ba) >= np.sum(bc)*pmaj:
                        b_cur[bc] = True

                ## Run GSA and perform identification
                X_cur = X_sel.loc[b_cur,:]
                X_pca_cur = X_pca_sel.loc[b_cur,:]
                y_clust_cur = y_clust_sel[b_cur]

                df_res_cur, df_score_cur, df_GSA_score_cur = \
                    run_gsa_and_clf(X_pca_cur, X_cur, cobj, y_clust_cur, mkr_lst, mkr_lst_neg, 
                                     method = None, N_components = N_components, pval_th = pval_th, 
                                     pct_margin_multi = pth_margin_mult, pct_min_rf = pct_th_min,
                                     Target_FPR = Target_FPR, pmaj = 0, minor_id_sep = minor_id_sep,
                                     min_logP_diff = min_logP_diff, thresholding = False,
                                     print_report = verbose)
                
                ######## Markercount_Ref ########
                bx = df_res_cur['>=th']
                if (ref_correct) & (len(list(df_res_cur.loc[bx, 'cell_type'].unique())) > 1):                
                    X_ref = X_cur.loc[bx,:].copy(deep = True)
                    cell_type_ref = df_res_cur.loc[bx, 'cell_type'].copy(deep = True)
                    X_test = X_cur.copy(deep = True)

                    if verbose: print('MarkerCount_Ref .. ', end = '')
                    cells_to_exclude = ['Unknown', 'Tumor']
                    df_res_ref, df_mkr_ref, df_param, df_score_ref = \
                        MkrCnt_Ref( X_ref = X_ref, cell_type_ref = cell_type_ref, X_test = X_test, \
                            df_mkr_mat = None, N_mkrs = 18, of_th = 0.85, \
                            cell_types_to_excl = [], log_transformed = False, \
                            N_clusters = None, cluster_label = None, p_maj = 1, \
                            X_pca = None, N_pca = 15, target_FPR = 0.05, file_to_save_marker = None, \
                            verbose = False )
                    if verbose: print('done.')

                    ## The identification results are in df_res['cell_type_pred']
                    # df_res_cur['cell_type'] = df_res_ref['cell_type_pred']  
                    df_res_ref = get_stat(df_score_ref)
                    df_res_cur['cell_type'] = df_res_ref['cell_type']  
                    df_res_cur['cell_type(2nd)'] = df_res_ref['cell_type(2nd)']  
                else:
                    df_score_ref = None
                    
                #################################

                y_pred_cur = df_res_cur['cell_type']
                ys_cur = df_res_cur['cell_type(rev)']

                if (N_neighbors_minor > 0) & (len(y_pred_cur) >= N_neighbors_minor):
                    if verbose:
                        print('Applying KNN rule to correct minor type .. ', end = '')

                    ## KNN correction
                    nbrs = NearestNeighbors(n_neighbors= int(N_neighbors_minor), algorithm='ball_tree').fit(X_pca_cur)
                    distances, nei_indices = nbrs.kneighbors(X_pca_cur)
                    nei_indices = pd.DataFrame( nei_indices, index = X_pca_cur.index.values )

                    gsa_score = df_res_cur['-logP']
                    y_pred_cur = apply_knn(nei_indices, y_pred_cur, gsa_score, pmaj = 0.5, nlog_pval_th = -np.log10(pval_th))        
                    if verbose: print('done. (%i) ' % (round(time.time() - start_time)))
                    
                    
                idx = idx_sel[b_cur]
                df_pred.loc[idx,'cell_type_minor'] = y_pred_cur
                y_pred_sel[b_cur] = y_pred_cur
                
                df_pred.loc[idx, 'cell_type_minor(1st)'] = df_res_cur['cell_type(1st)']
                df_pred.loc[idx, 'Confidence_minor(1st)'] = df_res_cur['-logP']
                df_pred.loc[idx, 'cell_type_minor(2nd)'] = df_res_cur['cell_type(2nd)']
                df_pred.loc[idx, 'Confidence_minor(2nd)'] = df_res_cur['-logP(2nd)']
                                
                if df_res_cur.shape[0] > 0:
                    dict_summary_res[tc + ' minor type'] = df_res_cur
                    dict_summary_score[tc + ' minor type'] = df_score_cur
                    dict_summary_score2[tc + ' minor type'] = df_GSA_score_cur
                    if df_score_ref is not None:
                        dict_summary_score3[tc + ' minor type'] = df_score_ref

                if cnt == 0:
                    df_res_sel = df_res_cur
                    # if df_score_cur is not None:
                    df_score_sel = df_score_cur
                    df_GSA_score_sel = df_GSA_score_cur
                else:
                    df_res_sel = pd.concat([df_res_sel, df_res_cur], axis = 0)
                    if (df_score_cur is not None) & (df_score_sel is not None):
                        df_score_sel = pd.concat([df_score_sel, df_score_cur], axis = 0)
                    df_GSA_score_sel = pd.concat([df_GSA_score_sel, df_GSA_score_cur], axis = 0)
                cnt += 1

                if verbose:
                    print('%s minor type identification done. (%i)' % \
                          (tc, round(time.time() - start_time)))
                    # plot_roc_result(df_score_sel, ys_sel, method)
                
            if np.sum(b_sel) != df_res_sel.shape[0]:
                print('WARNING: %i != %i. One or more cluster(s) is mixed.' % \
                      (np.sum(b_sel), df_res_sel.shape[0]))
            if verbose:
                print('%s type identification done. (%i)' % (tc, round(time.time() - start_time)))
                # plot_roc_result(df_score_sel, ys_sel, method)
            #'''

        '''
        if N_neighbors > 0:
            if verbose:
                print('Applying KNN rule to correct minor type .. ', end = '')

            ## KNN correction
            nbrs = NearestNeighbors(n_neighbors= int(N_neighbors), algorithm='ball_tree').fit(X_pca_sel)
            distances, nei_indices = nbrs.kneighbors(X_pca_sel)
            nei_indices = pd.DataFrame( nei_indices, index = X_pca_sel.index.values )

            gsa_score = df_res_sel['-logP']
            y_pred_sel = apply_knn(nei_indices, y_pred_sel, gsa_score, pmaj = 0.5, nlog_pval_th = -np.log10(pval_th))        

            if verbose: print('done. (%i) ' % (round(time.time() - start_time)))
        '''    
        df_pred.loc[b_sel, 'cell_type_minor'] = y_pred_sel
        
        # dict_summary_res['Minor_Type'] = df_res_sel
        # dict_summary_score['Minor_Type'] = df_score_sel
        # dict_summary_score2['Minor_Type'] = df_GSA_score_sel
            
        b_cur = y_pred_sel == 'unassigned'
        if verbose:
            print('Num of unassigned cells: %i among %i' % (np.sum(b_cur), len(b_cur)))

        ## make correction if major type and minor type does not match
        cell_type_map_dict = {}
        for tc in target_cell_lst.keys():
            tc_lst = target_cell_lst[tc]
            mkr_lst, mkr_lst_neg = get_markers_cell_type(marker_file, tc_lst, pnsh12 = mkr_selector, 
                                            rem_common = rem_cmn_mkr, verbose = False)

            cell_type_lst = list(mkr_lst.keys())
            cell_type_map_dict[tc] = cell_type_lst

        #'''
        map_dict = {}
        for key in list(cell_type_map_dict.keys()):
            for c in cell_type_map_dict[key]:
                map_dict[c] = key
        cell_types = list(map_dict.keys())

        y_pred_major = list(df_pred['cell_type_major'])
        y_pred_minor = list(df_pred['cell_type_minor'])

        y_pred_major_new = []
        y_pred_minor_new = []
        for yj, ym in zip(y_pred_major, y_pred_minor):
            if ym in cell_types:
                if map_dict[ym] == yj:
                    y_pred_major_new.append(yj)
                    y_pred_minor_new.append(ym)
                else:
                    # y_pred_major_new.append('unassigned')
                    # y_pred_minor_new.append('unassigned')
                    y_pred_major_new.append(map_dict[ym])
                    y_pred_minor_new.append(ym)
            else: # 'unassigned'
                y_pred_major_new.append('unassigned')
                y_pred_minor_new.append('unassigned')

        df_pred['cell_type_major'] = y_pred_major_new
        df_pred['cell_type_minor'] = y_pred_minor_new
        #'''
        '''
        if verbose:
            b_cur2 = y_pred_sel == 'unassigned'
            if np.sum(b_cur2) != np.sum(b_cur):
                print('%i cells among %i were filtered out.' % \
                      (np.sum(b_cur2)-np.sum(b_cur), len(b_cur)-np.sum(b_cur)))       
        #'''
        
    ### Subset identification
    if ident_level > 2:        

        df_pred['cell_type_subset(1st)'] = df_pred['cell_type_minor(1st)']
        df_pred['Confidence_subset(1st)'] = df_pred['Confidence_minor(1st)']
        df_pred['cell_type_subset(2nd)'] = df_pred['cell_type_minor(1st)']
        df_pred['Confidence_subset(2nd)'] = df_pred['Confidence_minor(1st)']
        
        cnt = 0
        for tc in target_cell_lst.keys():
            
            target_cell = target_cell_lst[tc]
            mkr_lst, mkr_lst_neg = get_markers_cell_type(marker_file, target_cell, 
                                            pnsh12 = mkr_selector, rem_common = rem_cmn_mkr,
                                                         verbose = False)
                
            if len(minor_types_to_exclude) > 0:
                mkr_lst = rem_minor( mkr_lst, minor_types_to_exclude )
                mkr_lst_neg = rem_minor( mkr_lst_neg, minor_types_to_exclude )
        
            cell_type_lst = list(mkr_lst.keys())
            b_cur = np.full(len(y_pred_sel), False)
            for c in cell_type_lst:
                b_cur = b_cur | (y_pred_sel == c)

            for c in cell_type_lst: # for each minor type
                b_cur = y_pred_sel == c
                mkr_lst, mkr_lst_neg = get_markers_minor_type2(marker_file, [c], 
                                    pnsh12 = mkr_selector, rem_common = rem_cmn_mkr, verbose = False)

                if len(minor_types_to_exclude) > 0:
                    mkr_lst = rem_minor( mkr_lst, minor_types_to_exclude )
                    mkr_lst_neg = rem_minor( mkr_lst_neg, minor_types_to_exclude )
        
                sm = ''
                for key in mkr_lst.keys():
                    sm = sm + '%s,' % key
                if verbose: 
                    print('%s subset identification: %s' % (c, sm[:-1]))

                if (len(list(mkr_lst.keys())) <= 1) | (np.sum(b_cur) < 10):
                    X_cur = X_sel.loc[b_cur,:]
                    idx = list(X_cur.index.values)
                    ctn = c
                    if len(list(mkr_lst.keys())) == 1:
                        ctn = list(mkr_lst.keys())[0]
                    df_pred.loc[idx,'cell_type_subset'] = ctn
                    
                else:
                    X_cur = X_sel.loc[b_cur,:]
                    X_pca_cur = X_pca_sel.loc[b_cur,:]
                    y_clust_cur = y_clust_sel[b_cur]

                    # print(mkr_lst)
                    df_res_cur, df_score_cur, df_GSA_score_cur = \
                        run_gsa_and_clf(X_pca_cur, X_cur, cobj, y_clust_cur, mkr_lst, mkr_lst_neg, 
                                         method = None, N_components = N_components, pval_th = pval_th, 
                                         pct_margin_multi = pth_margin_mult, pct_min_rf = pct_th_min,
                                         Target_FPR = Target_FPR, pmaj = 0, minor_id_sep = minor_id_sep,
                                         min_logP_diff = min_logP_diff, thresholding = False,
                                         print_report = verbose)

                    cts = list(mkr_lst.keys())
                    if len(cts) > 1:
                        bx = df_res_cur['-logP'] < -np.log10( pval_th2 )
                        df_res_cur.loc[bx, 'cell_type'] = 'unassigned'     
                        pass
                    
                    ######## Markercount_Ref ########
                    bx = df_res_cur['>=th']
                    if (ref_correct) & (len(list(df_res_cur.loc[bx, 'cell_type'].unique())) > 1):
                        X_ref = X_cur.loc[bx,:].copy(deep = True)
                        cell_type_ref = df_res_cur.loc[bx, 'cell_type'].copy(deep = True)
                        X_test = X_cur.copy(deep = True)

                        if verbose:
                            print('MarkerCount_Ref .. ', end = '')
                        cells_to_exclude = ['Unknown', 'Tumor']
                        df_res_ref, df_mkr_ref, df_param, df_score_ref = \
                            MkrCnt_Ref( X_ref = X_ref, cell_type_ref = cell_type_ref, \
                                X_test = X_test, \
                                df_mkr_mat = None, N_mkrs = 18, of_th = 0.85, \
                                cell_types_to_excl = [], log_transformed = False, \
                                N_clusters = None, cluster_label = None, p_maj = 0.2, \
                                X_pca = None, N_pca = 15, target_FPR = 0.05, \
                                file_to_save_marker = None, \
                                verbose = False )
                        if verbose:
                            print('done.')

                        ## The identification results are in df_res['cell_type_pred']
                        # df_res_cur['cell_type'] = df_res_ref['cell_type_pred']  
                        df_res_ref = get_stat(df_score_ref)
                        df_res_cur['cell_type'] = df_res_ref['cell_type']  
                        df_res_cur['cell_type(2nd)'] = df_res_ref['cell_type(2nd)']  
                    else:
                        df_score_ref = None

                    #################################
                
                    y_pred_cur = df_res_cur['cell_type']
                    ys_cur = df_res_cur['cell_type(rev)']

                    if (N_neighbors_subset > 2) & (len(y_pred_cur) >= N_neighbors_subset):
                        if verbose:
                            print('Applying KNN rule to correct cell type subset .. ', end = '')

                        ## KNN correction
                        nbrs = NearestNeighbors(n_neighbors= int(N_neighbors_subset), algorithm='ball_tree').fit(X_pca_cur)
                        distances, nei_indices = nbrs.kneighbors(X_pca_cur)
                        nei_indices = pd.DataFrame( nei_indices, index = X_pca_cur.index.values )

                        gsa_score = df_res_cur['-logP']
                        y_pred_cur = apply_knn(nei_indices, y_pred_cur, gsa_score, pmaj = 0.5, 
                                               nlog_pval_th = -np.log10(pval_th))        
                        if verbose: print('done. (%i) ' % (round(time.time() - start_time)))
                    
                    
                    if verbose:
                        print('%s subset identification done. (%i)' % \
                              (c, round(time.time() - start_time)))
                        # plot_roc_result(df_score_cur, ys_cur, method)

                    df_res_cur['cell_type'] = y_pred_cur

                    idx = list(X_cur.index.values)
                    df_pred.loc[idx,'cell_type_subset'] = y_pred_cur
                    
                    df_pred.loc[idx, 'cell_type_subset(1st)'] = df_res_cur['cell_type(1st)']
                    df_pred.loc[idx, 'Confidence_subset(1st)'] = df_res_cur['-logP']
                    df_pred.loc[idx, 'cell_type_subset(2nd)'] = df_res_cur['cell_type(2nd)']
                    df_pred.loc[idx, 'Confidence_subset(2nd)'] = df_res_cur['-logP(2nd)']
                
                    if df_res_cur.shape[0] > 0:
                        dict_summary_res['%s subset' % c] = df_res_cur
                        dict_summary_score['%s subset' % c] = df_score_cur
                        dict_summary_score2['%s subset' % c] = df_GSA_score_cur
                        if df_score_ref is not None:
                            dict_summary_score3['%s subset' % c] = df_score_ref

        b_cur = y_pred_correct == 'unassigned'
        if verbose:
            print('Num of unassigned cells: %i among %i' % (np.sum(b_cur), len(b_cur)))
    
    ## end if minor_subtype   
    #'''
    if cycling_cell & ('MKI67' in list(Xs.columns.values)):
        b = Xs['MKI67'] > 0
        
        y_new = []
        # y_pred_major = list(df_pred.loc[b, 'cell_type_minor'])
        y_pred_major = list(df_pred.loc[b, 'cell_type_major'])
        for y in y_pred_major:
            if y != 'unassigned':
                y_new.append('%s Cycling' % y)
            else:
                y_new.append(y)
                
        df_pred.loc[b, 'cell_type_minor'] = y_new
        df_pred.loc[b, 'cell_type_subset'] = y_new
    #'''
                
    if verbose: 
        print('MarkerCount_GSA done. (%i)' % round(time.time() - start_time))
    else:
        print(' done. (%i)' % round(time.time() - start_time))
        
    dict_summaries = {}
    dict_summaries['GSA_summary'] = dict_summary_res
    dict_summaries['GSA_scores'] = dict_summary_score2
    dict_summaries['Ref_scores'] = dict_summary_score3
    dict_summaries['Identification_model_scores'] = dict_summary_score
    dict_summaries['parameters'] = [pval_th, pth_margin_mult, pct_th_min]
    
       
    return df_pred, dict_summaries

