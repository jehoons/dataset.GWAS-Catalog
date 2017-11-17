import collections
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
from ipdb import set_trace
from os.path import exists,dirname,join
import pandas as pd 
import numpy as np 
import pickle,json

# ref. 
# https://github.com/dhimmel/gwas-catalog/blob/master/windows.ipynb
# SNAP proxy search for SNPs in LD with lead SNPs
# 1. Navigate to SNAP Proxy Search
# 2. Enter lead SNPs of interest under "Query SNPs"
# 3. Under "Search Options" set
#       SNP data set: 1000 Genomes Pilot 1
#       Population panel: CEU
#       r2 threshold: 0.5
#       Distance limit: 500
# 4. Under "Ouput Options"
#       Download to: File
#       Select "Include each query snp as a proxy for itself"
#       Select "Suppress warning messages in output"
# 5. Under "Filter By Array" make sure all are unselected
# 6. Under "Output Columns" select all besides "Associated Gene Annotations from GeneCruiser 
#       (decreases performance)"

test_root_dir = dirname(__file__)

# lead SNP - comes from where? GWAS-catalog? 
# in GWAS-catalog, diseases-SNP associated...  
inp_lead_snp = join(test_root_dir, 
    'data/snap/do-slim-lead-SNPs.txt')

inp_proxy_snp = join(test_root_dir, 
    'data/snap/do-slim-proxy-SNPs.tsv.gz')

chk_windows = join(test_root_dir, 
    'data/snap/windows.tsv')

chk_snap_errors = join(test_root_dir, 
    'data/snap/chk_snap_errors.pkl')

fig_lead_maf_vs_span_kb = join(test_root_dir, 
    'data/snap/fig-joint-dist-lead_maf-vs-span_kb.png')

fig_spankb_vs_span_cm = join(test_root_dir, 
    'data/snap/fig-span-kb-vs-span-cm.png')

fig_spankb_vs_span_cm_ver2 = join(test_root_dir, 
    'data/snap/fig-span-kb-vs-span-cm-ver2.png')

summary_json = './summary.json'

savefig_dpi = 150

def load_proxySNP():
    snap_df = pd.read_table(inp_proxy_snp, compression='gzip', 
        na_values=['N/A'])
    snap_df = snap_df.dropna(subset=['Coordinate_HG18'])
    snap_df.Coordinate_HG18 = snap_df.Coordinate_HG18.astype(int)    
    
    return snap_df

def load_leadSNP(): 
    with open(inp_lead_snp) as read_file:
        leads = {line.rstrip() for line in read_file}

    return leads 

def test_windows():
    summary = {}

    # Window boundaries are determined as the furthest upstream and
    # downstream SNPs with r-squared greater than or equal to this theshold 
    r2_theshold = 0.5
    summary['r2_theshold'] = r2_theshold 

    snap_df = load_proxySNP() # this data comes from SNP proxy search engine
    leads = load_leadSNP() # this is input data for SNP proxy search engine

    summary['#lead SNPs'] = len(leads)
    summary['#lead SNPs not in SNAP'] = len(leads - set(snap_df.SNP))

    # Condense SNAP proxy search by lead SNP and compute windows
    if exists(chk_windows) and exists(chk_snap_errors):
        window_df = pd.read_csv(chk_windows, sep='\t')
        with open(chk_snap_errors, 'rb') as fobj: 
            snap_errors = pickle.load(fobj)
    else: 
        from tqdm import tqdm
        snap_errors = list() 
        rows = list()
        for snp, group_df in tqdm(snap_df.groupby('SNP', as_index=False)):
            try:
                self_df = group_df[group_df.Proxy == snp]
                assert len(self_df) <= 1
                lead = self_df.iloc[0]
            except IndexError:
                snap_errors.append(snp)
                continue
            
            row = pd.Series()
            row['lead_snp'] = snp
            row['lead_chrom'] = lead['Chromosome']
            row['lead_coord'] = lead['Coordinate_HG18']
            row['lead_maf'] = lead['MAF']
            
            window_df = group_df[(group_df.Chromosome == row['lead_chrom']
                ) & (group_df.RSquared >= r2_theshold)]
            i_min = np.argmin(list(window_df.Coordinate_HG18))
            i_max = np.argmax(list(window_df.Coordinate_HG18))
            for side, i in ('lower', i_min), ('upper', i_max):
                window = window_df.iloc[i]
                row[side + '_snp'] = window.Proxy
                row[side + '_coord'] = window.Coordinate_HG18
                row[side + '_map_dist'] = window.GeneticMapDistance
                row[side + '_maf'] = window.MAF
            
            row['span_kb'] = (row['upper_coord'] - row['lower_coord']) / 1000
            row['span_cM'] = row['lower_map_dist'] + row['upper_map_dist']
            rows.append(row)
        # for
        with open(chk_snap_errors, 'wb') as fobj: 
            pickle.dump(snap_errors, fobj)

        window_df = pd.DataFrame(rows)
        window_df.to_csv(chk_windows, index=False, sep='\t', 
            float_format='%.5g')
    # if

    # Number of lead SNPs that returned problematic SNAP proxies
    summary['#SNAP errors'] = len(snap_errors)

    # import tabulate 
    seaborn.jointplot(window_df.lead_maf, np.sqrt(window_df.span_kb), 
        kind='hex'
        );
    plt.savefig(fig_lead_maf_vs_span_kb, dpi=savefig_dpi)

    # View window span (kilobases) versus span (centimorgans)
    seaborn.jointplot(np.log10(1 + window_df.span_kb), 
        np.log10(0.0001 + window_df.span_cM), kind='hex'
        );
    plt.savefig(fig_spankb_vs_span_cm, dpi=savefig_dpi)
    
    summary['#SpankB with zero size'] = int( sum(window_df.span_kb == 0) )
    summary['#SpancM with zero size'] = int( sum(window_df.span_cM == 0) )

    # after data cleaning, we draw figure again. 
    # Remove zero-span windows and view span (kilobases) versus span (centimorgans)
    nonzero_df = window_df[(window_df.span_kb > 0) & (
        window_df.span_cM > 0)]
    
    seaborn.jointplot(
        np.log10(nonzero_df.span_kb), 
        np.log10(nonzero_df.span_cM), kind='hex'
        );
    plt.savefig(fig_spankb_vs_span_cm_ver2, dpi=savefig_dpi)

    with open(summary_json, 'w') as fobj: 
        json.dump(summary, fobj, indent=4)


def test_coverage():
    # gwas-catalog로부터 데이터가 왔을 것으로 생각이 된다... 그런데 약간 차이가 나는구나. 

    from downloaded_201711 import test_gwas as gwas 
    dataset = gwas.load() 
    snps_from_gwas = set(dataset['catalog-v1.0']['SNPS'].unique().tolist())
    snps_from_dhimmel = load_leadSNP()
    # ipdb> snps_from_dhimmel - snps_from_gwas
    # {'rs2075800', 'rs6661053', 'rs473034', 'rs13420028', 'rs737387', 'rs28431436', 'rs10188442', 
    # 'rs6789378', 'rs6469823', 'rs492452', 'rs5751168', 'rs2456449', 'rs543686', 'rs9273012', 
    # 'rs6763848', 'rs1877455', 'rs3827735', 'rs2219968', 'rs17051310', 'rs10785581', 'rs10496288', 
    # 'rs921551', 'rs12317459', 'rs12907038', 'rs784411', 'rs10853029', 'rs2842346', 'rs12091564', 
    # 'rs10019279', 'rs7767084', 'rs4821941', 'rs11589568', 'rs7169431', 'rs13116936', 'rs16965039', 
    # 'rs11924705', 'rs2469997', 'rs1610677', 'rs200759', 'rs10853535', 'rs2902440', 'rs10277209',
    # 'rs6834555', 'rs2466032', 'rs6887846', 'rs12711517', 'rs17435444', 'rs1165669', 'rs2711721', 
    # 'rs12485321', 'rs7735940', 'rs17314229', 'rs10496289', 'rs10795917', 'rs5012808', 'rs4463179', 
    # 'rs35252396', 'rs11102807', 'rs10132579', 'rs2305016', 'rs12628051', 'rs998124', 'rs2842347', 
    # 'rs7789197', 'rs7694725', 'rs16875333', 'rs9757252', 'rs9350602', 'rs10218795', 'rs13192613', 
    # 'rs1916284', 'rs7556462', 'rs765899', 'rs17241910', 'rs13398206', 'rs1810320', 'rs1002979', 
    # 'rs7832443', 'rs12522034', 'rs28366298', 'rs4489787', 'rs7673097', 'rs12682851', 'rs735172', 
    # 'rs6470589', 'rs3752120', 'rs4238326', 'rs3789080', 'rs200752', 'rs3785856', 'rs1866967', 
    # 'rs2446581', 'rs7081208', 'rs2856683', 'rs9628987', 'rs3127599', 'rs10755578', 'rs11102800', 
    # 'rs11209003', 'rs10940579', 'rs10889676', 'rs1165668', 'rs2386661', 'rs3806872', 'rs2064689', 
    # 'rs2305795', 'rs10122902', 'rs13264970', 'rs6005451', 'rs7827545', 'rs10858047', 'rs13258681', 
    # 'rs10789230', 'rs11637980', 'rs11587400', 'rs3450', 'rs7847271', 'rs6766510', 'rs1372662', 
    # 'rs6089829', 'rs11720607', 'rs1125777', 'rs7572482', 'rs12568010', 'rs6669582', 'rs10489525', 
    # 'rs7960483', 'rs6507016', 'rs10885582', 'rs1343151', 'rs6835704', 'rs4915737', 'rs11582563', 
    # 'rs7535752', 'rs2409488', 'rs8453', 'rs731174', 'rs4455882', 'rs2261033', 'rs9331888', 'rs6452524', 
    # 'rs6917824L', 'rs2523395', 'rs17121983', 'rs11613298', 'rs7717572', 'rs6797852', 'rs76884941', 
    # 'rs111521887', 'rs2289731', 'rs7511633', 'rs6470588', 'rs2400997', 'rs1063635', 'rsA-DQB1*02:01, 
    # rs558702', 'rs11465802', 'rs4414128', 'rs7697839', 'rs9351730', 'rs3798440', 'rs9649213', 
    # 'rs290258', 'rs1004819', 'rs11209002', 'rs12567232', 'rs757369', 'rs1219414', 'rs9567349', 
    # 'rs11585926', 'rs4395927', 'rs2230754', 'rs16944141', 'rs10456809', 'rs11605083', 'rs1243647', 
    # 'rs17141741', 'rsA-DRB1*03:01, rs9275572', 'rs4455437', 'rs8057939', 'rs41453448', 'rs16867225'}
    
    # ipdb> len(snps_from_dhimmel - snps_from_gwas)
    # 181

    set_trace()
    pass 

