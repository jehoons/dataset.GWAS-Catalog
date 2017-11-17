import collections
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
from ipdb import set_trace
from os.path import exists,dirname,join
import pandas as pd 
import numpy as np 
import pickle,json

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

    summary = {}

    # Window boundaries are determined as the furthest upstream and
    # downstream SNPs with r-squared greater than or equal to this theshold 
    r2_theshold = 0.5
    summary['r2_theshold'] = r2_theshold 

    ######################################################
    # This is search result. 
    # Read SNAP proxy search
    # input data for this result is : 
    # 'data/snap/do-slim-lead-SNPs.txt'
    # snap_df = pd.read_table(inp_proxy_snp, compression='gzip', 
    #     na_values=['N/A'])
    # snap_df = snap_df.dropna(subset=['Coordinate_HG18'])
    # snap_df.Coordinate_HG18 = snap_df.Coordinate_HG18.astype(int)
    snap_df = load_proxySNP()
    ######################################################
    # This is input data for the SNAP proxy search engine. 
    # Read SNAP proxy search input: 
    # with open(inp_lead_snp) as read_file:
    #     leads = {line.rstrip() for line in read_file}
    leads = load_leadSNP()
    # check coverage: 
    # print('Of {0} lead SNPs, {1} not found in SNAP'.format(
    #    len(leads), len(leads - set(snap_df.SNP)))
    #    )
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

    # output format: 
    # 
    # lead_snp        lead_chrom      lead_coord      lead_maf        lower_snp       lower_coord     lower_map_dist  lower_maf       upper_snp       upper_coord     upp     er_map_dist  upper_maf       span_kb span_cM
    #    2 rs1000113       chr5    150220269       0.042   rs76767593      150150036       0.04616 0.042   rs73268545      150315475       0.02666 0.042   165.44  0.07282
    #    3 rs1000579       chr4    4770395 0.442   rs7693647       4735225 0.0338  0.492   rs4689926       4784293 0.0143  0.383   49.068  0.0481

    # Number of lead SNPs that returned problematic SNAP proxies
    # print('snap_errors:', len(snap_errors))
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
   
    # set_trace()

    with open(summary_json, 'w') as fobj: 
        json.dump(summary, fobj, indent=4)

    # set_trace()


def test_coverage():
    import downloaded_201711 as mycode 
    dataset = mycode.load() 
    set_trace()
    pass 


