import pandas as pd 
import os
import sys 
from os.path import join,dirname,exists
from ipdb import set_trace
from os import system
# from hetio.downloader import download, unzip
# from hetio.datasets import datasets_url
# from downloader import download
import json 
# pd.set_option("display.max_rows",1000)

datasets_url = 'http://192.168.0.97/share/StandigmDB/datasets'


def download(remote_dir=None, filename=None, local_dir=None, force=False):

    base_dir = dirname(__file__)

    if local_dir is None:
        local_dir = join(base_dir, 'scratch')

    os.system('mkdir -p %s' % local_dir)
    remotepath = join(remote_dir, filename)
    savepath = join(local_dir, filename)

    if not exists(savepath) or force==True: 
        os.system('wget -O %s %s' % (savepath, remotepath))

    return savepath


def load():
    # All data is downloaded from 
    # https://www.ebi.ac.uk/gwas/docs/file-downloads
    
    remote_dir = join(datasets_url, 'gwas', '2017-11')

    # All associations v1.0
    catalog_ver_1_0 = download(remote_dir, 
        'gwas_catalog_v1.0-associations_e90_r2017-11-13.tsv', 
        force=False
        )

    # All associations v1.0.1 - 
    # with added ontology annotations and GWAS Catalog study accession numbers
    catalog_associations_ver_1_0_1 = download(remote_dir, 
        'gwas_catalog_v1.0.1-associations_e90_r2017-11-13.tsv', 
        force=False
        )
    
    # All ancestry data (see FAQ A13 for more details)
    ancestry = download(remote_dir, 
        'gwas_catalog-ancestry_r2017-11-13.tsv', 
        force=False
        )

    # All associations v1.0.1 - 
    # with added ontology annotations and GWAS Catalog study accession numbers
    catalog_studies_ver_1_0_1  = download(remote_dir, 
        'gwas_catalog_v1.0.1-studies_r2017-11-13.tsv', 
        force=False
        ) 

    # All studies v1.0
    catalog_studies_ver_1_0 = download(remote_dir, 
        'gwas_catalog_v1.0-studies_r2017-11-13.tsv', 
        force=False
        )

    df_catalog_ver_1_0 = pd.read_csv(
        catalog_ver_1_0, sep='\t'
        )
    # ipdb> dataset['catalog-v1.0'].columns
    # Index(['DATE ADDED TO CATALOG', 'PUBMEDID', 'FIRST AUTHOR', 'DATE', 'JOURNAL',
    #    'LINK', 'STUDY', 'DISEASE/TRAIT', 'INITIAL SAMPLE SIZE',
    #    'REPLICATION SAMPLE SIZE', 'REGION', 'CHR_ID', 'CHR_POS',
    #    'REPORTED GENE(S)', 'MAPPED_GENE', 'UPSTREAM_GENE_ID',
    #    'DOWNSTREAM_GENE_ID', 'SNP_GENE_IDS', 'UPSTREAM_GENE_DISTANCE',
    #    'DOWNSTREAM_GENE_DISTANCE', 'STRONGEST SNP-RISK ALLELE', 'SNPS',
    #    'MERGED', 'SNP_ID_CURRENT', 'CONTEXT', 'INTERGENIC',
    #    'RISK ALLELE FREQUENCY', 'P-VALUE', 'PVALUE_MLOG', 'P-VALUE (TEXT)',
    #    'OR or BETA', '95% CI (TEXT)', 'PLATFORM [SNPS PASSING QC]', 'CNV'],
    #   dtype='object')
    df_catalog_associations_ver_1_0_1 = pd.read_csv(
        catalog_associations_ver_1_0_1, sep='\t'
        )
    # ipdb> dataset['catalog-associations-v1.0.1'].columns
    # Index(['DATE ADDED TO CATALOG', 'PUBMEDID', 'FIRST AUTHOR', 'DATE', 'JOURNAL',
    #    'LINK', 'STUDY', 'DISEASE/TRAIT', 'INITIAL SAMPLE SIZE',
    #    'REPLICATION SAMPLE SIZE', 'REGION', 'CHR_ID', 'CHR_POS',
    #    'REPORTED GENE(S)', 'MAPPED_GENE', 'UPSTREAM_GENE_ID',
    #    'DOWNSTREAM_GENE_ID', 'SNP_GENE_IDS', 'UPSTREAM_GENE_DISTANCE',
    #    'DOWNSTREAM_GENE_DISTANCE', 'STRONGEST SNP-RISK ALLELE', 'SNPS',
    #    'MERGED', 'SNP_ID_CURRENT', 'CONTEXT', 'INTERGENIC',
    #    'RISK ALLELE FREQUENCY', 'P-VALUE', 'PVALUE_MLOG', 'P-VALUE (TEXT)',
    #    'OR or BETA', '95% CI (TEXT)', 'PLATFORM [SNPS PASSING QC]', 'CNV',
    #    'MAPPED_TRAIT', 'MAPPED_TRAIT_URI', 'STUDY ACCESSION'],
    #   dtype='object')
    df_ancestry = pd.read_csv(
        ancestry, sep='\t'
        )
    # Index(['STUDY ACCCESSION', 'PUBMEDID', 'FIRST AUTHOR', 'DATE',
    #        'INITIAL SAMPLE DESCRIPTION', 'REPLICATION SAMPLE DESCRIPTION', 'STAGE',
    #        'NUMBER OF INDIVDUALS', 'BROAD ANCESTRAL CATEGORY', 'COUNTRY OF ORIGIN',
    #        'COUNTRY OF RECRUITMENT', 'ADDITONAL ANCESTRY DESCRIPTION'],
    #   dtype='object')    
    df_catalog_studies_ver_1_0_1 = pd.read_csv(
        catalog_studies_ver_1_0_1, sep='\t'
        )
    # ipdb> dataset['catalog-studies-v1.0.1'].columns
    # Index(['DATE ADDED TO CATALOG', 'PUBMEDID', 'FIRST AUTHOR', 'DATE', 'JOURNAL',
    #        'LINK', 'STUDY', 'DISEASE/TRAIT', 'INITIAL SAMPLE SIZE',
    #        'REPLICATION SAMPLE SIZE', 'PLATFORM [SNPS PASSING QC]',
    #        'ASSOCIATION COUNT', 'MAPPED_TRAIT', 'MAPPED_TRAIT_URI',
    #        'STUDY ACCESSION'],
    #       dtype='object')
    
    df_catalog_studies_ver_1_0 = pd.read_csv(
        catalog_studies_ver_1_0, sep='\t'
        )
    # ipdb> dataset['catalog-studies-v1.0'].columns
    # Index(['DATE ADDED TO CATALOG', 'PUBMEDID', 'FIRST AUTHOR', 'DATE', 'JOURNAL',
    #        'LINK', 'STUDY', 'DISEASE/TRAIT', 'INITIAL SAMPLE SIZE',
    #        'REPLICATION SAMPLE SIZE', 'PLATFORM [SNPS PASSING QC]',
    #        'ASSOCIATION COUNT'],
    #       dtype='object')
    dataset = {
        'catalog-v1.0': 
            df_catalog_ver_1_0,
        'catalog-associations-v1.0.1': 
            df_catalog_associations_ver_1_0_1, 
        'ancestry': 
            df_ancestry, 
        'catalog-studies-v1.0.1': 
            df_catalog_studies_ver_1_0_1, 
        'catalog-studies-v1.0': 
            df_catalog_studies_ver_1_0
        }

    # set_trace()

    example = {
        'catalog-v1.0': 
            dataset['catalog-v1.0'].head(1).to_json(orient='records'), 
        'catalog-associations-v1.0.1': 
            dataset['catalog-associations-v1.0.1'].head(1).to_json(orient='records'),
        'ancestry': 
            dataset['ancestry'].head(1).to_dict('index'),
        'catalog-studies-v1.0.1': 
            dataset['catalog-studies-v1.0.1'].head(1).to_json(orient='records'),
        'catalog-studies-v1.0': 
            dataset['catalog-studies-v1.0'].head(1).to_json(orient='records'),  
        }

    # set_trace()

    dataset['example'] = example
    # set_trace()

    # import json 
    with open('example.json', 'w') as f: 
        json.dump(example, f, indent=4) 

    return dataset 


def test_load():
    dataset = load() 

    # set_trace()

