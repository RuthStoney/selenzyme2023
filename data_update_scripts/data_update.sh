#!/usr/bin/env bash

# #### File requirements
# # raw_data_folder
# reac_prop.tsv		https://www.metanetx.org/mnxdoc/mnxref.html
# chem_prop.tsv		https://www.metanetx.org/mnxdoc/mnxref.html

# uniprot_sprot.fasta
# brenda_2023_1.txt	https://www.brenda-enzymes.org/download.php
# expasy_dat.txt		https://ftp.expasy.org/databases/enzyme/
# taxidlineage.dmp	https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/
# names.dmp			https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/

# # current dataset 
# reac_seqs.tsv

#### Software requirements
# RXNMapper  https://github.com/rxn4chemistry/rxnmapper


NEW_DATA_RAW=/home/ruth/code/update_selenzyme/run_folder_min_dec/raw_data_update/
OLD_DATA=/home/ruth/code/update_selenzyme/run_folder_min_dec/data/
NEW_DATA=/home/ruth/code/update_selenzyme/run_folder_min_dec/data_2023_t/


echo "\n     Filter_reactions run one"
python filter_reactions.py $NEW_DATA $NEW_DATA_RAW

echo "\n     Make reac_seq"
python make_reac_seq_from_brenda_expasy.py $NEW_DATA $NEW_DATA_RAW $OLD_DATA

echo "\n     Make fingerprints"
# requires RXNMapper
python make_fingerprint_atomMap.py $NEW_DATA $NEW_DATA_RAW

echo "\n     Make seq_org"
python make_seq_org_fasta_uniprotAPI.py $NEW_DATA $NEW_DATA_RAW $OLD_DATA

echo "\n     Make org_linage"
python make_org_lineage.py $NEW_DATA $NEW_DATA_RAW

echo "\n     Filter_reactions run two"
python filter_reactions.py $NEW_DATA $NEW_DATA_RAW

echo "\n     Make reac_smiles"
python make_reac_smiles.py $NEW_DATA $NEW_DATA_RAW

echo "\n     Make reac_xref"
python make_reac_xref.py $NEW_DATA $NEW_DATA_RAW

cp $NEW_DATA_RAW"uniprot_sprot.fasta" $NEW_DATA"seqs.fasta"
cp $NEW_DATA"Morgan/FP_Morg.npz" $NEW_DATA"FP_Morg.npz"
cp $NEW_DATA"Morgan/RF/FP_MorgRF.npz" $NEW_DATA"FP_MorgRF.npz"


echo "\n     Update complete!"
echo $NEW_DATA
