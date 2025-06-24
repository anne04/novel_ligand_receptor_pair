import pandas as pd
import json 
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict


def print_command(lrp_list, dict_gene_seq, prefix, path_to, AF_output_dir):
    for pair in lrp_list:
        ligand = pair[0]
        receptor = pair[1]
        lig_seq = dict_gene_seq[ligand]
        rec_seq = dict_gene_seq[receptor]
        f = open(path_to + prefix + ligand +'_'+receptor+".fasta", "w")
        f.write(">")
        f.write(ligand)
        f.write("\n")
        f.write(lig_seq)
        f.write("\n")
        f.write(">")
        f.write(receptor)
        f.write("\n")
        f.write(rec_seq)
        f.close()
        fasta = path_to + prefix + ligand + '_' + receptor + '.fasta'
        
        print('bash run_alphafold.sh -d ${DOWNLOAD_DIR} -o ' + AF_output_dir + ' -m model_1_multimer_v3  -p multimer -i ' + fasta + ' -t 2022-01-01 -r \'none\' -c reduced_dbs')
        print('mv ' + AF_output_dir + prefix + ligand+'_'+receptor+'/'+'ranking_debug.json '+ AF_output_dir + prefix + ligand +'_'+ receptor +'_score.json')
        print('rm -r '+ AF_output_dir + prefix + ligand+'_'+receptor)

if __name__ == "__main__":

    save_json_to = '/cluster/projects/schwartzgroup/fatema/LRbind/alphafold_input/'
    prefix = 'lrbind_'
    lr_list_file = '/cluster/home/t116508uhn/LRbind_output/without_elbow_cut/LRbind_V1_Human_Lymph_Node_spatial_1D_manualDB_geneLocalCorrKNN_bidir_prefiltered/model_LRbind_V1_Human_Lymph_Node_spatial_1D_manualDB_geneLocalCorrKNN_bidir_3L_prefiltered_down_up_deg_lr_list_sortedBy_totalScore_top_elbow_allLR.csv'
    #model_LRbind_V1_Human_Lymph_Node_spatial_1D_manualDB_geneLocalCorrKNN_bidir_3L_prefiltered_down_up_deg_novel_lr_list_sortedBy_totalScore_top_elbow_novelsOutOfallLR.csv'
    AF_output_dir = 'ParallelFold-main/output/'
    from_pair = 0
    to_pair = 30
    filter ='db_only' 
    marker = '+'
    import argparse
    parser = argparse.ArgumentParser()
    # ================ Specify data type firstly ===============
    parser.add_argument( '--file_name', type=str, default='uniprotkb_reviewed_true_AND_proteome_up_2025_02_27.tsv', help='The name of DB')
    parser.add_argument( '--database_path', type=str, default='database/NEST_database_no_predictedPPI.csv', help='The name of DB')
    parser.add_argument( '--result_path', type=str, default='result/')
    args = parser.parse_args()

    df = pd.read_csv(args.file_name, sep="\t")
    dict_gene_seq = dict()
    for i in range (0, df['Sequence'].shape[0]):
        if not isinstance(df['Gene Names'][i], str):
            continue            
        dict_gene_seq[(df['Gene Names'][i]).split(' ')[0]] = df['Sequence'][i]
        
    ##############################################################################


    probable_pairs = []
    df = pd.read_csv(lr_list_file, sep=",")
    to_pair = min(to_pair, len( df["Ligand-Receptor Pairs"]))
    for i in range (from_pair, to_pair):
        if filter == 'db_only' and df["Type"][i]!='From DB':
            continue
            
        ligand = df["Ligand-Receptor Pairs"][i].split(marker)[0]
        receptor = df["Ligand-Receptor Pairs"][i].split(marker)[1]       
        probable_pairs.append([ligand, receptor])

    
    ### see which of those are not calculated yet #############
    file_list = glob.glob(AF_output_dir+"*json")
    calculated_pairs = dict()
    for file_path in file_list:
        ligand = file_path.split('_')[1]
        receptor = file_path.split('_')[2]
        calculated_pairs[ligand + '+' + receptor] = ''


    lrp_list_to_run = []
    for pair in probable_pairs:
        ligand = pair[0]
        receptor = pair[1]
        if ligand + '+' + receptor not in calculated_pairs:
            lrp_list_to_run.append([ligand, receptor]) 

    
    print('Going to run AF on %d pairs'%len(lrp_list_to_run))
    print_command(lrp_list_to_run, dict_gene_seq, prefix, save_json_to, AF_output_dir)

