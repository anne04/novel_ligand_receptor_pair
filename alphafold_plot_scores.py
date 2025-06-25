import pandas as pd
import json 
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict


if __name__ == "__main__":

    save_json_to = '/cluster/projects/schwartzgroup/fatema/LRbind/alphafold_input/'
    prefix = 'lrbind_'
    lr_list_file = '/cluster/home/t116508uhn/LRbind_output/without_elbow_cut/LRbind_V1_Human_Lymph_Node_spatial_1D_manualDB_geneLocalCorrKNN_bidir_prefiltered/model_LRbind_V1_Human_Lymph_Node_spatial_1D_manualDB_geneLocalCorrKNN_bidir_3L_prefiltered_down_up_deg_lr_list_sortedBy_totalScore_top_elbow_allLR.csv'
    #model_LRbind_V1_Human_Lymph_Node_spatial_1D_manualDB_geneLocalCorrKNN_bidir_3L_prefiltered_down_up_deg_novel_lr_list_sortedBy_totalScore_top_elbow_novelsOutOfallLR.csv'
    from_dir = 'ParallelFold-main/output/'
    from_pair = 0
    to_pair = 30
    filter = 'all' #'db_only' 
    marker = '_to_'
    import argparse
    parser = argparse.ArgumentParser()
    # ================ Specify data type firstly ===============
    parser.add_argument( '--file_name', type=str, default='uniprotkb_reviewed_true_AND_proteome_up_2025_02_27.tsv', help='The name of DB')
    parser.add_argument( '--database_path', type=str, default='database/NEST_database_no_predictedPPI.csv', help='The name of DB')
    parser.add_argument( '--result_path', type=str, default='result/')
    args = parser.parse_args()

        
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
    file_list = glob.glob(from_dir+"*json")
    calculated_pairs = dict()
    for file_path in file_list:
        ligand = file_path.split('_')[1]
        receptor = file_path.split('_')[2]
        calculated_pairs[ligand + '+' + receptor] = ''


    lrp_list_to_run = dict()
    for pair in probable_pairs:
        ligand = pair[0]
        receptor = pair[1]
        if ligand + '+' + receptor in calculated_pairs:
            lrp_list_to_run[ligand +'+'+ receptor] = ''

    lrp_list_to_run = list(lrp_list_to_run.keys())
    print('Going to plot AF on %d pairs'%len(lrp_list_to_run))
    
    output_path = '/cluster/home/t116508uhn/LRbind_output/'
    plot_title = 'AlphaFold score distribution for predicted LRP lymph' #random  
    file_name = 'AF_score_distribution_predictedLRP_lymph_' #'AF_score_distribution_randomLRP_' #'AF_score_distribution_selfbindLRP_' #'AF_score_distribution_manualLRP_' # 'AF_score_distribution_falseLRP_' #
    lrpair_score_dict = defaultdict(list)
    score_list = []
    for file_path in file_list:
        with open(file_path, 'r') as file:
            data = json.load(file)
    
        AF_score = data['iptm+ptm']['model_1_multimer_v3_pred_0']
        
        ligand = file_path.split('_')[1]
        receptor = file_path.split('_')[2]
        if ligand + '+' + receptor in lrp_list_to_run:
            lrpair_score_dict['pair'].append(ligand + '_to_' + receptor)
            lrpair_score_dict['AF scores'].append(AF_score)
            score_list.append(AF_score)
    
    plt.clf()
    sns.histplot(score_list, bins=30, kde=True, color='skyblue') #kde adds a kernel density estimate
    plt.title(plot_title) #predicted #selfbind
    plt.xlabel('AlphaFold Score')
    plt.ylabel('Frequency')
    plt.savefig(output_path + file_name  +str(len(score_list))+'LRP.jpg')
    
    data_list_pd = pd.DataFrame(lrpair_score_dict)
    data_list_pd.to_csv(output_path + file_name +str(len(score_list))+'LRP.csv', index=False)
    

