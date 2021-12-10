import pandas as pd
import os
import datetime
import numpy as np
import seaborn as sns
import pandas as pd
from sklearn import preprocessing
import matplotlib.pyplot as plt
from plotly.offline import init_notebook_mode, iplot
from plotly.graph_objs import *
import plotly.graph_objects as go
from sklearn import metrics


def load_data_dca(path_data = "./data/data_dca_proteome.csv", protein = None, domain = None):
    df = pd.read_csv(path_data, sep = ',')
    if protein != None:
        df = df.loc[df['protein'] == protein]
    if domain != None:
        df = df.loc[df['domain'] == domain ]
    return df

def load_VOC(path_dir_VOC = "./data/VOC", rm_deletion = True, verbose = False):
    d_pange_who = {
        'B.1.1.7':"Alpha",
        'B.1.351': "Beta",
        'B.1.617.2': "Delta",
        'P.1':"Gamma",
        'BA.1':"Omicron"
        }
    #adapt cov-lineage protein names
    d_prot_name = {
        'E':'Envelope',
        'M':'Membrane',
        'N':'Nucleocapsid',
        'S':'Spike',
        'orf1a':'ORF1A',
        'orf1b':'ORF1B',
        'orf3a':'ORF3A',
        'Orf8':'ORF8',
        'orf8':'ORF8'
        }
    all_pango = []
    all_who = []
    all_mutation_type = []
    all_protein = []
    all_SNP = []
    all_pos_mutated = []
    for filename in os.listdir(path_dir_VOC):
        f = open(os.path.join(path_dir_VOC, filename))
        for line in f.readlines():
            for mut in line.split():
                try:
                    list_mut = mut.split(":")
                    if len(list_mut) == 2:
                        mut_type = list_mut[0]
                        if mut_type == 'nuc':
                            if verbose:
                                print('skipping {0}'.format(list_mut))
                            continue
                        mut_type = "aa"
                        mut_prot = list_mut[0]
                        mut_SNP = list_mut[1]
                    if len(list_mut) == 3:
                        mut_type = list_mut[0]
                        mut_prot = list_mut[1]
                        mut_SNP = list_mut[2]
                    if verbose:
                        print(filename, list_mut)
                    #append
                    pos_mutated = int(mut_SNP[1:-1])
                    all_pos_mutated.append(pos_mutated)
                    all_SNP.append(mut_SNP)
                    all_pango.append(filename)
                    all_who.append(d_pange_who[filename])
                    all_mutation_type.append(mut_type)
                    if mut_prot in d_prot_name.keys():
                        all_protein.append(d_prot_name[mut_prot])
                    else:
                        all_protein.append(mut_prot)
                except:
                    if verbose:
                        print("Error with {0}, pos={1}".format(filename, list_mut))
                    continue
    df_VOC = pd.DataFrame({'pango':all_pango, 'who':all_who,'mutation_type':all_mutation_type, 'protein':all_protein, 'mutation':all_SNP, 'position':all_pos_mutated})
    if rm_deletion:
        df_VOC_del = df_VOC.loc[df_VOC['mutation_type'] == 'del']
        df_VOC = df_VOC.loc[set(df_VOC.index)- set(df_VOC_del.index)]
    return df_VOC


def plot_roc(low_mut, high_mut, observed_mutability, list_score, df, add_obs_mut = False):
    #plot ROC
    df['tp'] = -1
    df.loc[df[observed_mutability].isin(low_mut), 'tp'] = 0
    df.loc[df[observed_mutability].isin(high_mut), 'tp'] = 1
    df = df.loc[df['tp'] != -1]
    for score in list_score:
        fpr, tpr, thresholds = metrics.roc_curve(df['tp'].values , df[score].values)
        AUC = metrics.auc(fpr, tpr)
        num_tp = np.sum(df['tp'])
        total = len(df['tp'])
        if add_obs_mut:
            lab = "{0}_AUC:{1} ({2})".format(score,str(round(AUC,2)), observed_mutability)
        else:
            lab = "{0}_AUC:{1}".format(score,str(round(AUC,2)))
        plt.plot(fpr,tpr, label = lab)
    plt.legend()

    return 0

def plot_dca_IEDB(df, score, list_pos = []):
    confidence_interval =  (df['IEDB_upperbound'].values -  df['IEDB_lowerbound'].values)
    size_scatter =  (1/confidence_interval)
    df['size_scatter'] = size_scatter
    fig = go.FigureWidget()
    #all data
    trace1 = fig.add_scattergl(
        x=df[score],
        y=df['IEDB_response_frequency'],
        text="wt: "+ df['aa_Wuhan-Hu-1'] + df['position_protein'].map(str),
        textposition='top right',
        textfont=dict(color='#E58606'),
        #mode='markers+text',
        mode='markers',
        marker=dict(color='#5D69B1', size=df['size_scatter']),
        #marker=dict(color=df['lineage'].map(str).map(len), size=df['size_scatter']),
        name='')
    if len(list_pos) > 0:
        df_tmp = df.loc[df['position_protein'].isin(list_pos)]
        #prot = df_tmp['protein'].values[0]
        #dom = df_tmp['domain'].values[0]
        list_mutated_prot_dom = df_tmp['position_protein'].values
        print("list mutated positions = {0}".format(list_mutated_prot_dom))
        fig.add_scattergl(
                x=df_tmp[score],
                y=df_tmp['IEDB_response_frequency'],
                text="wt: "+ df_tmp['aa_Wuhan-Hu-1'] + df_tmp['position_protein'].map(str),
                textposition='top right',
                textfont=dict(color='#E58606'),
                mode='markers',
                marker=dict(color='red', size=df_tmp['size_scatter']),
                showlegend=False
        )
    fig.layout = dict(
        plot_bgcolor="#FFF",
        legend=dict(
            # Adjust click behavior
            itemclick="toggleothers",
            itemdoubleclick="toggle",
        ),
        margin=dict(t=20, l=20, r=20, b=20),
        xaxis=dict(title=score, linecolor='#BCCCDC', showgrid=True, mirror=True),
        yaxis=dict(title='IEDB - Response Frequency', linecolor='#BCCCDC', showgrid=True, mirror=True),
        )
    iplot(fig)
    return 0




def plot_dca_IEDB_BTcell(df, score, list_pos = None, cell_type = "B_cell"):
    confidence_interval =  (df['upperbound_'+cell_type].values -  df['lowerbound_'+cell_type].values)
    size_scatter =  (1/confidence_interval)
    df['size_scatter'] = size_scatter
    fig = go.FigureWidget()
    #all data
    trace1 = fig.add_scattergl(
        x=df[score],
        y=df['rf_'+cell_type],
        text="wt: "+ df['aa_Wuhan-Hu-1'] + df['position_protein'].map(str),
        textposition='top right',
        textfont=dict(color='#E58606'),
        #mode='markers+text',
        mode='markers',
        marker=dict(color='#5D69B1', size=df['size_scatter']),
        #marker=dict(color=df['lineage'].map(str).map(len), size=df['size_scatter']),
        name='')
    if list_pos != None:
        df_tmp = df.loc[df['position_protein'].isin(list_pos)]
        fig.add_scattergl(
                x=df_tmp[score],
                y=df_tmp['rf_'+cell_type],
                text="wt: "+ df_tmp['aa_Wuhan-Hu-1'] + df_tmp['position_protein'].map(str),
                textposition='top right',
                textfont=dict(color='#E58606'),
                mode='markers',
                marker=dict(color='red', size=df_tmp['size_scatter']),
                showlegend=False
        )
    fig.layout = dict(
        plot_bgcolor="#FFF",
        legend=dict(
            # Adjust click behavior
            itemclick="toggleothers",
            itemdoubleclick="toggle",
        ),
        margin=dict(t=20, l=20, r=20, b=20),
        xaxis=dict(title=score, linecolor='#BCCCDC', showgrid=True, mirror=True),
        yaxis=dict(title='IEDB - Response Frequency', linecolor='#BCCCDC', showgrid=True, mirror=True),
        )
    iplot(fig)
    return 0


################################
# For updated IEDB data

def get_IEDB_versions(path_IEDB_epitope_data = "./data/IEDB_updated_data"):
    version_list = []
    for filename in os.listdir(path_IEDB_epitope_data):
        if filename[-3:] != "csv":
            continue
        version_date = filename.split("_")[-1][:-4]
        if version_date not in version_list:
            version_list.append(version_date)
    #sort by date
    version_list.sort(key=lambda date: datetime.datetime.strptime(date, "%d%b%Y"))
    print("IEDB available versions:", version_list)
    return version_list




def compute_RF(df, path_epitope):
    """ Compute mean Response Frequency for each position (from IEDB table with B/T cell epitopes)"""
    df_all_epi = pd.read_csv(path_epitope)
    df_non_linear = df_all_epi.loc[df_all_epi['Sequence'].apply(lambda x:len(x.split(","))) > 1]
    df_linear = df_all_epi.loc[df_all_epi['Epitope ID'].isin(df_non_linear['Epitope ID'].values) == False]
    #dictionary pos, tested, reactive subject
    d_pos_tested = {}
    d_pos_resp = {}
    #for linear epitope
    for i in df['position_protein'].values:
        beg = df_linear['Mapped Start Position']
        end = df_linear['Mapped End Position']
        df_tmp = df_linear.loc[(beg <= i) &(end >= i)]
        d_pos_tested[i] = (np.sum(df_tmp['Subjects Tested']))
        d_pos_resp[i] = (np.sum(df_tmp['Subjects Responded']))
    #for Non linear epitope
    for idx in df_non_linear.index:
        sub_tested = df_non_linear.loc[idx]['Subjects Tested']
        sub_resp = df_non_linear.loc[idx]['Subjects Responded']
        pos_epitope = df_non_linear.loc[idx]['Sequence'].split(",")
        for pos in pos_epitope:
            pos = pos.replace(' ',"")
            pos = int(pos[1:])
            if pos in d_pos_resp.keys():
                d_pos_resp[pos] += sub_resp
            else:
                d_pos_resp[pos] = sub_resp
            if pos in d_pos_tested.keys():
                d_pos_tested[pos] += sub_tested
            else:
                d_pos_tested[pos] = sub_resp
    #add subject reponded, test, and RF to df
    df['subj_tested'] = [d_pos_tested[pos] for pos in df['position_protein'].values]
    df['subj_responded'] = [d_pos_resp[pos] for pos in df['position_protein'].values]
    df['IEDB_response_frequency'] = df['subj_responded']/ df['subj_tested']
    return df


def get_updated_IEDB(df, version, path_IEDB_epitope_data = "./data/IEDB_updated_data"):
        print("Selecting *** IEDB {0} version ***".format(version))
        path_rf_lower_upper = os.path.join(path_IEDB_epitope_data, "response_frequency_"+str(version)+".csv")
        path_epitope = os.path.join(path_IEDB_epitope_data, "iedb_epitopes_"+str(version)+".csv")
        response_freq = pd.read_csv(path_rf_lower_upper)
        response_freq = response_freq.rename(columns = {'position': 'position_protein','upperbound':'IEDB_upperbound', 'lowerbound':'IEDB_lowerbound'})
        df = pd.merge(left = df, right = response_freq, on = 'position_protein')
        #compute also mean_rf (from epitope data)
        compute_RF(df, path_epitope)
        return df

