import pandas as pd
import numpy as np
import seaborn as sns
import pandas as pd
from sklearn import preprocessing
import matplotlib.pyplot as plt
from plotly.offline import init_notebook_mode, iplot
from plotly.graph_objs import *
import plotly.graph_objects as go
from sklearn import metrics


def read_dca_iedb_data(path_data = "./data/data_dca_proteome.csv", protein = 'Spike', domain = 'bCoV_S1_RBD'):
    df = pd.read_csv(path_data, sep = ',')
    if protein != None:
        df = df.loc[df['protein'] == protein]
    if domain != None:
        df = df.loc[df['domain'] == domain ]
    return df

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

def plot_dca_IEDB(df, score, list_pos = None):
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
    if list_pos != None:
        df_tmp = df.loc[df['position_protein'].isin(list_pos)]
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
