
'''
Script to download IEDB data for B/T cell epitope and then compute RF (upper/lowerbound)
output:
    - ./data/IEDB_updated_data/iedb_epitopes_[%d%b%Y].csv (downloaded) [current day]
    - ./data/IEDB_updated_data/response_frequency_[%d%b%Y].csv (generated)
'''

from selenium.webdriver.common.by import By
from selenium import webdriver
import webbrowser
import time
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import glob
import os
import shutil
from datetime import datetime
import pandas as pd
import numpy as np
import requests
from selenium.webdriver.chrome.service import Service



#from: https://stackoverflow.com/questions/13059011/is-there-any-python-function-library-for-calculate-binomial-confidence-intervals
def binP(N, p, x1, x2):
    p = float(p)
    q = p/(1-p)
    k = 0.0
    v = 1.0
    s = 0.0
    tot = 0.0
    while(k<=N):
            tot += v
            if(k >= x1 and k <= x2):
                    s += v
            if(tot > 10**30):
                    s = s/10**30
                    tot = tot/10**30
                    v = v/10**30
            k += 1
            v = v*q*(N+1-k)/k
    return s/tot

def calcBin(vx, vN, vCL = 95):
    '''
    Calculate the exact confidence interval for a binomial proportion
    Usage:
    calcBin(13,100)
    (0.07107391357421874, 0.21204372406005856)
    calcBin(4,7)
    (0.18405151367187494, 0.9010086059570312)
    '''
    vx = float(vx)
    vN = float(vN)
    #Set the confidence bounds
    vTU = (100 - float(vCL))/2
    vTL = vTU
    vP = vx/vN
    if(vx==0):
            dl = 0.0
    else:
            v = vP/2
            vsL = 0
            vsH = vP
            p = vTL/100
            while((vsH-vsL) > 10**-5):
                    if(binP(vN, v, vx, vN) > p):
                            vsH = v
                            v = (vsL+v)/2
                    else:
                            vsL = v
                            v = (v+vsH)/2
            dl = v
    if(vx==vN):
            ul = 1.0
    else:
            v = (1+vP)/2
            vsL =vP
            vsH = 1
            p = vTU/100
            while((vsH-vsL) > 10**-5):
                    if(binP(vN, v, 0, vx) < p):
                            vsH = v
                            v = (vsL+v)/2
                    else:
                            vsL = v
                            v = (v+vsH)/2
            ul = v
    return (dl, ul)

def compute_RF_upperlowerbound(df):
    '''
    Compute Response Frequency lower/upperbound for each position (from IEDB table with B/T cell epitopes)
    check: https://help.iedb.org/hc/en-us/articles/114094147751
    '''
    #beginning/end of the protein
    start_position = df_all_epi['Mapped Start Position'].values
    end_position = df_all_epi['Mapped End Position'].values
    all_positions = list(set(list(start_position) + list(end_position)))
    #divide between linear and non-linear epitopes
    df_non_linear = df_all_epi.loc[df_all_epi['Sequence'].apply(lambda x:len(x.split(","))) > 1]
    df_linear = df_all_epi.loc[df_all_epi['Epitope ID'].isin(df_non_linear['Epitope ID'].values) == False]
    #dictionary pos, tested, reactive subject
    d_pos_tested = {}
    d_pos_resp = {}
    #for linear epitope
    for i in all_positions:
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
    lowerbound = []
    upperbound = []
    for pos in all_positions:
        N = d_pos_tested[pos]
        R = d_pos_resp[pos]
        rf = np.round(d_pos_resp[pos]/d_pos_tested[pos],2)
        # Wilson score interval for N>=50
        if N>50:
            lower95 = np.round((((R/N) + 1.96*1.96/(2*N) - 1.96 * np.sqrt(((R/N)*(1-(R/N))+1.96*1.96/(4*N))/N))/(1+1.96*1.96/N)),2)
            upper95 = np.round((((R/N) + 1.96*1.96/(2*N) + 1.96 * np.sqrt(((R/N)*(1-(R/N))+1.96*1.96/(4*N))/N))/(1+1.96*1.96/N)),2)
        # Binomial proportion confidence interval for N<50
        if N<50:
            lower95, upper95 = calcBin(R, N, vCL = 95)
            lower95 = np.round(lower95,2)
            upper95 = np.round(upper95,2)
        #print(pos, rf, lower95, upper95)
        lowerbound.append(lower95)
        upperbound.append(upper95)
    df = pd.DataFrame( {"positions":all_positions,  "lowerbound":lowerbound, "upperbound":upperbound})
    return df



####################################################################################################
### SETTINGS ###
####################################################################################################


list_protein = ['Spike', 'Nucleocapsid', 'Membrane', 'ORF1a', 'ORF1b', 'Envelope', 'ORF3a', 'ORF8', 'ORF6', 'ORF7a', 'ORF10']

#N.B. ORF1a and ORF1b are Replicase polyprotein 1ab (UniProt:P0DTD1) in IEDB
list_path_iedb = [
    "https://www.iedb.org/immunomebrowser.php?cookie_id=e601a1&source_organism=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FNCBITaxon_2697049&source_organism_name=SARS-CoV2&source_antigen=http%3A%2F%2Fwww.uniprot.org%2Funiprot%2FP0DTC2&source_antigen_name=Spike+glycoprotein",
    "https://www.iedb.org/immunomebrowser.php?cookie_id=e609a1&source_organism=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FNCBITaxon_2697049&source_organism_name=SARS-CoV2&source_antigen=http%3A%2F%2Fwww.uniprot.org%2Funiprot%2FP0DTC2&source_antigen_name=Nucleoprotein",
    "https://www.iedb.org/immunomebrowser.php?cookie_id=e609a1&source_organism=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FNCBITaxon_2697049&source_organism_name=SARS-CoV2&source_antigen=http%3A%2F%2Fwww.uniprot.org%2Funiprot%2FP0DTC2&source_antigen_name=Membrane+protein",
    "https://www.iedb.org/immunomebrowser.php?cookie_id=e609a1&source_organism=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FNCBITaxon_2697049&source_organism_name=SARS-CoV2&source_antigen=http%3A%2F%2Fwww.uniprot.org%2Funiprot%2FP0DTD1&source_antigen_name=Replicase+polyprotein+1ab",
    "https://www.iedb.org/immunomebrowser.php?cookie_id=e609a1&source_organism=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FNCBITaxon_2697049&source_organism_name=SARS-CoV2&source_antigen=http%3A%2F%2Fwww.uniprot.org%2Funiprot%2FP0DTD1&source_antigen_name=Replicase+polyprotein+1ab",
    "https://www.iedb.org/immunomebrowser.php?cookie_id=e609a1&source_organism=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FNCBITaxon_2697049&source_organism_name=SARS-CoV2&source_antigen=http%3A%2F%2Fwww.uniprot.org%2Funiprot%2FP0DTC2&source_antigen_name=Envelope",
    "https://www.iedb.org/immunomebrowser.php?cookie_id=e609a1&source_organism=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FNCBITaxon_2697049&source_organism_name=SARS-CoV2&source_antigen=http%3A%2F%2Fwww.uniprot.org%2Funiprot%2FP0DTC3&source_antigen_name=ORF3a+protein",
    'https://www.iedb.org/immunomebrowser.php?cookie_id=e609a1&source_organism=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FNCBITaxon_2697049&source_organism_name=SARS-CoV2&source_antigen=http%3A%2F%2Fwww.uniprot.org%2Funiprot%2FP0DTC3&source_antigen_name=ORF8+protein',
    'https://www.iedb.org/immunomebrowser.php?cookie_id=e609a1&source_organism=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FNCBITaxon_2697049&source_organism_name=SARS-CoV2&source_antigen=http%3A%2F%2Fwww.uniprot.org%2Funiprot%2FP0DTC3&source_antigen_name=ORF6+protein',
    'https://www.iedb.org/immunomebrowser.php?cookie_id=e609a1&source_organism=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FNCBITaxon_2697049&source_organism_name=SARS-CoV2&source_antigen=http%3A%2F%2Fwww.uniprot.org%2Funiprot%2FP0DTC3&source_antigen_name=ORF7a+protein',
    'https://www.iedb.org/immunomebrowser.php?cookie_id=e609a1&source_organism=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FNCBITaxon_2697049&source_organism_name=SARS-CoV2&source_antigen=http%3A%2F%2Fwww.uniprot.org%2Funiprot%2FP0DTC3&source_antigen_name=ORF10+protein']


#download correct chromedriver version from https://chromedriver.chromium.org/downloads and set path
path_chromedriver = '/home/giancarlo/Documents/programs/chromedriver'
print("** Path chromedriver: ", path_chromedriver)
path_download_folder = '/home/giancarlo/Downloads/*'

for protein, url in zip(list_protein, list_path_iedb):
    ############################################################
    # 1. Download IEDB epitope data and move to ./data/IEDB_updated_data
    ############################################################
    os.environ["webdriver.chrome.driver"] = path_chromedriver
    driver = webdriver.Chrome(path_chromedriver)

    #s = Service(path_chromedriver)
    #driver = webdriver.Chrome(service=s)

    print(protein, url)
    driver.get(url)
    #wait until page is loaded.. may take a while
    try:
        element = WebDriverWait(driver, 500).until(EC.presence_of_element_located((By.CLASS_NAME, "txt")))
        time.sleep(500)
        element.click()
        time.sleep(500)
    finally:
        driver.quit()
    #mv IEDB data from ~/Downloads (latest file) to ./data/IEDB_updated_data/PROTEIN
    list_of_files = glob.glob(path_download_folder)
    latest_file = max(list_of_files, key=os.path.getctime)
    print(latest_file)
    if latest_file.split("/")[-1].split("_")[0] != "immunomebrowser":
        print("*****ERROR: file doesn't start with `immunomebrowser`*****")
    #date last update IEDB
    req = requests.get(url)
    for word in req.text.split("\n"):
        if "site_data:" in word:
            site_data = word
    date_str = site_data.split(": ")[-1][1:-2]
    date_last_update = datetime.strptime(date_str, "%B %d, %Y")
    #adapt to your format
    format = "%d%b%Y"
    time_file = date_last_update.strftime(format)
    #print("Formatted DateTime:", time_file)
    name_out = "iedb_epitopes_{0}.csv".format(time_file)
    name_folder = './data/IEDB_updated_data/{0}'.format(protein)
    name_out = "iedb_epitopes_{0}.csv".format(time_file)
    path_iedb_epitope = './data/IEDB_updated_data/{0}/{1}'.format(protein, name_out)
    #make folder if it doesn't exist
    isExist = os.path.exists(name_folder)
    if not isExist:
        os.makedirs(name_folder)
        print("Created {0}".format(name_folder))
    #if file already exists -> skip it
    isExist = os.path.exists(path_iedb_epitope)
    if isExist:
        continue
    #otherwise move is to the right folder
    shutil.move(latest_file, "{0}/{1}".format(name_folder, name_out))
    ############################################################
    # 2. from epitope data get upper/lower rf (you could also download them directly from IEDB webserver)
    ############################################################
    df_all_epi = pd.read_csv(path_iedb_epitope)
    df = compute_RF_upperlowerbound(df_all_epi)
    name_out_rf = "response_frequency_{0}.csv".format(time_file)
    path_iedb_epitope_rf = './data/IEDB_updated_data/{0}/{1}'.format(protein,name_out_rf)
    #adapt to IEDB format
    f = open(path_iedb_epitope_rf, "w")
    print("\"position\",\"lowerbound\",\"upperbound\"", file = f)
    for idx in df.index:
        n = int(df.loc[idx]['positions'])
        low = df.loc[idx]['lowerbound']
        upp = df.loc[idx]['upperbound']
        print(str(n) + ",\""+ str(low) +  "\",\"" + str(upp) +  "\"", file = f)
    f.close()
