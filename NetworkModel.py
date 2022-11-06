###################################################################
###########   Preprocessing of data_networks_raw.csv file  ##################
###################################################################
## Objective: generate a graph based on SCI
# import library packages
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import pygraphviz as pgv
import csv
import math
import os.path
import Gurobi_Interface
from scipy import stats
import numpy as np
import Graph
SCI_Cutoff = 0.05 #for women cohert
def network_construct_analysis():
    x = range(1, 27)
    y_phi = []
    y_sci = []
    num_patients = []
    for i in range(1, 27):
        df = pd.read_csv(str(i) + '_comorbidity_network.csv')
        d1 = df[df['p-values'] < 0.01]
        edges_phi = d1.shape[0]
        y_phi.append(edges_phi)
        d2 = df[df['SCI value'] > 0.05]
        edges_SCI = d2.shape[0]
        y_sci.append(edges_SCI)
    plt.figure(1)
    #plt.plot(x, y_phi, label='Pearson')
    plt.plot(x, y_sci, label='SCI', marker='x', markersize=1.5)
    for a, b in zip(x, y_phi):
        if a % 5 == 0:
            plt.text(a, b, str(b))
    for a, b in zip(x, y_sci):
        if a % 5 == 0:
            plt.text(a, b, str(b))
    plt.legend()
    plt.grid(True)
    plt.show()

def network_model(cluster_list):
    ## Load raw encounter table
    network_raw = pd.read_csv("input_files/Data_networks_raw.csv")
    # Build comorbidity networks and draw them for each cluster
    for i1 in cluster_list:
        network = network_raw[(network_raw.LOS_bin.isin(cluster_list[i1]))]
        ## Only keep Diagnose_ID and Encounter.
        network = network[['Diagnose_ID', 'ENCOUNTER_ID']]
        ## Drop duplicate rows.
        network = network.drop_duplicates()
        results = network.values
        len_result = len(results)
        # Like graph adjacency list, we find all all related encounter ID for each diagnosed ID;
        # eg: data_list[0] = {10000,100001,....} means that encounter 10000,100001,..have been diagnosed with disease 0
        data_list = {}
        # add diagnosed ID and its encounter ID
        for i in range(len_result):
            key1 = results[i][0]
            value1 = results[i][1]
            if key1 in data_list:
                data_list[key1].append(value1)
            else:
                data_list[key1] = []
                data_list[key1].append(value1)
        #G is the comorbidity network;
        G = nx.Graph()
        # calculate the total number of patients who are diagnosed both disease key1 and key2, and create an undirected graph
        csvfile = str(i1) + '_comorbidity_network' + '.csv'
        with open(os.path.join("Output_files/comorbidity_network", csvfile), mode='w') as disease_file:
            writer = csv.writer(disease_file)
            #writer.writerow(['Dignose_ID', 'Dignose_ID', '#commonPaitents', 'Phi', 'p-values', 'SCI value']) #this is for p-values
            #writer.writerow(['Dignose_ID', 'Dignose_ID', '#commonPaitents']) #this is for only generating a comborbidity graph
            writer.writerow(['Dignose_ID1', 'Dignose_ID 2', '#commonPaitents', 'Total_ID1','Total_ID2','SCI']) #this is for only generating a comborbidity graph
            for key1 in data_list:
                for key2 in data_list:
                    if key1 >= key2:
                        continue
                    commonPatient = list(set(data_list[key1]).intersection(data_list[key2]))
                    C = len(commonPatient)
                    SCI = C / math.sqrt(len(data_list[key1]) * len(data_list[key2]))
                    # calculate the Pearson coeffient and p-value; we can comment them if we do not calculate p-value for SCI CUT OFF
                    #*******************************************************************#
                    '''
                    N = len_result
                    P1 = len(data_list[key1])
                    P2 = len(data_list[key2])
                    phi = 0.0
                    pval = 0.0
                    phi = (C * N - P1 * P2) / math.sqrt(P1 * P2 * (N - P1) * (N - P2))
                    #calculate p value based on T
                    if phi >= 1 or C < 2:
                        continue
                    T = phi * math.sqrt(C - 2) / math.sqrt(1 - phi ** 2)
                    # p-value based the t
                    pval = stats.t.sf(np.abs(T), C - 1) * 2 #two-sided pvalue
                    #writer.writerow([key1, key2, C, phi, pval, SCI])
                    #G.add_edge(key1, key2)
                    '''
                    # *******************************************************************#
                    if SCI >= SCI_Cutoff:
                        writer.writerow([key1, key2, C,len(data_list[key1]),len(data_list[key2]),SCI])
                        G.add_edge(key1, key2)

        data_list.clear()
        #print("writing--", str(i1), "-file")
        # write Graph G to GraphML format
        graph_name = str(i1) + '_comorbidity_network' + '.graphml'
        nx.write_graphml(G,os.path.join("Output_files/comorbidity_network", graph_name))
        #method 1: directly find ego graph of node 0: c-diff
        ego_graph = nx.ego_graph(G, 0)
        Graph.save_figure(ego_graph,i1)
        # write Graph G to GraphML format
        graph_name_GraphML = str(i1) + '_egoNetwork' + '.graphml'
        nx.write_graphml(ego_graph, os.path.join("input_files/egoNetwork", graph_name_GraphML))
        #method 2: sove the model using Gurobi
        #Gurobi_Interface.clubs_model_solve(G, i1)

if __name__ == "__main__":
    cluster_list = {}
    for i in range(1,27):
        cluster_list[i] = [2*i-1,2*i]
    network_model(cluster_list)
