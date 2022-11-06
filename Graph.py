import networkx as nx
import pygraphviz as pgv
import statistics
import scipy.stats
import scipy.spatial.distance
import csv
import math
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import Gurobi_Interface
import silhouette_score
import NetworkModel
from netrd.distance import NetSimile
import numpy as np
import matplotlib.cm as cmap
from matplotlib.ticker import LinearLocator
import matplotlib.colors as colors

# read graphs from a folder
def read_graphs(input_directory):
    graphs = {}
    files_list = os.listdir(input_directory)
    n = 0  # num of files
    for file in files_list:
        file_name = input_directory + file
        if file_name == input_directory + "Icon\r":
            continue  # this is for skiping Icon file on Mac
        key = int(file.split('_')[0])
        graphs[key] = nx.read_graphml(file_name, node_type = int)
        n += 1
    return graphs

#This function is to aggreage cliques into a single node
def clique_signature(graph_list):
    #copy graph_list to G_list
    G_list = {}
    for k in graph_list:
        G_list[k] = graph_list[k].copy()
    n = len(graph_list) #num of files
    # clique signatures
    C = {} # used to store all cliques for each network
    node_size = {}
    for i in range(1, n+1):
        C[i] = []
    # while G_list is not empty
    while G_list:
        G_prime = {}
        for k1 in G_list:
            G_prime[k1] = G_list[k1].copy()
        D = []
        #find common nodes D
        isVisit = False
        for key1 in list(sorted(G_prime,key=lambda k: len(G_prime[k]), reverse=True)): #iterate dictionary by value length descending
            g = G_prime[key1]
            if not isVisit: #initialize the D using the graph with largest node size
                D = list(g.nodes).copy()
                isVisit = True
                continue
            if list(set(g.nodes).intersection(set(D))):
                D = list(set(g.nodes).intersection(set(D)))
            else:
                del G_prime[key1]
        # second while loop for finding a maximum clique in G_prime[D] in all networks
        D.sort()
        while D:
            D_prime = Gurobi_Interface.maximum_clique(G_prime,D) # find a subset D_prime such that it is maximum clique across all subgraph G_prime[D]
            for key2 in G_prime:
                g1 = G_prime[key2]
                C[key2].append(D_prime)
                g1.remove_nodes_from(D_prime)
                G_list[key2].remove_nodes_from(D_prime)
                if not g1.nodes:
                    del G_list[key2]
            D = list(set(D) - set(D_prime))
    sorted(C)
    return C
def ego_network_chart(graphs):
    #cliques = clique_signature(graphs)
    n = len(graphs)
    plt.figure()
    #graphs.keys must be sorted
    #for edges
    #plt.subplot(211)
    x_values = list(sorted(graphs.keys()))
    y_edges = [len(graphs[idx].edges) for idx in sorted(graphs.keys())]
    y_nodes = [len(graphs[idx].nodes) for idx in sorted(graphs.keys())]
    plt.plot(x_values,y_edges ,color='tab:blue',label='Edge',marker='x', linewidth=1, markersize=2)  # graph edges number
    d = {'time_window':x_values, 'nodes':y_nodes,'edges':y_edges}
    df = pd.DataFrame(data=d)
    df.to_csv('ego_node_edge.csv',index=True)
    #add data value on the plot
    '''
    for a, b in zip(x_values, y_edges):
        if a % 5 == 0:
            if b <= 100:
                plt.text(a-0.7, b + 9.5, str(b))
            elif 100 < b <= 1000:
                plt.text(a-0.5, b + 35, str(b))
            elif 1001 <b <= 2000:
                plt.text(a-0.5, b - 300, str(b))
            else:
                plt.text(a-0.5, b - 900, str(b))
    '''
    plt.yscale('log') # use log scale
    #plt.title('Number of edges of each ego network of C.diff')
    #plt.grid(True)
    plt.legend()
    #plt.subplot(212)  # for nodes
    plt.plot(x_values,y_nodes ,color='tab:orange',label='Node',marker='x', linewidth=1, markersize=2)  # graph nodes number;
    '''
    for a, b in zip(x_values, y_nodes):
        if a % 5 == 0:
            if b < 30:
                plt.text(a - 0.5, b + 2.5, str(b))
            else:
                plt.text(a - 0.5, b + 4.5, str(b))
    '''
    #plt.yscale('linear')
    #plt.title('Number of nodes and edges on ego network of C.diff')
    #plt.grid(True)
    plt.legend()
    plt.xlabel('Time window')
    #plt.plot(list(cliques.keys()), [len(cliques[idx]) for idx in graphs]) #clique number
    # Adjust the subplot layout, because the logit one may take more space
    # than usual, due to y-tick labels like "1 - 10^{-3}"
    plt.subplots_adjust(top=0.95, bottom=0.2, left=0.06, right=0.98, hspace=0.25,
                        wspace=0.45)
    plt.xticks(np.arange(2, 27, 2))  # x axis labels
    plt.savefig(os.path.join("output_files/figures", "egoNetwork.eps"), format="eps", dpi = 1000)
    plt.show()

def ego_network_chart_clustered(graphs,cliques):
    n = len(graphs)
    plt.figure()
    #graphs.keys must be sorted
    x_values = list(sorted(graphs.keys()))
    y_edges = [len(graphs[idx].edges) for idx in sorted(graphs.keys())]
    y_nodes = [len(graphs[idx].nodes) for idx in sorted(graphs.keys())]
    num_cliques = [len(cliques[i]) for i in sorted(cliques.keys())]
    plt.plot(x_values,y_edges ,color='tab:blue',label='Edge',marker='x',markersize=2, linewidth=0.5)  # graph edges number
    #add data value on the plot
    for a, b in zip(x_values, y_edges):
        if a % 1 == 0:
            if b <= 70:
                plt.text(a-0.05, b + 3.5, str(b))
            elif 70 < b <= 100:
                plt.text(a-0.08, b + 15, str(b))
            elif 100 <b <= 2000:
                plt.text(a-0.05, b - 100, str(b))
            else:
                plt.text(a-0.0, b - 900, str(b))
    plt.yscale('log') # use log scale
    #plt.title('Number of edges of each ego network of C.diff')
    plt.legend()
    plt.xticks(np.arange(1, 6, 1)) #x axis labels
    #plt.subplot(212)  # for nodes
    plt.plot(x_values,y_nodes ,color='tab:orange',label='Node',marker='x',markersize=2, linewidth=0.5)  # graph nodes number;
    for a, b in zip(x_values, y_nodes):
        if a % 1 == 0:
            if b < 13:
                plt.text(a - 0.05, b + 0.8, str(b))
            else:
                plt.text(a - 0.05, b + 2.5, str(b))
    #plt.yscale('linear')
    #plt.title('Number of nodes, edges, and cliques on each clustered network')
    plt.plot(x_values,num_cliques,color='tab:red',label='Clique',marker='x',markersize=2, linewidth=0.5)
    for a, b in zip(x_values, num_cliques):
        if a % 1 == 0:
            if b < 7:
                plt.text(a - 0.05, b + 0.5, str(b))
            else:
                plt.text(a - 0.09, b+1.1, str(b))
    #plt.grid(True)
    plt.legend()
    plt.xlabel('Clustered network')
    #plt.plot(list(cliques.keys()), [len(cliques[idx]) for idx in graphs]) #clique number
    # Adjust the subplot layout, because the logit one may take more space
    # than usual, due to y-tick labels like "1 - 10^{-3}"
    plt.subplots_adjust(top=0.95, bottom=0.2, left=0.06, right=0.98, hspace=0.25,
                        wspace=0.45)
    plt.savefig(os.path.join("output_files/figures", "clusterednetwork_stat.eps"), format="eps",dpi =1000)
    plt.show()

def similarity_figure(graphs,distances):
    n = len(graphs)
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    # Make data.
    x_values = list(sorted(graphs.keys()))
    x_len = len(x_values)
    y_values = list(sorted(graphs.keys()))
    y_len = len(y_values)
    xx, yy = np.meshgrid(x_values, y_values)
    z_values = []
    N = len(xx)
    for i1 in range(N):
        tmpList = []
        for i2 in range(N):
            i = xx[i1][i2]
            j = yy[i1][i2]
            tmpList.append(distances[i][j])
        z_values.append(tmpList)
    #z_values = [[distances[x_values[i1][i2]][y_values[i1][i2]] for i2 in range(x_len)]for i1 in range(x_len)]
    for i3 in range(N-1):
        print("The dissimilarity between graph %d and %d is: %f"%(i3+1,i3+2,z_values[i3][i3+1]))
    zz = np.array(z_values)
    #print(z_values)
    # Create an empty array of strings with the same shape as the meshgrid, and
    # populate it with two colors in a checkerboard pattern.
    #The following is to draw 3D figure
    '''
    colortuple = ('y', 'b')
    colors = np.empty(xx.shape, dtype=str)
    for y in range(y_len):
        for x in range(x_len):
            colors[x, y] = colortuple[(x + y) % len(colortuple)]
    # Plot the surface with face colors taken from the array we made.
    surf = ax.plot_surface(xx, yy, zz,rstride=1, facecolors=colors, linewidth=0)
    # Customize the z axis.
    ax.set_zlim(0, 1)
    ax.zaxis.set_major_locator(LinearLocator(6))
    '''
    #This is for a heatmap
    fig, ax = plt.subplots()
    im = ax.imshow(zz, cmap='Greens_r') #reverse the color Greens_r; otherwise Greens
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("Dissimilarity", rotation=-90, va="bottom")
    # We want to show all ticks...
    ax.set_xticks(np.arange(N))
    ax.set_yticks(np.arange(N))
    # ... and label them with the respective list entries
    ax.set_xticklabels(x_values,fontsize=8)
    ax.set_yticklabels(y_values,fontsize=8)
    # reverse the order of y label
    ax.invert_yaxis()

    # if needed, Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    '''
    for i in range(N):
        for j in range(N):
            text = ax.text(i, j, round(zz[i, j],2),
                           ha="center", va="center", color="w")
    '''
    #ax.set_title("Dissimilarity between Networks")
    plt.xlabel('Time window')
    plt.ylabel('Time window')
    fig.tight_layout()
    plt.subplots_adjust(top=0.95, bottom=0.16, left=-0.35, right=1, hspace=0.25,
                        wspace=0.25)
    plt.savefig(os.path.join('output_files/figures', 'Similarity_heatmap.eps'), format='eps', dpi=200)
    plt.show()

#This function is to generate new graphs by aggregating clique signatures into a single node
def aggregate_clique_into_node(graphs,cliques):
    # add all clique signatures of size >=2 into a set
    cliques_set = []
    for key in cliques:
        print("The graph", key, "has", len(cliques[key]), "cliques:", cliques[key], '\n')
        for val in cliques[key]:
            if val not in cliques_set and len(val) >= 2:
                cliques_set.append(val)
    cliques_set.sort()
    clique_nodes = {}  # for each clique, assign it to a node; if it includes 0 then call c0
    idx = 1
    for clique in cliques_set:
        if 0 in clique:
            clique_nodes['c0'] = clique
            continue
        else:
            clique_nodes['c' + str(idx)] = clique
        idx += 1
    new_graphs = {}  # store new graphs after agrregating cliques into a node
    for idx1 in graphs:
        g = graphs[idx1]
        new_graphs[idx1] = nx.Graph()
        tmp_cliques = cliques[idx1]
        node_size = len(tmp_cliques)
        for i1 in range(node_size):
            # add nodes
            nodes1 = tmp_cliques[i1]
            node_name1 = ''
            if len(nodes1) == 1:
                node_name1 = nodes1[0]
                new_graphs[idx1].add_node(node_name1, weight=1)
            else:
                node_name1 = list(clique_nodes.keys())[
                    list(clique_nodes.values()).index(nodes1)]  # find clique_nodes key by value
                new_graphs[idx1].add_node(node_name1, weight=len(nodes1))
            # add edges
            for j1 in range(i1 + 1, node_size):
                nodes2 = tmp_cliques[j1]
                node_name2 = ''
                if len(nodes2) == 1:
                    node_name2 = nodes2[0]
                else:
                    node_name2 = list(clique_nodes.keys())[
                        list(clique_nodes.values()).index(nodes2)]  # find clique_nodes key by value
                edge_weight = 0  # the total number of edges from nodes1 to nodes2 set
                for val1 in nodes1:
                    for val2 in nodes2:
                        if val2 in nx.neighbors(g, val1):
                            edge_weight += 1
                if edge_weight >= 1:
                    new_graphs[idx1].add_edge(node_name1, node_name2, weight=edge_weight)
    return new_graphs,clique_nodes #return a tuple with new graphs and cliques nodes

#save/show networkx as a figure file
def save_figure(G,index):
    # set node color
    val_map ={}
    if 0 in G.nodes:
        val_map = {'0': 'green'}  # set node 0 red
    else:
        val_map = {'c0': 'green'}  # set node c0 red
    values = [val_map.get(node, 'red') for node in G.nodes()]  # other nodes are green
    red_edges = []
    edge_colours = ['black' if not edge in red_edges else 'red'
                    for edge in G.edges()]
    black_edges = [edge for edge in G.edges() if edge not in red_edges]
    # Need to create a layout when doing
    # separate calls to draw nodes and edges
    pos = nx.spring_layout(G)
    plt.figure(index)
    nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'),
                           node_color=values, node_size=700)
    nx.draw_networkx_labels(G, pos)
    nx.draw_networkx_edges(G, pos, edgelist=red_edges, edge_color='r', arrows=True)
    nx.draw_networkx_edges(G, pos, edgelist=black_edges, arrows=False)
    figure_name = 'ego_network_LOS=' + str(index) + '.png'
    plt.savefig(os.path.join("output_files/egoNetwork", figure_name), format="PNG")
    # plt.show()

#generate csv files for Gephi visualization
#for general graph
def to_csv_egonetwork(graphs):
    for idx in sorted(graphs):
        n = len(graphs)
        g = graphs[idx]
        #for edge csv file:
        edge_csv_name = os.path.join("output_files/egoNetwork_26",'edges_ego_network_LOS=' + str(idx) + '.csv')
        f1 = open(edge_csv_name, 'w')
        edge_writer = csv.writer(f1)
        edge_writer.writerow(['Source', 'Target', 'Type', 'Id', 'Label', 'timeset'])
        count = 0
        for edge in g.edges():
            u = edge[0]
            v = edge[1]
            edge_writer.writerow([u,v, 'Undirected', count, '', ''])
            count += 1
        f1.close()

        #for nodes
        csv_name = os.path.join("output_files/egoNetwork_26", 'nodes_ego_network_LOS=' + str(idx) + '.csv')
        if idx == 1 and n >= 2:
            g_next = graphs[idx+1]
            red_nodes = []
            purple_nodes = []
            for node in g.nodes:
                if node in g_next.nodes:
                    purple_nodes.append(node)  # note that node is a string type
                else:
                    red_nodes.append(node)
            with open(csv_name, 'w') as f:
                writer = csv.writer(f)
                writer.writerow(['Id', 'Label', 'timeset', 'modularity_class'])
                for node in g.nodes():
                    if node == '0' or node == 0 or node == 'c0':
                        writer.writerow([node, node, '', 0])
                    elif node in purple_nodes:
                        writer.writerow([node, node, '', 3])
                    else:
                        writer.writerow([node, node, '', 4])
        elif 2 <= idx <= n-1:
            g_previous = graphs[idx-1]
            g_next = graphs[idx+1]
            blue_nodes = []
            aqua_nodes = []
            red_nodes = []
            purple_nodes = []
            for node in g.nodes:
                if (node in g_previous.nodes) and (node in g_next.nodes):
                    blue_nodes.append(node)  # note that node is a string type
                if (node in g_previous.nodes) and (node not in g_next.nodes):
                    aqua_nodes.append(node)
                if (node not in g_previous.nodes) and (node in g_next.nodes):
                    purple_nodes.append(node)
                if (node not in g_previous.nodes) and (node not in g_next.nodes):
                    red_nodes.append(node)
            with open(csv_name, 'w') as f:
                writer = csv.writer(f)
                writer.writerow(['Id', 'Label', 'timeset', 'modularity_class'])
                for node in g.nodes():
                    if node == '0' or node == 0 or node == 'c0':
                        writer.writerow([node, node, '', 0])
                    elif node in blue_nodes:
                        writer.writerow([node, node, '', 1])
                    elif node in aqua_nodes:
                        writer.writerow([node, node, '', 2])
                    elif node in purple_nodes:
                        writer.writerow([node, node, '', 3])
                    else:
                        writer.writerow([node, node, '', 4])
        elif idx == n:
            g_previous = graphs[idx - 1]
            aqua_nodes = []
            red_nodes = []
            for node in g.nodes:
                if node in g_previous.nodes:
                    aqua_nodes.append(node)
                else:
                    red_nodes.append(node)
            with open(csv_name, 'w') as f:
                writer = csv.writer(f)
                writer.writerow(['Id', 'Label', 'timeset', 'modularity_class'])
                for node in g.nodes:
                    if node == '0' or node == 0 or node == 'c0':
                        writer.writerow([node, node, '', 0])
                    elif node in aqua_nodes:
                        writer.writerow([node, node, '', 2])
                    else:
                        writer.writerow([node, node, '', 4])
#for graphs with cliques aggregation
def to_csv_aggreaged_network(graphs,cliques):
    # find the maximum weight in all graphs
    max_weight = 0
    for idx in graphs:
        g = graphs[idx]
        for edge in g.edges():
            u = edge[0]
            v = edge[1]
            if g[u][v]['weight'] > max_weight:
                max_weight = g[u][v]['weight']
    node_size = 50
    for idx in sorted(graphs):
        n = len(graphs)
        g = graphs[idx]
        #for edge csv file:
        edge_csv_name = os.path.join("output_files/clustered_network",'edges_ego_network_LOS=' + str(idx) + '.csv')
        f1 = open(edge_csv_name, 'w')
        edge_writer = csv.writer(f1)
        edge_writer.writerow(['Source', 'Target', 'Type', 'Id', 'Label', 'timeset', 'Weight'])
        count = 0
        for edge in g.edges():
            u = edge[0]
            v = edge[1]
            edge_writer.writerow([u,v, 'Undirected', count, '', '', node_size * g[u][v]['weight']/max_weight]) #edge weight is relative， = min node_size x weight/max_weight
            count += 1
        f1.close()

        #for nodes
        csv_name = os.path.join("output_files/clustered_network", 'nodes_ego_network_LOS=' + str(idx) + '.csv')
        if idx == 1 and n >= 2:
            g_next = graphs[idx+1]
            red_nodes = []
            purple_nodes = []
            for node in g.nodes:
                if node in g_next.nodes:
                    purple_nodes.append(node)  # note that node is a string type
                else:
                    red_nodes.append(node)
            with open(csv_name, 'w') as f:
                writer = csv.writer(f)
                writer.writerow(['Id', 'Label', 'timeset', 'modularity_class','weight'])
                for node in g.nodes():
                    if node == '0' or node == 0  or node == 'c0':
                        writer.writerow([node, node, '', 0, len(cliques[node]) if node in cliques else 1])
                    elif node in purple_nodes:
                        writer.writerow([node, node, '', 3, len(cliques[node]) if node in cliques else 1])
                    else:
                        writer.writerow([node, node, '', 4, len(cliques[node]) if node in cliques else 1])
        elif 2 <= idx <= n-1:
            g_previous = graphs[idx-1]
            g_next = graphs[idx+1]
            blue_nodes = []
            aqua_nodes = []
            red_nodes = []
            purple_nodes = []
            for node in g.nodes:
                if (node in g_previous.nodes) and (node in g_next.nodes):
                    blue_nodes.append(node)  # note that node is a string type
                if (node in g_previous.nodes) and (node not in g_next.nodes):
                    aqua_nodes.append(node)
                if (node not in g_previous.nodes) and (node in g_next.nodes):
                    purple_nodes.append(node)
                if (node not in g_previous.nodes) and (node not in g_next.nodes):
                    red_nodes.append(node)
            with open(csv_name, 'w') as f:
                writer = csv.writer(f)
                writer.writerow(['Id', 'Label', 'timeset', 'modularity_class','Weight'])
                for node in g.nodes():
                    if node == '0' or node == 0 or node == 'c0':
                        writer.writerow([node, node, '', 0, len(cliques[node]) if node in cliques else 1])
                    elif node in blue_nodes:
                        writer.writerow([node, node, '', 1, len(cliques[node]) if node in cliques else 1])
                    elif node in aqua_nodes:
                        writer.writerow([node, node, '', 2, len(cliques[node]) if node in cliques else 1])
                    elif node in purple_nodes:
                        writer.writerow([node, node, '', 3, len(cliques[node]) if node in cliques else 1])
                    else:
                        writer.writerow([node, node, '', 4, len(cliques[node]) if node in cliques else 1])
        elif idx == n:
            g_previous = graphs[idx - 1]
            aqua_nodes = []
            red_nodes = []
            for node in g.nodes:
                if node in g_previous.nodes:
                    aqua_nodes.append(node)
                else:
                    red_nodes.append(node)
            with open(csv_name, 'w') as f:
                writer = csv.writer(f)
                writer.writerow(['Id', 'Label', 'timeset', 'modularity_class','Weight'])
                for node in g.nodes:
                    if node == '0' or node == 0 or node == 'c0':
                        writer.writerow([node, node, '', 0, len(cliques[node]) if node in cliques else 1])
                    elif node in aqua_nodes:
                        writer.writerow([node, node, '', 2, len(cliques[node]) if node in cliques else 1])
                    else:
                        writer.writerow([node, node, '', 4, len(cliques[node]) if node in cliques else 1])

# This function is to generate ego network csv files for Gephi visulization
def ego_network_output_visualization(input_directory):
    #read graphs from a folder
    graphs = {}
    graphs = read_graphs(input_directory)
    n = len(graphs) #num of files
    cliques = {}
    cliques = clique_signature(graphs) #get all clique signatures from a collection of networks
    new_graphs,clique_node = aggregate_clique_into_node(graphs,cliques)
    for idx in new_graphs:
        save_figure(new_graphs[idx],idx)
    to_csv_aggreaged_network(new_graphs,clique_node)

def show_progress(state,idx, total):
    percentage = float(idx) / float(total) * 100
    sys.stdout.write('\r%s %d%% Complete' % (state, percentage))
    sys.stdout.flush()


def get_network_distance_KNC(graphs):
    n = len(graphs)
    distances = {i: {j: 0 for j in range(1,n+1)} for i in range(1,n+1)} #initialize the distances with 0
    for idx1 in range(1,n):
        show_progress("KNC Matrix Distance calculation",idx1,n)
        g1 = graphs[idx1]
        for idx2 in range(idx1+1, n+1):
            g2 = graphs[idx2]
            common_nodes = list(set(g1.nodes).intersection(set(g2.nodes)))
            common_nodes.sort()
            g1_subgraph = g1.subgraph(common_nodes)
            g2_subgraph = g2.subgraph(common_nodes)
            common_e1 = len(g1_subgraph.edges)
            common_e2 = len(g2_subgraph.edges)
            weight = (common_e1 + common_e2) / (len(g1.edges) + len(g2.edges)) #calculate weights in terms of edges overlap
            #Caution: netrd.distance.netsimile library only works for graphs with continous node numbers; for sparse graphs with nodes like 1,5,10...
            # we need to change one place in the source file netsimile.py: line 106
            # change egonets = [nx.ego_graph(G, n) for n in node_list] to egonets = {n:nx.ego_graph(G, n) for n in node_list}
            if idx1 == 2 and idx2 ==3:
                common_nodes1 = list(set(g1.nodes).intersection(set(g2.nodes)))

            matrix_distance = 0
            CN = len(common_nodes)
            for k1 in range(CN):
                node1 = common_nodes[k1]
                for k2 in range(k1+1,CN):
                    node2 = common_nodes[k2]
                    if (node1 in g1_subgraph.neighbors(node2) and node1 not in g2_subgraph.neighbors(node2)) or (node1 not in g1_subgraph.neighbors(node2) and node1 in g2_subgraph.neighbors(node2)):
                        matrix_distance += 1

            normalize_matrix_distance = math.sqrt(matrix_distance/(CN*(CN-1)/2)) #normalized distance
            distances[idx1][idx2] = normalize_matrix_distance
            distances[idx2][idx1] = normalize_matrix_distance
            #distances[idx1][idx2] = 1 - weight * (1 - normalize_matrix_distance)
            #distances[idx2][idx1] = 1 - weight * (1 - normalize_matrix_distance)
            if idx2 == idx1 +1:
                print('The size of common neighbors is:-',CN, 'and the normalize_matrix_distance is:',normalize_matrix_distance)
    for i1 in range(1,n):
        print('The modified distance between graph %d and %d is: %f'%(i1,i1+1,distances[i1][i1+1]))
    #convert distance to pandas data frame and then write to a csv file
    df = pd.DataFrame(data=distances)
    df.to_csv('input_files/distances_matrix_temp'+str(n) +'.csv',index=True)
    return distances
def get_network_Netdistance(graphs):
    n = len(graphs)
    distances = {i: {j: 0 for j in range(1,n+1)} for i in range(1,n+1)} #initialize the distances with 0
    for idx1 in range(1,n):
        show_progress("Modified NetSimile Distance calculation",idx1,n)
        g1 = graphs[idx1]
        for idx2 in range(idx1+1, n+1):
            g2 = graphs[idx2]
            common_nodes = list(set(g1.nodes).intersection(set(g2.nodes)))
            common_nodes.sort()
            g1_subgraph = g1.subgraph(common_nodes)
            g2_subgraph = g2.subgraph(common_nodes)
            common_e1 = len(g1_subgraph.edges)
            common_e2 = len(g2_subgraph.edges)
            weight = (common_e1 + common_e2) / (len(g1.edges) + len(g2.edges)) #calculate weights in terms of edges overlap
            #Caution: netrd.distance.netsimile library only works for graphs with continous node numbers; for sparse graphs with nodes like 1,5,10...
            # we need to change one place in the source file netsimile.py: line 106
            # change egonets = [nx.ego_graph(G, n) for n in node_list] to egonets = {n:nx.ego_graph(G, n) for n in node_list}
            netsimile = NetSimile()
            if idx1 == 2 and idx2 ==3:
                common_nodes1 = list(set(g1.nodes).intersection(set(g2.nodes)))
            netsimile_distance = netsimile.dist(g1_subgraph,g2_subgraph)
            distances[idx1][idx2] = 1 - weight * (1 - netsimile_distance/35)
            distances[idx2][idx1] = 1 - weight * (1 - netsimile_distance/35)
    for i1 in range(1,n):
        print('The modified distance between graph %d and %d is: %f'%(i1,i1+1,distances[i1][i1+1]))
    #convert distance to pandas data frame and then write to a csv file
    df = pd.DataFrame(data=distances)
    df.to_csv('input_files/distances_'+str(n) +'.csv',index=True)
    return distances
# find the cluster p value

#Currently this function is not used for our paper
def get_network_distance(graphs):
    n = len(graphs)
    distances = {i: {j: 0 for j in range(1,n+1)} for i in range(1,n+1)} #initialize the distances with 0
    for idx1 in range(1,n):
        show_progress("Modified NetSimile Distance calculation",idx1,n)
        g1 = graphs[idx1]
        for idx2 in range(idx1+1, n+1):
            g2 = graphs[idx2]
            common_nodes = list(set(g1.nodes).intersection(set(g2.nodes)))
            common_nodes.sort()
            common_e1 = len(g1.subgraph(common_nodes).edges)
            common_e2 = len(g2.subgraph(common_nodes).edges)
            weight = (common_e1 + common_e2) / (len(g1.edges) + len(g2.edges)) #calculate weights in terms of edges overlap
            #Caution: netrd.distance.netsimile library only works for graphs with continous node numbers; for sparse graphs with nodes like 1,5,10...
            # we need to change one place in the source file netsimile.py: line 106
            # change egonets = [nx.ego_graph(G, n) for n in node_list] to egonets = {n:nx.ego_graph(G, n) for n in node_list}
            netsimile = NetSimile()
            netsimile_distance = netsimile.dist(g1,g2)
            distances[idx1][idx2] = 1 - weight * (1 - netsimile_distance/35)
            distances[idx2][idx1] = 1 - weight * (1 - netsimile_distance/35)
    for i1 in range(1,n):
        print('The modified distance between graph %d and %d is: %f'%(i1,i1+1,distances[i1][i1+1]))
    #convert distance to pandas data frame and then write to a csv file
    df = pd.DataFrame(data=distances)
    df.to_csv('input_files/distances_'+str(n) +'.csv',index=True)
    return distances
# find the cluster p value
def read_distances(csv_filename):
    df = pd.read_csv(csv_filename, index_col=0)
    df.columns = df.columns.astype(int)
    distances = df.to_dict()
    return distances

def find_cluster_value_by_Silhouette(distances,tau):
    n = len(distances)
    x_axis_values = []
    y_axis_values = []
    for i in range(1,n+1):
        solution_dict = Gurobi_Interface.network_cluster(n, distances, i, tau)
        if not solution_dict:
            continue #if the solution is feasible, skip.
        x_axis_values.append(i)
        y_axis_values.append(silhouette_score.silhouette_score(distances,solution_dict,tau))
    plt.figure()
    plt.plot(x_axis_values, y_axis_values, marker='x', linewidth=1, markersize=3)
    plt.xticks(np.arange(2, 27, 2))
    plt.xlabel('Cluster value p')
    plt.ylabel('Silhouette index')
    #plt.title('Silhouette index chart for different values of p')
    #plt.grid(True)
    plt.subplots_adjust(top=0.95, bottom=0.16, left=0.12, right=0.97, hspace=0.25,
                        wspace=0.25)
    plt.savefig(os.path.join("output_files/figures", "Silhouette_method.eps"), format='eps', dpi=1000)
    plt.show()

def find_cluster_value_by_Elbow(distances,tau):
    n = len(distances)
    x_axis_values =[]
    y_axis_values =[]
    for i in range(1,n+1):
        objvalue = Gurobi_Interface.network_cluster(n, distances, i, tau)
        if objvalue == -1:
            continue
        x_axis_values.append(i)
        y_axis_values.append(objvalue)
    plt.plot(x_axis_values,y_axis_values, marker = 'x', linewidth=2, markersize=12)
    plt.savefig(os.path.join("output_files/figures", "Elbow_method.png"), format="PNG")
    plt.show()

if __name__ == "__main__":
    directory = "input_files/egoNetwork_26/"
    graphs = read_graphs(directory)
    n = len(graphs)
    distances = get_network_Netdistance(graphs)
    distances = read_distances("input_files/distances_" + str(n) + ".csv") #this is only 26 distances
    ego_network_chart(graphs)
    similarity_figure(graphs,distances)
    #to_csv_egonetwork(graphs) #generate 26 ego network csv files for Gephi visualization
    tau = 0.5
    find_cluster_value_by_Silhouette(distances, tau)
    #find_cluster_value_by_Elbow(distances, tau)
    #p value based on the figure: 3 is the best one
