#Original code was from Alexander Veremyev
#Modified by Yajun Lu for this project

import igraph as ig
from IPython.core.display import display, SVG
from random import randint
import networkx as nx
import os
import Graph
import colorsys
# graph g should be iGraph instance
def nxGraph_to_iGraph(nxGraph):
    # use iGraph to draw communities - create iG from G
    iG = ig.Graph()
    iG.add_vertices(list(nxGraph.nodes()))
    # since during transformuation nodes in iG are labeled from 0 to n-1, edge labels need to be converted as well
    # suppose nodes ID are not continous in networkX, we need a map between networkx and iGraph
    edges_ig = [(i, j) for i in iG.vs for j in iG.vs if (i['name'], j['name']) in nxGraph.edges() and i < j]
    iG.add_edges(edges_ig)
    return iG

def get_membership_index(node,mylist):
    #eg: mylist = [[1,2],[3,4,5],[9]], node 4 index is 1 (list index starting from 0)
    for idx in range(len(mylist)):
        list1 = mylist[idx]
        if node in list1:
            return idx

def get_membership(iGraph, cliques_set):
    # dict_mem = {1:0, 2:0, 3:0, 4:1, 5:1, 6:1,7:2,8:2}
    dict_mem = {}
    for vertex in iGraph.vs():
        dict_mem[vertex['name']] = get_membership_index(vertex['name'], cliques_set)
    membership = list()
    for vertex in iGraph.vs():
        membership.append(dict_mem[vertex['name']])
    return membership

def _plot(g, filename1,membership_colors,membership=None):
    #g: igraph type; membership: load membership - dictionary (key - node, value - membership)
    if membership is not None:
        gcopy = g.copy()
        edges = []
        edges_colors = []
        edges_widths = []
        for edge in g.es():
            if membership[edge.tuple[0]] != membership[edge.tuple[1]]:
                edges.append(edge)
                edges_colors.append("gray")
                edges_widths.append(1.5)
            else:
                edges_colors.append("black")
                edges_widths.append(3)
        gcopy.delete_edges(edges)
        layout = gcopy.layout("fr")
        g.es["color"] = edges_colors
        v = ig.VertexClustering(g, membership)
    else:
        layout = g.layout("fr")
        g.es["color"] = "gray"
    visual_style = {}
    visual_style["vertex_label_dist"] = 0
    visual_style["edge_color"] = g.es["color"]
    visual_style["vertex_size"] = 170
    visual_style['vertex_label_size'] = 80
    visual_style["edge_width"]= edges_widths
    visual_style["layout"] = layout
    if filename1 == "output_files/egoNetwork_5/"+'cluster_network_5.pdf' or filename1 == "output_files/egoNetwork_5/"+'cluster_network_4.pdf':
        visual_style["bbox"] = (6500, 4500)
    else:
        visual_style["bbox"] = (3000, 3000)
    visual_style["margin"] = 350
    #visual_style["edge_label"] = g.es["weight"]
    for vertex in g.vs():
        #vertex["label"] = vertex.index #this show new id in igraph
        vertex["label"] = vertex['name']
    if membership is not None:
        colors = []
        #if colors need to be defined randomly (for many communities)
        for i in range(0, max(membership)+1):
            colors.append('%06X' % randint(0, 0xFFFFFF))
        for vertex in g.vs():
            #vertex["color"] = str('#') + colors[membership[vertex.index]] #random color
            vertex["color"] = membership_colors[membership[vertex.index]] # assign a fix color for each membership
            vertex["shape"] = "circle"
            '''
            if membership[vertex.index]==0:
                vertex["shape"]="circle"
                vertex["color"]="red"
            #if membership[vertex.index]==1:
            else:
                vertex["shape"]="circle"
                vertex["color"] = "pink"
            '''
        visual_style["vertex_color"] = g.vs["color"]
        visual_style["vertex_shape"] = g.vs["shape"]
    visual_style["mark_groups"] = True # for each group, show shadow or not
    #ig.plot(v, **visual_style).save(filename1) #save as png file
    ig.plot(v,filename1, **visual_style) #save as a pdf file


if __name__ == "__main__":
    input_directory = "input_files/egoNetwork_5/"
    graphs = {}
    graphs = Graph.read_graphs(input_directory)
    n = len(graphs) #num of files
    cliques = {}
    cliques = Graph.clique_signature(graphs) #get all clique signatures from a collection of networks
    cliques_set = []
    for key in cliques:
        print("The graph", key, "has", len(cliques[key]), "cliques:", cliques[key], '\n')
        for val in cliques[key]:
            if val not in cliques_set:
                cliques_set.append(val)
    cliques_set.sort()
    N = len(cliques_set)
    membership_colors={}
    colors = ['#C4D67C', '#7BA505','#00ffc0',  '#EA36DE', '#97FBF7', '#AAB8B6', '#00bfff','#A5286E', '#EAF2D9']  # , '#A5286E', , '#327BEB' '#1D9862'
    for i in range(0, N):
        #color = str('#') + '%06X' % randint(0, 0xFFFFFF) #generate random colors
        color = colors[i] #use fix colors
        membership_colors[i] = color
        print('Memberhsip %d color and node list are: %s --'%(i,color),cliques_set[i])
    print(membership_colors)

    for idx in graphs:
        G = graphs[idx]
        iG = nxGraph_to_iGraph(G)
        membership = get_membership(iG,cliques_set)
        output_directory = "output_files/egoNetwork_5/"
        file_name = os.path.join(output_directory, "cluster_network_"+str(idx)+".pdf")
        _plot(iG, file_name, membership_colors,membership)

