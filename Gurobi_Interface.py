from gurobipy import *
import networkx as nx
import pygraphviz as pgv
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import csv

#Currently this function is not used for our paper
def clubs_model_solve(g,index_time_window):
    #This function is to find one-neighborhood (or egonet work)/2club/robust club/clique
    try:
        # Model
        m = Model("C-Diff")
        print('* Building the model of time window ' + str(index_time_window)+ '...\n')
        # Creating decision variables:
        x_vars = m.addVars(g.nodes(), vtype=GRB.BINARY, name="x")
        # constraints
        #if c-diff (node id=0) is in the node list, set x[0] = 1
        #find the subgraph induced by one-neighborhood of node 0; or called ego network
        for n1 in g.nodes():
            if n1 == 0:
                continue
            if n1 in g.neighbors(0):
                m.addConstr(x_vars[n1]== 1, "com_constr_" + str(n1))
            else:
                m.addConstr(x_vars[n1] == 0, "com_constr_" + str(n1))
        '''
        # Cliuqe Model constraints
        for n1 in g.nodes():
            for n2 in g.nodes():
                if n2 not in g.neighbors(n1) and n1 < n2:
                    m.addConstr(x_vars[n1] + x_vars[n2] <= 1, "com_constr_"+ str(n1)+','+str(n2))
        # 2-club neighborhood constraints
        for n1 in g.nodes():
            for n2 in g.nodes():
                if (n1<n2) and (n2 not in g.neighbors(n1)):
                    m.addConstr(x_vars[n1] + x_vars[n2] <= 1 + quicksum(x_vars[i] for i in nx.common_neighbors(g,n1,n2)),"com_constr_"+ str(n1)+','+str(n2))
        # r robust 2 club model constraints:
        r = 4 # try 4-robust 2-club
        for n1 in g.nodes():
            for n2 in g.nodes():
                if (n1 < n2):
                    if (n2 in g.neighbors(n1)):
                        m.addConstr((r-1)*(x_vars[n1] + x_vars[n2]-1) <= quicksum(x_vars[i] for i in nx.common_neighbors(g, n1, n2)), "com_constr_" + str(n1) + ',' + str(n2))
                    else:
                        m.addConstr(r*(x_vars[n1] + x_vars[n2]-1) <= quicksum(x_vars[i] for i in nx.common_neighbors(g, n1, n2)), "com_constr_" + str(n1) + ',' + str(n2))
        # k-plex model constraints
        k=2
        for n1 in g.nodes():
            tmpS = []
            for j in g.nodes():
                if j == n1:
                    continue
                if j not in g.neighbors(n1):
                    tmpS.append(j)
            m.addConstr(quicksum(x_vars[j] for j in tmpS) <= (k-1)*x_vars[n1] + (len(g.nodes())-g.degree(n1)-1)*(1-x_vars[n1]), "com_constr_" + str(n1))
        '''
        # Set optimization MAXIMIZE
        m.setObjective(x_vars.sum('*'), GRB.MAXIMIZE)
        m.setParam("OutputFlag", 0)  # turn off the log printing
        # Solve
        m.optimize()
        # Write model to a file
        #m.write('formulation.lp')
        print('*******************************************************************************')
        # Status checking
        status = m.Status
        solution_set = []
        if status == GRB.Status.OPTIMAL:
            print('\nOptmail value: %g' % m.objVal)
            print("num vars: ", m.NumVars)
            print("num constrs: ", m.NumConstrs)
            for node in g.nodes():
                if x_vars[node].x > 0.9:
                    solution_set.append(node);
             #visulize solutions
            B = pgv.AGraph()
            B.node_attr['color'] = 'blue'
            B.node_attr['shape'] = 'circle'
            B.node_attr['style'] = 'filled'
            B.node_attr['fixedsize'] = 'true'
            B.node_attr['width'] = 0.3
            B.node_attr['height'] = 0.3
            B.node_attr['fillcolor'] = 'red'
            B.edge_attr['color'] = 'blue'
            B.edge_attr['penwidth'] = '0.2'
            edge_csv_name = os.path.join("output_files/egoNetwork", 'edge_list_disease_networks_LOS=' + str(index_time_window) + '.csv')
            f1 = open(edge_csv_name,'w')
            edge_writer = csv.writer(f1)
            edge_writer.writerow(['Source','Target','Type','Id','Label','timeset','Weight'])
            id = 0
            g1 = nx.Graph()
            for n1 in g.nodes():
                for n2 in g.nodes():
                    if n1<n2 and n2 in g.neighbors(n1) and n1 in solution_set and n2 in solution_set:
                        B.add_edge(n1,n2)
                        g1.add_edge(n1,n2)
                        edge_writer.writerow([n1, n2, 'Undirected', id, '', '', 1])
                        id += 1
            f1.close()
            B.layout(prog='sfdp') # use sfdp; layout with default (neato)
            pngfile = 'solution_networks_LOS=' + str(index_time_window) + '.png'
            B.draw(os.path.join("output_files/egoNetwork", pngfile))  # draw png
            # write Graph G to GraphML format for later writing csv file use
            graph_name1 = str(index_time_window) + '_disease_networks_LOS=' + '.graphml'
            #this is for Gephi
            if index_time_window == 2:
                graph_name_current =str(index_time_window-1) + '_disease_networks_LOS=' + '.graphml'
                g_current = nx.read_graphml(os.path.join("output_files/egoNetwork", graph_name_current))
                red_nodes = []
                purple_nodes = []
                for node in g_current.nodes:
                    if int(node) in g1.nodes:
                        purple_nodes.append(node) #note that node is a string type
                    else:
                        red_nodes.append(node)
                csv_name = os.path.join("output_files/egoNetwork", 'disease_networks_LOS=' + str(index_time_window-1) + '.csv')
                with open(csv_name,'w') as f:
                    writer = csv.writer(f)
                    writer.writerow(['Id','Label','timeset','modularity_class'])
                    for node in g_current.nodes():
                        if node == '0':
                            writer.writerow([node, node, '', 0])
                        elif node in purple_nodes:
                            writer.writerow([node,node,'',3])
                        else:
                            writer.writerow([node, node, '', 4])
            elif index_time_window > 2 and index_time_window <= 26:
                graph_name_previous = str(index_time_window - 2) + '_disease_networks_LOS=' +  '.graphml'
                g_previous = nx.read_graphml(os.path.join("output_files/egoNetwork", graph_name_previous))
                graph_name_current = str(index_time_window-1) + '_disease_networks_LOS=' +  '.graphml'
                g_current = nx.read_graphml(os.path.join("output_files/egoNetwork", graph_name_current))
                blue_nodes = []
                aqua_nodes = []
                red_nodes = []
                purple_nodes = []
                for node in g_current.nodes:
                    if (node in g_previous.nodes) and (int(node) in g1.nodes):
                        blue_nodes.append(node) #note that node is a string type
                    if (node in g_previous.nodes) and (int(node) not in g1.nodes):
                        aqua_nodes.append(node) #note that node is a string type
                    if (node not in g_previous.nodes) and (int(node) in g1.nodes):
                        purple_nodes.append(node) #note that node is a string type
                    if (node not in g_previous.nodes) and (int(node) not in g1.nodes):
                        red_nodes.append(node) #note that node is a string type
                csv_name = os.path.join("output_files/egoNetwork", 'disease_networks_LOS=' + str(index_time_window-1) + '.csv')
                with open(csv_name,'w') as f:
                    writer = csv.writer(f)
                    writer.writerow(['Id','Label','timeset','modularity_class'])
                    for node in g_current.nodes():
                        if node == '0':
                            writer.writerow([node, node, '', 0])
                        elif node in blue_nodes:
                            writer.writerow([node, node, '', 1])
                        elif node in aqua_nodes:
                            writer.writerow([node, node, '', 2])
                        elif node in purple_nodes:
                            writer.writerow([node,node,'',3])
                        else:
                            writer.writerow([node, node, '', 4])
                if index_time_window == 26:
                    #this will generate csv for the 26th graph
                    graph_name_previous = str(index_time_window - 1) + '_disease_networks_LOS=' +  '.graphml'
                    g_previous = nx.read_graphml(os.path.join("output_files/egoNetwork", graph_name_previous))
                    aqua_nodes = []
                    red_nodes = []
                    for node in g1.nodes:
                        if str(node) in g_previous.nodes:
                            aqua_nodes.append(node)  # note that node is a string type
                        else:
                            red_nodes.append(node)  # note that node is a string type
                    csv_name = os.path.join("output_files/egoNetwork",
                                            'disease_networks_LOS=' + str(index_time_window) + '.csv')
                    with open(csv_name, 'w') as f:
                        writer = csv.writer(f)
                        writer.writerow(['Id', 'Label', 'timeset', 'modularity_class'])
                        for node in g1.nodes:
                            if node == 0:
                                writer.writerow([node, node, '', 0])
                            elif node in aqua_nodes:
                                writer.writerow([node, node, '', 2])
                            else:
                                writer.writerow([node, node, '', 4])
            nx.write_graphml(g1,os.path.join("output_files/egoNetwork", graph_name1))
            graph_name = str(index_time_window) + '_edge_disease_networks_LOS' + '.txt' #this file is for Network Simile code use
            nx.write_edgelist(g1, os.path.join("output_files/egoNetwork", graph_name), data=False)

        if status == GRB.Status.INF_OR_UNBD or \
                status == GRB.Status.INFEASIBLE or \
                status == GRB.Status.UNBOUNDED:
            m.computeIIS()
            m.write('formulation.ilp')
            print('The model cannot be solved because it is infeasible or unbounded')
            sys.exit(1)

        if status != GRB.Status.OPTIMAL:
            print('Optimization was stopped with status ' + str(status))
            sys.exit(1)
        print('*******************************************************************************')
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))
    except AttributeError:
        print('Encountered an attribute error')

#this function is modified p-median model for network clustering.
def network_cluster(n,d,p,tau):
    try:
        N = range(1,n+1) #Generate a list from 1 to n
        # Model
        m = Model("Network_Clustering")
        print('* Building the model...\n')
        # Creating decision variables:
        x_vars = m.addVars(N, N, vtype=GRB.BINARY, name="x")
        # constraints
        tmp_sum = 0
        for j in N:
            tmp_sum += x_vars[j,j]
            m.addConstrs((x_vars[i,j] + x_vars[i+1,j] <= 1 for i in range(1,n) if d[i][i+1] >= tau),'threshold contraints')
        m.addConstr(tmp_sum <= p,'p-median number constr')
        m.addConstrs((x_vars.sum(i,'*') == 1 for i in N),'C7 constr')
        for i in range(1,n-1):
            for j in range(i+2,n+1):
                for k in range(i+1,j):
                    m.addConstr(x_vars[i,j] <= x_vars[k,j],'C8-[%s,%s,%s]'%(str(i),str(j),str(k)))
        for i in range(3,n+1):
            for j in range(1,i-1):
                for k in range(j+1,i):
                    m.addConstr(x_vars[i,j] <= x_vars[k,j],'C9-[%s,%s,%s]'%(str(i),str(j),str(k)))
        '''
        m.addConstr(x_vars[1,1] + x_vars[1,2]==1, 'Network 1 constr')
        m.addConstr(x_vars[n,n-1]+x_vars[n,n] == 1, 'Network n constr')
        for i in range(2,n):
            m.addConstr(x_vars[i,i-1] + x_vars[i,i] + x_vars[i,i+1] == 1, 'Network-' + str(i) + '-constr')
        '''
        for i in N:
            for j in N:
                if i == j:
                    continue
                m.addConstr(x_vars[i,j] <= x_vars[j,j], 'Median constr %s,%s'%(i,j))

        # Set optimization MAXIMIZE
        m.setObjective(quicksum(d[i][j]*x_vars[i,j] for i in N for j in N), GRB.MINIMIZE)
        m.setParam("OutputFlag",0) #turn off the log printing
        # Solve
        m.optimize()
        # Write model to a file
        #m.write('formulation.lp')
        print('*******************************************************************************')
        # Status checking
        status = m.Status
        solution_dict = {}
        if status == GRB.Status.OPTIMAL:
            print('\nOptmail value: %g' % m.objVal)
            print("num vars: ", m.NumVars)
            print("num constrs: ", m.NumConstrs)
            for i in N:
                for j in N:
                    if x_vars[i,j].x>0.9:
                        if j in solution_dict:
                            solution_dict[j].append(i)
                        else:
                            solution_dict[j]= [i]
            print("There are %d clusters"%len(solution_dict))
            for key1 in solution_dict:
                print("The cluster centroid %d includes networks:"%key1, solution_dict[key1])
            return solution_dict
            #return m.getObjective().getValue()
        if status == GRB.Status.INF_OR_UNBD or \
                status == GRB.Status.INFEASIBLE or \
                status == GRB.Status.UNBOUNDED:
            m.computeIIS()
            m.write('formulation.ilp')
            print('The model cannot be solved because it is infeasible or unbounded')
            #sys.exit(1)
            return solution_dict
            #return -1
        if status != GRB.Status.OPTIMAL:
            print('Optimization was stopped with status ' + str(status))
            #sys.exit(1)
            return solution_dict
            #return -1
        print('*******************************************************************************')
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))
    except AttributeError:
        print('Encountered an attribute error')

def maximum_clique(graph_list, common_nodes):
    try:
        D = common_nodes
        # Model
        m = Model("maximum clique")
        print('* Building the model...\n')
        # Creating decision variables:
        x_vars = m.addVars(D, vtype=GRB.BINARY, name="x")
        # constraints
        for key in graph_list:
            g = graph_list[key]
            sub_Graph = g.subgraph(D) # the subgraph induced by set D
            for i in D:
                for j in D:
                    if int(i) < int(j) and i not in sub_Graph.neighbors(j):
                        m.addConstr(x_vars[i] + x_vars[j] <= 1, 'constraints %s,%s'%(i,j))
        # Set optimization MAXIMIZE
        m.setObjective(x_vars.sum('*'), GRB.MAXIMIZE)
        m.setParam("OutputFlag",0) #turn off the log printing
        # Solve
        m.optimize()
        # Write model to a file
        #m.write('formulation.lp')
        print('*******************************************************************************')
        # Status checking
        status = m.Status
        solution = []
        if status == GRB.Status.OPTIMAL:
            print('\nOptmail value: %g' % m.objVal)
            print("num vars: ", m.NumVars)
            print("num constrs: ", m.NumConstrs)
            for i in D:
                if x_vars[i].x > 0.9:
                    solution.append(i)
            print("The maximum clique is: ", solution)
            return solution
        if status == GRB.Status.INF_OR_UNBD or \
                status == GRB.Status.INFEASIBLE or \
                status == GRB.Status.UNBOUNDED:
            m.computeIIS()
            m.write('formulation.ilp')
            print('The model cannot be solved because it is infeasible or unbounded')
            sys.exit(1)
        if status != GRB.Status.OPTIMAL:
            print('Optimization was stopped with status ' + str(status))
            sys.exit(1)
        print('*******************************************************************************')
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))
    except AttributeError:
        print('Encountered an attribute error')