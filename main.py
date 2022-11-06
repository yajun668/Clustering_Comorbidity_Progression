#main function
import Graph
import Gurobi_Interface
import NetworkModel


def clustered_network_results():
    directory = "input_files/egoNetwork_26/"
    graphs = Graph.read_graphs(directory)
    n = len(graphs)
    distances = Graph.read_distances("input_files/distances_26.csv")
    tau = 0.5
    solutions = Gurobi_Interface.network_cluster(n, distances, 5, tau)
    # based on the clusters, we aggreage comorbidity and re-run network model
    new_clusters = {}
    num = 1
    for key in solutions:
        new_clusters[num] = solutions[key].copy()
        num += 1
    NetworkModel.network_model(new_clusters) #new ego network stored in the "input_files/egoNetwork/"
    directory2 = "input_files/egoNetwork/"
    Graph.ego_network_output_visualization(directory2)
def clusterNetwork_stat_figure():
    directory = "input_files/egoNetwork_5/"
    graphs = Graph.read_graphs(directory)
    n = len(graphs)
    cliques = Graph.clique_signature(graphs)
    Graph.ego_network_chart_clustered(graphs, cliques)

def generate_csv_for_egoNetwork26():
    directory = "input_files/egoNetwork_26/"
    graphs = Graph.read_graphs(directory)
    Graph.to_csv_egonetwork(graphs)

def generate_csv_for_aggregatedNetwork():
    directory = "input_files/egoNetwork_5/"
    graphs = Graph.read_graphs(directory)
    Graph.ego_network_output_visualization(directory)

def find_cluster_p():
    directory = "input_files/egoNetwork_26/"
    graphs = Graph.read_graphs(directory)
    #distances = Graph.get_network_distance(graphs) #only need to run one time and it generate a csv file; can be used next time using below function
    distances = Graph.read_distances("input_files/distances_26.csv")
    tau = 0.5
    Graph.find_cluster_value_by_Silhouette(distances, tau) # 5 is the optimal value

#In main(), when run one function, comment all other functions
if __name__ == "__main__":
    #find_cluster_p() #find the cluster p value based on Silhouette value
    clustered_network_results()
    #generate_csv_for_egoNetwork26()
    generate_csv_for_aggregatedNetwork()
    clusterNetwork_stat_figure()




