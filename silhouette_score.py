import numpy as np
from random import random
def silhouette_score(distance,cluster_dict,tau):
    silhouette_score = 0
    if len(cluster_dict)==0:
        return silhouette_score
    centroid_list = sorted(cluster_dict.keys())
    indices = len(centroid_list)
    s = {}
    a = {}
    b = {}
    if indices == 1:
        return silhouette_score
    for idx in range(indices):
        centroid = centroid_list[idx]
        for i in cluster_dict[centroid]:
            if len(cluster_dict[centroid]) == 1:
                if idx == 0:
                    next_G = cluster_dict[centroid_list[idx+1]][0]
                    if distance[i][next_G] > tau:
                        s[i] = 1
                    else:
                        s[i] = 0
                elif idx == indices - 1:
                    previous_G = cluster_dict[centroid_list[idx - 1]][-1]  # get the last element
                    if distance[i][previous_G] > tau:
                        s[i] = 1
                    else:
                        s[i] = 0
                else:
                    previous_G = cluster_dict[centroid_list[idx - 1]][-1]  # get the last element
                    next_G = cluster_dict[centroid_list[idx+1]][0]
                    if distance[i][previous_G] > tau and distance[i][next_G] > tau:
                        s[i] = 1
                    else:
                        s[i] = 0
                silhouette_score += s[i]
            elif len(cluster_dict[centroid]) >= 0.5*len(distance):
                s[i] = 0
                silhouette_score += s[i]
            else:
                a[i] = np.mean([distance[i][j] for j in cluster_dict[centroid] if not i == j])
                '''
                if idx == 0:
                    b[i] = np.mean([distance[i][j] for j in cluster_dict[centroid_list[idx+1]]])
                elif idx == indices - 1:
                    b[i] = np.mean([distance[i][j] for j in cluster_dict[centroid_list[idx - 1]]])
                else:
                    b[i] = np.min([np.mean([distance[i][j] for j in cluster_dict[centroid_list[idx - 1]]]), np.mean([distance[i][j] for j in cluster_dict[centroid_list[idx+1]]])])
                '''
                b[i] = np.min([np.mean([distance[i][j] for j in cluster_dict[centroid_list[k]]]) for k in range(indices) if not k == idx])
                s[i] = (b[i] - a[i])/np.max([a[i], b[i]])
                silhouette_score += s[i]
    print("s values are",s)
    print("distance len is",len(distance))
    return silhouette_score/len(distance)

if __name__ == "__main__":
    dist = {}
    for i in range(1,11):
        tmp = {}
        for j in range(1,11):
            tmp[j] = random()
        dist[i] = tmp
    cluster = {}
    cluster[1] = [1,2,3]
    cluster[4] = [4]
    cluster[6] =[5,6,7,8,9,10]
    print(silhouette_score(dist,cluster,0.5))


