##############################
# Author: Te-Yuan Liu
##############################

##############################
# Import Library
##############################
library(igraph)

##############################
# Define Function
##############################
plot_dispersion_embeddedness = function(g, core_id){
    p = induced.subgraph(g, c(which(V(g)$id==core_id), neighbors(g, which(V(g)$id==core_id))))
    V(p)$embeddedness = rep(-1, vcount(p))
    V(p)$dispersion = rep(-1, vcount(p))
    
    for(i in V(p)$id){
        if(i != core_id){
            print(i)
            mn_id_vec = intersect(neighbors(p, which(V(p)$id==i))$id, neighbors(p, which(V(p)$id==core_id))$id)
            if(length(mn_id_vec) == 0){
                V(p)[which(V(p)$id==i)]$embeddedness = 0
                V(p)[which(V(p)$id==i)]$dispersion = 0
            }
            else{
                V(p)[which(V(p)$id==i)]$embeddedness = length(mn_id_vec)
                p_del = delete.vertices(p, c(which(V(p)$id==i), which(V(p)$id==core_id)))
                dist_sum = 0
                if(length(mn_id_vec)==1){
                    dist_sum = 0  
                }
                else{
                    dist_all = shortest.paths(p_del, mode="all")
                    for(j in 1:(length(mn_id_vec)-1)){
                        v_j = mn_id_vec[j]
                        for(k in (j+1):length(mn_id_vec)){
                            v_k = mn_id_vec[k]
                            dist = dist_all[which(V(p_del)$id==v_j), which(V(p_del)$id==v_k)]
                            if(dist == Inf){
                                dist_sum = dist_sum + diameter(p)
                            }
                            else{
                                dist_sum = dist_sum + dist
                            }
                        }
                    }

                }
                V(p)[which(V(p)$id==i)]$dispersion = dist_sum
            }

        }
    }
    id_max_disp = V(p)[which.max(V(p)$dispersion)]$id

    print(id_max_disp)
}

##############################
# Main Function
##############################
main = function(){
    g = read.graph("facebook_combined.txt", format="edgelist", directed=FALSE)
    V(g)$id = 1:vcount(g)
    plot_dispersion_embeddedness(g, 1)
}
main()
