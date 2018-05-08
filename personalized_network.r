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
    V(p)$ratio = rep(-1, vcount(p))
    V(p)$product = rep(-1, vcount(p))
    
    for(i in V(p)$id){
        if(i != core_id){
            print(i)
            mn_id_vec = intersect(neighbors(p, which(V(p)$id==i))$id, neighbors(p, which(V(p)$id==core_id))$id)
            if(length(mn_id_vec) == 0){
                V(p)[V(p)$id==i]$embeddedness = 0
                V(p)[V(p)$id==i]$dispersion = 0
                V(p)[V(p)$id==i]$ratio = 0
                V(p)[V(p)$id==i]$product = 0
            }
            else{
                V(p)[V(p)$id==i]$embeddedness = length(mn_id_vec)
                p_del = delete.vertices(p, c(which(V(p)$id==i), which(V(p)$id==core_id)))
                dist_sum = 0
                inf_dist = FALSE
                if(length(mn_id_vec)==1){
                    dist_sum = 0  
                }
                else{
                    dist_all = shortest.paths(p_del, mode="all")
                    for(j in 1:(length(mn_id_vec)-1)){
                        if(inf_dist){
                            break
                        }
                        v_j = mn_id_vec[j]
                        for(k in (j+1):length(mn_id_vec)){
                            if(inf_dist){
                                break
                            }
                            v_k = mn_id_vec[k]
                            dist = dist_all[which(V(p_del)$id==v_j), which(V(p_del)$id==v_k)]
                            if(dist == Inf){
                                inf_dist = TRUE
                            }
                            else{
                                dist_sum = dist_sum + dist
                            }
                        }
                    }

                }
                if(inf_dist){
                    dist_sum = diameter(g) + 1
                }
                V(p)[V(p)$id==i]$dispersion = dist_sum
                V(p)[V(p)$id==i]$ratio = V(p)[V(p)$id==i]$dispersion/V(p)[V(p)$id==i]$embeddedness
                V(p)[V(p)$id==i]$product = V(p)[V(p)$id==i]$dispersion*V(p)[V(p)$id==i]$embeddedness

            }

        }
    }
    id_max_pro = which.max(V(p)$product)
    pro_edges = incident(p, V(p)[id_max_pro], mode="all")
    # plot community and highlight largest dispersion
    fg = fastgreedy.community(p)
    #lo = layout.fruchterman.reingold(p)
    ecol = rep("gray80", ecount(p))
    ecol[pro_edges] = "red"
    vcol = membership(fg)
    vcol[id_max_pro] = "red"
    vcol[V(p)$id==core_id] = "gold"
    vsize = rep(4, vcount(p))
    vsize[id_max_pro] = 20
    vsize[V(p)$id==core_id] = 10
    vlabelsize = rep(0.2, vcount(p))
    vlabelsize[id_max_pro] = 1
    vlabelsize[V(p)$id==core_id] = 0.5
    plot(fg, p, mark.groups=NULL, edge.color=ecol, col=vcol, vertex.size=vsize, vertex.label=V(p)$id, vertex.label.cex=vlabelsize, vertex.label.color="black")
    
    #id_max_disp = which.max(V(p)$dispersion)
    #disp_edges = incident(p, V(p)[id_max_disp], mode="all")
    # plot community and highlight largest dispersion
    #fg = fastgreedy.community(p)
    #lo = layout.fruchterman.reingold(p)
    #ecol = rep("gray80", ecount(p))
    #ecol[disp_edges] = "red"
    #vcol = membership(fg)
    #vcol[id_max_disp] = "red"
    #vcol[V(p)$id==core_id] = "gold"
    #vsize = rep(4, vcount(p))
    #vsize[id_max_disp] = 20
    #vsize[V(p)$id==core_id] = 10
    #vlabelsize = rep(0.2, vcount(p))
    #vlabelsize[id_max_disp] = 1
    #vlabelsize[V(p)$id==core_id] = 0.5
    #plot(fg, p, mark.groups=NULL, edge.color=ecol, col=vcol, vertex.size=vsize, vertex.label=V(p)$id, vertex.label.cex=vlabelsize, vertex.label.color="black")

    # plot community and highlight largest embeddedness and dispersion/embeddedness ratio
    #id_max_emb = which.max(V(p)$embeddedness)
    #id_max_ratio = which.max(V(p)$ratio)
    #ecol = rep("gray80", ecount(p))
    #vcol = membership(fg)
    #vsize = rep(4, vcount(p))
    #vlabelsize = rep(0.2, vcount(p))

    #if(id_max_emb == id_max_ratio){
    #    emb_edges = incident(p, V(p)[id_max_emb], mode="all")
    #    ecol[emb_edges] = "green"
    #    vcol[id_max_emb] = "green"
    #    vcol[V(p)$id==core_id] = "gold"
    #    vsize[id_max_emb] = 20
    #    vsize[V(p)$id==core_id] = 10
    #    vlabelsize[id_max_emb] = 1
    #    vlabelsize[V(p)$id==core_id] = 0.5
    #}
    #else{    
    #    emb_edges = incident(p, V(p)[id_max_emb], mode="all")
    #    ratio_edges = incident(p, V(p)[id_max_ratio], mode="all")
    #    ecol[emb_edges] = "red"
    #    ecol[ratio_edges] = "blue"
    #    vcol[id_max_emb] = "red"
    #    vcol[id_max_ratio] = "blue"
    #    vcol[V(p)$id==core_id] = "gold"
    #    vsize[id_max_emb] = 20
    #    vsize[id_max_ratio] = 20
    #    vsize[V(p)$id==core_id] = 10
    #    vlabelsize[id_max_emb] = 1
    #    vlabelsize[id_max_ratio] = 1
    #    vlabelsize[V(p)$id==core_id] = 0.5
    #}
    #plot(fg, p, mark.groups=NULL, edge.color=ecol, col=vcol, vertex.size=vsize, vertex.label=V(p)$id, vertex.label.cex=vlabelsize, vertex.label.color="black")
    
}

##############################
# Main Function
##############################
main = function(){
    g = read.graph("facebook_combined.txt", format="edgelist", directed=FALSE)
    V(g)$id = 1:vcount(g)
    plot_dispersion_embeddedness(g, 1)
    plot_dispersion_embeddedness(g, 108)
    plot_dispersion_embeddedness(g, 349)
    plot_dispersion_embeddedness(g, 484)
    plot_dispersion_embeddedness(g, 1087)

}
main()
