###########################
# Author: Te-Yuan Liu
###########################

###########################
# Import Library
###########################
library("igraph")

###########################
# Define Function
###########################
recommend = function(g, uid, numr, measure){
    neighbor_vec = neighbors(g, which(V(g)$id==uid))$id
    V(g)$measurement = rep(-1, vcount(g))
    recom_id_list = rep(-1, numr)
    for(i in V(g)$id){
        if(!(i %in% c(uid, neighbor_vec))){
            metric = -1
            tmp_vec = neighbors(g, which(V(g)$id==i))$id
            if(measure=="c"){
                metric = length(intersect(neighbor_vec, tmp_vec))
            }
            else if(measure=="j"){
                metric = length(intersect(neighbor_vec, tmp_vec))/length(union(neighbor_vec, tmp_vec))
            }
            else if(measure=="a"){
                j_set = intersect(neighbor_vec, tmp_vec)
                if(length(j_set) == 0){
                    metric = 0
                }
                else{
                    metric = 0
                    for(j in j_set){
                        j_vec = neighbors(g, which(V(g)$id==j))
                        if(length(j_vec) <= 1){
                            # exception case, skip it
                        }
                        else{
                            metric = metric + 1/log(length(j_vec))
                        }
                    }
                }
            }
            V(g)[V(g)$id==i]$measurement = metric
        }
    }
    n = sort(V(g)$measurement, decreasing=TRUE, index.return=TRUE)
    for(i in 1:length(recom_id_list)){
        recom_id_list[i] = V(g)[n$ix[1]]$id
    }
    return(recom_id_list)
}
accuracy_measure = function(g, measure){
    Nr = V(g)[V(g)$degree==4]$id
    print("************")
    print(Nr)
    acc1 = 0.0
    for(i in Nr){
        print("+++++++++++++++++++")
        acc2 = 0.0
        for(j in 1:10){
            print("///////////////////")
            p = induced_subgraph(g, V(g))
            neighbors_id = neighbors(p, which(V(p)$id==i))
            num_removal = ceiling(0.25*length(neighbors_id))
            if(num_removal==0){
                print("no removal error")
            }
            else{
                deleted_id = neighbors_id[sample(length(neighbors_id), num_removal, replace=FALSE)]
                print(deleted_id)
                num_e = ecount(p)
                for(k in deleted_id){
                    plot(p)
                    x = as.character(which(V(p)$id==i))
                    y = as.character(which(V(p)$id==k))
                    z = paste(x,y,sep="|")
                    p = delete_edges(p, z)
                    plot(p)
                    if(ecount(p)==num_e){
                        print("edge number error")
                    }
                }
                print(recommend(p,i,num_removal, measure))
                print(intersect(recommend(p, i, num_removal, measure), deleted_id))
                acc2 = acc2 + length(intersect(recommend(p, i, num_removal, measure), deleted_id))/num_removal
            }
        }
        acc2 = acc2/10
        acc1 = acc1 + acc2
    }
    return(acc1/length(Nr))
}
###########################
# Main Function
###########################
main = function(){
    g = read.graph("test.txt", format="edgelist", directed=FALSE)
    V(g)$id = 1:vcount(g)
    V(g)$degree = degree(g)
    #x = as.character(1)
    #y = as.character(4)
    #z = paste(x,y,sep="|")
    #g = delete_edges(g,z)
    #plot(g)
    print(accuracy_measure(g, "c"))
}
main()
