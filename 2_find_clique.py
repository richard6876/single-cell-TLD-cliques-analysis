import networkx as nx
import numpy as np
from function2 import  read_contactfile, get_dictionary, getKey, edges_nodes_list, graph, bk_initial_call, process_clique,read_celllist
import os

# get the TLD cliques per cell
if __name__ == '__main__':

    cell_list=read_celllist()

    for cell in cell_list:
        sig_list=read_contactfile(cell,"de2w")
        if sig_list == None:
           continue
        d=get_dictionary(sig_list)
        edges,nodes=edges_nodes_list(d,sig_list)
        fun_input=graph(nodes,edges)
        clique=bk_initial_call(fun_input)
        process_clique (d,clique,cell,"de2w")