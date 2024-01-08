import os
import networkx as nx
import numpy as np

#read cell name list
def read_celllist():
    f=open("/home/wanghy/deDoc2_clique/Nagano2017/cell_name.txt","r")
    cell_list=[]
    for line in f:
        cell_list.append(line.split('\n')[0])
    f.close()
    return cell_list


#read contact matrix
def read_contactfile(cell,software):
    path="/home/wanghy/deDoc2_clique/N_results/ted_sig_interaction/"+cell+"."+software+".sig"+".TAD"
    if os.path.getsize(path)>1024*2:
        f=open(path,"r")
        ls=[]
        for line in f:
            x=line.strip('\n').split('\t')
            ls.append([x[0:3],x[3:6]])
        f.close()
        ls=ls[4:]
        return ls
    else:
        return None


#Get a one-to-one dictionary
def get_dictionary(sig_list):
    d=dict()
    a=0
    for line in sig_list:
        if line[0] not in d.values():
            a=a+1
            d[a] = line[0]
        if line[1] not in d.values():
            a=a+1
            d[a] = line[1]
    return d


#Query the corresponding key based on the value
def getKey(dic,value):
    if value not in dic.values():
        return None
    result=[]
    for key in dic:
        if dic[key] == value:
            result.append(key)
    return result


#Get a list of edges and a list of points
def edges_nodes_list(d,sig_list):
    edges=[]
    for line in sig_list:
        key1=getKey(d,line[0])
        key2=getKey(d,line[1])
        edges.append([key1[0],key2[0]])
    nodes=list(d.keys()) 
    return edges,nodes
#e,n=edges_nodes_list(d,sig_list)


#Build the graph and convert it to dictionary list format
def graph(nodes,edges):
    G=nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    fun_input=nx.convert.to_dict_of_lists(G, nodelist=None)
    return fun_input


#Find clique
class Trace_Calls: 
#Credit goes to Richard Pattis for his Illustrate_Recursive class.
#Source:https://www.ics.uci.edu/~pattis/ICS-33/lectures/decoratorspackages.txt
    def __init__(self,f):
        self.f = f
        self.calls = 0
        self.trace = False
        self.record = []


    def illustrate(self,*args,**kargs):
        
        self.indent = 0
        self.trace = True
        answer = self.__call__(*args,**kargs)
        self.trace = False
        return answer
    
    # def __call__(self,*args,**kargs):  # bundle arbitrary arguments to this call
    #     self.calls += 1
    #     return self.f(*args,**kargs)  # unbundle arbitrary arguments to call f

    def display_records(self):
        return self.record

    def __call__(self,*args,**kargs):

        if self.trace:
            if self.indent == 0:
                print('Starting recursive illustration'+30*'-')
            print (self.indent*"."+"calling", self.f.__name__+str(args)+str(kargs))
            self.indent += 2
        self.calls += 1
        answer = self.f(*args,**kargs)
        if answer != None:
            self.record.append(answer)
        if self.trace:
            self.indent -= 2
            print (self.indent*"."+self.f.__name__+str(args)+str(kargs)+" returns", answer)
            if self.indent == 0:
                print('Ending recursive illustration'+30*'-')
        return answer
    def called(self):
        return self.calls

    def get_recursive_calls(self):
        return self.calls - 1
    
    def reset(self):
        self.calls = 0
        self.record = []


def trace(f): #Visualize recursive calls
    trace.recursive_calls = 0
    trace.depth = 0


    def _f(*args, **kwargs):


        print("  " * trace.depth, ">", f.__name__, args, kwargs)
        if trace.depth >= 1:
            trace.recursive_calls += 1
        trace.depth += 1
        res = f(*args, **kwargs)
        trace.depth -= 1
        print("  " * trace.depth, "<", res)
        print("recursive calls so far: {}".format(trace.recursive_calls))
        return res
    return _f
@Trace_Calls
def bron_kerbosch(R, P, X, graph, find_pivot=False):
    if len(P) == 0:
        if len(X) == 0:
            
            return R
    else:
        frontier = set(P)
        if find_pivot:
            #print("found_pivot")
            u = find_max_pivot(graph, P, X)
            #print(set(P), set(graph[u]))
            frontier = set(P) - set(graph[u])
        for v in frontier:
            # if DEBUG:
            #     print("BronKerbosch({}, {}, {})".format(
            #         R.union({v}),
            #         P.intersection(set(N(v,graph))),
            #         X.intersection(set(N(v,graph)))
            #         ))
            bron_kerbosch(
                R.union({v}),
                P.intersection(set(graph[v])),
                X.intersection(set(graph[v])),
                graph,
                find_pivot
             )

            P.remove(v)

            X = X.union({v})

def find_max_pivot(graph, P, X):
    nodes = list(P.union(X))
    u = nodes[0]
    max_intersection = len(set(graph[nodes[0]]).intersection(P))
    for n in nodes:
        if len(set(graph[n]).intersection(P)) > max_intersection:
            u = n
            max_intersection = len(set(graph[n]).intersection(P))

    return u


def bk_initial_call(graph, pivot=False, visualize=False):
    f = bron_kerbosch
    #print(f)

    if visualize:
        f.illustrate(set(), set(graph.keys()), set(), graph, pivot)
    else:
        f(set(), set(graph.keys()), set(), graph, pivot)
    #print(f.get_recursive_calls())
    #print(f.display_records())
    x=f.display_records()
    f.reset()
    return x
def N(v, g):
    # for i, n_v in enumerate(g[v]):
    #     print(i, n_v)
    #print("{}->{}".format(v,[n_v for i, n_v in enumerate(g[v]) if n_v]))

    return [n_v for i, n_v in enumerate(g[v]) if n_v]



#clique file processing
def process_clique (d,clique,cell,software):
    num_TAD=[]
    TAD_info=[]
    TAD_length=[]
    for item in clique:
        item=list(item)
        num_TAD.append(len(item))
        tad_info=[]
        tad_length=[]
        if len(item)>=3:
            for i in item:
                tad=d.get(i)
                tad_info.append(tad)
                tad_length.append(str(eval(tad[2])-eval(tad[1])))
            TAD_info.append(tad_info)
            TAD_length.append(tad_length)
    num_TAD1=np.unique(num_TAD,return_counts = True)
    
    #Location interval file
    f=open("/home/wanghy/deDoc2_clique/N_results/clique_results/"+cell+"."+software+".TAD_info"+".TAD","w")
    for item in TAD_info:
        x=[]
        i=0
        try:
            while item[i]:
                x.extend(item[i])
                i=i+1     
        except:
            pass 
        f.write("\t".join(x)+"\n")
    f.close() 

    #TLD length file
    f=open("/home/wanghy/deDoc2_clique/N_results/clique_results/"+cell+"."+software+".TAD_length"+".TAD","w")
    for item in TAD_length:
        f.write("\t".join(item)+"\n")
    f.close()

    #clique size statistics file
    np.savetxt("/home/wanghy/deDoc2_clique/N_results/clique_results/"+cell+"."+software+".num_TAD"+".TAD",num_TAD1,fmt='%d')
