import numpy as np
import networkx as nx
import random
from itertools import product
from itertools import combinations
from collections import Counter
import matplotlib.pyplot as plt
from sklearn import linear_model
import json
import copy
from scipy.stats import pareto
from scipy.stats import poisson




def list_split(items, n):
    #ivide the items evenly into n parts
    np.random.shuffle(items)
    items_random=list(items)
    d=int(len(items_random)/n)
    #print(d)
    A=[sorted(items_random[i * d:(i + 1) * d]) for i in range(n-1)]
    A.append(sorted(items_random[(n-1) * d:]))
    # A=[items_random[i * d:(i + 1) * d] for i in range(n-1)]
    # A.append(items_random[(n-1) * d:])
    return A

def num_split(num,n):
    #Split the integer num into n values randomly, with their sum equal to num
    num_n=random.sample(range(1,num),n-1)
    num_n.sort()
    num_n.append(num)
    num_n2=[]
    for i in range(len(num_n)):
        if i == 0:
            b=num_n[i]
        else:
            b=num_n[i] - num_n[i-1]
        num_n2.append(b)
    return num_n2

def list_split_name(key,items, n):#Divide the items into n parts. If the number of items is less than n, set n equal to the number of items.
    B={}
    if items<n:
        n=items
    A=num_split(items, n)
    for i in range(len(A)):
        B[key[i]]=A[i]
    return B


def func1(amount, num):
    #Generate a probability combination where the sum of probabilities is 1
    #random.seed(seed)
    list1 = []
    for i in range(0, num-1):
        a = random.random()   
        list1.append(a)
    list1.sort()                        
    list1.append(amount)                
    list2 = []
    for i in range(len(list1)):
        if i == 0:
            b = list1[i]                
        else:
            b = list1[i] - list1[i-1]   
        list2.append(b)
    return list2


def random_pick(some_list, prb):
    #Select an element with probabilit
    ran=random.uniform(0,1)
    #print("ran",ran)
    cumulative_prob=0.0
    for item, item_prob in zip(some_list,prb):
        cumulative_prob += item_prob
        if ran < cumulative_prob: break
    return item

def power_law_P(x_min,x_max,a):#Generate probabilities that follow a power-law distribution
    x = [float(i) for i in range(x_min, x_max)]
    p = pareto.pdf(x, a)  # x=a/x^(a+1)
    p_sum = sum(p)
    p_norm = p / p_sum
    #print("p_norm", p_norm, sum(p_norm))
    return p_norm

def permute(nums):
    #An array forming different permutations of order
    from itertools import permutations
    result = []
    for i in permutations(nums, len(nums)):
        result.append(list(i))
    return result



def generate_edage(nodelist,possion_mu):
    list_edges = []
    list_node =[]
    nodelist_copy=set(copy.deepcopy(nodelist))
    while len(nodelist_copy) != 0:
        size = poisson.rvs(possion_mu)
        if size >1:
            if len(nodelist_copy) <size:
                list_edges[-1].extend(list(nodelist_copy))
                list_node.extend(list(nodelist_copy))
                nodelist_copy = set()
            else:
                new_edge = np.random.choice(list(nodelist_copy),size=size,replace=False)
                nodelist_copy -= set(new_edge)
                list_edges.append(list(new_edge))
                list_node.extend(list(new_edge))
    #print('list_node',len(list_node))
    return list_edges


def combinations_fun(name, d_num):
    # An array named 'name', select d_num elements to form different order permutations
    result = []
    for com in combinations(name, d_num):
        p = np.zeros(C_num, dtype=int)
        for j in com:
            p[name == j]=1
        result.append(list(p))
    return result

def combinations_fun_all(name, d_num):
    #From the array of club names, select d_num clubs to form different ordered permutations
    result = []
    for num_i in range(d_num):
        for com in combinations(name, num_i+1):
            p = np.zeros(C_num, dtype=int)
            for j in com:
                p[name == j] = 1
            result.append(list(p))
    return result


def P_function(com_all_num, alpha):
    P=np.zeros(com_all_num,dtype=float)
    P[:C_num]=(1 - alpha)/C_num
    P[C_num:]=alpha/(com_all_num - C_num)
    return P


def calculate_prob(degrees,gamma):
    #The probability of node degree
    prob = dict()
    deg_sum = np.sum([(deg+gamma) for deg in degrees.values()])
    for node, degree in degrees.items():
        prob[node] = (degree+gamma) /deg_sum
    return prob


def Si_calculate(C_num,com_all,P):
    Si_all=[]
    com_all = np.array(com_all)
    for i in range(C_num):
        node_i_combine=com_all[com_all[:,i]== 1]
        #print("node_i_combine",node_i_combine)
        Si=sum(P[com_all[:, i] == 1])
        #print("Si",Si)
        Si_all.append(Si)
    return Si_all

class Hypergraph:
    def __init__(self):
        self.nodes_list = []#node
        self.edge ={}#edge
        self.degree={}#degree
        self.community={}#community
        self.node_i=0
        self.node_inedge = set()
        self.X=[]#Hyperedge length
        self.edge_between={}#community hyperedge

    #def hypergraph_PA(self,C0,p,M,p1, P,update_iter=10,Node_num_T=10):#X, P, r
    def hypergraph_PA(self,C0,p,M,p1, P,Node_num_T=10):
        self.community=C0
        for key,value in C0.items():
            self.nodes_list.extend(value)
            self.degree[key]=Counter(value)
        self.edge[0]=copy.deepcopy(self.nodes_list)
        self.nodes_list=set(self.nodes_list)
        self.node_i=len(self.nodes_list)
        self.some_list=list(range(len(P)))

        x_min = 2#Range of hyperedge lengths
        x_max = 21 
        p_norm = power_law_P(x_min, x_max, 1.25)
        edge_length = [i for i in range(x_min, x_max)]

        t=1
        #for t in range(1,update_iter+1):
        while len(self.nodes_list) < Node_num_T:
            # print("t",t)
            rand=random.uniform(0,1)
            #print("rand",rand)
            if rand <= p:#add node
                self.nodes_list.add(self.node_i)
                #print("self.node_i",self.node_i)
                item = random_pick(C_name, M)
                self.community[item].append(self.node_i)
                #print("self.degree[item]",self.degree[item])
                self.degree[item]=self.degree[item]+Counter([self.node_i])
                #print("self.degree[item]", self.degree[item])
                self.node_i += 1
                # print("self.community",self.community)
                #print("self.nodes_list", self.nodes_list)

            elif rand > p:#add edge

                # X_rad = random.uniform(0,1)
                # if X_rad <= p1:
                #     X_t = 2
                # else:
                #     X_t = 3
                X_t=random_pick(edge_length, p_norm)
                self.X.append(X_t)
                #print("X_t",X_t)
                P_id=random_pick(self.some_list, P)
                P_choice=np.array(com_all[P_id])
                #print("P_id", P_id,P_choice)
                C_name_choice=C_name[P_choice==1]
                #print("C_name_choice",len(C_name_choice),C_name_choice)

                split=list_split_name(C_name_choice,X_t,len(C_name_choice))
                #print("split",split)


                new_edge = []
                for C_i in split.keys():
                    p_C_i=calculate_prob(self.degree[C_i], gamma=gamma)
                    #print("p_C_i",p_C_i)
                    #list(random.choices(list(p_C_i.keys()), weights=list(p_C_i.values()), k=X_t[C_i_id]),)
                    z=split[C_i]
                    new_nodes=list(np.random.choice(list(p_C_i.keys()),p=list(p_C_i.values()),size=z,replace=False))
                    new_edge.extend(new_nodes)
                    self.degree[C_i] = self.degree[C_i] + Counter(new_nodes)
                    #print("self.degree[C_i]", new_nodes,self.degree[C_i])
                self.edge[t] = new_edge
                self.node_inedge |= set(new_edge)
                if  len(C_name_choice)>=2:
                    self.edge_between[t]=new_edge

            t+=1
        #print("t",t)
        return self.edge

    def hypergraph_family(self):
        edge=[]
        for nodelist in self.community.values():
            new_edges=generate_edage(nodelist, possion_mu=3)#possion_mu=2.3
            edge.extend(new_edges)
        self.edge[0]=edge
        return self.edge

    def hypergraph_PA_family(self,C0, p, M, p1, P, Node_num_T):
        #self.hypergraph_PA(C0, p, M, p1, P, update_iter, Node_num_T)
        self.hypergraph_PA(C0, p, M, p1, P, Node_num_T)
        self.hypergraph_family()
        edge=copy.deepcopy(self.edge)
        G_degree=copy.deepcopy(self.degree)
        G_community=copy.deepcopy(self.community)
        G_edge_between=copy.deepcopy(self.edge_between)
        #degree=copy.deepcopy(self.degree)
        return edge,G_degree,G_community,G_edge_between





if __name__ == "__main__":
    n_value=2
    C_num=n_value#Total number of community
    d_num=n_value#2#n_value#The number from different community
    #alpha = 0.2 #
    seed_all=[1,2,3,4,5,6,7,8,9,10]
    for seed in seed_all:
        print("seed",seed)
        #alpha_all=[1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05,0.01]
        #alpha_all=[0.09,0.07,0.05,0.03,0.01]
        #alpha_all=[0.009,0.007,0.005,0.003,0.001]
        #alpha_all=[0.08,0.06,0.04,0.008,0.006,0.004,0.002]
        alpha_all=[0.00]
        for alpha in alpha_all:
            #seed=5
            random.seed(seed)
            np.random.seed(seed)
            Node_num_0=20*C_num
            gamma=2
            p=0.5#
            p1=0.5
            #update_iter=10000#
            Node_num_T=1000#
            C_name=np.array(['C%d'%(i+1) for i in range(C_num)])#
            #print(C_name)
            Node=[j for j in range(Node_num_0)]#
            Node_split=list_split(Node,C_num)
            #print("Node_split",Node_split)
        
            C0={}#初始社团
            for key,value in zip(C_name,Node_split):
                C0[key]=value
            #print("C0",C0)
            #M=func1(1,C_num)#
            M=[1/C_num for i in range(C_num)]#
            #print(M)
        
            com_all=combinations_fun_all(C_name, d_num)#
            #print("com_all",len(com_all),com_all)
            P=P_function(len(com_all), alpha)#
            # X = X_time_sequence(C_num, update_iter,p1=p1)
        
            # G=Hypergraph()#
            # #node_edge=G.hypergraph_PA(C0,p=p,M=M,p1=p1,P=P,update_iter=update_iter,Node_num_T=Node_num_T)
            # node_edge = G.hypergraph_PA(C0, p=p, M=M, p1=p1, P=P, Node_num_T=Node_num_T)
            # print("node_edge",len(node_edge)-1,node_edge)
            # G.community
            # print("G.community",G.community)
        
            G = Hypergraph()
            node_edge,G_degree,G_community,G_edge_between= G.hypergraph_PA_family(C0, p=p, M=M, p1=p1, P=P, Node_num_T=Node_num_T)
            #("node_edge",len(node_edge),node_edge)
            #print("G_community",G_community)
         
            #print("self.edge_between", G_edge_between)
        
            f=open('./dataset5/noloop_network_seed%d_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f.txt'%(seed,Node_num_T, C_num,d_num,Node_num_0,gamma,p,alpha,p1),'w')
            f.write(str(node_edge))
            f.close()
        
            f2 = open('./dataset5/noloop_seed%d_degree_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f.txt'%(seed,Node_num_T, C_num,d_num,Node_num_0,gamma,p,alpha,p1), 'w')
            f2.write(str(G_degree))
            f2.close()
        
            f3 = open('./dataset5/noloop_seed%d_community_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f.txt'%(seed,Node_num_T, C_num,d_num, Node_num_0, gamma, p ,alpha, p1), 'w')
            f3.write(str(G_community))
            f3.close()
        
            f4 = open('./dataset5/noloop_seed%d_network_between_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f.txt'%(seed,Node_num_T, C_num,d_num,Node_num_0,gamma,p,alpha,p1),'w')
            f4.write(str(G_edge_between))
            f4.close()
        
        




