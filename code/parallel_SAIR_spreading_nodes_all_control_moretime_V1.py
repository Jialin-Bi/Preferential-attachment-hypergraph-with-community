import numpy as np
import random
import math
from collections import Counter
import SAIR_on_hypergraph_choice_E_jianhua
import csv
import copy
from parallel import parallel
import time
from spreading_method import fixed_length_edges,delete_fixed_edges,delete_fixed_length_edges,delete_edges,direct_edge_comnine

seed_all=[3, 5, 6, 14, 526, 18, 21, 533, 536, 544, 40, 59, 573, 81, 83, 598, 600, 95, 611, 625, 117, 119, 638, 648, 136, 652,
     655, 145, 660, 151, 157, 675, 678, 173, 689, 177, 179, 693, 193, 709, 713, 203, 204, 719, 209, 723, 734, 235, 240,
     241, 758, 764, 783, 786, 274, 281, 799, 287, 809, 299, 816, 817, 820, 314, 832, 324, 332, 342, 344, 345, 346, 347,
     350, 869, 875, 892, 390, 905, 910, 409, 412, 417, 422, 423, 426, 939, 945, 450, 966, 968, 970, 973, 466, 491, 503,1,2,4,7,9]
'''
def SAIR_function(x,beta1,beta2,mu,lamda,SAIR,seed):
    SAIR.initial_setup(I_percentage=None, I_num=None, fixed_nodes_to_infect=[x],seed=seed, Gcommunity=None)
    s, e, i, r, slist, elist, ilist, rlist = SAIR.SAIRmodel(t_max=60, beta1=beta1, beta2=beta2,lamda=lamda, mu=mu)
    srho, srho_t = SAIR.get_stationary_rho_s(normed=True, last_k_values=1)
    rho, rho_t = SAIR.get_stationary_rho_r(normed=True, last_k_values=1)
    erho, erho_t = SAIR.get_stationary_rho_e(normed=True, last_k_values=1)
    irho, irho_t = SAIR.get_stationary_rho_i(normed=True, last_k_values=1)


    return (x,srho, srho_t,erho, erho_t,irho, irho_t ,rho, rho_t)
'''
def SAIR_function(x,beta1,beta2,mu,lamda,SAIR,I_num,Gcommunityi):
    infected_this_setup =SAIR.initial_setup(I_percentage=None, I_num=I_num, seed=x, Gcommunity=Gcommunityi)
    #print("infected_this_setup", infected_this_setup)
    s, e, i, r, slist, elist, ilist, rlist = SAIR.SAIRmodel(t_max=60, beta1=beta1, beta2=beta2,lamda=lamda, mu=mu)
    srho, srho_t = SAIR.get_stationary_rho_s(normed=True, last_k_values=1)
    rho, rho_t = SAIR.get_stationary_rho_r(normed=True, last_k_values=1)
    erho, erho_t = SAIR.get_stationary_rho_e(normed=True, last_k_values=1)
    irho, irho_t = SAIR.get_stationary_rho_i(normed=True, last_k_values=1)


    return (x,srho, srho_t,erho, erho_t,irho, irho_t ,rho, rho_t)

if __name__ == '__main__':
    random.seed(5)
    np.random.seed(5)

    time_start = time.time()
    C_num=2#
    d_num=2#
    Node_num_0=20*C_num
    gamma=2
    p=0.5#
    alpha = 0.2 #
    p1=0.5
    Node_num_T=1000

    #读取
    # beta1 = 0.01
    # beta2 = 0.05
    beta1 = 0.05
    beta2 = 0.1
    lamda=0.1 #
    mu = 1 #
    fixed_nodes_to_infect=None#
    I_num=10#

    method=0 #Randomly delete edges within the community
    #method=1#elete the edges between the community
    #method = 2  # Randomly delete edges 5, max
    ##method=3#Randomly delete edges
    #method=4#Randomly delete edge 2, max

    C_name = np.array(['C%d' % (i + 1) for i in range(C_num)])  # 
    #result='result3_more'
    result='result52'

    f=open('./dataset5/noloop_network_seed5_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f.txt'%(Node_num_T, C_num,d_num,Node_num_0,gamma,p,alpha,p1),'r')
    a=f.read()
    direct_edges=eval(a)
    print("dict_data",direct_edges)
    edges1=list(direct_edges.values())
    print("edges1", edges1)
    edges=edges1[0]+edges1[1:]
    print("edges", edges)
    f.close()

    f2 = open('./dataset5/noloop_seed5_degree_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f.txt'%(Node_num_T, C_num,d_num,Node_num_0,gamma,p,alpha,p1), 'r')
    b = f2.read()
    Gdegree = eval(b)
    f2.close()

    f3 = open('./dataset5/noloop_seed5_community_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f.txt'%(Node_num_T, C_num,d_num, Node_num_0, gamma, p ,alpha, p1), 'r')
    c = f3.read()
    Gcommunity = eval(c)
    f3.close()

    f4 = open('./dataset5/noloop_seed5_network_between_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f.txt'%(Node_num_T, C_num,d_num,Node_num_0,gamma,p,alpha,p1), 'r')
    d = f4.read()
    edges_between = eval(d)
    f4.close()
    #print('edges_between',len(edges_between),edges_between)

    nodes=[]
    for node in Gcommunity.values():
        nodes.extend(node)
    nodes=set(nodes)
    #print("nodes",len(nodes),nodes)
    
    
    if Node_num_T != None:
        #
        resultname1 = (
        './'+result+'/choice_noloop_network_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f I_num%d beta1%0.2f beta2%0.2f lamda%0.2f mu%0.2f SAIR_s_range_time_method%d repeat more.csv' % (
        Node_num_T, C_num, d_num, Node_num_0, gamma, p, alpha, p1, I_num, beta1,beta2,lamda, mu,method))
        File1 = open(resultname1, 'w+', newline='')
        Writing_File1 = csv.writer(File1)

        resultname2 = (
        './'+result+'/choice_noloop_network_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f I_num%d beta1%0.2f beta2%0.2f lamda%0.2f mu%0.2f SAIR_r_range_time_method%d repeat more.csv' % (
        Node_num_T, C_num, d_num, Node_num_0, gamma, p, alpha, p1, I_num, beta1,beta2,lamda, mu,method))
        File2 = open(resultname2, 'w+', newline='')
        Writing_File2 = csv.writer(File2)

        resultname3 = (
        './'+result+'/choice_noloop_network_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f I_num%d beta1%0.2f beta2%0.2f lamda%0.2f mu%0.2f SAIR_e_range_time_method%d repeat more.csv' % (
        Node_num_T, C_num, d_num, Node_num_0, gamma, p, alpha, p1, I_num, beta1, beta2, lamda, mu,method))
        File3 = open(resultname3, 'w+', newline='')
        Writing_File3 = csv.writer(File3)

        resultname4 = (
        './'+result+'/choice_noloop_network_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f I_num%d beta1%0.2f beta2%0.2f lamda%0.2f mu%0.2f SAIR_i_range_time_method%d repeat more.csv' % (
        Node_num_T, C_num, d_num, Node_num_0, gamma, p, alpha, p1, I_num, beta1, beta2, lamda, mu,method))
        File4 = open(resultname4, 'w+', newline='')
        Writing_File4 = csv.writer(File4)

        resultname5 = (
        './'+result+'/choice_noloop_network_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f I_num%d beta1%0.2f beta2%0.2f lamda%0.2f mu%0.2f SAIR_r_range_method%d repeat more.csv' % (
        Node_num_T, C_num, d_num, Node_num_0, gamma, p, alpha, p1, I_num, beta1, beta2, lamda, mu,method))
        File5 = open(resultname5, 'w+', newline='')
        Writing_File5 = csv.writer(File5)
        
  
    flag = 1
    degree=0
    for key, value in Gdegree.items():
        if flag == 1:
            degree = value
            flag = 0
        else:
            degree += value
    #print("degree", degree)
    degree_count = Counter(degree.values())
    #print(Counter(degree.values()))
    degree_values = sorted(degree_count.keys())
    #print("degree_values", degree_values)
    degree_values_time = []
    for key in degree_values:
        degree_values_time.append(degree_count[key])
    #print(degree_values_time)

    degree_sorted = sorted(degree.items(), key=lambda x: x[1], reverse=True)
    n_node_sorted=[]
    n_degree_sorted = []
    n_node_sorted_2=[]
    n_degree_sorted_2 = []
    for n in degree_sorted:
        n_node=n[0]
        n_degree=n[1]
        n_node_sorted.append(n_node)
        n_degree_sorted.append(n_degree)
        if n_degree >1:
            n_node_sorted_2.append(n_node)
            n_degree_sorted_2.append(n_degree)


    direct_edges_one = direct_edge_comnine(direct_edges=direct_edges)#
    #print("direct_edges_one", len(direct_edges_one))
    direct_edges_one_value=direct_edges_one.values()



    edge_del2 = [len(edges_between), len(edges_between) - 1, len(edges_between) - 2, len(edges_between) - 3, len(edges_between) - 4,len(edges_between) - 5]
    edge_del1 = [int(float(i)/100*len(edges_between)) for i in range(100, 0, -5)]
    edge_all = edge_del2 + edge_del1[1:] + list([0])#

    edge_have1=[float(i) for i in range(0, 100, 5)]#
    edge_have=[0,100/len(edges_between),200/len(edges_between),300/len(edges_between),400/len(edges_between),500/len(edges_between)]+edge_have1[1:]+ list([100])#


    edge_all_flag=1
    for edges_num in edge_all:
        for times in range(5):#
            print("edges_num",edges_num)
            #
            if method==1:
                # 
                direct_edges_del,edges_del_key = delete_fixed_edges(direct_edges_one, direct_edges_inlength=edges_between, edges_num=edges_num,seed=times)
                #print("direct_edges", len(direct_edges_del.keys()),direct_edges_del.keys())
    
            if method==2:
                length_min = 5
                # direct_edges 
                #direct_edges_inlength = fixed_length_edges(direct_edges_one,length_min=10,length_max=10 )
                direct_edges_inlength = fixed_length_edges(direct_edges_one, length_min=length_min)
                #print("direct_edges_inlength", direct_edges_inlength)
                direct_edges_del,edges_del_key = delete_fixed_edges(direct_edges_one, direct_edges_inlength, edges_num=edges_num,seed=times)
                #print("direct_edges_del", len(direct_edges_del))
            if method==4:
                length_min = 2
                direct_edges_inlength = fixed_length_edges(direct_edges_one, length_min=length_min)
                #print("direct_edges_inlength", direct_edges_inlength)
                direct_edges_del,edges_del_key = delete_fixed_edges(direct_edges_one, direct_edges_inlength, edges_num=edges_num,seed=times)
                #print("direct_edges_del", len(direct_edges_del))
            if method==3:
                #
                direct_edges_del,edges_del_key = delete_fixed_edges(direct_edges_one, direct_edges_inlength=direct_edges_one, edges_num=edges_num,seed=times)
            if method==5:  
                edges_incom=dict()
                for key in direct_edges_one.keys():
                    if key not in edges_between.keys():
                        edges_incom[key]=direct_edges_one[key]
                direct_edges_del,edges_del_key = delete_fixed_edges(direct_edges_one, direct_edges_inlength=edges_incom, edges_num=edges_num,seed=times)        
                    
                #
                
                
            nodes_del=[]
            edges=list(direct_edges_del.values())#
            for key_id in edges_del_key:#edges_del_key
                nodes_del.extend(direct_edges_one[key_id])
            nodes_del_set=set(nodes_del)
            nodes_del_len=len(nodes_del_set)
            nodes_del_all_len=len(nodes_del)
            
            print(nodes_del_len)

    
    
    
            SAIR = SAIR_on_hypergraph_choice_E_jianhua.SAIRModel(edges,nodes)
            #for seed in seed_all[:5]:
                #print("seed",seed)
            for ci in C_name:
                print("seed",ci)
                #
                node_list=[]
                srhos_list = []
                erhos_list = []
                rhos_list=[]
                irhos_list = []
                #RList= parallel(SAIR_function, n_node_sorted, beta1=beta1,beta2=beta2, mu=mu, SAIR=SAIR, lamda=lamda,seed=seed)
                #print("RList",RList)
                RList= parallel(SAIR_function, seed_all[:100], beta1=beta1,beta2=beta2, mu=mu, SAIR=SAIR, lamda=lamda,I_num=I_num,Gcommunityi=Gcommunity[ci])
                for R in RList:
                    Writing_File1.writerow([R[0]] + list(R[2]))
                    Writing_File3.writerow([R[0]] + list(R[4]))
                    Writing_File4.writerow([R[0]] + list(R[6]))
                    Writing_File2.writerow([R[0]] + list(R[8]))
    
    
                    node_list.append(R[0])
                    srhos_list.append(R[1])
                    erhos_list.append(R[3])
                    irhos_list.append(R[5])
                    rhos_list.append(R[7])
    
    
                if edge_all_flag==1:
                    Writing_File5.writerow(['nodes'] + node_list)
                    edge_all_flag=0
                #Writing_File5.writerow(['%0.3f_%d_%d' %(edge_have[edge_all.index(edges_num)],nodes_del_len,nodes_del_all_len)] + rhos_list)
                Writing_File5.writerow(['%0.3f' %(edge_have[edge_all.index(edges_num)])] + rhos_list)


    File1.close()
    File2.close()
    File3.close()
    File4.close()
    File5.close()
    
    time_end=time.time()
    print("totally time",time_end-time_start)

