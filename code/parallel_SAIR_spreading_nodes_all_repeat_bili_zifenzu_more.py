import numpy as np
import random
import math
from collections import Counter
import SAIR_on_hypergraph_choice_E_jianhua_zifenzu
import csv
import copy
from parallel import parallel
import time
from spreading_method import fixed_length_edges,delete_fixed_edges,delete_fixed_length_edges,delete_edges,direct_edge_comnine

seed_all=[3, 5, 6, 14, 526, 18, 21, 533, 536, 544, 40, 59, 573, 81, 83, 598, 600, 95, 611, 625, 117, 119, 638, 648, 136, 652,
     655, 145, 660, 151, 157, 675, 678, 173, 689, 177, 179, 693, 193, 709, 713, 203, 204, 719, 209, 723, 734, 235, 240,
     241, 758, 764, 783, 786, 274, 281, 799, 287, 809, 299, 816, 817, 820, 314, 832, 324, 332, 342, 344, 345, 346, 347,
     350, 869, 875, 892, 390, 905, 910, 409, 412, 417, 422, 423, 426, 939, 945, 450, 966, 968, 970, 973, 466, 491, 503,1,2,4,7,9]



def SAIR_function(x,beta1,beta2,mu,lamda,SAIR,I_num,Gcommunity_ci,Gcommunity):
    infected_this_setup =SAIR.initial_setup(I_percentage=None, I_num=I_num, seed=x, Gcommunity_ci=Gcommunity_ci,Gcommunity=Gcommunity)
    #print("infected_this_setup", infected_this_setup)
    s, e, i, r, slist, elist, ilist, rlist,slists, elists, ilists, rlists = SAIR.SAIRmodel(t_max=60, beta1=beta1, beta2=beta2,lamda=lamda, mu=mu,Gcommunity=Gcommunity)
    srho, srho_t = SAIR.get_stationary_rho_s(normed=True, last_k_values=1)
    rho, rho_t = SAIR.get_stationary_rho_r(normed=True, last_k_values=1)
    erho, erho_t = SAIR.get_stationary_rho_e(normed=True, last_k_values=1)
    irho, irho_t = SAIR.get_stationary_rho_i(normed=True, last_k_values=1)
    

    return (x,srho, srho_t,erho, erho_t,irho, irho_t ,rho, rho_t,slists, elists, ilists, rlists)
    #return (x,srho, srho_t,erho, erho_t,irho, irho_t ,rho, rho_t)
# def SAIR_function(x,beta1,beta2,mu,lamda,SAIR,seed):
#     SAIR.initial_setup(I_percentage=None, I_num=None, fixed_nodes_to_infect=[x],seed=seed, Gcommunity=None)
#     s, e, i, r, slist, elist, ilist, rlist = SAIR.SAIRmodel(t_max=60, beta1=beta1, beta2=beta2,lamda=lamda, mu=mu)
#     srho, srho_t = SAIR.get_stationary_rho_s(normed=True, last_k_values=1)
#     rho, rho_t = SAIR.get_stationary_rho_r(normed=True, last_k_values=1)
#     erho, erho_t = SAIR.get_stationary_rho_e(normed=True, last_k_values=1)
#     irho, irho_t = SAIR.get_stationary_rho_i(normed=True, last_k_values=1)
#
#
#     return (x,srho, srho_t,erho, erho_t,irho, irho_t ,rho, rho_t)

if __name__ == '__main__':
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    random.seed(5)
    np.random.seed(5)

    time_start = time.time()
    C_num=10#
    d_num=10#
    Node_num_0=20*C_num
    #update_iter=10000
    gamma=2
    p=0.5#
    alpha = 0.2#
    p1=0.5
    Node_num_T=1000

    #read data
    beta=0.005
    lamda=0.1 #
    mu =1 #
    #fixed_nodes_to_infect=True#
    I_num=10#
    C_name = np.array(['C%d' % (i + 1) for i in range(C_num)])  # 
    '''
    #alpha_all=[1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01]
    #alpha_all=[0.09,0.07,0.05,0.03,0.01]
    #alpha_all=[0.009,0.007,0.005,0.003,0.001]
    #alpha_all=[0.08,0.06,0.04,0.008,0.006,0.004,0.002]
    #alpha_all=[1.0,0.8,0.6,0.4,0.2,0.1,0.08,0.06,0.04,0.02,0.01,0.008,0.006,0.004,0.002,0.001]
    #alpha_all=[0.09,0.07,0.05,0.03,0.3,0.5,0.7,0.9]
    alpha_all=[0.009,0.007,0.005,0.003]
    print("alpha_all",alpha_all)
    for alpha in alpha_all:#
    '''
    alpha_all=[1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05,0.01]
    for alpha in alpha_all:#
    
    #beta_all=[0.001,0.003,0.005,0.007,0.009,0.01,0.03,0.05,0.07,0.09,0.10,0.13,0.15,0.17,0.19,0.2]
    #for beta in beta_all:#
    
    
        beta1 = beta
        beta2 = beta
        print("beta",beta)

        f=open('./dataset5/noloop_network_seed5_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f.txt'%(Node_num_T, C_num,d_num,Node_num_0,gamma,p,alpha,p1),'r')
        a=f.read()
        direct_edges=eval(a)
        #print("dict_data",direct_edges)
        edges1=list(direct_edges.values())
        #print("edges1", edges1)
        edges=edges1[0]+edges1[1:]
        #print("edges", edges)
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
    
    
    
        if I_num != None:
            #save
    
            resultname1 = (
            './result51/choice_noloop_network_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f I_num%d beta1%0.3f beta2%0.3f lamda%0.1f mu%0.1f SAIR_s_range_time.csv' % (
            Node_num_T, C_num, d_num, Node_num_0, gamma, p, alpha, p1, I_num, beta1, beta2,lamda, mu))
            File1 = open(resultname1, 'w+', newline='')
            Writing_File1 = csv.writer(File1)
    
            resultname2 = (
            './result51/choice_noloop_network_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f I_num%d beta1%0.3f beta2%0.3f lamda%0.1f mu%0.1f SAIR_r_range_time.csv' % (
            Node_num_T, C_num, d_num, Node_num_0, gamma, p, alpha, p1, I_num, beta1,beta2,lamda, mu))
            File2 = open(resultname2, 'w+', newline='')
            Writing_File2 = csv.writer(File2)
    
            resultname3 = (
            './result51/choice_noloop_network_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f I_num%d beta1%0.3f beta2%0.3f lamda%0.1f mu%0.1f SAIR_e_range_time.csv' % (
            Node_num_T, C_num, d_num, Node_num_0, gamma, p, alpha, p1, I_num, beta1, beta2, lamda, mu))
            File3 = open(resultname3, 'w+', newline='')
            Writing_File3 = csv.writer(File3)
    
            resultname4 = (
            './result51/choice_noloop_network_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f I_num%d beta1%0.3f beta2%0.3f lamda%0.1f mu%0.1f SAIR_i_range_time.csv' % (
            Node_num_T, C_num, d_num, Node_num_0, gamma, p, alpha, p1, I_num, beta1, beta2, lamda, mu))
            File4 = open(resultname4, 'w+', newline='')
            Writing_File4 = csv.writer(File4)
    
            resultname5 = (
            './result51/choice_noloop_network_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f I_num%d beta1%0.3f beta2%0.3f lamda%0.1f mu%0.1f SAIR_r_range.csv' % (
            Node_num_T, C_num, d_num, Node_num_0, gamma, p, alpha, p1, I_num, beta1, beta2, lamda, mu))
            File5 = open(resultname5, 'w+', newline='')
            Writing_File5 = csv.writer(File5)
            
            resultname6 = (
            './result51/choice_noloop_network_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f I_num%d beta1%0.3f beta2%0.3f lamda%0.1f mu%0.1f SAIR_s_range_time_fenzu.csv' % (
            Node_num_T, C_num, d_num, Node_num_0, gamma, p, alpha, p1, I_num, beta1, beta2,lamda, mu))
            File6 = open(resultname6, 'w+', newline='')
            Writing_File6 = csv.writer(File6)
    
            resultname7 = (
            './result51/choice_noloop_network_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f I_num%d beta1%0.3f beta2%0.3f lamda%0.1f mu%0.1f SAIR_e_range_time_fenzu.csv' % (
            Node_num_T, C_num, d_num, Node_num_0, gamma, p, alpha, p1, I_num, beta1,beta2,lamda, mu))
            File7 = open(resultname7, 'w+', newline='')
            Writing_File7 = csv.writer(File7)
    
            resultname8 = (
            './result51/choice_noloop_network_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f I_num%d beta1%0.3f beta2%0.3f lamda%0.1f mu%0.1f SAIR_i_range_time_fenzu.csv' % (
            Node_num_T, C_num, d_num, Node_num_0, gamma, p, alpha, p1, I_num, beta1, beta2, lamda, mu))
            File8 = open(resultname8, 'w+', newline='')
            Writing_File8 = csv.writer(File8)
    
            resultname9 = (
            './result51/choice_noloop_network_node%d C_num%d d_num%d Node%d gamm%d p%0.1f alpha%0.4f p1%0.1f I_num%d beta1%0.3f beta2%0.3f lamda%0.1f mu%0.1f SAIR_r_range_time_fenzu.csv' % (
            Node_num_T, C_num, d_num, Node_num_0, gamma, p, alpha, p1, I_num, beta1, beta2, lamda, mu))
            File9 = open(resultname9, 'w+', newline='')
            Writing_File9 = csv.writer(File9)
    
    
    
        flag = 1
        degree=0
        for key, value in Gdegree.items():
            if flag == 1:
                degree = value
                flag = 0
            else:
                degree += value
        degree_count = Counter(degree.values())
        degree_values = sorted(degree_count.keys())
        degree_values_time = []
        for key in degree_values:
            degree_values_time.append(degree_count[key])
    
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

        direct_edges_one_value=direct_edges_one.values()
    
        edge_all_flag = 1
    
        edges=list(direct_edges_one_value)

    
    
        SAIR = SAIR_on_hypergraph_choice_E_jianhua_zifenzu.SAIRModel(edges,nodes)

    
        for key, value in Gcommunity.items():
            Gcommunity[key]=set(Gcommunity[key])
            Gcommunity[key] &= copy.deepcopy(SAIR.nodes)
        #for ci in C_name[:5]:
        for ci in C_name:
            print("seed",ci)
            #for ci in C_name:
            node_list=[]
            srhos_list = []
            erhos_list = []
            rhos_list=[]
            irhos_list = []
            RList= parallel(SAIR_function, seed_all[:100], beta1=beta1,beta2=beta2, mu=mu, SAIR=SAIR, lamda=lamda,I_num=I_num,Gcommunity_ci=Gcommunity[ci],Gcommunity=Gcommunity)
    
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
                
                #
                for row in R[9].values():
                    Writing_File6.writerow([R[0]] + list(row))
                for row in R[10].values():
                    Writing_File7.writerow([R[0]] + list(row))
                for row in R[11].values():
                    Writing_File8.writerow([R[0]] + list(row))
                for row in R[12].values():
                    Writing_File9.writerow([R[0]] + list(row))
    
    
            #Writing_File5.writerow(['Community=%s' % ci] + rhos_list)
            if edge_all_flag==1:
                Writing_File5.writerow(['nodes'] + node_list)
                edge_all_flag=0
            #Writing_File5.writerow(['nodes'] + node_list)
            Writing_File5.writerow(['0'] + rhos_list)
            
    
                
    
        File1.close()
        File2.close()
        File3.close()
        File4.close()
        File5.close()
        File6.close()
        File7.close()
        File8.close()
        File9.close()
        time_end=time.time()
        print("totally time",time_end-time_start)

