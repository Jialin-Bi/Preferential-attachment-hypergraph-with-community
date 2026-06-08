import numpy as np
import random
import math
from collections import Counter
import matplotlib.pyplot as plt
import copy

class SAIRModel():
    def __init__(self,edge,nodes):
        self.nodes = set(nodes)
        self.edges=edge
        self.N=len(self.nodes)#
        print("self.N",self.N)


    def initial_setup(self,I_num=None,I_percentage=None,fixed_nodes_to_infect=None,seed=0,Gcommunity=None):
        self.sAgentSet = set()
        self.eAgentSet = set()
        self.iAgentSet = set()
        self.rAgentSet = set()
        self.sList = []
        self.eList = []
        self.iList = []
        self.rList = []
        self.t = 0

        random.seed(seed)
        np.random.seed(seed)
        infected_this_setup = []

        if I_num != None:#
            self.I = I_num
            self.sAgentSet = copy.deepcopy(self.nodes)
            Gcommunity_ci = copy.deepcopy(Gcommunity)
            to_infect=np.random.choice(list(Gcommunity_ci),size= I_num,replace=False)
            infected_this_setup.extend(to_infect)
            self.sAgentSet -=set(infected_this_setup)
            self.eAgentSet =set(infected_this_setup)


        if fixed_nodes_to_infect != None: #
            fixed_nodes_to_infect=list(fixed_nodes_to_infect)
            self.sAgentSet = copy.deepcopy(self.nodes)
            self.sAgentSet -=set(fixed_nodes_to_infect)
            self.eAgentSet =set(fixed_nodes_to_infect)
            infected_this_setup = fixed_nodes_to_infect
            # for to_infect in fixed_nodes_to_infect:
            #     self.infectAgent(to_infect)
            #     infected_this_setup.append(to_infect)
            #return infected_this_setup

        elif I_percentage != None:
            # infect nodes
            C_N=len(Gcommunity)
            self.I = int(I_percentage * self.N/100)
            #self.I = int(I_percentage * C_N /100)
            self.sAgentSet = copy.deepcopy(self.nodes)
            if fixed_nodes_to_infect == None:  # the first time I create the model (the instance __init__)
                Gcommunity_ci = copy.deepcopy(Gcommunity)
                for ite in range(self.I):  # we will infect E agents
                    # select one to infect among the supsceptibles
                    to_infect = random.choice(list(Gcommunity_ci))
                    Gcommunity_ci.remove(to_infect)
                    self.infectAgent(to_infect)
                    infected_this_setup.append(to_infect)
            else:  # I already have run the model and this is not the first run, I want to infect the same nodes

                for to_infect in fixed_nodes_to_infect:
                    self.infectAgent(to_infect)
                    infected_this_setup.append(to_infect)

        self.sList.append(len(self.sAgentSet))
        self.eList.append(len(self.eAgentSet))
        self.iList.append(0)
        self.rList.append(0)

        return infected_this_setup


    def infectAgent(self, agent):
        self.eAgentSet.add(agent)
        self.sAgentSet.remove(agent)
        return 1


    def SAIRmodel(self,t_max, beta1, beta2, lamda, mu,theta1=None):
        self.t_max = t_max
        # while len(self.iAgentSet) > 0 and len(self.sAgentSet) !=0 and self.t < self.t_max:
        while len(self.iAgentSet) > 0 or len(self.eAgentSet) > 0:
            #print("t",self.t)
            newElist = set()#
            newIlist = set()#
            newRlist = set()#

           
            #
            for e in self.edges:
            # 
                ecopy = e.copy()
                e_set = set(ecopy)
                e_s = e_set.intersection(self.sAgentSet)  # 
                if len(e_s)>0 and len(e_s)<len(e_set):
                    e_i = e_set.intersection(self.iAgentSet)#
                    e_e = e_set.intersection(self.eAgentSet)# 
                    #e_i=[node for node in e if node in self.iAgentSet ] #
                    #e_e = [node for node in e if node in self.eAgentSet]  # 
                    #beta = (beta1*np.sqrt(len(e_e))+beta2*np.sqrt(len(e_i)))#
                    beta = (beta1 * len(e_e) + beta2 * len(e_i))  # 
                    #beta = (beta1 * np.log(1+len(e_e)) + beta2 * np.log(1+len(e_i)))  # 
                    for n1 in e_s:
                        if(random.random() <= beta):
                            newElist.add(n1)
                            self.sAgentSet.remove(n1)


            #print("newElist",newElist)
            # 
            if len(self.eAgentSet)!=0:
                for ne in self.eAgentSet:
                    if(random.random()<= lamda):#由E转I
                        newIlist.add(ne)


            #
            if len(self.iAgentSet) != 0:
                for ni in self.iAgentSet:
                    if (random.random() <= mu):
                        newRlist.add(ni)

            #
            #self.sAgentSet -= newElist
            self.eAgentSet |= newElist
            self.eAgentSet -= newIlist
            self.iAgentSet |= newIlist
            self.iAgentSet -= newRlist
            self.rAgentSet |= newRlist
            st=len(self.sAgentSet)
            et=len(self.eAgentSet)
            it=len(self.iAgentSet)
            rt=len(self.rAgentSet)
            self.sList.append(st)
            self.eList.append(et)
            self.iList.append(it)
            self.rList.append(rt)
            self.t += 1
        s_range = len(self.sAgentSet)
        e_range = len(self.eAgentSet)
        i_range = len(self.iAgentSet) #
        r_range = len(self.rAgentSet)
        return s_range,e_range,i_range,r_range,self.sList,self.eList,self.iList,self.rList


    def get_stationary_rho_s(self, last_k_values,normed=True,):
        s = self.sList
        if len(s)==0:
            return 0
        if normed:
            s = (1.*np.array(s))/self.N
        if s[-1]==1:
            return 1,s
        elif s[-1]==0:
            return 0,s
        else:
            avg_s = np.mean(s[-last_k_values:])
            avg_s = np.nan_to_num(avg_s) #if there are no infected left nan->0
            s=np.nan_to_num(s)
            return avg_s,s#


    def get_stationary_rho_e(self, last_k_values,normed=True,):
        e = self.eList
        if len(e)==0:
            return 0
        if normed:
            e = (1.*np.array(e))/self.N
        if e[-1]==1:
            return 1,e
        elif e[-1]==0:
            return 0,e
        else:
            avg_e = np.mean(e[-last_k_values:])
            #print("i[-last_k_values:]",i[-last_k_values:])
            #avg_e = np.nan_to_num(avg_e) #if there are no infected left nan->0
            #e=np.nan_to_num(e)
            return avg_e,e


    def get_stationary_rho_i(self, last_k_values,normed=True,):
        i = self.iList
        if len(i)==0:
            return 0
        if normed:
            i = (1.*np.array(i))/self.N
        if i[-1]==1:
            return 1,i
        elif i[-1]==0:
            return 0,i
        else:
            avg_i = np.mean(i[-last_k_values:])
            #print("i[-last_k_values:]",i[-last_k_values:])
            #avg_i = np.nan_to_num(avg_i) #if there are no infected left nan->0
            #i=np.nan_to_num(i)
            return avg_i,i

    def get_stationary_rho_r(self, last_k_values,normed=True,):
        r = self.rList
        if len(r)==0:
            return 0
        if normed:
            r = (1.*np.array(r))/self.N
        if r[-1]==1:
            return 1,r
        elif r[-1]==0:
            return 0,r
        else:
            avg_r = np.mean(r[-last_k_values:])
            #avg_r = np.nan_to_num(avg_r) #if there are no infected left nan->0
            #r=np.nan_to_num(r)
            return avg_r,r

if __name__ == '__main__':
    C_num = 5  # 