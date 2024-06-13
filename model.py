import numpy as np
import matplotlib.pyplot as plt
import matplotlib
%matplotlib inline
import networkx as nx
import itertools
import random
import warnings
import seaborn as sns
from scipy.stats import entropy
import pandas as pd
import pickle

fit_inc=1.05 #Increase in fitness for next intermediate level
fit_cross=2 #Increase in fitness for next level (i.e. cross)

### Generate Fitness Vector
rrr=50

starting_fitness=np.array([1,1,1,1])
fitnessvec=[]
fitnessvec.append(starting_fitness)
fitnessvec.append(starting_fitness*fit_inc)
fitnessvec.append(starting_fitness*fit_inc*fit_cross)
for i in range(rrr-1):
    fitnessvec.append(fitnessvec[len(fitnessvec)-1]*fit_inc)
    fitnessvec.append(fitnessvec[len(fitnessvec)-1]*fit_inc)
    fitnessvec.append(fitnessvec[len(fitnessvec)-1]*fit_cross)
fitnessvec = [item for sublist in fitnessvec for item in sublist]

### Generate Discoveries Dictionary
discoveries_dict={}
for k in range(rrr):
    c=16*k
    #A
    discoveries_dict[c,c+4]=c+8
    discoveries_dict[c+4,c]=c+8
    
    discoveries_dict[c+4,c+8]=c+12
    discoveries_dict[c+8,c+4]=c+12
    
    discoveries_dict[c+8,c+12]=c+16
    discoveries_dict[c+12,c+8]=c+16
    

    discoveries_dict[c+16,c+17]=c+20
    discoveries_dict[c+16,c+18]=c+20
    discoveries_dict[c+16,c+19]=c+20
    
    #B
    discoveries_dict[c+1,c+4+1]=c+8+1
    discoveries_dict[c+4+1,c+1]=c+8+1
    
    discoveries_dict[c+4+1,c+8+1]=c+12+1
    discoveries_dict[c+8+1,c+4+1]=c+12+1
    
    discoveries_dict[c+8+1,c+12+1]=c+16+1
    discoveries_dict[c+12+1,c+8+1]=c+16+1
    
    discoveries_dict[c+17,c+16]=c+20+1
    discoveries_dict[c+17,c+18]=c+20+1
    discoveries_dict[c+17,c+19]=c+20+1
    
    #C
    discoveries_dict[c+2,c+4+2]=c+8+2
    discoveries_dict[c+4+2,c+2]=c+8+2
    
    discoveries_dict[c+4+2,c+8+2]=c+12+2
    discoveries_dict[c+8+2,c+4+2]=c+12+2
    
    discoveries_dict[c+8+2,c+12+2]=c+16+2
    discoveries_dict[c+12+2,c+8+2]=c+16+2
    
    discoveries_dict[c+18,c+16]=c+20+2
    discoveries_dict[c+18,c+17]=c+20+2
    discoveries_dict[c+18,c+19]=c+20+2
    
    #D
    discoveries_dict[c+3,c+7]=c+11
    discoveries_dict[c+7,c+3]=c+11
    
    discoveries_dict[c+7,c+11]=c+15
    discoveries_dict[c+11,c+7]=c+15
    
    discoveries_dict[c+11,c+15]=c+19
    discoveries_dict[c+15,c+11]=c+19
    
    discoveries_dict[c+19,c+16]=c+20+3
    discoveries_dict[c+19,c+17]=c+20+3
    discoveries_dict[c+19,c+18]=c+20+3
    
### Generate Level Vectors, Line Vectors, Color Vectors
levelvec=[[0,0,0,0],[0,0,0,0]]
for i in range(rrr):
    for ii in range(3):
        levelvec.append([i+1,i+1,i+1,i+1])
levelvec = [item for sublist in levelvec for item in sublist]


levelvec2=[[0,0,0,0],[0,0,0,0]]
for i in range(rrr*3):
    levelvec2.append([i+1,i+1,i+1,i+1])
levelvec2 = [item for sublist in levelvec2 for item in sublist]


linevec=[['A0.1','B0.1','C0.1','D0.1'],
         ['A0.2','B0.2','C0.2','D0.2']]
for i in range(rrr):
    for ii in range(3):
        linevec.append(['A'+str(i+1)+'.'+str(ii+1),'B'+str(i+1)+'.'+str(ii+1),
                        'C'+str(i+1)+'.'+str(ii+1),'D'+str(i+1)+'.'+str(ii+1)])
linevec = [item for sublist in linevec for item in sublist]


colorvec=[['tab:red','tab:blue','tab:green','tab:olive'],
          ['tab:red','tab:blue','tab:green','tab:olive']]
for i in range(rrr):
    for ii in range(3):
        colorvec.append(['tab:red','tab:blue','tab:green','tab:olive'])
colorvec = [item for sublist in colorvec for item in sublist]


linevec2=[['A','B','C','D'],
          ['A','B','C','D']]
for i in range(rrr):
    for ii in range(3):
        linevec2.append(['A','B','C','D'])
linevec2 = [item for sublist in linevec2 for item in sublist]


###FUNZIONI PER SIMULAZIONE
def generate_network(n,M,parent_r,couple_r,sibl_r):
    """
    This function generate the starting network
    """
    N=M*n
    n_fam=n//5
    if n%5 != 0:
        raise ValueError("Families should be divisible by 5")
    G=nx.Graph()
    communities={}
    families={}
    tools={}
    cult_line={}
    cult_level={}
    a=0
    for i in range(M):
        e=list(itertools.combinations(range(n*i,n*i+n),2)) #Create complete graph in each community
        G.add_edges_from(e)
        for n_f in range(n_fam):
            family=[(a,a+1,couple_r),(a,a+2,parent_r),(a,a+3,parent_r),(a,a+4,parent_r),(a+1,a+2,parent_r),
                    (a+1,a+3,parent_r),(a+1,a+4,parent_r),(a+2,a+3,sibl_r),(a+2,a+4,sibl_r),(a+3,a+4,sibl_r)]
            G.add_weighted_edges_from(family,'relatedness')
            families[a]=n_f+i*n_fam
            families[a+1]=n_f+i*n_fam
            families[a+2]=n_f+i*n_fam
            families[a+3]=n_f+i*n_fam
            families[a+4]=n_f+i*n_fam
            a=a+5
        for ii in range(n):
            communities[n*i+ii]=i #Add community attribute to each node
            tools[n*i+ii]=[0,1,2,3,4,5,6,7]
            cult_line[n*i+ii]=0
    nx.set_node_attributes(G, communities, name='community')
    nx.set_node_attributes(G, families, name='family')
    nx.set_node_attributes(G, tools, name='tools')
    nx.set_node_attributes(G, tools, name='tools_hist')
    nx.set_node_attributes(G, cult_line, name='cult_line')
    return(G)

def map_fitness(tools):
    """
    This function receive tools vector and gives back fitness vector for tools
    """
    tools_fit=[]
    for tool in tools:
        tools_fit.append(fitnessvec[tool])
    return tools_fit
def min_fit_value(tools):
    """
    This function returns minimum fitness from a tool vector
    """
    return min(map_fitness(tools))
def min_fit_idxs(tools):
    """
    This function returns indexes of tools with lowest fitness
    """
    return [i for i in range(len(tools)) if map_fitness(tools)[i] == min_fit_value(tools)]
def min_fit_idx_rand(tools):
    """
    This function selects at random the index of one of the tools with lowest fitness
    """
    return random.choice(min_fit_idxs(tools))
def update_tools(tools,discovery,memory):
    """
    This function checks if discovery already in memory.
    If new discovery has better fitness and replace tool if yes.
    This function uses the memory parameter.
    """
    if discovery not in tools:
        if len(tools)==memory:
            if fitnessvec[discovery] >= min_fit_value(tools):
                tools[min_fit_idx_rand(tools)]=discovery
        if len(tools)<memory:
            tools.append(discovery)
    return sorted(tools)
def update_tools_hist(tools,discovery):
    """
    This function checks if discovery already in memory.
    This function tracks historical discoveries.
    """
    if discovery not in tools:
        tools.append(discovery)
    return sorted(tools)
def distance_dict(G):
    """
    This function generate the distance dictionary (input of dictionary = tuple)
    """
    dic_dist={}
    nodes=G.nodes()
    couples=list(itertools.combinations(G.nodes(),2))
    for cc in couples:
        tool1=set(nx.get_node_attributes(G,'tools')[cc[0]])
        tool2=set(nx.get_node_attributes(G,'tools')[cc[1]])
        dic_dist[cc]=1-len(set.intersection(tool1,tool2))/max(len(tool1),len(tool2))
    return dic_dist
def average_distance(node,list_nodes,dic_dist):
    """
    This function generate the average distance between a node and a list of nodes according to tools
    """
    touple_list=[]
    for i in list_nodes:
        touple_list.append(tuple(sorted([node,i])))
    distance_list=[]
    for i in touple_list:
        distance_list.append(dic_dist[i])
    return sum(distance_list)/len(distance_list)
def successful_share(node1,node2,rel_dic,proximity_p):
    """
    Check if share among two nodes is successful
    """
    edge=tuple(sorted([node1,node2]))
    if random.random() < rel_dic.get(edge,proximity_p):
        outcome = True
    else:
        outcome = False
    return outcome

#Functions for tracking purposes
def map_level(tools):
    """
    This function receive tools vector and gives back level vector for tools (A1.1=1,A1.2=1,...,A.3.3=2)
    """
    tools_fit=[]
    for tool in tools:
        tools_fit.append(levelvec[tool])
    return tools_fit
def map_level2(tools):
    """
    This function receive tools vector and gives back level2 vector for tools (A1.1=1,A1.2=2,...)
    """
    tools_fit=[]
    for tool in tools:
        tools_fit.append(levelvec2[tool])
    return tools_fit
def max_level_value(tools):
    """
    This function returns maximum level from a tool vector
    """
    return max(map_level(tools))
def max_level_idxs(tools):
    """
    This function returns a random index of tools with highest level
    """
    return random.choice([i for i in range(len(tools)) if map_level(tools)[i] == max_level_value(tools)])
def give_line(tools):
    """
    This function returns the cultural line at the highest level (at random if more than 1)
    """
    return linevec[tools[max_level_idxs(tools)]]
def give_level(tools):
    """
    This function returns the level of highest cultural
    """
    return levelvec[tools[max_level_idxs(tools)]]
def max_fit_value(tools):
    """
    This function returns minimum fitness from a tool vector
    """
    return max(map_fitness(tools))
def give_color(tools):
    """
    This function returns color of node depending on highest cultural level
    """
    return colorvec[tools[max_level_idxs(tools)]]
def map_line(tools):
    """
    This function receive tools vector and gives back the cultural line vector
    """
    tools_fit=[]
    for tool in tools:
        tools_fit.append(linevec[tool])
    return tools_fit
def map_line2(tools):
    """
    This function receive tools vector and gives back the cultural line vector
    """
    tools_fit=[]
    for tool in tools:
        tools_fit.append(linevec2[tool])
    return tools_fit
def get_inclusive_fit(G,fitness_i,rel_dic):
    """
    Compute inclusive fitness for all nodes
    """
    inclusive_fitness=[]
    inclusive_fitness_f=[]
    for node in G.nodes():
        temp_fit=fitness_i[node]
        list_family=[x for x,y in G.nodes(data=True)
                               if y['family']==nx.get_node_attributes(G, 'family')[node]]
        list_family.remove(node)
        for fam in list_family:
            pair=tuple(sorted([node,fam]))
            temp_fit=temp_fit+rel_dic[pair]*fitness_i[fam]
        inclusive_fitness.append(temp_fit)
    return inclusive_fitness
def gini(list_of_values):
    """
    Compute Gini coefficient
    """
    sorted_list = sorted(list_of_values)
    height, area = 0, 0
    for value in sorted_list:
        height += value
        area += height - value / 2.
    fair_area = height * len(list_of_values) / 2.
    return (fair_area - area) / fair_area
def count_tot_tools(G):
    tools=list(nx.get_node_attributes(G,'tools').values())
    tools_union=set([])
    for i in tools:
        tools_union=set.union(tools_union,set(i))
    return len(tools_union)
def count_tot_tools_fit(G):
    tools=list(nx.get_node_attributes(G,'tools').values())
    tools_union=set([])
    for i in tools:
        tools_union=set.union(tools_union,set(i))
    return sum(map_fitness(tools_union))
def count_av_tools(G):
    tools=list(nx.get_node_attributes(G,'tools').values())
    n_tools=[]
    for i in tools:
        n_tools.append(len(i))
    return sum(n_tools)/len(n_tools)
def count_av_tools_fit(G):
    tools=list(nx.get_node_attributes(G,'tools').values())
    n_tools=[]
    for i in tools:
        n_tools.append(sum(map_fitness(i)))
    return sum(n_tools)/len(n_tools)
def count_av_frac_tools(G):
    tools=list(nx.get_node_attributes(G,'tools').values())
    total_tools=count_tot_tools(G)
    f_tools=[]
    for i in tools:
        f_tools.append(len(i)/total_tools)
    return sum(f_tools)/len(f_tools)
def count_av_frac_tools_fit(G):
    tools=list(nx.get_node_attributes(G,'tools').values())
    return count_av_tools_fit(G)/count_tot_tools_fit(G)
def gini_av_frac_tools(G):
    tools=list(nx.get_node_attributes(G,'tools').values())
    total_tools=count_tot_tools(G)
    f_tools=[]
    for i in tools:
        f_tools.append(len(i)/total_tools)
    return gini(f_tools)
def fitness_all_tools(G):
    tools=list(nx.get_node_attributes(G,'tools').values())
    av_fit_i=[]
    for i in tools:
        av_fit_i.append(sum(map_fitness(i))/len(i))
    return sum(av_fit_i)/len(av_fit_i)
def count_specialization(G):
    """
    This function compute the vector specialization according to the number in each cultural line
    """
    tools=list(nx.get_node_attributes(G,'tools').values())
    specialization=[]
    for i in tools:
        A=map_line2(i).count('A')
        B=map_line2(i).count('B')
        C=map_line2(i).count('C')
        D=map_line2(i).count('D')
        tot=A+B+C+D
        specialization.append([A/tot,B/tot,C/tot,D/tot])
    return specialization
def count_specialization_fit(G):
    """
    This function compute the specialization vector using the fitness weights
    """
    tools=list(nx.get_node_attributes(G,'tools').values()) 
    specialization_fit=[]
    for i in tools:
        tools_fit=[]
        A=map_line2(i).count('A')
        B=map_line2(i).count('B')
        C=map_line2(i).count('C')
        D=map_line2(i).count('D')
        
        weight_fit=np.array(map_fitness(i))
        tot_fit=np.sum(weight_fit)
        A_fit = np.sum(weight_fit[0:A])
        B_fit = np.sum(weight_fit[A:A+B])
        C_fit = np.sum(weight_fit[A+B:A+B+C])
        D_fit = np.sum(weight_fit[A+B+C:A+B+C+D])
        
        specialization_fit.append([A_fit/tot_fit,B_fit/tot_fit,C_fit/tot_fit,D_fit/tot_fit])
    return specialization_fit
def gini_specialization(G):
    specialization=count_specialization(G)
    gini_i=[]
    for i in specialization:
        gini_i.append(gini(i))
    return sum(gini_i)/len(gini_i)
def entropy_specialization_pdf(G):
    """
    This function returns the entropy values for the specialization (# version)
    """
    specialization=count_specialization(G)
    entropy_i=[]
    for i in specialization:
        entropy_i.append(entropy(i))
    return entropy_i
def entropy_specialization_fit_pdf(G):
    """
    This function returns the entropy values for the specialization (fit version)
    """
    specialization=count_specialization_fit(G)
    entropy_i=[]
    for i in specialization:
        entropy_i.append(entropy(i))
    return entropy_i
def entropy_specialization_camp_pdf(G):
    """
    This function returns the entropy values for camp specialization (# version)
    """
    specialization=count_specialization(G)
    entropy_i=[]
    for i in specialization:
        entropy_i.append(entropy(i))
    n_camps=len(set(nx.get_node_attributes(G,'community').values()))
    n_per_camp=len(G.nodes())//n_camps
    entropy_camp=[]
    for i in range(n_camps):
        entropy_camp.append(np.mean(entropy_i[i*n_per_camp:(i+1)*n_per_camp]))
    return entropy_camp
def entropy_specialization_fit_camp_pdf(G):
    """
    This function returns the entropy values for camp specialization (fit version)
    """
    specialization=count_specialization_fit(G)
    entropy_i=[]
    for i in specialization:
        entropy_i.append(entropy(i))
    n_camps=len(set(nx.get_node_attributes(G,'community').values()))
    n_per_camp=len(G.nodes())//n_camps
    entropy_camp=[]
    for i in range(n_camps):
        entropy_camp.append(np.mean(entropy_i[i*n_per_camp:(i+1)*n_per_camp]))
    return entropy_camp
def threshold_specialization(G,thresh):
    specialization=count_specialization(G)
    count=[]
    for i in specialization:
        if max(i)>=thresh:
            count.append(1)
        else:
            count.append(0)
    return sum(count)/len(count)
def get_family_fit(G,fitness_i):
    """
    Compute family fitness for all nodes
    """
    family_fitness=[]
    for node in G.nodes():
        list_family=[x for x,y in G.nodes(data=True)
                               if y['family']==nx.get_node_attributes(G, 'family')[node]]
        temp_fam=[]
        for fam in list_family:
            temp_fam.append(fitness_i[fam])
        family_fitness.append(max(temp_fam))
    return sum(family_fitness)/len(family_fitness)
def count_tools_hist(G):
    tools=list(nx.get_node_attributes(G,'tools_hist').values())
    tot_tools=set([])
    for i in tools:
        iii=set(i)
        tot_tools.update(iii)
    return len(tot_tools)
def jaccard_dist(list1, list2):
    intersection = len(list(set(list1).intersection(set(list2))))
    union = (len(set(list1)) + len(set(list2))) - intersection
    return 1 - float(intersection) / union
def distance_dict_j(G):
    """
    This function generate the Jaccard distance dictionary (input of dictionary = tuple)
    """
    dic_dist_j={}
    nodes=G.nodes()
    couples=list(itertools.combinations(G.nodes(),2))
    for cc in couples:
        tool1=set(nx.get_node_attributes(G,'tools')[cc[0]])
        tool2=set(nx.get_node_attributes(G,'tools')[cc[1]])
        dic_dist_j[cc]=jaccard_dist(tool1,tool2)
    return dic_dist_j


def modello(n=20, #Number of members per group !!!! MUST BE DIVISIBLE BY 5
            M=15,  #Number of groups
            Epochs_max=100, #Number of epochs
            memory=10000, #Maximum tools in memory
            parent_r=0.5, #Relatedness parent-son
            couple_r=0.5, #Relatedness couple
            sibl_r=0.25, #Relatedness siblings
            proximity_p=0.05, #Sharing probability within same group non-kin
            fit_inc=1.05, #Increase in fitness for next intermediate level
            fit_cross=2, #Increase in fitness for next level (i.e. cross)
            sharing_discoveries=True): #Whether individuals share discoveries (True-False)

    N=M*n
    n_fam=n//5

    ###GENERATE NETWORK
    G=generate_network(n,M,parent_r,couple_r,sibl_r)

    max_fit_t=[]
    gini_fit_t=[]
    max_fit_av_t=[]
    family_fit_av_t=[]
    inclusive_fit_av_t=[]
    prob_out_t=[]
    max_level_t=[]
    gini_level_t=[]
    max_level_av_t=[]
    fraction_knowledge_av_t=[]
    gini_fraction_knowledge_t=[]
    n_tools_av_t=[]
    fitness_all_av_t=[]
    gini_specialization_t=[]
    specialized_individuals_av_t=[]
    pairwise_dist_t=[]
    fraction_knowledge_fit_av_t = []
    

    ### SIMULAZIONE
    stop_sim=0
    matrix_node_tools=[[0]*8+[-1]*(len(fitnessvec)-8)][0] #tracking vector
    rel_dic=nx.get_edge_attributes(G,'relatedness')
    for epoch in range(Epochs_max):
        
        #Compute dictionary of distances
        dic_dist=distance_dict_j(G)
 
        ###TRACKING MEASURES
        prob_out_t_i=[]
        fitness_i=map_fitness([max(p) for p in list(nx.get_node_attributes(G,'tools').values())])
        max_fit_t.append(max(fitness_i))
        gini_fit_t.append(gini(fitness_i))
        max_fit_av_t.append(sum(fitness_i)/N)
        family_fit_av_t.append(get_family_fit(G,fitness_i))
        inclusive_fit_av_t.append(sum(get_inclusive_fit(G,fitness_i,rel_dic))/N)
        level_i=map_level2([max(p) for p in list(nx.get_node_attributes(G,'tools').values())])
        max_level_t.append(max(level_i))
        gini_level_t.append(gini(np.array(level_i)+1))
        max_level_av_t.append(sum(level_i)/N)
        n_tools_av_t.append(count_av_tools(G))
        fraction_knowledge_av_t.append(count_av_frac_tools(G))
        gini_fraction_knowledge_t.append(gini_av_frac_tools(G))
        fitness_all_av_t.append(fitness_all_tools(G))
        gini_specialization_t.append(gini_specialization(G))
        specialized_individuals_av_t.append(threshold_specialization(G,0.5))
        pairwise_dist_t.append(np.mean(list(dic_dist.values())))
        fraction_knowledge_fit_av_t.append(count_av_frac_tools_fit(G)) 

        #### FINISH SIMULATION IF MAXIMUM EPOCH REACHED #####
        if epoch==Epochs_max:
            stop_sim=1
        if (stop_sim==1):
            break
        ####################################################

        ###select the nodes in a random order
        random_order_list=np.random.choice(range(len(G)),size=len(G),replace=False)
        ###from list of size len(G), select len(G) items with no replacement
        for node in random_order_list:
        ###select node

            #compute average distance inside community
            list_insiders=[x for x,y in G.nodes(data=True)
                           if y['community']==nx.get_node_attributes(G, 'community')[node]]
            list_insiders.remove(node)
            av_distance_in=average_distance(node,list_insiders,dic_dist)

            #Compute average distance outside community
            list_outsiders=[x for x,y in G.nodes(data=True)
                            if y['community']!=nx.get_node_attributes(G, 'community')[node]]
            av_distance_out=average_distance(node,list_outsiders,dic_dist)


            #Compute probability of picking node outside community

            prob_out = (av_distance_out - av_distance_in)/2
            if prob_out < 0:
                prob_out = 0
            prob_out_t_i.append(prob_out)

            ### Do we select neighbor in our group?

            if random.random()>prob_out:
                select_in_group=True
                selected_neigh=random.choice(list_insiders)
            else:
                select_in_group=False
                #select node from different community
                selected_neigh=random.choice(list_outsiders)


            #find ingredients for node and for neighbor
            ingredients_node=nx.get_node_attributes(G, 'tools')[node]
            ingredients_neigh=nx.get_node_attributes(G, 'tools')[selected_neigh]

            #select ingredient for node
            node_tool=random.choices(ingredients_node,weights=map_fitness(ingredients_node))[0]

            #select ingredient for neighbor
            neigh_tool=random.choices(ingredients_neigh,weights=map_fitness(ingredients_neigh))[0]

            #Check discovery
            pair=(node_tool,neigh_tool)
            new_tool=discoveries_dict.get(pair,-1) #-1 if not discovered anything

            #Was the discovery successful for the node?
            #First check if unsuccessful
            if new_tool==-1:
                node_success=False
            #Check if already existent tool
            elif new_tool in ingredients_node:
                node_success=False
            #Now check if discovery is better then the worst one in memory
            elif fitnessvec[new_tool]>=min_fit_value(ingredients_node):
                node_success=True
            else:
                node_success=False

            #If successful discovery outside community, create link
            if node_success==True and select_in_group==False:
                G.add_edge(node, selected_neigh)

            #Now check if the discovery was successful for the neighbor
            if new_tool==-1:
                neigh_success=False
            #Check if already existent tool
            elif new_tool in ingredients_neigh:
                neigh_success=False
            #Now check if discovery is better then the worst one in memory
            elif fitnessvec[new_tool]>=min_fit_value(ingredients_neigh):
                neigh_success=True
            else:
                neigh_success=False


            #Update discovery for node and share with some nodes in its community according to relatedness
            if node_success==True:
                G.nodes[node]["tools"]=update_tools(G.nodes[node]["tools"],new_tool,memory)
                G.nodes[node]["tools_hist"]=update_tools_hist(G.nodes[node]["tools_hist"],new_tool)
                if sharing_discoveries==True:
                    for node_comm in list_insiders:
                        if successful_share(node,node_comm,rel_dic,proximity_p) == True:
                            G.nodes[node_comm]["tools"]=update_tools(G.nodes[node_comm]["tools"],new_tool,memory)
                            G.nodes[node_comm]["tools_hist"]=update_tools_hist(G.nodes[node_comm]["tools_hist"],new_tool)


            #Update discovery for neigh and share with some nodes in its community according to relatedness
            if neigh_success==True:
                G.nodes[selected_neigh]["tools"]=update_tools(G.nodes[selected_neigh]["tools"],new_tool,memory)
                G.nodes[selected_neigh]["tools_hist"]=update_tools_hist(G.nodes[selected_neigh]["tools_hist"],new_tool)
                list_insiders_neigh=[x for x,y in G.nodes(data=True)
                                if y['community']==nx.get_node_attributes(G, 'community')[selected_neigh]]
                list_insiders_neigh.remove(selected_neigh)
                if sharing_discoveries==True:
                    for node_comm in list_insiders_neigh:
                        if successful_share(selected_neigh,node_comm,rel_dic,proximity_p) == True:
                            G.nodes[node_comm]["tools"]=update_tools(G.nodes[node_comm]["tools"],new_tool,memory)
                            G.nodes[node_comm]["tools_hist"]=update_tools_hist(G.nodes[node_comm]["tools_hist"],new_tool)

            #Track temporal dynamics
            if (new_tool>-1):
                if (matrix_node_tools[new_tool]<0):
                    matrix_node_tools[new_tool]=epoch+1


        prob_out_t.append(sum(prob_out_t_i)/N)
        print('time:'+str(epoch+1)+'/'+str(Epochs_max))

    ###TRACKING MEASURES
    for node in G.nodes():
        G.nodes[node]["cult_line"]=give_line(G.nodes[node]["tools"])
        G.nodes[node]["cult_level"]=give_level(G.nodes[node]["tools"])
        G.nodes[node]["max_fitness"]=max_fit_value(G.nodes[node]["tools"])

    levels=list(nx.get_node_attributes(G, 'cult_level').values())
    av_max_fitness=np.mean(list(nx.get_node_attributes(G, 'max_fitness').values()))
    av_extra_degree=np.mean(list(dict(G.degree()).values()))-n+1
    
    return (G, levels, av_max_fitness, av_extra_degree, max_fit_t, gini_fit_t, max_fit_av_t, family_fit_av_t,
            inclusive_fit_av_t, prob_out_t, max_level_t, gini_level_t, max_level_av_t, fraction_knowledge_av_t,
            gini_fraction_knowledge_t, n_tools_av_t, fitness_all_av_t, gini_specialization_t, specialized_individuals_av_t,
           pairwise_dist_t,fraction_knowledge_fit_av_t, dic_dist)


memories= [8,12,16,20]
sharing_discoveriess= [False,True]
runs=50
for counting in range(50):
    for shar in sharing_discoveriess:
        for memo in memories:
            nomeee='modelrun'+'_'+str(shar)+'_'+str(memo)+'_'+str(counting)+'.pkl'
            model_run = modello(Epochs_max=150, memory=memo, sharing_discoveries=shar)
            print(nomeee)
            with open(nomeee, 'wb') as pickle_file:
                pickle.dump(model_run, pickle_file)