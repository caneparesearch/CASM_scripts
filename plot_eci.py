import os,json
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.patches as patches
from matplotlib.ticker import MultipleLocator

class ClusterFunctionGroup: # this is a collection of ClusterFunction object
    def __init__(self,cluster_function_list):
        self.cluster_function_list = cluster_function_list
        self.num_cluster_function = len(self.cluster_function_list)
        self.eff_cluster_list = self.get_eff_cluster_list(tol=0.1)
        self.num_eci = len(self.eff_cluster_list)

    def construct_from_file(fname):
        try:
            with open(fname,'r') as f:
                data = json.load(f)
        except IOError:
            raise Exception('Cannot find',fname,' in this directory: ',os.getcwd())
        clex_funcs_json = data['orbits']
        cluster_function_list = []
        for i in clex_funcs_json:
            cluster_function_list.append(ClusterFunction(i))
        return ClusterFunctionGroup(cluster_function_list)

    def append(self,new_clex_function): # append a cluster fuction to this group
        self.cluster_function_list.append(new_clex_function)
        self.num_cluster_function+=1
        
    def get_property_list(self,prop):
        return [i.get_property(prop) for i in self.cluster_function_list]

    def get_eff_property_list(self,prop):
        return [i.get_property(prop) for i in self.eff_cluster_list]
    
    def select_cluster_functions(self,prop,cond):# return another ClusterFunctionGroup with prop==cond
        return ClusterFunctionGroup([i for i in self.cluster_function_list if i.get_property(prop)==cond])
    
    def get_eff_cluster_list(self,tol=0.1): # calculate the number of ECI (exclude them with abs(normalized_eci) smaller than tol (default: 0.1 meV/multiplicity))
        return [i for i in self.cluster_function_list if abs(i.get_property('normalized_eci'))>tol]

    def print_info(self):
        print('There are',self.num_cluster_function,'cluster functions in this group. N_ECI = ',self.num_eci)

class ClusterFunction: # this is used to deal with cluster_function section in eci.json
    def __init__(self,clex_func_json):
        self.linear_function_index = clex_func_json["cluster_functions"][0]['linear_function_index']
        self.index = self.linear_function_index+1 # index start from 1 instead of 0
        self.mult = clex_func_json['mult']
        #self.orbit = clex_func_json['orbit']
        self.prototype = clex_func_json['prototype']
        #self.prototype_function = clex_func_json['prototype_function']
        if len(self.prototype['sites']) == 0:
            self.category = 'Empty'
        elif len(self.prototype['sites']) == 1:
            self.category = 'Point'
        elif len(self.prototype['sites']) == 2:
            self.category = 'Pair'
        elif len(self.prototype['sites']) == 3:
            self.category = 'Triplet'
        elif len(self.prototype['sites']) == 4:
            self.category = 'Quadruplet'
        else: 
            raise Exception('Error occurred in number of sites!')
        if 'eci' in clex_func_json["cluster_functions"][0]:
            self.eci = clex_func_json["cluster_functions"][0]['eci']
            self.normalized_eci = self.eci/self.mult*1000 # convert into meV normalized by multiplicity
        else:
            self.eci = 0.0
            self.normalized_eci = 0.0
    def get_property(self,prop):
        return self.__dict__[prop]
    def print_info(self):
        print(self.linear_function_index,self.mult,self.prototype_function,self.eci,self.category)
        
        
def draw_stem(ax,group,x,y,color,marker,label):# do a stem plot for "group" with "x" vs "y" using "color" and 'marker' with "label" on "ax"
    markerline, stemlines, baseline=ax.stem(group.get_eff_property_list(x),group.get_eff_property_list(y),markerfmt='.',use_line_collection=True,label=label)
    plt.setp(stemlines,'color',color,'lw',1)
    plt.setp(markerline,'color',color,'lw',1,'marker',marker,'ms',3)
    plt.setp(baseline,'visible',False)   


# initialize ClusterFunctionGroup from eci.json file
clust_funcs = ClusterFunctionGroup.construct_from_file('eci.json')
clust_funcs.print_info()
# slice cluster_funcs based on their category
# point=clust_funcs.select_cluster_functions('category','Point')
pair=clust_funcs.select_cluster_functions('category','Pair')
triplet=clust_funcs.select_cluster_functions('category','Triplet')
quadruplet=clust_funcs.select_cluster_functions('category','Quadruplet')
print(pair.get_eff_property_list('index'))
print(pair.get_eff_property_list('eci'))
print(pair.get_eff_property_list('normalized_eci'))
print(quadruplet.get_eff_property_list('index'))
print(quadruplet.get_eff_property_list('eci'))
print(quadruplet.get_eff_property_list('normalized_eci'))

fig, ax = plt.subplots(1,figsize=(4,4))
# draw_stem(ax,point,'index','normalized_eci','#384259','Point')
draw_stem(ax,pair,'index','normalized_eci','#e74c3c','o','Pair')
draw_stem(ax,triplet,'index','normalized_eci','#3498db','^','Triplet')
draw_stem(ax,quadruplet,'index','normalized_eci','#95e1d3','s','Quadruplet')

ax.legend() 
ax.axhline(y=0, xmin=0.0, xmax=300, marker='', linestyle='--', linewidth=0.5, color="black", antialiased=True, label="")
ax.set_ylabel('ECI/multiplicity (meV)')
ax.set_xlabel('Cluster Function Index')
ax.set_xlim(0,25)
#ax.set_ylim(-100,200)
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_major_locator(MultipleLocator(5))
fig.tight_layout()
fig.savefig('eci.pdf')