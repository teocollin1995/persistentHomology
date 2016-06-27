#!/usr/bin/env python
import numpy as np
import sympy as sy
import itertools 
import networkx
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
import copy
#import networkx as nx
#plt.ion() #for interactive ploting

def unzip(z): return(zip(*z))

def concat(tlist):
    print(tlist)
    xs = []
    ys = []
    for x in range(0,len(tlist)):

        xs.extend(list(tlist[x][0]))
        ys.extend(list(tlist[x][1]))
    return([xs,ys])
        

def eucliddistance(p1,p2,dim):
    #assume len(p1) == len(p2) == dim
    return((sum(map(lambda x: (p1[x] - p2[x])**2, range(0,dim))))**(0.5))

#https://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-in-python-whilst-preserving-order
#def list_eq_if_shares_all_elements(l1,l2):
#    return(set(l1) == set(l2))
#any(map(lambda y: set(y) == set(x), seen))

def f7(seq):
    seen = list()
    seen_add = seen.append
    return [x for x in seq if not (any(map(lambda y: set(y) == set(x), seen)) or seen_add(x))]


class Data():

    def __init__(self,dim,data_list,distance, N):
        self.dim = dim
        self.points = data_list #should be list of lists where each list is of size dim and made of floats
        self.datasize = len(data_list)
        self.distance = distance #distance function for two points
        self.pdict = dict()
        self.max = None
        self.min = None
        self.distancecount = N
        self.distances = None
        

    def generate_dict(self):
        tuples = (itertools.combinations_with_replacement(range(0,self.datasize),2))
        for x in tuples:

            self.pdict[(self.points[x[0]],self.points[x[1]])] = self.distance(self.points[x[0]],self.points[x[1]])
            if x[0] != x[1]:
                self.pdict[(self.points[x[1]],self.points[x[0]])] = self.distance(self.points[x[0]],self.points[x[1]])

        l = self.pdict.values()
        self.max = max(l) * 2.0
        self.min = min(l) / 2.0

    def generate_distancelist(self):
        increment = (self.max - self.min)/ float(self.distancecount)
        self.distances = range(int(self.min),int(self.max)+1,int(increment))



class Filtration():
    def __init__(self,data, d = True):
        #data is a member of the data class,
        self.data = data
        self.data.generate_dict()
        self.data.generate_distancelist()
        self.simplicies = [] #this is a list of simplicies. A simplex is a list of points so [p1] is a 0 simple, [p1,p2] is a 2 and so one
        self.display = d
        #if d and data.dim == 2:
        #    self.fig = plt.figure()

    def build_epsilon(self,epsilon):
        avaliable_edges = (filter(lambda x: self.data.pdict[x] <= epsilon, self.data.pdict.keys()))
        #        print(avaliable_edges)
        temp_sim = []
        for x in range(0,self.data.dim+1):
            #print(x)
            if x == 0: 
                temp_sim += map(lambda x: [x], list(set(map(lambda x: x[0], avaliable_edges))))

            else:
                lastdim = filter(lambda y: len(y) == x, temp_sim)
                #print(temp_sim)
                #print(lastdim)
                newdim = []
                for y in lastdim:
                    extensions = map(lambda z: [z], list(set(map(lambda z: z[1], filter(lambda x: x[0] in y and x[1] not in y , avaliable_edges))))) #still flawed...
                    #extensions = map(lambda x: [x[1]], filter(lambda x: all(map(lambda z: (z == x[0]) and z != x, y)), avaliable_edges))
                    #print("o:", y,"ext:",extensions)
                    #for all points in the simplex, the point in the data is within epsilon of every point
                    ychildren = filter(lambda z: len(set(z)) == x+1, map(lambda q: y + q, extensions))
                    #print("child:",ychildren)
                    newdim +=  ychildren #erg
                    
                    
                temp_sim += f7(newdim)
        #horribly inefficient 
        
        for x in self.simplicies:
            try:
                temp_sim.remove(x)
            except:
                pass 
 #        lines = []
 #        points = [[],[]]
 #        for x in f7(avaliable_edges):
 #            if x[0] != x[1]:
 # #               lines.extend([[x[0][0],x[1][0]],[x[0][1],x[1][1]], 'b'])
 #                plt.plot([x[0][0],x[1][0]],[x[0][1],x[1][1]], 'b')
 #                plt.draw()
 #                plt.pause(0.1)
 #                print("erg")

#        plt.plot(*lines)
#        plt.draw()
        # plotting region 
        # for z in temp_sim:
        #     line = [[],[]]
        #     for y in z:
        #         line[0].append(y[0])
        #         line[1].append(y[1])
        #     plt.plot(line)
        #     plt.draw()
        #     plt.pause(0.05)
        # showing = concat(map(unzip,temp_sim))
        # t = len(showing[0])
        # if t > 0:        
        #     plt.plot(showing,'ro')
        #     plt.draw()
        #     plt.pause(0.05)
        self.simplicies += temp_sim
            
    def build_simplex_complex(self):
        for x in self.data.distances:
            # if self.display:
            #     for p in self.data.points:
            #         c = plt.Circle(p,x,color='b')
            #         self.fig.gca().add_artist(c)
            #                plt.draw()
            #            plt.pause(5)
            self.build_epsilon(x)

    def build_boundary_matrix(self): #from here on out is an implementation of algorthims described in p2.pdf
        self.simplicies_count = len(self.simplicies)
        self.boundary_matrix = sy.Matrix([[0.0 for x in range(0,self.simplicies_count)] for x in range(0,self.simplicies_count)])
        for i in range(0,self.simplicies_count):
            s = self.simplicies[i]
            current_dim = len(s)
            codiemsnion_one_simplicies = filter(lambda x: len(x) == current_dim+1 and all(map(lambda q: q in x, s)) , self.simplicies)
            if codiemsnion_one_simplicies == []:
                pass
            else:
                for j in codiemsnion_one_simplicies:
                    self.boundary_matrix[i,self.simplicies.index(j)] = 1

    def low(self,j,m):
        col = m.col(j)
        onesandindexes = [x for x in range(0,self.simplicies_count) if col[x] == 1]
        if onesandindexes == []:
           return(None) #undefined
        else:
            return(max(onesandindexes))
        

    def boundary_to_barcodes(self):
        for i in range(0,self.simplicies_count):
            ilow = self.low(i,self.boundary_matrix)
            if ilow != None:
                for j in range(i+1,self.simplicies_count):
                    jlow = self.low(j,self.boundary_matrix)
                    if ilow == jlow:
                        for x in range(0,self.simplicies_count):
                            self.boundary_matrix[x,j] = ((self.boundary_matrix[x,j] + self.boundary_matrix[x,i]) % 2) #over the field F_2
                            
            # while any([self.low(i,self.bardcode_matrix) == self.low(j,self.bardcode_matrix) and self.low(j,self.bardcode_matrix) != None for j in range(i,self.simplicies_count)]):
            #     for x in range(0,self.simplicies_count):
            #         print(x,j)
                    
            #         self.bardcode_matrix[x,j] = self.bardcode_matrix[x,j] + self.bardcode_matrix[x,i]
                
            
                
        
        
        




            
        

    


        


