#!/usr/bin/env python
import numpy as np
import sympy as sy
import itertools 
import networkx
import turtle as t

def draw_circle(x,y,r,t):
    t.pu()
    t.goto(x+r,y)
    t.setheading(90)
    t.pd()
    t.circle(r)
    t.pu()
    

#helper functions:
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



def f7(seq):
    seen = list()
    seen_add = seen.append
    return [x for x in seq if not (any(map(lambda y: set(y) == set(x), seen)) or seen_add(x))]


#data storage:

class Data():

    def __init__(self,dim,data_list,distance, N):
        self.dim = dim
        self.points = data_list #should be list of lists where each list is of size dim and made of floats
        self.datasize = len(data_list)
        self.distance = distance #distance function for two points
        self.pdict = dict() #dictionary of (point1,point2) -> distance
        self.max = None
        self.min = None
        self.distancecount = N
        self.distances = None
        

    def generate_dict(self): #fills pdict
        tuples = (itertools.combinations_with_replacement(range(0,self.datasize),2))
        for x in tuples:

            self.pdict[(self.points[x[0]],self.points[x[1]])] = self.distance(self.points[x[0]],self.points[x[1]])
            if x[0] != x[1]:
                self.pdict[(self.points[x[1]],self.points[x[0]])] = self.distance(self.points[x[0]],self.points[x[1]])

        l = self.pdict.values()
        self.max = max(l) * 2.0
        self.min = min(l) / 2.0


    def generate_distancelist(self): #prepares the distances we will go through to build our complex
        increment = (self.max - self.min)/ float(self.distancecount)
        self.distances = range(int(self.min),int(self.max)+1,int(increment))


#d = [(0,0),(3,3),(6,0),(3,8),(7,9),(10,15),(13,9),(15,7),(19,2),(17,0),(26,0),(18,-5),(10,-7),(8,-11),(12,-10),(4,-3)]
class Filtration():
    def __init__(self,data, d = False,p =False, computed_matrix = None):
        self.data = data        #data is a member of the data class,
        self.data.generate_dict()
        self.data.generate_distancelist()
        self.simplicies = [] 
        #this is a list of simplicies. A simplex is a list of points so [p1] is a 0 simple, [p1,p2] is a 2 and so one
        self.simplicies_index = dict() 
        #dict where the index of a simplix in self.simplices corresponds to the epsilon at which it was created
        self.d = d #display
        self.p = p #present
        if p:
            self.precomputed_boundary_matrix = computed_matrix
        if d and data.dim == 2:
            self.r = None
            self.expand_factor = 10 #figure out appropraite way to calculate this.
            self.turtle = t
            self.turtle.hideturtle()
            self.turtle.speed(0)#fastest possible speed
            for x in data.points:
                self.turtle.pu()
                self.turtle.goto((x[0]*self.expand_factor,x[1]*self.expand_factor))
                self.turtle.pd()
                self.turtle.dot()
                

    def build_epsilon(self,epsilon):
        if self.d == True:
            r = epsilon * self.expand_factor
            for x in self.data.points:
                self.turtle.pen(pencolor="red")
                draw_circle(self.expand_factor*x[0],self.expand_factor*x[1],r,self.turtle)
                if self.r != None:
                     self.turtle.pen(pencolor="white")
                     draw_circle(self.expand_factor*x[0],self.expand_factor*x[1],self.r,self.turtle)
                     self.turtle.pen(pencolor="red")
            
        self.r = r
        avaliable_edges = (filter(lambda x: self.data.pdict[x] <= epsilon, self.data.pdict.keys()))
        temp_sim = []
        for x in range(0,self.data.dim+1):
            if x == 0: 
                temp_sim += map(lambda x: [x], list(set(map(lambda x: x[0], avaliable_edges))))

            else:
                lastdim = filter(lambda y: len(y) == x, temp_sim)
                newdim = []
                for y in lastdim:
                    extensions = [[z[1]] for z in avaliable_edges if all(map(lambda x: (x,z[1]) in avaliable_edges,y))]
                    ychildren = filter(lambda z: len(set(z)) == x+1, map(lambda q: y + q, extensions))
                    newdim +=  ychildren #erg
                new = f7(newdim)


                        
                temp_sim += new


        #horribly inefficient 
        for x in self.simplicies:
            try:
                temp_sim.remove(x)
            except:
                pass 

        if self.simplicies == []:
            for i in range(0,len(temp_sim)):

                self.simplicies_index[i] = (epsilon)
        else:
            s = len(self.simplicies)

            e = s + len(temp_sim)

            last_epsilon = self.simplicies_index[s-1]
            for i in range(s,e):

                self.simplicies_index[i] = (epsilon+last_epsilon)

            if d and temp_sim != []:
                for sim in temp_sim:

                    dim = len(sim)
                    if dim == 1:
                        pass
                    else:
                        for i in range(0,dim-1):
                            self.turtle.pu()
                            self.turtle.goto(sim[i][0]*self.expand_factor,sim[i][1]*self.expand_factor)
                            if i == 0:
                                self.turtle.begin_fill()
                            self.turtle.pd()
                            self.turtle.goto(sim[i+1][0]*self.expand_factor,sim[i+1][1]*self.expand_factor)
                            if i+1 == dim-1:
                                self.turtle.end_fill()

        self.simplicies += temp_sim
        self.simplicies = f7(self.simplicies)
            
    def build_simplex_complex(self): #O(2^|S|-1)
        try:
            for x in self.data.distances:
                self.build_epsilon(x)
        except KeyboardInterrupt:
            return
    def build_boundary_matrix(self): 
        #from here on out is an implementation of algorthims described in p2.pdf
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

    def low(self,j,m): #defined in p2.pdf
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
                            self.boundary_matrix[x,j] = ((self.boundary_matrix[x,j] + self.boundary_matrix[x,i]) % 2) #adding cols over the field F_2

    def read_barcodes(self):
        barcode = []
        lows = map(lambda j: self.low(j,self.boundary_matrix),range(0,self.simplicies_count))
        for x in range(0,self.simplicies_count):
            l = lows[x]
            if l != None:
                barcode.append(((l,self.simplicies_index[l]),(x,self.simplicies_index[x]))) #something is born with l and lies with x

            else:
                if x in lows:
                    inverse_x = lows.index(x)
                    barcode.append(((x,self.simplicies_index[x]),(inverse_x,self.simplicies_index[inverse_x])))
                else:
                    barcode.append(((x,self.simplicies_index[x]),(None,None))) #no death of the feature
        return(list(set(barcode)))

    def run(self):
        print("Building Complex")
        self.build_simplex_complex()
        print("Complex built")
        print("Generating Boundary Matrix")
        self.build_boundary_matrix()
        print("Generated Boundary Matrix")
        print("Generating barcodes")
        self.boundary_to_barcodes()
        print("Barcodes ready:")
        return(self.read_barcodes())
        
                
  
        
        
        




            
        

    


        


