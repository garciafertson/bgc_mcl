import numpy as np
from collections import defaultdict
class sim_index:
    def __init__(self,index_type,a,b):
        self.a=a
        self.b=b
        self.index_type=index_type
    def jaccard_index(self):
        #calculate jaccard index for pairs of bgc
        a=set(self.a)
        b=set(self.b)
        ji=len(a.intersection(b))/len(a.union(b))
        return(ji)
    def overlap_coef(self):
        a=set(self.a)
        b=set(self.b)
        oc=len(a.intersection(b))/min(len(a),len(b))
        return(oc)
    def domdupsim(self): #Measures the amount of duplication of a shared domain
        #get elements in  union
        union=set(self.a).union(set(self.b))
        #count the number of times each domain exists in each set
        count_a=defaultdict(int)
        for element in self.a: #count number of times each element appears in a
            count_a[element]+=1
        count_b=defaultdict(int)
        for element in self.b:
            count_b[element]+=1
        #calculate S and ref Kui Lin 2006
        S=0
        for dom in union:
            S+= max(count_a[dom], count_b[dom])
        #calculate D index [0.36 ,1]
        D=0
        x=0
        if S>0:
            for dom in union:
                x+= abs(count_a[dom]-count_b[dom])/S
            D=np.e**-x
            #D= (np.e**-x) - (x*(np.e**-x))
            #D=1-sm
        return(D)
    def ji_ds(self):
        coef_ji=0.365
        coef_ds=0.635
        j_d= coef_ji*self.jaccard_index()+ coef_ds*self.domdupsim()
        return(j_d)

    def main(self):
        if self.index_type==1:
            index=self.overlap_coef()
        elif self.index_type==2:
            index=self.jaccard_index()
        elif self.index_type==3:
            index=self.domdupsim()
        elif self.index_type==4:
            index=self.ji_ds()
        else:
            sys.exit("Error: --index, Index type not definded")
        return(index)

