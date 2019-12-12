class sim_index:
    def _init_(self,index_type,a,b):
        self.a=a
        self.b=b
        self.index_type=index_type
    def jaccard_index(a,b):
        #calculate jaccard index for pairs of bgc
        a=set(a)
        b=set(b)
        ji=len(a.intersection(b))/len(a.union(b))
        return(ji)
    def overlap_coef(a,b):
        a=set(a)
        b=set(b)
        oc=len(a.intersection(b))/min(len(a),len(b))
        return(oc)
    def domdupsim(a,b): #Measures the amount of duplication of a shared domain
        #get elements in  union
        union=set(a).union(set(b))
        #count the number of times each domain exists in each set
        count_a=defaultdict(int)
        for element in a: #count number of times each element appears in a
            count_a[element]+=1
        count_b=defaultdict(int)
        for element in b:
            count_b[element]+=1
        #calculate S and ref Kui Lin 2006
        S=0
        for dom in union:
            S+= max(count_a[dom], count_b[dom])
        #calculate D index [0.36 ,1]
        exp=0
        D=0
        if S>0:
            for dom in union:
                x+= abs(count_a[dom]-count_b[dom])/S
            D=np.e**exp
            #D= (np.e**-x) - (x*(np.e**-x))
            #D=1-sm
        return(D)
    def ji_ds(a,b):
    def main(self):
        if self.index_type==1:
            index=overlap_coef(self.a,self.b)
        elif self.index_type==2:
            index=jaccard_index(self.a,self.b)
        elif sefl.index==3:
            index=domdupsim(self.a,self.b)
        elif index==4:
            index=domdupsim(self.a,self.b)
        else:
            sys.exit("Error: --index, Index type not definded")
        return(index)

