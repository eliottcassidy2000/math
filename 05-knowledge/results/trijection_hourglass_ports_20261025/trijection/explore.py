from itertools import combinations

def edge(u,v): return min(u,v),max(u,v)
def edges(q): return {edge(u,v) for u,v in zip(q,q[1:])}

class Parity:
    def __init__(self,n):
        self.pairs=list(combinations(range(n),2));self.ids={p:i for i,p in enumerate(self.pairs)}
        self.c=len(self.pairs);self.p=list(range(self.c+1));self.x=[0]*(self.c+1);self.ok=True
    def find(self,a):
        if self.p[a]!=a:
            r,t=self.find(self.p[a]);self.x[a]^=t;self.p[a]=r
        return self.p[a],self.x[a]
    def join(self,a,b,v):
        a,x=self.find(a);b,y=self.find(b)
        if a==b:self.ok&=(x^y)==v
        else:self.p[a]=b;self.x[a]=x^y^v
    def eq(self,u,v,a,b,value=0):
        self.join(self.ids[edge(u,v)],self.ids[edge(a,b)],value^(u>v)^(a>b))
    def force(self,u,v,value=1):self.join(self.ids[edge(u,v)],self.c,value^(u>v))
    def anti(self,q,flips=set()):
        r=dict(zip(q,q[::-1]))
        for u,v in self.pairs:
            a,b=r[v],r[u]
            self.eq(u,v,a,b,(edge(u,v) in flips)^(edge(a,b) in flips))
    def bits(self):
        rc,xc=self.find(self.c)
        return [self.find(i)[1]^(xc if self.find(i)[0]==rc else 0) for i in range(self.c)]

def source(n,h,ell):
    d=Parity(n);P=tuple(range(n));d.anti(P)
    for u,v in zip(P,P[1:]):d.force(u,v)
    I={h,n-3,n-2}
    d.force(n-1,0)
    for v in range(1,n-1):d.force(v,0,int(v in I))
    d.force(ell+3,2)
    return d

def words(n,ell):
    D=tuple(range(3,ell+3));U=tuple(range(ell+3,n-1))
    return ((0,2)+D+U+(n-1,1),D+(0,2)+U+(n-1,1),D+(0,1,2)+U+(n-1,))

def feasible(n,h,ell,j):
    d=source(n,h,ell)
    if not d.ok:return None
    qs=words(n,ell)
    flips=set()
    if j==1:flips=set(d.pairs)^edges(qs[0])
    if j==2:flips=set(d.pairs)^edges(tuple(range(n)))
    d.anti(qs[j],flips)
    return d.ok

if __name__=='__main__':
    for n in range(15,35):
        out=[]
        for ell in range(1,n-3,2):
            bad=[]
            for h in range(2,n-3):
                if h in {2,ell+2}:continue
                f=[feasible(n,h,ell,j) for j in (1,2)]
                if any(f):bad.append((h,[j+2 for j,x in enumerate(f) if x]))
            if bad:out.append((ell,bad))
        print(n,out)
