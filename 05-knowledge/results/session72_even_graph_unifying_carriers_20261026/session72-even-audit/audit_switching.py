#!/usr/bin/env python3
"""Independent root-gauge audit of switching quotient Fourier/gap statements.
Every source switching class is gauge-fixed with vertex zero a source.
Every dual Eulerian graph is reconstructed from non-root chord coordinates.
No repository implementation is imported. All gates survive python -O.
"""
import ast
import itertools
import json
from pathlib import Path

CHECKS=0

def check(ok,message):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(message)

def parity(x): return x.bit_count()&1

def perm_mask(x,images):
    out=0
    while x:
        b=x&-x
        out^=images[b.bit_length()-1]
        x^=b
    return out

def det_exact(matrix):
    a=[r[:] for r in matrix]
    if not a: return 1
    prev=1
    sign=1
    for k in range(len(a)-1):
        j=next((j for j in range(k,len(a)) if a[j][k]),None)
        if j is None: return 0
        if j!=k:
            a[j],a[k]=a[k],a[j]
            sign=-sign
        pivot=a[k][k]
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):
                numerator=pivot*a[i][j]-a[i][k]*a[k][j]
                check(numerator%prev==0,'integer determinant division')
                a[i][j]=numerator//prev
            a[i][k]=0
        prev=pivot
    return sign*a[-1][-1]

class Switching:
    def __init__(self,n):
        self.n=n
        self.edges=list(itertools.combinations(range(n),2))
        self.m=len(self.edges)
        idx={e:i for i,e in enumerate(self.edges)}
        self.chords=[i for i,(u,v) in enumerate(self.edges) if u>0]
        self.size=1<<len(self.chords)
        self.expand=[sum(1<<self.chords[j] for j in range(len(self.chords)) if q>>j&1) for q in range(self.size)]
        self.cuts=[]
        for s in range(1<<max(0,n-1)):
            self.cuts.append(sum(1<<i for i,(u,v) in enumerate(self.edges) if ((s>>(u-1)&1) if u else 0)^((s>>(v-1)&1) if v else 0)))
        self.euler=[]
        for x in self.expand:
            y=x
            for v in range(1,n):
                if parity(x&sum(1<<i for i,e in enumerate(self.edges) if v in e)):
                    y^=1<<idx[(0,v)]
            self.euler.append(y)
        self.perms=[]
        for p in itertools.permutations(range(n)):
            images=[1<<idx[tuple(sorted((p[u],p[v])))] for u,v in self.edges]
            invsrc=sum(1<<i for i,(u,v) in enumerate(self.edges) if p[u]>p[v])
            flip=perm_mask(invsrc,images)
            self.perms.append((images,invsrc,flip))
        self.sorbits,self.scanon=self.orbits(True)
        self.eorbits,self.ecanon=self.orbits(False)
    def compress(self,x):
        return sum(1<<j for j,i in enumerate(self.chords) if x>>i&1)
    def normalize(self,x):
        s=sum(1<<(v-1) for v in range(1,self.n) if x>>(v-1)&1)
        return self.compress(x^self.cuts[s])
    def act(self,q,p,affine):
        images,invsrc,flip=self.perms[p]
        if affine: return self.normalize(perm_mask(self.expand[q],images)^flip)
        return self.compress(perm_mask(self.euler[q],images))
    def orbits(self,affine):
        orbits={}
        canon=[None]*self.size
        for q in range(self.size):
            if canon[q] is not None: continue
            orb={self.act(q,p,affine) for p in range(len(self.perms))}
            check(min(orb)==q,'orbit minimality')
            for r in orb: canon[r]=q
            orbits[q]=sorted(orb)
        return orbits,canon
    def quartic(self,q):
        x=self.expand[q]
        K=[[0]*self.n for _ in range(self.n)]
        for i,(u,v) in enumerate(self.edges):
            K[u][v]=-1 if x>>i&1 else 1
            K[v][u]=-K[u][v]
        K2=[[sum(K[i][k]*K[k][j] for k in range(self.n)) for j in range(self.n)] for i in range(self.n)]
        return sum(a*a for row in K2 for a in row)-self.n*(self.n-1)*(2*self.n-3)
    def run(self):
        even=[]
        for q in self.eorbits:
            x=self.euler[q]
            if all(not parity(x&self.perms[p][1]) for p in range(len(self.perms)) if self.act(q,p,False)==q): even.append(q)
        source=list(self.sorbits)
        check(len(even)==len(source),'switching/even-Euler equinumerosity')
        cols=[]
        for q in even:
            x=self.euler[q]
            signed={}
            for p,(images,invsrc,flip) in enumerate(self.perms):
                y=perm_mask(x,images)
                sign=-1 if parity(x&invsrc) else 1
                check(y not in signed or signed[y]==sign,'sign transport well-defined')
                signed[y]=sign
            vals=[sum(s*(-1 if parity(a&y) else 1) for y,s in signed.items()) for a in self.expand]
            lam=self.m-2*x.bit_count()
            for a in range(self.size):
                check(vals[a]==vals[self.scanon[a]],'switching Fourier invariant')
                check(sum(vals[self.normalize(self.expand[a]^(1<<e))] for e in range(self.m))==lam*vals[a],'native quotient eigenfunction')
            cols.append([vals[a] for a in source])
        matrix=[list(row) for row in zip(*cols)]
        check(det_exact(matrix)!=0,'full eigenbasis')
        quartic=[self.quartic(q) for q in range(self.size)]
        for q in range(self.size):
            check(quartic[q]==quartic[self.scanon[q]],'quartic orbit invariant')
            check(sum(quartic[self.normalize(self.expand[q]^(1<<e))] for e in range(self.m))==(self.m-8)*quartic[q],'quartic eigen8')
        positive=[2*self.euler[q].bit_count() for q in even if self.euler[q]]
        if self.n>=4:
            check(min(positive)==8 and positive.count(8)==1,'gap8 simple')
            check(any(quartic),'quartic nonzero')
            four=next(i for i,q in enumerate(even) if self.euler[q].bit_count()==4)
            check(all(quartic[a] in (8*cols[four][i],-8*cols[four][i]) for i,a in enumerate(source)),'quartic Fourier scaling')
            nonzero=next(i for i,a in enumerate(source) if cols[four][i])
            sign=quartic[source[nonzero]]//(8*cols[four][nonzero])
            check(sign in (-1,1) and all(quartic[a]==sign*8*cols[four][i] for i,a in enumerate(source)),'single global quartic phase')
        else:
            check(not positive and len(source)==1 and not any(quartic),'small n one-state boundary')
        idx={q:i for i,q in enumerate(source)}
        adjacency=[[0]*len(source) for _ in source]
        for i,q in enumerate(source):
            for e in range(self.m): adjacency[i][idx[self.scanon[self.normalize(self.expand[q]^(1<<e))]]]+=1
        return {'n':self.n,'labelled_switching_classes':self.size,'source_type_count':len(source),'ordinary_eulerian_type_count':len(self.eorbits),'sign_even_eulerian_type_count':len(even),'laplacian_eigenvalues':[2*self.euler[q].bit_count() for q in even],'quartic_by_type':[quartic[q] for q in source],'adjacency':adjacency,'weighted_rowsums':self.m}

def main():
    rows=[Switching(n).run() for n in range(1,7)]
    tree=ast.parse(Path(__file__).read_text())
    check(not any(isinstance(x,ast.Assert) for x in ast.walk(tree)),'no assert nodes')
    check(not any(isinstance(x,ast.Constant) and isinstance(x.value,float) for x in ast.walk(tree)),'no floating constants')
    print(json.dumps({'status':'PASS','checks':CHECKS,'universe':'Every labelled switching class and Eulerian graph n=1..6; all vertex permutations; root gauge exact','rows':rows},indent=2,sort_keys=True))

if __name__=='__main__': main()
