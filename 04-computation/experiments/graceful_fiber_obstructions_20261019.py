"""Exact small controls for fixed-upper-tail graceful tree exchange fibers."""
from itertools import product, combinations
from collections import Counter

def check(x,msg):
    if not x: raise RuntimeError(msg)

def components(n,edges):
    p=list(range(n))
    def root(x):
        while p[x]!=x:
            p[x]=p[p[x]]; x=p[x]
        return x
    for a,b in edges:
        a,b=root(a),root(b)
        if a==b:return None
        p[a]=b
    return tuple(root(i) for i in range(n))

def classes(states):
    states=set(states); out=[]
    while states:
        seed=min(states);states.remove(seed);todo=[seed]
        for a in todo:
            near=[b for b in states if sum(x!=y for x,y in zip(a,b))==1]
            for b in near:states.remove(b);todo.append(b)
        out.append(tuple(sorted(todo)))
    return tuple(out)

def support(n,c,d):
    return {tuple(sorted((c[a],c[a+d]))) for a in range(n-d) if c[a]!=c[a+d]}

def two_fiber(n,starts):
    fixed=[(a,a+d) for d,a in zip(range(3,n),starts)]
    c=components(n,fixed)
    if c is None:return None
    check(len(set(c))==3,'three components')
    q1,q2=support(n,c,1),support(n,c,2)
    check(len(q1)>=2 and len(q2)>=2,'connected three-vertex supports')
    states=[]
    for a,b in product(range(n-1),range(n-2)):
        if components(n,fixed+[(a,a+1),(b,b+2)]) is not None:states.append((a,b))
    cc=classes(states)
    blocked=q1==q2 and len(q1)==2
    check(len(cc)==(2 if blocked else 1),'two-color classification')
    return blocked,cc

def certificate(k,m,d):
    fixed=[(1,m) if j==m-1 else (0,j) for j in range(k+1,m+1)]
    b=m-1
    free=[(b-j,b) if j==d else (1,j+1) if j<d else (0,j) for j in range(1,k+1)]
    tree=fixed+free
    check(len(tree)==m and components(m+1,tree) is not None,'certificate tree')
    check(sorted(abs(a-b) for a,b in tree)==list(range(1,m+1)),'graceful spectrum')
    comp=components(m+1,fixed); h=comp[0]
    check(comp[b]!=h,'isolated distinguished vertex')
    for j in range(1,k+1):
        for a in (b-j,b+j):
            if 0<=a<=m:check(comp[a]==h,'all short neighbors in hub')
    return tuple(a for a,b in free)

def main():
    check(len(classes([(0,0),(1,0)]))==1, "n3 two-color boundary")
    total=0;bad=0
    for n in range(4,11):
        counts=Counter()
        for starts in product(*(range(n-d) for d in range(3,n))):
            r=two_fiber(n,starts)
            if r is not None:counts[r[0]]+=1
        total+=sum(counts.values());bad+=counts[True]
        print('two-free-color fibers n=%d connected=%d disconnected=%d'%(n,counts[False],counts[True]))
    r=two_fiber(7,(0,0,1,0))
    check(r is not None and r[0],'order-seven hostile')
    expected=(((1,3),(2,3)),((4,0),(4,2),(5,0),(5,2)))
    check(r[1]==expected,'exact six-state fiber')
    print('n=7 exact components:',r[1])
    print('classification checked on',total,'acyclic tails;',bad,'disconnected fibers')
    count=0
    for k in range(2,21):
        for m in (2*k+2,2*k+3,3*k+5):
            for d in range(1,k+1):certificate(k,m,d);count+=1
    print('infinite-family positive certificates checked:',count,'(the all-k proof is in the note)')
    # Independent enumeration of all completions for three family examples.
    for k,m in ((2,6),(3,8),(4,10)):
        fixed=[(1,m) if j==m-1 else (0,j) for j in range(k+1,m+1)]
        states=[]; colors={}
        for a in product(*(range(m+1-j) for j in range(1,k+1))):
            free=[(x,x+j) for j,x in enumerate(a,1)]
            if components(m+1,fixed+free) is None:continue
            incident=[j for j,e in enumerate(free,1) if m-1 in e]
            check(len(incident)==1,'single hub attachment')
            states.append(a);colors[a]=incident[0]
        cc=classes(states)
        check(len(cc)>=k,'at least k components')
        for group in cc:check(len({colors[a] for a in group})==1,'attachment color invariant')
        print('fixed-band family k=%d n=%d states=%d components=%d'%(k,m+1,len(states),len(cc)))
    print('ALL CHECKS PASSED')
if __name__=='__main__':main()
