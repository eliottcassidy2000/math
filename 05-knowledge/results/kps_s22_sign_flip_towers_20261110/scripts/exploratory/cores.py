"""cores.py - strong cores of unit sets. Input us_k.txt lines: "order_of_L codeB Omask H(B) M" (B of order k-1).
For each unit set find the strong component Q containing O (lemma: all other components are single vertices),
and print statistics: core order, feedback vertex number of Q, |O|, H(Q), odd part c of M."""
import sys, itertools, collections
def decode(code,m):
    A=[[0]*m for _ in range(m)]; b=0
    for i in range(m):
        for j in range(i+1,m):
            if code>>b&1: A[i][j]=1
            else: A[j][i]=1
            b+=1
    return A
def comps(A):
    n=len(A); reach=[[A[i][j] or i==j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    if reach[k][j]: reach[i][j]=1
    comp=[min(j for j in range(n) if reach[i][j] and reach[j][i]) for i in range(n)]
    return comp
def transitive(A,keep):
    sc=sorted(sum(A[i][j] for j in keep) for i in keep)
    return sc==list(range(len(keep)))
def fvs(A):
    n=len(A)
    for k in range(0,n+1):
        for rem in itertools.combinations(range(n),k):
            keep=[i for i in range(n) if i not in rem]
            if transitive(A,keep): return k
stats=collections.Counter(); viol=0; byc=collections.defaultdict(set)
for fn in sys.argv[1:]:
    for line in open(fn):
        k,code,O,H,M=map(int,line.split()); m=k-1
        A=decode(code,m); comp=comps(A); Os=[i for i in range(m) if O>>i&1]
        cs=set(comp[i] for i in Os)
        sizes=collections.Counter(comp)
        if len(cs)!=1 or any(sizes[c]>1 for c in sizes if c not in cs): viol+=1; continue
        core=[i for i in range(m) if comp[i]==comp[Os[0]]]
        Q=[[A[i][j] for j in core] for i in core]
        c=M
        while c%2==0: c//=2
        f=fvs(Q)
        stats[(len(core),f,len(Os))]+=1; byc[(len(core),f)].add(c)
print("violations of the one-component lemma:",viol)
print("(core order, fvs, |O|): count")
for key in sorted(stats): print(key,stats[key])
print("distinct odd c by (core order, fvs):",{k:len(v) for k,v in sorted(byc.items())})
print("c from one-apex cores (fvs=1):",sorted(set().union(*[v for k,v in byc.items() if k[1]==1]))[:60])
