#!/usr/bin/env python3
"""Finite exact checks; no third-party dependencies, no assert statements."""
import itertools as it, json, math, sys, random

def require(p,msg):
    if not p: raise RuntimeError(msg)
def pairs(n): return list(it.combinations(range(n),2))
def adj(n,mask):
    a=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(pairs(n)): a[i][j]=(mask>>k)&1;a[j][i]=1-a[i][j]
    return a

def encode(a,p): return sum(a[p[i]][p[j]]<<k for k,(i,j) in enumerate(pairs(len(a))))
def canonical(a): return min(encode(a,p) for p in it.permutations(range(len(a))))
def hp(a): return [p for p in it.permutations(range(len(a))) if all(a[p[i]][p[i+1]] for i in range(len(a)-1))]
def normalize(a,p):
    free=[(i,j) for i,j in pairs(len(a)) if j>i+1]
    return sum((1-a[p[i]][p[j]])<<k for k,(i,j) in enumerate(free))
def observer(n,x):
    a=[[int(i<j) for j in range(n)] for i in range(n)]
    free=[(i,j) for i,j in pairs(n) if j>i+1]
    for k,(i,j) in enumerate(free): a[i][j]=1-((x>>k)&1);a[j][i]=1-a[i][j]
    return a

def census(n):
    free=[(i,j) for i,j in pairs(n) if j>i+1];m=len(free)
    def rev(x): return sum(((x>>free.index((n-1-j,n-1-i)))&1)<<k for k,(i,j) in enumerate(free))
    fibres={}; reps={}; symmetric={}
    for x in range(1<<m):
        a=observer(n,x);c=canonical(a);fibres[c]=fibres.get(c,0)+1;reps[c]=a
        symmetric[c]=symmetric.get(c,0)+(rev(x)==x)
        converse=[[a[j][i] for j in range(n)] for i in range(n)]
        require(normalize(converse,tuple(reversed(range(n))))==rev(x),'converse reflection')
        require(rev(x^((1<<m)-1))==rev(x)^((1<<m)-1),'commuting involutions')
    records=[]
    for c in sorted(fibres):
        a=reps[c];paths=hp(a);auts=[p for p in it.permutations(range(n)) if encode(a,p)==encode(a,range(n))]
        seen={normalize(a,p) for p in paths}
        require(len(seen)==fibres[c],'normalized fibre coverage')
        require(len(paths)==len(auts)*fibres[c],'free action fibre')
        records.append({'canonical_mask':c,'H':len(paths),'Aut':len(auts),'observers':fibres[c],'reflection_fixed':symmetric[c]})
    fixed=sum(symmetric.values());f=(n-1)//2
    require(fixed==2**((m+f)//2),'reflection count')
    return {'n':n,'free_bits':m,'observers':1<<m,'blue_edges':fixed//2,'black_edges':((1<<m)-fixed)//2,'classes':records}

def hostile():
    # A_0,A_1,A_2,B_0,B_1,B_2; cyclic order inside each cell.
    a=[[0]*6 for _ in range(6)]
    for off in (0,3):
        for u in range(3):a[off+u][off+(u+1)%3]=1
    for u in range(3):
        for v in range(3):a[u][3+v]=int((v-u)%3 in (1,2));a[3+v][u]=1-a[u][3+v]
    valid=[];literal=[]
    for x,y in it.product(range(3),repeat=2):
        p=tuple((x+i)%3 for i in range(3))+tuple(3+(y+i)%3 for i in range(3))
        ok=all(a[p[i]][p[i+1]] for i in range(5))
        require(ok==((y-x)%3 in (0,1)),'native boundary criterion')
        if ok:valid.append((x,y));literal.append(p)
    require(len(set(literal))==6,'six different actual paths')
    def heap(x,y,z):return tuple((xx-yy+zz)%3 for xx,yy,zz in zip(x,y,z))
    witness=((0,1),(0,0),(0,1));out=heap(*witness)
    require(all(t in valid for t in witness) and out not in valid,'heap failure')
    preserving=[]
    for x,y in it.product(range(3),repeat=2):
        if all(a[(u+x)%3][3+(v+y)%3]==a[u][3+v] for u,v in it.product(range(3),repeat=2)):preserving.append((x,y))
    require(preserving==[(0,0),(1,1),(2,2)],'whole-block group')
    opposite=[[a[j][i] for j in range(6)] for i in range(6)]
    require(all(all(opposite[p[i+1]][p[i]] for i in range(5)) for p in literal),'literal converse paths')
    guarded=[(x,y) for x,y in valid if (x-y)%8 in (0,7)]
    require(len(guarded)==5,'synthetic mod8 guard')
    return {'scores':[sum(row) for row in a],'valid_chart_pairs':valid,'actual_paths':literal,'whole_block_stabilizer':preserving,'heap_input':witness,'heap_output':out,'guarded_pairs':guarded,'all_Hamilton_paths':len(hp(a))}

def native_cycle_tests():
    rng=random.Random(430031);total=0;nonempty=0;cycle_rejected=0
    subsets=[(0,),(1,),(2,),(0,1),(0,2),(1,2)]
    for trial in range(72):
        blocks={(i,j):set(rng.choice(subsets)) for i,j in pairs(3)}
        template=[(i,u) for u in range(3) for i in range(3)]
        if trial%2:rng.shuffle(template)
        def arc(i,u,j,v):
            if i==j:return (v-u)%3==1
            if i<j:return (v-u)%3 in blocks[i,j]
            return not arc(j,v,i,u)
        edges=[];internal=True
        for (i,u),(j,v) in zip(template,template[1:]):
            if i==j:internal &= arc(i,u,j,v)
            else:edges.append((i,j,[d for d in range(3) if arc(i,u,j,(v+d)%3)]))
        actual=[]
        for a in it.product(range(3),repeat=3):
            path=[(i,(u+a[i])%3) for i,u in template]
            if all(arc(i,u,j,v) for (i,u),(j,v) in zip(path,path[1:])):actual.append(a)
        expanded=0
        for ds in it.product(*(e[2] for e in edges)):
            potentials={};components=0;ok=True
            graph=[[] for _ in range(3)]
            for (i,j,_),d in zip(edges,ds):graph[i].append((j,d));graph[j].append((i,-d))
            for root in range(3):
                if root in potentials:continue
                components+=1;potentials[root]=0;stack=[root]
                while stack:
                    i=stack.pop()
                    for j,d in graph[i]:
                        v=(potentials[i]+d)%3
                        if j in potentials:ok &= (potentials[j]==v)
                        else:potentials[j]=v;stack.append(j)
            if ok:expanded+=3**components
            else:cycle_rejected+=1
        if not internal:expanded=0
        require(expanded==len(actual),'cycle-potential native compiler')
        if actual:nonempty+=1
        total+=len(actual)
    return {'seed':430031,'cases':72,'accepted_chart_tuples':total,'nonempty_cases':nonempty,'inconsistent_difference_branches':cycle_rejected}

r={'scope':'all Hamilton observer masks at n=3,4,5; one actual six-vertex two-cell tournament; no all-size native completeness claim','observer_censuses':[census(n) for n in (3,4,5)],'native_union_hostile':hostile(),'native_cycle_tests':native_cycle_tests()}
print(json.dumps(r,sort_keys=True,indent=2))
