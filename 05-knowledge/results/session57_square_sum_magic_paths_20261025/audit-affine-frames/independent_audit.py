"""Independent F3-vector-space audit; no source verifier imports."""
import itertools as it, json
from pathlib import Path
P=Path(__file__).parent

def require(condition, message):
    if not condition: raise RuntimeError(message)

def add(a,b): return tuple((x+y)%3 for x,y in zip(a,b))
def neg(a): return tuple(-x%3 for x in a)
def sub(a,b): return add(a,neg(b))
def canonical(B,D): return tuple(sorted((tuple(sorted(B)),tuple(sorted(D)))))

def census(n,C):
    V=list(it.product(range(3),repeat=n)); zero=(0,)*n
    U={v for v in V if not any(v[2:])}; A={add(C,u) for u in U}
    directions=sorted({tuple(sorted((zero,u,neg(u)))) for u in U-{zero}})
    frames=set()
    for L,M in it.combinations(directions,2):
        baseline={}
        for g in V:
            B={add(g,l) for l in L}; D={add(sub(C,g),m) for m in M}
            require({add(b,d) for b in B for d in D}==A,'cross sums')
            k=canonical(B,D); require(k not in baseline.values(),'injectivity')
            baseline[g]=k
            if not B&D: frames.add(k)
        for h in U:
            lh=next(l for l in L if sub(h,l) in M)
            for g in V:
                gp=add(g,lh); Cp=add(C,h)
                require(canonical({add(gp,l) for l in L},{add(sub(Cp,gp),m) for m in M})==baseline[g],'all direction C changes')
    require(len(frames)==6*(3**n-9),'frame formula')
    u=next(x for x in directions[0] if x!=zero);v=next(x for x in directions[1] if x!=zero)
    B0=set(directions[0]);D0={C,add(C,v),add(C,sub(u,v))}
    require(add(add(*list(D0)[:2]),list(D0)[2])!=zero,'non AP')
    require(sum(not ({add(g,b) for b in B0}&{sub(d,g) for d in D0}) for g in V)==3**n-9,'non AP all gauges')
    return frames

counts=[]
for n in (2,3,4):
    for C in [(0,)*n]+([(0,0)+(1,)*(n-2)] if n>2 else []):
        frames=census(n,C); counts.append({'n':n,'C':C,'frames':len(frames)})
# Independent complete 9! search, testing all eight magic lines.
A=list(it.product(range(3),repeat=2)); zero=(0,0)
lines=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
magic=set()
for p in it.permutations(A):
    sums=[add(add(p[i],p[j]),p[k]) for i,j,k in lines]
    if len(set(sums))==1: magic.add(p)
require(len(magic)==432,'all magic arrangements, not merely explicit family')
require(all(add(add(p[0],p[1]),p[2])==zero for p in magic),'magic constants necessarily zero')
# Every path fixes its unordered alternating bipartition. Explicitly verify injectivity.
paths=set();incidences=0
for B,D in frames:
    for b in it.permutations(B):
        for d in it.permutations(D):
            path=tuple(x for pair in zip(b,d) for x in pair)
            path=min(path,path[::-1]); incidences+=1
            require(path not in paths,'cross frame or within frame duplicated path')
            paths.add(path)
require(len(paths)==incidences==36*len(frames),'incidence injection')
result={'status':'PASS','coordinate_counts':counts,'complete_9_factorial_magic_count':len(magic),'distinct_paths_last_census':len(paths),'path_incidences_last_census':incidences}
(P/'independent-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
