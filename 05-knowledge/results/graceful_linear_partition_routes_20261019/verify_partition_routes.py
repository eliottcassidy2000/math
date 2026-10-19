"""Explicit cycle-arborescence routing for full-rank staircase partition matroids."""
from collections import deque,defaultdict
from itertools import product,combinations
from random import Random
import json

CHECKS=MOVES=0

def check(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(msg)

def support(A):return {p for p in A if p is not None}
def block_bound(b):return 1 if b<=2 else 2*b-3

def route(classes,source,target):
    """Each element is its actual index in its color list; list values are parts."""
    global MOVES
    r=len(classes);current=list(source);path=[tuple(source)]
    S=[support(A) for A in classes]
    representatives=[]
    for A in classes:
        rr={}
        for j,p in enumerate(A):
            if p is not None:rr.setdefault(p,j)
        representatives.append(rr)
    check(all(len(A)>=r-i for i,A in enumerate(S)),'staircase')
    check(len(set().union(*S))==r,'fully tight rank')
    def valid(st):
        pp=[classes[i][a] for i,a in enumerate(st)]
        return None not in pp and len(set(pp))==r
    check(valid(source) and valid(target),'endpoint bases')
    def add(new):
        old=tuple(current);new=tuple(new)
        if new==old:return
        check(sum(a!=b for a,b in zip(old,new))<=2,'width two')
        check(valid(new),'all parts distinct')
        current[:]=new;path.append(new)
    if r<=2:
        add(target);MOVES+=len(path)-1;return path
    suffix=[set() for _ in range(r+1)]
    for i in range(r-1,-1,-1):suffix[i]=suffix[i+1]|S[i]
    tight=[i for i in range(r+1) if len(suffix[i])==r-i]
    check(tight[0]==0 and tight[-1]==r,'complete tight chain')
    total_bound=0
    for lo,hi in reversed(list(zip(tight,tight[1:]))):
        b=hi-lo;total_bound+=block_bound(b)
        parts=suffix[lo]-suffix[hi]
        ids=list(range(lo,hi));labels={classes[i][target[i]]:i-lo for i in ids}
        check(set(labels)==parts,'target block labels')
        avail=[{labels[p] for p in S[i]&parts} for i in ids]
        check(len(avail[0])==b,'full buffer class')
        for mask in range(1,(1<<(b-1)) if b<=8 else 1):
            J=[i+1 for i in range(b-1) if mask>>i&1]
            check(len(set().union(*(avail[j] for j in J)))>len(J),'bounded strict Hall audit')
        sigma=[labels[classes[i][current[i]]] for i in ids]
        before=len(path)
        if sigma==list(range(b)):
            wrong=[i for i in ids if current[i]!=target[i]]
            for z in range(0,len(wrong),2):
                new=current[:]
                for i in wrong[z:z+2]:new[i]=target[i]
                add(new)
        else:
            cycles=[];which=[None]*b
            for j in range(b):
                if which[j] is not None:continue
                cyc=[];k=j
                while which[k] is None:
                    which[k]=len(cycles);cyc.append(k);k=sigma[k]
                cycles.append(cyc)
            root=which[0];children=defaultdict(list);seen={root};queue=deque([root])
            while queue:
                C=queue.popleft()
                for a in cycles[C]:
                    for j in range(1,b):
                        D=which[j]
                        if D not in seen and a in avail[j]:
                            seen.add(D);queue.append(D);children[a].append((D,j))
            check(len(seen)==len(cycles),'cycle arborescence reachable')
            inv={v:k for k,v in labels.items()}
            def move(j):
                global MOVES
                i=lo+j;hole=classes[lo][current[lo]];oldpart=classes[i][current[i]]
                new=current[:]
                new[i]=target[i] if hole==classes[i][target[i]] else representatives[i][hole]
                new[lo]=target[lo] if oldpart==classes[lo][target[lo]] else representatives[lo][oldpart]
                add(new)
            events=[('enter',root,0)]
            while events:
                kind,u,v=events.pop()
                if kind=='move':
                    move(u);continue
                if kind=='check':
                    check(labels[classes[lo][current[lo]]]==u,'hole at scheduled position');continue
                C,entry=u,v
                if C!=root:move(entry)
                sequence=[];a=sigma[entry]
                while True:
                    sequence.append(a)
                    if a==entry:break
                    a=sigma[a]
                for a in reversed(sequence):
                    if a!=entry or C!=root:events.append(('move',a,None))
                    for D,j in reversed(children[a]):events.append(('enter',D,j))
                    events.append(('check',a,None))
            exact=b-1+len(cycles)-1
            check(len(path)-before==exact,'exact arborescence route cost')
        check(current[lo:hi]==list(target[lo:hi]),'exact target realizations')
        check(len(path)-before<=block_bound(b),'linear block bound')
    check(tuple(current)==tuple(target),'global target')
    check(len(path)-1<=total_bound<=2*r-3,'uniform linear bound')
    MOVES+=len(path)-1
    return path

def states(classes):
    ans=[]
    def rec(i,used,st):
        if i<0:ans.append(tuple(st));return
        for a,p in enumerate(classes[i]):
            if p is not None and p not in used:
                st[i]=a;rec(i-1,used|{p},st)
    rec(len(classes)-1,set(),[0]*len(classes));return ans

def random_state(classes,rng):
    used=set();st=[None]*len(classes)
    for i in range(len(classes)-1,-1,-1):
        choices=[a for a,p in enumerate(classes[i]) if p is not None and p not in used]
        st[i]=rng.choice(choices);used.add(classes[i][st[i]])
    return tuple(st)

def main():
    rng=Random(20261019);systems=routes=0;maxratio=0
    for r in range(1,17):
        for trial in range(60):
            A=[]
            for i in range(r):
                v=rng.sample(range(r),rng.randrange(r-i,r+1))
                if trial%3==0:v+=[None]
                if trial%4==0:v+=[rng.choice(v[:len(v)-(trial%3==0)])]
                A.append(tuple(v))
            pairs=[(random_state(A,rng),random_state(A,rng)) for _ in range(6)]
            if r<=5:
                st=states(A)
                if len(st)<=60:pairs.extend((x,y) for x in st for y in st)
            for x,y in pairs:
                pp=route(A,x,y);routes+=1
                maxratio=max(maxratio,(len(pp)-1)/r)
            systems+=1
    hostile=((0,1,2,3),(0,1,3),(1,2),(2,0))
    source=(3,0,0,0);target=(3,1,1,1)
    hp=route(hostile,source,target)
    check(len(hp)-1==4,'known five-state hostile route')
    print(json.dumps({'status':'PASS','checks':CHECKS,'systems':systems,'routes':routes,'validated_moves':MOVES,'largest_rank':16,'maximum_route_to_rank_ratio':maxratio,'target_preservation_hostile_length':len(hp)-1,'scope':'Explicit polynomial cycle-arborescence routing, not state-graph BFS. Random finite controls plus complete endpoint pairs in bounded small systems; universal bounds are proved separately.'},indent=2))

if __name__=='__main__':main()
