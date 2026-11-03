#!/usr/bin/env python3
"""Named adaptive-root updates, actual reset, and minimal profile nonclosure.
Complete ordinary-class augmentation through six; no repository imports.
"""
from collections import Counter,deque
from functools import lru_cache
from itertools import combinations, permutations, product
from math import factorial
import json


def need(ok,why):
    if not ok:raise RuntimeError(why)


def add(T,v):
    n=len(T)
    return tuple(row|((1<<n) if i==v else 0) for i,row in enumerate(T))+(((1<<n)-1)^(1<<v),)


def induced(T,keep):
    return tuple(sum(((T[u]>>v)&1)<<j for j,v in enumerate(keep)) for u in keep)


def strong(T):
    n=len(T)
    if not n:return False
    R=[x|1<<i for i,x in enumerate(T)]
    for k in range(n):
        for i in range(n):
            if R[i]>>k&1:R[i]|=R[k]
    return all(x==(1<<n)-1 for x in R)


@lru_cache(None)
def table(T):
    n=len(T);dp=[[0]*n for _ in range(1<<n)];H=[0]*(1<<n);H[0]=1
    for i in range(n):dp[1<<i][i]=1
    for M in range(1,1<<n):
        for v,c in enumerate(dp[M]):
            if c:
                choices=T[v]&~M
                while choices:
                    bit=choices&-choices;choices-=bit
                    dp[M|bit][bit.bit_length()-1]+=c
        H[M]=sum(dp[M])
    return H,dp


@lru_cache(None)
def profiles(T):
    n=len(T);full=(1<<n)-1
    H,ends=table(T);_,starts=table(tuple(full^r^(1<<i) for i,r in enumerate(T)))
    return tuple((H[-1],starts[-1][v],H[full^(1<<v)],
                  sum(ends[M][v]*H[full^M] for M in range(1,full+1))) for v in range(n))


def keep(z):
    h,s,q,F=z;return h+F-s,q,q,2*F-q


def new(z):
    h,s,q,F=z;return h+F-s,h-s,h,h+F


@lru_cache(None)
def from_start(T,u):
    n=len(T);dp=[[0]*n for _ in range(1<<n)];dp[1<<u][u]=1
    for M in range(1,1<<n):
        for v,c in enumerate(dp[M]):
            if c:
                choices=T[v]&~M
                while choices:
                    bit=choices&-choices;choices-=bit
                    dp[M|bit][bit.bit_length()-1]+=c
    return dp


@lru_cache(None)
def joint(T,u,v):
    """A_uv: P starts u ends v. B_uv: P ends v,Q starts u.
    G_uv: P ends u,Q ends v,third slot arbitrary. All cover T.
    """
    need(u!=v,'distinct joint ports')
    n=len(T);full=(1<<n)-1;H,ends=table(T)
    _,starts=table(tuple(full^r^(1<<i) for i,r in enumerate(T)))
    fixed=from_start(T,u)
    A=sum(fixed[M][v]*H[full^M] for M in range(1,full+1))
    B=sum(ends[M][v]*starts[full^M][u] for M in range(1,full+1))
    G=0
    for M in range(1,full+1):
        if not ends[M][u]:continue
        available=full^M;Q=available
        while Q:
            G+=ends[M][u]*ends[Q][v]*H[available^Q]
            Q=(Q-1)&available
    return A,B,G


def encode(T,p):
    return sum(((T[p[i]]>>p[j])&1)<<k for k,(i,j) in enumerate(combinations(range(len(T)),2)))


def decode(n,code):
    T=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        a,b=(i,j) if code>>k&1 else (j,i);T[a]|=1<<b
    return tuple(T)


def orders(T):
    cells={}
    for i,r in enumerate(T):cells.setdefault(r.bit_count(),[]).append(i)
    for blocks in product(*(permutations(cells[d]) for d in sorted(cells))):
        yield sum(blocks,())


def canonical(T):return min(encode(T,p) for p in orders(T))


def inventories(N):
    layers=[{()}]
    for n in range(1,N+1):
        codes=set()
        for T in layers[-1]:
            for cut in range(1<<(n-1)):
                R=tuple(row|((1<<(n-1)) if cut>>i&1 else 0) for i,row in enumerate(T))+(((1<<(n-1))-1)^cut,)
                codes.add(canonical(R))
        layers.append({decode(n,c) for c in codes})
    need([len(x) for x in layers]==[1,1,1,2,4,12,56][:N+1],'ordinary class universe')
    return layers


def legal_last(T):
    n=len(T)
    return [v for v,row in enumerate(T) if n-1-row.bit_count()==1
            and strong(induced(T,tuple(i for i in range(n) if i!=v)))]


def accepted(T):
    if not strong(T):return False
    while len(T)>3:
        xs=legal_last(T);need(len(xs)<=1,'unique inverse')
        if not xs:return False
        T=induced(T,tuple(i for i in range(len(T)) if i!=xs[0]))
    return len(T)==3 and strong(T)


def hist(T,v):
    return Counter(h+F-s for h,s,q,F in profiles(add(T,v)))


def history(word):
    T=(2,4,1)
    for root in word:T=add(T,root)
    return T


def limited(seeds,bound):
    todo=deque(seeds);seen={z:(i,()) for i,z in enumerate(seeds)}
    while todo:
        z=todo.popleft()
        for letter,fun in (('R',keep),('N',new)):
            nz=fun(z)
            need(nz[0]>z[0],'strict growth in actual restricted language')
            if nz[0]<=bound and nz not in seen:
                i,w=seen[z];seen[nz]=(i,w+(letter,));todo.append(nz)
    return seen


def main():
    layers=inventories(6);controls=Counter();minimal=[]
    for n in range(3,7):
        seen={};collisions=[];accepted_count=0
        for T in sorted(layers[n]):
            if not strong(T):continue
            need(len(legal_last(T))<=1,'all-strong unique inverse')
            accepted_count+=accepted(T)
            for v,z in enumerate(profiles(T)):
                h,s,q,F=z;delta=F-h-q+s;indegree=n-1-T[v].bit_count()
                need(delta>=indegree*(indegree-1),'positive mixed-endpoint defect')
                need((delta==0)==(indegree==1),'strong degree-one profile detector')
                need((delta==0 and q>s)==(v in legal_last(T)),'profile legal-inverse detector')
                child=add(T,v)
                need(profiles(child)[v]==keep(z),'same root update')
                need(profiles(child)[n]==new(z),'new root update')
                reset=keep(new(z));h,s,q,F=z
                rz=keep(z);nz=new(z)
                need(rz[3]-rz[0]-rz[2]+rz[1]==delta,'keep preserves defect')
                need(nz[3]-nz[0]-nz[2]+nz[1]==0,'new annihilates defect')
                need(reset==(h+2*F,h,h,h+2*F),'exact balanced reset')
                need(reset[0]>z[0],'reset growth')
                controls['all_strong_root_update_checks']+=1
                for w in range(n):
                    if w==v:continue
                    A,B,G=joint(T,v,w);Awv=joint(T,w,v)[0]
                    oldw=profiles(T)[w];neww=profiles(child)[w]
                    need(neww[1]==Awv,'old root start update')
                    need(neww[3]==2*oldw[3]-A-B+2*G,'old root cover update')
                    keep_vertices=tuple(i for i in range(n) if i!=w)
                    deleted=induced(T,keep_vertices);v_deleted=keep_vertices.index(v)
                    dz=profiles(deleted)[v_deleted]
                    need(neww[2]==keep(dz)[0],'old root deletion update')
                    controls['all_strong_ordered_old_root_pair_checks']+=1
                signature=tuple(sorted(hist(T,v).items()))
                if z in seen and seen[z][0]!=signature:
                    collisions.append((z,seen[z],(signature,T,v)))
                seen.setdefault(z,(signature,T,v))
        minimal.append({'n':n,'profile_collision':bool(collisions),'accepted_classes':accepted_count})
        if n<6:need(not collisions,'earlier profile nonclosure')
    need([r['accepted_classes'] for r in minimal]==[1,1,4,20],'adaptive class counts')
    A=(40,17,3,6,13,30);B=(34,24,3,5,13,30)
    z=(37,5,13,71)
    need(profiles(A)[2]==profiles(B)[3]==z,'profile witness')
    need(accepted(A) and accepted(B),'actual adaptive membership')
    need(canonical(A)!=canonical(B),'ordinary source distinction')
    values_A=[h+F-s for h,s,q,F in profiles(add(A,2))]
    values_B=[h+F-s for h,s,q,F in profiles(add(B,3))]
    need(values_A==[287,321,219,329,275,233,179] and values_B==[289,307,339,219,271,227,179],'root map')
    # Explicit histories and their named transport into the displayed gauge.
    for word,transport,target,root in (((0,1,2),(3,1,0,2,4,5),A,3),
                                        ((0,2,1),(3,0,1,2,4,5),B,0)):
        T=history(word);R=[0]*6
        for i in range(6):
            for j in range(6):
                if T[i]>>j&1:R[transport[i]]|=1<<transport[j]
        need(tuple(R)==target,'named adaptive transport')
        need(profiles(T)[root]==z,'transported root')
    # Independent full vertex-word counts on final eight-vertex witnesses.
    literal={}
    for name,T,v,w in (('A233',A,2,5),('A287',A,2,0),('B307',B,3,1)):
        U=add(add(T,v),w);h=table(U)[0][-1]
        direct=sum(all(U[p[i]]>>p[i+1]&1 for i in range(7)) for p in permutations(range(8)))
        need(direct==h,'literal final target')
        literal[name]=h
    # Closed two-operation language with actual T4 seeds, not arbitrary profiles.
    T4=history((0,));seeds=profiles(T4)
    c3_states=limited(((3,1,1,3),),199)
    t4_states=limited(seeds,199)
    c3_hits={z[0] for z in c3_states};t4_hits={3}|{z[0] for z in t4_states}
    need(11 not in c3_hits and 13 not in c3_hits,'C3 keep/new obstruction')
    need(39 not in t4_hits and 67 not in t4_hits,'all-root T4 keep/new obstruction')
    # Replay every retained finite-state witness as an actual named graph.
    for z,(seed,w) in t4_states.items():
        T=T4;root=seed
        for letter in w:
            n=len(T);T=add(T,root)
            if letter=='N':root=n
        need(profiles(T)[root]==z,'actual finite-state recipe')
        controls['restricted_actual_recipe_replays']+=1
    # Both nonadaptive strong order-five classes have H15, with different stopping mechanisms.
    stop=[]
    for T in sorted(layers[5]):
        if strong(T) and not accepted(T):
            stop.append({'H':table(T)[0][-1], 'indegree_one_vertices':[i for i,r in enumerate(T) if 4-r.bit_count()==1], 'rows':T})
    need(len(stop)==2 and all(s['H']==15 for s in stop),'earliest class hostiles')
    return {'status':'PROVED exact rooted updates/reset; finite exact minimal profile hostile and restricted-language exclusions',
            'controls':dict(controls),'minimal_profile_scan':minimal,
            'same_profile':(37,5,13,71),
            'A':{'rows':A,'root':2,'next_root_values':values_A,'history':[0,1,2]},
            'B':{'rows':B,'root':3,'next_root_values':values_B,'history':[0,2,1]},
            'literal_targets':literal,'order5_nonadaptive_strong':stop,
            'C3_keep_new_first_missing':[h for h in range(3,200,2) if h not in c3_hits][:15],
            'T4_all_roots_keep_new_missing_to199':[h for h in range(3,200,2) if h not in t4_hits],
            'T4_restricted_states':len(t4_states)}


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
