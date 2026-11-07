#!/usr/bin/env python3
"""Unpruned permutation canonicalization + literal SCC-word quotient audit."""
import argparse,itertools,json,math
from collections import Counter
from pathlib import Path

def need(p,msg):
    if not p:raise RuntimeError(msg)

def hp(r):
    n=len(r);d=[[0]*n for _ in range(1<<n)]
    for u in range(n):d[1<<u][u]=1
    for s in range(1,1<<n):
        for u in range(n):
            if not s>>u&1:continue
            t=s^(1<<u)
            if t:d[s][u]=sum(d[t][v] for v in range(n) if t>>v&1 and r[v]>>u&1)
    return sum(d[-1])

def canonical(r):
    # Every permutation is tested; no degree partition or producer import.
    n=len(r);pairs=list(itertools.combinations(range(n),2));best=None;aut=0
    for p in itertools.permutations(range(n)):
        c=sum(((r[p[u]]>>p[v])&1)<<j for j,(u,v) in enumerate(pairs))
        if best is None or c<best:best,aut=c,1
        elif c==best:aut+=1
    return best,aut

def decode(n,c):
    r=[0]*n
    for k,(u,v) in enumerate(itertools.combinations(range(n),2)):
        if c>>k&1:r[u]|=1<<v
        else:r[v]|=1<<u
    return tuple(r)

def converse(r):return tuple(((1<<len(r))-1)^x^(1<<u) for u,x in enumerate(r))

def strong(r):
    for start in range(len(r)):
        seen=1<<start
        while True:
            nxt=seen
            for u in range(len(r)):
                if seen>>u&1:nxt|=r[u]
            if nxt==seen:break
            seen=nxt
        if seen!=(1<<len(r))-1:return False
    return True

def compositions(total,length):
    if length==1:
        yield (total,);return
    for a in range(total+1):
        for b in compositions(total-a,length-1):yield (a,)+b

parser=argparse.ArgumentParser();parser.add_argument('--packet',required=True)
base=Path(parser.parse_args().packet)
raw=[tuple(map(int,line.split())) for line in (base/'n7_representatives.txt').read_text().splitlines() if line.strip()]
need(len(raw)==456,'atlas row count')
atoms={};subsets=0
for n in range(3,8):
    mats=set()
    for r in raw:
        for vs in itertools.combinations(range(7),n):
            mats.add(tuple(sum(((r[u]>>v)&1)<<j for j,v in enumerate(vs)) for u in vs))
    for r in mats:
        subsets+=1
        h=hp(r)
        if h not in (3,5,9,15,27) or not strong(r):continue
        code,aut=canonical(r);key=(n,code)
        if key not in atoms:
            dual=canonical(converse(r))[0]
            atoms[key]=(h,aut,(n,dual))
inventory=Counter((n,h,aut,key==dual) for key,(h,aut,dual) in atoms.items() for n in [key[0]])
expected=Counter({(3,3,3,True):1,(4,5,1,True):1,(5,9,1,True):2,
 (5,15,3,True):1,(5,15,5,True):1,(6,15,1,False):2,
 (6,27,3,False):2,(7,27,1,False):2,(7,27,1,True):1})
need(inventory==expected,('atom inventory',inventory))

def cores(h):
    if h==1:yield ();return
    for key,(hh,_,_) in atoms.items():
        if h%hh==0:
            for w in cores(h//hh):yield (key,)+w

results={}
singleton=(1,0)
for h in (15,27):
    words={}
    for core in cores(h):
        used=sum(k[0] for k in core)
        if used>14:continue
        for spaces in compositions(14-used,len(core)+1):
            w=[]
            for j,key in enumerate(core):w += [singleton]*spaces[j]+[key]
            w += [singleton]*spaces[-1]
            words[tuple(w)]=math.prod(atoms[key][1] for key in core)
    merged={};sc=0
    for w,aut in words.items():
        dual=tuple(singleton if key==singleton else atoms[key][2] for key in w[::-1])
        need(dual in words,('converse closure',h,w))
        sc+=w==dual
        merged[min(w,dual)]=h//aut
    histogram=Counter(merged.values())
    need(sc==0,('self-converse order14',h))
    need(histogram==({3:5,5:41,15:9} if h==15 else {1:28,9:65,27:12}),('degree histogram',h,histogram))
    results[h]={'ordinary':len(words),'merged':len(merged),'self_converse':sc,'degrees':dict(sorted(histogram.items()))}

r=decode(14,292187853813044375413825338)
need(strong(r),'large witness strong')
need(hp(r)==11390625==15**6,'large witness reverse DP')
print(json.dumps({'status':'PASS','induced_matrices_examined':subsets,
 'canonicalization':'all vertex permutations, no degree pruning',
 'atom_count':len(atoms),'n14':results,'order14_H':hp(r)},indent=2,sort_keys=True))
