"""Exact end-rooted forward growth and independent presentation controls."""
from pathlib import Path
from itertools import permutations, product
from collections import Counter, defaultdict
from functools import lru_cache
from math import factorial
import json, gzip, sys, time
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'modular_centers_collatz_pell_20261111'))
import content_checks as cc
import tournament_checks as tc

def need(ok,message):
    if not ok: raise RuntimeError(message)

def save(name,data):
    tag='' if __debug__ else '.optimized'
    (HERE/(name+tag+'.json.gz')).write_bytes(gzip.compress(json.dumps(data,separators=(',',':')).encode(),mtime=0))

def relabel(rows,order):
    return tuple(sum(((rows[a]>>b)&1)<<j for j,b in enumerate(order)) for a in order)

@lru_cache(None)
def rooted(rows,r):
    n=len(rows)
    signatures={i:(rows[i].bit_count(),(rows[i]>>r)&1) for i in range(n) if i!=r}
    cells=[[i for i,s in signatures.items() if s==sig] for sig in sorted(set(signatures.values()))]
    best=None; aut=0
    for blocks in product(*(permutations(cell) for cell in cells)):
        order=sum(blocks,())+(r,)
        candidate=relabel(rows,order)
        if best is None or candidate<best:best=candidate;aut=1
        elif candidate==best:aut+=1
    return best,aut

@lru_cache(None)
def ordinary(rows):return cc.canonical(rows)

def endpoints(rows):
    n=len(rows);end=[[0]*n for _ in range(1<<n)]
    for v in range(n):end[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not end[mask][v]:continue
            available=rows[v]&~mask
            while available:
                bit=available&-available;available-=bit
                end[mask|bit][bit.bit_length()-1]+=end[mask][v]
    return end[-1]

def tiling(n,bits):
    rows=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            a,b=(i,j) if j==i+1 or bits>>k&1 else (j,i)
            rows[a]|=1<<b
            if j>i+1:k+=1
    return tuple(rows)

def reflection(n,bits):
    pairs=[(i,j) for i in range(n) for j in range(i+2,n)]
    index={p:k for k,p in enumerate(pairs)}
    return sum(((bits>>k)&1)<<index[n-1-j,n-1-i] for k,(i,j) in enumerate(pairs))

def run(max_n=7):
    weights={(0,):1};levels=[];records=[];transitions=[];hostile=None
    expect=[1,1,2,4,12,56,456,6880]
    for n in range(1,max_n+1):
        t0=time.time();codes={ordinary(T) for T in weights}
        need(len(codes)==expect[n-1],'ordinary census')
        grouped=defaultdict(int);root_counts=Counter();strong_roots=0
        for T,mass in sorted(weights.items()):
            code=ordinary(T);key,aut=rooted(T,n-1);ends=endpoints(T)
            need(key==T and ends[-1]==aut*mass,'rooted orbit mass from independent subset paths')
            if n<=5:
                direct=sum(p[-1]==n-1 and all(T[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(range(n)))
                need(direct==ends[-1],'independent permutation path count')
            grouped[code]+=mass;root_counts[code]+=1
            # Every vertex can end a Hamiltonian path exactly when T is strong.
            strong_roots+=int(all(ends))
            records.append(dict(n=n,rows=T,ordinary_code=code,aut_root=aut,end_paths=ends[-1],mass=mass))
        for code,mass in grouped.items():
            rows=cc.decode(n,code);aut=tc.automorphism_count(rows)
            need(sum(endpoints(rows))==aut*mass,'ordinary fiber mass')
        need(sum(weights.values())==1<<((n-1)*(n-2)//2),'all tiling mass')
        converse={code:ordinary(tuple(((1<<n)-1)^(1<<i)^row for i,row in enumerate(cc.decode(n,code)))) for code in codes}
        selfc=sum(code==opp for code,opp in converse.items())
        level=dict(n=n,ordinary_types=len(codes),end_rooted_types=len(weights),strong_rooted_types=strong_roots,
                   merged_types=(len(codes)+selfc)//2,self_converse_types=selfc,tilings=sum(weights.values()))
        if n<=6:
            literal=Counter();color=defaultdict(dict);profiles=defaultdict(list)
            for bits in range(1<<((n-1)*(n-2)//2)):
                T=tiling(n,bits);key=rooted(T,n-1)[0];literal[key]+=1
                blue=reflection(n,bits)==bits
                color[key].setdefault(blue,bits)
                if hostile is None and len(color[key])==2:
                    hostile=dict(n=n,rooted_rows=key,blue_bits=color[key][True],black_bits=color[key][False])
            need(dict(literal)==weights,'complete tiling enumeration versus rooted transfer')
            level['literal_tilings_checked']=sum(literal.values())
            level['rooted_states_with_both_colors']=sum(len(c)==2 for c in color.values())
        levels.append(level);print(level,'seconds',round(time.time()-t0,2),flush=True)
        if n<max_n:
            following=Counter()
            for T,mass in weights.items():
                counts=Counter(rooted(tuple(cc.child(T,mask)),n)[0] for mask in range(1<<(n-1)))
                need(sum(counts.values())==1<<(n-1),'constant lift degree')
                for U,multiplicity in counts.items():
                    following[U]+=mass*multiplicity
                    transitions.append(dict(n=n,parent=T,child=U,multiplicity=multiplicity))
            weights=dict(following)
            rooted.cache_clear();ordinary.cache_clear()
    save('rooted_growth',dict(scope=f'Every end-rooted tournament type through {max_n}; all forward attachments; complete tiling controls through6',levels=levels,records=records,transitions=transitions,color_hostile=hostile))
    print('color hostile',hostile,flush=True)

if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 7)
