#!/usr/bin/env python3
"""Focused H=63 edge surgery; literal cores, no census assumption in controls."""
from itertools import combinations, permutations
from collections import Counter
import json
from pathlib import Path

CORES=((16,5,17,7,74,31,175,63),(2,16,3,23,165,79,31,111))

def inv(rows):
    n=len(rows); full=(1<<n)-1
    end=[[0]*n for _ in range(full+1)]
    H=[0]*(full+1); H[0]=1
    for v in range(n): end[1<<v][v]=1
    for mask in range(1,full+1):
        for v in range(n):
            if not(mask>>v&1): continue
            z=mask^(1<<v)
            if z: end[mask][v]=sum(end[z][u] for u in range(n) if z>>u&1 and rows[u]>>v&1)
        H[mask]=sum(end[mask])
    starts=[0]*n
    for v in range(n):
        z=full^(1<<v)
        # reversal of paths identifies starts through the converse DP.
        starts[v]=sum(path_count_from(rows,z,w) for w in range(n) if z>>w&1 and rows[v]>>w&1)
    rooted=[]
    for v in range(n):
        F=sum(end[mask][v]*H[full^mask] for mask in range(1,full+1) if mask>>v&1)
        rooted.append([H[full],starts[v],H[full^(1<<v)],F])
    return H[full],rooted

def path_count_from(rows,mask,start):
    # Only used for eight tiny rooted controls, memoized per rows.
    if mask==(1<<start): return 1
    z=mask^(1<<start)
    return sum(path_count_from(rows,z,w) for w in range(len(rows)) if z>>w&1 and rows[start]>>w&1)

def hp(rows):
    n=len(rows); full=(1<<n)-1; dp=[[0]*n for _ in range(full+1)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,full+1):
        for v in range(n):
            c=dp[mask][v]
            if not c: continue
            avail=rows[v]&(full^mask)
            while avail:
                b=avail&-avail; avail-=b
                dp[mask|b][b.bit_length()-1]+=c
    return sum(dp[full])

def flip(rows,edges):
    r=list(rows)
    for a,b in edges: r[a]^=1<<b; r[b]^=1<<a
    return tuple(r)

def main():
    pairs=list(combinations(range(8),2)); records=[]
    for index,T in enumerate(CORES):
        assert hp(T)==63
        H,states=inv(T)
        direct=sum(all(T[p[i]]>>p[i+1]&1 for i in range(7)) for p in permutations(range(8)))
        assert direct==63
        sig=sorted(r.bit_count()for r in T)
        bysize={}; witnesses=[]
        for size in (1,2,3):
            count=Counter(); got=[]
            for E in combinations(pairs,size):
                U=flip(T,E)
                if hp(U)!=63: continue
                newclass=sorted(r.bit_count()for r in U)!=sig
                count['all']+=1; count['opposite_class' if newclass else 'same_class']+=1
                got.append({'edges':E,'rows':U,'opposite_class':newclass})
            bysize[size]=dict(count)
            witnesses.extend(got)
        records.append({'index':index,'rows':T,'states':states,'neutral':bysize,'witnesses':witnesses})
    target=Path(__file__).with_name('probe_receipt.json')
    target.write_text(json.dumps({'universe':'all subsets of1,2,3pairs on each literal8core','records':records},indent=2)+'\n')
    print(json.dumps([{k:v for k,v in R.items()if k!='witnesses'}for R in records],indent=2))

if __name__=='__main__': main()
