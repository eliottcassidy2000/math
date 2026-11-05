#!/usr/bin/env python3
"""Independent exhaustive labelled-tournament audit, with no repo imports."""
import argparse
import itertools
import json
import math
from collections import Counter
from pathlib import Path


def adjacency(n, code):
    out = [0] * n
    for bit, (a,b) in enumerate(itertools.combinations(range(n), 2)):
        u,v = (a,b) if (code >> bit) & 1 else (b,a)
        out[u] |= 1 << v
    return out


def hamilton_paths(out):
    n=len(out)
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        remaining=((1<<n)-1)^mask
        for i, count in enumerate(dp[mask]):
            if not count: continue
            todo=out[i]&remaining
            while todo:
                bit=todo&-todo
                todo-=bit
                dp[mask|bit][bit.bit_length()-1]+=count
    return sum(dp[-1])


def triangle_count(out):
    return sum(bool(((out[a]>>b)&1) and ((out[b]>>c)&1) and ((out[c]>>a)&1)) +
               bool(((out[b]>>a)&1) and ((out[a]>>c)&1) and ((out[c]>>b)&1))
               for a,b,c in itertools.combinations(range(len(out)),3))


def automorphisms(out):
    # A permutation preserves score classes; enumerate only those possibilities.
    by_score={}
    for i,row in enumerate(out): by_score.setdefault(row.bit_count(),[]).append(i)
    classes=list(by_score.values())
    count=0
    for choices in itertools.product(*(itertools.permutations(c) for c in classes)):
        p=list(range(len(out)))
        for c,q in zip(classes,choices):
            for i,j in zip(c,q): p[i]=j
        if all(bool((out[i]>>j)&1)==bool((out[p[i]]>>p[j])&1)
               for i in range(len(out)) for j in range(i+1,len(out))):
            count+=1
    return count


def isprime(k):
    return k>=2 and all(k%d for d in range(2, math.isqrt(k)+1))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--max-n', type=int, default=6)
    ap.add_argument('--output', type=Path)
    args=ap.parse_args()
    result={'method':'All labelled tournaments, exact subset DP, score-preserving permutation automorphism count; independent of repository code.','orders':[]}
    for n in range(1,args.max_n+1):
        spectrum=Counter()
        primes=Counter()
        prime_aut=Counter()
        selected_aut=Counter()
        tested=0
        for code in range(1<<(n*(n-1)//2)):
            out=adjacency(n,code)
            h=hamilton_paths(out)
            c3=triangle_count(out)
            assert h&1, (n,code,h)
            assert h>=1+2*c3, (n,code,h,c3)
            spectrum[h]+=1
            if isprime(h) or h in (1,21):
                aut=automorphisms(out)
                assert h%aut==0, (n,code,h,aut)
                tested+=1
                if isprime(h):
                    primes[h]+=1
                    prime_aut[h,aut]+=1
                    if h>=5: assert aut==1,(n,code,h,aut)
                if h in (1,3,7,21):
                    selected_aut[h,aut]+=1
                    if h==1: assert aut==1,(n,code,h,aut)
                    if h==3: assert aut==3,(n,code,h,aut)
                    if h==21: assert aut in (1,3),(n,code,h,aut)
            if n>=3 and len(set(row.bit_count() for row in out))==1:
                assert c3==n*(n*n-1)//24, (n,code,c3)
        row={'n':n,'labelled_tournaments':sum(spectrum.values()),'full_spectrum':dict(sorted(spectrum.items())),
             'automorphism_audited':tested,'prime_H_labelled_counts':dict(sorted(primes.items())),
             'prime_H_automorphism_orders':{f'H={h},Aut={a}':v for (h,a),v in sorted(prime_aut.items())},
             'selected_H_automorphism_orders':{f'H={h},Aut={a}':v for (h,a),v in sorted(selected_aut.items())}}
        result['orders'].append(row)
        print(json.dumps(row),flush=True)
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__': main()
