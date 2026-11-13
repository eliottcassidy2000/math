"""Exact word/address controls for the proved physical retained-chain grammar."""
from collections import Counter, defaultdict
from math import gcd
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def forward(pair,k):
    q,s=pair; S=(1<<(k-1))*(3*q+2*s)
    return q+3*S,S


def inverse(pair):
    Q,S=pair
    v=S&-S; k=v.bit_length()
    q=Q-3*S
    numerator=S//v-3*q
    if q<=0 or numerator<=0 or numerator%4:
        return None
    return (q,numerator//2),k


def decode(pair,seed):
    reverse=[]
    while pair!=seed:
        old=inverse(pair)
        if old is None:
            return None
        previous,k=old
        need(sum(previous)<sum(pair),'inverse strict descent')
        reverse.append(k);pair=previous
    return tuple(reversed(reverse))


def main():
    cap=24;expected=[0]*(cap+1);expected[0]=1
    for n in range(4,cap+1):
        expected[n]=expected[n-1]+expected[n-4]
    results=[]
    for seed in ((135,38),(215,74),(87,20)):
        by_size=Counter();seen={};totals=defaultdict(list)
        def walk(pair,word,cost):
            need(pair not in seen,'two words share a rooted profile')
            seen[pair]=word;by_size[cost]+=1
            need(decode(pair,seed)==word,'exact address decoding')
            q,s=pair
            need(q%2==1 and s%2==0 and gcd(q,s)==1,'primitive chart guard')
            totals[q+3*s].append((word,pair,cost))
            for k in range(2,cap-cost-1):
                walk(forward(pair,k),word+(k,),cost+k+2)
        walk(seed,(),0)
        need([by_size[n] for n in range(cap+1)]==expected,'exact rational series')
        results.append(dict(seed=seed,words=len(seen),counts_by_added_vertices=expected,
                            scalar_collisions=[dict(H=H,addresses=rows) for H,rows in sorted(totals.items()) if len(rows)>1]))
    out=dict(status='PASS: recurrence-level exact controls of an independently proved physical grammar',
             universe=dict(seeds=3,max_added_vertices=cap,letters='every integer k>=2 with cost k+2'),results=results)
    (HERE/'chain_language_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'PASS','words_per_seed':len(seen),
                      'scalar_collision_totals':[len(r['scalar_collisions']) for r in results]}))


if __name__=='__main__':
    main()
