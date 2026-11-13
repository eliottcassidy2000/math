"""Compile the inherited complete unique-start core grammar, retaining addresses.

No enumeration of binary outer words is used by the target inverse. This small
core bank is a stated input language, not all rooted tournament classes.
"""
from functools import lru_cache
from itertools import product
from pathlib import Path
import argparse,json
HERE=Path(__file__).resolve().parent
SINGLE=(1,1,1,2,1,1)
@lru_cache(None)
def coeff(word,m):
    t=len(word)
    l=(m-1)**t+sum(m**i*(m-1)**(t-i-1) for i,b in enumerate(word) if b)
    d=(m-1)**t+sum((m-1)**i*m**(t-i-1) for i,b in enumerate(word) if not b)
    state={(0,0):1}
    for bit in word:
        new={}
        for (last,started),weight in state.items():
            moves=[((last,started),m-1),((1 if bit else 2,started),1)]
            if started or not bit:moves.append(((last,1),1))
            for key,mult in moves:new[key]=new.get(key,0)+weight*mult
        state=new
    return l,d,sum(v for (last,_),v in state.items() if last!=2)
def advance(x,word):
    h,e,s,p,f,t=x;r=len(word)
    (l1,d1,b1),(l2,d2,b2),(l3,d3,_)=tuple(coeff(word,m) for m in (1,2,3))
    return (l2*s+b1*h-d1*e,l1*h-(0 if r else 1)*e,
            2**r*s+d2*h,2*(l3*t+b2*p-d2*f),l2*p-f,3**r*t+d3*p)
def bank(maxn):
    out={1:[(SINGLE,[])]}
    for n in range(3,maxn+1):
        rows=[]
        for m in [1]+list(range(3,n)):
            t=n-m-1
            words=[()] if not t else [tuple((j>>i)&1 for i in range(t-1))+(1,) for j in range(1<<(t-1))]
            for w in words:
                rows.extend((advance(s,w),address+[''.join(map(str,w))]) for s,address in out[m])
        out[n]=rows
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument('--order',type=int,default=10);a=p.parse_args()
    rows=[]
    for n,entries in bank(a.order).items():
        for rank,(state,address) in enumerate(entries):rows.append(dict(order=n,rank=rank,state=state,address=address))
    (HERE/'core_bank.json').write_text(json.dumps({'max_core_order':a.order,'cores':rows},indent=2)+'\n')
    (HERE/'core_bank.txt').write_text(''.join(' '.join(map(str,(r['order'],r['rank'],*r['state'])))+'\n' for r in rows))
    print('compiled',len(rows),'complete core addresses through',a.order)
if __name__=='__main__':main()
