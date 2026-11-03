"""Smallest source-block extension scouts; actual strong five-vertex charts."""
from pathlib import Path
from itertools import combinations
from collections import Counter
import json
from transfer import compile_descriptor,SMALL,LETTERS,need,apply3
from reach import closure,word_for,minimum_orders,primes
HERE=Path(__file__).resolve().parent
def chart(n,bits):
    rows=[0]*n;k=0
    for i,j in combinations(range(n),2):
        forward=j==i+1 or (bits>>k&1)
        if j!=i+1:k+=1
        rows[i if forward else j]|=1<<(j if forward else i)
    return tuple(rows)
def strong(rows):
    n=len(rows);reach=[r|1<<i for i,r in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    return all(r==(1<<n)-1 for r in reach)
def descriptor(rows):
    n=len(rows);answer=[]
    for i in reversed(range(n)):
        mask=1
        for j in range(i+2,n):
            if not (rows[i]>>j&1):mask|=1<<(j-i-1)
        answer.append(mask)
    return tuple(answer)
def main():
    bank={};details={}
    for bits in range(64):
        source=chart(5,bits)
        if not strong(source):continue
        name='p5_'+str(bits);desc=descriptor(source);mat=compile_descriptor(desc)
        need(all(x>=0 for row in mat for x in row),'positive matrix')
        need(all(row[2]>=row[0]+row[1] for row in mat),'cone invariant')
        need(all(row[0]+row[1]>=1 for row in mat),'monotone H')
        bank[name]=mat;details[name]={'bits':bits,'source_rows':source,'descriptor':desc,'matrix':mat}
    need(len(bank)==50,'strong chart count')
    fixes={}
    for name,mat in bank.items():
        mats=dict(SMALL,**{name:mat});states,parents,seen=closure(343,mats)
        hs={sum(v[:2]):i for i,v in enumerate(states)}
        got={str(h):word_for(hs[h],parents) for h in (37,43,47,63,133,343) if h in hs}
        if got:fixes[name]=got
    mats=dict(SMALL,**bank)
    states,parents,seen=closure(10001,mats)
    hs={sum(v[:2]):i for i,v in enumerate(states)}
    coefficient_bound=max(2,max((row[2]+row[0]+row[1]-1)//(row[0]+row[1]) for m in bank.values() for row in m))
    for v in states:need(v==(1,0,0) or sum(v[:2])<=v[2]<=coefficient_bound*sum(v[:2]),'expanded cone')
    cert47,_,_=closure(47,mats)
    result={'strong5_charts':details,'single_symbol_repairs':fixes,'C_over_H_bound':coefficient_bound,
            'closed_certificate_47':sorted(cert47),
            'full_bank_10001':{'states':len(states),'values':len(hs),
            'missing':[h for h in range(1,10002,2) if h not in hs],
            'prime_lane_missing':[p for p in primes(10001) if p!=7 and p not in hs],
            'seven_prime_lane_missing_p':[p for p in primes(10001//7) if p!=3 and 7*p not in hs],
            'selected':{str(h):word_for(hs[h],parents) if h in hs else None for h in (37,43,47,63,133,343)}}}
    (HERE/'extension.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'symbol_count':len(bank),'repairs37':[(x,v['37']) for x,v in fixes.items() if '37' in v],
                      'full_bank':result['full_bank_10001']}))
if __name__=='__main__':main()
