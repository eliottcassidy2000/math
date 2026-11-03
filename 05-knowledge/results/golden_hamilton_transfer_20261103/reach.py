"""Complete bounded-H reachability, not a finite-order scan."""
from collections import deque
from pathlib import Path
from itertools import combinations
import argparse, hashlib, json, heapq
from transfer import SMALL,LETTERS,apply3,need,vector3,targets,subset_h
HERE=Path(__file__).resolve().parent
def primes(n):
    sieve=bytearray(b'\1')*(n+1)
    if n>=0:sieve[0]=0
    if n>=1:sieve[1]=0
    for p in range(2,int(n**.5)+1):
        if sieve[p]:sieve[p*p::p]=b'\0'*(((n-p*p)//p)+1)
    return [p for p in range(3,n+1,2) if sieve[p]]
def closure(limit,bank=SMALL):
    start=(1,0,0);seen={start:0};states=[start];parents=[(-1,"")];todo=deque([0])
    while todo:
        idx=todo.popleft();v=states[idx]
        for letter,mat in bank.items():
            w=apply3(v,mat);h=w[0]+w[1]
            if h>limit or w in seen:continue
            need(h>=v[0]+v[1],('monotonicity',letter,v,w))
            if bank is SMALL:need(h<=w[2]<=8*h,('cone',letter,w))
            seen[w]=len(states);parents.append((idx,letter));states.append(w);todo.append(len(states)-1)
    return states,parents,seen
def word_for(idx,parents):
    word=[]
    while parents[idx][0]>=0:
        idx0,letter=parents[idx];word.append(letter);idx=idx0
    return tuple(reversed(word))
def minimum_orders(states,seen,bank=SMALL,costs=None):
    if costs is None:costs={l:len(d) for l,d in LETTERS.items()}
    dist=[10**99]*len(states);dist[0]=0;prev=[(-1,"")]*len(states);heap=[(0,0)]
    while heap:
        n,i=heapq.heappop(heap)
        if n!=dist[i]:continue
        for l,m in bank.items():
            w=apply3(states[i],m);j=seen.get(w)
            if j is None:continue
            if n+costs[l]<dist[j]:
                dist[j]=n+costs[l];prev[j]=(i,l);heapq.heappush(heap,(dist[j],j))
    return dist,prev
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--limit',type=int,default=100000);args=ap.parse_args()
    states,parents,seen=closure(args.limit)
    dist,parents=minimum_orders(states,seen)
    best={}
    for i,(a,b,c) in enumerate(states):
        h=a+b
        if h not in best or dist[i]<dist[best[h]]:best[h]=i
    witnesses={}
    for h in (1,3,5,7,11,21,31,35,37,43,49,55,56,57,63,77,91,111,119,133,343):
        i=best.get(h)
        if i is None:witnesses[str(h)]=None;continue
        w=word_for(i,parents);rows=targets(w)
        need(subset_h(rows)==h,('independent selected witness',h,w))
        witnesses[str(h)]={'word':w,'order':len(rows),'state':states[i],'rows':rows}
    pp=[p for p in primes(args.limit) if p!=7]
    seven=[p for p in primes(args.limit//7) if p!=3]
    missing_p=[p for p in pp if p not in best]
    missing_7p=[p for p in seven if 7*p not in best]
    certificates={}
    for cap in (37,133,343):
        ss,_,_ = closure(cap)
        certificates[str(cap)]=sorted(ss)
    digest=hashlib.sha256('\n'.join(','.join(map(str,v)) for v in sorted(states)).encode()).hexdigest()
    r={'status':'PASS','bound':args.limit,'reachable_states':len(states),
       'attained_H_values':len(best),'state_sha256':digest,
       'prime_lane':{'tested':len(pp),'attained':len(pp)-len(missing_p),'first_missing':missing_p[:30]},
       'seven_prime_lane':{'tested':len(seven),'attained':len(seven)-len(missing_7p),'first_missing_p':missing_7p[:30]},
       'selected_witnesses':witnesses,
       'H55_sources':[{'state':v,'word':word_for(i,parents),'after_s':sum(apply3(v,SMALL['s'])[:2])}
                      for i,v in enumerate(states) if sum(v[:2])==55],
       'closed_certificates':certificates,
       'max_minimum_order':max(dist),'order_histogram':{str(n):sum(x==n for x in dist) for n in sorted(set(dist))}}
    (HERE/'reachability.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:r[k] for k in ('bound','reachable_states','attained_H_values','prime_lane','seven_prime_lane','H55_sources','max_minimum_order')}))
if __name__=='__main__':main()
