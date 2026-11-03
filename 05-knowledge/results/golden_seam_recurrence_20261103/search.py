"""Complete bounded-H reachability; cost is original vertex count."""
from heapq import heappop, heappush
from pathlib import Path
from time import perf_counter
import argparse
import hashlib
import json
from transfer import POSITIVE, COST, fast_vector, need

HERE=Path(__file__).resolve().parent


def transition(v,m):
    h,x,y=v
    return (h*m[0][0]+x*m[1][0]+y*m[2][0],
            h*m[0][1]+x*m[1][1]+y*m[2][1],
            h*m[0][2]+x*m[1][2]+y*m[2][2])


def explore(bound, matrices=POSITIVE, costs=COST, cone=4):
    root=(1,0,0)
    best={root:0}
    parent={root:None}
    heap=[(0,root)]
    witnesses={}
    transitions=0
    ordered=tuple(sorted(matrices.items()))
    while heap:
        cost,v=heappop(heap)
        if best[v]!=cost: continue
        witnesses.setdefault(v[0],v)
        for letter,m in ordered:
            transitions+=1
            next_h=v[0]*m[0][0]+v[1]*m[1][0]+v[2]*m[2][0]
            if next_h>bound: continue
            w=transition(v,m)
            need(w[0]>=v[0] and 0<=w[1]<=cone*w[0] and 0<=w[2]<=cone*w[0], 'monotone bounded cone')
            newcost=cost+costs[letter]
            if w not in best or newcost<best[w]:
                best[w]=newcost
                parent[w]=(v,letter)
                heappush(heap,(newcost,w))
    def recover(v):
        word=[]
        while parent[v] is not None:
            v,letter=parent[v]
            word.append(letter)
        return tuple(word)
    proof=[]
    for h,v in sorted(witnesses.items()):
        word=recover(v)
        reconstructed=root
        for letter in reversed(word): reconstructed=transition(reconstructed,matrices[letter])
        need(reconstructed==v,('literal witness',h))
        proof.append({'H':h,'vertices':best[v],'word':word,'vector':v})
    digest=hashlib.sha256()
    for h,x,y in sorted(best): digest.update(f'{h},{x},{y}\n'.encode('ascii'))
    return {'bound':bound,'states':len(best),'transitions':transitions,
            'state_sha256':digest.hexdigest(),
            'largest_minimum_source_order':max(best.values()),
            'values':proof,'missing_odds':[h for h in range(1,bound+1,2) if h not in witnesses]}


def primes(bound):
    sieve=bytearray(b'\x01')*(bound+1)
    sieve[:2]=b'\x00\x00'
    for p in range(2,int(bound**.5)+1):
        if sieve[p]:
            sieve[p*p::p]=b'\x00'*(((bound-p*p)//p)+1)
    return [p for p in range(3,bound+1,2) if sieve[p]]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--bound',type=int,default=1000)
    parser.add_argument('--output',default='coverage.json')
    args=parser.parse_args()
    started=perf_counter()
    result=explore(args.bound)
    present={r['H'] for r in result['values']}
    p=primes(args.bound)
    result['prime_lane']={'eligible':len([a for a in p if a!=7]),
                          'missing':[a for a in p if a!=7 and a not in present]}
    p7=primes(args.bound//7)
    result['seven_prime_lane']={'eligible':len([a for a in p7 if a!=3]),
                                'missing_p':[a for a in p7 if a!=3 and 7*a not in present]}
    (HERE/args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('values','missing_odds','prime_lane','seven_prime_lane')}|{
        'prime_lane':{'eligible':result['prime_lane']['eligible'],'missing':len(result['prime_lane']['missing']),
                      'first_missing':result['prime_lane']['missing'][:10]},
        'seven_prime_lane':{'eligible':result['seven_prime_lane']['eligible'],'missing':len(result['seven_prime_lane']['missing_p']),
                            'first_missing_p':result['seven_prime_lane']['missing_p'][:10]},
        'attained_count':len(present),'missing_odd_count':len(result['missing_odds']),
        'seconds':round(perf_counter()-started,3)}))


if __name__=='__main__': main()
