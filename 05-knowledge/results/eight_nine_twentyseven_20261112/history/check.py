#!/usr/bin/env python3
"""Small exact layer and retained-clock controls; no external claims audit."""
import argparse
import hashlib
import itertools
import json
import math
import random
from pathlib import Path


def need(p,message):
    if not p:
        raise RuntimeError(message)


def layers(rows,c):
    n=len(rows)
    full=(1<<n)-1
    dp=[{} for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][(v,int(v<c))]=1
    for mask in range(1,1<<n):
        for (last,k),count in dp[mask].items():
            allowed=rows[last] & (full^mask)
            while allowed:
                bit=allowed & -allowed
                v=bit.bit_length()-1
                kk=k+int(v<c and last>=c)
                key=(v,kk)
                nxt=dp[mask|bit]
                nxt[key]=nxt.get(key,0)+count
                allowed ^= bit
    out=[0]*(c+1)
    for (_,k),count in dp[full].items():
        out[k]+=count
    return out


def random_context(c,e,r,rng):
    n=c+e
    rows=[0]*n
    def arc(x,y):
        rows[x] |= 1<<y
    for i,j in itertools.combinations(range(c),2):
        arc(i,j) if rng.randrange(2) else arc(j,i)
    for i,j in itertools.combinations(range(c,n),2):
        arc(i,j) if rng.randrange(2) else arc(j,i)
    for j in range(c,n):
        mask=rng.randrange(1,(1<<c)-1) if j<c+r else rng.choice([0,(1<<c)-1])
        for i in range(c):
            arc(j,i) if mask>>i & 1 else arc(i,j)
    return rows


def top_probe(k,s):
    n=3*k-3
    rows=[0]*n
    for i,j in itertools.combinations(range(n),2):
        rows[i] |= 1<<j
    # Exterior labels: x0,y0,x1,y1,...,z. Set C->x, y->C.
    for j in range(k,n-1):
        if (j-k)%2:
            for i in range(k):
                rows[i] &= ~(1<<j)
                rows[j] |= 1<<i
    z=n-1
    for i in range(s):
        rows[i] &= ~(1<<z)
        rows[z] |= 1<<i
    return rows


def run():
    rng=random.Random(89127)
    controls=[]
    for c,e in [(2,2),(3,3),(4,4),(5,5)]:
        for r in range(min(2,e)+1):
            for repeat in range(6):
                rows=random_context(c,e,r,rng)
                h=layers(rows,c)
                for k,count in enumerate(h[1:],1):
                    div=math.factorial(max(0,k-r))*math.factorial(max(0,k-r-1))
                    need(count%div==0,'shifted factorial divisor')
                if len(rows)<=8:
                    literal=[0]*(c+1)
                    for p in itertools.permutations(range(len(rows))):
                        if all(rows[x]>>y & 1 for x,y in zip(p,p[1:])):
                            k=sum(v<c and (i==0 or p[i-1]>=c) for i,v in enumerate(p))
                            literal[k]+=1
                    need(literal==h,'independent literal path layers')
                controls.append({'core_order':c,'exterior_order':e,'selective_vertices':r,
                                 'rows':rows,'layers':h})
    top=[]
    for k in range(2,6):
        for s in range(1,k):
            rows=top_probe(k,s)
            h=layers(rows,k)
            expected=(k-1)*s*(k-s)*math.factorial(k-2)**2
            need(h[k]==expected,'exact marked top-layer family')
            top.append({'k':k,'s':s,'top_layer':h[k],'whole_count':sum(h)})

    path=[8,1,15,10,6,3,13,12,4,5,11,14,2,7,9]
    need(set(path)==set(range(1,16)),'square path host')
    clock=[x+y for x,y in zip(path,path[1:])]
    need(set(clock)=={9,16,25},'three actual square clocks')
    def degree(v):
        return sum(1<=t-v<=15 and t-v!=v for t in [9,16,25])
    endpoints=[v for v in range(1,16) if degree(v)==1]
    need(endpoints==[8,9],'deleted fixed point and unmatched boundary')
    starts=[]
    for t in range(1,16):
        p=[t]
        for target in clock:
            p.append(target-p[-1])
        if len(set(p))==15 and set(p)==set(range(1,16)):
            starts.append(t)
    need(starts==[8],'indexed square-clock reconstruction')

    menus={1:[2,4],2:[3,7],4:[2,6],5:[3,5],7:[4,6],8:[5,7]}
    carry=[]
    for n in [7,25,43]:
        child=(16*n-1)//3
        need(n%9==7 and (16*n-1)%3==0,'common first clock branch')
        carry.append({'source':n,'source_mod27':n%27,'first_menu':menus[n%9],
                      'exponent':4,'child':child,'child_mod9':child%9,
                      'next_menu':menus[child%9]})
    need(len({tuple(x['next_menu']) for x in carry})==3,'three distinct retained ternary lifts')
    return {'status':'FINITE-EXACT controls for proved shifted factorial and inherited clock maps',
            'context_cases':controls,'top_layer_family':top,
            'square_path':path,'square_clock':clock,'square_endpoints':endpoints,
            'compatible_indexed_starts':starts,'collatz_carry_control':carry}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    data=json.dumps(run(),indent=2,sort_keys=True)+'\n'
    p=Path(__file__).with_name('receipt.json')
    if args.write:
        p.write_text(data)
    else:
        need(p.read_text()==data,'receipt drift')
    print(json.dumps({'passed':True,'receipt_sha256':hashlib.sha256(data.encode()).hexdigest()}))
