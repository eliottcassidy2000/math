#!/usr/bin/env python3
"""Independent targeted audit of bridge kernels and symmetric class inventory.

No producer mathematical function is imported. Kernels use a product of
forced/optional root-gap polynomials for each exterior permutation.
"""
from collections import Counter
from itertools import permutations
from math import comb,isqrt
from pathlib import Path
import json,hashlib
HERE=Path(__file__).resolve().parent;BRIDGE=HERE.parent/'bridge'

def need(ok,message):
    if not ok: raise RuntimeError(message)

def matrix_frame(n,word):
    a=[[False]*n for _ in range(n)];pos=0
    for i in range(n):
        for j in range(i+1,n):
            bit=True if j==i+1 else bool(word>>pos&1)
            if j>i+1:pos+=1
            a[i][j]=bit;a[j][i]=not bit
    return a

def is_strong(a):
    n=len(a)
    for reverse in (False,True):
        seen={0};todo=[0]
        while todo:
            v=todo.pop()
            for w in range(n):
                if w not in seen and (a[w][v] if reverse else a[v][w]):
                    seen.add(w);todo.append(w)
        if len(seen)!=n:return False
    return True

def gap_kernel(a):
    n=len(a);result=[0]*(n+1)
    if n==1:return [1,1]
    for p in permutations(range(1,n)):
        required=0;optional=int(a[0][p[0]])+int(a[p[-1]][0]);valid=True
        for u,v in zip(p,p[1:]):
            bridge=a[u][0] and a[0][v]
            if not a[u][v]:
                if not bridge:valid=False;break
                required+=1
            elif bridge:optional+=1
        if valid:
            for j in range(optional+1):
                result[required+j]+=comb(optional,j)
    return result

def rooted_canonical(a):
    n=len(a)
    return min(tuple(int(a[p[i]][p[j]]) for i in range(n) for j in range(i+1,n))
               for tail in permutations(range(1,n)) for p in [(0,)+tail])

def is_prime(x):
    return x>=2 and all(x%d for d in range(2,isqrt(x)+1))

def hamilton_predecessors(rows):
    n=len(rows);pred=[0]*n
    for v in range(n):
        pred[v]=sum(1<<u for u in range(n) if rows[u]>>v&1)
    dp=[[0]*n for _ in range(1<<n)]
    for s in range(1,1<<n):
        bits=s
        while bits:
            bit=bits&-bits;bits-=bit;v=bit.bit_length()-1;old=s^bit
            if old==0:dp[s][v]=1;continue
            options=old&pred[v];total=0
            while options:
                b=options&-options;options-=b;total+=dp[old][b.bit_length()-1]
            dp[s][v]=total
    return sum(dp[-1])

def main():
    data=json.loads((BRIDGE/'check.json').read_text())
    C=data['module_ordered_covers']
    need(C==[0,175,1638,6678,15120,20160,15120,5040],'module vector source')
    inventory={};tested=0
    for n in range(3,7):
        for word in range(1<<((n-1)*(n-2)//2)):
            a=matrix_frame(n,word)
            if not is_strong(a):continue
            k=gap_kernel(a);tested+=1
            value=sum(x*y for x,y in zip(k,C))
            if value<19467 and value%7==0 and is_prime(value//7) and value//7>=11:
                inventory.setdefault((value,n+6),set()).add(rooted_canonical(a))
    observed={f'{value}:{n}':len(classes) for (value,n),classes in inventory.items()}
    expected={'4151:10':2,'5789:10':2,'14791:12':2,'16079:12':4,'17717:12':2}
    need(observed==expected,'independent marked quotient inventory')
    basischecks=[]
    for item in data['frame_universes']:
        n=item['q']
        for row in item['pareto_kernels']:
            a=matrix_frame(n,row['frame']);got=gap_kernel(a)
            need(got[1:]==row['coefficients'],'independent Pareto source kernel')
            basischecks.append({'q':n,'frame':row['frame'],'coefficients':got[1:]})
    witnesschecks=[]
    for row in data['complete_symmetric_inventory']:
        rows=row['rows'];value=hamilton_predecessors(rows)
        need(value==row['H'],'independent actual target path count')
        p=[(v+1)%7 if v<7 else v for v in range(len(rows))]
        need(all(bool(rows[i]>>j&1)==bool(rows[p[i]]>>p[j]&1)
                 for i in range(len(rows)) for j in range(len(rows))),'actual C7 symmetry')
        witnesschecks.append([value,len(rows)])
    print(json.dumps({'status':'PASS','scope':'Targeted independent audit, not a second full q7 scan',
        'strong_root_frames_q3_through6':tested,'inventory':observed,
        'pareto_sources':basischecks,'actual_witnesses':witnesschecks,
        'producer_receipt_sha256':hashlib.sha256((BRIDGE/'check.json').read_bytes()).hexdigest()},indent=2,sort_keys=True))
if __name__=='__main__':main()
