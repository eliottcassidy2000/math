"""A000568 from actual fresh-path reset and independent odd-cycle Burnside."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
import json
from math import factorial, gcd
from pathlib import Path
from fields import need


def odd_partitions(n,minimum=1):
    if n==0:
        yield ()
    for i in range(minimum,n+1,2):
        for tail in odd_partitions(n-i,i): yield (i,)+tail


def burnside(n):
    result=Fraction(0); partitions=0
    for lengths in odd_partitions(n):
        ex=sum((l-1)//2 for l in lengths)+sum(gcd(a,b) for a,b in combinations(lengths,2))
        z=1
        for l,k in Counter(lengths).items(): z *= (l**k)*factorial(k)
        result += Fraction(1<<ex,z); partitions+=1
    need(result.denominator==1,'Burnside integrality')
    return int(result),partitions


def reset_census(n):
    free=[(i,j) for i in range(n) for j in range(i+2,n)]
    index={e:k for k,e in enumerate(free)}; M=len(free); N=1<<M
    maps=[]
    for p in permutations(range(n)):
        ones=zeros=0; impossible=False
        for u,v in zip(p,p[1:]):
            low,high=sorted((u,v))
            if high==low+1:
                if u>v: impossible=True; break
            elif u<v: ones |= 1<<index[low,high]
            else: zeros |= 1<<index[low,high]
        if impossible: continue
        base=0; channels=[]
        for k,(i,j) in enumerate(free):
            u,v=p[i],p[j]; low,high=sorted((u,v))
            if high==low+1:
                if u<v: base |= 1<<k
            else:
                if u>v: base |= 1<<k
                channels.append((index[low,high],k))
        maps.append((ones,zeros,base,channels))
    signatures=[]; diagonal=[]; class_id={}; trace=Fraction(0); hp_total=0
    for x in range(N):
        row=Counter()
        for ones,zeros,base,channels in maps:
            if x&ones!=ones or x&zeros: continue
            y=base
            for old,new in channels:
                if x>>old&1: y ^= 1<<new
            row[y]+=1
        need(x in row,'identity path')
        aut=row[x]; H=sum(row.values()); keys=tuple(sorted(row))
        need(set(row.values())=={aut},'uniform orbit fibres')
        need(H==aut*len(keys),'correct presentation mass')
        signatures.append(keys); diagonal.append(Fraction(aut,H)); trace+=diagonal[-1]; hp_total+=H
        class_id.setdefault(keys,len(class_id))
    for x,keys in enumerate(signatures):
        need(all(signatures[y]==keys for y in keys),'identical reset rows / projector')
    need(trace==len(class_id),'reset trace counts classes')
    need(trace==burnside(n)[0],'reset versus independent Burnside')
    # R*tau is the actual complement after reset. Its trace need not count classes.
    tau_trace=sum((diagonal[x] if x^(N-1) in signatures[x] else Fraction(0)) for x in range(N))
    return dict(n=n,presentations=N,ordinary_classes=int(trace),reset_rank=int(trace),
                reset_trace=str(trace),reset_complement_trace=str(tau_trace),
                path_instances=hp_total,class_size_histogram=dict(sorted(Counter(map(len,class_id)).items())))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output');parser.add_argument('--max-order',type=int,default=50);args=parser.parse_args()
    sequence=[]; partition_counts=[]
    for n in range(args.max_order+1):
        value,count=burnside(n); sequence.append(value); partition_counts.append(count)
    strong=[0]*(args.max_order+1)
    for n in range(1,args.max_order+1): strong[n]=sequence[n]-sum(strong[k]*sequence[n-k] for k in range(1,n))
    need(sequence[:11]==[1,1,1,2,4,12,56,456,6880,191536,9733056],'published initial values')
    need(strong[:9]==[0,1,0,1,1,6,35,353,6008],'strong composition head')
    result=dict(ordinary_counts=sequence,strong_counts=strong,odd_partition_counts=partition_counts,
                reset_censuses=[reset_census(n) for n in range(1,7)],empty_class_convention=1)
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output: Path(args.output).write_text(text,encoding='utf-8',newline='\n')
    else: print(text,end='')


if __name__=='__main__':main()
