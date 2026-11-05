#!/usr/bin/env python3
"""Independent controls for sharp H=4151 equality, beyond the source verifier."""
import itertools
import json
import hashlib
from pathlib import Path
from verify_prime_rigidity import adjacency,automorphisms
from verify_symmetry_threshold import strong,inflate


def all_subset_H(out):
    n=len(out); full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        unused=full^mask
        for i,count in enumerate(dp[mask]):
            if count:
                todo=out[i]&unused
                while todo:
                    b=todo&-todo; todo-=b
                    dp[mask|b][b.bit_length()-1]+=count
    return [sum(row) if m else 1 for m,row in enumerate(dp)]


def strong_mask(out,mask):
    for v in range(len(out)):
        if not mask>>v&1: continue
        seen=1<<v
        while True:
            new=seen
            for i in range(len(out)):
                if seen>>i&1: new|=out[i]&mask
            if new==seen: break
            seen=new
        if seen!=mask: return False
    return True


def isomorphism(a,b):
    ca={};cb={}
    for i,row in enumerate(a): ca.setdefault(row.bit_count(),[]).append(i)
    for i,row in enumerate(b): cb.setdefault(row.bit_count(),[]).append(i)
    if {k:len(v) for k,v in ca.items()}!={k:len(v) for k,v in cb.items()}:return None
    keys=sorted(ca)
    for options in itertools.product(*(itertools.permutations(cb[k]) for k in keys)):
        p=[0]*len(a)
        for k,opt in zip(keys,options):
            for i,j in zip(ca[k],opt):p[i]=j
        if all(bool(a[i]>>j&1)==bool(b[p[i]]>>p[j]&1)
               for i in range(len(a)) for j in range(i+1,len(a))):return p
    return None


def main():
    total=0;strict_cases=0;rooted_cases=0;strong_counts={}
    for n in range(1,7):
        strong_counts[n]=0
        for code in range(1<<(n*(n-1)//2)):
            out=adjacency(n,code); total+=1
            if not strong(out):continue
            strong_counts[n]+=1
            hs=all_subset_H(out);full=(1<<n)-1;covered=0
            for mask in range(1,full+1):
                if mask.bit_count()<3 or not strong_mask(out,mask):continue
                if mask<full:
                    if not hs[mask]<hs[full]:raise RuntimeError(('strictness',n,code,mask))
                    strict_cases+=1
                if mask.bit_count()==4:covered|=mask
            if n>=4:
                if covered!=full:raise RuntimeError(('rooted strong4',n,code))
                rooted_cases+=n
    q=[2,8,3,5]
    c=[sum(1<<((i+j)%7) for j in (1,2,3)) for i in range(7)]
    cores=[inflate(q,r,c) for r in (0,2)]
    rev=lambda a:[(((1<<len(a))-1)^(1<<i))^row for i,row in enumerate(a)]
    if isomorphism(cores[0],cores[1]) is not None:raise RuntimeError('cores isomorphic')
    if isomorphism(cores[0],rev(cores[0])) is not None:raise RuntimeError('first self-converse')
    if isomorphism(cores[1],rev(cores[1])) is not None:raise RuntimeError('second self-converse')
    converse_map=isomorphism(cores[0],rev(cores[1]))
    if converse_map is None:raise RuntimeError('cores not converse')
    core_data=[]
    for core in cores:
        hs=all_subset_H(core);aut=automorphisms(core)
        if hs[-1]!=4151 or aut!=7 or not strong(core):raise RuntimeError('core invariant')
        core_data.append({'rows':core,'H':hs[-1],'Aut':aut,'strong':True})
    result={'status':'PASS','labelled_universe_n1_through_n6':total,
            'strong_labelled_counts':strong_counts,'strict_strong_inclusion_cases':strict_cases,
            'rooted_strong4_vertex_cases':rooted_cases,'cores':core_data,
            'cores_distinct':True,'neither_self_converse':True,'first_to_reverse_second_isomorphism':converse_map}
    path=Path(__file__).with_name('equality_independent_receipt.json')
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,sort_keys=True))

if __name__=='__main__':main()
