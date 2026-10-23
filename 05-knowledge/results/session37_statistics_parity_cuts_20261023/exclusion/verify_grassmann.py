#!/usr/bin/env python3
"""Finite tests for the explicitly defined exterior embeddings and hostile."""
from itertools import combinations, permutations
from pathlib import Path
import json

def wedge(a,b):
    if a&b:return 0,0
    inv=sum((b&((1<<i)-1)).bit_count() for i in range(a.bit_length()) if a>>i&1)
    return a|b,(-1)**inv

def images(n,even=True):
    return [3<<(2*i) for i in range(n)] if even else [1]+[3<<(2*i-1) for i in range(1,n)]

def basis_images(n,even):
    out={}
    gen=images(n,even)
    for s in range(1<<n):
        z=0;sgn=1
        for i in range(n):
            if s>>i&1:
                z,t=wedge(z,gen[i]);sgn*=t
        out[s]=(z,sgn)
    assert len(set(out.values()))==1<<n
    for a in range(1<<n):
        for b in range(1<<n):
            z,sign=wedge(out[a][0],out[b][0]);sign*=out[a][1]*out[b][1]
            if a&b:assert sign==0
            else:assert (z,sign)==out[a|b]
    return len(out)

def signed_count(n,mask):
    a=[[0]*n for _ in range(n)]
    for b,(i,j) in enumerate(combinations(range(n),2)):
        a[i][j]=mask>>b&1;a[j][i]=1-a[i][j]
    h=hs=0
    for p in permutations(range(n)):
        if all(a[p[i]][p[i+1]] for i in range(n-1)):
            h+=1
            hs+=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
    return h,hs

if __name__=='__main__':
    result={'faithful_embeddings':[]}
    for n in range(1,8):
        for even in [True,False]:
            size=basis_images(n,even)
            result['faithful_embeddings'].append({'n':n,'even':even,'odd_generators':2*n if even else 2*n-1,'basis_size':size})
    for n in range(1,4):
        for mask in range(1<<(n*(n-1)//2)):
            h,hs=signed_count(n,mask);assert h==abs(hs)
    h,hs=signed_count(4,4);assert (h,hs)==(5,-1)
    result['minimal_absolute_signed_count_hostile']={'n':4,'mask':4,'H':h,'H_alt':hs}
    pack=json.loads(Path(__file__).with_name('witness_verification.json').read_text())
    pol=[]
    for r in pack:
        terms={}
        for key,c in r['packings'].items():
            sizes=tuple(map(int,key.split(','))) if key else ()
            k=(sum(sizes),sum(j%2 for j in sizes),sum(j%2==0 for j in sizes))
            terms[k]=terms.get(k,0)+c
        pol.append(terms)
    diff={k:pol[0].get(k,0)-pol[1].get(k,0) for k in pol[0]|pol[1]}
    diff={k:v for k,v in diff.items() if v}
    assert diff=={(9,3,0):1,(9,1,0):-1,(8,0,2):1,(8,2,0):-1}
    result['realized_kernel_terms']=[{'occupied':k[0],'odd_cycles':k[1],'even_cycles':k[2],'coefficient':v} for k,v in sorted(diff.items())]
    responses=[]
    for terms in pol:
        response=[sum(c*(2**u) for (z,u,v),c in terms.items() if v==j) for j in range(3)]
        responses.append(response)
    assert responses==[[2843,1567,96],[2841,1567,95]]
    result['even_response_coefficients']=responses
    Path(__file__).with_name('grassmann_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
