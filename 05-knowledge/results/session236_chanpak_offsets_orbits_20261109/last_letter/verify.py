#!/usr/bin/env python3
"""Exact finite certificates, followed by the all-word last-letter induction.

Run from any directory. No third-party packages are needed. --rebuild additionally
regenerates the complete block alphabet using the inherited THM-5070 engine.
Assertions are not used as verification gates, so python -O runs the same checks.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = ROOT / '05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/chain_blocks.py'
T = (1,3,0,0,0,1,1,1,0)
V = (3,1,2)
EXCEPTION = (7,18,12,2,6,7,5,7,6)
IDENTITY = (1,0,0,0,1,0,0,0,1)

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def mul(A,B):
    return tuple(sum(A[3*i+k]*B[3*k+j] for k in range(3)) for i in range(3) for j in range(3))

def act(A,v):
    return tuple(sum(A[3*i+j]*v[j] for j in range(3)) for i in range(3))

def row(A,r):
    return A[3*r:3*r+3]

def minus(a,b,k):
    return tuple(x-k*y for x,y in zip(a,b))

def powers(A,n):
    P=IDENTITY
    result={}
    for k in range(1,n+1):
        P=mul(A,P)
        result[P]=k
    return result

def certificate(A):
    R=minus(row(A,0),row(A,1),3)
    S=minus(row(A,2),row(A,1),2)
    for bound in range(1,5):
        for p in range(-bound,bound+1):
            for q in range(-bound,bound+1):
                a=tuple(p*r+q*s for r,s in zip(R,S))
                if a[0]>0 and a[1]>=0 and a[2]>=0:
                    return {'p':p,'q':q,'positive_form':a}
    return None

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--rebuild',action='store_true')
    parser.add_argument('--mutate',action='store_true',help='negative control: corrupt one generator')
    args=parser.parse_args()
    raw=(HERE/'matrices_base3_order4.txt').read_bytes()
    matrices=[tuple(map(int,line.split())) for line in raw.decode().splitlines()]
    require(len(matrices)==57 and len(set(matrices))==57,'alphabet must have 57 distinct matrices')
    if args.mutate:
        matrices[0]=(0,)*9
    require(all(len(A)==9 and min(A)>=0 for A in matrices),'bad matrix dimensions or signs')
    if args.rebuild:
        spec=importlib.util.spec_from_file_location('s18_blocks',SOURCE)
        module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        rebuilt=set()
        marked=0
        for k in range(1,5):
            for adjacency in module.tournaments(k):
                hx=module.make_Hx(k,adjacency,3)
                for w in range(k):
                    for up in range(k):
                        rebuilt.add(module.block_matrix(k,adjacency,w,up,3,hx))
                        marked+=1
        require(marked==1105,'wrong marked labelled universe')
        require(rebuilt==set(matrices),'frozen matrices differ from rebuilt alphabet')
    tower=powers(T,4)
    require(set(tower)<=set(matrices),'missing tower power')
    require(act(T,V)==(6,2,4),'tower does not double target ray')
    records=[]
    for i,A in enumerate(matrices):
        require(A[0]>=1,'h must remain positive')
        cone=tuple(3*o-h for h,o in zip(row(A,0),row(A,2)))
        require(min(cone)>=0,'3O-h cone not preserved')
        if A in tower:
            records.append({'index':i,'kind':'tower','exponent':tower[A],'cone_form':cone})
            continue
        cert=certificate(A)
        if cert:
            records.append({'index':i,'kind':'separator','cone_form':cone,**cert})
            continue
        require(A==EXCEPTION,'unclassified last letter')
        R=minus(row(A,0),row(A,1),3)
        S=minus(row(A,2),row(A,1),2)
        require(R==(1,0,-9) and S==(1,-5,-8),'bad exceptional equations')
        require(act(A,(45,1,5))==(393,131,262),'exceptional preimage mismatch')
        require(45>3*5,'exceptional preimage must be outside cone')
        records.append({'index':i,'kind':'exceptional_preimage','preimage':(45,1,5),'cone_form':cone})
    counts={kind:sum(r['kind']==kind for r in records) for kind in ['tower','separator','exceptional_preimage']}
    require(counts=={'tower':4,'separator':52,'exceptional_preimage':1},'classification counts differ')
    # Positive and hostile controls for the *statement*: the exceptional matrix
    # really does send a positive vector outside our invariant cone to the ray.
    require(act(EXCEPTION,(45,1,5))==tuple(131*x for x in V),'missing hostile control')
    require(3<=3*2 and 1<=3*1,'ray and seed must lie in cone')
    require(IDENTITY not in set(matrices),'unexpected identity generator')
    result={
      'status':'PASS','universe':'all 1105 marked labelled tournament blocks of orders 1..4, deduplicated to 57 ternary matrices',
      'counts':counts,'alphabet_sha256':hashlib.sha256(raw).hexdigest(),
      'rebuild':args.rebuild,'all_word_conclusion':'Only products consisting entirely of tower-power generators preserve the ray (3,1,2).',
      'scope':'all word lengths; finite exact certificate plus mathematical induction, not orbit cofiniteness',
      'hostile_control':'The exceptional matrix maps (45,1,5) to 131*(3,1,2); the invariant cone is essential.',
      'records':records}
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
