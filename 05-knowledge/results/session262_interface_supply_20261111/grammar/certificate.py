#!/usr/bin/env python3
"""Cheap new controls only: nonlinear reuse and finite inverse-cover certificates."""
from pathlib import Path
from math import lcm
import importlib.util,json,random,sys
sys.dont_write_bytecode=True
import sympy as s
ROOT=Path(__file__).resolve().parents[4]
BASE=Path(__file__).resolve().parent
p=ROOT/'05-knowledge/results/session261_creative_chanpak_20261111/wildcard/check_automaton.py'
spec=importlib.util.spec_from_file_location('spine',p); spine=importlib.util.module_from_spec(spec);spec.loader.exec_module(spine)

def inverse_cover(B,M,residues,edges):
    # Edges are (source,destination,a,b), with a >= 2.
    C=max([B+1]+[a*B+b for _,_,a,b in edges]+[(-b)//(a-1)+1 for _,_,a,b in edges])
    L=lcm(M,*[a*M for _,_,a,b in edges])
    missing=[]
    for q,R in residues.items():
        for n in range(L):
            if n%M not in R: continue
            if not any(dst==q and (n-b)%a==0 and ((n-b)//a)%M in residues[src] for src,dst,a,b in edges):missing.append([q,n])
    seeds={q:[n for n in range(B,C) if n%M in R] for q,R in residues.items()}
    return {'B':B,'M':M,'C':C,'L':L,'required_initial_values':seeds,'uncovered_residues':missing,'arithmetic_cover':not missing}

def mm(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def mv(A,v):return [sum(a*b for a,b in zip(row,v)) for row in A]
def ident(d):return [[int(i==j) for j in range(d)] for i in range(d)]

def main():
    x,h=s.symbols('x h');g=h-15*x+1
    for e in [0,1]:
        Fh=2**e*(h+1)**2/15-1
        residual=s.expand(Fh-15*2**e*x*x+1)
        factor=2**e*(h+1+15*x)/15
        if s.expand(residual-g*factor)!=0:raise RuntimeError('ideal membership failure')
    # Physical source and copy operations, not merely scalar values.
    arcs=[[False,True],[False,False]]
    iv,acc,tr=spine.compile_automaton(2,arcs);d=len(iv)
    mats=[]
    for transitions in tr:
        A=[[0]*d for _ in range(d)]
        for i,ds in enumerate(transitions):
            for j in ds:A[j][i]+=1
        mats.append(A)
    checked=[]
    for m in range(1,10):
        word=[1,2,0]+[2]*m
        H=spine.automaton_count(2,arcs,word)
        if H!=15*2**m-1:raise RuntimeError('ray failure')
        for e in [0,1]:
            target=word[:3]+word[3:]+word[3:]+([2] if e else [])
            actual=spine.automaton_count(2,arcs,target)
            expected=2**e*(H+1)**2//15-1
            if actual!=expected:raise RuntimeError('physical copy failure')
        checked.append([m,H])
    # New all-word square-context identity on 85 tiny words.
    from itertools import product
    square_checks=0
    for n in range(4):
        for w in product(range(4),repeat=n):
            X=ident(d)
            for b in w:X=mm(mats[b],X)
            v=mv(mm(X,X),iv)
            H=sum(a*b for a,b in zip(acc,v))
            if H!=spine.automaton_count(2,arcs,w+w):raise RuntimeError('context square failure')
            square_checks+=1
    # Independent direct path DP for six new copied ray instances.
    literal=0
    for m in [1,2,3]:
        for e in [0,1]:
            w=[1,2,0]+[2]*(2*m+e)
            if spine.hp(spine.tournament(2,arcs,w))!=15*2**(2*m+e)-1:raise RuntimeError('literal DP failure')
            literal+=1
    binary=inverse_cover(23,2,{'q':{1}},[('q','q',2,-1),('q','q',2,1)])
    one_digit=inverse_cover(23,2,{'q':{1}},[('q','q',2,1)])
    if not binary['arithmetic_cover'] or one_digit['arithmetic_cover']:raise RuntimeError('coverage hostile failed')
    if binary['required_initial_values']['q']!=list(range(23,46,2)):raise RuntimeError('initial case mismatch')
    # Compact code m in binary: initial m=1, each subsequent digit applies m->2m+e.
    for m in range(1,1001):
        state=1
        for e in bin(m)[3:]:state=2*state+int(e)
        if state!=m:raise RuntimeError('initial-case decoder')
    report={'status':'FINITE-EXACT new controls; all-size proofs in PROOF.md','ideal_identity':'g(F_e)=2^e*(h+1+15x)*g/15','physical_ray_values':checked,'square_context_cases':square_checks,'independent_literal_path_counts':literal,'binary_index_codes':1000,'two_digit_arithmetic_certificate':binary,'plus_only_hostile':one_digit,'phantom_zero_module':{'formal_x':1,'formal_h':14,'actual_empty_block_count':spine.automaton_count(2,arcs,[1,2,0])}}
    (BASE/'receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
