#!/usr/bin/env python3
"""Independent odd-partition comparisons for the production reverse engine."""
from importlib.util import module_from_spec,spec_from_file_location
from math import factorial,gcd
from pathlib import Path
import hashlib,json


def load(name,path):
    spec=spec_from_file_location(name,path);module=module_from_spec(spec);spec.loader.exec_module(module);return module
def v2(n):return (n&-n).bit_length()-1
def partitions(n,lo=1):
    if n==0:yield ();return
    for d in range(lo,n+1,2):
        for tail in partitions(n-d,d):yield (d,)+tail
def partition_terms(n):
    out=[]
    for p in partitions(n):
        E=sum((d-1)//2 for d in p)+sum(gcd(p[i],p[j]) for i in range(len(p)) for j in range(i+1,len(p)))
        den=1
        for d in set(p):den*=d**p.count(d)*factorial(p.count(d))
        valuation=E-v2(den)
        out.append((valuation,den>>v2(den)))
    return out
def main():
    base=Path(__file__).resolve().parent
    reverse=load('reverse_production',base/'reverse_ordinary_check.py')
    modular=load('absolute_production',base.parent/'port_residual_integral_20261025/counting/modular.py')
    comparisons=edges=0;samples=[];zero_seen=False;cuts_seen=False
    for n in range(61):
        graph=reverse.build_graph(n);terms=partition_terms(n)
        B=min(v for v,d in terms)
        if B!=graph['baseline'] or B!=(n//2 if n>=3 else 0):raise RuntimeError('independent baseline')
        edges+=sum(map(len,graph['edges'].values()))
        for K in range(1,9):
            mod=1<<K
            expected=sum((pow(2,v-B,mod)*pow(d,-1,mod)) for v,d in terms)%mod
            result,stats=reverse.normalized_target(n,K,graph)
            if result!=expected:raise RuntimeError(('partition disagreement',n,K,result,expected))
            comparisons+=1;zero_seen|=stats['retained_zero_residue_nodes']>0;cuts_seen|=stats['active_nodes']<stats['structural_nodes']
            if n in (14,30,50,60) and K in (1,2,8):samples.append(stats|{'residue':result})
        if n<=40 or n==50:
            absolute,_=modular.modular_target(n,B+3,False)
            result,_=reverse.normalized_target(n,3,graph)
            if absolute%(1<<B) or absolute>>B!=result:raise RuntimeError('absolute engine disagreement')
    hostile,_=reverse.normalized_target(50,2)
    if hostile!=2 or not zero_seen or not cuts_seen:raise RuntimeError('precision/support hostile not exercised')
    for n,K in [(-1,2),(2,0),(True,2),(3,1.5)]:
        try:reverse.normalized_target(n,K)
        except ValueError:pass
        else:raise RuntimeError('invalid input accepted')
    print(json.dumps({'status':'FINITE-EXACT','orders':[0,60],'normalized_bits':[1,8],
                     'independent_partition_comparisons':comparisons,'structural_edges_checked':edges,
                     'retained_zero_residue_nodes_exercised':zero_seen,'structural_precision_pruning_exercised':cuts_seen,
                     'a50_normalized_mod4':hostile,'samples':samples,
                     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))
if __name__=='__main__':main()
