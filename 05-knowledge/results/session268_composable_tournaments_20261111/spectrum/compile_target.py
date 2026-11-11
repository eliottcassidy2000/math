#!/usr/bin/env python3
"""Minimum vertex witness inside an explicit infinite compositional language.

Usage: python3 compile_target.py INTEGER [--output witness.json]
No result outside this language is implied by a negative answer.
"""
from functools import lru_cache
from pathlib import Path
import argparse,json
from replay import transitive,cyclic,atom,mask_graph,seam,cover_vector,count

def ordinal(parts):
    rows=[];total=sum(map(len,parts));offset=0
    for P in parts:
        after=total-offset-len(P);cross=((1<<after)-1)<<(offset+len(P))
        rows.extend((r<<offset)|cross for r in P);offset+=len(P)
    return tuple(rows)

def seeds():
    C3=cyclic((0,),(0,),(0,))
    return (ordinal([(0,),cyclic((0,),(0,),transitive(3))]),ordinal([C3,C3]))

def bank(bound,modes=('atom','seam','lift')):
    best={}
    def put(h,size,kind,pars):
        if h<=bound and (h not in best or size<best[h]['size']):
            best[h]={'h':h,'size':size,'kind':kind,'parameters':pars}
    bits=bound.bit_length()
    for b in range(1,bits+1):
        for c in range(1,bits+1):
            h=atom(b,c)
            if 'atom' in modes and b<=c:put(h,b+c+1,'atom',[b,c])
            if 'seam' in modes:
                for i,a,r in [(0,24,15),(1,20,11)]:
                    put(a*h+r*((1<<b)-2)-1,b+c+6,'seam',[i,b,c])
    if 'lift' in modes:
        for i,(h,Q) in enumerate([(9,150),(9,126)]):
            for b in range(1,bits+1):put((h+Q)*(1<<b)+h-Q,b+7,'lift',[i,b])
    return best

def realize(spec):
    p=spec['parameters'];kind=spec['kind']
    if kind=='atom':return cyclic((0,),transitive(p[0]),transitive(p[1]))
    if kind=='seam':
        i,b,c=p
        return seam(cyclic((0,),transitive(b),transitive(c)),mask_graph(5,[8,10][i]))
    if kind=='lift':
        i,b=p
        return cyclic((0,),transitive(b),seeds()[i])
    raise ValueError(kind)

def compile_target(N,modes=('atom','seam','lift'),materialize=True):
    if N<1 or N%2==0:raise ValueError('target must be positive odd')
    if N==1:return {'target':1,'vertices':1,'factors':[],'rows':[0] if materialize else None,'states':1}
    generators=bank(N,modes);states=0
    @lru_cache(None)
    def solve(m):
        nonlocal states
        states+=1
        if m==1:return 0,()
        best=None
        for g,spec in generators.items():
            if m%g:continue
            rest=solve(m//g)
            if rest is None:continue
            cand=(rest[0]+spec['size'],rest[1]+(g,))
            if best is None or cand<best:best=cand
        return best
    result=solve(N)
    if result is None:return {'target':N,'attained_in_language':False,'states':states}
    cost,factors=result;specs=[generators[g] for g in factors]
    return {'target':N,'attained_in_language':True,'vertices':cost,'factors':specs,
            'rows':list(ordinal([realize(s) for s in specs])) if materialize else None,'states':states}

def measurements():
    expected=[(9,150,756,1704,1800,720),(9,126,684,1656,1800,720)]
    for T,q in zip(seeds(),expected):
        if cover_vector(T)!=q:raise AssertionError('literal hidden-state seed')
    targets=[59,71,153,177,393,495,59*47,153*177,1000001]
    data=[]
    for N in targets:
        row={'target':N}
        for name,modes in [('atom',('atom',)),('atom_seam',('atom','seam')),('full',('atom','seam','lift'))]:
            ans=compile_target(N,modes,materialize=False)
            row[name+'_vertices']=ans.get('vertices');row[name+'_states']=ans['states']
        full=compile_target(N)
        if full.get('vertices',1000)<=16 and count(tuple(full['rows']))!=N:
            raise AssertionError('compiled physical count')
        data.append(row)
    return {'scope':'minimum within declared generator bank plus ordered SCC product language, not all tournaments',
            'targets':data,'generator_bank_to_1e6':{name:len(bank(10**6,modes)) for name,modes in
             [('atom',('atom',)),('atom_seam',('atom','seam')),('full',('atom','seam','lift'))]}}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('target',type=int,nargs='?');p.add_argument('--output');p.add_argument('--measure',action='store_true');a=p.parse_args()
    result=measurements() if a.measure else compile_target(a.target)
    data=json.dumps(result,indent=2)+'\n'
    if a.output:Path(a.output).write_text(data)
    print(data,end='')
