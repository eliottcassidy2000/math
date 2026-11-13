#!/usr/bin/env python3
"""Exact small controls and literal interfaces for the wrapped neutral refill."""
from functools import lru_cache
from itertools import permutations
from math import comb, gcd
from pathlib import Path
import json

A = (2,28,25,81,97,79,7)
B = (2,80,3,7,173,143,189,15)
PI = (462,508,120,496,481,449,384,260,4)
AP = [0,14,0,60,0,20,0,1]

def require(ok, text):
    if not ok: raise RuntimeError(text)

def polyadd(*terms):
    out=[0]*max(map(len,terms))
    for term in terms:
        for i,a in enumerate(term): out[i]+=a
    while len(out)>1 and out[-1]==0: out.pop()
    return out

def polymul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return out

def shift(f,s):
    out=[0]*len(f)
    for i,a in enumerate(f):
        for j in range(i+1): out[j]+=a*comb(i,j)*s**(i-j)
    return out

def wrap(f):
    return polyadd([0]+shift(f,1),[0,0]+f,[0]+[-x for x in shift(f,-1)])

def value(f,x):
    out=0
    for a in reversed(f): out=out*x+a
    return out

def count(j,f):
    h,p,w=(value(f,k) for k in (1,2,3))
    return 120*3**j*w+(240*3**j+1)*p+(1-120*3**j)*h

def named(rows,prefix):
    labels=[f'{prefix}{i}' for i in range(len(rows))]
    return {x:{labels[j] for j in range(len(rows)) if rows[i]>>j&1}
            for i,x in enumerate(labels)}

def single(label): return {label:set()}

def substitute(rows,blocks):
    out={v:set(E) for block in blocks for v,E in block.items()}
    require(len(out)==sum(map(len,blocks)), 'overlapping labels')
    for i,block in enumerate(blocks):
        for j,other in enumerate(blocks):
            if rows[i]>>j&1:
                for v in block: out[v].update(other)
    return out

def ordinal(*blocks):
    n=len(blocks)
    return substitute(tuple(sum(1<<j for j in range(i+1,n)) for i in range(n)),blocks)

def cyclic(T,x='x',y='y'):
    return substitute((2,4,1),(T,single(x),single(y)))

def outer(M,T):
    return substitute(PI,(single('a'),cyclic(T),single('c'),M,single('d0'),
                          single('d1'),single('x0'),single('f0'),single('f1')))

def masks(G):
    labels=list(G)
    return tuple(sum(1<<j for j,w in enumerate(labels) if w in G[v]) for v in labels)

def subset_h(rows):
    n=len(rows); dp=[[0]*n for _ in range(1<<n)]; hs=[0]*(1<<n); hs[0]=1
    for v in range(n): dp[1<<v][v]=1
    full=(1<<n)-1
    for S in range(1,1<<n):
        hs[S]=sum(dp[S])
        for v,a in enumerate(dp[S]):
            if not a: continue
            rest=rows[v]&~S&full
            while rest:
                bit=rest&-rest; rest-=bit
                dp[S|bit][bit.bit_length()-1]+=a
    return hs

def slot(rows,k):
    h=subset_h(rows); full=len(h)-1; old=[0]*len(h); old[0]=1
    for _ in range(k):
        new=[0]*len(h)
        for S in range(len(h)):
            sub=S
            while True:
                new[S]+=old[S^sub]*h[sub]
                if not sub: break
                sub=(sub-1)&S
        old=new
    return old[full]

def run_signature(j):
    M=ordinal(*(single(f'l{i}') for i in range(j)))
    G=outer(M,single('t')); outside=[v for v in G if v!='t']; full=(1<<len(outside))-1
    values=[]
    for runs in (1,2,3):
        @lru_cache(None)
        def f(S,k,last):
            if S==full and k==runs: return 1
            z=0
            for i,v in enumerate(outside):
                if not S>>i&1 and (last=='' or v in G[last]):
                    z+=f(S|1<<i,k,v)
            if k<runs and last!='t' and (last=='' or 't' in G[last]):
                z+=f(S,k+1,'t')
            return z
        values.append(f(0,0,''))
    require({v for v in outside if 't' in G[v]}=={'a','y'},'run support guard')
    require(values==[720*3**j+3,600*3**j+1,120*3**j],'run signature')
    return values

def cycle(M,T,f,j,d,step):
    old=outer(M,T); U=T; uf=f
    for k in range(d):
        U=cyclic(U,f'i{step}_{k}x',f'i{step}_{k}y'); uf=wrap(uf)
    fresh=named(B,f'b{step}_'); w=f'b{step}_7'
    prepared=outer(M,ordinal(U,fresh))
    replacement=named(A,f'b{step}_')
    neutral=outer(M,ordinal(U,replacement,single(w)))
    moved={v:set(E) for v,E in neutral.items()}
    for z in ('c','x'):
        require(z in moved[w], 'wrong transfer direction')
        moved[w].remove(z); moved[z].add(w)
    newT=ordinal(U,replacement); newM=ordinal(single(w),M)
    target=outer(newM,newT)
    require(moved==target,'two-flip target mismatch')
    for G in (prepared,neutral,target):
        require(all((G[v]&old.keys())==E for v,E in old.items()),'old induced graph changed')
    nf=polymul(uf,AP)
    # Q_B=t Q_A, certified in inherited S310; protected outside readouts agree.
    prepared_f=polymul(uf,[0]+AP)
    require(count(j,prepared_f)==count(j,[0]+nf),'neutral count')
    h,p=(value(nf,k) for k in (1,2))
    drift=(240*3**j-1)*p-240*3**j*h
    require(count(j+1,nf)-count(j,[0]+nf)==drift,'transfer drift')
    return newM,newT,nf,j+1,{
        'd':d,'old_order':len(old),'new_order':len(target),
        'prepared_count':count(j,prepared_f),'neutral_count':count(j,[0]+nf),
        'final_count':count(j+1,nf),'two_flips':[[w,'c'],[w,'x']],
        'old_vertices_retained':len(old),'inner_order':len(newT),
        'mark': 'seed0', 'mark_retained': 'seed0' in newT}

def modular_cycle(values,d,mod):
    for _ in range(d):
        values=[(t*values[(t+1)%mod]+t*t*values[t]-t*values[(t-1)%mod])%mod
                for t in range(mod)]
    return [value(AP,t)*x%mod for t,x in enumerate(values)]

def select_target(target,j0=1):
    """Complete inverse for this fixed seed/rule language; no module oracle."""
    K0=240*3**j0*2**7
    answers=[]; examined=0
    def visit(word,f,j,weight):
        nonlocal examined
        if K0*weight>target: return
        examined+=1
        if count(j,f)==target: answers.append(list(word))
        child_weight=384*weight; wf=f; d=0
        while K0*child_weight<=target:
            visit(word+(d,),polymul(AP,wf),j+1,child_weight)
            wf=wrap(wf); d+=1; child_weight*=4
    visit((),AP,j0,1)
    return {'target':target,'addresses':answers,'bounded_states_examined':examined}

def main():
    controls=[]
    small=[(0,),(2,0),(2,4,1),(2,4,9,3)]
    for rows in small:
        T=named(rows,'t'); C=cyclic(T)
        for k in (1,2,3):
            rhs=k*slot(rows,k+1)+k*k*slot(rows,k)-k*slot(rows,k-1)
            require(slot(masks(C),k)==rhs,'cyclic profile transform')
        if len(rows)<=3:
            for j in (1,2):
                M=ordinal(*(single(f'l{i}') for i in range(j)))
                h,p,w=(slot(rows,k) for k in (1,2,3))
                expected=120*3**j*w+(240*3**j+1)*p+(1-120*3**j)*h
                actual=subset_h(masks(outer(M,T)))[-1]
                require(actual==expected,'outer count')
                controls.append({'inner_rows':rows,'j':j,'H':actual})
    # Fixed rigid prime quotient of the outer frame.
    quotient=(110,124,56,112,97,64,4)
    modules=[]
    for S in range(1,127):
        if S.bit_count()<2: continue
        if all((quotient[v]&S) in (0,S) for v in range(7) if not S>>v&1): modules.append(S)
    aut=[p for p in permutations(range(7)) if all(
        ((quotient[i]>>j)&1)==((quotient[p[i]]>>p[j])&1)
        for i in range(7) for j in range(7))]
    require(not modules and len(aut)==1,'prime rigid quotient')
    paths=[]
    for word in ((0,0),(0,1,0),(0,2,0),(1,0)):
        M=single('l0'); T=named(A,'seed'); f=AP; j=1; path=[]
        residues={q:[value(f,t)%q for t in range(q)] for q in (3,7,49)}
        for step,d in enumerate(word):
            M,T,f,j,record=cycle(M,T,f,j,d,step); path.append(record)
            for q in residues:
                residues[q]=modular_cycle(residues[q],d,q)
                require(residues[q]==[value(f,t)%q for t in range(q)],'modular observer')
        paths.append({'word':word,'steps':path,'Hmod3':count(j,f)%3,
                      'Hmod7':count(j,f)%7,'Hmod49':count(j,f)%49})
    inverse_controls=[select_target(count(1,AP)),select_target(count(1,AP)+2),
                      select_target(count(2,polymul(AP,AP)))]
    require([] in inverse_controls[0]['addresses'],'empty-word inverse missed')
    require(not inverse_controls[1]['addresses'],'nearby hostile inverse changed')
    require([0] in inverse_controls[2]['addresses'],'one-cycle inverse missed')
    first_counts=[inverse_controls[0]['target'],inverse_controls[2]['target']]
    require(gcd(*first_counts)==1,'common-divisor control')
    receipt={'status':'FINITE-EXACT; arbitrary-word theorem is proved separately',
             'small_direct_controls':controls,'run_signatures':{str(j):run_signature(j) for j in (1,2)},
             'prime_quotient':{'rows':quotient,'proper_modules':modules,
                               'automorphism_count':len(aut),'permutations_tested':5040},
             'literal_cycles':paths,'exact_inverse_controls':inverse_controls,
             'consecutive_coprime_counts':first_counts,
             'scope':'No census, scalar injectivity, primality, or all-target claim.',
             'reproduce':'python3 -B 05-knowledge/results/symmetry_port_return_20261113/clock/refill.py'}
    Path(__file__).with_name('receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','literal_words':len(paths),'direct_controls':len(controls),
                      'signature':receipt['run_signatures'],'rigid_quotient':len(aut)},indent=2))

if __name__=='__main__': main()
