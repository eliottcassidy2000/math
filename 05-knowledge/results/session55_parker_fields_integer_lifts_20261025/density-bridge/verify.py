#!/usr/bin/env python3
"""Exact finite controls; standard library, portable source-pinned arithmetic."""
import argparse, hashlib, importlib.util, itertools, json
from fractions import Fraction as Q
from pathlib import Path
D=Path(__file__).resolve().parent

def check(ok,msg):
    if not ok: raise RuntimeError(msg)

def module(name,file):
    spec=importlib.util.spec_from_file_location(name,D/'inputs'/file)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

Fmod=module('field_source','verify.py')
Rmod=module('ring_source','ring_verify.py')

def run():
    pins=json.loads((D/'source-pins.json').read_text())
    for pin in pins:
        check(hashlib.sha256((D/pin['copy']).read_bytes()).hexdigest()==pin['sha256'],'source hash '+pin['copy'])
    known={r['q']:r for r in json.loads((D/'inputs/small_counts.json').read_text())}
    fields=[];slice_count=0
    for q,p,n in Fmod.prime_powers(126):
        F=Fmod.Field(p,n);square=F.squares-{0}
        chi=lambda a: 0 if a==0 else (1 if a in square else -1)
        A={v for v in range(q) if F.add(1,v) in square and F.sub(1,v) in square}
        check(len(A)==Q(q-chi(F.neg(1)),4)-Q(1+chi(F.add(1,1)),2),'exact A count')
        E={0,1,F.neg(1),F.add(1,1),F.neg(F.add(1,1))};G=A-E
        check(abs(Q(len(G))-Q(q,4))<=Q(25,4),'G error')
        good=0
        for v in G:
            roots={1,F.neg(1),F.add(1,v),F.neg(F.add(1,v)),F.sub(1,v),F.neg(F.sub(1,v))}
            check(len(roots)==6,'six distinct roots')
            N=sum(Fmod.valid(F,g) and 0 not in g for u in range(q) for g in [Fmod.magic(F,1,u,v)])
            # |64N-q| <= 655+114sqrt(q), without floating point.
            check(max(0,abs(64*N-q)-655)**2 <= 114**2*q,'two-sided slice error')
            good+=N;slice_count+=1
        norm=Fmod.normalized(F,True);C=norm['center_one_count'];Z=(q-1)//2*norm['center_zero_u1_count']
        check(C==known[q]['center_one'],'source center-one count')
        check(norm['center_zero_u1_count']==known[q]['center_zero_u1'],'source zero-center count')
        check(0<=C-good<=13*q,'exceptional/zero center-one remainder')
        check(0<=Z<=q*q,'zero-center bound')
        # |C-q²/256| <= 2q^(3/2)+24q, again exactly squared.
        error=max(Q(0),abs(Q(C)-Q(q*q,256))-24*q)
        check(error*error<=4*q**3,'center-one error')
        total=(q-1)//2*C+Z
        error=max(Q(0),abs(Q(total)-Q(q**3,512))-14*q*q)
        check(error*error<=q**5,'full field error')
        fields.append({'q':q,'p':p,'degree':n,'admissible_v':len(A),'good_v':len(G),'good_slice_arrays':good,'center_one':C,'center_zero':Z,'total':total})
    rings=[]
    for p,K in [(2,9),(3,5),(5,3),(7,2),(29,1)]:
        counts={-1:0,0:0}
        for k in range(1,K+1):
            m=p**k;square={x*x%m for x in range(m)}
            N=P=U=0
            for c,a,b in itertools.product(square,repeat=3):
                g=Rmod.grid(c,a,b,m)
                if len(set(g))!=9 or not all(x in square for x in g):continue
                check(Rmod.magic(g,m),'eight direct ring line sums')
                N+=1;P+=any(x%p for x in g);U+=all(x%p for x in g)
                if not any(x%p for x in g):
                    check(k>=2 and all(x%(p*p)==0 for x in g),'square-content divisibility')
                    h=m//(p*p);reduced=tuple(x//(p*p) for x in g)
                    check(len(set(reduced))==9 and Rmod.magic(reduced,h),'descent injectivity and magic')
                    small={x*x%h for x in range(h)}
                    check(all(x in small for x in reduced),'descent square predicate')
            check(N==P+counts[k-2],'general two-step recurrence')
            if p in (2,3):
                check(P==U and N==Rmod.count(p,k),'special all-unit branch and exact source formula')
            counts[k]=N
            rings.append({'p':p,'k':k,'modulus':m,'N':N,'primitive_at_least_one_unit':P,'all_unit':U,'nonunit':N-P})
    witness=(1,22,6,5,0,24,23,7,28)
    check(Rmod.magic(witness,29) and len(set(witness))==9 and set(witness)<={x*x%29 for x in range(29)},'mixed primitive hostile witness')
    check(any(r['p']==29 and r['primitive_at_least_one_unit']>r['all_unit'] for r in rings),'primitive does not imply all-unit')
    constants={'field':Q(1,512),'dyadic_primitive':Q(1,512),'ternary_primitive':Q(1,27)}
    constants['dyadic_total']=constants['dyadic_primitive']/(1-Q(1,2**6))
    constants['ternary_total']=constants['ternary_primitive']/(1-Q(1,3**6))
    check(constants['dyadic_total']==Q(1,504) and constants['ternary_total']==Q(27,728),'ring limits')
    constants['dyadic_field_ratio']=constants['dyadic_total']/constants['field']
    constants['ternary_field_ratio']=constants['ternary_total']/constants['field']
    check(constants['dyadic_field_ratio']==Q(64,63) and constants['ternary_field_ratio']==Q(1728,91),'comparison constants')
    convergence=[]
    for p,limit in [(2,Q(1,504)),(3,Q(27,728))]:
        previous=None
        for k in (16,24,40):
            density=Q(Rmod.count(p,k),p**(3*k));error=abs(density-limit)
            check(previous is None or error<previous,'finite convergence control')
            previous=error
            convergence.append({'p':p,'k':k,'density':str(density),'error':str(error)})
    return {'scope':'Labelled distinct square VALUES, zero allowed. Exact finite controls support separate all-level proofs; no novelty claim.','source_pins_checked':len(pins),'field_slices_checked':slice_count,'field_counts':fields,'ring_counts':rings,'constants':{k:str(v) for k,v in constants.items()},'finite_convergence_controls':convergence,'mixed_primitive_mod29':witness}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    out=json.dumps(run(),sort_keys=True,indent=2)+'\n'
    if args.write:(D/'certificate.json').write_text(out)
    else:check((D/'certificate.json').read_text()==out,'certificate replay')
    print('PASS: source pins, field slices and counts, general content recurrence, hostile controls, exact density constants')
