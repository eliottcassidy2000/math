#!/usr/bin/env python3
"""Exact entry-coset maxima; inherited rate formula, independent finite orbits."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import importlib.util,json
ROOT=Path(__file__).resolve().parents[3]
pin=ROOT/'05-knowledge/results/session114_continuation_symmetry_20261030/collatz/check_periodic_rates.py'
spec=importlib.util.spec_from_file_location('inherited_rates',pin);r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
def need(b,msg):
    if not b:raise RuntimeError(str(msg))
def matrix(w):
    W=r.coeff(w);t,h=r.params(W)
    prefixes=[r.coeff(w[:j]) for j in range(1,len(w)+1)]
    suffixes=[r.coeff(w[j:]) for j in range(1,len(w)+1)]
    D=[[r.cutoff(W,r.compose(v,u)) for v in prefixes] for u in suffixes]
    rates=[sum((r.mass(t,h,x) for x in row),F(0)) for row in D]
    return D,rates
def main():
    count=stall=0;maximal=[];worst=F(0);worst_word=None
    for w in r.words(14):
        D,rs=matrix(w);best=max(rs);ratio=best/sum(w);count+=1
        need(ratio<=F(5,8),('entry-uniform finite controller failure',w,rs))
        if ratio==F(5,8):maximal.append(w)
        if set(w)!={1} and ratio>worst:worst=ratio;worst_word=w
        if r.params(r.coeff(r.primitive_root(w)))[1]>1:stall+=1
        need(all(D[i][j]==D[j][i] for i in range(len(w)) for j in range(len(w))),'symmetric coset depth')
    dominance=0;ultrametric=0
    for w in r.words(8):
        D,rs=matrix(w)
        for i,j,k in product(range(len(w)),repeat=3):
            need(D[i][k]>=min(D[i][j],D[j][k]),'ultrametric depth');ultrametric+=1
        for u in [()]+list(r.words(5)):
            depths=[r.phase_rate(w,u,w[:j])[1] for j in range(1,len(w)+1)]
            j=max(range(len(w)),key=depths.__getitem__)
            need(all(y>=x for x,y in zip(depths,D[j])),('componentwise dominance',w,u));dominance+=1
    # Single-outlier theorem: prefix powers give floor-subadditive row sums.
    outlier=0
    for a,b,L in product(range(1,9),range(1,9),range(1,21)):
        if a==b:continue
        w=(a,)*L+(b,);D,rs=matrix(w)
        need(max(rs)==r.rate(w),'single-outlier endpoint maximizes every entry')
        outlier+=1
    # The general repeated-letter equality grammar, with literal intersection moduli.
    equality_checks=0
    for a,b,L in product(range(1,5),range(1,5),range(1,5)):
        if a==b:continue
        w=(a,)*L+(b,);D,rs=matrix(w);maximum=max(rs);ells=[]
        for E in range(1,7):
            q=5**E;H=r.orbit(r.coeff(w),q);gm,gb=r.pair(r.coeff((a,)),q)
            x,y=1,0;k=0
            while True:
                x,y=gm*x%q,(gm*y+gb)%q;k+=1
                if (x,y) in H:break
            ells.append(k)
            if k>L:break
        need(ells[-1]>L,'equality grammar finite separating depth')
        accepted=[j for j in range(L+1) if all(j%e<=L%e for e in ells)]
        need(all((rs[j-1] if j else rs[-1])==maximum for j in accepted),'maximal phase classes')
        need(sum(x==maximum for x in rs)==len(accepted),'complete maximal phase count')
        entries=[()]+list(r.words(6))+[w[j:]+w for j in range(L+1)]
        for u in entries:
            expected=False
            for j in accepted:
                joined=(a,)*j+u;B0,rem=divmod(len(joined),len(w))
                if rem==0 and joined==w*B0:expected=True
            need((r.rate(w,u)==maximum)==expected,('full equality grammar',w,u,ells));equality_checks+=1
    # Independently enumerate every affine entry and literal orbit at finite precision.
    finite=0
    def valmod(x,E):
        h=0
        while h<E and x%5==0:h+=1;x//=5
        return h
    for w in r.words(5):
        W=r.coeff(w);D,rs=matrix(w)
        for E in (1,2):
            q=5**E;orb=r.orbit(W,q);T=len(orb)
            predicted=max(sum((r.mass(*r.params(W),min(E,x)) for x in row),F(0)) for row in D)
            best=-1
            for m,b in product(range(1,q),range(q)):
                if m%5==0:continue
                initial=(m,b);x,y=initial;total=0
                for k in range(T):
                    for a in w:
                        inv=pow(2**a,-1,q);x,y=3*x*inv%q,(3*y+1)*inv%q
                        total+=min(valmod(x-1,E),valmod(y,E))
                need((x,y)==initial,'literal group clock');best=max(best,total);finite+=1
            need(F(best,T)==predicted,('full finite entry maximum',w,E))
    # Literal integer-prefix charge and the explicit original-entry-clock error.
    discrepancy=0;max_excess=F(0)
    for w0 in r.words(5):
        _,means=matrix(w0);Rmax=max(means);A=sum(w0);d=len(w0)
        for u0 in [()]+list(r.words(4)):
            U=sum(u0)
            for N in (1,2,5,25):
                P,Q,B=1,1,0
                for aa in u0:P,Q,B=3*P,(2**aa)*Q,3*B+Q
                S=0
                for aa in w0*N:
                    P,Q,B=3*P,(2**aa)*Q,3*B+Q
                    from math import gcd
                    g=gcd(abs(P-Q),B);h=0
                    while g%5==0:g//=5;h+=1
                    S+=h
                target=2*3**(U+4*A)*(U+N*A+1);K=0;power=1
                while power<target:power*=5;K+=1
                need(S<=N*Rmax+d*K,'explicit original-clock discrepancy')
                max_excess=max(max_excess,S-N*Rmax);discrepancy+=1
    w=(1,3,1,1,3);D,rs=matrix(w)
    need(r.rate(w)==F(5,8) and max(rs)==F(369,200),'entry changes mean hostile')
    # Pointwise dominance is false even when the infinite depths agree.
    need(r.val(3**2-2**2)==1 and r.val(3-2)==0,'one-block time-origin hostile')
    out=dict(status='finite exact controls for proved entry-dominance theorem',
             controllers_clock_at_most_14=count,order_stall_controllers=stall,
             equality_controllers=len(maximal),equality_only_all_ones=all(set(w)=={1} for w in maximal),
             maximum_nonones_clock_rate=str(worst),maximum_nonones_word=worst_word,
             componentwise_entry_checks=dominance,ultrametric_checks=ultrametric,
             single_outlier_controls=outlier,literal_finite_entries=finite,
             maximizing_entry_grammar_checks=equality_checks,integer_clock_discrepancy_checks=discrepancy,maximum_tested_excess=str(max_excess),
             variable_entry_example=dict(controller=w,empty_mean='5/8',suffix_means=list(map(str,rs)),
                                         depths=[[str(x) for x in row] for row in D]))
    Path(__file__).with_name('entry_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
