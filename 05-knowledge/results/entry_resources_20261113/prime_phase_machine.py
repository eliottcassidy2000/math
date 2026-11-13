"""Exact exponent sieve as a finite graph cellular automaton and unary DFA."""
from pathlib import Path
from math import gcd,lcm,isqrt
from fractions import Fraction
from itertools import combinations
import json,gzip,sys
HERE=Path(__file__).resolve().parent
def need(p,m):
    if not p:raise RuntimeError(m)
def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def factor(n):
    out={};p=2
    while p*p<=n:
        while n%p==0:out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=out.get(n,0)+1
    return out
def merge(a,m,b,n):
    g=gcd(m,n)
    if (b-a)%g:return None
    t=((b-a)//g)*pow(m//g,-1,n//g)%(n//g)
    return (a+m*t)%lcm(m,n),lcm(m,n)
def run():
    primes=[p for p in range(2,100) if prime(p)];active=[];inactive=[]
    for p in primes:
        if p in [2,3,11,23]:inactive.append(p);continue
        orbit=[];x=1
        while x not in orbit:orbit.append(x);x=2*x%p
        zeros=[k for k,v in enumerate(orbit) if (99*v+23)%p==0]
        if zeros:need(len(zeros)==1,'single discrete-log fiber');active.append({'p':p,'period':len(orbit),'residue':zeros[0]})
        else:inactive.append(p)
    inherited_period=lcm(*(r['period'] for r in active))
    redundant=[r for r in active if any(r is not s and r['period']%s['period']==0 and r['residue']%s['period']==s['residue'] for s in active)]
    essential=[r for r in active if r not in redundant]
    L=lcm(*(r['period'] for r in essential));fac=factor(L);mods=[p**e for p,e in fac.items()]
    def allowed(k):return all(k%r['period']!=r['residue'] for r in active)
    separators=[]
    for p in fac:
        shift=L//p;k=next(k for k in range(1,100001) if allowed(k)!=allowed(k+shift))
        separators.append({'prime':p,'proper_period_candidate':shift,'k':k,'accept_k':allowed(k),'accept_shifted':allowed(k+shift)})
    density=Fraction(1)
    for size in range(1,len(active)+1):
        for subset in combinations(active,size):
            state=(0,1)
            for r in subset:
                state=merge(*state,r['residue'],r['period'])
                if state is None:break
            if state is not None:density+=Fraction((-1)**size,state[1])
    need((L*density).denominator==1,'integer complete-period survivor count')
    clauses=[]
    for r in essential:
        terms=[]
        for i,m in enumerate(mods):
            d=gcd(m,r['period'])
            if d>1:terms.append({'cycle':i,'mask':sum(1<<v for v in range(m) if v%d==r['residue']%d)})
        clauses.append(terms)
    state=[1]*len(mods);checks=0
    for k in range(50001):
        output=not any(all(state[t['cycle']]&t['mask'] for t in clause) for clause in clauses)
        literal=all((99*pow(2,k,p)+23)%p for p in primes)
        need(output==allowed(k)==literal,'CA, exponent classes and direct divisibility agree');checks+=1
        state=[((s<<1)&((1<<m)-1))|(s>>(m-1)) for s,m in zip(state,mods)]
        need(all(s.bit_count()==1 for s in state),'one marker on every local cycle')
    M=lcm(L,418)
    need(prime(419) and pow(2,418,419)==1 and all(pow(2,d,419)!=1 for d in range(1,418) if 418%d==0),'419 order control')
    need(allowed(2) and M%L==0 and M%418==0,'composite-survivor progression')
    old=json.loads((HERE.parent/'endpoint_coverage_compiler_20261111/wheel.json').read_text(encoding='utf8'))
    need(old['active']==active and old['wheel_period']==inherited_period and old['survivor_residues']==int(inherited_period*density),'independent recovery of inherited sieve')
    out={'scope':'Inherited exact25-prime divisibility wheel independently reconstructed; redundant clauses removed, minimal-period separators, binary permutation CA on disjoint directed cycles. Not a primality test or one-dimensional Rule30 claim.','active':active,'inactive':inactive,'redundant_clauses':redundant,'essential_clauses':essential,'inherited_valid_period':inherited_period,'minimal_unary_DFA_states':L,'density':[density.numerator,density.denominator],'survivors':int(L*density),'cycle_lengths':mods,'binary_CA_cells':sum(mods),'readout_clauses':clauses,'minimal_period_separators':separators,'literal_controls':checks,'composite_progression':{'start':2,'period':M,'proper_divisor':419,'t_min':1}}
    tag='' if __debug__ else '.optimized';(HERE/('prime_phase_machine'+tag+'.json.gz')).write_bytes(gzip.compress((json.dumps(out,separators=(',',':'))+'\n').encode(),mtime=0))
    print({k:v for k,v in out.items() if k not in ['scope','readout_clauses']},flush=True)
if __name__=='__main__':run()
