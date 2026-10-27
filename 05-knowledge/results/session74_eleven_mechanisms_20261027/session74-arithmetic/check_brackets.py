#!/usr/bin/env python3
"""Exact Session74 bracket arithmetic; stdlib, exceptions remain under python -O."""
from fractions import Fraction
from math import gcd, isqrt
import json

checks = 0
def eq(got, want, label):
    global checks
    checks += 1
    if got != want:
        raise RuntimeError(f'{label}: {got!r} != {want!r}')

def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n)+1))

def bracket(n):
    if n < 2:
        raise ValueError('the singleton {1} is excluded')
    return (isqrt(n-1)+1)//2

def same(n, r):
    return r.numerator*n <= r.denominator*(2*bracket(n)+1)**2

def intervals(q):
    return [(m, (2*m-1)**2+1, (q-1)*(2*m+1)**2//q)
            for m in range(1, 2*q-1)]

def interval_count(q):
    return sum(hi-lo+1 for _,lo,hi in intervals(q))

def residue_count(q):
    tax = sum((-(2*n+1)**2) % q for n in range(q))
    return Fraction(16*q*q-24*q+8,3)-Fraction(2*tax,q)

def legendre(a,p):
    if a % p == 0:
        return 0
    z=pow(a,(p-1)//2,p)
    if z not in (1,p-1):
        raise RuntimeError('nonprime Legendre modulus')
    return 1 if z==1 else -1

def character_sum(p):
    return sum(a*legendre(a,p) for a in range(1,p))

def pp_count(p,e):
    q=p**e
    z=p**(e//2)
    geom=(p**((e+1)//2)-1)//(p-1)
    return Fraction(16*q*q-27*q+8,3)+z+Fraction(2*character_sum(p)*geom,p)

# Baseline exact classification; proof in README establishes the all-height cutoff.
original=[n for n in range(2,10001) if same(n,Fraction(2))]
eq(original,[2,3,4,10,11,12],'doubling integer set')
eq([n for n in original if prime(n)],[2,3,11],'doubling prime set')
pairs=[(k,p) for p in range(2,100) if prime(p) for k in range(2,10)
       if same(p,Fraction(k))]
eq(pairs,[(2,2),(3,2),(4,2),(2,3),(3,3),(2,11)],'prime multiplier pairs')
plateau=[]
for r in [Fraction(25,13),Fraction(2),Fraction(25,11),Fraction(25,11)+Fraction(1,10000)]:
    found=[n for n in range(2,10000) if prime(n) and same(n,r)]
    plateau.append({'r':str(r),'primes':found})
eq(plateau[0]['primes'],[2,3,11,13],'lower endpoint not in plateau')
eq(plateau[1]['primes'],[2,3,11],'interior plateau')
eq(plateau[2]['primes'],[2,3,11],'upper endpoint in plateau')
eq(plateau[3]['primes'],[2,3],'above upper endpoint')

# Two exact paths for all q through512; direct element-membership is a third path through64.
for q in range(2,513):
    rows=intervals(q)
    eq(all(hi>=lo for _,lo,hi in rows),True,f'nonempty q={q}')
    m=2*q-1
    eq((q-1)*(2*m+1)**2//q < (2*m-1)**2+1,True,f'cutoff q={q}')
    eq(Fraction(interval_count(q)),residue_count(q),f'floor/residue q={q}')
    if q<=64:
        brute=[n for n in range(2,(4*q-3)**2+1) if same(n,Fraction(q,q-1))]
        literal=[n for _,lo,hi in rows for n in range(lo,hi+1)]
        eq(brute,literal,f'element/interval q={q}')

# Prime-power formula and independent unit/nonunit square-residue recurrence.
pp_rows=[]
for p in range(3,102,2):
    if not prime(p):
        continue
    qr=sorted({x*x%p for x in range(1,p)})
    S=character_sum(p)
    eq(S,2*sum(qr)-p*(p-1)//2,f'QR/character p={p}')
    for e in range(1,7):
        q=p**e
        if q>20000:
            break
        R=sum(x*x%q for x in range(q))
        Z=sum(x*x%q==0 for x in range(q))
        Ru=sum(x*x%q for x in range(q) if x%p)
        eq(Z,p**(e//2),f'zero stratum p={p},e={e}')
        eq(Fraction(Ru),Fraction(p-1,2)*p**(2*e-1)+p**(e-1)*S,
           f'unit lift p={p},e={e}')
        if e>=2:
            Rsmall=sum(x*x%(p**(e-2)) for x in range(p**(e-2)))
            eq(R,Ru+p**3*Rsmall,f'nonunit recurrence p={p},e={e}')
        geom=(p**((e+1)//2)-1)//(p-1)
        eq(Fraction(R,q),Fraction(q-Z,2)+Fraction(S*geom,p),
           f'weighted lift p={p},e={e}')
        eq(Fraction(interval_count(q)),pp_count(p,e),f'prime-power count p={p},e={e}')
        eq(residue_count(q),pp_count(p,e),f'residue/prime-power p={p},e={e}')
        pp_rows.append({'p':p,'e':e,'q':q,'S_p':S,'zero_roots':Z,'count':int(pp_count(p,e))})

# q=11 source arithmetic and prime census are distinct claims.
q11rows=intervals(11)
q11primes=[n for _,lo,hi in q11rows for n in range(lo,hi+1) if prime(n)]
eq(len(q11primes),100,'finite q11 prime census')
eq(character_sum(11),-11,'11 weighted defect')
eq(character_sum(23),-69,'hostile to universal minus-two correction')
eq(pp_count(11,1),548,'q11 count')
eq(pp_count(11,2),77008,'q121 count')
eq(pp_count(11,3),9436336,'q1331 count')

# CRT image is the actual inherited K, with all source multiplicities and weights.
q=341
K={pow(4,a,q)*pow(81,b,q)%q for a in range(5) for b in range(15)}
unit_image={x*x%q for x in range(q) if gcd(x,q)==1}
eq(K,unit_image,'joint-clock/squaremap support')
eq(len(K),75,'joint-clock group order')
strata=[]
for name,condition in [
    ('zero',lambda x:x==0),
    ('order11_axis',lambda x:x!=0 and x%31==0),
    ('order31_axis',lambda x:x!=0 and x%11==0),
    ('units',lambda x:gcd(x,q)==1)]:
    xs=[x for x in range(q) if condition(x)]
    hist={}
    for x in xs:
        hist[x*x%q]=hist.get(x*x%q,0)+1
    strata.append({'stratum':name,'root_count':len(xs),'image_count':len(hist),
                   'fiber_sizes':sorted(set(hist.values())),
                   'square_residue_sum':sum(x*x%q for x in xs),
                   'ceiling_tax':sum((-x*x)%q for x in xs)})
eq([r['root_count'] for r in strata],[1,10,30,300],'341 source strata')
eq([r['image_count'] for r in strata],[1,5,15,75],'341 image strata')
eq([r['fiber_sizes'] for r in strata],[[1],[2],[2],[4]],'341 multiplicities')
eq([r['ceiling_tax'] for r in strata],[0,2046,4092,53196],'341 weighted taxes')
eq(interval_count(341),617092,'341 exact count')
units_only=Fraction(16*q*q-24*q+8,3)-Fraction(2*strata[-1]['ceiling_tax'],q)
eq(units_only,617128,'units-only hostile')
eq(units_only-interval_count(341),36,'units-only overcount')

# Quotient failures: neither squareclass nor prime tag alone sees bracket escape.
eq(same(11,Fraction(2)),True,'11 positive control')
eq(same(44,Fraction(2)),False,'squareclass hostile44')
eq(same(5,Fraction(2)),False,'generic basis hostile5')

out={
 'status':'PROVED elementary formulas; FINITE-EXACT checks, not a prime-count theorem',
 'checks':checks,
 'universes':{'all_integer_q':[2,512],'direct_element_q':[2,64],
              'prime_power_base_primes':[3,101],'prime_power_max_exponent':6,
              'prime_power_max_modulus':20000,'prime_power_cases':len(pp_rows)},
 'original_doubling':{'integers':original,'primes':[2,3,11],'pairs':pairs,'plateau_controls':plateau},
 'q11':{'count':548,'prime_count_finite':100,'prime_list':q11primes,
         'quadratic_residues':[1,3,4,5,9],'weighted_character_sum':-11,
         'terminal_bracket':q11rows[-1]},
 'prime_power_rows':pp_rows,
 'q341':{'count':617092,'K_order':75,'strata':strata,'units_only_wrong_count':int(units_only),
         'missing_axis_correction':36},
 'hostiles':{'squareclass':[11,44],'generic_basis_replacement':[2,3,5],
             'p23_weighted_character_sum':-69,'p23_polynomial_correction':-6}
}
print(json.dumps(out,indent=2,sort_keys=True))
