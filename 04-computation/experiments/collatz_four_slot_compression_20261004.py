"""Exact four-slot bridges, carry repair, and a new paid Collatz cylinder.

No convergence oracle is used. Infinite claims are proved in the companion
note; enumerations below keep their own finite scopes. Standard library only.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[2]
CHECKS = 0


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ArithmeticError(label)


def coeff(word):
    p, d, b = 1, 1, 0
    for a in word:
        p, d, b = 3*p, d*2**a, 3*b+d
    return p, d, b


def literal(n, word):
    for a in word:
        z = 3*n+1
        actual = (z & -z).bit_length()-1
        check(actual == a, 'literal exact valuation')
        n = z//2**a
    return n


def rank(n):
    if n == 1:
        return 0, 0
    z = n-1
    k = (z & -z).bit_length()-1
    return 3**k*(z//2**k)**2, k


def affine_family(n, period, word):
    """An exact symbolic guard for every nonnegative integer parameter."""
    for a in word:
        check((3*n+1) % 2**(a+1) == 2**a, 'least member exact valuation')
        check(3*period % 2**(a+1) == 0, 'valuation constant on whole progression')
        n, period = (3*n+1)//2**a, 3*period//2**a
    return n, period


def profiles_and_four_core():
    # Four events: P0 precedes P1; Q and R are independent.
    masks = range(16)
    ideals = [s for s in masks if not(s & 2) or s & 1]
    occupancy = lambda s: ((s & 1 != 0)+(s & 2 != 0), int(bool(s & 4)), int(bool(s & 8)))
    box = set(product(range(3), range(2), range(2)))
    check(len(ideals) == 12 and {occupancy(s) for s in ideals} == box, 'ideal/divisor bijection')
    for s,t in product(ideals, repeat=2):
        x,y = occupancy(s), occupancy(t)
        check(((s & t) == s) == all(a <= b for a,b in zip(x,y)), 'order isomorphism')
        check(occupancy(s & t) == tuple(min(a,b) for a,b in zip(x,y)), 'ideal meet')
        check(occupancy(s | t) == tuple(max(a,b) for a,b in zip(x,y)), 'ideal join')
    fibres = Counter(occupancy(s) for s in masks)
    check(Counter(fibres.values()) == {1:8, 2:4}, 'raw-subset quotient keeps collision weights')
    check(occupancy(1) == occupancy(2) and occupancy(1 & 2) != occupancy(1),
          'raw occupancy is not an intersection homomorphism')
    proper = box-{(0,0,0),(2,1,1)}
    squarefree = {x for x in proper if max(x) <= 1}
    prime = {x for x in proper if sum(x) == 1}
    check((len(proper),len(squarefree),len(prime)) == (10,7,3), 'proper divisor balance')
    check(prime < squarefree and len(proper-squarefree) == len(prime), 'nested sets versus count balance')
    # Q4 = C3[TT2,1,1], with P={0,1}, Q={2}, R={3}.
    edges = {(0,1),(0,2),(1,2),(2,3),(3,0),(3,1)}
    paths = [w for w in permutations(range(4)) if all(e in edges for e in zip(w,w[1:]))]
    check(len(paths) == 5, 'strong four-core Hamiltonian count')
    check(sorted(sum((u,v) in edges for v in range(4)) for u in range(4)) == [1,1,2,2], 'four-core scores')
    modules = [pair for pair in combinations(range(4),2)
               if all(((v,pair[0]) in edges) == ((v,pair[1]) in edges)
                      for v in range(4) if v not in pair)]
    check(modules == [(0,1)], 'intrinsic repeated pair')
    check(not any(tuple(0 if v in (0,1) else v-1 for v in w) == (0,1,0,2)
                  for w in paths), '121a order does not falsely become a Q4 path')
    profile_cases = 0
    for q,r in product(range(1,17), repeat=2):
        k = q+r
        all_exponents = set(product(range(k+1),range(q+1),range(2)))
        pp = all_exponents-{(0,0,0),(k,q,1)}
        sk = sum(max(x) < k for x in pp)
        uu = sum(sum(x) == 1 for x in pp)
        check(len(pp) == 2*(k+1)*(q+1)-2, 'general F')
        check(sk == 2*k*(q+1)-1 and uu == 3, 'general k-free and prime counts')
        check(len(pp)-sk-uu == 2*q-2, 'transverse-width defect')
        check(len([1,2]*q+[1]*r+[7]) == 2*q+r+1, 'word event count')
        check(Counter([1,2]*q+[1]*r+[7]) == {1:q+r,2:q,7:1}, 'abelianization profile')
        profile_cases += 1
    return dict(ideal_count=12, raw_subset_count=16, fibre_sizes=dict(Counter(fibres.values())),
                balance=[10,7,3], four_core_paths=paths, profile_cases=profile_cases)


def trace_and_carry():
    trace_cases = 0
    for u,v in product([F(i,j) for i in range(1,9) for j in range(1,6)], repeat=2):
        slots=(u*u,u*v,v*u,v*v)
        z,d=u+v,u*v
        check(sum(slots) == z*z and u*u+v*v == z*z-2*d, 'four-slot trace identity')
        check((u*u)*(v*v) == d*d, 'determinant transport')
        trace_cases += 1
    for i,j in product(range(1,15), repeat=2):
        lam=F(i,j);z=lam+1/lam
        check(lam**2+lam**-2 == z*z-2, 'reciprocal trace semiconjugacy')
        check((lam+2)+1/(lam+2) == z+2-2/(lam*(lam+2)), 'plus-two correction retained')
    swap_cases = 0
    for prefix,suffix in product([(),(1,),(2,1),(3,2,1)],repeat=2):
        for a,b in product(range(1,7),repeat=2):
            w1=prefix+(a,b)+suffix;w2=prefix+(b,a)+suffix
            p1,d1,c1=coeff(w1);p2,d2,c2=coeff(w2)
            check((p1,d1)==(p2,d2), 'swapping keeps two clocks')
            check(c1-c2 == 3**len(suffix)*2**sum(prefix)*(2**a-2**b), 'exact ordered-carry defect')
            swap_cases += 1
    clock_cases=0
    for q,r,a in product(range(1,17),range(3,17),range(4,13)):
        w=[1,2]*q+[1]*r+[a]
        w2=[1,2]*(q+1)+[1]*(r-2)+[a-1]
        p,d,b=coeff(w);p2,d2,b2=coeff(w2)
        check((p,d)==(p2,d2), 'three-to-two clock kernel')
        check(b2-b == 2**(3*q+2)*(3**(r-1)-2**(r-1)) > 0, 'clock fibre retains carry coordinate')
        clock_cases+=1
    check(coeff((1,2,1,1,1,5)) == (729,2048,925), 'first same-clock control')
    check(coeff((1,2,1,2,1,4)) == (729,2048,1085), 'second same-clock control')
    return dict(trace_cases=trace_cases, swap_cases=swap_cases, clock_fibre_cases=clock_cases,
                same_clock_carries=[925,1085])


def paid_four_slots():
    rows=[]
    for a in range(3,65):
        words=sorted(set(permutations((1,1,2,a))))
        start_one=[w for w in words if w[0] == 1]
        check(len(words)==12 and len(start_one)==6, 'four-event distinct orderings')
        carries=[coeff(w)[2] for w in start_one]
        check(max(carries)==45+14*2**a, 'max carry by descending suffix')
        check(7*(2**(a+4)-81) > max(carries), 'uniform n-at-least-seven payment bound')
        for w in start_one:
            p,d,b=coeff(w);modulus=2*d
            source=((d-b)*pow(p,-1,modulus))%modulus
            check(source>=7 and source%4==3, 'nonvacuous legal source and original-rank stratum')
            child,step=affine_family(source,modulus,w)
            check(source>child>0 and modulus>step>0, 'all-height paid family')
            for t in (0,1,2,17,10**6):
                n=source+modulus*t;m=child+step*t
                check(literal(n,w)==m and rank(m)<rank(n), 'literal family and inherited-rank controls')
            if a<=4:rows.append(dict(word=w,source=source,period=modulus,child=child,child_period=step))
    check(literal(123,(1,2,1,2))==157, 'hostile: lowering exit threshold to two is unpaid')
    check(literal(151,(1,1,10))==1 and literal(1,(2,))==1,
          'general permutation arithmetic may require first-hit truncation')
    check(affine_family(187,256,(1,2,1,3))==(119,162), 'new entire weak-exit cylinder')
    check(affine_family(7,256,(1,1,2,3))==(5,162), 'permuted entire weak-exit cylinder')
    for t in range(4096):
        for source,child,w in ((187,119,(1,2,1,3)),(7,5,(1,1,2,3))):
            n=source+256*t;m=child+162*t
            check(literal(n,w)==m and 0<m<n and rank(m)<rank(n), 'new source-cell census')
    phases=[k for k in range(64) if pow(3,k,256)==187]
    check(phases==[27] and pow(3,32,256)!=1 and pow(3,64,256)==1, 'unique exponent phase')
    for t in range(101):
        k=27+64*t
        check(pow(3,k,256)==187, 'modular infinite-family control')
        if t<8:
            n=3**k
            check(literal(n,(1,2,1,3))==(81*n+85)//128<n, 'expanded power control')
    return dict(exit_range=[3,64], orderings_per_exit=6, rows=rows,
                full_motif_odd_density='3/32', new_sources=[[187,256],[7,256]], new_children=[[119,162],[5,162]],
                new_odd_density='1/64', new_exponent=[27,64], added_relative_exponent_density='1/8')


def permuted_bank(H=64):
    rows=[]; lower=F(0)
    for h in range(2,H+1):
        a=3
        if h!=2:
            while 2**(h+2+a)<=2*3**(h+2):a+=1
        word=[1]*h+[2,a]
        p,d,b=coeff(word)
        check(b==3**(h+2)-2**(h+1), 'ordered all-ones-first carry')
        pp,dd,bb=coeff([1,2]+[1]*(h-1)+[a])
        check((pp,dd)==(p,d) and bb-b==12*(3**(h-1)-2**(h-1)),
              'same balanced profile transports clocks and changes carry explicitly')
        check(h==2 or 2*p<d, 'safe threshold beyond sharp four-slot base')
        source=(-b*pow(p,-1,d))%d
        check(source%8==7, 'permuted bank has no powers of three')
        holes=[]; density=F(2,d)
        if a<=6:
            pp,dd,bb=coeff([1]*h+[2,6])
            hole=((dd-bb)*pow(pp,-1,2*dd))%(2*dd)
            check(hole%d==source and (2*dd)%d==0, 'exact old-bank overlap contained')
            holes=[dict(residue=hole,modulus=2*dd)]
            density-=F(1,dd)
        for t in (0,1,2,17,10**6):
            n=source+d*t;x=n
            x=literal(x,[1]*h+[2])
            z=3*x+1;actual=(z & -z).bit_length()-1;m=z//2**actual
            check(actual>=a and 0<m<n and rank(m)<rank(n), 'permuted coarse cylinder literal payment')
        check((d-p)*source>b, 'formal least-source inequality for saved permuted cells')
        rows.append(dict(h=h,L=a,residue=source,modulus=d,holes=holes,odd_density=str(density)))
        lower+=density
    check([x['h'] for x in rows if x['holes']]==[2,3,4,5,6], 'all exact-six old-bank overlaps')
    return dict(H=H,rows=rows,lower=str(lower),tail_bound=str(F(1,2**(H+4))),decimal=float(lower))


def coverage(permuted):
    oldpath=ROOT/'05-knowledge/results/collatz_paid_portrait_controllers_20261004.json'
    old=json.loads(oldpath.read_text());c=old['coverage']
    cells=old['rational_bank']['rows']+old['switched_bank']['rows']
    cells += [x for x in old['safe_cylinders'] if x['c']==17]
    cells += [old['switched_bank']['checkpoint']]
    for source in (187,7):
        for x in cells:
            check((source-x['residue']) % math.gcd(256,x['modulus']) != 0, 'disjoint from saved prior binary cells')
    check((187-7)%256 != 0, 'new cells mutually disjoint')
    check(all(pow(3,k,256)!=7 for k in range(64)), 'second cell has no power-three phase')
    for row in permuted['rows']:
        for x in cells:
            check((row['residue']-x['residue']) % math.gcd(row['modulus'],x['modulus']) != 0,
                  'permuted bank disjoint from preceding additions')
        check((row['residue']-187)%math.gcd(row['modulus'],256)!=0, 'permuted bank disjoint from1213')
    for i,row in enumerate(permuted['rows']):
        for other in permuted['rows'][:i]:
            check((row['residue']-other['residue'])%math.gcd(row['modulus'],other['modulus'])!=0,
                  'permuted rows distinguished by initial ones')
    frozen=json.loads((ROOT/'05-knowledge/results/collatz_binary_ternary_guard_fusion_20261004.json').read_text())
    check([(x['s'],x['source_word'][:-1]) for x in frozen['debt_rows'] if x['e']==1]
          == [(1,[1,6]),(2,[6])], 'inherited binary16 disjointness prefixes')
    # Reconstruct the finite strengthened-origin prefix and tail independently.
    origin=[(1,2)]
    for k in range(1,32):
        d=1
        while 3**d<=2**(d+2*k-1):d+=1
        residue=-((4**k+2)//3)*pow(4**k,-1,3**d)%3**d
        if not any(residue%3**j==rr for j,rr in origin):origin.append((d,residue))
    tl=sum(F(1,3**d) for d,_ in origin)
    d=1
    while 3**d<=2**(d+63):d+=1
    tu=tl+F(27,26*3**d)
    check([str(tl),str(tu)]==c['strengthened_origin_density_interval'], 'frozen ternary interval reproduced')
    oldres=list(map(F,c['strengthened_fused_residual_interval']))
    gainlo=F(1,128)+F(permuted['lower'])
    gainhi=gainlo+F(permuted['tail_bound'])
    residual=(oldres[0]-gainhi*(1-tl),oldres[1]-gainlo*(1-tu))
    exponents=tuple(F(x)+F(1,8) for x in c['repaired_fraction_of_exponents_3_mod8_interval'])
    binary=(F(c['new_binary_density_interval'][0])+gainlo,F(c['new_binary_density_interval'][1])+gainhi)
    return dict(previous_artifact_sha256=hashlib.sha256(oldpath.read_bytes()).hexdigest(),
                compared_prior_cells=len(cells), new_binary_interval=list(map(str,binary)),
                incremental_binary_interval=list(map(str,(gainlo,gainhi))),incremental_binary_decimal=float(gainlo),
                named_residual_interval=list(map(str,residual)), named_residual_decimal=float(residual[0]),
                exponent_3_mod8_coverage_interval=list(map(str,exponents)),
                exponent_3_mod8_coverage_decimal=float(exponents[0]),
                qualification='Paid dependencies in the specified banks; universal root coverage OPEN.')


def weaker_exit_probe():
    """A finite conjecture probe, deliberately separate from the all-height proof."""
    fail=[]
    for q,r in product(range(1,129),repeat=2):
        p=3**(2*q+r+1)
        a=max(2,p.bit_length()-3*q-r)
        d=2**(3*q+r+a)
        b=5*p-12*3**r*2**(3*q)-2**(3*q+r+1)
        check(coeff([1,2]*q+[1]*r+[a])==(p,d,b), 'independent closed switched coefficients')
        check(d>p and (a==2 or d//2<=p), 'minimum contracting exit')
        source=(-b*pow(p,-1,d))%d
        check(source%2==1, 'coarse source is odd')
        if (d-p)*source<=b:fail.append([q,r,a,source])
    return dict(status='FINITE-EXACT; all-parameter minimum-exit claim OPEN',
                q_range=[1,128],r_range=[1,128],cases=128**2,least_source_payment_failures=fail)


def completed_family():
    controls=[]
    for h,t in product(range(2,9),(0,1,3)):
        a=3**(h+1)*(2*t+1)-1
        denominator=3**(h+2)
        numerator=2**(h+1)*(2**(a+1)+1)
        check(numerator%denominator==0, 'completed family integrality')
        n=numerator//denominator-1
        check(n>1 and n%8==7, 'completed family positive source')
        x=n
        for letter in [1]*h+[2,a]:
            check(x>1, 'completed family first-hit boundary')
            x=literal(x,(letter,))
        check(x==1, 'completed family reaches root')
        controls.append(dict(h=h,parameter=t,exit=a,source_bits=n.bit_length(),
                             source=n if h<=3 else None))
    check(controls[0]['source']==13256071, 'smallest displayed completed-family source')
    return dict(status='All-height constructed completed family, not arbitrary-source coverage',controls=controls)


def main():
    permuted=permuted_bank()
    out=dict(status='PROVED scoped bridges and new paid family; universal Collatz OPEN',
             profiles=profiles_and_four_core(),carry=trace_and_carry(),paid=paid_four_slots(),
             permuted_bank=permuted,coverage=coverage(permuted),weaker_exit_probe=weaker_exit_probe(),
             completed_family=completed_family())
    out['checks']=CHECKS
    base=ROOT/'05-knowledge/results/collatz_four_slot_compression_20261004'
    base.with_suffix('.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    report='\n'.join([
        out['status'],
        'PROFILE: 12 ideals / 16 raw subsets; (F,S,U)=(10,7,3); width defect=2q-2.',
        'CARRY: exact adjacent-swap law; same clocks retain carries925 and1085.',
        'PAID: all six initial-one orderings of {1,1,2,a}, every a>=3 (proved in note).',
        'NEW CELL: 187+256t ->119+162t; exponents27+64t; added relative exponent density1/8.',
        'SECOND CELL: 7+256t ->5+162t; two disjoint cells add odd-relative binary density1/64.',
        f"ALL-DEPTH PERMUTATION BANK: density {permuted['decimal']:.16f}, exact-six holes removed.",
        f"TOTAL BINARY INCREMENT: {out['coverage']['incremental_binary_decimal']:.16f}",
        f"COMBINED EXPONENT DOMAIN: {out['coverage']['exponent_3_mod8_coverage_decimal']:.16f}",
        f"NAMED BANK RESIDUAL: {out['coverage']['named_residual_decimal']:.16f}",
        'COMPLETED AT EVERY DEPTH: a=3^(h+1)(2t+1)-1, n=2^(h+1)(2^(a+1)+1)/3^(h+2)-1; word1^h2a ->1.',
        'FINITE MINIMUM-EXIT PROBE: '+json.dumps(out['weaker_exit_probe'],sort_keys=True),
        f'PASS: {CHECKS} explicit checks.',
    ])+'\n'
    base.with_suffix('.out').write_text(report)
    print(report,end='')


if __name__=='__main__':
    main()
