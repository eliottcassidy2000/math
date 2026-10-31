#!/usr/bin/env python3
"""Session 131: exact arithmetic, typed prime sets, and FRACTRAN relabeling.
Python 3 standard library only. No network, git, input files, or assertions.
All accepted checks remain active under python -O. Prints deterministic JSON.
Class-number scan: primes 5 <= p < 1000. Point-count scan: p <= 199.
Level-two enumeration: thirteen Ogg primes >= 5, no global extrapolation.
"""
from fractions import Fraction
from functools import lru_cache
from math import gcd, isqrt, prod
import json

CHECKS = 0

def require(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise ValueError('CHECK FAILED: ' + label)


def primes_below(n):
    return [x for x in range(2, n) if all(x % d for d in range(2, isqrt(x) + 1))]


def valuation(n, p):
    e = 0
    while n % p == 0:
        e += 1
        n //= p
    return e


def factor(n):
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


@lru_cache(None)
def forms(D):
    if D >= 0 or D % 4 not in (0, 1):
        raise ValueError('Not a negative quadratic-order discriminant')
    out = []
    for a in range(1, isqrt((-D) // 3) + 1):
        for b in range(-a, a + 1):
            if (b*b - D) % (4*a):
                continue
            c = (b*b - D) // (4*a)
            if c < a or gcd(gcd(a, b), c) != 1:
                continue
            if b < 0 and (abs(b) == a or a == c):
                continue
            out.append((a, b, c))
    return tuple(out)


def genus_row(p):
    """Cited Deuring/class-number formulas, evaluated exactly for p >= 5."""
    h4 = len(forms(-4*p))
    hp = len(forms(-p)) if p % 4 == 3 else 0
    total = p // 12 + {1: 0, 5: 1, 7: 1, 11: 2}[p % 12]
    require((h4 + hp) % 2 == 0, f'integral rational count p={p}')
    rational = (h4 + hp) // 2
    require(total >= rational and (total-rational) % 2 == 0,
            f'nonnegative conjugate-pair count p={p}')
    gplus = (total-rational) // 2
    # Direct X0(p) genus formula, independent of the floor expression.
    legendre = lambda a: 0 if a % p == 0 else (1 if pow(a % p, (p-1)//2, p) == 1 else -1)
    g0 = Fraction(p+1, 12) - Fraction(1+legendre(-1), 4) - Fraction(1+legendre(-3), 3)
    require(g0.denominator == 1 and g0 == total-1, f'X0 genus p={p}')
    require(Fraction(2*g0+2-h4-hp, 4) == gplus, f'Fricke genus p={p}')
    return dict(p=p, h_minus_4p=h4, h_minus_p=hp if hp else None,
                total_ss_j=total, rational_ss_j=rational,
                frobenius_plus_dim=(total+rational)//2,
                frobenius_minus_dim=gplus, integral_eigensum_index=2**gplus, genus_X0=int(g0), genus_X0_plus=gplus)


def rational_ss_j(p):
    """Enumerate j in Fp; each model is nonsingular with that j; trace zero iff ss."""
    chi = [-1] * p
    chi[0] = 0
    for x in range(1, p):
        chi[x*x % p] = 1
    cubic = [(x*x*x) % p for x in range(p)]
    out = []
    for j in range(p):
        if j == 0:
            a, b = 0, 1
        elif j == 1728 % p:
            a, b = 1, 0
        else:
            k = j * pow((1728-j) % p, -1, p) % p
            a, b = 3*k % p, 2*k % p
        disc = (4*a*a*a + 27*b*b) % p
        require(disc != 0, f'nonsingular p={p},j={j}')
        require(1728 * 4*a*a*a * pow(disc, -1, p) % p == j,
                f'model j invariant p={p},j={j}')
        trace = -sum(chi[(cubic[x]+a*x+b) % p] for x in range(p))
        if trace == 0:
            out.append(j)
    return out


def level_two(p, ss):
    """Direct Fp2 u-enumeration; s=u+4096/u is the Fricke quotient coordinate."""
    d = next(x for x in range(2, p) if pow(x, (p-1)//2, p) == p-1)
    def add(x, y):
        return ((x[0]+y[0]) % p, (x[1]+y[1]) % p)
    def mul(x, y):
        return ((x[0]*y[0]+d*x[1]*y[1]) % p,
                (x[0]*y[1]+x[1]*y[0]) % p)
    def inv(x):
        z = pow((x[0]*x[0]-d*x[1]*x[1]) % p, -1, p)
        return (x[0]*z % p, -x[1]*z % p)
    def jmap(u):
        v = add(u, (256 % p, 0))
        return mul(mul(mul(v, v), v), mul(inv(u), inv(u)))
    supers = set(ss)
    svalues = set()
    edges = {}
    for a in range(p):
        for b in range(p):
            u = (a, b)
            if u == (0, 0):
                continue
            ui = inv(u)
            require(mul(u, ui) == (1, 0), f'Fp2 inversion p={p},u={u}')
            j = jmap(u)
            if j[1] or j[0] not in supers:
                continue
            dual = mul((4096 % p, 0), ui)
            s = add(u, dual)
            svalues.add(s)
            jdual = jmap(dual)
            require(jdual[1] == 0 and jdual[0] in supers, f'dual ss endpoint p={p}')
            edges.setdefault(s, set()).add(tuple(sorted((j[0], jdual[0]))))
    rational = sorted(a for a, b in svalues if b == 0)
    nonrational = sorted((a, b) for a, b in svalues if b != 0)
    quadratics = set()
    endpoints = set()
    for a, b in svalues:
        require((a, -b % p) in svalues, f'Frobenius closure p={p}')
        if b:
            quadratics.add((1, -2*a % p, (a*a-d*b*b) % p))
            endpoints.update(edges[(a, b)])
    return dict(p=p, field_nonresidue=d, total_points=len(svalues), rational_points=rational,
                nonrational_points=nonrational, minus_dimension=len(nonrational)//2, integral_eigensum_index=2**(len(nonrational)//2),
                irreducible_quadratics=sorted(quadratics), endpoints=sorted(endpoints))


def relabel_check():
    program = [(17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),
               (77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1)]
    support = sorted(set().union(*(set(factor(a)) | set(factor(b)) for a,b in program)))
    require(support == primes_below(30), 'printed program first-ten support')
    def rename(n):
        e = valuation(n, 29)
        return (n // (29**e)) * (67**e)
    relabeled = [(rename(a),rename(b)) for a,b in program]
    require(relabeled[4] == (67,33) and relabeled[5] == (77,67), '29-to-67 fractions')
    n = m = 2
    seen_changed_register = False
    output = []
    trace = []
    for step in range(1, 10001):
        original_applicable = [i for i,(a,b) in enumerate(program) if n*a % b == 0]
        relabeled_applicable = [i for i,(a,b) in enumerate(relabeled) if m*a % b == 0]
        require(original_applicable == relabeled_applicable, f'all guards at step={step}')
        i = original_applicable[0]
        a,b = program[i]
        aa,bb = relabeled[i]
        n = n*a//b
        m = m*aa//bb
        require(m == rename(n), f'conjugacy at step={step}')
        seen_changed_register |= n % 29 == 0
        trace.append(i+1)
        if n & (n-1) == 0:
            require(m == n, f'unchanged output at step={step}')
            output.append([step,n.bit_length()-1])
    require(seen_changed_register, 'relabel witness actually reaches changed register')
    # Structural control graph uses one-hot flag primes; state-free loops removed.
    control = {11,13,17,19,23,29}
    colors = {11:0,17:0,23:0,13:1,19:1,29:1,None:1}
    arcs=[]
    for a,b in program:
        ins = set(factor(b)) & control
        outs = set(factor(a)) & control
        require(len(ins) <= 1 and len(outs) <= 1, 'one input/output flag per fraction')
        src = next(iter(ins),None)
        dst = next(iter(outs),None)
        if src is None and dst is None:
            continue
        require(colors[src] != colors[dst], 'control graph bipartition')
        arcs.append([src,dst])
    return dict(relabel={29:67}, steps=10000, changed_fractions=[[67,33],[77,67]],
                pure_two_outputs=output, control_arcs=arcs,
                relabeled_non_Ogg_prime=67, original_residue_mod8=29%8,
                relabeled_residue_mod8=67%8)


def main():
    A={2,3,5,7,11,13}; Btail={17,19,23,31,47}; C={29,41,59,71}; D={37,43,53,61,71}
    mexp={2:46,3:20,5:9,7:6,11:2,13:3,17:1,19:1,23:1,29:1,31:1,41:1,47:1,59:1,71:1}
    bexp={2:41,3:13,5:6,7:2,11:1,13:1,17:1,19:1,23:1,31:1,47:1}
    fexp={2:17,3:9,5:2,7:1,11:1,13:1}
    M=prod(p**e for p,e in mexp.items()); B=prod(p**e for p,e in bexp.items())
    Fi22=prod(p**e for p,e in fexp.items())
    require(M==808017424794512875886459904961710757005754368000000000, 'Monster decimal order')
    require(B==4154781481226426191177580544000000, 'Baby decimal order')
    require(Fi22==64561751654400, 'Fi22 decimal order')
    H=2*B; require(M % H==0, 'centralizer divides group order'); I=M//H
    require(I==97239461142009186000, '2A index decimal')
    iexp={p:mexp[p]-bexp.get(p,0)-(1 if p==2 else 0) for p in mexp}
    iexp={p:e for p,e in iexp.items() if e}
    require(prod(p**e for p,e in iexp.items())==I, '2A index factorization')
    require(set(bexp)&set(iexp)==A, 'support intersection = A')
    require(set(bexp)-set(iexp)==Btail, 'support stabilizer only = Btail')
    require(set(iexp)-set(bexp)==C, 'support index only = C')
    require({p for p,e in mexp.items() if e>=2}==A, 'Monster repeated-prime support')
    require({p for p,e in bexp.items() if e>=2}=={2,3,5,7}, 'Baby repeated-prime support')
    require(set(fexp)==A, 'Fi22 support')
    require(set(mexp)==A|Btail|C and not (A&Btail or A&C or Btail&C), 'three disjoint shells')
    omitted=set(primes_below(72))-set(mexp)
    require(omitted=={37,43,53,61,67}, 'actual bounded complement')
    require(D&set(mexp)=={71} and omitted-D=={67}, 'supplied fourth set incidence')
    for disc,expected in [(-3,1),(-4,1),(-7,1),(-8,1),(-11,1),(-12,1),(-15,2),(-20,2),(-23,3),(-163,1)]:
        require(len(forms(disc))==expected, f'class number control D={disc}')
    rows=[genus_row(p) for p in primes_below(1000) if p>=5]
    require({r['p'] for r in rows if r['genus_X0_plus']==0}==set(mexp)-{2,3}, 'bounded Ogg classification below 1000')
    pointsets={}
    for row in rows:
        p=row['p']
        if p<=199:
            pointsets[p]=rational_ss_j(p)
            require(len(pointsets[p])==row['rational_ss_j'],f'direct point-count comparison p={p}')
    levels=[level_two(p,pointsets[p]) for p in sorted(mexp) if p>=5]
    expected={29:([1,20,16],[[2,25]]),41:([1,8,33],[[3,32]]),
              59:([1,49,53],[[28,48]]),71:([1,49,53],[[41,66]])}
    for row in levels:
        p=row['p']
        require(row['minus_dimension']==int(p in C), f'level-two odd sector p={p}')
        if p in C:
            quadratic,endpoint=expected[p]
            require(row['irreducible_quadratics']==[tuple(quadratic)], f'known irreducible quadratic p={p}')
            require(row['endpoints']==[tuple(x) for x in endpoint], f'parallel-edge endpoints p={p}')
    relabel=relabel_check()
    genus_D={r['p']:r['genus_X0_plus'] for r in rows if r['p'] in D|{67}}
    require(genus_D=={37:1,43:1,53:1,61:1,67:2,71:0}, 'fourth-set genus fingerprint')
    result=dict(status='PASS',check_count=CHECKS,
                scope=dict(genus_primes='5 <= p < 1000',point_count_primes='5 <= p <= 199',
                           level_two_primes=sorted(set(mexp)-{2,3}),small_primes_2_3='cited, not point-counted'),
                group_orders=dict(M=M,B=B,Fi22=Fi22,centralizer_double_cover_B=H,index_2A=I,index_factorization=iexp),
                sets=dict(A=sorted(A),Btail=sorted(Btail),C=sorted(C),supplied_D=sorted(D),
                          complement_below_72=sorted(omitted),D_genus=genus_D),
                genus_rows=rows,direct_pointcount_ss_j=pointsets,level_two=levels,primegame=relabel)
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
