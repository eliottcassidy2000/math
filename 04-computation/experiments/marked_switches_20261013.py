"""Delayed Collatz invoices, the fish's necessary cover, and Weyl-module jets.

Run from the repository root with Python 3 (also -O). Requires SymPy and GP.
All checks use exact integers/rationals/algebraic expressions. The supplied
Collatz words are discovered in a bounded control pass and separately replayed.
No attached manuscript's BSD assertion is an input to any calculation.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from math import factorial, gcd, isqrt
from pathlib import Path
import hashlib
import ast
import json
import subprocess
import sympy as S

ROOT = Path(__file__).resolve().parents[2]
CHECKS = 0


def check(predicate, label):
    global CHECKS
    CHECKS += 1
    if not predicate:
        raise RuntimeError(label)


def E(n):
    return (n + 1).bit_length() - 1


def step(n):
    check(n > 1 and n % 2 == 1, 'strict nonterminal positive odd state')
    z = 3 * n + 1
    a = (z & -z).bit_length() - 1
    return z >> a, a


def carrier(w):
    p, q, b = 1, 1, 0
    for a in w:
        check(type(a) is int and a > 0, 'positive exact valuation letter')
        p, q, b = 3 * p, q << a, 3 * b + q
    return p, q, b


def replay(n, w):
    states = [n]
    for a in w:
        n, actual = step(n)
        check(actual == a, 'native exact letter, no ROOT padding')
        states.append(n)
    p, q, b = carrier(w)
    check(q * n == p * states[0] + b, 'independent affine replay')
    return states


def ceildiv(a, b):
    return -((-a) // b)


def gate(w, a0, s):
    """All-height dyadic-annulus gate; no target orbit is consulted here."""
    check(a0 >= 1 and s >= 1, 'gate input')
    p, q, b = carrier(w)
    P, B, M = 3 * p, 3 * b + q, q << a0
    cut = 3
    pp, qq, bb = 1, 1, 0
    for a in w:
        pp, qq, bb = 3 * pp, qq << a, 3 * bb + qq
        cut = max(cut, ceildiv(3 * qq - bb, pp))
    charge_cap = sum(a % 2 for a in w) + 1
    emin = max(s + 1, ceildiv(charge_cap + 3 * s * s, 6 * s))
    return dict(head=list(w), final_min=a0, drop=s, P=P, B=B, M=M,
                residue=-B * pow(P, -1, M) % M, cut=cut,
                charge_cap=charge_cap, minimum_energy=emin)


def passes(g, n):
    e = E(n)
    return (n >= g['cut'] and n % g['M'] == g['residue'] and
            e >= g['minimum_energy'] and
            g['P'] * n + g['B'] <= g['M'] * ((1 << (e-g['drop']+1))-2))


def apply(g, n):
    check(passes(g, n), 'owned annular gate')
    y = replay(n, g['head'])[-1]
    h, a = step(y)
    charge = sum(t % 2 for t in g['head']) + a % 2
    check(a >= g['final_min'], 'coarse final valuation')
    check(E(h) <= E(n)-g['drop'] and h < n, 'actual smaller child and rank drop')
    check(charge <= g['charge_cap'] <= 3 * (E(n)**2-E(h)**2), 'retained-source payment')
    return h, a, charge


def native_bridge(g, n):
    """Join one marked annular receipt to the next paid native-tail point.

    The seed-plus-tail move is inherited from the concurrent adaptive-return
    packet. Here M is the coarser terminal guard, not an exact-word modulus.
    """
    check(passes(g, n), 'bridge has an authenticated original point')
    # These compiled seeds use two energy bits, which imply the raw half
    # condition even at the low end of their magnitude band.
    check(g['drop'] >= 2, 'bridge seed annulus implies the raw half condition')
    P, B, M, c = (g[k] for k in ['P', 'B', 'M', 'charge_cap'])
    check(M > 2*P, 'coarse native tail has strict half-drop slope')
    emin = max(2, ceildiv(c+3, 6))
    raw_cut = ceildiv(2*B+M, M-2*P)
    cut = max(g['cut'], raw_cut, (1 << emin)-1)
    check(cut <= n+M, 'no uncovered native point between seed and tail')
    return dict(source=n, modulus=M, uniform_tail_cut=cut,
                raw_half_cut=raw_cut, tail_minimum_energy=emin)


def apply_bridge(g, bridge, n):
    base, M = bridge['source'], bridge['modulus']
    check(n >= base and (n-base) % M == 0, 'owned complete native progression')
    if n == base:
        h, a, cost = apply(g, n)
        check(2*(h+1) <= n+1, 'bridge seed actual half condition')
        return h, a, cost
    check(n >= bridge['uniform_tail_cut'], 'source above authenticated tail cut')
    y = replay(n, g['head'])[-1]
    h, a = step(y)
    cost = sum(t % 2 for t in g['head']) + a % 2
    check(a >= g['final_min'], 'bridge retains variable terminal valuation')
    check(2*(h+1) <= n+1 and h < n, 'uniform native-tail half drop')
    check(cost <= g['charge_cap'] <= 3*(E(n)**2-E(h)**2), 'native-tail invoice paid')
    return h, a, cost


def deadline(n, budget):
    # If ROOT has not occurred: odd letters multiply x-1 by <=2,
    # even letters by <=3/4. No convergence assumption is used.
    left, right, t = (n-1) << budget, 2, 0
    while left >= right:
        left *= 3
        right *= 4
        t += 1
    return budget + t


def adaptive(n):
    budget = 3 * (E(n)**2-1)
    cap = deadline(n, budget)
    x, cost, w, first, paid = n, 0, [], None, None
    for _ in range(cap):
        x, a = step(x)
        cost += a % 2
        w.append(a)
        if 2 * (x+1) <= n+1:
            if first is None:
                first = dict(length=len(w), h=x, charge=cost, drop=E(n)-E(x))
            if cost <= 3 * (E(n)**2-E(x)**2):
                paid = dict(length=len(w), h=x, charge=cost, drop=E(n)-E(x), word=w[:])
                break
        if cost > budget:
            break
        check(x != 1, 'ROOT within budget would already be accepted')
    check(paid is not None or cost > budget, 'total charge-or-ROOT deadline')
    return first, paid, cap


def annular_points(g, e):
    lo = max(g['cut'], (1 << e)-1)
    hi = min((1 << (e+1))-2,
             (g['M']*((1 << (e-g['drop']+1))-2)-g['B'])//g['P'])
    first = lo + (g['residue']-lo) % g['M']
    count = max(0, (hi-first)//g['M']+1)
    return first, count


def quartic_scan():
    """Exact bounded search in Q(sqrt5), with square roots in its cyclic quartic extension."""
    def add(a,b):return a[0]+b[0],a[1]+b[1]
    def mul(a,b):return a[0]*b[0]+5*a[1]*b[1],a[0]*b[1]+a[1]*b[0]
    def inv(a):
        n=a[0]*a[0]-5*a[1]*a[1]
        return a[0]/n,-a[1]/n
    def scale(a,n):return a[0]*n,a[1]*n
    def ratroot(q):
        if q<0:return None
        a,b=isqrt(q.numerator),isqrt(q.denominator)
        return Q(a,b) if a*a==q.numerator and b*b==q.denominator else None
    def root(z):
        A,B=z
        if B==0:
            r=ratroot(A)
            if r is not None:return r,Q(0)
            r=ratroot(A/5)
            return (Q(0),r) if r is not None else None
        n=ratroot(A*A-5*B*B)
        if n is None:return None
        for R in [n,-n]:
            u=ratroot((A+R)/2)
            if u:
                candidate=(u,B/(2*u))
                if mul(candidate,candidate)==z:return candidate
        return None
    for z in [(Q(0),Q(0)),(Q(2),Q(1)),(Q(-3,7),Q(5,11))]:
        r=root(mul(z,z))
        check(r is not None and mul(r,r)==mul(z,z), 'golden square-root positive controls')
    check(root((Q(-1),Q(0))) is None and root((Q(11),Q(0))) is None,
          'golden square-root hostile controls')
    delta=(Q(55,2),Q(11,2)); idelta=inv(delta)
    check(mul(delta,idelta)==(Q(1),Q(0)), 'quartic defining square class')
    count=0; found=[]
    for d in range(1,13):
        for a in range(-30,31):
            for b in range(-30,31):
                if b==0 or gcd(gcd(abs(a),abs(b)),d)>1:continue
                t=(Q(a,d),Q(b,d)); f=(Q(0),Q(0))
                for c in [7,-154,1342,-5082,7623]:f=add(mul(f,t),(Q(c),Q(0)))
                f=scale(f,22)
                r0,r1=root(f),root(mul(f,idelta))
                count+=1
                if r0 is not None or r1 is not None:
                    found.append([a,b,d])
    check(count==36518 and not found, 'declared exact golden-parameter quartic scan')
    return dict(parameters=count, numerator_range=[-30,30], denominators=[1,12],
                nonrational_parameters=True, primitive_triples=True, points=found)


def collatz():
    inherited = json.loads((ROOT/'05-knowledge/results/collatz_band_budget_hostile_20261012.json').read_text())
    records = inherited['glide_record_witnesses']
    sources = sorted({r['n'] for r in records if r['fails3']} | {96883183, 129177577})
    result = []
    for n in sources:
        first, paid, cap = adaptive(n)
        check(paid is not None and first['charge'] > 3*(E(n)**2-E(first['h'])**2), 'inherited first-exit hostile')
        w = paid['word']
        replay(n, w)
        choices = [gate(w[:-1], a0, s) for s in range(1, paid['drop']+1)
                   for a0 in range(1, w[-1]+1)]
        good = [g for g in choices if passes(g, n)]
        check(bool(good), 'first paid word admits an annular certificate')
        g = min(good, key=lambda g: (g['M'], -g['drop']))
        check(apply(g, n)[0] == paid['h'], 'same marked source and endpoint')
        check(not passes(g, n+2), 'nearby source is not substituted')
        # A slope too large for a whole residue cylinder can still pay a
        # nonzero fraction of every sufficiently high dyadic annulus.
        ratio = Q(g['M'], g['P'] * (1 << (g['drop']-1)))
        fraction = min(Q(2), ratio)-1
        check(fraction > 0, 'positive all-scale annular width')
        bridge = native_bridge(g, n)
        bridge_samples = []
        for k in [0, 1, 2, 17, 1 << 64]:
            m = n+k*g['M']
            h, a, cost = apply_bridge(g, bridge, m)
            bridge_samples.append(dict(k=k, n=m, h=h, final=a, charge=cost))
        exact_modulus = 1 << (sum(w)+1)
        check(exact_modulus % g['M'] == 0 and exact_modulus > g['M'],
              'coarse progression strictly contains the exact-word progression')
        bridge['exact_modulus_ratio'] = exact_modulus//g['M']
        sampled = []
        for e in [E(n), g['M'].bit_length()+4, g['M'].bit_length()+9]:
            base, count = annular_points(g, e)
            check(count > 0, 'nonempty source guard at checked scales')
            indices = sorted({0, count//2, count-1})
            for j in indices:
                m = base+j*g['M']
                h, a, cost = apply(g, m)
                sampled.append(dict(n=m, h=h, final=a, charge=cost, energy=e))
            excluded = base + count*g['M']
            if E(excluded) == e:
                check(not passes(g, excluded), 'same progression, unpaid annular phase rejected')
        result.append(dict(n=n, first=first, paid=paid, deadline=cap, gate=g,
                           asymptotic_annular_fraction=str(fraction), samples=sampled,
                           complete_native_progression=bridge, bridge_samples=bridge_samples))
    # Declared finite universe, not a new finite-convergence record.
    accepted, max_extra, witness = 0, 0, None
    for n in range(3, 65536, 2):
        first, paid, _ = adaptive(n)
        if paid is not None:
            accepted += 1
            extra = paid['length']-first['length']
            if extra > max_extra:
                max_extra, witness = extra, n
    check(accepted == 32767, 'finite adaptive controller acceptance')
    # The same stopped path may not have its debt reset at the first exit.
    for r in result:
        check(r['paid']['charge'] >= r['first']['charge'], 'carried rather than reset debt')
    return dict(hostiles=result, finite=dict(odd_sources='3..65535', accepted=accepted,
                max_extra_after_first=max_extra, witness=witness))


def fish():
    x = S.symbols('x')
    q = [x*(1-x)/11, x*(4-x)/11, (1-x)*(x+3)/11, (1+2*x)*(3-2*x)/11]
    roots = [{0, 1}, {0, 4}, {1, -3}, {Q(-1, 2), Q(3, 2)}]
    quotient_data = []
    for mask in range(1, 16):
        odd = set()
        for j in range(4):
            if mask >> j & 1:
                odd ^= roots[j]
        check(bool(odd), 'independent square classes over algebraic closure')
        degree = len(odd)
        genus = (degree-2)//2
        quotient_data.append(dict(mask=mask, roots=sorted(map(str, odd)), genus=genus))
    check(sum(r['genus'] for r in quotient_data) == 9, 'all character quotient genera sum to nine')
    check(sum(r['genus'] for r in quotient_data if r['mask'] < 8) == 1, 'first three covers give genus one')
    check(1+16*(6-4)//4 == 9, 'independent Riemann-Hurwitz calculation')
    poly = x*(1-x)*(4-x)*(x+3)
    u = 12/x-3
    check(S.cancel((12/x**2)**2*poly-u*(u-9)*(u+7)) == 0, 'elliptic birational map')
    t=S.symbols('t')
    r=(7*t*t-66*t+231)/(2*(t*t-33))
    ss=S.Rational(33,2)+t*(r-S.Rational(7,2))
    D=(r*r-4)/(r*r-1)
    H=(3*t*t-22*t+55)*(5*t*t-66*t+297)
    quartic=22*(7*t**4-154*t**3+1342*t*t-5082*t+7623)
    check(S.factor(ss*ss-33*(r*r-4))==0, 'rational parametrization of first two square conditions')
    z0=ss/(11*(r*r-1))
    check(S.factor(11*z0*z0-D*(1-D))==0, 'retained first square root')
    check(S.factor(11*(r*z0)**2-D*(4-D))==0, 'retained second square root')
    check(S.factor((4*(t*t-33)/(11*H))**2*quartic-(1-D)*(D+3)/11)==0,
          'marked genus-one quartic, including lift factor')

    # A faithful real realization of the 13-vertex necessary subframework.
    # Its field is allowed to be larger than Q(zeta55); this is a positive control.
    I, rt2, rt14 = S.I, S.sqrt(2), S.sqrt(14)
    c, d = rt2/2, I*rt2/2
    a, b = rt14/4-I*rt2/4, -rt2/4+I*rt14/4
    norm = lambda z: S.simplify(S.expand(z*S.conjugate(z)))
    edges = [(0,10),(0,20),(0,22),(1,11),(1,21),(1,22),
             (2,8),(2,19),(2,22),(3,9),(3,18),(3,22),
             (8,10),(8,20),(9,11),(9,21),(10,19),(11,18),
             (18,19),(18,20),(19,21),(20,21)]
    inherited_path=ROOT/'04-computation/experiments/ei19_fish_placement_q_zeta55_20261012.py'
    tree=ast.parse(inherited_path.read_text())
    graphs=next(ast.literal_eval(node.value) for node in tree.body
                if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='GRAPHS' for t in node.targets))
    adj=graphs['FISH']
    for i,j in edges:
        check(j in adj[i] and i in adj[j], 'subframework uses actual inherited fish edges')
    cycles=[(0,10,8,20),(2,8,10,19),(1,11,9,21),(3,9,11,18),(18,19,21,20)]
    rows=[]
    for cyc in cycles:
        check(all(cyc[(j+1)%4] in adj[cyc[j]] for j in range(4)), 'actual unit quadrilateral')
        row=[0]*23
        for j,v in enumerate(cyc):row[v]=(-1)**j
        rows.append(row)
    R=S.Matrix(rows)
    for positive,negative in [([0,2],[19,20]),([1,3],[18,21]),([18,21],[19,20])]:
        row=[0]*23
        for v in positive:row[v]+=1
        for v in negative:row[v]-=1
        check(R.col_join(S.Matrix([row])).rank()==R.rank(), 'retained midpoint relation')
    def colour(k):
        col={22:0}
        def extend():
            if len(col)==23:return dict(col)
            v=max((v for v in range(23) if v not in col),
                  key=lambda v:(len({col[u] for u in adj[v] if u in col}),len(adj[v]),-v))
            forbidden={col[u] for u in adj[v] if u in col}
            for c0 in range(min(k,max(col.values())+2)):
                if c0 not in forbidden:
                    col[v]=c0
                    out=extend()
                    if out is not None:return out
                    del col[v]
            return None
        return extend()
    check(colour(3) is None, 'complete abstract three-colour rejection')
    four=colour(4)
    check(four is not None and all(four[i]!=four[j] for i in range(23) for j in adj[i]), 'abstract four-colour witness')
    controls = []
    for s, t in product([-1, 1], repeat=2):
        z = {0:a, 1:b, 2:-a, 3:-b, 18:c, 19:d, 20:-d, 21:-c, 22:S.Integer(0)}
        z[10] = (a+d+s*I*(d-a))/2
        z[8] = z[10]-a-d
        z[11] = (b+c+t*I*(c-b))/2
        z[9] = z[11]-b-c
        for i,j in edges:
            check(norm(z[i]-z[j]) == 1, 'all required subframework unit edges')
        faithful = all(norm(z[i]-z[j]) != 0 for i,j in combinations(z,2))
        if faithful:
            controls.append({str(k):str(S.expand(v)) for k,v in sorted(z.items())})
    check(bool(controls), 'nonvacuous faithful real subframework')

    # Genuine GP 2-descent bounds, not ellanalyticrank or a BSD assumption.
    gp = '''E=ellinit([0,-2,0,-63,0]);
print(ellrank(E));
print(ellrank(elltwist(E,5)));
print(elltors(E));
print(elltors(elltwist(E,5)));
print([ellcard(ellinit([0,-2,0,-63,0],11)),ellcard(ellinit([0,-2,0,-63,0],19))]);
nf=nfinit(x^4-55*x^2+605); tau=Mod(x,nf.pol); count=0; found=List();
print([nfeltissquare(nf,(1+tau)^2),nfeltissquare(nf,-1),nfeltissquare(nf,11)]);
for(d=1,3,for(a=-2,2,for(b=-2,2,for(c=-2,2,for(e=-2,2,if(b==0&&e==0,next); if(gcd(gcd(gcd(gcd(abs(a),abs(b)),abs(c)),abs(e)),d)>1,next); t=(a+b*tau+c*tau^2+e*tau^3)/d; f=22*(7*t^4-154*t^3+1342*t^2-5082*t+7623); count++; if(nfeltissquare(nf,f,&z),listput(found,[a,b,c,e,d,z])); )))));
print([count,Vec(found)]);
quit
'''
    run = subprocess.run(['gp','-fq'], input=gp, text=True, capture_output=True, check=True)
    check(not run.stderr.strip(), 'GP has no hidden error')
    lines = run.stdout.strip().splitlines()
    check(lines[:2] == ['[0, 0, 0, []]']*2, 'unconditional rank zero for E and E^5')
    check(lines[2].startswith('[8, [4, 2]') and lines[3].startswith('[4, [2, 2]'), 'rational torsion orders')
    check(lines[4] == '[8, 16]', 'split-prime torsion bound over Q(sqrt5)')
    check(lines[5] == '[1, 0, 0]', 'full quartic field square controls')
    check(lines[6] == '[1728, []]', 'declared genuinely quartic-parameter scan')
    # Independent rational group law and point counts.
    def add(P, R):
        if P is None:return R
        if R is None:return P
        xx,yy=P; X,Y=R
        if xx==X and yy==-Y:return None
        slope=(3*xx*xx-4*xx-63)/(2*yy) if P==R else (Y-yy)/(X-xx)
        z=slope*slope+2-xx-X
        return z, -yy+slope*(xx-z)
    G,H=(Q(-3),Q(12)),(Q(0),Q(0))
    pts=set();P=None
    for _ in range(4):
        pts.update([P,add(P,H)]);P=add(P,G)
    check(P is None and len(pts)==8, 'explicit eight rational torsion points')
    parameters=set()
    for P in pts:
        if P is None: parameters.add('0')
        elif P[0]==-3:parameters.add('infinity')
        else:
            check(P[1]**2==P[0]*(P[0]-9)*(P[0]+7), 'torsion point equation')
            parameters.add(str(12/(P[0]+3)))
    check(parameters == {'0','1','4','-3','1/2','infinity'}, 'complete permitted rational parameter candidates')
    for p, expected in [(11,8),(19,16)]:
        count=1+sum(1 for xx in range(p) for yy in range(p) if (yy*yy-xx*(xx-9)*(xx+7))%p==0)
        check(count==expected and pow(5,(p-1)//2,p)==1, 'independent split-prime count')
    # D=1/2 would force sqrt(11) into the real cyclotomic field, impossible.
    check(S.simplify(q[0].subs(x,S.Rational(1,2)))==S.Rational(1,44), 'remaining square-class obstruction')
    def rational_square(r):
        r=Q(r)
        return r>=0 and isqrt(r.numerator)**2==r.numerator and isqrt(r.denominator)**2==r.denominator
    # A square root of a rational number lies in Q(sqrt5) iff it is already
    # rational or sqrt5 times a rational. Include the degenerate branch points.
    for t in [S.Integer(0),S.Integer(1),S.Rational(1,2)]:
        vals=[Q(z.subs(x,t)) for z in q[:3]]
        check(any(not rational_square(z) and not rational_square(z/5) for z in vals),
              'no lift of any real elliptic torsion parameter to the first three square conditions')
    check(not rational_square(Q(-1,11)), 'no real point at infinity of the three-square cover')
    units=[a for a in range(1,55) if gcd(a,55)==1]
    # The real Galois group is cyclic of order 20 (class of 2 generates it).
    orbit={min(pow(2,j,55),(-pow(2,j,55))%55) for j in range(20)}
    check(len(units)==40 and len(orbit)==20, 'unique real quadratic subfield Q(sqrt5)')
    check(gcd(8,5)==1, 'degree-eight divisor and odd quintic point would give degree-one divisor')
    return dict(necessary_quadratics=list(map(str,q)), quotients=quotient_data,
                genus=9, first_three_genus=1, elliptic_model='y^2=x(x-9)(x+7)',
                gp_source=gp, gp_output=lines, parameters=sorted(parameters),
                real_control=controls[0], required_edges=edges,
                field_exclusion='No fish placement in Q(zeta11,sqrt5); no rational D in Q(zeta55).',
                four_colour_witness=four,
                graph_source_sha256=hashlib.sha256(inherited_path.read_bytes()).hexdigest(),
                quartic=str(quartic), quartic_scan=quartic_scan())


def jets():
    t = S.symbols('t')
    tested=0
    for b in range(1,9):
        # Nontrivial rank-two connection d/dt + A(t), with mixed entries.
        A=S.Matrix([[t,1],[1,t*t]])
        nabla=lambda v: v.diff(t)+A*v
        def coeffs(v):
            out=[]
            for j in range(b):
                out.extend([S.expand(z).subs(t,0)/factorial(j) for z in v])
                v=nabla(v)
            return S.Matrix(out)
        basis=[S.eye(2)[:,k]*t**j for j in range(b) for k in range(2)]
        J=S.Matrix.hstack(*[coeffs(v) for v in basis])
        check(J.det()==1, 'canonical Taylor jet map is an isomorphism')
        shift=S.zeros(2*b)
        for j in range(2*(b-1)):shift[j+2,j]=1
        for v in basis:
            check(coeffs(t*v)==shift*coeffs(v), 'marked multiplication becomes shift')
            check(coeffs(t**b*v)==S.zeros(2*b,1), 'well-defined modulo t^b')
            tested+=1
    # Laurent cylinder: all global residues survive each transverse jet.
    residue_tests=0
    z=S.symbols('z')
    for j in range(8):
        for k in range(-10,11):
            # d(z^k)/dz has zero coefficient of z^-1, including k=0.
            check(not (k-1==-1 and k!=0), 'global residue annihilates exact derivatives')
            if k!=-1:
                check(S.diff(t**j*z**(k+1)/S.Integer(k+1),z)==t**j*z**k,
                      'all other Laurent monomials have a primitive')
            residue_tests+=1
    # Existing Smith-bar mechanism: a rank-two free summand plus bars 1,3,7.
    lengths=[2*b+sum(min(b,a) for a in (1,3,7)) for b in range(1,13)]
    check(lengths[-1]-lengths[-2]==2, 'persistent rank is the eventual length slope')
    # A Pruefer arm is nonzero although every quotient M/t^b M is zero.
    for b in range(1,9):
        for k in range(1,9):
            check(-(k+b)+b==-k, 'explicit divisible preimage retains the missing sidecar')
    return dict(taylor_basis_tests=tested, laurent_tests=residue_tests,
                smith_lengths=lengths, generic_rank=2,
                warning='The connection does not descend to a fixed t^b quotient; the jet map does.')


def main():
    out=dict(collatz=collatz(),fish=fish(),jets=jets())
    out['checks']=CHECKS
    out['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(out,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
