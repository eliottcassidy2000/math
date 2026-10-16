#!/usr/bin/env python3
"""Exact controls for marked recursion, Sturm certificates, and geometric events.

Run from the repository root; stdout is a deterministic JSON receipt.
Only the explicit finite universes below are enumerated. No orbit conjecture
or attached paper's formalization is assumed by these computations.
"""
import json
import os
from math import gcd, isqrt
from fractions import Fraction as F
from itertools import combinations, product
# Some bundled SymPy/Flint combinations try to order nmod coefficients
# during factor sorting. Use SymPy's exact Python ground domain consistently.
os.environ['SYMPY_GROUND_TYPES'] = 'python'
import sympy as s

x = s.Symbol('x')


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def sturm(poly):
    a = s.Poly(poly, x, domain=s.QQ)
    b = a.diff()
    chain = [a, b]
    while not b.is_zero:
        c = -a.rem(b)
        if c.is_zero:
            break
        chain.append(c)
        a, b = b, c
    return chain


def signs(chain, at):
    result = []
    for p in chain:
        if at == '+inf':
            y = p.LC()
        elif at == '-inf':
            y = p.LC() * (-1)**p.degree()
        else:
            y = p.eval(at)
        result.append(int(s.sign(y)))
    return result


def variation(values):
    v = [a for a in values if a]
    return sum(a != b for a, b in zip(v, v[1:]))


def roots(poly, lo, hi):
    chain = sturm(poly)
    return variation(signs(chain, lo)) - variation(signs(chain, hi))


def inertia(matrix):
    """Exact congruence elimination, including zero-diagonal 2x2 pivots."""
    a = [[F(v) for v in row] for row in matrix]
    pos = neg = zero = 0
    while a:
        n = len(a)
        k = next((i for i in range(n) if a[i][i]), None)
        if k is not None:
            order = [k] + [i for i in range(n) if i != k]
            a = [[a[i][j] for j in order] for i in order]
            d = a[0][0]
            pos += d > 0
            neg += d < 0
            a = [[a[i][j] - a[i][0]*a[0][j]/d
                  for j in range(1, n)] for i in range(1, n)]
        else:
            pair = next(((i, j) for i in range(n) for j in range(i+1, n)
                         if a[i][j]), None)
            if pair is None:
                zero += n
                break
            i, j = pair
            order = [i, j] + [k for k in range(n) if k not in pair]
            a = [[a[i][j] for j in order] for i in order]
            b = a[0][1]
            pos += 1
            neg += 1
            a = [[a[i][j] - (a[i][0]*a[1][j]+a[i][1]*a[0][j])/b
                  for j in range(2, n)] for i in range(2, n)]
    return [pos, neg, zero]


def octic_control():
    f = s.Poly(x**8-20*x**7+178*x**6-842*x**5+1923*x**4-496*x**3
               -6754*x**2+12420*x-6865, x)
    chain = sturm(f)
    records = {a: signs(chain, a) for a in ['-inf', '0', '+inf']}
    # Newton sums computed by companion multiplication, separate from Sturm.
    coeff = list(reversed(f.all_coeffs()[1:]))
    companion = s.zeros(8)
    for i in range(7):
        companion[i+1, i] = 1
    for i in range(8):
        companion[i, 7] = -coeff[i]
    power = s.eye(8)
    moments = []
    for k in range(16):
        moments.append(int(s.trace(power)))
        power = power * companion
    H = [[moments[i+j] for j in range(8)] for i in range(8)]
    weighted = lambda a: [[moments[i+j+1]-F(a)*moments[i+j]
                           for j in range(8)] for i in range(8)]
    lo, hi = F(3877,1000), F(1939,500)
    inertias = {'trace': inertia(H), 'times_x': inertia(weighted(0)),
                'x_minus_lo': inertia(weighted(lo)),
                'x_minus_hi': inertia(weighted(hi))}
    check(inertias == {'trace':[5,3,0], 'times_x':[4,4,0],
                      'x_minus_lo':[4,4,0], 'x_minus_hi':[3,5,0]},
          'unexpected trace signatures')
    check([variation(records[a]) for a in ['-inf','0','+inf']] == [5,4,3],
          'paper Sturm signs')
    check(roots(f, s.Rational(lo), s.Rational(hi)) == 1, 'positive root isolation')
    P = 5*x**8-10*x**7-2*x**6+14*x**5+12*x**4-6*x**3+2*x**2+2*x-1
    C = -1373*x**8+4318*x**7+15386*x**6+29302*x**5+16044*x**4+4914*x**3+1414*x**2+186*x-31
    d = 1+2*x-x*x
    check(s.cancel(d**8*f.as_expr().subs(x,(6*x+4)/d)-P*C) == 0,
          'marked parameter eliminant')
    check(roots(P,s.Rational(9,25),s.Rational(37,100)) == 1, 'tilt root')
    factors = {}
    def xp_mod(power, modulus):
        answer = s.Poly(1,x,modulus=modulus.get_modulus())
        base = s.Poly(x,x,modulus=modulus.get_modulus())
        while power:
            if power & 1:
                answer = (answer*base).rem(modulus)
            base = (base*base).rem(modulus)
            power //= 2
        return answer
    # Explicit factors from Appendix C, checked by the finite-field
    # irreducibility criterion, rather than an opaque factorization result.
    for p, linear, degree, proper in [(7,x-1,7,1),(29,s.Integer(1),8,4),
                                    (131,(x+30)*(x+55)*(x-55),5,1)]:
        fp = f.set_modulus(p)
        lp = s.Poly(linear,x,modulus=p)
        g = fp.exquo(lp)
        xx = s.Poly(x,x,modulus=p)
        check((xp_mod(p**degree,g)-xx).rem(g).is_zero,'Frobenius return')
        check(s.gcd(xp_mod(p**proper,g)-xx,g).degree()==0,'proper factors excluded')
        check(s.gcd(fp,fp.diff()).degree()==0,'squarefree reduction')
        factors[str(p)] = [[1,1]]*lp.degree()+[[degree,1]]
    # Direct exact sign boundary hostile: a polynomial need not be zero
    # in a reducible quotient to vanish at the selected real root.
    reducible = s.Poly((x*x-1)*(x*x-2), x)
    g = s.Poly(x*x-2, x)
    check(not g.rem(reducible).is_zero, 'hostile nonzero remainder')
    check(roots(s.gcd(reducible,g), s.Rational(5,4), s.Rational(3,2)) == 1,
          'selected-root gcd')
    return {'degrees':[p.degree() for p in chain], 'endpoint_signs':records,
            'inertias':inertias, 'positive_root_interval':[str(lo),str(hi)],
            'modular_factor_degrees_and_multiplicities':factors,
            'parameter_identity':True, 'reducible_embedding_hostile':True}


def hinge_control():
    """Seven actual plane points; univariate parameter, exact Q(sqrt3)."""
    r = s.sqrt(3)
    c, sn = (1-x*x)/(1+x*x), 2*x/(1+x*x)
    base = [(0,0),(r,0),(r/2,s.Rational(1,2)),(r/2,-s.Rational(1,2))]
    rotate = lambda a,b:(s.expand(c*a-sn*b),s.expand(sn*a+c*b))
    points = base+[rotate(*p) for p in base[1:]]
    lo,hi = s.Rational(3,10),s.Rational(31,100)
    event = s.Poly(11*x*x-1,x)
    persistent, isolated, all_events = [], [], []
    for i,j in combinations(range(7),2):
        dx,dy = [s.expand(points[i][k]-points[j][k]) for k in range(2)]
        dist = s.cancel(s.expand(dx*dx+dy*dy))
        check(dist != 0, 'distinct rational point functions')
        for kind,expr in [('edge',dist-1),('collision',dist)]:
            num = s.cancel(expr).as_numer_denom()[0].expand()
            if num == 0:
                check(kind == 'edge', 'coincident points')
                persistent.append([i,j])
                continue
            norm = s.Poly(s.expand(num*num.xreplace({r:-r})),x,domain=s.QQ)
            squarefree = norm.sqf_part()
            check(squarefree.eval(lo)!=0 and squarefree.eval(hi)!=0,
                  'hinge interval endpoint event')
            count = roots(squarefree,lo,hi)
            if count:
                all_events.append([i,j,kind,count])
                check(kind == 'edge' and [i,j] == [1,4], 'unexpected hinge event')
                check(s.rem(squarefree,event).is_zero, 'event equation')
                check(s.simplify(expr.subs(x,1/s.sqrt(11)))==0,'actual selected-root edge')
                isolated.append([i,j])
    check(len(persistent) == 10 and isolated == [[1,4]], 'hinge edge list')
    def ncolour(edges,q):
        return sum(all(v[i]!=v[j] for i,j in edges)
                   for v in product(range(q),repeat=7))
    counts = {'generic_3':ncolour(persistent,3),
              'event_3':ncolour(persistent+isolated,3),
              'event_4':ncolour(persistent+isolated,4)}
    check(counts['generic_3'] > 0 and counts['event_3'] == 0 and counts['event_4'] > 0,
          'hinge chromatic jump')
    return {'interval':['3/10','31/100'],'persistent_edges':persistent,
            'event_edges':isolated,'event_equation':'11*x^2-1',
            'all_norm_events':all_events,'colouring_counts':counts}


def paley_control():
    H = {i*i%11 for i in range(1,11)}
    A = s.Matrix(11,11,lambda i,j:int((j-i)%11 in H))
    J,I = s.ones(11),s.eye(11)
    K = A-A.T
    check(A+A.T==J-I and A*A.T==3*I+2*J, 'Paley incidence')
    check(K*K==J-11*I, 'Paley oriented complex structure')
    pairs=list(combinations(range(11),2))
    orbit={tuple(sorted((b,(a+b)%11))) for a in H for b in range(11)}
    check(orbit==set(pairs), 'regular pair action')
    # It is the actual difference predicate, not a cosmetic orientation.
    for a,b in product(H,range(11)):
        check(all(A[i,j]==A[(a*i+b)%11,(a*j+b)%11]
                  for i in range(11) for j in range(11)), 'orientation invariance')
    return {'residues':sorted(H),'group_order':55,'pair_orbit':len(orbit),
            'AA_transpose':'3I+2J','skew_square':'J-11I',
            'adjacency_characteristic':str(A.charpoly(x).as_expr().factor()),
            'vertex_pair_orbitals':3}


def recursion_control():
    # Read complete future ROOT observations of powers of two.
    # No termination beyond this explicitly known family is used.
    horizon=32
    words=[]
    for k in range(horizon+1):
        n=2**k
        word=[]
        for clock in range(horizon+1):
            word.append(int(n==1))
            if n!=1:
                n=n//2 if n%2==0 else (3*n+1)//2
        words.append(tuple(word))
    check(len(set(words))==33, 'exact clock future languages')
    # Relational gluing: identical unary marginals do not retain closure.
    eq={(i,i) for i in range(6)}
    shift={(i,(i+1)%6) for i in range(6)}
    check({a for a,b in eq}=={a for a,b in shift}==set(range(6)), 'left marginals')
    check({b for a,b in eq}=={b for a,b in shift}==set(range(6)), 'right marginals')
    check(len(eq & shift)==0, 'gluing hostile')
    # Expanding and deleting a two-ended staircase have different ranks.
    ranks=[]
    for n in range(4,20):
        d=(n-1)*(n-2)//2
        inner=(n-3)*(n-4)//2
        check(d==inner+2*(n-3)+1, 'marked staircase decomposition')
        ranks.append([n,d,inner,2*n-5])
    check(inertia([[0,1],[1,0]])==[1,1,0], 'zero pivot hostile')
    return {'powers_of_two_horizon':horizon,'distinct_future_words':33,
            'full_six_colour_marginal_relations':[len(eq),len(shift)],
            'compatible_gluings_same':len(eq),'compatible_gluings_shifted':0,
            'staircase_rank_controls':ranks,'zero_diagonal_inertia_hostile':True}


def cyclotomic_hinge_control():
    z=s.Symbol('z')
    modulus=s.Poly(s.cyclotomic_poly(55,z),z,domain=s.QQ)
    def p(expr):
        return s.Poly(expr,z,domain=s.QQ).rem(modulus)
    def mul(a,b):
        return (a*b).rem(modulus)
    def bar(a):
        return p(sum(co*z**((-mon[0])%55) for mon,co in a.terms()))
    one=p(1)
    residues={i*i%11 for i in range(1,11)}
    g11=p(sum((1 if a in residues else -1)*z**(5*a) for a in range(1,11)))
    g5=p(1+2*(z**11+z**44))
    check(mul(g11,g11)==p(-11) and mul(g5,g5)==p(5),'Gauss square roots')
    g55=mul(g5,g11)
    cases=[]
    for d,g in [(11,g11),(55,g55)]:
        a=(one+g)*s.Rational(1,2)
        radius=s.Rational(d+1,4)
        gamma=(p(d-1)+2*g)*s.Rational(1,d+1)
        check(mul(a,bar(a))==p(radius),'marked radius')
        check(mul(gamma,bar(gamma))==one,'rotation norm')
        check(mul(gamma,a)==a-one,'difference-one hinge')
        points=[a]
        for _ in range(11):
            points.append(mul(gamma,points[-1]))
        edges=[]
        for i,j in combinations(range(12),2):
            delta=points[i]-points[j]
            check(not delta.is_zero,'finite orbit injectivity')
            if mul(delta,bar(delta))==one:
                edges.append([i,j])
        check(edges==[[i,i+1] for i in range(11)],'orbit is an induced path')
        cases.append({'d':d,'squared_radius':str(radius),
                      'vertices':12,'induced_edges':edges})
    a=(one+g11)*s.Rational(1,2)
    b=(p(3)+g11)*s.Rational(1,5)
    gamma55=(p(3)+g55)*s.Rational(1,8)
    check(mul(b,bar(b))==p(s.Rational(4,5)),'mixed shell radius')
    check(mul(a-b,bar(a-b))==one,'mixed shell unit cross edge')
    check(mul(gamma55,bar(gamma55))==one,'mixed shell rotation')
    delta=b-mul(gamma55,b)
    check(mul(delta,bar(delta))==one,'mixed shell unit chord')
    quartic=s.Poly(-275*x**4+150*x*x-3,x)
    check(s.gcd(quartic,quartic.diff()).degree()==0,'genus one nonsingularity')
    check(quartic.eval(s.Rational(1,5))==s.Rational(8,5)**2,'rational quartic point')
    hits=[]
    for den in range(1,101):
        for num in range(1,den):
            if gcd(num,den)!=1:
                continue
            r=F(num,den)
            radius=(1+55*r*r)/4
            disc=radius*radius-8*radius+4
            if disc>=0:
                continue
            for d in [11,55]:
                value=-disc/d
                if (isqrt(value.numerator)**2==value.numerator and
                    isqrt(value.denominator)**2==value.denominator):
                    hits.append([str(r),str(radius),d,str(value)])
    check(hits==[['1/5','4/5',11,'4/25']],'bounded mixed-shell census')
    # The quartic maps to Y^2=X^3+150X^2+825X by
    # (X,Y)=(-275*r^2,-275*r*y). Good reductions bound rational torsion.
    elliptic_counts={}
    cubic=s.Poly(x**3+150*x*x+825*x,x)
    for prime in [7,17]:
        reduced=cubic.set_modulus(prime)
        check(s.gcd(reduced,reduced.diff()).degree()==0,'elliptic good reduction')
        elliptic_counts[str(prime)]=1+sum(
            (yy*yy-(xx**3+150*xx*xx+825*xx))%prime==0
            for xx in range(prime) for yy in range(prime))
    check(elliptic_counts=={'7':6,'17':20},'torsion reduction counts')
    def add(P,Q):
        X,Y=P; U,V=Q
        slope=(3*X*X+300*X+825)/(2*Y) if P==Q else (V-Y)/(U-X)
        xx=slope*slope-150-X-U
        return xx,-Y+slope*(X-xx)
    P=(F(-11),F(-88))
    twice=add(P,P)
    triple=add(P,twice)
    r=F(91,135)
    y=triple[1]/(-275*r)
    check(triple[0]==-275*r*r and y==F(10552,3645),'tripled point lifts')
    check(y*y==-275*r**4+150*r*r-3,'new quartic point')
    return {'field':'Q[z]/Phi_55(z)', 'degree':40, 'gauss_square_roots':[-11,5,-55],
            'difference_one_controls':cases,
            'mixed_shell':{'squared_radii':['3','4/5'],'unit_cross_edge':True,
                           'quartic':'y^2=-275*r^4+150*r^2-3',
                           'rational_point':['1/5','8/5'],
                           'squarefree_quartic':True,
                           'reduced_parameter_universe':'0<r<1, denominator<=100',
                           'admissible_hits':hits,
                           'elliptic_good_reduction_counts':elliptic_counts,
                           'elliptic_point':['-11','-88'],
                           'double_point':[str(t) for t in twice],
                           'triple_point':[str(t) for t in triple],
                           'new_quartic_point':[str(r),str(y)],
                           'new_squared_radius':str((1+55*r*r)/4)}}


def mixed_patch_control():
    """Independent degree-four arithmetic: basis 1,j,k,jk, j^2=-11,k^2=-55."""
    def c(a=0,b=0,d=0,e=0):
        return tuple(map(F,(a,b,d,e)))
    def sub(a,b):
        return tuple(x-y for x,y in zip(a,b))
    def scale(a,q):
        return tuple(x*q for x in a)
    def mul(a,b):
        out=[F(0)]*4
        for i,ai in enumerate(a):
            for j,bj in enumerate(b):
                v=ai*bj
                if i&j&1:
                    v*=-11
                if i&j&2:
                    v*=-55
                out[i^j]+=v
        return tuple(out)
    def bar(a):
        return (a[0],-a[1],-a[2],a[3])
    def power(a,k):
        if k<0:
            return power(bar(a),-k)
        out=c(1)
        for _ in range(k):
            out=mul(out,a)
        return out
    alpha=c(F(1,2),F(1,2))
    gamma=c(F(5,6),F(1,6))
    points=set()
    for i in range(-3,4):
        g=power(gamma,i)
        points.add(mul(g,alpha))
        for r,y in [(F(1,5),F(8,5)),(F(91,135),F(10552,3645))]:
            radius=(1+55*r*r)/4
            b=scale(mul(c(radius+2,y/4),alpha),F(1,6))
            delta=c((55*r*r-1)/(1+55*r*r),0,2*r/(1+55*r*r))
            check(mul(delta,bar(delta))==c(1),'patch generator norm')
            for j in range(-1,2):
                points.add(mul(g,mul(power(delta,j),b)))
    points=sorted(points)
    edges=[]
    adjacency=[set() for _ in points]
    for i,j in combinations(range(len(points)),2):
        delta=sub(points[i],points[j])
        if mul(delta,bar(delta))==c(1):
            edges.append([i,j]); adjacency[i].add(j); adjacency[j].add(i)
    colours=[None]*len(points)
    colours[0]=0
    queue=[0]
    for i in queue:
        for j in sorted(adjacency[i]):
            if colours[j] is None:
                colours[j]=1-colours[i]; queue.append(j)
            check(colours[j]!=colours[i],'patch bipartite')
    check(len(points)==49 and len(edges)==48 and len(queue)==49,'patch is a tree')
    return {'universe':'i=-3..3; r=1/5,91/135; radial steps j=-1..1',
            'basis':['1','sqrt(-11)','sqrt(-55)','sqrt(-11)*sqrt(-55)'],
            'coordinates':[[str(v) for v in point] for point in points],
            'vertices':len(points),'edges':edges,'two_colouring':colours,
            'connected':True,'status':'FINITE-EXACT induced tree; not forcing'}


def main():
    data={'status':'FINITE-EXACT; explicit controls only',
          'octic':octic_control(),'hinge':hinge_control(),
          'paley':paley_control(),'recursion':recursion_control(),
          'cyclotomic_hinge':cyclotomic_hinge_control(),
          'mixed_patch':mixed_patch_control()}
    print(json.dumps(data,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
