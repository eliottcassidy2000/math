"""Exact controls for marked ROOT completions and cyclotomic contexts.

Python + sympy; no removable assertions. The infinite statements are proved
in the companion notes. No finite census here proves universal Collatz.
"""
from fractions import Fraction as Q
from itertools import product, combinations
from math import prod
import json
import sympy as s

checks = 0


def check(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)


def nu(m):
    return Q(1, (m+1)**6)-Q(1, (m+2)**6)


def replay(n, aa):
    states = [n]
    for a in aa:
        check(n > 1 and n % 2 == 1, "strict positive nonroot state")
        y = 3*n+1
        check((y & -y).bit_length()-1 == a, "literal valuation")
        n = y >> a
        states.append(n)
    check(n == 1, "literal ROOT endpoint")
    return states


def terminal_parameters(n):
    """Exact family membership; failure says nothing about convergence."""
    if not isinstance(n,int) or isinstance(n,bool) or n<=1 or n%2==0:
        return None
    d = 243*n+207
    v = (d & -d).bit_length()-1
    a = v-3
    if a<18 or a%18:
        return None
    z = ((d>>v)+5)//2
    if z<=0 or z & (z-1):
        return None
    exponent = z.bit_length()-1
    if exponent%2:
        return None
    k = exponent//2
    if k<242 or (k+1)%243:
        return None
    return a//18,(k+1)//243


def collatz():
    check([k for k in range(243) if pow(4,k,1458) == 1276] == [242],
          "complete terminal phase")
    check(s.n_order(4,729) == 243, "independent order")
    rows = []
    for j,r in product(range(1,7), range(1,7)):
        k, a = 243*r-1, 18*j
        t, rem = divmod(4**k-1276,1458)
        check(rem == 0 and t >= 0, "nonnegative integral terminal parameter")
        n, rem = divmod(8*2**a*(2*4**k-5)-207,243)
        check(rem == 0 and n % 2 == 1, "integral source")
        check(terminal_parameters(n)==(j,r),"source-only parameter decoder")
        for offset in (-4,-2,2,4):
            check(terminal_parameters(n+offset) is None,"neighbour cannot inherit source receipt")
        child = 223+256*t
        common = (4**k-1)//3
        check(n == (2264*2**a-23)//27+96*2**a*t, "inherited marked bridge")
        check(425+486*t == common, "terminal common future")
        aa = (1,2,a,1,2*k)
        states = replay(n,aa)
        child_states = replay(child,(1,1,1,1,3,2*k))
        check(states[-2] == child_states[-2] == common, "same marked endpoint")
        K = sum((b-1)//2 for b in aa)
        check(K == 9*j+k-2, "counter formula")
        f = Q(1,4**K)*Q(3,4)**4
        check(f == Q(81,4**(9*j+k+2)), "quarter price formula")
        W = Q(240,prod(range(K+1,K+7)))
        check(W >= Q(9,16)*f*f, "inherited beta lower transport")
        Z = Q(3*n+1,4)*f
        check(Z == Q(1,4)-Q(5,8*4**k)-Q(63,32*2**a*4**k),
              "normalized defect identity")
        check(Z > Q(1,5), "uniform family floor")
        m, q1 = (n-1)//2, (3*n+1)//2
        check(nu(m)*f > Q(12,5*q1**8), "independent source mass floor")
        for d in range(1,5):
            q = q1*(4**d-1)//3
            A = nu(m)+Q(1,4**d)*nu(m+q)
            E = Q(9,448*q**8)
            # M >= A*f follows from the source and its actual sibling.
            lower = f-E/A
            check(lower > Q(1777,1792)*f, "floor survives fixed-source refinement")
            check(A*f-E > 0, "positive signed measurement margin")
        if j == r == 1:
            rows.append({"j":j,"r":r,"k":k,"source":str(n),
                         "source_bits":n.bit_length(),"word":aa,"K":K,
                         "f":str(f),"W":str(W),
                         "refinement_relative_floor_lower":"1777/1792"})
    # Cheap hostile: the smallest neighbouring exponent does not even meet
    # the terminal congruence. No absent child is silently certified.
    for k in (1,2,240,241,243,244,484,486,487):
        check((4**k-1276) % 1458 != 0, "wrong terminal phase rejected")
    check(pow(4,485,1458)==1276, "second legal phase retained as positive control")
    for n in (0,1,2,3,27,223,True,Q(3,2)):
        check(terminal_parameters(n) is None,"outside-family decoder boundary")
    return {"universe":"j,r=1..6; source-preserving depths d=1..4",
            "families_tested":36,"examples":rows,"universal_coverage":"OPEN"}


# Exact Q(zeta5, sqrt(-55)) arithmetic, a subfield of Q(zeta55).
# A vector is a+b*w, w^2=-55, each coefficient in Q[z]/Phi5.
ZERO = (Q(0),)*8
ONE = (Q(1),)+(Q(0),)*7


def add(a,b):
    return tuple(x+y for x,y in zip(a,b))


def neg(a):
    return tuple(-x for x in a)


def scale(a,c):
    return tuple(c*x for x in a)


def mul4(a,b):
    c = [Q(0)]*7
    for i,j in product(range(4),repeat=2):
        c[i+j] += a[i]*b[j]
    for d in range(6,3,-1):
        for j in range(1,5):
            c[d-j] -= c[d]
    return tuple(c[:4])


def mul(a,b):
    ac,bd = mul4(a[:4],b[:4]),mul4(a[4:],b[4:])
    ad,bc = mul4(a[:4],b[4:]),mul4(a[4:],b[:4])
    return tuple(x-55*y for x,y in zip(ac,bd))+tuple(x+y for x,y in zip(ad,bc))


def power(a,k):
    v = ONE
    for _ in range(k):
        v = mul(v,a)
    return v


ALPHA = (Q(0),Q(1),Q(0),Q(0))+(Q(0),)*4
AP = [power(ALPHA,i) for i in range(5)]


def conj(a):
    v = ZERO
    for i in range(4):
        v = add(v,scale(AP[(-i)%5],a[i]))
        p = AP[(-i)%5]
        v = add(v,tuple(Q(0) for _ in range(4))+tuple(-a[i+4]*x for x in p[:4]))
    return v


def norm(a):
    return mul(a,conj(a))


def geometric():
    eta = (Q(-3,8),Q(0),Q(0),Q(0),Q(1,8),Q(0),Q(0),Q(0))
    check(norm(eta) == ONE,"quadratic nonintegral unit")
    walk = [ZERO,ONE,add(ONE,eta),scale(ONE,Q(1,4))]
    for i,j in combinations(range(4),2):
        check((norm(add(walk[i],neg(walk[j]))) == ONE) == (j == i+1),
              "induced three-edge rational-quarter path")
    x = s.symbols('x')
    P = s.cyclotomic_poly(5,x)
    u = s.rem(x*(x-3)*s.invert(1-3*x,P,x),P,x)
    U = tuple(Q(u.coeff(x,i)) for i in range(4))+(Q(0),)*4
    check(norm(U) == ONE,"conjugate unit direction")
    us = []
    for a in range(1,5):
        v = ZERO
        for i in range(4):
            v = add(v,scale(AP[(a*i)%5],U[i]))
        check(norm(v) == ONE,"all four Galois unit directions")
        us.append(v)
    trace = ZERO
    for v in us:
        trace = add(trace,v)
    check(trace == scale(ONE,Q(419,121)),"exact unit trace")
    check(add(scale(trace,11),scale(ONE,-38)) == scale(ONE,Q(1,11)),
          "82-edge rational-eleventh walk")
    # Integral 3-colouring on the tensor basis alpha^i beta^j.
    a_values,b_values = [1]*4+[-1],[1]*10+[-1]
    check(sum(a_values)%3 == sum(b_values)%3 == 0,"cyclotomic relations in F3")
    roots = {(i,j) for i,j in product(range(5),range(11))}
    for sign,(i,j) in product((-1,1),roots):
        check((sign*a_values[i]*b_values[j])%3 != 0,"all 110 integral directions")
    P55 = s.cyclotomic_poly(55,x)
    ell_power = [a_values[k%5]*b_values[(-2*k)%11] for k in range(40)]
    for k in range(110):
        reduced = s.Poly(s.rem(x**k,P55,x),x)
        observed = sum(int(reduced.nth(i))*ell_power[i] for i in range(40))%3
        check(observed == (a_values[k%5]*b_values[(-2*k)%11])%3,
              "independent power-basis colour reduction")
    check(s.totient(165) == 2*s.totient(55),"no primitive cube root in field")
    # Complete induced Cartesian incidence, not only listed edges.
    pentagon = [ZERO]
    for i in range(4):
        pentagon.append(add(pentagon[-1],AP[i]))
    chosen = None
    for exponent in range(1,9):
        rotation = power(eta,exponent)
        pts = [add(v,scale(rotation,j)) for v in pentagon for j in range(3)]
        if len(set(pts)) != 15:
            continue
        ok, edges = True,[]
        for a,b in combinations(range(15),2):
            i,j = divmod(a,3); k,l = divmod(b,3)
            wanted = (j == l and (i-k)%5 in (1,4)) or (i == k and abs(j-l)==1)
            actual = norm(add(pts[a],neg(pts[b]))) == ONE
            if actual != wanted:
                ok = False
                break
            if actual:
                edges.append([a,b])
        if ok:
            chosen = exponent
            break
    check(chosen is not None,"exact field-preserving Cartesian realization")
    for a,b in combinations(range(15),2):
        i,j = divmod(a,3); k,l = divmod(b,3)
        wanted = (j == l and (i-k)%5 in (1,4)) or (i == k and abs(j-l)==1)
        check((norm(add(pts[a],neg(pts[b]))) == ONE) == wanted,
              "all Cartesian pairs including nonedges")
    # Every labelled proper three-colouring on the original C5 extends;
    # the formula is checked on every edge, not just an orbit representative.
    seed_count = 0
    for c in product(range(3),repeat=5):
        if any(c[i] == c[(i+1)%5] for i in range(5)):
            continue
        seed_count += 1
        extension = [(c[i]+j)%3 for i in range(5) for j in range(3)]
        check(all(extension[a]!=extension[b] for a,b in edges),"pinned colouring lift")
        check(extension[::3] == list(c),"original marked fibre unchanged")
    check(seed_count == 30,"all labelled C5 three-colourings")
    # A deliberate seam with a known field solution, retaining its square
    # class. Such a solution alone does not rule out vertex collisions.
    A1,B1,rotation = add(ONE,ALPHA),ONE,neg(ONE)
    check(norm(add(A1,mul(rotation,B1)))==ONE,"positive deliberate-seam control")
    C = add(ONE,neg(add(norm(A1),norm(B1))))
    D = mul(norm(A1),norm(B1))
    w = mul(mul(conj(A1),B1),rotation)
    root55 = (Q(0),)*4+(Q(1),)+(Q(0),)*3
    h = mul(add(w,neg(conj(w))),scale(root55,Q(-1,55)))
    check(h == conj(h),"seam square-root lies in real subfield")
    check(add(scale(D,4),neg(mul(C,C))) == scale(mul(h,h),55),
          "exact square-class seam criterion")
    # All ordered pentagon displacements are roots of unity times 1 or phi.
    phi = add(ONE,add(ALPHA,conj(ALPHA)))
    mu10 = AP+[neg(v) for v in AP]
    differences = {add(a,neg(b)) for a,b in product(pentagon,repeat=2) if a!=b}
    check(differences == set(mu10+[mul(phi,v) for v in mu10]),
          "complete pentagon displacement spectrum")
    for eps,difference in product(mu10,differences):
        check(mul(eps,difference) in differences,"pentagon difference root symmetry")
    # Every non-(1,1) chord-length pair demands the same fifth-root trace.
    for r,tau in ((ONE,phi),(phi,ONE),(phi,phi)):
        check(add(ONE,neg(add(mul(r,r),mul(tau,tau)))) == neg(mul(phi,mul(r,tau))),
              "pentagon cross-edge trace is minus phi")
    golden_discriminant = add(scale(phi,7),neg(ONE))
    # Norm from Q(sqrt5) of 7phi-1 is -55, so it has mixed real signs.
    phi_conjugate = add(ONE,neg(phi))
    check(mul(golden_discriminant,add(scale(phi_conjugate,7),neg(ONE))) == scale(ONE,-55),
          "mixed-sign golden norm blocks the length-two seam")
    # Positive nonintegral seam: two unit circles with rational centres3/4
    # apart. Their two intersection points form a faithful C4.
    rhombus = [ZERO,neg(eta),scale(ONE,Q(3,4)),neg(conj(eta))]
    for i,j in combinations(range(4),2):
        check((norm(add(rhombus[i],neg(rhombus[j])))==ONE) == ((i-j)%4 in (1,3)),
              "all six nonintegral rhombus incidences")
    return {"quarter_path":[[str(c) for c in p] for p in walk],
            "trace_unit":str(u),"unit_trace":"419/121","eleventh_walk_length":82,
            "cartesian":{"rotation_power":chosen,"vertices":15,"edges":edges,
                         "all_seed_colourings_extended":seed_count},
            "pentagon_difference_count":len(differences),
            "golden_length_two_discriminant_norm":-55,
            "rational_rhombus_diagonal":"3/4",
            "six_colour_forcing":"OPEN"}


def primes():
    fixed,allowed = [],[]
    for p in s.primerange(2,501):
        if p in (5,11):
            actual = False
        else:
            order = int(s.n_order(p,55))
            actual = any(pow(p,k,55)==54 for k in range(order))
            criterion = p%5 == 4 and s.legendre_symbol(p,11) == -1
            check(actual == criterion,"conjugation-fixed prime criterion")
        if actual:
            fixed.append([int(p),int(p)**(int(s.n_order(p,55))//2)])
        else:
            allowed.append(int(p))
    check(min(q for p,q in fixed)==109,"smallest fixed-prime norm base control")
    return {"tested_primes_through":500,"inverted_rational_primes_through_43":
            [p for p in allowed if p<=43],"fixed_primes_through_109":
            [[p,q] for p,q in fixed if p<=109]}


def timeform():
    x,t,E = s.symbols('x t E')
    H = x+x**4*t**5
    multiplier = s.diff(H,x)
    check(s.expand((1-4*x**3*t**5)*s.diff(H,x)+Q(16,5)*x**2*t**6*s.diff(H,t))==1,
          "inherited gradient Bezout control")
    check(s.cancel(multiplier.subs(t**5,(E-x)/x**4)-(4*E/x-3))==0,
          "exact clock multiplier on generic fibre")
    resultant = s.factor(s.resultant(H-E,multiplier,x))
    affine_zero = 256*E**3*t**5+27
    check(resultant == -t**15*affine_zero,"full resultant including degree-drop factor")
    check(multiplier.subs(t,0)==1,"degree-drop factor has no affine clock zero")
    # Orders at the three branch points of the cyclic degree-five x-cover.
    # At infinity x has pole order5 and t has zero order3.
    divisor = {"t":{"x=0":-4,"x=E":1,"x=infinity":3},
               "dt":{"x=0":-5,"x=E":0,"x=infinity":2,
                     "five_affine_simple_zeros":5}}
    check(sum(divisor['t'].values())==0,"function divisor degree zero")
    check(sum(divisor['dt'].values())==2,"genus-two canonical degree")
    check(2*2-2 == -2*5+3*(5-1),"degree-five Riemann-Hurwitz")
    # A finite truncation can only move the defect, never remove it.
    q, coeff = t,s.Integer(1)
    defects=[]
    for K in range(7):
        J = s.expand(s.diff(H,x)*s.diff(q,t)-s.diff(H,t)*s.diff(q,x))
        expected = coeff*(5*K+4)*x**(3*K+3)*t**(5*K+5)
        check(s.expand(J-1-expected)==0,"inherited character-one truncation defect")
        defects.append(str(expected))
        coeff *= -s.Rational(5*K+4,5*K+6)
        q += coeff*x**(3*K+3)*t**(5*K+6)
    return {"carrier":"H=x+x^4*t^5","exact_clock":"t", "multiplier":str(multiplier),
            "full_resultant":str(resultant),"affine_zero_equation":str(affine_zero),"divisor":divisor,
            "truncation_defects":defects,"standard_area_mate":"excluded by THM-3418/3716"}


if __name__ == '__main__':
    output = {"collatz":collatz(),"geometry":geometric(),"localization":primes(),
              "planar":timeform()}
    output['checks'] = checks
    print(json.dumps(output,sort_keys=True,indent=2))
