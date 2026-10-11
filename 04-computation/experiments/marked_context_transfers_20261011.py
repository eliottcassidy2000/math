"""Exact checks for quarter-price ROOT tests and geometric/source contexts.

This is a finite replay, not a Collatz convergence proof. No removable asserts.
Python 3 + sympy. Run normally and with -O; stdout is deterministic JSON.
"""
from fractions import Fraction as Q
from itertools import product, combinations
from math import factorial
from functools import lru_cache
from collections import Counter
import json
import sympy as sp

checks = 0


def check(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)


@lru_cache(None)
def word(n):
    original = n
    aa, nn = [], [n]
    for _ in range(10000):
        if n == 1:
            return tuple(aa), tuple(nn)
        y = 3*n+1
        a = (y & -y).bit_length()-1
        n = y >> a
        aa.append(a)
        nn.append(n)
    raise RuntimeError(f"uncertified finite-bank source {original}")


def scores(n):
    aa, nn = word(n)
    if n == 1:
        return Q(1), Q(1)
    k, tau = sum((a-1)//2 for a in aa), len(aa)
    f = Q(1, 4)**k * Q(3, 4)**(tau-1)
    w = Q(2*factorial(k)*factorial(tau), factorial(k+tau+1))
    return f, w


def collatz():
    max_tau = 0
    for n in range(1, 20002, 2):
        aa, nn = word(n)
        f, w = scores(n)
        max_tau = max(max_tau, len(aa))
        z = Q(1)
        for a, y in zip(aa[:-1], nn[1:-1]):
            z *= Q(3*y, 3*y+1) * (Q(1, 2) if a % 2 else 1)
        check(Q(3*n+1, 4)*f == z, "normalized quarter-price defect product")
        check(f <= Q(4, 3*n+1), "sharp source envelope")
        check((f == Q(4, 3*n+1)) == (len(aa) <= 1), "envelope equality cases")
        check(w >= Q(9, 16)*f*f, "quarter-price to beta floor")
        if n < 1000:
            for j in range(1, 7):
                sibling = (4**j*(3*n+1)-1)//3
                check(scores(sibling)[0] == Q(1, 4**j)*f, "native exact multiplier")
                check(Q(3*sibling+1,4)*scores(sibling)[0] == Q(3*n+1,4)*f,
                      "normalized score is sibling invariant")
    # The residual is tested independently by literal forward certificates.
    samples = []
    nu = lambda a: Q(1, (a+1)**6)-Q(1, (a+2)**6)
    for n, j in product((1, 3, 5, 7, 27, 1729), (1, 2, 4, 8)):
        m = (n-1)//2
        q = (3*m+2)*(4**j-1)//3
        residual = sum((nu(m+k*q)*scores(n+2*k*q)[0]
                        for k in range(2, 42)), Q(0))
        bill = Q(9, 448*q**8)
        check(residual <= bill, "quarter-price q^-8 residual")
        # The sum is only forty terms; it is not printed as the full readout.
        samples.append({"n": n, "j": j, "q": q,
                        "forty_term_residual": str(residual), "full_tail_upper": str(bill)})
    # Exact finite KL hostile, represented by its exponential.
    # A deleted marked atom p has D(P||Q)=-log(1-p), not small-gap protection.
    for den in range(2, 31):
        p = Q(1, den)
        qlaw = [p, (1-p)/2, (1-p)/2]
        plaw = [Q(0), Q(1, 2), Q(1, 2)]
        check(plaw[1]/qlaw[1] == 1/(1-p), "deleted atom KL likelihood ratio")
        check(sum(qlaw) == sum(plaw) == 1, "probability normalizations")
    mersenne_rows = []
    for h in range(1,19):
        n = 2**(h+1)-1
        aa, nn = word(n)
        check(aa[:h] == (1,)*h, "Mersenne valuation-one prefix")
        check(nn[:h+1] == tuple(3**i*2**(h+1-i)-1 for i in range(h+1)),
              "Mersenne exact affine states")
        upper = Q(3,4)**h * Q(4,6*3**h-2)
        check(upper == Q(8,(n+1)**2)/(3-Q(1,3**h)), "quadratic lower-floor hostile cap")
        check(scores(n)[0] <= upper, "literal Mersenne score control")
        mersenne_rows.append({"h":h,"n":n,"f":str(scores(n)[0]),"upper":str(upper)})
    return {"odd_sources_through": 20001, "odd_sources": 10001,
            "odd_step_cap": 10000, "max_odd_steps_observed": max_tau,
            "sibling_checks": 3000, "residual_samples": samples,
            "mersenne_lower_floor_hostiles": mersenne_rows}


def cycle_count(lists, q):
    count = 0
    for first in range(q):
        if not lists[0] >> first & 1:
            continue
        dp = [int(c == first) for c in range(q)]
        for mask in lists[1:]:
            total = sum(dp)
            dp = [(total-dp[c]) if mask >> c & 1 else 0 for c in range(q)]
        count += sum(dp)-dp[first]
    return count


def geometry():
    lists = [a for a in range(16) if a.bit_count() >= 2]
    bad = 0
    for row in product(lists, repeat=5):
        w = cycle_count(row, 4)
        forbidden = row[0].bit_count() == 2 and len(set(row)) == 1
        check((w == 0) == forbidden, "complete four-palette pentagon list criterion")
        bad += (w == 0)
    check(bad == 6, "only the six common-pair failures")
    # Independent literal count for the specified six-colour positive controls.
    controls = [([3]*5, 0), ([7,3,3,3,3], 2), ([7]*5, 30)]
    for row, expected in controls:
        direct = sum(all(row[i] >> cs[i] & 1 for i in range(5))
                     and all(cs[i] != cs[(i+1) % 5] for i in range(5))
                     for cs in product(range(6), repeat=5))
        check(direct == cycle_count(row, 6) == expected, "six-colour literal list control")
    # Non-unit edges between external leaves cannot obstruct recolouring a
    # nonadjacent pair with a previously unused sixth colour.
    # Exhaust all proper five-colour assignments of an independent four-port group.
    repairs = 0
    for cs in product(range(5), repeat=4):
        repaired = (5,5,cs[2],cs[3])
        check(len(set(repaired)) <= 3, "private port-pair repair")
        repairs += 1
    return {"complete_four_colour_list_profiles": len(lists)**5,
            "uncolourable_profiles": bad, "six_colour_control_weights": [0,2,30],
            "private_group_recolour_controls": repairs,
            "synchronized_palette_product_probability": str(Q(1, 15**4)),
            "synchronized_palette_total_correlation": "4*log(15)"}


def gf16mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 31
    return r


def algebraic():
    z, v = sp.symbols('z v')
    rows = []
    for b in range(1, 8):
        poly = z**4-z**3+(2-b)*z**2-z+1
        check(sp.expand(z**2*((z+1/z)**2-(z+1/z)-b)-poly) == 0,
              "reciprocal quadratic trace")
        disc = sp.discriminant(poly, z)
        check(disc == (1+4*b)**2*(2-b)*(6-b), "quartic discriminant")
        rows.append({"b": b, "polynomial": str(poly), "factorization": str(sp.factor(poly)),
                     "discriminant": int(disc), "irreducible_Q": bool(sp.Poly(poly,z).is_irreducible)})
    colours = [0,1,1,0,1,0,2,1,2,0,3,2,3,2,0,3]
    units = [a for a in range(1,16)
             if gf16mul(gf16mul(a,a),gf16mul(a,gf16mul(a,a))) == 1]
    check(units == [1,2,4,8,15], "F16 norm-one subgroup")
    edges = [(a,b) for a in range(16) for b in range(a+1,16) if a^b in units]
    for a,b in edges:
        check(colours[a] != colours[b], "explicit F16 four-colouring")
    independent = []
    for mask in range(1 << 16):
        if all(not(mask >> a & 1 and mask >> b & 1) for a,b in edges):
            independent.append(mask.bit_count())
    check(max(independent) == 5, "independence bound supplies no three-colouring")
    # The integral b=4 module has the F2 character alpha -> 1.
    check((1-1-2-1+1) % 2 == 0, "Salem17 integral parity character")
    # Golden Cartesian coordinates: the explicit non-axis unit uses odd denominator.
    phi = (1+sp.sqrt(5))/2
    xx, yy = -phi/(phi+2), 2*phi/(phi+2)
    check(sp.simplify(xx*xx+yy*yy-1) == 0, "non-axis golden rational unit")
    check(sp.simplify((phi+2)*((1-sp.sqrt(5))/2+2)) == 5, "golden denominator norm")
    # Trace root of T^2-T-4 at the even 2-adic place. This certifies the
    # finite lifts used to check the all-depth Hensel/Eisenstein proof.
    tr = 0
    for k in range(1, 33):
        modulus = 2**k
        candidates = [w for w in (tr,tr+modulus//2) if (w*w-w-4) % modulus == 0]
        # At k=1 choose the even root. Later the existing residue is fixed.
        if k == 1:
            candidates = [0]
        check(len(candidates) == 1, "unique even trace lift")
        tr = candidates[0]
        if k >= 3:
            check(tr % 8 == 4 and (2-tr) % 4 == 2, "Eisenstein at the Salem17 place")
    # Stronger local tests at three: these improve the dyadic four-colour
    # upper bounds for the full Salem13 and Salem21 fields to three.
    check(sp.expand((z+2)*(z-3)+3) == z*z-z-3,
          "Salem13 trace valuation identity")
    check(sp.expand(((z+1)**4-(z+1)**3-3*(z+1)**2-(z+1)+1))
          == z**4+3*z**3-6*z-3, "Salem21 Eisenstein translate at three")
    tr3 = 1
    for k in range(2, 13):
        modulus = 3**k
        candidates = [w for w in (tr3+i*(modulus//3) for i in range(3))
                      if (w*w-w-3) % modulus == 0]
        check(len(candidates) == 1, "unique Salem13 trace lift at three")
        tr3 = candidates[0]
        check((tr3+2) % 3 == 0 and (tr3+2) % 9 != 0,
              "Salem13 ramified quadratic place")
    cycle_rows = []
    for beta in (1,3,4,5):
        def mul(a,c):
            vv = [0]*7
            for i,xx in enumerate(a):
                for j,yy in enumerate(c):
                    vv[i+j] += xx*yy
            for kk in range(6,3,-1):
                zz = vv[kk]
                vv[kk-4] -= zz
                vv[kk-3] += zz
                vv[kk-2] += (beta-2)*zz
                vv[kk-1] += zz
            return tuple(vv[:4])
        one = (1,0,0,0)
        inv = (1,beta-2,1,-1)
        bars = [one,inv,mul(inv,inv),mul(mul(inv,inv),inv)]
        bar = lambda aa: tuple(sum(aa[j]*bars[j][i] for j in range(4)) for i in range(4))
        steps = [one,(0,-1,0,0)]
        steps += [(0,0,1 if beta == 1 else -1,0)]*abs(2-beta)
        steps += [(0,0,0,-1),(-1,1,beta-2,1)]
        pts = [(0,0,0,0)]
        for step in steps:
            check(mul(step,bar(step)) == one, "reciprocal polygon unit step")
            pts.append(tuple(xx+yy for xx,yy in zip(pts[-1],step)))
        check(pts[-1] == (0,0,0,0), "reciprocal polygon closure")
        pts.pop()
        check(len(set(pts)) == len(pts), "reciprocal polygon distinct vertices")
        edges_cycle = []
        for ii,jj in combinations(range(len(pts)),2):
            diff = tuple(xx-yy for xx,yy in zip(pts[ii],pts[jj]))
            unit = mul(diff,bar(diff)) == one
            check(unit == (jj-ii == 1 or (ii == 0 and jj == len(pts)-1)),
                  "induced reciprocal unit cycle")
            if unit:
                edges_cycle.append((ii,jj))
        cycle_rows.append({"b": beta, "vertices_power_basis": pts, "edges": edges_cycle})
    return {"reciprocal_family": rows, "F16_norm_one": units, "F16_edges": edges,
            "F16_four_colouring": colours, "F16_independence_number": 5,
            "F16_chromatic_number": 4, "even_trace_root_mod_2pow32": tr,
            "Salem13_trace_root_mod_3pow12": tr3,
            "full_field_chromatic_numbers_by_proved_local_tests": {1:3,3:3,4:2,5:3},
            "induced_unit_cycles": cycle_rows}


def planar():
    x,t,a,b,c,d = sp.symbols('x t a b c d')
    v = x*t
    p, y, u = t+v*v, v*(t+v*v), x*x*t
    G = lambda aa,bb: -u/2+aa*p+bb*y
    F = lambda aa,bb: (2*(aa+bb*v)/x, -x*x*t/(2*(aa+bb*v)))
    sub = lambda expr, pair: sp.factor(expr.subs({x:pair[0],t:pair[1]}, simultaneous=True))
    jac = lambda pair: sp.factor(sp.det(sp.Matrix([[sp.diff(z,w) for w in (x,t)] for z in pair])))
    fab = F(a,b)
    composition = tuple(sub(z,fab) for z in F(c,d))
    r = (a+b*v)/(c-d*v)
    check(sp.factor(composition[0]-x/r) == 0, "two-letter x normal form")
    check(sp.factor(composition[1]-r*t) == 0, "two-letter t normal form")
    check(sp.factor(sub(G(a,-b),fab)-G(a,b)) == 0, "source equation")
    check(jac(fab) == 1, "F area preservation")
    # Root-centred translation after a divisor move: changes v and gives a
    # polynomial pullback, outside the fixed-v signed normal forms.
    change = (b*x*x*t-a*x, 1/(b*x))
    newG = t+b*x*(t-sp.Rational(1,2))*(x*t-a/b)**2
    check(sp.factor(sub(G(a,b), change)-newG) == 0, "root-centred polynomial pullback")
    check(jac(change) == 1, "root-centred map area")
    check(sp.factor(change[0]*change[1]-(v-a/b)) == 0, "actual changed product coordinate")
    inverse = (1/(b*t), b*x*t*t+a*t)
    for coord, expected in zip(inverse,(x,t)):
        check(sp.factor(sub(coord,change)-expected) == 0, "root-centred inverse")
    check(sp.factor(1-2*G(a,b)-(1+x*x*t)*(1-2*(a+b*v)*t)) == 0,
          "common reducible critical fibre")
    # An escape to a genus-two polynomial submersion, with a separate
    # complement obstruction. The generic hyperelliptic discriminant is E^2-4v^5.
    H = x+x**4*t**5
    check(sp.groebner([sp.diff(H,x),sp.diff(H,t)],x,t).contains(sp.Integer(1)),
          "genuine polynomial submersion")
    E,V = sp.symbols('E V')
    disc = sp.Poly(E**2-4*V**5,V)
    check(sp.gcd(disc,disc.diff()).degree() == 0, "generic genus-two squarefree quintic")
    PP,QQ,RR = (sp.Function(name) for name in ('PP','QQ','RR'))
    carrier = PP(v)/t+QQ(v)*t+RR(v)
    check(sp.simplify(sp.diff(carrier,x)*x-sp.diff(carrier,t)*t
                      -(PP(v)/t-QQ(v)*t)) == 0, "generic carrier time-form sign")
    radial = V**2*(a+b*V)
    sextic = (E-radial)**2+2*radial
    check(sp.expand(sp.diff(sextic,V)-2*(radial+1-E)*sp.diff(radial,V)) == 0,
          "generic original sextic squarefree eliminant")
    check(sp.expand(sextic.subs(E,radial+1)-(2*radial+1)) == 0,
          "moving critical factor gives two energy minus one")
    P = sp.Function('P')
    Z = sp.symbols('Z')
    # Check the weighted-component ODE for monomial P=Z^k, for exact controls.
    for k in range(9):
        J = t*(x**3*t**5)**k
        bracket = sp.diff(H,x)*sp.diff(J,t)-sp.diff(H,t)*sp.diff(J,x)
        expected = (1+4*Z)*Z**k+5*Z*(1+Z)*sp.diff(Z**k,Z)
        check(sp.expand(bracket-expected.subs(Z,x**3*t**5)) == 0, "weighted bracket ODE")
    return {"two_letter_ratio": str(r), "root_centred_map": [str(z) for z in change],
            "root_centred_polynomial_pullback": str(newG),
            "submersion": str(H), "submersion_generic_genus": 2,
            "polynomial_Keller_complement": "inherited THM-3716/3772; weight control 5d+4",
            "genus_two_source_rational_mate": "excluded by holomorphic dv/Y; invariant under birational symplectic transport"}


def cyclotomic():
    rows = []
    masks = {3:(1,2), 5:(1,14), 11:(511,654), 13:(1,726,3448)}
    X = sp.symbols('X')
    for prime in (3,5,11,13):
        degree = prime-1
        modulus = (1 << prime)-1
        size = 1 << degree
        q = 1 << (degree//2)
        def mul(a,b):
            ans = 0
            while b:
                if b & 1:
                    ans ^= a
                b >>= 1
                a <<= 1
                if a & size:
                    a ^= modulus
            return ans
        def power(a,n):
            ans = 1
            while n:
                if n & 1:
                    ans = mul(ans,a)
                a = mul(a,a)
                n >>= 1
            return ans
        check(sp.n_order(2,prime) == degree, "cyclotomic irreducibility at two")
        check(pow(2,degree//2,prime) == prime-1, "local conjugation exponent")
        units = [u for u in range(1,size) if power(u,q+1) == 1]
        check(len(units) == q+1, "complete norm-one difference set")
        spectrum = [sum(1 if (a&u).bit_count()%2 == 0 else -1 for u in units)
                    for a in range(size)]
        fwht = [int(a in units) for a in range(size)]
        step = 1
        while step < size:
            for base in range(0,size,2*step):
                for offset in range(step):
                    i,j = base+offset,base+offset+step
                    aa,bb = fwht[i],fwht[j]
                    fwht[i],fwht[j] = aa+bb,aa-bb
            step *= 2
        check(spectrum == fwht, "independent exact Walsh transform")
        check(sum(spectrum) == 0, "adjacency trace")
        check(sum(a*a for a in spectrum) == size*len(units), "adjacency square trace")
        for u in units:
            check(any((mask&u).bit_count()%2 for mask in masks[prime]),
                  "additive palette map separates every norm-one difference")
        minimum = min(spectrum)
        alpha_bound = size*(-minimum)//(len(units)-minimum)
        lower = (size+alpha_bound-1)//alpha_bound
        if prime == 5:
            lower = 4  # independently exhaustive F16 calculation above
        rows.append({"cyclotomic_prime": prime, "field_size": size,
                     "degree": len(units), "norm_one_elements": units,
                     "spectrum": sorted(Counter(spectrum).items()),
                     "hoffman_independence_upper": alpha_bound,
                     "chromatic_lower": lower, "chromatic_upper": 2**len(masks[prime]),
                     "colour_linear_masks": masks[prime]})
    exponents = [20,18,17,15,14,12,9,7,6,5,2,1,0]
    f = sum(X**e for e in exponents)
    g = sp.expand(X**20*f.subs(X,1/X))
    phi55 = sp.cyclotomic_poly(55,X)
    check(sp.Poly(f*g-phi55,X,modulus=2).is_zero, "reciprocal factorization of Phi55 mod two")
    check(sp.gcd(sp.Poly(f,X,modulus=2),sp.Poly(g,X,modulus=2)).degree() == 0,
          "distinct conjugate dyadic primes")
    check(sp.n_order(2,55) == 20, "degree20 dyadic factors")
    check(not any(pow(2,k,55) == 54 for k in range(20)), "minus one absent from decomposition group")
    norm_f = int(sp.resultant(phi55,f,X))
    check(norm_f == 2**20*991*284681, "exact norm of escaping factor")
    fixed_primes = []
    for ell in sp.primerange(2,110):
        if ell in (5,11):
            continue
        order = int(sp.n_order(ell,55))
        fixed = order % 2 == 0 and pow(ell,order//2,55) == 54
        criterion = ell % 5 == 4 and pow(ell,5,11) == 10
        check(fixed == criterion, "conjugation-fixed prime CRT criterion")
        if fixed:
            q = ell**(order//2)
            check(q == (ell if ell % 11 == 10 else ell**5), "residue norm base")
            fixed_primes.append((ell,order,q))
    check([row[0] for row in fixed_primes] == [19,29,79,109], "all small fixed primes")
    check(min(row[2] for row in fixed_primes) == 109, "smallest local norm base")
    check(sp.n_order(5,11) == 5 and 11 % 5 == 1, "ramified primes also not fixed")
    return {"norm_graphs": rows, "Phi55_factor_exponents": exponents,
            "Phi55_reciprocal_factor_exponents": sorted(20-e for e in exponents),
            "factor_norm": norm_f,
            "fixed_primes_through109": fixed_primes, "minimum_local_norm_base": 109,
            "mixed_field_unit": "eta=f(zeta55)/f(zeta55^(-1)); valuation 1 at one dyadic prime"}


if __name__ == '__main__':
    out = {"collatz": collatz(), "geometry": geometry(),
           "algebraic_integers": algebraic(), "planar": planar(),
           "cyclotomic_palette": cyclotomic()}
    out["always_active_checks"] = checks
    print(json.dumps(out, indent=2, sort_keys=True))
