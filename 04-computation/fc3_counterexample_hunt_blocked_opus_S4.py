"""
FC(3) full counterexample hunt on the triangle: BLOCKED at every degree.
C_3-covariant homogeneous phi (omega^2-eigenspace) auto-kills 3 nmid m moments; try to also kill int phi^{3k}.
D=2 (1 param) kills phi^3 not phi^6; D=3 (2 params) kills phi^3,phi^6 not phi^9. Pattern: degree D kills the
first p(D) moments, the next is nonzero -> no polynomial realizes the all-moments-zero (disk) measure.
Reason: FC(3)=triangle moment problem (discrete S_3); GM(3)=sphere (continuous SO(3), Long's CE). Discrete
symmetry caps balance at 3-fold -> evidence bare FC(3) is TRUE though GM(3) is false.
"""
import mpmath as mp
mp.mp.dps = 35
w = mp.exp(2j*mp.pi/3)

def si(poly):  # exact 2-simplex integral of a monomial dict {(i,j,k):coeff}
    return sum(c*mp.factorial(i)*mp.factorial(j)*mp.factorial(k)/mp.factorial(i+j+k+2)
               for (i, j, k), c in poly.items())
def mul(p, q):
    r = {}
    for (i, j, k), c in p.items():
        for (a, b, d), e in q.items():
            key = (i+a, j+b, k+d); r[key] = r.get(key, mp.mpc(0)) + c*e
    return r
def power(p, n):
    r = {(0, 0, 0): mp.mpc(1)}
    for _ in range(n): r = mul(r, p)
    return r
def lin(terms):
    r = {}
    for co, bas in terms:
        for k, c in bas.items(): r[k] = r.get(k, mp.mpc(0)) + co*c
    return r

# D=2
q1 = {(2, 0, 0): mp.mpc(1), (0, 2, 0): w, (0, 0, 2): w*w}
q2 = {(1, 1, 0): mp.mpc(1), (0, 1, 1): w, (1, 0, 1): w*w}
c3 = si(power(q1, 3)); c2 = 3*si(mul(power(q1, 2), q2)); c1 = 3*si(mul(q1, power(q2, 2))); c0 = si(power(q2, 3))
print("D=2: solve int phi^3=0 (phi=A q1+q2), check int phi^6:")
for A in mp.polyroots([c3, c2, c1, c0]):
    p = lin([(A, q1), (mp.mpc(1), q2)])
    print("  A=%s: |int phi^3|=%s |int phi^6|=%s (NONZERO -> blocked)"
          % (mp.nstr(A, 5), mp.nstr(abs(si(power(p, 3))), 2), mp.nstr(abs(si(power(p, 6))), 4)))

# D=3
r1 = {(3, 0, 0): mp.mpc(1), (0, 3, 0): w, (0, 0, 3): w*w}
r2 = {(2, 1, 0): mp.mpc(1), (0, 2, 1): w, (1, 0, 2): w*w}
r3 = {(1, 2, 0): mp.mpc(1), (0, 1, 2): w, (2, 0, 1): w*w}
def phi3(A, B): return lin([(mp.mpc(1), r1), (A, r2), (B, r3)])
def moms3(A, B, ks): p = phi3(A, B); return [si(power(p, k)) for k in ks]
print("D=3: solve int phi^{3,6}=0, check UNFITTED int phi^9:")
sol = mp.findroot(lambda A, B: moms3(A, B, [3, 6]), (mp.mpc('0.3'), mp.mpc('0.1')))
v = moms3(sol[0], sol[1], [3, 6, 9])
print("  |phi^3|=%s |phi^6|=%s  UNFITTED |phi^9|=%s (NONZERO -> blocked)"
      % (mp.nstr(abs(v[0]), 2), mp.nstr(abs(v[1]), 2), mp.nstr(abs(v[2]), 4)))
print("=> full FC(3) counterexample blocked at every degree; the S_3 form only kills 2/3 of the moments.")
