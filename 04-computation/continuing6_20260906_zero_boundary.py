"""Exact all-phase certificate on the anchored zero-root boundary.

No producer imports. Reconstruct the original binomial-14 Hadamard carriers,
derive the two C/B Gram constraints, then verify complete positive polynomial
coefficient certificates. No floating arithmetic or numerical optimization.
"""
from math import comb
from pathlib import Path
import hashlib
import json
import sys
import sympy as S

sys.stdout.reconfigure(newline="\n")
HERE = Path(__file__).resolve().parent
DEST = HERE.parent / "05-knowledge" / "results" if HERE.name == "04-computation" else HERE
GATES = 0


def require(ok, label):
    global GATES
    if not ok:
        raise RuntimeError(label)
    GATES += 1


def same(a, b, label):
    require(S.expand(a-b) == 0, label)


x, y, z, s, u, v, X = S.symbols("x y z s u v X")


def product(a, b):
    out = {}
    for i, ai in a.items():
        for j, bj in b.items():
            out[i+j] = out.get(i+j, 0)+ai*bj
    return {i:S.expand(ai) for i, ai in out.items()}


def add(a, b):
    return {i:S.expand(a.get(i, 0)+b.get(i, 0)) for i in set(a)|set(b)}


def shift(a, d, c=1):
    return {i+d:c*ai for i, ai in a.items()}


O = {j:comb(14, 2*j+1) for j in range(7)}
E = {j:comb(14, 2*j) for j in range(8)}
beta = dict(zip(range(-1, 5), [1, 13, 55, x, y, z]))
cr = dict(zip(range(-1, 4), [1, 12, 45, 2*x/3, 3*y/7]))
dr = dict(zip(range(-1, 4), [1, 11, 36, 5*x/12, y/7]))
P = {i:S.expand(O[i]*beta[i]) for i in O.keys() & beta.keys()}
qa = add(product(O, O), shift(product(E, E), -1))
qb = add(product(beta, beta), shift(product(cr, dr), 1, 2))
Q = {i:S.expand(qa[i]*qb[i]) for i in qa.keys() & qb.keys()}
same(Q[-1], 28, "original inverse carry")
p = S.expand(sum(pi*(-s)**i for i, pi in P.items()).subs(z, 0)/2002)
same(p, -S.Rational(12, 7)*y*s**3+x*s*s-10*s+S.Rational(1, 11), "original zero-boundary phase")
T = S.expand(sum(qi*(-s)**i*s for i, qi in Q.items()).subs(z, 0))
require(S.Poly(T, s).degree() == 8, "zero-boundary response degree with mixed carry")
same(S.diff(T, x, 2), 397670*s**5*(66-85*s), "raw separate convexity in x")
same(S.diff(T, y, 2), S.Rational(4011660, 7)*s**7*(140-13*s), "raw separate convexity in y")

# Exact residue moments; the packet is necessary for C interlacing alone.
den = [1, -13, 55, -x, y]
numC = [1, -12, 45, -2*x/3, 3*y/7]
mu = []
for k in range(9):
    mu.append(S.expand((numC[k] if k < 5 else 0)-sum(den[j]*mu[k-j] for j in range(1, min(k, 4)+1))))
    same(sum(den[j]*mu[k-j] for j in range(min(k, 4)+1)),
         numC[k] if k < 5 else 0, "formal C/B division coefficient")
shifted2 = S.Matrix(2, 2, lambda i, j:mu[i+j+1]).det()
ordinary3 = S.Matrix(3, 3, lambda i, j:mu[i+j]).det()
cap = S.Rational(7, 72)*(x-75)*(135-x)
same(shifted2, (x-75)/3, "C shifted H2")
same(ordinary3, (x-75)*(135-x)/9-8*y/7, "C ordinary H3")
same(cap, S.Rational(175, 2)-S.Rational(7, 72)*(x-105)**2, "global packet rectangle")

# Solve exactly the original phase, without changing its variable or carriers.
Y = 7*x/(12*s)-S.Rational(35, 6)/s**2+S.Rational(7, 132)/s**3
require(S.cancel(p.subs(y, Y)) == 0, "original phase elimination")
F = S.Poly(S.cancel(-968*T.subs(y, Y)), x, s).as_expr()
expectedF = (5182210385*s**6*x*x+8439169200*s**5*x*x
             -235681249320*s**5*x-585418008010*s**4*x
             +2361406947900*s**4+38085709200*s**3*x
             +6069283306240*s**3-385825440*s*s*x
             -498310898015*s*s+8357151600*s-38651536)
same(F, expectedF, "complete eleven-term restricted response")
require(len(S.Poly(F, x, s).terms()) == 11, "complete restricted response support")
A = S.expand(F).coeff(x, 2)
same(A, 264385*s**5*(19601*s+31920), "strict quadratic convexity of F")
J = s**3*(x*x-210*x+10125)+6*x*s*s-60*s+S.Rational(6, 11)
require(S.cancel(S.Rational(72, 7)*s**3*(Y-cap)-J) == 0, "necessary phase-packet inequality J<=0")
print("Original P/Q rebuilt with binomial14 and inverse carry28; sQ=-F/968 on the original zero phase.")
print("Packet: x>=75, 0<=y<=(7/72)(x-75)(135-x); hence x<=135,y<=175/2.")

# Uniform root coverage from the packet; all endpoint inequalities are exact.
left1, right1 = S.Rational(1, 110), S.Rational(1, 90)
left2, right2 = S.Rational(63, 1000), S.Rational(1, 8)
left3, tail = S.Rational(9, 16), S.Rational(3, 5)
lo1 = p.subs({x:75, y:S.Rational(175, 2), s:left1})
hi1 = p.subs({x:135, y:0, s:right1})
hi2 = p.subs({x:135, y:0, s:left2})
require(lo1 > 0 and hi1 < 0 and hi2 < 0, "uniform rectangular phase endpoints")
phase_bounds = [str(lo1), str(hi1), str(hi2)]
for ss, wanted in ((right2, S.Rational(3, 2816)), (left3, S.Rational(3335, 22528))):
    lower = S.expand(p.subs({y:cap, s:ss}))
    centre = 105-3/ss
    minimum = lower.subs(x, centre)
    same(minimum, wanted, "packet phase minimum")
    same(lower, minimum+ss**3*(x-centre)**2/6, "complete square for phase endpoint")
    require(minimum > 0, "strict positive phase endpoint")
    phase_bounds.append(str(minimum))
print("Uniform phases: (1/110,1/90), (63/1000,1/8), and (9/16,infinity) for y>0; first two for y=0.")

certificates = []
coefficient_count = 0


def interval_positive(label, poly, a, b, degree=None):
    global coefficient_count
    pp = S.Poly(poly, s)
    degree = pp.degree() if degree is None else degree
    require(pp.degree() <= degree, "declared homogeneous degree")
    transformed = S.Poly(S.expand(sum(pp.nth(i)*(a+b*u)**i*(1+u)**(degree-i)
                                    for i in range(degree+1))), u)
    cc = [transformed.nth(i) for i in range(degree+1)]
    for c in cc:
        require(c > 0, "strict positive interval coefficient: "+label)
    coefficient_count += len(cc)
    certificates.append({"label":label, "interval":[str(a), str(b)], "degree":degree,
                         "power_coefficients":[str(pp.nth(i)) for i in range(degree+1)],
                         "homogeneous_coefficients":[str(c) for c in cc]})
    print(f"{label}: {len(cc)} positive coefficients, minimum {min(cc)}.")


def quadrant_positive(label, x0, s0):
    global coefficient_count
    shifted = S.Poly(S.expand(F.subs({x:X+x0, s:u+s0})), X, u)
    coefficients = []
    for i in range(3):
        for j in range(7):
            c = shifted.coeff_monomial(X**i*u**j)
            require(c > 0, "complete strictly positive quadrant coefficient: "+label)
            coefficients.append({"x_degree":i, "s_degree":j, "coefficient":str(c)})
    coefficient_count += 21
    certificates.append({"label":label, "quadrant_base":[str(x0), str(s0)], "coefficients":coefficients})
    print(f"{label}: all21 shifted coefficients positive.")


# First branch: raw T is separately convex throughout the packet rectangle.
require(66-85*right1 > 0 and 140-13*right1 > 0, "first-branch separate convexity signs")
for xx in (75, 135):
    for yy in (0, S.Rational(175, 2)):
        interval_positive(f"raw first corner x={xx},y={yy}", -T.subs({x:xx, y:yy}), left1, right1, 8)

# Second branch: a discriminant region, then an admissible-x barrier.
disc = S.discriminant(F, x)
disc_reduced = S.cancel(disc/s**4)
require(S.Poly(disc_reduced, s).degree() == 6, "complete restricted discriminant degree")
mid2 = S.Rational(11, 100)
interval_positive("negative discriminant", -disc_reduced, left2, mid2)
interval_positive("J at x=105", J.subs(x, 105), mid2, right2)
interval_positive("F at x=105", F.subs(x, 105), mid2, right2)
interval_positive("negative Fx at x=105", -S.diff(F, x).subs(x, 105), mid2, right2)
same(S.diff(J, x), s*s*(2*s*(x-105)+6), "packet derivative orientation")
same(S.diff(F, x), 2*A*(x-105)+S.diff(F, x).subs(x, 105), "response derivative orientation")

# Third branch: the phase packet forces x>84 before the x>=75 tail applies.
interval_positive("J at x=84", J.subs(x, 84), left3, tail)
require(6-42*left3 < 0, "J decreases for x<=84 on third overlap")
quadrant_positive("F quadrant x>=84,s>=9/16", S.Integer(84), left3)
quadrant_positive("F quadrant x>=75,s>=3/5", S.Integer(75), tail)
require(coefficient_count == 107, "complete coefficient certificate count")

# Positive full-model and hostile controls, all kept in exact original types.
G = v**4-13*v**3+55*v*v-x*v+y
require(S.Poly(G.subs({x:84, y:35}), v).count_roots(0, S.oo) == 4, "nonvacuous zero-beta full geometry")
for name, numerator in (("C", numC), ("D", [1, -11, 36, -5*x/12, y/7])):
    seq = []
    for k in range(9):
        seq.append(S.expand((numerator[k] if k < 5 else 0)-sum(den[j]*seq[k-j] for j in range(1, min(k, 4)+1))))
    H = S.Matrix(5, 5, lambda i, j:seq[i+j].subs({x:84, y:35}))
    minors = [H[:k, :k].det() for k in range(1, 6)]
    for minor in minors:
        require(minor > 0, "full-model positive residue control")
    print(f"(84,35,0), {name}/B full Gram leading minors: {minors}.")
same(cap.subs(x, 75), 0, "C-only repeated hostile lies in packet")
hostile = {x:84, y:S.Rational(1050, 11), s:S.Rational(1, 6)}
same(p.subs(hostile), 0, "packet-deletion hostile uses original phase")
same(T.subs(hostile), S.Rational(228261457669, 313632), "positive response without C packet")
same(ordinary3.subs(hostile), -S.Rational(639, 11), "exact first failed packet condition")
off = {x:84, y:35, s:S.Rational(1, 4)}
same(T.subs(off), S.Rational(412444713751, 32768), "positive off-phase response at full-model point")
same(p.subs(off), S.Rational(335, 176), "original-phase firewall")

payload = {"response_F":str(F), "response_T":str(T), "phase_p":str(p),
           "necessary_J":str(J), "phase_endpoint_bounds":phase_bounds,
           "positive_coefficient_count":coefficient_count, "certificates":certificates}
path = DEST / (Path(__file__).stem+"_certificate.json")
path.write_text(json.dumps(payload, sort_keys=True, indent=2)+"\n", encoding="utf-8", newline="\n")
print("Coefficient certificate SHA256: "+hashlib.sha256(path.read_bytes()).hexdigest())
print("Entire zero-beta boundary is closed, even with C alone. Shared-root boundaries and general sign remain OPEN.")
print(f"PASS: {GATES} always-active exact gates.")
