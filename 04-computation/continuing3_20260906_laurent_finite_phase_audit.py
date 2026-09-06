"""Independent exact referee: literal carriers, corner Sturm, inverse differences.

Does not import or execute the producer. Requires Python 3 and SymPy.
The inherited e4>1/100 proof is an explicitly external analytical dependency.
"""
from pathlib import Path
from fractions import Fraction
from itertools import combinations, product
from math import comb
import hashlib
import json
import sys
import sympy as S

sys.stdout.reconfigure(newline="\n")
HERE = Path(__file__).resolve().parent
BASE = "continuing3_20260906_laurent_finite_phase"
PINS = {
    ".py": "8c49a16864d2b17c3fc888df0313e241b6fd5cdfb47c55a06920d709902269d5",
    ".out": "85de04cb1e32eae863d12c321e0defe35fddd68ddc19182164f2400c27ad3133",
    "_certificate.json": "afe514688067ab0fdc38f33ca650fbecdb0ddb40e08466b6395c7c5c29659385",
}
COUNT = 0


def need(ok, description):
    global COUNT
    if not bool(ok):
        raise ArithmeticError(description)
    COUNT += 1


def zero(expr, description):
    need(S.expand(expr) == 0, description)


def existing(suffix):
    name = BASE + suffix
    candidate = HERE / name
    if candidate.exists():
        return candidate
    if suffix == ".out":
        return HERE.parent / "05-knowledge/results" / name
    return candidate


for suffix, digest in PINS.items():
    need(hashlib.sha256(existing(suffix).read_bytes()).hexdigest() == digest,
         "frozen producer byte pin " + suffix)
cert = json.loads(existing("_certificate.json").read_text(encoding="utf-8"))

x, y, z, s, t, inv = S.symbols("x y z s t inv")
GB = 1 + 13*t + 55*t**2 + x*t**3 + y*t**4 + z*t**5
GC = 1 + 12*t + 45*t**2 + S.Rational(2, 3)*x*t**3 + S.Rational(3, 7)*y*t**4
GD = 1 + 11*t + 36*t**2 + S.Rational(5, 12)*x*t**3 + S.Rational(1, 7)*y*t**4
odd = sum(S.binomial(14, 2*j + 1)*t**j for j in range(7))
even = sum(S.binomial(14, 2*j)*t**j for j in range(8))
# Ordinary shifts before any coefficient multiplication:
# alpha=t^-1*(t O^2+E^2), beta-packet=t^-2*(GB^2+2t GC GD).
alpha = S.Poly(t*odd**2 + even**2, t)
packet = S.Poly(GB**2 + 2*t*GC*GD, t)
q = {}
for j in range(-1, 14):
    need(alpha.nth(j + 1) == S.binomial(28, 2*j + 2),
         "literal even-odd square multiplier")
    if j <= 8:
        q[j] = S.expand(alpha.nth(j + 1)*packet.nth(j + 2))
need(alpha.nth(-1) == 0, "unavailable carry has zero alpha multiplier")
need(q[-1] == 28, "lower carry exactly 28/t")
zero(q[8] - S.binomial(28, 18)*z**2, "top square coefficient")
zero(q[7] - S.binomial(28, 16)*(2*y*z + S.Rational(6, 49)*y**2),
     "next square-plus-crossing coefficient")
F = S.expand(sum(q[j]*S.Integer(-1)**j*s**(j + 1) for j in q))
first = S.expand(sum(S.Poly(odd, t).nth(j)*S.Poly(GB, t).nth(j+1)*
                    (-s)**j for j in range(5)))
zero(first - (182 - 20020*s + 2002*x*s**2 - 3432*y*s**3 + 2002*z*s**4),
     "literal first polynomial")
zero(first.subs({x:84, y:35, z:1}) - 2002*(s**4 - 60*s**3 + 84*s**2 - 10*s + S.Rational(1, 11)),
     "inherited actual endpoint-27 monic first row normalization")


def unpack(rows, variables):
    return S.expand(sum(S.Rational(value)*S.prod(v**k for v, k in zip(variables, exponent))
                        for exponent, value in rows))


zero(F - unpack(cert["full_T"], (x,y,z,s)), "entire certificate is original sQ(-s)")
for monomial, coefficient in S.Poly(F, x,y,z,s).terms():
    need(coefficient*S.Integer(-1)**(monomial[3] - 1) > 0,
         "all raw carried coefficients are positive")

# Verify the elementary root bounds in a different algebra system.
roots = S.symbols("a0:5")


def es(seq, k):
    return sum(S.prod(packet) for packet in combinations(seq, k))


e1, e2, e3, e4 = (es(roots, k) for k in range(1, 5))
zero(e1**3 - 3*e1*sum(a*a for a in roots) + 2*sum(a**3 for a in roots) - 6*e3,
     "formal third-power identity")
sos1 = sum((roots[i]-roots[j])**2 * es([roots[k] for k in range(5) if k not in (i,j)], 2)
           for i,j in combinations(range(5), 2))
sos2 = sum((a*b-c*d)**2 + (a*c-b*d)**2 + (a*d-b*c)**2
           for a,b,c,d in combinations(roots,4))
zero(e2**2 - 20*e4 - sos1 - sos2/3, "fourth-coefficient sum of squares")
zero(e1*e3 - 4*e4 - sum(a*a*es([roots[j] for j in range(5) if j != i], 2)
                        for i,a in enumerate(roots)), "relative fourth/third bound")
need(13**2 > 2*59, "at least three strictly positive roots")
need(5*S.Rational(71,10)**2 - 26*S.Rational(71,10) - 67 > 0,
     "maximum root less than 71/10")
need((S.Rational(71,10)*59-52)/3 < 123 < 130, "third-coefficient bound")
need(S.Rational(55**2,20) < 152, "fourth-coefficient bound")
need(S.Rational(55,10)**5 < 72**2, "fifth-coefficient pair-product bound")
need(2002-S.Rational(11154,110) > 0, "first equation positive before lower endpoint")
upper = first.subs({x:130,y:0,z:72,s:S.Rational(1,90)})
slope = -20020+S.Rational(4004*130,90)+S.Rational(8008*72,90**3)
need(upper < 0 and slope < 0, "unique simple first branch upper endpoint")

# An independent full-box sign proof: coordinatewise convexity reduces
# the whole coefficient box to eight univariate corner polynomials.
expected_quadratic = [198835*s**5*(66-85*s),
                      S.Rational(2005830,7)*s**7*(140-13*s),
                      13123110*s**9]
for variable, expected in zip((x,y,z), expected_quadratic):
    zero(S.diff(F,variable,2)/2 - expected, "separate-convexity coefficient")
need(66-S.Rational(85,80) > 0 and 140-S.Rational(13,80) > 0,
     "strict separate convexity on whole phase interval")


def changes(values):
    signs = [S.sign(v) for v in values if v != 0]
    return sum(a != b for a,b in zip(signs, signs[1:]))


corner_records = []
low, high, middle = S.Rational(1,120), S.Rational(1,80), S.Rational(1,100)
for xx, yy, zz in product((0,130),(0,152),(0,72)):
    pol = S.Poly(F.subs({x:xx,y:yy,z:zz})+400, s, domain=S.QQ)
    chain = pol.sturm()
    left = changes([member.eval(low) for member in chain])
    right = changes([member.eval(high) for member in chain])
    need(left == right, "corner has zero closed-interval Sturm roots")
    need(pol.eval(low) < 0 and pol.eval(high) < 0 and pol.eval(middle) < 0,
         "corner sign is strictly negative with margin 400")
    need(pol.count_roots(low,high) == 0, "independent library root-count agreement")
    corner_records.append([xx,yy,zz,int(left),int(right)])
print("EIGHT CORNER STURM RECORDS", corner_records)
print("Separate convexity proves sQ(-s)<-400 throughout the complete coefficient/phase box.")

# Validate every frozen Bernstein coefficient independently by an inverse
# finite-difference transform, rather than producer forward summation.
degrees = (2,2,2,9)
bank = {tuple(exponent): S.Rational(value) for exponent,value in cert["bernstein"]}
universe = set(product(*(range(d+1) for d in degrees)))
need(set(bank) == universe and len(bank) == 270, "complete tensor universe")
need(max(bank.values()) == -S.Rational(22645374245632441,52254720000000),
     "published Bernstein maximum")
X,Y,Z,U = S.symbols("X Y Z U")
cube = S.Poly(F.subs({x:130*X,y:152*Y,z:72*Z,s:low+(high-low)*U}, simultaneous=True), X,Y,Z,U)
for alpha_index in sorted(universe):
    need(bank[alpha_index] < -400, "each frozen Bernstein sign")
    difference = sum(S.Integer(-1)**sum(a-b for a,b in zip(alpha_index,beta_index)) *
                     S.prod(S.binomial(a,b) for a,b in zip(alpha_index,beta_index))*bank[beta_index]
                     for beta_index in product(*(range(a+1) for a in alpha_index)))
    rebuilt = S.prod(S.binomial(d,a) for d,a in zip(degrees,alpha_index))*difference
    need(rebuilt == cube.coeff_monomial(S.prod(v**a for v,a in zip((X,Y,Z,U),alpha_index))),
         "inverse differences reconstruct every power coefficient")

# Exact symbolic elimination from the independently reconstructed raw row.
substitute_z = S.Rational(12,7)*y*inv-x*inv**2+10*inv**3-inv**4/11
R = S.Poly(S.cancel((F/s**8).subs(s,1/inv).subs(z,substitute_z)), x,y,inv)
zero(R.as_expr()-unpack(cert["eliminated_Q_over_s7"], (x,y,inv)),
     "all eliminated coefficients agree with original-root substitution")
zero(first.subs({s:1/inv,z:substitute_z}, simultaneous=True),
     "elimination lies exactly on original first root")
envelope = {int(k):S.Rational(value) for k,value in cert["tail_envelope"]}
rebuilt = {}
for (a,b,k), coefficient in R.terms():
    if coefficient > 0:
        need(a <= 1 and b <= 1 and k > 0, "positive tail monomial tariff scope")
        rebuilt[k] = rebuilt.get(k,0) + coefficient*123**a*100**(2-b)
need(rebuilt == envelope, "all and only positive eliminated monomials retained")
zero(R.nth(0,2,0)+S.Rational(26075790,7), "strict negative retained constant")
tail = -S.Rational(26075790,7)+sum(v/S.Integer(75000)**k for k,v in envelope.items())
need(tail == S.Rational(cert["tail_at_75000"]) and tail < -120000,
     "entire tail envelope at exact cutoff")
need(all(v>0 and k>0 for k,v in envelope.items()), "monotone tail envelope")
hessian = S.det(S.hessian(R.as_expr(),(x,y)))
zero(hessian + S.Rational(1031232600,49)*inv**2*(144097056*inv**2+72692884*inv+10966105),
     "quoted indefinite Hessian factor")
print("TAIL AT 75000", tail)

# Literal exact valid-shape control and two deliberately invalid relaxations.
actual = [S.Rational(v,5) for v in (1,3,9,22,30)]
ev = [es(actual,k) for k in range(1,6)]
need(ev[:2] == [13,55], "valid model root anchors")
xp,yp,zp = ev[2:]
Cpoly = s**4-12*s**3+45*s**2-S.Rational(2,3)*xp*s+S.Rational(3,7)*yp
Dpoly = s**4-11*s**3+36*s**2-S.Rational(5,12)*xp*s+yp/7
for index, root in enumerate(actual):
    need((-1)**index*Cpoly.subs(s,root)>0 and (-1)**index*Dpoly.subs(s,root)>0,
         "both strict model interlacings")
need(F.subs({x:xp,y:yp,z:zp,s:75000})>0 and first.subs({x:xp,y:yp,z:zp,s:75000}) != 0,
     "large phase without original-root predicate is hostile")
need(F.subs({x:130,y:0,z:0,s:S.Rational(1,60)}) == S.Rational(34734566083,93312000)>0,
     "larger coefficient-box phase hostile")

hx,hy,hz,hs = S.Integer(104),S.Integer(50),S.Rational(37435088,3898125),S.Rational(15,2)
hostile_first = S.Poly(first.subs({x:hx,y:hy,z:hz}),s)
zero(hostile_first.as_expr() - S.Rational(26,50625)*(2*s-15)*
     (18717544*s**3-26680920*s**2+2595600*s-23625), "hostile exact first factorization")
need(hostile_first.eval(hs)==0, "hostile is on original root")
need(hostile_first.count_roots(0,S.oo)==4 and hostile_first.discriminant()!=0,
     "hostile first row has four distinct positive roots")
for lo,hi in ((S.Rational(1,99),S.Rational(1,98)),(S.Rational(3,32),S.Rational(5,53)),
              (S.Rational(70,53),S.Rational(37,28))):
    need(hostile_first.eval(lo)*hostile_first.eval(hi)<0, "three remaining phase brackets")
hostile_value = F.subs({x:hx,y:hy,z:hz,s:hs})/hs
need(hostile_value == S.Rational(78541969368658673,18480)>0, "exact positive hostile response")
margins = [13**2-S.Rational(5,2)*55,55**2-2*13*hx,hx**2-2*55*hy,hy**2-S.Rational(5,2)*hx*hz]
need(margins == [S.Rational(63,2),321,5316,S.Rational(2437924,779625)] and min(margins)>0,
     "all four Newton inequalities remain strict")
beta = S.Poly(s**5-13*s**4+55*s**3-hx*s**2+hy*s-hz,s)
need(beta.count_roots(-S.oo,S.oo)==beta.count_roots(0,S.oo)==1,
     "hostile beta has exactly one real root")
need(beta.count_roots(-S.oo,0)==0 and beta.eval(0)!=0,
     "no nonpositive beta root restores model membership")
print("NEWTON/REAL-FIRST-ROW HOSTILE Q", hostile_value)
print("PASS", COUNT, "always-active exact audit gates; producer never imported or executed.")
print("Accepted scope: first branch without interlacers; tail requires both inherited interlacers and the original root.")
