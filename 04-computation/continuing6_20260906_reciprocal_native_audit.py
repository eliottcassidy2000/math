"""Independent standard-library referee for native reciprocal repair packets.

No producer imports. Recover packets from every literal integer triangle,
then compare with exhaustive modular-root scans in the permitted slope box.
Infinite affine templates are checked by degree-bounded interpolation and
coefficientwise native bounds, independently of the producer's convolution.
"""
from fractions import Fraction
from itertools import combinations
from math import comb, gcd, isqrt
from pathlib import Path
import argparse
import hashlib
import json
import sys

sys.stdout.reconfigure(newline="\n")
HERE = Path(__file__).resolve().parent
STEM = "continuing6_20260906_reciprocal_native"
parser = argparse.ArgumentParser()
parser.add_argument("--certificate", type=Path)
args = parser.parse_args()
path = args.certificate or HERE / (STEM+"_certificate.json")
if not path.exists():
    path = HERE.parent / "05-knowledge" / "results" / path.name
GATES = 0


def require(ok, label):
    global GATES
    if not ok:
        raise RuntimeError(label)
    GATES += 1


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


require(sha(HERE / (STEM+".py")) == "d0981913961389520200a8215b1c2f8611ea7090bd379ba2024ceb9a6b345c11", "frozen producer source")
require(sha(path) == "1da5e84836de4aee72e0636e792a43715c9cb85eac9e621203e215ac7e228380", "frozen certificate")
data = json.loads(path.read_text(encoding="utf-8"))


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p)+1))


def determinant(pts):
    (x0,y0),(x1,y1),(x2,y2) = pts
    return x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1)


def packet_ok(p, packet):
    a,b,h,k,l,c = packet
    return (a >= 1 and h >= 5 and a*h == p-1 and 2 <= k < l and k+l == h
            and 1 <= c < k and b >= 1 and b*k*l == h+c*p
            and gcd(a,b) == 1 and b*l <= p-1)


small = [p for p in range(3,44,2) if prime(p)]
primes = small+[61,97,113,197]
require([row["p"] for row in data["controls"]] == primes, "complete declared prime universe")
triangle_universe = 0
eligible_slopes = 0
accepted = 0
rows = {}
for rec in data["controls"]:
    p = rec["p"]
    # Modular exponentiation is independent of the producer's inverse recurrence.
    pts = [(0,1),(1,0)]+[(i,pow(i,p-2,p)) for i in range(2,p)]
    require(len({y for x,y in pts}) == p, "one point per column")
    require({(y,x) for x,y in pts} == set(pts), "literal transpose invariance")
    bad = set()
    for tri in combinations(pts,3):
        triangle_universe += 1
        if determinant(tri) == 0:
            bad.add(tuple(sorted(tri)))
    through_P = [tri for tri in bad if (0,1) in tri]
    through_Q = [tri for tri in bad if (1,0) in tri]
    require(set(through_P).isdisjoint(through_Q) and len(through_P)+len(through_Q) == len(bad), "every native event has exactly one moved point")
    require({tuple(sorted((y,x) for x,y in tri)) for tri in through_P} == set(through_Q), "disjoint transpose pairing")
    native = set()
    for tri in through_P:
        (rx,ry),(sx,sy) = sorted(pt for pt in tri if pt != (0,1))
        g = gcd(rx, ry-1)
        a, b = (ry-1)//g, rx//g
        require(rx % b == sx % b == 0, "native line spacing")
        k,l = rx//b,sx//b
        h = k+l
        require((b*k*l-h) % p == 0, "integer carry extracted from native line")
        c = (b*k*l-h)//p
        packet = a,b,h,k,l,c
        require(packet_ok(p,packet), "all primitive arithmetic conditions recovered")
        require((rx,ry) == (b*k,a*k+1) and (sx,sy) == (b*l,a*l+1), "native coordinate recovery")
        require(5*a <= p-1 and 1 <= b < 2*a, "sharp slope enumeration bounds")
        native.add(packet)
    require(len(native) == len(through_P), "packet-to-triple injection")
    require(native == {tuple(packet) for packet in rec["packets"]}, "all saved primitive packets from literal triangles")
    require(rec["triples"] == len(bad) == 2*len(native), "exact doubled defect count")

    # No quadratic formula, square-root choice, or producer packet loops.
    arithmetic = set()
    squares = {j*j % p for j in range(1,p)}
    for a in range(1,(p-1)//5+1):
        if (p-1) % a:
            continue
        h = (p-1)//a
        for b in range(1,2*a):
            if gcd(a,b) != 1:
                continue
            eligible_slopes += 1
            H = min(h-1,(p-1)//b)
            roots = [k for k in range(p) if (a*b*k*k+b*k-1) % p == 0]
            delta = b*(b+4*a) % p
            require(len(roots) == (1 if delta == 0 else 2 if delta in squares else 0), "exact discriminant classification")
            bounded = len(roots) == 2 and all(1 <= k <= H for k in roots)
            one_bounded = delta != 0 and delta in squares and any(h-H <= k <= H for k in roots)
            require(bounded == one_bounded, "one-root interval equivalence in standard lifts")
            if bounded:
                k,l = roots
                require(k+l == h and (b*k*l-h) % p == 0, "native unwrapped root sum and carry")
                packet = a,b,h,k,l,(b*k*l-h)//p
                require(packet_ok(p,packet), "bounded roots imply full primitive packet")
                arithmetic.add(packet)
    require(arithmetic == native, "complete independent slope/root enumeration")
    chi = 0 if p == 5 else 1 if 5 % p in squares else -1
    require(chi == rec["chi5"], "literal quadratic character of five")
    if p != 5:
        require((chi == 1) == (p % 5 in (1,4)), "finite special-reciprocity controls")
    if chi == 1:
        require(any(a == b == 1 for a,b,h,k,l,c in native), "character-five native obstruction")
    rows[p] = native,len(bad),chi
    accepted += len(native)
    print(f"p={p}: {len(bad)} native triples, {len(native)} primitive packets, chi5={chi}.")
require(triangle_universe == sum(comb(p,3) for p in primes), "literal complete triangle universe")
require(min(p for p in small if rows[p][1] and rows[p][2] == -1) == 37, "first nonresidue failure among all smaller odd primes")
require(all(rows[p][1] == 0 for p in (3,5,7,13,17,23)), "genuine successful controls")
require((12,11,5,2,3,1) in rows[61][0], "cofactor-five sharpness")
require((4,3,28,8,20,4) in rows[113][0] and rows[113][1] == 4, "primitive-gauge and extra-defect boundary")
require([k for k in range(13) if (2*k*k+k-1) % 13 == 0] == [7,12], "square-discriminant lift hostile")
require(determinant(((0,1),(7,2),(12,12))) == 65, "modular determinant need not vanish over integers")
require([k for k in range(5) if (k*k+k-1) % 5 == 0] == [2], "ramified character-five control")

# All-parameter template proof by interpolation of degree<=2 identities.
def lin(a, t):
    return a[0]+a[1]*t


for rec in data["templates"]:
    r0,step,h,k,l,c = (rec[key] for key in ("residue","modulus","h","k","l","c"))
    require((r0,step) in ((37,360),(43,70),(97,120)), "exact template list")
    require(gcd(r0,step) == 1, "Dirichlet coprimality hypothesis")
    require(step % 5 == 0 and r0 % 5 in (2,3), "fixed nonresidue progression")
    require(h == k+l and 2 <= k < l and 1 <= c < k, "template parameters")
    a = [Fraction(r0-1,h),Fraction(step,h)]
    b = [Fraction(h+c*r0,k*l),Fraction(c*step,k*l)]
    require(a == rec["a"] and b == rec["b"], "independent affine slope coefficients")
    for j,coords,quotient in zip((k,l),rec["points"],rec["product_quotients"]):
        xx,yy = coords
        require(xx == [j*b[0],j*b[1]] and yy == [j*a[0]+1,j*a[1]], "template coordinate formulas")
        for coord in (xx,yy):
            require(coord[0] >= 2 and coord[1] >= 0, "all-parameter coordinate lower bound")
            require(r0-1-coord[0] >= 0 and step-coord[1] >= 0, "all-parameter native upper bound")
        # The product and line identities have degree<=2. Three evaluations
        # prove them as polynomials, rather than sampling an infinite family.
        for t in (0,1,2):
            p = r0+step*t
            aa,bb = lin(a,t),lin(b,t)
            xp,yp = lin(xx,t),lin(yy,t)
            require(xp*yp-1 == p*lin(quotient,t), "degree-two reciprocal product identity")
            require(bb*(yp-1) == aa*xp, "degree-two anchored-line identity")
    # Strictly distinct points and positive slope persist for every t>=0.
    require(all(ai.denominator == 1 for ai in a+b), "integral affine construction")
    require(a[0] > 0 and a[1] >= 0 and b[0] > 0 and b[1] >= 0, "positive all-parameter slope")
    print(f"p={r0} mod {step}: all-parameter native identities and box inequalities pass; infinitude uses CITED Dirichlet.")

print(f"Independent universe: {triangle_universe} literal integer triples, {eligible_slopes} eligible primitive slopes, {accepted} accepted packets.")
print("Prime-family identities, native intervals, transpose multiplicity, and primitive gauge all pass.")
print("No successful infinite family, arbitrary-swap classification, or 2p-point construction inferred.")
print(f"PASS: {GATES} always-active exact gates.")
