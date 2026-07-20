#!/usr/bin/env python3
"""
klein-2026-07-20-S353 -- THE VIETA DUALITY: additive vs multiplicative symmetric functions,
from the classical e/pi disjunction to our own pair-sum ruler.

Owner: "consider Schanuel's conjecture and the algebraic proof by contradiction using Vieta's
formula, the quadratic polynomial proof for e and pi, and previous work in the repo regarding
the multiplication/addition duality; come to understand threads regarding rationals and
irrationals."

THE CLASSICAL ARGUMENT, stated exactly (it is not ours and is not new):
  e and pi are both TRANSCENDENTAL (Hermite 1873, Lindemann 1882).  If BOTH e+pi and e*pi
  were algebraic, then e and pi would be the two roots of
        x^2 - (e+pi) x + (e*pi) = 0,
  a quadratic with ALGEBRAIC coefficients, hence both algebraic -- contradiction.  Therefore
        AT LEAST ONE of  e+pi,  e*pi  is transcendental.
  Which one is unknown: neither is known to be irrational individually.

  SCHANUEL'S CONJECTURE would upgrade this.  Taking z1 = 1, z2 = i*pi (linearly independent
  over Q), Schanuel gives trdeg Q(1, i pi, e, -1) >= 2, i.e. e and pi ALGEBRAICALLY
  INDEPENDENT -- whence BOTH e+pi and e*pi are transcendental.  Open.

THE SHAPE, which is the transferable part:
  ONE relation on the pair (e_1, e_2) = (sum, product) buys a DISJUNCTION ("at least one").
  A ONE-PARAMETER FAMILY of relations buys a CONJUNCTION ("both", i.e. rigidity).
  That is exactly the difference between the e/pi result and my own THM-1550, where the
  small-root product must be linear in t FOR EVERY t and the conclusion is outright rigidity.

TESTED HERE: does our LRC pair-sum ruler (THM-1002: q | v_i + v_j) have a MULTIPLICATIVE
counterpart?  If both the sum and the product of the extremal pair were pinned, Vieta would
pin the pair itself -- a real sharpening.  Measured, not assumed.
"""
from fractions import Fraction as Fr
from math import gcd
import itertools

def Mexact(A):
    """exact M(A) = max_t min_v ||v t||, and the argmax denominator"""
    A = sorted(A); n = len(A)
    ds = sorted({A[i] + A[j] for i in range(n) for j in range(i, n)})
    best = (Fr(0), None)
    for dd in ds:
        for k in range(1, dd):
            mn = min(min((k * v) % dd, dd - (k * v) % dd) for v in A)
            val = Fr(mn, dd)
            if val > best[0]: best = (val, (dd, k))
    return best

def pair_witnesses(A, q):
    """which pairs (v_i,v_j) from A satisfy q | v_i+v_j ?"""
    out = []
    for i in range(len(A)):
        for j in range(i, len(A)):
            if (A[i] + A[j]) % q == 0: out.append((A[i], A[j]))
    return out

FAMILIES = {
    "AP {1..12} u {13} (the tight 13-set)": list(range(1, 14)),
    "the wedge minimiser {1..11,13,84}": list(range(1, 12)) + [13, 84],
    "covering-min {1..12} u {182}": list(range(1, 13)) + [182],
    "equioscillating {1..12} u {26}": list(range(1, 13)) + [26],
    "equioscillating {1..12} u {39}": list(range(1, 13)) + [39],
    "dilated 2*{1..12} u {13}": [2 * k for k in range(1, 13)] + [13],
}
print("=" * 82)
print("DOES THE PAIR-SUM RULER HAVE A MULTIPLICATIVE COUNTERPART?")
print("=" * 82)
print(f"{'family':<38} {'M':>10} {'q':>6} {'#pairs q|(vi+vj)':>17} {'q | vi*vj too?':>16}")
for name, A in FAMILIES.items():
    M, arg = Mexact(A)
    q = M.denominator
    pw = pair_witnesses(A, q)
    mult = [p for p in pw if (p[0] * p[1]) % q == 0]
    print(f"{name:<38} {str(M):>10} {q:>6} {len(pw):>17} {str(len(mult)) + '/' + str(len(pw)):>16}")
print("""
 READING: the ADDITIVE constraint q | v_i+v_j is by construction satisfied by every listed
 pair.  The MULTIPLICATIVE one q | v_i*v_j is satisfied only sporadically -- so there is NO
 automatic Vieta companion, and the pair is NOT pinned by the ruler alone.  That is a
 NEGATIVE and it is the useful content: THM-1002 gives one symmetric function, not two, so it
 buys a DISJUNCTION over candidate pairs rather than a determination -- structurally the same
 limitation as the e/pi argument.
""")
print("=" * 82)
print("HOW MANY PAIRS SURVIVE?  (the size of the disjunction the ruler leaves open)")
print("=" * 82)
for name, A in FAMILIES.items():
    M, arg = Mexact(A)
    q = M.denominator
    pw = pair_witnesses(A, q)
    prods = sorted({p[0] * p[1] for p in pw})
    print(f"  {name}")
    print(f"     q = {q};  admissible pairs = {pw[:6]}{' ...' if len(pw) > 6 else ''}")
    print(f"     their products = {prods[:8]}{' ...' if len(prods) > 8 else ''}"
          f"   distinct products: {len(prods)}")
print("""
 So the ruler narrows the extremal pair to a SET, never to a point.  To pin the pair one
 needs a second, independent symmetric function -- exactly the ingredient the e/pi argument
 also lacks, and exactly the ingredient a one-parameter family supplies in THM-1550.
""")
