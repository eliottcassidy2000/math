#!/usr/bin/env python3
"""elementary_abelian_escape_kps_S128c107.py -- kind-pasteur-2026-07-20-S128c107

DOES THE INDEX-1 COLLAPSE SURVIVE A BIGGER GROUP?

mac-mini's THM-1385: for a FREE Z/2 action on an ASPHERICAL space, ind = 1.  So on
the resonance k-torus Borsuk-Ulam collapses to a single equation, for every k and
every free involution, with an explicit zero-free odd map T^k -> R^2 as the
sharpness witness.

TWO QUESTIONS THIS SCRIPT SETTLES.

(A) CAN A SUBSPACE ESCAPE?  No, and for a soft reason: ind is MONOTONE under
    Z/2-maps, and if A is an invariant subspace with free action then the inclusion
    A -> T^k is a Z/2-map, so ind(A) <= ind(T^k) = 1.  Corollary: T^k contains NO
    invariant free S^m for m >= 2.  Nothing to compute -- recorded because it is the
    first escape anyone would try, and it is closed.

(B) CAN A BIGGER GROUP ESCAPE?  YES -- and this is the point.  THM-1385 is a
    statement about Z/2, i.e. about a SINGLE character.  Take G = (Z/2)^n acting on
    T^k by n half-translations h_1..h_n that are independent over F_2 (so n <= k);
    the action is still free, the quotient is still the aspherical T^k, and Gamma is
    still torsion-free -- mac-mini's hypotheses all hold.  But the relevant index is
    now the Fadell-Husseini ideal

        Index_G(T^k) = ker( F_2[x_1..x_n] -> H^*(T^k/G; F_2) = Lambda(e_1..e_k) ),
        x_i |-> w_i.

    The exterior algebra kills SQUARES, so x_i^2 is in the kernel -- that is exactly
    the Z/2 collapse.  But it does NOT kill PRODUCTS OF DISTINCT generators:
    w_1...w_n != 0 whenever the h_i are independent.  Since Index_G(S(V)) = (x_1...x_n)
    for V the sum of the n DISTINCT characters, there is no G-map T^k -> S(V), so

        every G-equivariant f : T^k -> V has a ZERO,   codimension n, not 1.

    The collapse is a property of the GROUP, not of the torus.  The sharpness witness
    (cos 2*pi*t_1, sin 2*pi*t_1) escapes precisely because BOTH coordinates carry the
    SAME character; distinct characters cannot be escaped.

WHAT IS VERIFIED HERE
  1. mac-mini's zero-free witness, reproduced at several k (confirming ind <= 1).
  2. The algebraic heart: w_1...w_n != 0 in Lambda(F_2^k) iff the h_i are independent
     -- an exact F_2 rank computation.
  3. The topological claim at n = k = 2: MANY random (Z/2)^2-equivariant maps
     T^2 -> R^2 with distinct characters, each certified to have a zero by an exact
     winding number around a grid cell.
  4. The CONTROL that makes (3) meaningful: the same random search with REPEATED
     characters, where zero-free maps must and do exist.  A test that cannot produce
     the escape when the escape is available proves nothing (MISTAKE-196).
"""
import sys
import math
import random
from itertools import combinations, product

random.seed(20260720)
NSAMP = int(sys.argv[1]) if len(sys.argv) > 1 else 400
GRID = int(sys.argv[2]) if len(sys.argv) > 2 else 120


# ---------------------------------------------------------------- 1. witness
print("=" * 78)
print("(1) mac-mini's SHARPNESS WITNESS reproduced: a zero-free ODD map T^k -> R^2")
print("=" * 78)
print("  f(t) = (cos 2*pi*t_1, sin 2*pi*t_1),  involution t -> t + h, h = (1/2,0,..,0)")
for k in (1, 2, 3, 5, 8):
    worst = 1e9
    odd_err = 0.0
    for _ in range(4000):
        t = [random.random() for _ in range(k)]
        f = (math.cos(2 * math.pi * t[0]), math.sin(2 * math.pi * t[0]))
        worst = min(worst, math.hypot(*f))
        t2 = list(t)
        t2[0] = (t2[0] + 0.5) % 1.0
        g = (math.cos(2 * math.pi * t2[0]), math.sin(2 * math.pi * t2[0]))
        odd_err = max(odd_err, abs(g[0] + f[0]), abs(g[1] + f[1]))
    print("  k = %-2d  min |f| = %.6f  (zero-free)   max oddness error = %.2e"
          % (k, worst, odd_err))
print("  -> confirms ind(T^k) <= 1 for the single Z/2: BOTH coordinates carry the")
print("     SAME character, so the target is really S^1 and an equivariant map exists.")


# ---------------------------------------------------------------- 2. algebra
def f2_rank(vecs, k):
    rows = [int(v) for v in vecs]
    r = 0
    piv = []
    for row in rows:
        cur = row
        for p in piv:
            cur = min(cur, cur ^ p)
        if cur:
            piv.append(cur)
            piv.sort(reverse=True)
            r += 1
    return r


print()
print("=" * 78)
print("(2) THE ALGEBRAIC HEART:  w_1...w_n != 0 in Lambda(F_2^k)  iff  h_i independent")
print("=" * 78)
print("  In an exterior algebra a product of degree-1 classes is non-zero exactly")
print("  when they are linearly independent.  So the (Z/2)^n index survives to")
print("  degree n = rank of the half-period subgroup, and dies at n = rank + 1.")
for k in (2, 3, 4):
    print("  k = %d :" % k)
    for n in range(1, k + 2):
        # try to pick n independent half-periods in (Z/2)^k
        best = 0
        for combo in combinations(range(1, 1 << k), n):
            best = max(best, f2_rank(combo, k))
            if best == n:
                break
        ok = (best == n)
        print("     n = %-2d  independent set exists : %-5s  ->  w_1..w_n %s"
              % (n, ok, "!= 0  (index reaches degree %d)" % n if ok else "= 0  (collapse)"))
print("  -> the reachable codimension is exactly k, the rank of the 2-torsion.")


# ---------------------------------------------------------------- 3. topology
def winding_cell(fs, x0, y0, hstep, steps=40):
    """Winding number of f = (f1,f2) around the boundary of a grid cell."""
    pts = []
    for i in range(steps):
        pts.append((x0 + hstep * i / steps, y0))
    for i in range(steps):
        pts.append((x0 + hstep, y0 + hstep * i / steps))
    for i in range(steps):
        pts.append((x0 + hstep - hstep * i / steps, y0 + hstep))
    for i in range(steps):
        pts.append((x0, y0 + hstep - hstep * i / steps))
    tot = 0.0
    prev = None
    for (x, y) in pts + [pts[0]]:
        a, b = fs(x, y)
        if a == 0 and b == 0:
            return None
        th = math.atan2(b, a)
        if prev is not None:
            d = th - prev
            while d > math.pi:
                d -= 2 * math.pi
            while d < -math.pi:
                d += 2 * math.pi
            tot += d
        prev = th
    return int(round(tot / (2 * math.pi)))


def make_map(chars, nterm=3):
    """Random equivariant f: T^2 -> R^2. chars[i] = (p_i, q_i) in {0,1}^2 says the
    parity of (m_1, m_2) forced on the Fourier support of f_i by the character."""
    comps = []
    for (p, q) in chars:
        terms = []
        for _ in range(nterm):
            m1 = 2 * random.randint(-2, 2) + p
            m2 = 2 * random.randint(-2, 2) + q
            terms.append((m1, m2, random.uniform(-1, 1), random.uniform(-1, 1)))
        comps.append(terms)

    def f(x, y):
        out = []
        for terms in comps:
            s = 0.0
            for (m1, m2, a, b) in terms:
                ang = 2 * math.pi * (m1 * x + m2 * y)
                s += a * math.cos(ang) + b * math.sin(ang)
            out.append(s)
        return tuple(out)
    return f


def has_zero(f, grid=GRID):
    h = 1.0 / grid
    for i in range(grid):
        for j in range(grid):
            w = winding_cell(f, i * h, j * h, h, steps=16)
            if w is None or w != 0:
                return True
    return False


print()
print("=" * 78)
print("(3) THE TOPOLOGICAL CLAIM at n = k = 2: DISTINCT characters force a zero")
print("=" * 78)
print("  h_1 = (1/2, 0), h_2 = (0, 1/2).  f_1 odd under h_1 / even under h_2")
print("  (Fourier support m_1 odd, m_2 even); f_2 the mirror.  Sampling %d maps."
      % NSAMP)
distinct = [(1, 0), (0, 1)]
nz = 0
for _ in range(NSAMP):
    f = make_map(distinct)
    if has_zero(f):
        nz += 1
print("  maps with a CERTIFIED zero (non-zero winding on some cell) : %d / %d"
      % (nz, NSAMP))
print("  -> %s" % ("ALL have zeros, as the index argument predicts"
                   if nz == NSAMP else "COUNTEREXAMPLE FOUND -- the claim is FALSE"))

print()
print("=" * 78)
print("(4) THE CONTROL -- can this search FIND a zero-free map when one exists?")
print("=" * 78)
print("  Same machinery, REPEATED character: f_1 and f_2 both odd under h_1 and even")
print("  under h_2.  Then the target is a single character's R^2 = S^1 and mac-mini's")
print("  witness applies, so zero-free maps MUST exist and the search must find them.")
same = [(1, 0), (1, 0)]
free = 0
for _ in range(NSAMP):
    f = make_map(same)
    if not has_zero(f):
        free += 1
print("  zero-free maps found with REPEATED characters : %d / %d" % (free, NSAMP))
print("  explicit witness (cos 2pi x, sin 2pi x) is zero-free : %s"
      % (not has_zero(lambda x, y: (math.cos(2 * math.pi * x),
                                    math.sin(2 * math.pi * x)))))
print()
if free > 0 and nz == NSAMP:
    print("  INSTRUMENT VALIDATED: the search DOES find zero-free maps when they exist")
    print("  (repeated characters) and finds NONE when the index argument forbids them")
    print("  (distinct characters).  The contrast is the result:")
    print("     Z/2, one character   -> ind = 1, codimension 1   (THM-1385)")
    print("     (Z/2)^n, n distinct  -> codimension n <= k")
    print("  The index-1 collapse is a property of the GROUP, not of the torus.")
elif free == 0:
    print("  INSTRUMENT NOT VALIDATED: the search failed to find any zero-free map even")
    print("  where one provably exists, so its negative results carry no weight.")
