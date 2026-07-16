# opus-2026-07-16-S332 -- HYP-7135 / THM-925 (B): T2 LAYER 2 -- THE ENGINE.
# err(A,B,C) := mu(W cap D_C) - (2/13) mu(W),  W = D_A cap D_B (coprime A<B).
# EXACT IDENTITY: err = sum over components [psi(right) - psi(left)] where
#   psi(t) = mu([0,t] cap D_C) - (2/13) t  (continuous, period 1/C).
# STRUCTURE (Layer 1, verified): components indexed by beat m = iB - jA;
#   per regime the LEFT endpoints form an exact AP with step B*/A or B~/B
#   (Bezout: B B* - A B~ = 1), and the RIGHT endpoints with the OTHER step.
# ENGINE: fold by period: psi(t) = psi1({C t}); each edge sequence becomes a
#   rotation orbit y_m = y0 + m*beta mod 1, beta = {C * delta_f} exact
#   rational; DENJOY-KOKSMA: |sum_{m<N} psi~(y_m)| <= (floor(N/q)+1) * V1
#   for any continued-fraction denominator q of beta, V1 = variation of psi
#   over one period (~ 4/(13 C)).  RESONANCE: q large iff the triple forms
#   |qC - rA - sB| stay large (e_A(q) = |q C B* mod+- A| etc).
# This script VERIFIES, exactly: (V1) the identity; (V2) per-edge-family
#   AP-ness + the DK inequality for every CF denominator; (V3) the assembled
#   T2 bound vs exact err on a grid incl. resonant triples; (V4) q* vs the
#   triple small form y3*.
from fractions import Fraction
from math import gcd
import random, bisect

F = Fraction

def teeth(x):
    w = F(1, 13*x); out = []
    for j in range(x):
        a, b = F(j, x) - w, F(j, x) + w
        if a < 0: out += [(F(0), b), (a + 1, F(1))]
        elif b > 1: out += [(a, F(1)), (F(0), b - 1)]
        else: out.append((a, b))
    out.sort()
    m = []
    for iv in out:
        if m and iv[0] <= m[-1][1]: m[-1] = (m[-1][0], max(m[-1][1], iv[1]))
        else: m.append(list(iv))
    return [(a, b) for a, b in m]

def isect(u, v):
    out, i, j = [], 0, 0
    while i < len(u) and j < len(v):
        a, b = max(u[i][0], v[j][0]), min(u[i][1], v[j][1])
        if a < b: out.append((a, b))
        if u[i][1] < v[j][1]: i += 1
        else: j += 1
    return out

def mu(iv): return sum(b - a for a, b in iv)

class Psi:
    """psi(t) = mu([0,t] cap D_C) - (2/13) t, exact; plus mean and per-period
    variation V1."""
    def __init__(self, C):
        self.C = C
        tt = teeth(C)
        self.starts = [a for a, b in tt]; self.ends = [b for a, b in tt]
        self.cum = [F(0)]
        for a, b in tt: self.cum.append(self.cum[-1] + (b - a))
        # mean of psi over [0,1): integrate piecewise-linear psi exactly
        # psi' = 1_{D_C} - 2/13; breakpoints = tooth edges
        pts = sorted({F(0), F(1)} | set(self.starts) | set(self.ends))
        tot = F(0); prev_t, prev_v = F(0), F(0)
        for t in pts[1:]:
            v = self.val(t)
            tot += (prev_v + v) * (t - prev_t) / 2
            prev_t, prev_v = t, v
        self.mean = tot
        # V1: variation over one period [-1/(13C) .. shifted]: psi rises on
        # teeth (slope 11/13), falls off-teeth (slope -2/13); one period has
        # one rise 2/(13C)*(11/13) and one fall of equal size
        self.V1 = 2 * F(2, 13*C) * F(11, 13)
    def val(self, t):
        t = t % 1
        k = bisect.bisect_right(self.starts, t) - 1
        s = F(0)
        if k >= 0:
            s = self.cum[k] + min(max(t - self.starts[k], F(0)),
                                  self.ends[k] - self.starts[k])
        return s - F(2, 13) * t
    def tilde(self, t): return self.val(t) - self.mean

def cf_denoms(fr):
    """continued-fraction convergent denominators of fr in (0,1)."""
    p, q = fr.numerator, fr.denominator
    ds = []; h0, h1 = 0, 1
    while q:
        a = p // q
        h0, h1 = h1, a * h1 + h0
        ds.append(h1)
        p, q = q, p - a * q
    return [d for d in ds if d >= 1]

def component_tags(W, A, B):
    """exact (m, left, right) per component; i,j recovered and verified."""
    tags = []
    wA, wB = F(1, 13*A), F(1, 13*B)
    for (a, b) in W:
        c = (a + b) / 2
        i = round(c * A)
        if not (F(i, A) - wA <= a and b <= F(i, A) + wA):
            for di in (-1, 1):
                if F(i+di, A) - wA <= a and b <= F(i+di, A) + wA: i += di; break
        j = round(c * B)
        if not (F(j, B) - wB <= a and b <= F(j, B) + wB):
            for dj in (-1, 1):
                if F(j+dj, B) - wB <= a and b <= F(j+dj, B) + wB: j += dj; break
        assert F(i, A) - wA <= a and b <= F(i, A) + wA, (A, B, a, b, i)
        assert F(j, B) - wB <= a and b <= F(j, B) + wB
        tags.append((i * B - j * A, a, b))
    tags.sort()
    return tags

def ap_families(tags):
    """split (m, left, right) list into maximal runs where BOTH the left and
    right endpoint sequences are exact APs (consecutive m)."""
    fams = []
    cur = [tags[0]]
    dl = dr = None
    for k in range(1, len(tags)):
        m0, l0, r0 = cur[-1]; m1, l1, r1 = tags[k]
        ndl, ndr = (l1 - l0) % 1, (r1 - r0) % 1
        if m1 == m0 + 1 and (dl is None or (ndl == dl and ndr == dr)):
            if dl is None and len(cur) == 1: dl, dr = ndl, ndr
            cur.append(tags[k])
        else:
            fams.append((cur, dl, dr)); cur = [tags[k]]; dl = dr = None
    fams.append((cur, dl, dr))
    return fams

def t2_bound_and_check(A, B, C, psi=None, verbose=False):
    """returns (exact err, assembled T2 bound, q*, per-family DK all-ok)."""
    W = isect(teeth(A), teeth(B))
    if psi is None: psi = Psi(C)
    WC = isect(W, teeth(C))
    err = mu(WC) - F(2, 13) * mu(W)
    # (V1) identity
    ident = sum(psi.val(b) - psi.val(a) for a, b in W)
    assert ident == err, (A, B, C)
    tags = component_tags(W, A, B)
    fams = ap_families(tags)
    bound = F(0); dk_ok = True; qstar = None
    for comps, dl, dr in fams:
        N = len(comps)
        if N <= 2:
            bound += N * psi.V1  # junction: |psi(r)-psi(l)| <= V1/2 each*2
            continue
        for (edge_idx, dstep) in ((1, dl), (2, dr)):
            beta = (C * dstep) % 1
            S = sum(psi.tilde(c[edge_idx]) for c in comps)
            if beta == 0:
                q = 1
            else:
                qs = [q for q in cf_denoms(beta) if q <= N] or [1]
                q = max(qs)
            dkb = (N // q + 1) * psi.V1
            if abs(S) > dkb:
                dk_ok = False
                if verbose:
                    print(f"    DK FAIL A={A} B={B} C={C} fam N={N} "
                          f"edge{edge_idx} q={q}: |S|={float(abs(S)):.3e} "
                          f"> {float(dkb):.3e}")
            bound += dkb
            qstar = q if qstar is None else min(qstar, q)
    return err, bound, (qstar or 1), dk_ok

print("=" * 72)
print("T2 ENGINE VERIFICATION (identity exact-asserted on every triple)")
print("=" * 72)
random.seed(3320)
worst = []
n_ok = n_tot = 0
dk_all = True
rows = []
for (A, B) in [(19, 24), (50, 63), (87, 101), (55, 89), (120, 143), (73, 90)]:
    if gcd(A, B) != 1: continue
    Cs = {A + B, A + B + 1, A + B - 1, 2*A + B, A + 2*B, 2*(A + B) + 1,
          3*A + 2*B + 2, 13*A + 7}
    while len(Cs) < 18:
        Cs.add(random.randint(B + 1, 13 * A))
    for C in sorted(Cs):
        if C in (A, B) or C <= B: continue
        err, bnd, qst, ok = t2_bound_and_check(A, B, C, verbose=True)
        dk_all &= ok
        n_tot += 1
        if abs(err) <= bnd: n_ok += 1
        y3 = min(abs(q*C - r*A - s*B) for q in range(1, 14)
                 for r in range(-13, 14) for s in range(-13, 14))
        rows.append((A, B, C, float(err), float(bnd), qst, y3))
print(f"triples: {n_tot}; identity: all exact; DK inequality all families: "
      f"{dk_all}; assembled bound covers err: {n_ok}/{n_tot}")
print()
print("q* vs y3* vs err/bound (sample, sorted by q*):")
rows.sort(key=lambda r: r[5])
for r in rows[:8] + rows[-8:]:
    A, B, C, err, bnd, qst, y3 = r
    print(f"  A={A:3d} B={B:3d} C={C:4d}: err={err:+.3e} bound={bnd:.3e} "
          f"ratio={abs(err)/bnd if bnd else 0:5.2f}  q*={qst:4d}  y3*={y3:3d}")
print()
# correlation summary
import statistics
lo = [r for r in rows if r[6] == 0]
hi = [r for r in rows if r[6] >= 4]
print(f"resonant (y3*=0):   n={len(lo)}, max|err|={max(abs(r[3]) for r in lo):.3e}, "
      f"median q*={statistics.median(r[5] for r in lo)}")
if hi:
    print(f"clean (y3*>=4):     n={len(hi)}, max|err|={max(abs(r[3]) for r in hi):.3e}, "
          f"median q*={statistics.median(r[5] for r in hi)}")
mid = [r for r in rows if 1 <= r[6] <= 3]
if mid:
    print(f"near-res (1<=y3*<=3): n={len(mid)}, max|err|={max(abs(r[3]) for r in mid):.3e}, "
          f"median q*={statistics.median(r[5] for r in mid)}")
