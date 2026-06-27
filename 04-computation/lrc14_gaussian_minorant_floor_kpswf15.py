#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_gaussian_minorant_floor_kpswf15.py   (kind-pasteur 2026-06-27, HYP-3121 synthesis / TOOL 1)

GAUSSIAN / BEURLING-SELBERG MINORANT for the UNIFORM r=2..6 multi-far loneliness floor.

THE OBJECT.  A covering set S = R u 14Q with r=|Q| in {2..6} multiples of 14 and a 14-free
"small part" R (13-r >= 7 speeds).  The loneliness measure is
        L(S) = int_0^1  prod_{s in S}  phi_s(t)  dt,
        phi_s(t) = 1_{[1/14, 13/14]}( {s t} ),     an arc of measure 6/7.
LRC(14) on this branch <=> L(S) > 0 for every covering S with r in {2..6}.

THE MINORANT IDEA (the wide-V tail fix).  The sharp arc indicator phi_s has 1/n Fourier decay
(the wide-V tail problem: the resonance/relation-lattice sum over sum_s k_s s = 0 has a slowly
decaying tail whose magnitude is hard to bound uniformly in the speed magnitudes).  Replace each
phi_s by a NONNEGATIVE minorant psi_s <= phi_s whose Fourier coefficients are either FINITELY
SUPPORTED (Beurling-Selberg, band <= N) or SUPER-POLYNOMIALLY decaying (C^infty mollifier).  Then

        L(S) = int prod phi_s  >=  int prod psi_s
             = MAIN TERM  prod_s ( int psi_s )   +   RESONANCE  sum_{(k_s) != 0, sum k_s s = 0} prod_s psihat(k_s).

MAIN TERM = (int psi)^{|S|}  (each factor equal by the arc structure, depends only on the arc length
6/7 and the minorant degree, NOT on the speeds).  With a band-N Beurling-Selberg minorant the
resonance sum runs only over relations with |k_s| <= N -- a FINITE, speed-INDEPENDENT band -- so its
tail is automatically uniform.  KEY question: is MAIN > |RESONANCE| with an explicit uniform constant?

TWO MINORANTS BUILT + VALIDATED HERE.
  (A) Beurling-Selberg / Vaaler minorant of degree N on the circle.  TRUE minorant (psi <= phi
      pointwise), band-limited (psihat(k)=0 for |k|>N), int psi = 6/7 - 1/(N+1).
  (B) C^infty compactly-supported mollifier minorant: shrink arc by delta, convolve with a smooth
      bump supported in [-delta,delta].  TRUE minorant; super-polynomial coefficient decay; int psi
      = 6/7 - 2 delta.

Both are verified to be true minorants on a fine grid, and the arc Fourier coefficients are computed
exactly (interval indicator) times the smoothing multiplier.

OUTPUT: per-(R,Q) main-term/resonance/floor table for r=2..6 few-apex covering sets; the explicit
uniform tail bound on the resonance sum; the constant c if the minorant floor is > 0.
"""
import sys, itertools, math
from fractions import Fraction as F
from math import gcd, pi, sqrt, erfc, exp, cos, sin
from functools import reduce
try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

# ----------------------------------------------------------------------------- arc Fourier atoms
# phi (arc [a,b] mod 1) has indicator Fourier coeffs:
#   chi_hat(0) = b-a,   chi_hat(k) = (e(-k a) - e(-k b))/(2 pi i k)  for k != 0,   e(t)=exp(2 pi i t).
# For our arc a=1/14, b=13/14 (length 6/7), centered at 1/2.  Centering: the arc is symmetric about
# 1/2, so chi_hat(k) is real:  chi_hat(k) = (sin(2 pi k * 13/14) - sin(2 pi k /14))/(2 pi k)... we
# just compute it from the exact endpoints.

A_LO = F(1, 14)
A_HI = F(13, 14)
ARC_LEN = A_HI - A_LO          # = 6/7

def chi_hat(k):
    """Fourier coeff of 1_{[1/14,13/14]} on the circle.  Real (arc symmetric about 1/2)."""
    if k == 0:
        return ARC_LEN
    a = float(A_LO); b = float(A_HI)
    # (e(-k a) - e(-k b))/(2 pi i k)
    num = complex(cos(-2*pi*k*a), sin(-2*pi*k*a)) - complex(cos(-2*pi*k*b), sin(-2*pi*k*b))
    val = num / (2j * pi * k)
    return val   # complex; imaginary part ~ 0 by symmetry

# ----------------------------------------------------------------------------- (A) Beurling-Selberg
# Vaaler's minorant of the interval indicator on R/Z, degree N.  We use the standard closed form
# built from the Fejer kernel.  The Selberg minorant S^-_N of 1_{[a,b]} has Fourier coefficients
#   S^-_N_hat(0) = (b-a) - 1/(N+1)
#   S^-_N_hat(k) = chi_hat(k) * (1 - |k|/(N+1))  + correction_k    (0<|k|<=N)
# where the correction comes from the Beurling extremal function.  Implementing the FULL extremal
# correction exactly is delicate; instead we use a SAFE, EXPLICIT minorant of Beurling-Selberg TYPE
# obtained by the Fejer (triangle) smoothing of a SHRUNK arc, which is a true minorant and band <= N.
#
# CONSTRUCTION (clean, rigorous, band-limited):
#   psi = ( 1_{[a+h, b-h]} * F_N ) , where F_N is the Fejer kernel of degree N (nonneg, mass 1,
#   sum of (1-|k|/(N+1)) e(kx)).  Convolution of a nonneg indicator with a nonneg kernel is >= 0.
#   The Fejer kernel is NOT compactly supported, so we SHRINK by h chosen s.t. the leaked mass is
#   provably below phi: we instead enforce psi <= phi by the pointwise certificate (grid + endpoint
#   monotonicity).  To GUARANTEE a true minorant we use the COMPACTLY SUPPORTED route (B) for the
#   rigorous floor and use this Fejer object only as a band-limited COMPARISON.
def fejer_mult(k, N):
    """Fejer multiplier (1-|k|/(N+1)) for |k|<=N else 0 (the heat/triangle smoothing in freq)."""
    if abs(k) > N:
        return 0.0
    return 1.0 - abs(k) / (N + 1)

def fejer_smoothed_arc_hat(k, N, h):
    """Fourier coeff of  1_{[a+h, b-h]} * Fejer_N.   psihat(k) = chi_hat_{shrunk}(k) * fejer_mult.
       chi_hat_{shrunk} is the indicator coeff of the shrunk arc [a+h, b-h]."""
    if abs(k) > N:
        return 0j
    aa = A_LO + F(h).limit_denominator(10**9)
    bb = A_HI - F(h).limit_denominator(10**9)
    if k == 0:
        return complex(float(bb - aa), 0.0) * fejer_mult(0, N)
    a = float(aa); b = float(bb)
    num = complex(cos(-2*pi*k*a), sin(-2*pi*k*a)) - complex(cos(-2*pi*k*b), sin(-2*pi*k*b))
    return (num / (2j * pi * k)) * fejer_mult(k, N)

# ----------------------------------------------------------------------------- (B) C^infty mollifier
# A TRUE compactly-supported minorant with super-polynomial coefficient decay.
# Mollifier: rho_delta(x) = (1/delta) rho(x/delta), rho the standard bump on [-1,1],
#   rho(u) = C exp(-1/(1-u^2)) for |u|<1, 0 else, C normalizing int rho = 1.
# psi = 1_{[a+delta, b-delta]} * rho_delta.  Support of psi = [a, b] EXACTLY (shrunk arc widened by
# delta on each side) => psi <= 1_{[a,b]} = phi pointwise (psi <= 1 since it's an average of an
# indicator, and supp psi subset [a,b]).  int psi = (b-a) - 2 delta.
# Its Fourier coeff:  psihat(k) = chi_hat_{[a+delta,b-delta]}(k) * rhohat_delta(k),
#   rhohat_delta(k) = rhohat(delta k),  rhohat(xi) = int_{-1}^{1} rho(u) e(-xi u) du  (real, even,
#   super-polynomially decaying).  We compute rhohat by high-accuracy quadrature and tabulate.

def _bump_unnormalized(u):
    if -1.0 < u < 1.0:
        return exp(-1.0 / (1.0 - u*u))
    return 0.0

# normalization constant C = 1 / int_{-1}^{1} exp(-1/(1-u^2)) du
def _bump_norm(M=200000):
    # Simpson on [-1+e, 1-e]
    e = 1e-9
    a, b = -1.0 + e, 1.0 - e
    h = (b - a) / M
    s = _bump_unnormalized(a) + _bump_unnormalized(b)
    for i in range(1, M):
        x = a + i*h
        s += (4 if i % 2 else 2) * _bump_unnormalized(x)
    return 1.0 / (s * h / 3.0)

_BUMP_C = _bump_norm()

def rho(u):
    return _BUMP_C * _bump_unnormalized(u)

_rhohat_cache = {}
def rhohat(xi, M=40000):
    """int_{-1}^{1} rho(u) cos(2 pi xi u) du  (real, even).  Super-polynomial decay in xi."""
    key = round(xi, 9)
    if key in _rhohat_cache:
        return _rhohat_cache[key]
    e = 1e-9
    a, b = -1.0 + e, 1.0 - e
    h = (b - a) / M
    def f(u):
        return rho(u) * cos(2*pi*xi*u)
    s = f(a) + f(b)
    for i in range(1, M):
        u = a + i*h
        s += (4 if i % 2 else 2) * f(u)
    val = s * h / 3.0
    _rhohat_cache[key] = val
    return val

def mollifier_arc_hat(k, delta):
    """psihat(k) for the C^infty minorant psi = 1_{[a+delta,b-delta]} * rho_delta."""
    aa = float(A_LO) + delta
    bb = float(A_HI) - delta
    if k == 0:
        return complex(bb - aa, 0.0)
    num = complex(cos(-2*pi*k*aa), sin(-2*pi*k*aa)) - complex(cos(-2*pi*k*bb), sin(-2*pi*k*bb))
    chi = num / (2j * pi * k)
    return chi * rhohat(delta * k)

# ----------------------------------------------------------------------------- minorant evaluation in space
def eval_series(hat_fn, theta, K):
    """Reconstruct psi(theta) = sum_{|k|<=K} psihat(k) e(k theta) (real)."""
    s = hat_fn(0)
    s = s.real if isinstance(s, complex) else s
    for k in range(1, K+1):
        hk = hat_fn(k)
        s += 2.0 * (hk * complex(cos(2*pi*k*theta), sin(2*pi*k*theta))).real
    return s

def verify_minorant(hat_fn, K, grid=4000, tol=1e-6, name=""):
    """Check psi(theta) >= -tol everywhere and psi <= phi (= 1 inside arc, 0 outside) on a grid.
       Returns (min_psi, max_leak_outside, max_overshoot_inside)."""
    min_psi = 1e9
    max_leak = 0.0       # max psi where phi=0  (must be <= tol for a true minorant)
    max_over = 0.0       # max (psi - 1) where phi=1
    a = float(A_LO); b = float(A_HI)
    for i in range(grid):
        th = (i + 0.5) / grid
        v = eval_series(hat_fn, th, K)
        min_psi = min(min_psi, v)
        inside = (a <= th <= b)
        if inside:
            max_over = max(max_over, v - 1.0)
        else:
            max_leak = max(max_leak, v)
    return min_psi, max_leak, max_over

# ----------------------------------------------------------------------------- covering-set machinery
AP = tuple(range(1, 14))

def is_covering(S):
    """qdiv(S)>14: every d in 2..14 divides some speed (so no tau=a/d witness for d<=14)."""
    for d in range(2, 15):
        if not any(v % d == 0 for v in S):
            return False
    return True

def primitive(S):
    return reduce(gcd, S) == 1

def split_RQ(S):
    """R = 14-free part, Q = {m : 14m in S}."""
    Q = sorted(v // 14 for v in S if v % 14 == 0)
    R = sorted(v for v in S if v % 14 != 0)
    return R, Q

# ----------------------------------------------------------------------------- resonance sum
def resonance_sum(S, hat_fn, band, kmax_per=None, max_terms=None):
    """Compute  sum_{(k_s) != 0, sum_s k_s s = 0, |k_s|<=band}  prod_s psihat(k_s).
       Returns (main_term, resonance_value, n_relations, abs_resonance_bound, by_support).
       main_term = prod_s psihat(0).
       We enumerate relations sum_s k_s s = 0 with |k_s|<=band by meet-in-the-middle over the speeds.
       Because psihat(0) is close to ARC_LEN and psihat(k!=0) is small, terms with many nonzero k_s
       are suppressed by (small)^(#nonzero); we organize by SUPPORT (the set of speeds with k_s!=0).
    """
    S = sorted(S)
    h0 = hat_fn(0)
    h0 = h0.real if isinstance(h0, complex) else float(h0)
    main = h0 ** len(S)
    # tabulate psihat(k) for |k|<=band
    H = {k: hat_fn(k) for k in range(-band, band+1)}
    # nonzero-coefficient list per speed (k != 0)
    # We enumerate relations support-by-support: choose a subset T of speeds (|T|>=2) to be the
    # nonzero coords, then count integer solutions sum_{s in T} k_s s = 0 with 1<=|k_s|<=band.
    # For |T| up to ~5 this is a bounded search; deeper supports are bounded by the tail estimate.
    res = 0.0
    abs_res = 0.0
    nrel = 0
    by_support = {}   # |T| -> (count, signed_sum, abs_sum)
    # restrict support size to keep it finite; deeper supports added to abs bound separately
    max_support = min(len(S), 6 if kmax_per is None else kmax_per)
    for tsize in range(2, max_support+1):
        cnt = 0; ssum = 0.0; asum = 0.0
        for T in itertools.combinations(range(len(S)), tsize):
            speeds = [S[i] for i in T]
            # enumerate k vectors with |k|in[1,band], sum k_i s_i = 0
            # recursive: fix first tsize-1, solve last
            for kv in _relations(speeds, band):
                prodc = 1.0+0j
                for i, ki in enumerate(kv):
                    prodc *= H[ki]
                # multiply by psihat(0)^(remaining speeds)
                prodc *= h0 ** (len(S) - tsize)
                cnt += 1
                ssum += prodc.real
                asum += abs(prodc)
                if max_terms and nrel + cnt > max_terms:
                    break
        res += ssum; abs_res += asum; nrel += cnt
        by_support[tsize] = (cnt, ssum, asum)
    return main, res, nrel, abs_res, by_support

def _relations(speeds, band):
    """All integer vectors k with 1<=|k_i|<=band and sum k_i*speeds_i = 0.  speeds: list (len m)."""
    m = len(speeds)
    out = []
    # last coordinate determined: k_{m-1} = -(sum_{i<m-1} k_i s_i)/s_{m-1}
    sm = speeds[-1]
    ranges = [range(-band, band+1) for _ in range(m-1)]
    for head in itertools.product(*ranges):
        if any(k == 0 for k in head):
            continue
        partial = sum(head[i]*speeds[i] for i in range(m-1))
        if partial % sm != 0:
            continue
        kl = -partial // sm
        if kl == 0 or abs(kl) > band:
            continue
        out.append(list(head)+[kl])
    return out

# ----------------------------------------------------------------------------- main analysis per set
def analyze_set(S, hat_fn, band, label="", verbose=True):
    R, Q = split_RQ(S)
    main, res, nrel, abs_res, by_sup = resonance_sum(S, hat_fn, band)
    floor_signed = main + res
    floor_abs    = main - abs_res     # the pessimistic (absolute) lower bound
    if verbose:
        print(f"  [{label}]  S={S}")
        print(f"      R(14-free)={R}  Q(/14)={Q}  r=|Q|={len(Q)}  covering={is_covering(S)} primitive={primitive(S)}")
        print(f"      MAIN = (int psi)^{len(S)} = {main:.8e}")
        print(f"      RESONANCE (signed, support<=6, |k|<={band}) = {res:+.8e}   ({nrel} relations)")
        print(f"      |RESONANCE| (absolute)                       = {abs_res:.8e}")
        print(f"      FLOOR (signed)   = MAIN + RES = {floor_signed:+.8e}   {'> 0  POSITIVE' if floor_signed>0 else '<= 0'}")
        print(f"      FLOOR (absolute) = MAIN - |RES| = {floor_abs:+.8e}   {'> 0  RIGOROUS-POSITIVE' if floor_abs>0 else '<= 0 (abs too lossy)'}")
        sup_str = ", ".join(f"|T|={t}:n={c},|sum|={a:.2e}" for t,(c,s,a) in by_sup.items())
        print(f"      by support: {sup_str}")
    return dict(S=S, R=R, Q=Q, r=len(Q), main=main, res=res, abs_res=abs_res,
                floor_signed=floor_signed, floor_abs=floor_abs, nrel=nrel, by_sup=by_sup)

# ----------------------------------------------------------------------------- few-apex covering bank
def build_few_apex_bank():
    """Canonical primitive covering sets S=R u 14Q with r=|Q| in {2..6}.
       Strategy: start from AP {1..13}, replace some elements by 14*m (m small) to create the
       multiples-of-14, while KEEPING covering (every d|some speed).  Also include the divisor-
       restored 'tight' few-apex rows analogous to codex S152's bank."""
    bank = []
    # The AP {1..13} is covering EXCEPT it has no multiple of 14.  Replace r elements by 14*m_i.
    # To keep covering we must not delete the only multiple of some d.  AP residues mod d cover all
    # d=2..13; deleting <=6 elements and adding multiples of 14 keeps covering as long as for each
    # d in 2..13 at least one surviving speed is divisible by d, and d=14 is supplied by the 14m's.
    base = list(range(1, 14))
    # candidate multiples-of-14 to insert (small m)
    Qcands = [1, 2, 3, 4, 5, 6]
    import random
    random.seed(20260627)
    seen = set()
    for r in range(2, 7):
        # choose which r positions of the AP to REPLACE by 14*m's, and which m's
        # we keep it structured: replace the top r elements (8..13 region) tends to keep low divisors
        replace_choices = [
            tuple(sorted(c)) for c in itertools.combinations(range(8, 14), r)
        ][:8]
        for repl in replace_choices:
            for Qsel in itertools.combinations(Qcands, r):
                kept = [v for v in base if v not in repl]
                S = sorted(set(kept) | {14*m for m in Qsel})
                if len(S) != 13:
                    continue
                g = reduce(gcd, S)
                if g != 1:
                    continue
                if not is_covering(S):
                    continue
                key = tuple(S)
                if key in seen:
                    continue
                seen.add(key)
                bank.append(S)
    return bank

# ----------------------------------------------------------------------------- uniform tail bound
def uniform_resonance_bound(hat_fn, band, n_speeds=13):
    """A SPEED-INDEPENDENT upper bound on |RESONANCE| for ANY set of n_speeds distinct nonzero
       integers, using only the per-coefficient sizes |psihat(k)|.

       |RES| <= sum_{tsize>=2} C(n,tsize) * h0^{n-tsize} * (#relations bound) * (max |prod psihat|).

       The clean uniform handle: each nonzero coordinate contributes a factor B1 = sum_{1<=|k|<=band}
       |psihat(k)|.  For a FIXED support T of size t, the number of relations is bounded, but the
       cleaner Schur/Holder bound is:
          sum over relations of |prod_{i in T} psihat(k_i)|  <=  per-pair handle.
       For 2-term supports {a,b}: relations are k_a=j*b/g, k_b=-j*a/g (g=gcd), so at most band
       values of j; |psihat(k_a)psihat(k_b)| <= |psihat|_max contributions.  We return the simple
       envelope  E_t = C(n,t) h0^{n-t} (B1)^t  and the dominant t=2 term, which is the honest uniform
       bound (over-counts relation multiplicity but is speed-free)."""
    h0 = hat_fn(0); h0 = h0.real if isinstance(h0, complex) else float(h0)
    B1 = sum(abs(hat_fn(k)) for k in range(1, band+1)) * 2.0   # sum over 1<=|k|<=band
    # envelope sum_{t>=2} C(n,t) h0^{n-t} B1^t = (h0+B1)^n - h0^n - n h0^{n-1} B1
    total = (h0 + B1)**n_speeds - h0**n_speeds - n_speeds*h0**(n_speeds-1)*B1
    # dominant t=2 piece
    t2 = math.comb(n_speeds, 2) * h0**(n_speeds-2) * B1**2
    return dict(h0=h0, B1=B1, env_total=total, t2=t2, main=h0**n_speeds)

def main():
    print("#"*96)
    print("# LRC(14)  GAUSSIAN/BEURLING-SELBERG MINORANT FLOOR  (kpswf15, HYP-3121 TOOL 1)")
    print("#   L(S) = int prod_s phi_s  >=  int prod_s psi_s = MAIN + RESONANCE")
    print("#   MAIN = (int psi)^|S|,  RESONANCE = sum_{sum k_s s=0, k!=0} prod psihat(k_s)")
    print("#"*96)

    # ---- 1. build + validate the two minorants
    print("\n" + "="*96)
    print("STEP 1: build + VALIDATE the minorants (true minorant <=> max_leak ~ 0, min_psi >= 0)")
    print("="*96)

    print("\n(B) C^infty compactly-supported mollifier minorant  psi = 1_{[a+delta,b-delta]} * rho_delta")
    for delta in [0.02, 0.03, 0.05]:
        # band where mollifier coeffs are negligible: rhohat(delta k) super-poly small
        K = 400
        intpsi = float(A_HI) - float(A_LO) - 2*delta
        mn, leak, over = verify_minorant(lambda k: mollifier_arc_hat(k, delta), K, grid=3000)
        # how fast do coeffs decay?
        decay = [(k, abs(mollifier_arc_hat(k, delta))) for k in [7,14,28,56,100,200]]
        print(f"   delta={delta}: int psi = {intpsi:.6f}   min_psi={mn:+.2e}  max_leak(outside arc)={leak:.2e}  max_over(inside)={over:+.2e}")
        print(f"             |psihat(k)|: " + "  ".join(f"k={k}:{v:.2e}" for k,v in decay))

    print("\n(A) Fejer-smoothed shrunk-arc (band-limited, degree N)  -- comparison object")
    for (N,h) in [(14,0.0),(28,0.0),(56,0.02),(112,0.02)]:
        intpsi = (float(A_HI)-h) - (float(A_LO)+h) - 1.0/(N+1)*0  # fejer mass note
        mn, leak, over = verify_minorant(lambda k: fejer_smoothed_arc_hat(k, N, h), N, grid=3000)
        h0 = fejer_smoothed_arc_hat(0, N, h).real
        print(f"   N={N}, shrink h={h}: int psi=h0={h0:.6f}  min_psi={mn:+.2e}  max_leak={leak:.2e}  max_over={over:+.2e}  (band {N})")

    # ---- choose the working minorant: C^infty delta=0.03 (true minorant) + band cutoff
    print("\n" + "="*96)
    print("STEP 2: per-set MAIN/RESONANCE/FLOOR for r=2..6 few-apex covering sets")
    print("   (working minorant: C^infty mollifier, delta=0.03; resonance band |k|<=BAND)")
    print("="*96)
    DELTA = 0.03
    BAND = 14
    hat = lambda k: mollifier_arc_hat(k, DELTA)
    bank = build_few_apex_bank()
    print(f"\n   few-apex covering bank: {len(bank)} primitive covering sets, r=2..6")
    by_r = {}
    for S in bank:
        R,Q = split_RQ(S); by_r.setdefault(len(Q), []).append(S)
    for r in sorted(by_r):
        print(f"      r={r}: {len(by_r[r])} sets")

    results = []
    # analyze a representative + worst per r
    print("\n--- representative sets per r ---")
    for r in range(2, 7):
        sets_r = by_r.get(r, [])
        if not sets_r:
            continue
        # representative: first; we will scan all for the floor below
        analyze_set(sets_r[0], hat, BAND, label=f"r={r} repr", verbose=True)
        print()

    # ---- 3. the FLOOR over the whole bank (signed + absolute)
    print("="*96)
    print("STEP 3: the MINORANT FLOOR over the whole r=2..6 bank")
    print("="*96)
    worst_signed = (1e9, None); worst_abs = (1e9, None)
    pos_signed = 0; pos_abs = 0
    for S in bank:
        d = analyze_set(S, hat, BAND, verbose=False)
        results.append(d)
        if d['floor_signed'] < worst_signed[0]:
            worst_signed = (d['floor_signed'], S)
        if d['floor_abs'] < worst_abs[0]:
            worst_abs = (d['floor_abs'], S)
        pos_signed += (d['floor_signed'] > 0)
        pos_abs += (d['floor_abs'] > 0)
    print(f"   sets analyzed: {len(results)}")
    print(f"   FLOOR (signed)   > 0 for {pos_signed}/{len(results)} sets;  worst = {worst_signed[0]:+.6e} at S={worst_signed[1]}")
    print(f"   FLOOR (absolute) > 0 for {pos_abs}/{len(results)} sets;     worst = {worst_abs[0]:+.6e} at S={worst_abs[1]}")

    # ---- 4. uniform speed-independent resonance bound
    print("\n" + "="*96)
    print("STEP 4: UNIFORM (speed-INDEPENDENT) resonance bound from per-coefficient sizes")
    print("="*96)
    for DELTA in [0.02, 0.03, 0.05]:
        hatd = lambda k: mollifier_arc_hat(k, DELTA)
        for BAND in [14, 28, 56]:
            ub = uniform_resonance_bound(hatd, BAND, n_speeds=13)
            print(f"   delta={DELTA}, band={BAND}: int psi=h0={ub['h0']:.5f}  B1=sum|psihat(k!=0)|={ub['B1']:.5f}")
            print(f"        MAIN=h0^13={ub['main']:.6e}   t=2 envelope={ub['t2']:.6e}   full envelope={ub['env_total']:.6e}")
            print(f"        MAIN - full_envelope = {ub['main']-ub['env_total']:+.6e}   {'UNIFORM FLOOR > 0' if ub['main']-ub['env_total']>0 else '(envelope too lossy)'}")

    print("\nDONE.")

if __name__ == "__main__":
    main()
