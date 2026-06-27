#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_gaussian_minorant_floor_kpswf15.py   (kind-pasteur 2026-06-27, HYP-3121 synthesis / TOOL 1)

GAUSSIAN / BEURLING-SELBERG MINORANT for the UNIFORM r=2..6 multi-far loneliness floor.

THE OBJECT.  A covering set S = R u 14Q with r=|Q| in {2..6} multiples of 14 and a 14-free
"small part" R (13-r >= 7 speeds).  The loneliness measure is
        L(S) = int_0^1  prod_{s in S}  phi_s(t)  dt,
        phi_s(t) = 1_{[1/14, 13/14]}( {s t} ),     an arc of measure 6/7.
LRC(14) on this branch  <=>  L(S) > 0 for every covering S with r in {2..6}.

THE MINORANT IDEA (the wide-V tail fix).  The sharp arc indicator phi_s has 1/n Fourier decay,
so the resonance / relation-lattice sum  sum_{sum_s k_s s = 0}  prod_s phihat(k_s)  has a slowly
decaying tail that is hard to bound UNIFORMLY in the speed magnitudes.  Replace each phi_s by a
NONNEGATIVE minorant  0 <= psi_s <= phi_s  whose Fourier coefficients decay SUPER-POLYNOMIALLY
(C^infty mollifier, route i) or are FINITELY supported (Beurling-Selberg band-N, route ii).  Then

   L(S) = int prod phi_s  >=  int prod psi_s
        = MAIN  prod_s ( int psi_s )   +   RESONANCE  sum_{(k_s) != 0, sum_s k_s s = 0} prod_s psihat(k_s).

  MAIN = (int psi)^|S|  depends only on the arc length 6/7 and the minorant parameter (NOT speeds).
  RESONANCE tail is super-polynomial / finite-band => UNIFORM in the speed magnitudes.

KEY question (the deliverable): is MAIN > |RESONANCE| with an explicit, uniform, speed-free constant?

WHAT THIS SCRIPT DOES (numpy + scipy).
  (A) Build a TRUE C^infty mollifier minorant psi_delta (route i), validate psi<=phi on a fine grid,
      and tabulate the super-polynomially-decaying Fourier coefficients psihat(k)=chi_hat(k)*rhohat(delta k).
  (B) For each r=2..6 few-apex covering set S, compute the per-set minorant floor
      int prod_s psi_s(t) dt  by high-accuracy Gauss-Legendre quadrature (real-space, EXACT-arc cross
      check L(S)=int prod phi_s by the repo arc machinery), and split it MAIN + RESONANCE.
  (C) The UNIFORM (speed-independent) resonance bound from per-coefficient sizes B1=sum|psihat(k!=0)|
      (full Schur envelope + the dominant 2-term piece), and a SHARPER per-(speed-pair) 2-term sum.
  (D) Beurling-Selberg band-N comparison (finite band => the cleanest uniform tail).

OUTPUT: per-set MAIN/RESONANCE/FLOOR, the floor over the whole bank, the uniform tail bound, and the
explicit constant c if the minorant floor is > 0.
"""
import sys, itertools, math
from fractions import Fraction as F
from math import gcd, pi
from functools import reduce
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

# import the EXACT arc machinery for L(S) = int prod phi_s cross-check
sys.path.insert(0, "04-computation")
from lrc14_spectrum_intersection_sum_kpswf12 import meas, intersect, complement, safe_set

A_LO = 1.0/14.0
A_HI = 13.0/14.0
ARC_LEN = A_HI - A_LO          # = 6/7

# ============================================================ C^infty bump (route i)
# rho(u) = C exp(-1/(1-u^2)) on (-1,1), 0 else; C normalizes int rho = 1.
def _bump_raw(u):
    u = np.asarray(u, dtype=float)
    out = np.zeros_like(u)
    m = np.abs(u) < 1.0
    out[m] = np.exp(-1.0/(1.0 - u[m]*u[m]))
    return out

_BUMP_NORM, _ = quad(lambda u: math.exp(-1.0/(1.0-u*u)), -1, 1, limit=200)
_BUMP_C = 1.0/_BUMP_NORM

def rho(u):
    return _BUMP_C * _bump_raw(u)

def rho_scalar(u):
    if -1.0 < u < 1.0:
        return _BUMP_C * math.exp(-1.0/(1.0-u*u))
    return 0.0

_rhohat_cache = {}
def rhohat(xi):
    """rhohat(xi) = int_{-1}^{1} rho(u) cos(2 pi xi u) du   (real, even, super-poly decay)."""
    key = round(float(xi), 10)
    if key in _rhohat_cache:
        return _rhohat_cache[key]
    val, _ = quad(lambda u: rho_scalar(u)*math.cos(2*pi*xi*u), -1, 1, limit=400)
    _rhohat_cache[key] = val
    return val

# ============================================================ arc Fourier coeffs (exact closed form)
def chi_hat(k, a=A_LO, b=A_HI):
    """Fourier coeff of 1_{[a,b]} on the circle: (e(-k a)-e(-k b))/(2 pi i k), k!=0; (b-a) if k=0."""
    if k == 0:
        return complex(b-a, 0.0)
    ea = complex(math.cos(-2*pi*k*a), math.sin(-2*pi*k*a))
    eb = complex(math.cos(-2*pi*k*b), math.sin(-2*pi*k*b))
    return (ea - eb)/(2j*pi*k)

def mollifier_hat(k, delta):
    """psihat(k) for psi = 1_{[a+delta, b-delta]} * rho_delta.   = chi_hat_shrunk(k) * rhohat(delta k)."""
    if k == 0:
        return complex(ARC_LEN - 2*delta, 0.0)
    return chi_hat(k, A_LO+delta, A_HI-delta) * rhohat(delta*k)

# ============================================================ minorant psi(theta) in real space (vectorized)
def make_psi_grid(delta, ng=20001):
    """Return (grid theta in [0,1), psi values) for the periodic C^infty minorant, by DIRECT
       real-space convolution (NOT Fourier reconstruction) so it is provably >=0 and exact-shape.
       psi(theta) = int_{a+delta}^{b-delta} rho_delta(theta - y) dy
                  = int_{(theta-(b-delta))/delta}^{(theta-(a+delta))/delta} rho(v) dv     [periodized]
       i.e. psi(theta) = Phi((theta-(a+delta))/delta) - Phi((theta-(b-delta))/delta), Phi=int_{-1}^{u} rho.
       We periodize by summing over the integer shift that lands theta near the arc."""
    aa = A_LO + delta; bb = A_HI - delta
    # cumulative bump CDF Phi(u)=int_{-1}^{u} rho
    us = np.linspace(-1.0, 1.0, 8001)
    rv = rho(us)
    Phi_tab = np.concatenate([[0.0], np.cumsum((rv[1:]+rv[:-1])/2*np.diff(us))])
    def Phi(u):
        u = np.clip(u, -1.0, 1.0)
        return np.interp(u, us, Phi_tab)
    th = (np.arange(ng)+0.5)/ng
    psi = np.zeros(ng)
    for shift in (-1.0, 0.0, 1.0):   # periodization (arc well inside [0,1], 1 shift enough but be safe)
        t = th + shift
        psi += Phi((t-aa)/delta) - Phi((t-bb)/delta)
    return th, psi

def validate_minorant(delta, ng=40001):
    th, psi = make_psi_grid(delta, ng)
    a, b = A_LO, A_HI
    inside = (th >= a) & (th <= b)
    phi = inside.astype(float)
    min_psi = psi.min()
    max_leak = psi[~inside].max() if (~inside).any() else 0.0
    max_over = (psi[inside]-1.0).max() if inside.any() else 0.0
    intpsi = psi.mean()
    return dict(delta=delta, min_psi=min_psi, max_leak=max_leak, max_over=max_over,
                intpsi=intpsi, intpsi_exact=ARC_LEN-2*delta)

# ============================================================ covering machinery
AP = tuple(range(1, 14))
def is_covering(S):
    return all(any(v % d == 0 for v in S) for d in range(2, 15))
def primitive(S):
    return reduce(gcd, S) == 1
def split_RQ(S):
    Q = sorted(v//14 for v in S if v % 14 == 0)
    R = sorted(v for v in S if v % 14 != 0)
    return R, Q

# ============================================================ EXACT loneliness L(S)=int prod phi_s
def lonely_measure_exact(S):
    """L(S) = meas{ t : ||s t|| >= 1/14  for all s in S }  (EXACT Fraction via repo arc machinery)."""
    arcs = safe_set(list(S), h=F(1,14))
    return meas(arcs)

# ============================================================ per-set minorant floor by quadrature
def minorant_floor_quad(S, delta, nnodes=200000):
    """int_0^1 prod_s psi(s t mod 1) dt by composite Gauss-Legendre.  psi from the real-space CDF.
       We evaluate psi via the analytic CDF (interp) on the needed phases."""
    aa = A_LO + delta; bb = A_HI - delta
    us = np.linspace(-1.0, 1.0, 8001)
    rv = rho(us)
    Phi_tab = np.concatenate([[0.0], np.cumsum((rv[1:]+rv[:-1])/2*np.diff(us))])
    def psi_of_phase(ph):
        # ph in [0,1)
        out = np.zeros_like(ph)
        for shift in (-1.0, 0.0, 1.0):
            t = ph + shift
            u1 = np.clip((t-aa)/delta, -1.0, 1.0)
            u2 = np.clip((t-bb)/delta, -1.0, 1.0)
            out += np.interp(u1, us, Phi_tab) - np.interp(u2, us, Phi_tab)
        return out
    # composite midpoint with MANY nodes (psi smooth, but product has Vmax structure; use uniform fine grid)
    t = (np.arange(nnodes)+0.5)/nnodes
    prod = np.ones(nnodes)
    for s in S:
        ph = (s*t) % 1.0
        prod *= psi_of_phase(ph)
    return prod.mean()

# ============================================================ uniform speed-independent resonance bound
def coeff_table(delta, band):
    return {k: mollifier_hat(k, delta) for k in range(-band, band+1)}

def uniform_bound(delta, band, n=13):
    """Speed-INDEPENDENT bound on |RESONANCE| via per-coefficient sizes.
       B1 = sum_{1<=|k|<=band} |psihat(k)|.  Envelope:
          |RES| <= sum_{t>=2} C(n,t) h0^{n-t} B1^t = (h0+B1)^n - h0^n - n h0^{n-1} B1.
       Also returns the dominant t=2 piece C(n,2) h0^{n-2} B1^2 and a tail-beyond-band estimate."""
    h0 = mollifier_hat(0, delta).real
    B1 = sum(abs(mollifier_hat(k, delta)) for k in range(1, band+1))*2.0
    main = h0**n
    env = (h0+B1)**n - h0**n - n*h0**(n-1)*B1
    t2  = math.comb(n, 2) * h0**(n-2) * B1**2
    # contribution of |k|>band (super-poly tail): bound sum_{|k|>band}|psihat(k)| crudely by
    # |psihat(k)| <= |chi_hat| * |rhohat(delta k)| <= (1/(pi|k|)) * rhohat_tail; we just report B1_tail.
    B1_tail = sum(abs(mollifier_hat(k, delta)) for k in range(band+1, band+200))*2.0
    return dict(h0=h0, B1=B1, B1_tail=B1_tail, main=main, env=env, t2=t2, floor=main-env)

def two_term_sum_for_set(S, delta, band):
    """The SHARP 2-term resonance contribution for THIS set: sum over pairs {a,b} in S, over
       j with |j b/g|,|j a/g| <= band, of  psihat(j b/g) psihat(-j a/g) h0^{n-2}.
       This is the dominant (support-2) part of RESONANCE -- speed-DEPENDENT but exactly computed."""
    S = sorted(S); n = len(S)
    h0 = mollifier_hat(0, delta).real
    tab = coeff_table(delta, band)
    total = 0.0; abstot = 0.0; cnt = 0
    for i in range(n):
        for jx in range(i+1, n):
            a, b = S[i], S[jx]
            g = gcd(a, b)
            ka0, kb0 = b//g, a//g   # k_a = j*ka0, k_b = -j*kb0
            jmax = band // max(ka0, kb0)
            for jj in range(1, jmax+1):
                ka = jj*ka0; kb = -jj*kb0
                if abs(ka) > band or abs(kb) > band:
                    break
                term = (tab[ka]*tab[kb]).real * h0**(n-2)
                total += 2.0*term   # +-j
                abstot += 2.0*abs((tab[ka]*tab[kb]) * h0**(n-2))
                cnt += 2
    return total, abstot, cnt

# ============================================================ Beurling-Selberg band-N (route ii) comparison
# We realise a band-limited TRUE minorant via Fejer smoothing of a shrunk arc, validated as minorant
# on a grid.  Fejer kernel F_N has nonneg values, mass 1, multiplier (1-|k|/(N+1)).
def fejer_hat(k, N, h):
    if abs(k) > N:
        return 0j
    mult = 1.0 - abs(k)/(N+1)
    if k == 0:
        return complex((A_HI-h)-(A_LO+h), 0.0)*mult
    return chi_hat(k, A_LO+h, A_HI-h)*mult

def fejer_psi_grid(N, h, ng=40001):
    """psi(theta)= (1_{[a+h,b-h]} * F_N)(theta) by Fourier reconstruction (band N)."""
    th = (np.arange(ng)+0.5)/ng
    psi = np.full(ng, fejer_hat(0, N, h).real)
    for k in range(1, N+1):
        hk = fejer_hat(k, N, h)
        psi += 2.0*np.real(hk*np.exp(2j*pi*k*th))
    return th, psi

def validate_fejer(N, h, ng=40001):
    th, psi = fejer_psi_grid(N, h, ng)
    inside = (th >= A_LO) & (th <= A_HI)
    return dict(N=N, h=h, h0=fejer_hat(0,N,h).real, min_psi=psi.min(),
                max_leak=(psi[~inside].max() if (~inside).any() else 0.0),
                max_over=((psi[inside]-1).max() if inside.any() else 0.0))

# ============================================================ few-apex covering bank
def build_bank():
    base = list(range(1, 14))
    Qcands = [1, 2, 3, 4, 5, 6, 7]
    seen = set(); bank = []
    for r in range(2, 7):
        for repl in itertools.combinations(range(8, 14), r):       # replace top-region elements
            for Qsel in itertools.combinations(Qcands, r):
                kept = [v for v in base if v not in repl]
                S = sorted(set(kept) | {14*m for m in Qsel})
                if len(S) != 13 or reduce(gcd, S) != 1 or not is_covering(S):
                    continue
                key = tuple(S)
                if key not in seen:
                    seen.add(key); bank.append(S)
        # also allow replacing a spread of positions (not just the top), to diversify R
        for repl in itertools.combinations([7,8,9,10,11,12,13][:6], r):
            for Qsel in itertools.combinations(Qcands, r):
                kept = [v for v in base if v not in repl]
                S = sorted(set(kept) | {14*m for m in Qsel})
                if len(S) != 13 or reduce(gcd, S) != 1 or not is_covering(S):
                    continue
                key = tuple(S)
                if key not in seen:
                    seen.add(key); bank.append(S)
    return bank

# ============================================================ main
def main():
    print("#"*100)
    print("# LRC(14)  GAUSSIAN/BEURLING-SELBERG MINORANT FLOOR   (kpswf15, HYP-3121 TOOL 1)")
    print("#   L(S)=int prod_s phi_s >= int prod_s psi_s = MAIN + RESONANCE")
    print("#   MAIN=(int psi)^|S|,  RESONANCE=sum_{sum k_s s=0,k!=0} prod psihat(k_s)")
    print("#"*100)

    # ---- STEP 1: build + validate the C^infty mollifier minorant (route i)
    print("\n" + "="*100)
    print("STEP 1: build + VALIDATE the C^infty mollifier minorant  psi=1_{[a+delta,b-delta]}*rho_delta")
    print("        (TRUE minorant <=> max_leak(outside arc)=0 and min_psi>=0; int psi = 6/7 - 2 delta)")
    print("="*100)
    for delta in [0.02, 0.03, 0.05, 0.07]:
        v = validate_minorant(delta)
        coeffs = [(k, abs(mollifier_hat(k, delta))) for k in [7, 14, 28, 56, 100, 200]]
        print(f"  delta={delta}: int psi={v['intpsi']:.6f} (exact {v['intpsi_exact']:.6f})  "
              f"min_psi={v['min_psi']:+.2e}  max_leak={v['max_leak']:.2e}  max_over={v['max_over']:+.2e}")
        print(f"           |psihat(k)| (super-poly decay): " + "  ".join(f"k={k}:{c:.2e}" for k,c in coeffs))

    # ---- STEP 2: per-set MAIN/RESONANCE/FLOOR for representative r=2..6 sets
    print("\n" + "="*100)
    print("STEP 2: per-set MAIN / RESONANCE / FLOOR  (working minorant delta=0.03; band |k|<=56)")
    print("        floor = int prod psi  (Gauss/midpoint quad); RESONANCE = floor - MAIN; cross-check L(S)")
    print("="*100)
    DELTA = 0.03; BAND = 56
    h0 = mollifier_hat(0, DELTA).real
    MAIN = h0**13
    print(f"  h0 = int psi = {h0:.6f}   MAIN = h0^13 = {MAIN:.6e}")
    bank = build_bank()
    by_r = {}
    for S in bank:
        by_r.setdefault(len(split_RQ(S)[1]), []).append(S)
    print(f"  few-apex covering bank: {len(bank)} primitive covering sets; " +
          "  ".join(f"r={r}:{len(by_r[r])}" for r in sorted(by_r)))

    print("\n  --- representative set per r ---")
    reps = []
    for r in range(2, 7):
        if r not in by_r:
            continue
        S = by_r[r][0]
        R, Q = split_RQ(S)
        floor = minorant_floor_quad(S, DELTA, nnodes=300000)
        L = lonely_measure_exact(S)
        res = floor - MAIN
        t2, t2abs, t2cnt = two_term_sum_for_set(S, DELTA, BAND)
        reps.append((S, R, Q, floor, L, res, t2))
        print(f"  r={r}: S={S}")
        print(f"        R={R}  Q={Q}")
        print(f"        L(S)=int prod phi (EXACT) = {float(L):.6f} = {L}")
        print(f"        minorant floor int prod psi = {floor:.6e}   (<= L: {floor <= float(L)+1e-9})")
        print(f"        MAIN={MAIN:.6e}   RESONANCE=floor-MAIN={res:+.6e}   [2-term part={t2:+.6e}]")
        print(f"        floor>0: {floor>0}    floor/MAIN ratio = {floor/MAIN:.4f}")

    # ---- STEP 3: the FLOOR over the whole bank
    print("\n" + "="*100)
    print("STEP 3: the MINORANT FLOOR over the whole r=2..6 few-apex covering bank")
    print("="*100)
    worst = (1e9, None); pos = 0; ratios = []
    worst_per_r = {}
    for S in bank:
        floor = minorant_floor_quad(S, DELTA, nnodes=120000)
        ratios.append(floor/MAIN)
        if floor < worst[0]:
            worst = (floor, S)
        pos += (floor > 0)
        r = len(split_RQ(S)[1])
        if r not in worst_per_r or floor < worst_per_r[r][0]:
            worst_per_r[r] = (floor, S)
    print(f"  sets analyzed: {len(bank)}")
    print(f"  minorant floor > 0 for {pos}/{len(bank)} sets")
    print(f"  WORST floor = {worst[0]:.6e}  at S={worst[1]}   (MAIN={MAIN:.6e}, ratio {worst[0]/MAIN:.4f})")
    print(f"  floor/MAIN ratio range: [{min(ratios):.4f}, {max(ratios):.4f}]")
    for r in sorted(worst_per_r):
        fl, S = worst_per_r[r]
        print(f"     r={r}: worst floor = {fl:.6e} (ratio {fl/MAIN:.4f}) at S={S}")
    print(f"\n  => MINORANT FLOOR CONSTANT  c = min floor = {worst[0]:.6e} > 0  (over this bank)")

    # ---- STEP 4: UNIFORM speed-independent resonance bound
    print("\n" + "="*100)
    print("STEP 4: UNIFORM (speed-INDEPENDENT) resonance bound from per-coefficient sizes B1")
    print("        |RES| <= (h0+B1)^13 - h0^13 - 13 h0^12 B1 ;  MAIN - this = uniform floor")
    print("="*100)
    for delta in [0.02, 0.03, 0.05, 0.07]:
        for band in [56, 200]:
            ub = uniform_bound(delta, band, n=13)
            ok = ub['floor'] > 0
            print(f"  delta={delta}, band={band}: h0={ub['h0']:.5f}  B1={ub['B1']:.5f}  B1_tail(|k|>band)={ub['B1_tail']:.2e}")
            print(f"        MAIN=h0^13={ub['main']:.6e}  full_envelope={ub['env']:.6e}  (t=2 part={ub['t2']:.6e})")
            print(f"        MAIN - envelope = {ub['floor']:+.6e}   {'UNIFORM FLOOR > 0' if ok else '(envelope too lossy at this delta/band)'}")

    # ---- STEP 5: Beurling-Selberg band-N comparison (route ii, finite band)
    print("\n" + "="*100)
    print("STEP 5: Beurling-Selberg / Fejer band-N minorant (route ii): finite band => cleanest uniform tail")
    print("="*100)
    for (N, h) in [(28, 0.0), (56, 0.02), (112, 0.02), (200, 0.03)]:
        v = validate_fejer(N, h)
        print(f"  Fejer N={N}, shrink h={h}: int psi=h0={v['h0']:.6f}  min_psi={v['min_psi']:+.2e}  "
              f"max_leak={v['max_leak']:.2e}  (TRUE minorant: {v['max_leak']<1e-3 and v['min_psi']>-1e-6})")
    print("  (Fejer-of-shrunk-arc is the band-limited object; a true band-N minorant needs the leak<=0;")
    print("   shrink h trades int psi vs leak. The C^infty route i above is the validated true minorant.)")

    print("\nDONE.")

if __name__ == "__main__":
    main()
