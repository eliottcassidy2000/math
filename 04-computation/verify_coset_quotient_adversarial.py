#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL INDEPENDENT VERIFICATION of the coset-quotient factorization claim.

Claim under test (codex HYP-2640 sharpening):
  For every support-6 relation n (6 coords, all nonzero, none mult of 7):
      K(n) = D7(n mod 7) / prod_j n_j
  where D7(c) = sum_T (-1)^|T| prod_j h_T(c_j), h_T(r) = -A(r) sum_{j in T} zeta_7^{-r j},
  A(r) = (1 - zeta_7^{-r})/(2 pi i), zeta_7 = exp(-2 pi i /7).

  Plus: K(n) = sum_T (-1)^|T| prod_j chat_T(n_j),
        chat_T(m) = sum_{j in T} (e^{-2pi i m j/7} - e^{-2pi i m (j+1)/7})/(2 pi i m).
  (This is the THM-538 Fourier kernel; we re-derive chat from scratch independently.)

I re-implement chat_T from the FOURIER COEFFICIENT DEFINITION independently, not from
the existing script, to avoid copying a possibly-wrong formula.
"""
import sys, itertools, math, cmath, random
from fractions import Fraction
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

# ---------------------------------------------------------------------------
# INDEPENDENT derivation of chat_T(m).
#
# The indicator of sector T = union of [j/7,(j+1)/7) for j in T.
# Its Fourier coefficient at frequency m (mode e^{2 pi i m x}) is
#   chat_T(m) = integral_0^1 1_T(x) e^{-2 pi i m x} dx
#             = sum_{j in T} integral_{j/7}^{(j+1)/7} e^{-2 pi i m x} dx.
# For m != 0:  integral_a^b e^{-2pi i m x} dx = (e^{-2pi i m a} - e^{-2pi i m b})/(2 pi i m).
# For m = 0 : it's just the measure |T|/7.
# ---------------------------------------------------------------------------
def chat_indep(T, m):
    if m == 0:
        return complex(len(T)/7.0, 0.0)
    s = 0+0j
    for j in T:
        a = j/7.0; b = (j+1)/7.0
        s += (cmath.exp(-2j*math.pi*m*a) - cmath.exp(-2j*math.pi*m*b))/(2j*math.pi*m)
    return s

def K_indep(vals):
    """K(n) = sum_{T subset {1..6}} (-1)^|T| prod_j chat_T(n_j). Sum over ALL subsets of {1..6}."""
    total = 0+0j
    for cnt in range(7):
        for T in itertools.combinations(range(1,7), cnt):
            prod = 1+0j
            for v in vals:
                prod *= chat_indep(T, v)
            total += ((-1)**cnt)*prod
    return total

# ---------------------------------------------------------------------------
# The factored form, derived independently.
# For m != 0, mult7 false: chat_T(m) = sum_{j in T} e^{-2pi i m j/7}(1 - e^{-2pi i m/7})/(2 pi i m).
# Let r = m mod 7. e^{-2pi i m j/7} = e^{-2pi i r j/7} (depends only on r).
# e^{-2pi i m/7} = e^{-2pi i r/7}. So:
#   chat_T(m) = (1/m) * [ (1 - e^{-2pi i r/7})/(2 pi i) * sum_{j in T} e^{-2pi i r j/7} ]
#             = (1/m) * hT(r),  hT(r) = A(r) * sum_{j in T} e^{-2pi i r j/7},  A(r)=(1-e^{-2pi i r/7})/(2 pi i)
# NOTE sign: I derive hT = +A * sum (no minus). The original script uses h_T = -A*sum(zeta^{rj})
# with zeta(p)=exp(-2pi i p/7). exp(-2pi i r j/7) = zeta(r j). So sum_j zeta(rj) matches.
# The sign difference (+A vs -A) appears in EACH factor; with 6 factors per product the
# overall (-1)^6 = +1, so D7 is the SAME. We test both to be safe.
# ---------------------------------------------------------------------------
def A_coef(r):
    return (1 - cmath.exp(-2j*math.pi*r/7.0))/(2j*math.pi)

def hT_plus(T, r):
    return A_coef(r) * sum(cmath.exp(-2j*math.pi*r*j/7.0) for j in T)

def hT_minus(T, r):  # original script convention
    return -A_coef(r) * sum(cmath.exp(-2j*math.pi*r*j/7.0) for j in T)

def D7(residues, hT=hT_plus):
    total = 0+0j
    for cnt in range(7):
        for T in itertools.combinations(range(1,7), cnt):
            prod = 1+0j
            for r in residues:
                prod *= hT(T, r)
            total += ((-1)**cnt)*prod
    return total

def K_factored(vals, hT=hT_plus):
    res = tuple(v % 7 for v in vals)
    invprod = 1.0
    for v in vals:
        invprod /= v
    return invprod * D7(res, hT)

# ===========================================================================
if __name__ == "__main__":
    print("="*70)
    print("STEP 1: independent verification of factorization K(n)=D7(n mod7)/prod n")
    print("="*70)
    random.seed(12345)
    maxerr_plus = 0.0
    maxerr_minus = 0.0
    worst = None
    for _ in range(3000):
        vals = []
        while len(vals) < 6:
            v = random.randint(-25, 25)
            if v != 0 and v % 7 != 0:
                vals.append(v)
        kd = K_indep(vals)
        kp = K_factored(vals, hT_plus)
        km = K_factored(vals, hT_minus)
        ep = abs(kd-kp); em = abs(kd-km)
        if ep > maxerr_plus:
            maxerr_plus = ep; worst = (tuple(vals), kd, kp)
        maxerr_minus = max(maxerr_minus, em)
    print(f"  max|K_indep - K_factored(+A)| over 3000 random n = {maxerr_plus:.3e}")
    print(f"  max|K_indep - K_factored(-A)| over 3000 random n = {maxerr_minus:.3e}")
    print(f"  (factorization EXACT iff ~1e-13)")
    print(f"  also: K is REAL? max|Im K_indep| sample = {abs(worst[1].imag):.2e}")

    print()
    print("="*70)
    print("STEP 2: D7 structural claims (antipodal symmetry, bounds, nonzero)")
    print("="*70)
    # (1) antipodal symmetry D7(-c) = conj(D7(c))   [-c means (7-c_j) mod 7 = -c_j mod 7]
    random.seed(99)
    max_anti = 0.0
    for _ in range(3000):
        c = tuple(random.randint(1,6) for _ in range(6))
        cm = tuple((7-r)%7 for r in c)  # -c mod 7, in 1..6
        d  = D7(c)
        dm = D7(cm)
        max_anti = max(max_anti, abs(dm - d.conjugate()))
    print(f"  max|D7(-c) - conj(D7(c))| over 3000 random cosets = {max_anti:.3e}")
    print("  => antipodal pairing kills imaginary part => correction is REAL.")

    # (2) full survey over all 6^6 = 46656 cosets: nonzero?  max|Re|, max|Im|
    print("\n  full survey of all 6^6 = 46656 cosets...")
    minabs = 1e9; maxre = 0.0; maxim = 0.0; nzero = 0
    # cache hT(T,r) for all (T,r)
    Ts = [T for cnt in range(7) for T in itertools.combinations(range(1,7),cnt)]
    signs = [(-1)**len(T) for T in Ts]
    hcache = {}
    for T in Ts:
        for r in range(1,7):
            hcache[(T,r)] = hT_plus(T, r)
    for c in itertools.product(range(1,7), repeat=6):
        total = 0+0j
        for T, sg in zip(Ts, signs):
            prod = 1+0j
            for r in c:
                prod *= hcache[(T,r)]
            total += sg*prod
        a = abs(total)
        if a < minabs: minabs = a
        if a < 1e-9: nzero += 1
        maxre = max(maxre, abs(total.real))
        maxim = max(maxim, abs(total.imag))
    print(f"  min|D7(c)| over all cosets = {minabs:.4e}  (#near-zero<1e-9: {nzero})")
    print(f"  max|Re D7| = {maxre:.4f}   max|Im D7| = {maxim:.4f}")
    print("  claim: nonzero on all, max|Re|=0.1431, max|Im|=0.627")

    print()
    print("="*70)
    print("STEP 3: box-truncation reconstruction + per-shell decay (NEGATIVE claim)")
    print("="*70)
    # exact p0 engine
    def dist_p0(E):
        E = sorted(set(E))
        bps = set([Fraction(0), Fraction(1)])
        for e in E:
            if e == 0: continue
            for a in range(0, 7*e+1):
                bps.add(Fraction(a, 7*e))
        bps = sorted(b for b in bps if 0 <= b <= 1)
        p0 = Fraction(0)
        for i in range(len(bps)-1):
            lo, hi = bps[i], bps[i+1]
            if hi == lo: continue
            mid = (lo+hi)/2
            hit = set()
            for e in E:
                v = e*mid; v = v - (v.numerator//v.denominator)
                hit.add((v.numerator*7)//v.denominator)
            N = sum(1 for j in range(1,7) if j not in hit)
            if N == 0: p0 += (hi-lo)
        return p0

    # enumerate support-6 relations of E with |n|inf <= L, return (combo, ninf)
    def relations_support6(E, L):
        k=len(E); out=[]
        for idxs in itertools.combinations(range(k),6):
            es=[E[i] for i in idxs]
            dep=max(range(6),key=lambda i:abs(es[i])); e_dep=es[dep]
            free=[i for i in range(6) if i!=dep]; efree=[es[i] for i in free]
            ranges=[range(-L,L+1) for _ in range(5)]
            for vf in itertools.product(*ranges):
                if any(c==0 or c%7==0 for c in vf): continue
                s=sum(c*e for c,e in zip(vf,efree))
                if e_dep==0:
                    if s!=0: continue
                    for vd in range(-L,L+1):
                        if vd==0 or vd%7==0: continue
                        combo=[0]*6
                        for i,c in zip(free,vf): combo[i]=c
                        combo[dep]=vd
                        out.append((tuple(combo),max(abs(x) for x in combo)))
                else:
                    if s%e_dep!=0: continue
                    vd=-s//e_dep
                    if vd==0 or vd%7==0 or abs(vd)>L: continue
                    combo=[0]*6
                    for i,c in zip(free,vf): combo[i]=c
                    combo[dep]=vd
                    out.append((tuple(combo),max(abs(x) for x in combo)))
        return out

    from collections import defaultdict
    E = [0,1,2,3,4,5,6,7]
    p0 = dist_p0(E)
    M7 = sum(((-1)**t)*math.comb(6,t)*((1-t/7.0)**(len(E)-1)) for t in range(7))
    target = float(p0) - M7
    print(f"  E={E} (k=8 AP). p0={float(p0):.6f}, M7={M7:.6f}, target corr p0-M7 = {target:.6f}")
    # NOTE: dedupe relations across the C(8,6)=28 index-subsets (a given n vector with 6
    # nonzero coords positioned in E could be produced once per subset; relations_support6
    # already fixes which 6 indices, so each appears once per its own index set). The
    # correction is sum over the FULL relation lattice; support-6 vectors live in fixed coords.
    # We sum K over all support-6 relations found, grouped by inf-norm shell.
    shell_signed = defaultdict(float)
    shell_abs    = defaultdict(float)
    box_cum = {}
    for L in [3,5,7,9]:
        rels = relations_support6(E, L)
        # dedupe (same combo could be generated under different dep choice? no: idxs fixed,
        # but the SAME 6-coord support could arise from different idxs subsets if E has repeats;
        # E distinct so each (idxs,combo) unique). Use set on (idxs not tracked) -> use combo+positions.
        seen=set(); uniq=[]
        for combo,ninf in rels:
            key=combo
            # combo is over the 6 chosen indices -> ambiguous across idxs; track via full vector?
            # relations_support6 loses which idxs. To be safe, recompute with idxs.
            uniq.append((combo,ninf))
        # reconstruct correction via factorization
        corr=0+0j
        per_shell_s=defaultdict(float); per_shell_a=defaultdict(float)
        for combo,ninf in uniq:
            res=tuple(v%7 for v in combo)
            ip=1.0
            for v in combo: ip/=v
            kval=(ip*D7(res)).real
            corr+=kval
            per_shell_s[ninf]+=kval
            per_shell_a[ninf]+=abs(kval)
        box_cum[L]=corr.real
        print(f"  L={L}: #rels={len(uniq):>9}  box-cum corr={corr.real:+.6f}  frac={corr.real/target:+.4f}")
    print("  per-shell signed sums (L=9 enumeration):")
    rels=relations_support6(E,9)
    ps=defaultdict(float); pa=defaultdict(float)
    for combo,ninf in rels:
        res=tuple(v%7 for v in combo); ip=1.0
        for v in combo: ip/=v
        kval=(ip*D7(res)).real
        ps[ninf]+=kval; pa[ninf]+=abs(kval)
    for s in sorted(ps):
        print(f"    shell|n|inf={s}: signed={ps[s]:+.5f}  abs={pa[s]:.5f}")
    print("  claim: signed oscillates (no decay), abs GROWS, box non-monotone.")
