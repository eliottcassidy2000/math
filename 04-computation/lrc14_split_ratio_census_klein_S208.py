#!/usr/bin/env python3
"""
lrc14_split_ratio_census_klein_S208.py

HYP-5691: the P∪L-split ratio census of the k>=8 covering family, and the
honest per-instrument firing map for the realization (THM-527-A) chain.

CONTEXT (klein-2026-07-09-S208). THM-665 (monad-explorer-S1) consequence 2
states "the a-priori existence route never fires on covering clusters
(V/spread ~ 1 there)" using the ALL-RUNNER co-offset convention
(spread = Vmax - v_min ~ Vmax whenever small speeds are present -- always,
for covering sets). But THM-527's architecture splits S = P ∪ L with
P = S ∩ [1,13] handled MEASURE-side by G_P; the cluster co-offsets are
E = {Vmax - u : u in L}, so the relevant spread is Vmax - min(L), which the
covering constraint does NOT force to be proportional: L's covering duty
(multiples of the q's P misses) can be discharged inside [0.643*Vmax, Vmax].

QUESTIONS:
 (Q1) Stratify the k>=8 covering family by split-ratio r = Vmax/spread(E).
      Both strata (confined r>2.8 and proportional r<=2.8) should be
      populated => THM-665 cons.2 needs the split-convention refinement.
 (Q2) Per stratum, which realization instrument HONESTLY fires?
      (i)  aliasing existence (THM-665): V > V0 = sqrt(TV(W')/(12*Int W))
      (ii) drift embed (klein-S205 LRCDriftEmbed.minReach_ge_of_driftGap):
           exists grid point x_j = j/V in G_P^slack whose tooth-gap (a, a+g)
           satisfies 1/7 + 2*spread*(a+g/2)/V < g  [the EXACT Lean hypothesis]
      (iii) j=1 fat period (LEM-010): fires only if P = emptyset
      (iv) direct M(S) >= 1/14 (ground truth, exact rational local-maxima)
 (Q3) The ADAPTIVE-SPLIT ladder: moving s slow cluster members into P~:
      |P~| = |P| + s, k~ = k - s. If k~ <= 7 the mu=1 pigeonhole applies and
      the floor becomes meas(G_P~) alone (Bonferroni > 0 iff |P~| <= 6).
      Count the strata (s, |P~|, k~) to see where the ladder closes the
      MEASURE floor and where it hits the apex-7/Fraenkel wall.

All W/TV computations are EXACT (Fraction arithmetic on the piecewise-linear
structure); M(S) is exact via the p/(v_i+v_j) local-maxima method (S206).
"""

from fractions import Fraction
from math import gcd, isqrt
import itertools, random, sys

random.seed(20260709)

# ----------------------------------------------------------------------
# exact cluster functionals: gaps of {frac(e_i x)}, W(x), Int W, TV(W')
# ----------------------------------------------------------------------

def breakpoints(E, extra_denoms=()):
    """All breakpoints in [0,1) of the piecewise-linear gap structure of
    {frac(e x) : e in E}: collisions m/|e_i-e_j| and wraps m/e_i.
    (theta-crossings of W are NOT breakpoints of the phase order; W's own
    kinks are handled by refining with theta-crossing points per cell later
    -- for Int W and TV(W') we integrate the max/sum structure per cell by
    subdividing cells at all gap-theta crossings, computed exactly.)"""
    dens = set()
    E = sorted(E)
    for i in range(len(E)):
        if E[i] > 0:
            dens.add(E[i])
        for j in range(i + 1, len(E)):
            d = E[j] - E[i]
            if d > 0:
                dens.add(d)
    for d in extra_denoms:
        if d > 0:
            dens.add(d)
    pts = set()
    for d in dens:
        for m in range(d):
            pts.add(Fraction(m, d))
    return sorted(pts)


def gaps_at(E, x):
    """circular gaps of the multiset {frac(e x)} (list, unsorted)."""
    ph = sorted(Fraction(e) * x - int(Fraction(e) * x) for e in E)
    k = len(ph)
    gs = [ph[i + 1] - ph[i] for i in range(k - 1)]
    gs.append(1 - ph[-1] + ph[0])
    return gs


THETA = Fraction(1, 7)


def W_at(E, x):
    return sum((g - THETA) for g in gaps_at(E, x) if g > THETA)


def maxgap_at(E, x):
    return max(gaps_at(E, x))


def exact_intW_and_TV(E, cap_breaks=250000):
    """Exact Int_0^1 W and TV(W') by cell decomposition.
    Within a cell (between consecutive phase breakpoints) each gap is linear
    with integer slope; W = sum (g_i - theta)_+ is piecewise linear inside
    the cell with kinks where some g_i crosses theta. We subdivide at those
    crossings (exact rationals) and accumulate the integral of W and the
    jump-variation of W'."""
    bps = breakpoints(E)
    n = len(bps)
    if n > cap_breaks:
        return None, None  # too big for exact treatment here
    intW = Fraction(0)
    slopes_seq = []  # W' value on each (sub)cell, in order
    for idx in range(n):
        a = bps[idx]
        b = bps[idx + 1] if idx + 1 < n else Fraction(1)
        if b <= a:
            continue
        mid = (a + b) / 2
        # gap structure at the cell midpoint: identify adjacency & slopes
        ph = sorted((Fraction(e) * mid - int(Fraction(e) * mid), e) for e in E)
        k = len(ph)
        pairs = [(ph[i][1], ph[i + 1][1]) for i in range(k - 1)] + [(ph[-1][1], ph[0][1])]
        # gap g between tooth (e_i at pos p_i) and next tooth e_j:
        # g(x) = (e_j - e_i) x + const on the cell (mod-1 constants fixed inside)
        # compute exact g(a+)->g(b-) endpoints via values at a and b using
        # continuity: evaluate gap value at mid and use slope to extend.
        cell_kinks = {a, b}
        gap_data = []
        for (ei, ej) in pairs:
            slope = ej - ei if pairs[-1] != (ei, ej) or k == 1 else ej - ei
            # slope of frac(ej x) - frac(ei x) inside cell = ej - ei;
            # for the wraparound pair (last->first) g = 1 - span, slope = e_first - e_last
            pass
        # simpler: recompute each adjacent gap as linear function from value
        # at mid and slope = (e_next - e_this) [wrap: e_first - e_last]
        gmid = []
        for t in range(k - 1):
            gmid.append((ph[t + 1][0] - ph[t][0], ph[t + 1][1] - ph[t][1]))
        gmid.append((1 - ph[-1][0] + ph[0][0], ph[0][1] - ph[-1][1]))
        # kinks: g(x) = gm + s*(x - mid) crosses theta
        for (gm, s) in gmid:
            if s != 0:
                xc = mid + (THETA - gm) / s
                if a < xc < b:
                    cell_kinks.add(xc)
        pts = sorted(cell_kinks)
        for t in range(len(pts) - 1):
            u, v = pts[t], pts[t + 1]
            m2 = (u + v) / 2
            Wm, Sm = Fraction(0), 0
            for (gm, s) in gmid:
                gv = gm + s * (m2 - mid)
                if gv > THETA:
                    Wm += gv - THETA
                    Sm += s
            intW += Wm * (v - u)
            slopes_seq.append(Sm)
    # TV(W') = sum of |slope jumps| around the circle
    tv = 0
    m = len(slopes_seq)
    for i in range(m):
        tv += abs(slopes_seq[(i + 1) % m] - slopes_seq[i])
    return intW, tv


# ----------------------------------------------------------------------
# exact M(S) via local maxima of min_i ||v_i t|| at t = m/(v_i+v_j) and m/(2 v_i)
# ----------------------------------------------------------------------

def dist_num(v, num, den):
    """||v * num/den|| as Fraction, exact."""
    r = (v * num) % den
    return Fraction(min(r, den - r), den)


def M_exact(S, oneover14=Fraction(1, 14)):
    """exact M(S) = max_t min_i ||v_i t||; candidates at meeting points
    m/(v_i+v_j) (up-down crossings) and m/(2 v_i) (peaks). Early exit if a
    candidate certifies >= 1/14."""
    S = sorted(S)
    best = Fraction(0)
    dens = set()
    for i in range(len(S)):
        dens.add(2 * S[i])
        for j in range(i, len(S)):
            dens.add(S[i] + S[j])
    for D in sorted(dens):
        for m in range(1, D // 2 + 1):
            if gcd(m, D) != 1:
                continue
            val = min(dist_num(v, m, D) for v in S)
            if val > best:
                best = val
    return best


# ----------------------------------------------------------------------
# covering machinery
# ----------------------------------------------------------------------

QS = list(range(2, 15))


def is_covering(S):
    return all(any(s % q == 0 for s in S) for q in QS)


def missed_qs(P):
    return [q for q in QS if not any(p % q == 0 for p in P)]


# ----------------------------------------------------------------------
# the drift-embed exact hypothesis check on a grid point (Lean-faithful)
# ----------------------------------------------------------------------

def gp_slack_ok(P, x, V):
    """x in G_P with the snap/embed slack: nearInt(p x) >= 1/14 + p*(1)/V
    (drift of a P-runner across the full fast period)."""
    for p in P:
        xr = Fraction(p) * x
        d = xr - int(xr)
        d = min(d, 1 - d)
        if d < Fraction(1, 14) + Fraction(p, V):
            return False
    return True


def drift_embed_fires_at(E, P, V, j, spread):
    """EXACT hypothesis of minReach_ge_of_driftGap at grid point j/V, over all
    circular gaps of the teeth {frac(e j/V)}: exists gap (a, a+g) with
    1/7 + 2*spread*(a+g/2)/V < g, AND the P-part safe with slack at the
    realized tau = (j + a+g/2)/V (checked directly)."""
    x = Fraction(j, V)
    ph = sorted(Fraction(e) * x - int(Fraction(e) * x) for e in E)
    k = len(ph)
    cands = []
    for t in range(k - 1):
        a, g = ph[t], ph[t + 1] - ph[t]
        cands.append((a, g))
    # wraparound gap: from ph[-1] to 1+ph[0]; represent a = ph[-1], g = 1-ph[-1]+ph[0]
    cands.append((ph[-1], 1 - ph[-1] + ph[0]))
    for (a, g) in cands:
        if Fraction(1, 7) + 2 * (Fraction(spread) * (a + g / 2) / V) < g:
            # check P at realized tau
            tau = (Fraction(j) + a + g / 2) / V
            ok = True
            for p in P:
                d = Fraction(p) * tau
                d = d - int(d)
                d = min(d, 1 - d)
                if d < Fraction(1, 14):
                    ok = False
                    break
            if ok:
                return True
    return False


# ----------------------------------------------------------------------
# instance generation: k>=8 covering sets, three L-placement styles
# ----------------------------------------------------------------------

def gen_instance(P, V, style, k):
    """Build L (size k) covering Q_missed(P), speeds in (13, V], max = V.
    style: 'confined' -> members in [ceil(0.66 V), V]
           'uniform'  -> anywhere in (13, V]
           'slowheavy'-> force ~2 members in [14, V//3]"""
    Qm = missed_qs(P)
    lo_map = {'confined': max(14, (66 * V) // 100), 'uniform': 14, 'slowheavy': 14}
    lo = lo_map[style]
    L = {V}
    # covering duty first: for each missed q pick a multiple in [lo, V]
    for q in Qm:
        if any(u % q == 0 for u in L):
            continue
        mlo = (lo + q - 1) // q
        mhi = V // q
        if mlo > mhi:
            return None
        L.add(q * random.randint(mlo, mhi))
    # slow-heavy: inject slow members
    if style == 'slowheavy':
        for _ in range(2):
            if len(L) < k:
                L.add(random.randint(14, max(15, V // 3)))
    # fill to k
    tries = 0
    while len(L) < k and tries < 400:
        L.add(random.randint(lo, V))
        tries += 1
    if len(L) != k:
        return None
    S = sorted(set(P) | L)
    if len(S) != len(P) + k:
        return None
    if not is_covering(S):
        return None
    return sorted(L)


def analyze(P, L, do_exact_M=False, do_tv=True):
    V = max(L)
    E = sorted(V - u for u in L)
    spread = E[-1]  # = V - min L  (e=0 present)
    r = V / spread if spread > 0 else float('inf')
    slow = [u for u in L if u < Fraction(9, 14) * V]  # minL < (1-1/2.8)V ~ 0.643V
    s = len(slow)
    kt = len(L) - s
    Pt = len(P) + s
    # adaptive-split diagnosis
    if kt <= 7 and Pt <= 6:
        ladder = 'CLOSED(mu=1,Bonf)'
    elif kt <= 7:
        ladder = 'apex7-wall(G_Ptilde>0 needed, |Pt|=%d)' % Pt
    else:
        ladder = 'bars(k~=%d)+G_Ptilde(|Pt|=%d)' % (kt, Pt)
    out = dict(V=V, k=len(L), P=tuple(P), minL=min(L), spread=spread, r=r,
               s=s, ladder=ladder)
    if do_tv:
        intW, tv = exact_intW_and_TV(E)
        if intW is not None:
            out['intW'] = float(intW)
            out['TV'] = tv
            v0sq = Fraction(tv, 12) / intW if intW > 0 else None
            out['V0'] = float(v0sq) ** 0.5 if v0sq is not None else float('inf')
            out['aliasing_fires'] = (intW > 0) and (V * V > v0sq)
        else:
            out['aliasing_fires'] = None
    # drift-embed scan over grid points in G_P-with-slack
    fired = False
    for j in range(1, V):
        if not gp_slack_ok(P, Fraction(j, V), V):
            continue
        if drift_embed_fires_at(E, P, V, j, spread):
            fired = True
            break
    out['embed_fires'] = fired
    if do_exact_M:
        out['M'] = M_exact(sorted(set(P) | set(L)))
    return out


def main():
    print("=" * 78)
    print("HYP-5691 census: P∪L-split ratios + instrument firing map, k>=8 covering")
    print("=" * 78)

    # --- (Q1) both strata populated: explicit constructions -------------
    print("\n[Q1] Explicit k>=8 covering instances, both strata (exact):")
    demos = []
    # confined: P covers 2..6,8..13 densely; L covers {7,13,14,...} high up
    P1 = (8, 9, 10, 12)          # covers 2,3,4,5,6,8,9,10,12
    # missed: 7,11,13,14
    V1 = 1000
    L1 = [994, 990, 988, 986, 980, 973, 966, 952, 1000]  # 994=14*71,990=11*90,988=13*76...
    # ensure covering + k>=8: check
    S1 = sorted(set(P1) | set(L1))
    # slow-member instance:
    P2 = (8, 9, 10, 12)
    V2 = 1000
    L2 = [22, 994, 988, 700, 850, 920, 960, 980, 1000]   # 22 slow member (11|22), 994=14*71, 988=13*76
    for (P, L, name) in [(P1, L1, 'confined-L'), (P2, L2, 'slow-member-L')]:
        S = sorted(set(P) | set(L))
        cov = is_covering(S)
        E = sorted(max(L) - u for u in L)
        spread = E[-1]
        print(f"  {name}: P={P}, |L|={len(L)}, V={max(L)}, minL={min(L)}, "
              f"spread(E)={spread}, r={max(L)/spread:.2f}, covering={cov}")

    # --- (Q1/Q2) random census ------------------------------------------
    print("\n[Q2] Census (random, stratified by style); exact intW/TV where |E| breakpoints allow:")
    Ps = [(8, 9, 10, 12), (7, 9, 10, 11, 12), (11, 12, 13), (10, 11, 12, 13), (9, 11, 13)]
    rows = []
    for P in Ps:
        for style in ('confined', 'uniform', 'slowheavy'):
            got = 0
            attempts = 0
            while got < 4 and attempts < 60:
                attempts += 1
                V = random.choice([120, 150, 200, 260])
                k = random.choice([8, 9])
                L = gen_instance(list(P), V, style, k)
                if L is None:
                    continue
                # keep exact-TV tractable: spread cap for the tv computation
                E = sorted(max(L) - u for u in L)
                do_tv = (E[-1] <= 130)
                try:
                    row = analyze(list(P), L, do_exact_M=False, do_tv=do_tv)
                except Exception as ex:
                    continue
                row['style'] = style
                rows.append(row)
                got += 1
    # print stratified summary
    def strat(rows, pred):
        return [z for z in rows if pred(z)]
    conf = strat(rows, lambda z: z['r'] > 2.8)
    mid = strat(rows, lambda z: 1.41 < z['r'] <= 2.8)
    prop = strat(rows, lambda z: z['r'] <= 1.41)
    print(f"  total instances: {len(rows)}")
    for nm, grp in (('r>2.8 (confined)', conf), ('1.41<r<=2.8 (window)', mid), ('r<=1.41 (proportional)', prop)):
        if not grp:
            print(f"  {nm}: 0 instances")
            continue
        al = [z for z in grp if z.get('aliasing_fires') is True]
        alx = [z for z in grp if z.get('aliasing_fires') is not None]
        em = [z for z in grp if z['embed_fires']]
        print(f"  {nm}: n={len(grp)}, aliasing fires {len(al)}/{len(alx)} (exact-TV subset), "
              f"drift-embed fires {len(em)}/{len(grp)}")
    print("\n  per-instance detail (first 20):")
    for z in rows[:20]:
        print(f"   P={z['P']} k={z['k']} V={z['V']} minL={z['minL']} r={z['r']:.2f} s={z['s']} "
              f"style={z['style']} intW={z.get('intW','-')} TV={z.get('TV','-')} "
              f"V0={z.get('V0','-') if not isinstance(z.get('V0'),float) else round(z['V0'],1)} "
              f"alias={z.get('aliasing_fires','-')} embed={z['embed_fires']} ladder={z['ladder']}")

    # --- (Q3) adaptive-split ladder counts --------------------------------
    print("\n[Q3] Adaptive-split ladder strata over the census:")
    from collections import Counter
    cnt = Counter(z['ladder'] for z in rows)
    for lad, c in sorted(cnt.items()):
        print(f"   {lad}: {c}")

    # --- ground truth on a few small instances ---------------------------
    print("\n[Q4] Ground truth M(S) on small slow-member instances (exact):")
    checked = 0
    for z in rows:
        if z['s'] >= 1 and z['V'] <= 160 and checked < 4:
            pass
    # build fresh small ones for M check
    smalls = []
    for _ in range(200):
        P = random.choice(Ps)
        L = gen_instance(list(P), random.choice([60, 80, 100]), 'slowheavy', 8)
        if L is None:
            continue
        S = sorted(set(P) | set(L))
        if len(S) == len(P) + 8:
            smalls.append((P, L))
        if len(smalls) >= 4:
            break
    for (P, L) in smalls:
        S = sorted(set(P) | set(L))
        M = M_exact(S)
        E = sorted(max(L) - u for u in L)
        r = max(L) / E[-1]
        print(f"   S={S}: r={r:.2f} M={M} = {float(M):.5f} >= 1/14? {M >= Fraction(1,14)}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
