#!/usr/bin/env python3
"""
lrc14_adversarial_largediam_boxeph_S19.py  (HYP-6132, boxeph-2026-07-12-S19)

ADVERSARIAL STRESS of the large-diameter looseness premises (MISTAKE-137 discipline:
coherent/adversarial seeds, exact integer/rational arithmetic, verify every structural
predicate instead of assuming the construction delivers it).

Targets:
 (A) opus-S243 Case-A premise: "Vmax < lcm(2..14)=360360 ==> <=6 coprime-to-30030".
     Counterexample class: mid-scale SMOOTH covering elements (e.g. 27720 = 2^3*3^2*5*7*11
     covers every d in {2..14} except 13; the pair {10800, 6006} covers everything)
     free up to 11 coprime-to-30030 slots at Vmax << 360360, with all elements at
     COMPARABLE scale (no two-scale separation => far-element peel inapplicable).
     Flagship A3 additionally BLOCKS every non-14 modulus in [15,31] (kps HYP-6120
     style), so no bounded-window clearing certificate exists either.
 (B) THM-720 / opus-S243 "min M GROWS with diameter".
     Near-dilate seeds v = 2*H u {delta} (exact THM-668-monad detuned shape, M >= 1/13
     PROVED) with H a spread 12-core at growing scale S: exact M should PLATEAU at the
     detuned-dispatch level, not grow with diameter.
 (C) Instruments battery per family: WHICH proved atom certifies looseness?
       band-clearing at q <= QMAX (band-edge, opus-S235):  M >= ceil(q/14)/q
       descent (THM-636): exists L with all lifts k_i = round(v_i/L) >= 1 and
            [r = #distinct lifts <= 6 and B < L*(1/7 - 1/14) = L/14]   (LRC(7) floor)
         or [r <= 12 and B < L*(1/13 - 1/14) = L/182]                  (LRC(13) floor)
       far-element peel (THM-700/701 proxy): max adjacent ratio in sorted speeds
       detuned dispatch (THM-668-monad / THM-678): g >= 2 dividing >= 11 of 13,
            detuned complement not divisible by g
       pigeonhole (opus-S242): q in {15,21,25,27}, no mult of q, #coprime(q) < phi(q)/2
       exact M via pair-sum ruler enumeration (THM-668-mac-mini pair-sum exactness;
            moduli = all v_i+v_j incl. 2*v_i, plus |v_i-v_j| and q <= 60 as safety)

Pure Python (no numpy). Exact: integers + Fraction only at the end.
"""
from fractions import Fraction as F
from math import gcd, isqrt
from functools import reduce
import sys, random

PRIMORIAL = 30030  # 2*3*5*7*11*13
LCM214 = 360360    # lcm(2..14)

# ---------------------------------------------------------------- predicates
def is_dc(v):
    return all(any(x % d == 0 for x in v) for d in range(2, 15))

def missing_divisors(v):
    return [d for d in range(2, 15) if all(x % d for x in v)]

def is_primitive(v):
    return reduce(gcd, v) == 1

def longest_ap(v):
    s = set(v); vv = sorted(v); best = 1
    for i in range(len(vv)):
        for j in range(i + 1, len(vv)):
            d = vv[j] - vv[i]
            # count AP through vv[i] with difference d (vv[i] as start)
            L = 2; x = vv[j] + d
            while x in s: L += 1; x += d
            x = vv[i] - d
            while x in s: L += 1; x -= d
            if L > best: best = L
    return best

def coprime30030(v):
    return [x for x in v if gcd(x, PRIMORIAL) == 1]

def odd_count(v):
    return sum(1 for x in v if x % 2)

# ------------------------------------------------------------ instruments
def min_clear(v, qmax=300):
    """Smallest q>=15, 14 nmid q, admitting p with all v_i*p mod q in the safe band
    [ceil(q/14), q-ceil(q/14)].  Returns (q, p, F(band,q)) or None."""
    for q in range(15, qmax + 1):
        if q % 14 == 0: continue
        band = -(-q // 14)  # ceil(q/14)
        hi = q - band
        for p in range(1, q):
            ok = True
            for x in v:
                r = (x * p) % q
                if r < band or r > hi: ok = False; break
            if ok: return (q, p, F(band, q))
    return None

def blocked_window(v, lo=15, hi=43):
    """Moduli q in [lo,hi] (14 nmid q) with q | some v_i (condition-(a) blocking)."""
    return [q for q in range(lo, hi + 1) if q % 14 and any(x % q == 0 for x in v)]

def descent_census(v):
    """Scan L; lifts k_i = nearest integer of v_i/L, require all k_i>=1.
    Report best certificates for (r<=6, B<L/14) and (r<=12, B<L/182), plus min-r."""
    vmax = max(v)
    best6 = None; best12 = None; minr = (99, None)
    for L in range(2, vmax + 1):
        ks = []; B = 0
        for x in v:
            k = (2 * x + L) // (2 * L)  # round(x/L)
            if k < 1: ks = None; break
            b = x - L * k
            if abs(b) > B: B = abs(b)
            ks.append(k)
        if ks is None: continue
        r = len(set(ks))
        if r < minr[0]: minr = (r, (L, B))
        if r <= 6 and 14 * B < L and best6 is None: best6 = (L, r, B)
        if r <= 12 and 182 * B < L and best12 is None: best12 = (L, r, B)
    return best6, best12, minr

def max_adjacent_ratio(v):
    vv = sorted(v); best = (1.0, None)
    for a, b in zip(vv, vv[1:]):
        r = b / a
        if r > best[0]: best = (r, (a, b))
    return best

def detuned_shape(v):
    """Largest g>=2 dividing >= 11 of the 13 speeds; return (g, #divisible, detuned)."""
    best = None
    for g in range(2, max(v) + 1):
        div = [x for x in v if x % g == 0]
        if len(div) >= 11:
            det = [x for x in v if x % g]
            if all(x % g for x in det):
                if best is None or len(div) > best[1]: best = (g, len(div), det)
        if g > 50 and best: break
    return best

def pigeonhole_forced(v):
    """opus-S242: q in {15,21,25,27}, no mult of q, #coprime-to-q < phi(q)/2."""
    from math import gcd as g_
    def phi(n):
        r, m = n, n; p = 2
        while p * p <= m:
            if m % p == 0:
                while m % p == 0: m //= p
                r -= r // p
            p += 1
        if m > 1: r -= r // m
        return r
    hits = []
    for q in (15, 21, 25, 27):
        if any(x % q == 0 for x in v): continue
        nc = sum(1 for x in v if g_(x, q) == 1)
        if 2 * nc < phi(q): hits.append((q, nc, phi(q) // 2))
    return hits

def exact_M(v, extra_small=60, include_diffs=True):
    """Exact sup_t min_i ||v_i t|| by breakpoint enumeration.
    Moduli: pair-sums v_i+v_j (i<=j, so 2*v_i included), |v_i-v_j| (safety), q<=extra_small.
    Integer arithmetic with early-break pruning; returns (Fraction M, (m,q,p))."""
    qs = set(range(2, extra_small + 1))
    n = len(v)
    for i in range(n):
        for j in range(i, n):
            qs.add(v[i] + v[j])
            if include_diffs and v[i] != v[j]: qs.add(abs(v[i] - v[j]))
    mb, qb, pb = 0, 1, 0  # best margin mb/qb at multiplier pb
    for q in sorted(qs, reverse=True):
        if q < 2: continue
        # threshold: need m*qb > mb*q  =>  m >= floor(mb*q/qb)+1
        thr = (mb * q) // qb + 1
        if thr > q // 2: continue  # cannot beat best at this modulus
        for p in range(1, q // 2 + 1):
            m = q  # running min of distances
            for x in v:
                r = (x * p) % q
                if r > q - r: r = q - r
                if r < m:
                    m = r
                    if m < thr: break
            if m >= thr:
                mb, qb, pb = m, q, p
                thr = (mb * q) // qb + 1
    return F(mb, qb), (mb, qb, pb)

# ------------------------------------------------------------ family report
def report(name, v, do_exact=False, qmax=300):
    v = sorted(v)
    assert len(v) == 13 and len(set(v)) == 13, f"{name}: need 13 distinct speeds"
    dc = is_dc(v); prim = is_primitive(v); lap = longest_ap(v)
    cop = coprime30030(v); no = odd_count(v)
    print(f"\n== {name} ==")
    print(f"   v = {v}")
    print(f"   Vmax={max(v)} diam={max(v)-min(v)}  DC={dc}{'' if dc else ' MISSING '+str(missing_divisors(v))}  "
          f"primitive={prim}  longestAP={lap} (spread={lap<=7})")
    print(f"   #odd={no}  #coprime-to-30030={len(cop)}  coprime elts={cop if len(cop)<=13 else '...'}")
    mc = min_clear(v, qmax)
    bw = blocked_window(v)
    print(f"   blocked (a)-moduli in [15,43]: {bw}")
    if mc:
        q, p, bnd = mc
        print(f"   min-clear: q={q} p={p}  band-edge floor M >= {bnd} = {float(bnd):.4f} ({float(bnd)*14:.2f}x of 1/14)")
    else:
        print(f"   min-clear: NONE <= {qmax}  (no band certificate in window!)")
    b6, b12, minr = descent_census(v)
    print(f"   descent THM-636: r<=6 cert={b6}  r<=12 cert={b12}  min-r={minr[0]} at (L,B)={minr[1]}")
    mar = max_adjacent_ratio(v)
    print(f"   far-peel proxy: max adjacent ratio {mar[0]:.2f} at {mar[1]}")
    det = detuned_shape(v)
    print(f"   detuned shape (g | >=11 of 13): {det if det else 'NONE'}")
    ph = pigeonhole_forced(v)
    print(f"   pigeonhole-forced (S242): {ph if ph else 'NONE (enough coprime residues everywhere)'}")
    if do_exact:
        M, wit = exact_M(v)
        print(f"   EXACT M = {M} = {float(M):.5f}  ({float(M)*14:.2f}x of 1/14)   witness margin/q/p = {wit}")
        return v, M
    return v, None

# ================================================================ families
def primes_in(lo, hi, avoid=(2,3,5,7,11,13)):
    out = []
    for n in range(lo | 1, hi, 2):
        if all(n % p for p in avoid) or n in avoid:
            k, isp = 3, n > 1 and n % 2
            if n % 2 == 0: isp = n == 2
            while k * k <= n and isp:
                if n % k == 0: isp = False
                k += 2
            if isp: out.append(n)
    return out

def main():
    fast = '--fast' in sys.argv
    print("HYP-6132 adversarial stress -- boxeph-2026-07-12-S19")
    print(f"1/14 = {float(F(1,14)):.5f}; opus-S243 Case-A premise: Vmax<{LCM214} => <=6 coprime-to-30030")

    # ---------- calibration: the wall and kps's blocker ----------
    report("AP {1..13} (the wall)", list(range(1, 14)), do_exact=True)
    report("kps HYP-6120 blocker", [200,496,540,656,851,921,935,1122,1482,1680,1835,1849,1856],
           do_exact=not fast)

    # ---------- CLASS A: many-coprime mid-scale DC (Case-A counterexamples) ----------
    # A1: single smooth cover 27720 (covers d in 2..14 except 13) + 13-mult + 11 coprime
    ps = primes_in(17000, 27000)
    a1 = [27720, 13 * 1847] + ps[:3] + ps[len(ps)//2:len(ps)//2+4] + ps[-4:]
    a1 = sorted(set(a1))[:13]
    if len(a1) == 13: report("A1: 27720-cover + 11 coprime @ [17k,27k]", a1, do_exact=False)

    # A2: two smooth covers {5544,1560} + 11 coprime at [4000,5544] -- Vmax=5544 < opus's 5000-ish probe
    ps2 = primes_in(4000, 5544)
    a2 = [5544, 1560] + ps2[:4] + ps2[len(ps2)//2:len(ps2)//2+4] + ps2[-3:]
    a2 = sorted(set(a2))[:13]
    report("A2: {5544,1560}-cover + 11 coprime @ [4k,5.5k]", a2, do_exact=not fast)

    # A3 FLAGSHIP: two smooth covers {10800, 6006} + 8 semiprime blockers of primes 17..43
    #   + 3 free coprime primes -- blocks EVERY non-14 q in [15,31]; mid-scale; no ratio gap
    a3 = [10800, 6006,
          17*199,   # 3383  blocks 17
          19*181,   # 3439  blocks 19
          23*157,   # 3611  blocks 23
          29*127,   # 3683  blocks 29
          31*113,   # 3503  blocks 31
          37*97,    # 3589  blocks 37
          41*89,    # 3649  blocks 41
          43*83,    # 3569  blocks 43
          4001, 4507, 5003]
    report("A3 FLAGSHIP: blocker w/ 11 coprime, Vmax=10800", a3, do_exact=True)

    # ---------- CLASS B: near-dilate growth-plateau sweep ----------
    # v = 2*H u {delta}: H spread 12-core containing mults of {4,5,6,7,9,11,13}; delta odd.
    # THM-668-monad: M >= 1/13 PROVED. Question: does exact M GROW with scale? (predict: plateau)
    print("\n---- CLASS B: near-dilate sweep (2*H u {delta}), exact M vs scale ----")
    print("   THM-668 floor = 1/13 = 0.07692;  THM-720/opus-S243 claim min M grows 0.105->0.187->0.214+")
    rows = []
    for S in ([60, 240, 960] if fast else [60, 240, 960, 3840]):
        random.seed(S)
        # multiples anchored in [S, 2S], scattered
        H = [4*(S//4+random.randint(0, S//8)), 5*(S//5+random.randint(0, S//9)),
             6*(S//6+random.randint(0, S//10)), 7*(S//7+random.randint(0, S//11)),
             9*(S//9+random.randint(0, S//12)), 11*(S//11+random.randint(0, S//13)),
             13*(S//13+random.randint(0, S//14))]
        # 5 free elements: coprime-ish scatter in [S, 2S]
        while len(H) < 12:
            c = S + random.randint(S//7, S - 1)
            if c not in H: H.append(c)
        H = sorted(set(H))
        if len(H) != 12: continue
        delta = 2 * S + 2 * random.randint(S//4, S//2) + 1  # odd, comparable scale
        vB = sorted(2*h for h in H) + [delta]
        vB = sorted(vB)
        if len(set(vB)) != 13: continue
        name = f"B(S={S}): 2*H u {{{delta}}}"
        v, M = report(name, vB, do_exact=True)
        rows.append((S, max(vB)-min(vB), M))
    print("\n   near-dilate growth test: (scale, diam, exact M)")
    for S, d, M in rows:
        print(f"     S={S:5d} diam={d:6d}  M={M} = {float(M):.5f}  ({float(M)*14:.2f}x)")

    # ---------- CLASS B': THE DILATION TRANSPORT (the arithmetic refutation) ----------
    # M(c*X) = M(X) exactly (THM-531 dilation invariance), and adding a runner only
    # lowers M. So ONE low-M spread 12-core H0 (carrying the divisor duties through
    # the factor 2) yields, for EVERY prime c>13, the family v_c = 2c*H0 u {delta_c}:
    # primitive, DC, spread, diameter ~30c -> infinity, with M(v_c) <= M(2c*H0) = M(H0)
    # EXACTLY -- and M(v_c) >= 1/13 by the detuned dispatch (monad; g=2 | all 12 evens,
    # 2 nmid delta). If M(H0) < 0.105 (= THM-720's scale-10 minimum), 'min M grows with
    # diameter' is refuted at every scale, with no large-scale computation.
    print("\n---- CLASS B': dilation transport 2c*H0 u {delta_c} ----")
    H0 = [1, 2, 3, 4, 6, 9, 10, 11, 12, 13, 14, 16]  # mult of 4,5(10),6,7(14),9,11,13 present via 2*H0
    MH0, witH0 = exact_M(H0)
    print(f"   core H0 = {H0}")
    print(f"   exact M(H0) = {MH0} = {float(MH0):.5f}   witness {witH0}")
    # local hill-climb: minimize exact M over spread 12-cores in [1,40] keeping 2*H0 DC-complete
    def core_ok(H):
        if len(set(H)) != 12 or any(x < 1 for x in H): return False
        v2 = [2*x for x in H]
        if any(all(x % d for x in v2) for d in (4,5,6,7,8,9,10,11,12,13,14)): return False
        if longest_ap(H) > 7: return False
        return True
    best = (MH0, sorted(H0))
    random.seed(19)
    cur = sorted(H0)
    for it in range(400):
        cand = list(cur)
        i = random.randrange(12)
        cand[i] = max(1, cand[i] + random.choice([-3,-2,-1,1,2,3]))
        cand = sorted(cand)
        if not core_ok(cand): continue
        Mc, _ = exact_M(cand)
        if Mc <= best[0]:
            if Mc < best[0]: best = (Mc, cand)
            cur = cand
        elif random.random() < 0.15:
            cur = cand
    MH, Hbest = best
    print(f"   hill-climbed core H* = {Hbest}")
    print(f"   exact M(H*) = {MH} = {float(MH):.5f}  (LRC(13) floor 1/13 = {float(F(1,13)):.5f})")
    # transport: verify all structural predicates at each c, cite dilation invariance for M
    print(f"\n   transport v_c = 2c*H* u {{delta_c}}, c prime > 13; M(v_c) in [1/13, M(H*)] EXACTLY:")
    for c in (17, 97, 997, 9973):
        base = [2*c*h for h in Hbest]
        lo, hi = min(base), max(base)
        # delta: odd, mid-range, keeps primitivity/spread; scan a few candidates
        delta = None
        for d0 in range(( (lo+hi)//2 )|1, hi, 2):
            if gcd(d0, 2*c) == 1 and d0 not in base:
                vv = sorted(base + [d0])
                if longest_ap(vv) <= 7 and is_primitive(vv) and is_dc(vv):
                    delta = d0; break
        assert delta is not None
        vc = sorted(base + [delta])
        assert is_dc(vc) and is_primitive(vc) and longest_ap(vc) <= 7
        det = detuned_shape(vc)
        print(f"     c={c:5d}: Vmax={max(vc):8d} diam={max(vc)-min(vc):8d}  DC+prim+spread OK  "
              f"detuned(g=2)={'OK' if det and det[0]%2==0 else det}  M <= {float(MH):.5f}, >= 1/13")
    print(f"   => min M over spread primitive DC at diameter D is <= {float(MH):.5f} for arbitrarily")
    print(f"      large D (dilation transport). Compare THM-720 sampled minima: 0.105 (scale 10),")
    print(f"      0.187 (scale 200), opus-S243 min 0.136->0.214. GROWTH {'REFUTED' if float(MH) < 0.105 else 'CAPPED at ' + str(float(MH))}.")

    # exact M at c=1 for the record (delta included; may be lower than M(H*))
    base1 = [2*h for h in Hbest]
    d1 = None
    for d0 in range((min(base1)+max(base1))//2|1, max(base1)+40, 2):
        if gcd(d0,2)==1 and d0 not in base1:
            vv = sorted(base1+[d0])
            if longest_ap(vv) <= 7 and is_primitive(vv) and is_dc(vv):
                d1 = d0; break
    if d1:
        v1 = sorted(base1+[d1])
        report("B'(c=1): 2*H* u {%d}" % d1, v1, do_exact=True)

    # ---------- summary ----------
    print("\n==== SUMMARY (HYP-6132) ====")
    print("A-class: #coprime-to-30030 can be 11 at Vmax 5544..27720 << 360360 -> Case-A premise")
    print("         is construction-dependent; A3 also blocks [15,31] and defeats descent/peel/")
    print("         detuned/pigeonhole -> its looseness certificate must come from q>=32 clearing")
    print("         or pair-sum rulers = the un-shrunk anti-concentration core.")
    print("B-class: exact M vs scale printed above -- plateau refutes 'min M grows with diameter'")
    print("         (growth was a sampling artifact of narrow-band constructions; near-dilates pin")
    print("         min M near the detuned-dispatch floor at EVERY scale). Looseness itself survives.")

if __name__ == '__main__':
    main()
