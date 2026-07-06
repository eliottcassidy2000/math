#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S3 -- HYP-4252: THE UNIFORM CELL LEMMA + the next-four-cells
anchored probe.

PART A (verification of the lemma, 0 violations expected): for ANY finite speed
set W with M(W) = c/q (reduced) attained at t* = m/q* (reduced):
  (0) q | q*   [gcd(c,q)=1: binder grid distance q*c/q is an integer]
  (i) every binder v (d_{q*}(vm) = c q*/q = ck) satisfies k | v, k = q*/q
  (ii) scaling: d_{qk}(k x) = k d_q(x)  => quotient binders sit at exactly c mod q
  (iii) both binder signs occur (kink), and q | v' + w' for every up/down
        quotient binder pair (v = w allowed: the single-runner peak edge case)
  (iv) NEW -- BINDER UNITS: gcd(v', q) = 1 for every quotient binder
       [p | gcd(v',q) => p | c, contradicting gcd(c,q) = 1]
  (v) NEW -- WITNESS DETERMINISM: m == +- c * (v')^{-1} (mod q) per binder.
The lemma never uses in-gapness: every attained M of every family is a test.

PART B (the anchored realization probe, next four cells + 3/38 cross-check):
cells (3,38), (4,51), (5,63), (5,64), (6,77).  For each unit pair shape
{a, q-a} (gcd(a,q)=1, a < q/2 -- phi(q)/2 shapes): families containing the pair
as actual elements, a covering core (multiples of 2..12, chosen to clear level
c mod q at the FORCED witness dilation m = c a^{-1}), free slots filled with
random clearers.  Exact M on near-gap candidates.  Expected: ZERO families with
M in (1/13, 2/25); recorded: the quantization gap min |M - c/q| per cell.
"""
import random, itertools, math
from fractions import Fraction as F
from math import gcd

# ---------------------------------------------------------------- exact M + argmax
def dist_int(x, q):
    """distance of x to 0 in Z/q"""
    r = x % q
    return min(r, q - r)

def exact_M_argmax(W):
    """Exact M(W) = max_t min_v ||vt|| and one maximizer, via the critical grid:
    t = j/(v+w), j/(v-w) (v!=w), and half-integer peaks (2v)."""
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = F(0); best_t = F(0)
    seen = set()
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = F(j, s)          # auto-reduces; non-reduced j/s are VALID candidates
            if t in seen:
                continue
            seen.add(t)
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best, best_t = mv, t
    return best, best_t

def float_M(W, grid_mult=6):
    """Float lower-envelope max over the same critical grid (fast pre-pass)."""
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    best = 0.0
    for s in dens:
        for j in range(1, s):
            t = j / s
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best

# ---------------------------------------------------------------- PART A
def check_lemma(W, verbose=False):
    """Returns list of violation strings for the uniform cell lemma on W."""
    M, t = exact_M_argmax(W)
    if M == 0:
        return []
    c, q = M.numerator, M.denominator
    m, qs = t.numerator, t.denominator
    viol = []
    # (0) q | q*
    if qs % q != 0:
        return [f"(0) q={q} does not divide q*={qs} (M={M}, t={t}, W={W})"]
    k = qs // q
    ck = c * k
    binders = [v for v in W if dist_int(v * m, qs) == ck]
    if not binders:
        return [f"(binder-empty) M={M} t={t} W={W}"]
    ups   = [v for v in binders if (v * m) % qs == ck]
    downs = [v for v in binders if (v * m) % qs == qs - ck]
    # (iii) both signs (v=w allowed when ck == qs-ck i.e. 2ck=qs)
    if 2 * ck == qs:
        ups = downs = binders
    if not ups or not downs:
        viol.append(f"(iii-kink) missing a binder sign: M={M} t={t} W={W}")
    for v in binders:
        # (i) k | v
        if v % k != 0:
            viol.append(f"(i) binder {v} not divisible by k={k} (M={M}, t={t}, W={W})")
            continue
        vp = v // k
        # (ii) quotient sits at exactly c mod q
        if dist_int(vp * m, q) != c:
            viol.append(f"(ii) quotient binder {vp} not at level c={c} mod q={q} (W={W})")
        # (iv) binder units
        if gcd(vp, q) != 1:
            viol.append(f"(iv) quotient binder {vp} not a unit mod q={q} (M={M}, W={W})")
        else:
            # (v) witness determinism: m == +- c * vp^{-1} mod q
            inv = pow(vp, -1, q)
            if (m - c * inv) % q != 0 and (m + c * inv) % q != 0:
                viol.append(f"(v) m={m} != +-c*inv({vp}) mod {q} (M={M}, W={W})")
    # (iii) pair divisibility on quotients
    if ups and downs and k > 0:
        for v in ups:
            for w in downs:
                if v % k or w % k:
                    continue
                if (v // k + w // k) % q != 0:
                    viol.append(f"(iii) quotient pair {v//k}+{w//k} not 0 mod q={q} (W={W})")
    return viol

def part_A():
    print("=" * 78)
    print("PART A: uniform cell lemma vs real attainers (0 violations expected)")
    print("=" * 78)
    random.seed(3)
    total, viols, cells_seen = 0, [], {}
    # structured + random families, sizes 3..12, heights <= 60
    fams = []
    for n in (3, 4, 5, 6, 8, 10, 12):
        for _ in range(60):
            W = sorted(random.sample(range(1, 61), n))
            g = math.gcd(*W) if n > 1 else W[0]
            W = [v // g for v in W] if g > 1 else W
            fams.append(sorted(set(W)))
    # structured: APs, lifts, dilates, pair-anchored
    for d in range(1, 6):
        fams.append([d * i for i in range(1, 13)])
    fams.append([1,2,3,5,7,8,9,10,11,12,17,19])   # the 2/25 attainer
    fams.append(list(range(1, 12)) + [168])        # the deep well
    for a in (1, 3, 5, 7, 9, 11, 13, 15, 17):      # 38-pair anchored smalls
        base = [a, 38 - a] + random.sample([v for v in range(1, 40) if v not in (a, 38-a)], 6)
        fams.append(sorted(set(base)))
    for W in fams:
        if len(W) < 2:
            continue
        total += 1
        vs = check_lemma(W)
        if vs:
            viols.extend(vs)
        M, t = exact_M_argmax(W)
        if M > 0:
            cells_seen[(M.numerator, M.denominator)] = cells_seen.get((M.numerator, M.denominator), 0) + 1
    print(f"  families tested: {total}; DISTINCT attained cells (c,q): {len(cells_seen)}")
    ks = {}
    for W in fams[:0]:
        pass
    if viols:
        print(f"  VIOLATIONS ({len(viols)}):")
        for v in viols[:20]:
            print("   ", v)
    else:
        print("  ZERO violations of (0),(i),(ii),(iii),(iv),(v) -- lemma verified")
    # how often is k > 1?  (the stratification is real, not vacuous)
    kcount = {}
    for W in fams:
        if len(W) < 2:
            continue
        M, t = exact_M_argmax(W)
        if M == 0:
            continue
        k = t.denominator // M.denominator
        kcount[k] = kcount.get(k, 0) + 1
    print(f"  k-spectrum over attainers: {dict(sorted(kcount.items()))}")
    return len(viols)

# ---------------------------------------------------------------- PART B
CELLS = [(3, 38), (4, 51), (5, 63), (5, 64), (6, 77)]
LO, HI = F(1, 13), F(2, 25)

def clearers(q, m, c, lo, hi):
    """values v in [lo,hi] with d_q(v m) >= c"""
    return [v for v in range(lo, hi + 1) if dist_int(v * m, q) >= c]

def probe_cell(c, q, n_per_shape=120, hmax=90, seed=7):
    random.seed(seed + q)
    shapes = [a for a in range(1, q // 2 + (q % 2)) if gcd(a, q) == 1]
    in_gap_hits = []
    best_gap = (float('inf'), None, None)   # (|M - c/q| float, W, M)
    below = at_or_above = under_lo = 0
    tested = 0
    target = c / q
    for a in shapes:
        pair = (a, q - a)
        for sgn in (1, -1):
            m = (sgn * c * pow(a, -1, q)) % q
            # covering core: for each d in 2..12 a multiple of d clearing level c at m
            cores = []
            ok = True
            for d in range(2, 13):
                cand = [d * h for h in range(1, hmax // d + 1) if dist_int(d * h * m, q) >= c]
                if not cand:
                    ok = False
                    break
                cores.append(cand)
            if not ok:
                continue
            pool = clearers(q, m, c, 1, hmax)
            pool = [v for v in pool if v not in pair]
            if len(pool) < 10:
                continue
            for _ in range(n_per_shape // 2):
                W = set(pair)
                for cand in cores:
                    W.add(random.choice(cand))
                while len(W) < 12:
                    W.add(random.choice(pool))
                if len(W) != 12:
                    continue
                W = sorted(W)
                tested += 1
                fM = float_M(W)
                if fM < float(LO) - 2e-9:
                    under_lo += 1
                elif fM > float(HI) + 2e-9:
                    at_or_above += 1
                else:
                    # near/inside gap at float precision -> exact check
                    M, t = exact_M_argmax(W)
                    if LO < M < HI:
                        in_gap_hits.append((W, M, t))
                    elif M <= LO:
                        under_lo += 1
                    else:
                        at_or_above += 1
                d = abs(fM - target)
                if d < best_gap[0]:
                    best_gap = (d, W, fM)
    return dict(shapes=len(shapes), tested=tested, in_gap=in_gap_hits,
                best=best_gap, under_lo=under_lo, above=at_or_above)

def part_B():
    print()
    print("=" * 78)
    print("PART B: anchored realization probe -- cells (3,38),(4,51),(5,63),(5,64),(6,77)")
    print("=" * 78)
    grand_hits = 0
    for c, q in CELLS:
        r = probe_cell(c, q)
        tgt = F(c, q)
        print(f"  cell {c}/{q} (= {float(tgt):.6f}): shapes phi(q)/2 = {r['shapes']}, "
              f"families tested = {r['tested']}")
        print(f"    M < 1/13: {r['under_lo']};  M >= 2/25: {r['above']};  "
              f"IN-GAP: {len(r['in_gap'])}")
        d, W, fM = r['best']
        print(f"    nearest approach to {c}/{q}: |M - c/q| ~ {d:.6f} at M ~ {fM:.6f}")
        if r['in_gap']:
            grand_hits += len(r['in_gap'])
            for W, M, t in r['in_gap'][:5]:
                print(f"    *** IN-GAP HIT: W={W}, M={M}, t={t}")
    print()
    if grand_hits == 0:
        print("  ZERO in-gap realizations across all five cells -- the void holds "
              "on every anchored template probed.")
    else:
        print(f"  *** {grand_hits} IN-GAP HITS -- (G) IS FALSE OR A BUG; INVESTIGATE")
    return grand_hits

if __name__ == "__main__":
    va = part_A()
    hb = part_B()
    print()
    print(f"SUMMARY: partA violations = {va}; partB in-gap hits = {hb}")
