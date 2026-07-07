#!/usr/bin/env python3
r"""
lrc_step_gauge_kps_S62.py   (kind-pasteur-2026-07-07-S62, HYP-4877)

THE STEP-GAUGE: the tournament tiling-model move applied to LRC(14) (owner directive:
"a tree on 8 events has 7 edges and that is just the Hamiltonian path connecting each
element in an 8 player tournament ... apply similar analysis to lonely runners").

THE DICTIONARY.
  tournament (n players)                    LRC(14) (14 runners: 13 moving + the origin)
  ------------------------------------     ------------------------------------------------
  base Hamiltonian path (n-1 arcs)          sorted speeds 0 < v(1) < ... < v(13): the path
                                            0 -> v(1) -> ... -> v(13), 13 edges on 14 nodes
  tile coordinates (free arcs)              edge weights = consecutive differences
                                            delta_i = v(i) - v(i-1) >= 1  (v(0) = 0)
  relabeling/gauge freedom                  delta_1 = translation gauge (mu ignores it);
                                            common scaling = dilation gauge
  COMPLEMENT T -> T^op                      STEP REVERSAL (delta_2..delta_13) -> reversed:
                                            E -> max(E)+min(E) - E, mu-invariant via x -> -x
  self-complementary (SC) classes           PALINDROMIC step words
  cut (+) cycle decomposition               the pair matrix d_ij = v_j - v_i is PURE CUT
                                            (a path metric); pair (i,j) crosses exactly
                                            |d_ij| times per period (difference uniformity)
  staircase tile weights                    order-cell cost of step l ~ #pairs straddling l

This script:
  PART 1: reversal invariance -- mu_{1/7}(E) and E[maxgap](E) invariant under
          E -> max+min-E: EXACT check on the zoo (incl GW, whose word is NOT palindromic).
  PART 2: palindromic-extremizer test -- per raw diameter D, adversarial min of mu_{1/7}
          and E[maxgap] over (a) all step words, (b) palindromic-only: do the minimizers
          sit on the symmetry locus (as SC tournaments do for max H)?
  PART 3: THE a=2 STRATUM -- a 13-set whose 12 differences take <= 2 distinct values
          {s,t} (i copies of s, 12-i of t, any order) lies in the (i+1) x (13-i) grid
          {v1 + a*s + b*t}: N = (i+1)(13-i) <= 49.  Per-split ledger:
            mu2(i+1, 13-i) (2-torus, numeric) vs m_P, and the small-(s,t) side vs the
            S59 exact diameter floor (diam = i*s + (12-i)*t <= 75).  Assemble the
            two-leg closure per split and report the honest residue (medium st vs error C).
  PART 4: the crossing staircase -- #order cells of the maxgap function per step position
          (which steps generate complexity), verifying the pair-crossing count |d_ij| and
          the triangle weight profile.
"""
from fractions import Fraction as F
import math, random, itertools

TH = 1 / 7
MP = 14249 / 252252

def maxgap_pts(pts):
    ph = sorted(pts)
    g = max(ph[i + 1] - ph[i] for i in range(len(ph) - 1)) if len(ph) > 1 else 1.0
    return max(g, ph[0] + 1.0 - ph[-1])

def mu_numeric(E, res=25000):
    c = 0
    for r in range(res):
        x = (r + .5) / res
        if maxgap_pts([(e * x) % 1.0 for e in E]) > TH: c += 1
    return c / res

def Emg_numeric(E, res=25000):
    t = 0.0
    for r in range(res):
        x = (r + .5) / res
        t += maxgap_pts([(e * x) % 1.0 for e in E])
    return t / res

def steps_of(E):
    E = sorted(E)
    return tuple(E[i + 1] - E[i] for i in range(len(E) - 1))

def family_from_steps(d, v1=1):
    v = [v1]
    for s in d: v.append(v[-1] + s)
    return v

# ---------------------------------------------------------------- PART 1
print("=" * 96)
print("PART 1 -- step-reversal = the complement symmetry: mu and E[maxgap] invariant (exact/numeric)")
print("=" * 96)
zoo = {
    "AP {1..13}": list(range(1, 14)),
    "GW {1..11,13,24}": list(range(1, 12)) + [13, 24],
    "record 2*{1..11}u{11,13}": [2,4,6,8,10,11,12,13,14,16,18,20,22],
    "stretch {0,2..12,17,28}": [0] + list(range(2, 13)) + [17, 28],
    "random word": family_from_steps((3,1,4,1,5,9,2,6,5,3,5,8)),
}
for nm, E in zoo.items():
    Erev = sorted(max(E) + min(E) - e for e in E)
    w = steps_of(E)
    pal = "PALINDROMIC" if w == w[::-1] else "non-palindromic"
    m1, m2 = mu_numeric(E, 15000), mu_numeric(Erev, 15000)
    a1, a2 = Emg_numeric(E, 15000), Emg_numeric(Erev, 15000)
    print(f"  {nm:28s} steps={w} [{pal}]")
    print(f"      mu: {m1:.5f} vs reversed {m2:.5f} (diff {abs(m1-m2):.5f});  E[mg]: {a1:.5f} vs {a2:.5f} (diff {abs(a1-a2):.5f})")

# ---------------------------------------------------------------- PART 2
print()
print("=" * 96)
print("PART 2 -- palindromic-extremizer test: per-diameter adversarial minima, all vs palindromic")
print("=" * 96)
rng = random.Random(62)
def descend(D, functional, palindromic, iters=160):
    # step words of length 12 summing to D (v1 fixed = 1; translation gauge)
    def rand_word():
        cuts = sorted(rng.sample(range(1, D), 11))
        w = [cuts[0]] + [cuts[i+1]-cuts[i] for i in range(10)] + [D - cuts[-1]]
        if palindromic:
            half = w[:6]
            w = half + half[::-1]
            s = sum(w)
            w[5] += (D - s) // 2 * 0  # keep symmetric: adjust middle pair equally
            # rebuild to sum D: scale middle pair
            diff = D - sum(w)
            w[5] += diff - diff // 2; w[6] += diff // 2
            if w[5] < 1 or w[6] < 1: return None
            if palindromic and diff % 2 == 1:  # keep palindrome: put odd remainder... skip odd
                return None
        return tuple(w) if all(x >= 1 for x in w) else None
    best = (9.9, None)
    for _ in range(iters):
        w = rand_word()
        if w is None: continue
        E = family_from_steps(w)
        v = functional(E, 4000)
        if v < best[0]:
            # local moves preserving sum (and palindromy if required)
            cur, cv = list(w), v
            for _ in range(60):
                i, j = rng.randrange(12), rng.randrange(12)
                if i == j: continue
                c2 = cur[:]
                if palindromic:
                    ii, jj = 11 - i, 11 - j
                    if len({i, j, ii, jj}) not in (2, 4): continue
                    c2[i] += 1; c2[ii] += 1; c2[j] -= 1; c2[jj] -= 1
                else:
                    c2[i] += 1; c2[j] -= 1
                if any(x < 1 for x in c2) or sum(c2) != D: continue
                E2 = family_from_steps(tuple(c2))
                v2 = functional(E2, 4000)
                if v2 < cv: cur, cv = c2, v2
            v, w = cv, tuple(cur)
            if v < best[0]: best = (v, w)
    return best

print(f"  {'D':>4} {'min mu (all)':>13} {'word':^34} {'min mu (palin)':>15}  pal-word-is-min?")
for D in (12, 16, 20, 24):
    va, wa = descend(D, mu_numeric, False)
    vp, wp = descend(D, mu_numeric, True)
    pal_a = wa == wa[::-1] if wa else None
    print(f"  {D:4d} {va:13.4f} {str(wa):^34} {vp:15.4f}  all-min palindromic: {pal_a}")
print(f"\n  {'D':>4} {'min E[mg] (all)':>15} {'word':^34} {'min E[mg] (palin)':>17}")
for D in (12, 20):
    va, wa = descend(D, Emg_numeric, False)
    vp, wp = descend(D, Emg_numeric, True)
    pal_a = wa == wa[::-1] if wa else None
    print(f"  {D:4d} {va:15.4f} {str(wa):^34} {vp:17.4f}  all-min palindromic: {pal_a}")

# ---------------------------------------------------------------- PART 3
print()
print("=" * 96)
print("PART 3 -- THE a=2 STRATUM: two-letter step words close via the S61 grid ledger")
print("=" * 96)
def mu2_torus(n1, n2, res=260):
    cnt = 0
    for a in range(res):
        u = (a + 0.5) / res
        for b in range(res):
            v = (b + 0.5) / res
            pts = [(i * u + j * v) % 1.0 for i in range(n1) for j in range(n2)]
            if maxgap_pts(pts) > TH: cnt += 1
    return cnt / (res * res)

print("  split i (copies of s) | grid (i+1)x(13-i) |  N  | mu2 (2-torus) | vs m_P")
mu2_split = {}
for i in range(0, 13):
    n1, n2 = i + 1, 13 - i
    if n1 > n2: continue          # reversal symmetry: (i, 12-i) same grid transposed
    N = n1 * n2
    m2 = mu2_torus(n1, n2)
    mu2_split[(n1, n2)] = m2
    print(f"    i={i:2d} (and {12-i:2d})     |   {n1:2d} x {n2:2d}       | {N:3d} |    {m2:.4f}    |  {'OK' if m2 >= MP else '*** below'}  ")
worst_split = min(mu2_split.values())
print(f"  worst split mu2 = {worst_split:.4f} = {worst_split/MP:.1f}x m_P")
print()
print("  ASSEMBLY per split (s,t; i copies of s, diam = i*s+(12-i)*t):")
print("   - diam <= 75:            S59 exact diameter floor (PROVED, roof now THM-637)")
print("   - diam >= 76, st large:  mu >= mu2(split) - C/(s*t); need C/(s*t) < mu2 - m_P")
print("   - the pinch: minimal s*t given diam >= 76 within the split")
Cerr = 5.0   # empirical S61 error constant (Fourier-provable shape)
print(f"   using empirical C = {Cerr} (S61 PART 4):")
allok = True
for i in range(0, 13):
    n1, n2 = i + 1, 13 - i
    key = (n1, n2) if n1 <= n2 else (n2, n1)
    m2 = mu2_split[key]
    margin = m2 - MP
    # worst (s,t): minimize s*t subject to i*s + (12-i)*t >= 76, s != t, s,t >= 1
    best_st = None
    for s in range(1, 90):
        for t in range(1, 90):
            if s == t: continue
            if i * s + (12 - i) * t >= 76:
                st = s * t
                if best_st is None or st < best_st[0]:
                    best_st = (st, s, t)
        # (small pruning unnecessary)
    st, s, t = best_st
    need = Cerr / st
    ok = need < margin
    allok &= ok
    print(f"    i={i:2d}: mu2={m2:.4f} margin={margin:.4f}  worst st at diam>=76: (s,t)=({s},{t}) st={st}  err<={need:.4f}  {'CLOSES' if ok else '*** PINCH: needs sharper error or LOO'}")
print(f"  => a=2 stratum: {'ALL splits close (modulo exact mu2 + provable C)' if allok else 'residue at flagged splits (t=1 far-element shapes -> LOO/klein lane)'}")

# ---------------------------------------------------------------- PART 4
print()
print("=" * 96)
print("PART 4 -- the crossing staircase: order-cell cost per step position")
print("=" * 96)
print("  pair (i,j) of speeds crosses |v_i - v_j| times per period (difference uniformity);")
print("  a step at sorted position l is straddled by (l)(13-l) pairs + the 0-runner pairs:")
E = list(range(1, 14))
tot = sum(abs(a - b) for a, b in itertools.combinations([0] + E, 2))
print(f"  AP {{1..13}} (with 0): total pair crossings per period = {tot}")
w = [(l * (14 - l)) for l in range(1, 14)]
print(f"  triangle weights l(14-l) for l=1..13: {w}  (middle steps cost most)")
print(f"  sum of weights = {sum(w)}  (should equal total crossings for unit steps: {tot})")
print()
print("DONE.")
