#!/usr/bin/env python3
"""
LRC(14) A3 - SDR / SECTION-SHUFFLING CENSUS.  kind-pasteur S2-wf.

How do shared sections index the hard families?

Setup (matches lrc14_sdr_hall_symmetry / regions_grounding):
  At grid time tau=a/14 (a in (Z/14)*={1,3,5,9,11,13}), runner i sits in
  SECTION r_i = v_i*a mod 14.  Observer LONELY at a/14  <=>  no r_i == 0.
  A 13-runner config's ON-GRID picture is a function of its RESIDUE MULTISET
  R = { v_i mod 14 }.  (Z/14)* acts on R by entrywise multiplication.

PARTS
  1. Census orbits of 13-element residue-multisets over Z/14 under (Z/14)*.
     Classify by SECTION-SHARING PATTERN = the multiset of multiplicities on
     the NONZERO sections {1..13} plus the count in section 0.  Perfect SDR
     {1..13} = the maximally-spread fixed shape (pattern all-singletons, 0 empty).
  2. EASY->HARD degradation. For real speed-sets (not just residues): how the
     section-sharing pattern tracks hardness.  Is "more sharing" monotone "harder"?
     Compute exact M and meas(G_S) (the lonely set measure) on a graded family.
  3. EASY/HARD correspondence via covering families {C U {84m}}.  Does (Z/14)*
     send covering->covering and core->core?  Does each hard family correspond
     to a fixed easy section pattern of its core?
  4. Hall: trivial mod 14 (distinct&nonzero).  Does a FINER grid restore a
     non-trivial matching/Hall obstruction for hard cores?

Stdlib only. Exact rationals.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
from collections import Counter, defaultdict

N = 14

# ---------------- exact M tool (verbatim, validated) ----------------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def g(S, t):
    return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2):
            C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at
def binding(S, t):
    return sorted(v for v in S if nrm(v*t) == g(S, t))

units = [a for a in range(1, N) if gcd(a, N) == 1]   # {1,3,5,9,11,13}

# ============================================================
print("="*74)
print("PART 1 - CENSUS of 13-elt residue-multisets over Z/14 under (Z/14)*")
print("="*74)
print("A 13-runner config's on-grid picture depends only on its residue MULTISET")
print("R = multiset of (v_i mod 14), |R|=13.  (Z/14)* acts by entrywise *a.")
print("SECTION-SHARING PATTERN of R: count in section 0 (=z) + the partition of")
print("the remaining 13-z runners by which nonzero section they share.")
print()

def share_pattern(R):
    """R: tuple of 13 residues in 0..13.  Return (z, sorted multiplicity partition
       on nonzero sections, #distinct-nonzero-sections-used, #empty-nonzero-sections)."""
    c = Counter(R)
    z = c.get(0, 0)
    nonzero_mults = sorted((m for s, m in c.items() if s != 0), reverse=True)
    used = len([s for s in c if s != 0])
    empty = 13 - used    # nonzero sections {1..13} not occupied
    return (z, tuple(nonzero_mults), used, empty)

def mult_by(R, a):
    return tuple(sorted((r*a) % N for r in R))

def orbit(R):
    return frozenset(mult_by(R, a) for a in units)

# Full census of ALL 13-multisets over Z/14 is C(13+13,13)=10400600 - too big to
# enumerate naively but we don't need every multiset; we need the ORBIT TYPES,
# which are determined by the share-pattern up to the unit action.  KEY FACT:
# multiplying by a unit a is a BIJECTION of Z/14 fixing 0, so it permutes the
# nonzero sections and FIXES the multiplicity-partition and the count z.
# => the SECTION-SHARING PATTERN is a (Z/14)*-INVARIANT.  We verify this, then
# census the patterns (= partitions of 13 into z + parts), and for each ask:
# is it grid-loneable (z==0)? is it perfect SDR (z==0 & all parts ==1)?

print("[1a] CLAIM: the section-sharing pattern (z, mult-partition) is (Z/14)*-INVARIANT.")
print("     (mult by a unit fixes 0, permutes {1..13}, preserves multiplicities.)")
import random
random.seed(3)
ok = True
for _ in range(20000):
    R = tuple(random.randint(0, 13) for _ in range(13))
    p0 = share_pattern(R)
    for a in units:
        if share_pattern(mult_by(R, a)) != p0:
            ok = False; break
    if not ok: break
print("     verified invariant on 20000 random 13-multisets:", ok)
print()

# Census the share-PATTERNS (these are the orbit *types*).  A pattern is:
#   z in 0..13 ; a partition of (13 - z) into parts, with #parts <= 13 (only 13
#   nonzero sections available).  Loneable <=> z==0.  Perfect SDR <=> z==0 and
#   partition = all 1's (13 singletons fill sections 1..13 exactly).
def partitions(n, maxpart=None, maxlen=None):
    if maxpart is None: maxpart = n
    if maxlen is None: maxlen = n
    if n == 0:
        yield (); return
    if maxlen == 0: return
    for first in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - first, first, maxlen - 1):
            yield (first,) + rest

patterns = []
for z in range(0, 14):
    rem = 13 - z
    for part in partitions(rem, maxpart=rem, maxlen=13):
        patterns.append((z, part))
total = len(patterns)
loneable = [p for p in patterns if p[0] == 0]
sdr = [p for p in patterns if p[0] == 0 and all(x == 1 for x in p[1])]
print("[1b] Number of distinct section-sharing PATTERNS (orbit types):", total)
print("     of these, GRID-LONEABLE (z=0, section 0 empty):", len(loneable))
print("     of those, PERFECT SDR (z=0 & all-singletons):", len(sdr),
      "-> UNIQUE pattern:", sdr[0] if sdr else None)
print()
print("     The perfect SDR is the SINGLE maximally-spread loneable pattern:")
print("       z=0, partition=(1^13) - every nonzero section used exactly once.")
print("     ALL OTHER loneable patterns SHARE at least one section (some part>=2)")
print("     and leave at least one nonzero section EMPTY.")
print()
# show the loneable patterns graded by "sharing", measured by max multiplicity
print("[1c] Loneable patterns graded by sharing (max nonzero multiplicity):")
by_maxmult = defaultdict(list)
for (z, part) in loneable:
    by_maxmult[max(part)].append(part)
for mm in sorted(by_maxmult):
    ps = by_maxmult[mm]
    egs = ps[:3]
    print("     maxmult=%2d : %4d patterns   e.g. %s" % (mm, len(ps), egs))
print("     SDR is exactly maxmult=1 (the unique least-shared loneable pattern).")
print()

# ============================================================
print("="*74)
print("PART 2 - EASY->HARD degradation: does more section-sharing => harder?")
print("="*74)
print("We need ACTUAL speed sets (off-grid M sees integer speeds, not just residues).")
print("Take a graded family that degrades {1..13} (perfect SDR) toward a covering set")
print("by pushing residues to collide / land on section 0, and measure M + meas(G_S).")
print()

# meas(G_S): Lebesgue measure of {tau in [0,1) : ||v tau||>=1/14 for all v in S}.
# G_S is a union of open intervals; its endpoints are among the binding-pair
# crossings where some ||v tau|| == 1/14 exactly, i.e. v*tau = k +/- 1/14.
# Collect all such breakpoints, sort, and on each subinterval test the midpoint.
def measG(S, thr=F(1, 14)):
    S = sorted(set(S))
    bps = {F(0), F(1)}
    for v in S:
        # ||v tau|| = thr  <=>  v tau = k + thr or k - thr, for integer k
        # tau in [0,1): k ranges so that tau in [0,1)
        kmax = v  # generous
        for k in range(-1, kmax + 2):
            for off in (thr, -thr):
                tau = F(k, 1) + off
                tau = tau / v
                if F(0) <= tau < F(1):
                    bps.add(tau)
    bp = sorted(bps)
    meas = F(0)
    for i in range(len(bp) - 1):
        lo, hi = bp[i], bp[i+1]
        if hi <= lo: continue
        mid = (lo + hi) / 2
        if g(S, mid) >= thr:
            meas += (hi - lo)
    return meas

def section_share_score(S, n=N):
    """At a representative unit (a=1), how shared are the sections?
       Return (z=count in sec0, used distinct nonzero sections, empty nonzero,
       sharing-entropy)."""
    R = tuple(v % n for v in S)
    z, part, used, empty = share_pattern(R)
    # entropy of the occupancy distribution over nonzero sections (base-2)
    import math
    tot = len(R) - z
    H = 0.0
    if tot > 0:
        for m in part:
            p = m / tot
            H -= p * math.log2(p)
    return z, used, empty, part, H

print("Graded family (13 speeds), degrading the AP {1..13}:")
fam = []
fam.append(("perfect SDR (AP)        ", list(range(1, 14))))
fam.append(("one collision (swap 13->1 dup)", [1,2,3,4,5,6,7,8,9,10,11,12,15]))  # 15==1 mod14 collides w/1
fam.append(("two collisions          ", [1,2,3,4,5,6,7,8,9,10,11,15,16]))         # 15==1,16==2
fam.append(("Goddyn-Wong {1..11,13,24}", list(range(1,12))+[13,24]))               # 24==10 mod14 (dup of 10)
fam.append(("covering core {1..11,13,84}", list(range(1,12))+[13,84]))            # 84==0 mod14 (sec 0!)
fam.append(("covering {1..11,13,168}  ", list(range(1,12))+[13,168]))             # 168==0 mod14 (m=2)
for name, S in fam:
    z, used, empty, part, Hent = section_share_score(S)
    m, at = M(S)
    mg = measG(S)
    nonzero = all(v % N != 0 for v in S)
    print("  %s" % name)
    print("     residues mod14 = %s" % sorted(v % N for v in S))
    print("     z(sec0)=%d  distinct-nonzero-sections=%d  empty=%d  share-part=%s  entropy=%.3f"
          % (z, used, empty, part, Hent))
    print("     EXACT M = %s = %.6f  (>1/14? %s)   meas(G_S)=%s=%.6f"
          % (m, float(m), m > F(1,14), mg, float(mg)))
    print("     binding at tau*=%s : %s" % (at, binding(S, at)))
    print()

print("Reading (CORRECTED - the naive 'more sharing = harder' is FALSE):")
print(" * Within z=0 (no runner parked in section 0): MORE section-sharing makes")
print("   M LARGER (EASIER), not smaller.  Perfect SDR {1..13} = the AP has the")
print("   SMALLEST gap M=1/14 (it IS the tight extremal); one collision -> 1/13,")
print("   two collisions -> 1/12.  So the maximally-spread SDR is the HARDEST z=0")
print("   config, and sharing sections RELAXES it.  Section-sharing is NOT a")
print("   monotone proxy for 'harder' inside the easy (z=0) regime - it is the")
print("   REVERSE.  The AP/SDR sits at the bottom of the well.")
print(" * The covering hard families live at z=1 (one runner in section 0): there")
print("   grid loneliness DIES (gridM=0) and the true M is OFF-grid.  Crucially")
print("   M(S_m)=7m/(84m+5) is INCREASING in m: 7/89 < 14/173 < 21/257 -> 1/12.")
print("   So m=1 (the LEAST-degraded covering set, parked runner closest to 0 in")
print("   absolute speed) is the HARDEST, approaching 1/14 from above but never")
print("   reaching it.  z=1 with the SMALLEST off-grid parked runner = the worst.")
print(" * NET: the relevant 'hardness coordinate' is NOT the share-partition; it is")
print("   (a) whether section 0 is occupied (z>=1 forces off-grid), and (b) HOW the")
print("   parked runner's true position degrades.  Section-sharing on the nonzero")
print("   sections is essentially orthogonal to M in this family.")
print()

# ---- Part 2b: robustness of the reversal + meas(G_S) trend in the covering family
print("="*74)
print("PART 2b - ROBUSTNESS: is the 'sharing relaxes M' reversal a fluke? + meas trend")
print("="*74)
print("Sweep collisions 0..6 in the AP (replace the top residues by +14 duplicates),")
print("staying z=0, and confirm M is NON-DECREASING as sharing increases:")
prevM = None
mono = True
for k in range(0, 7):
    # take {1..13}, push the top k speeds up by +14 so they duplicate residues 1..k
    base = list(range(1, 14 - k))            # 1..(13-k)
    dups = [14 + i for i in range(1, k + 1)] # 15,16,... duplicate residues 1,2,...
    S = base + dups
    while len(S) < 13:                        # pad if needed (k small keeps len 13)
        S.append(S[-1] + 14)
    S = S[:13]
    z, used, empty, part, Hent = section_share_score(S)
    m, at = M(S)
    mg = measG(S)
    flag = ""
    if prevM is not None:
        if m < prevM: mono = False; flag = " <-- M DECREASED"
    print("  collisions=%d : residues=%s  z=%d empty=%d  M=%s=%.6f  meas=%.6f%s"
          % (k, sorted(v % N for v in S), z, empty, m, float(m), float(mg), flag))
    prevM = m
print("  M non-decreasing across the WHOLE sweep?", mono, "(not strictly - see k=6)")
print("  But M(SDR)=1/14 is the GLOBAL MINIMUM of the sweep, and every sharing step")
print("  k=1..5 RAISES M.  So: spreading (SDR) is the tight floor; sharing generally")
print("  relaxes M but not strictly monotonically in collision-count (k=6 dips to")
print("  2/21=0.0952 from 1/10).  HONEST: SDR=min, sharing helps, not a clean monotone.")
print()
print("meas(G_S) across the covering family S_m (the gap-1/14 lonely-set measure):")
Ccore = list(range(1, 12)) + [13]
for m in (1, 2, 3, 4):
    S = Ccore + [84*m]
    mm, at = M(S)
    mg = measG(S)
    print("  m=%d: M=%s=%.6f  meas(G_S)=%s=%.6f  (M>1/14? %s)"
          % (m, mm, float(mm), mg, float(mg), mm > F(1,14)))
print("  meas(G_S) GROWS with m (less-tight as the parked runner gets large/slow),")
print("  while M also grows toward 1/12: the m=1 set is the unique hardest, both in")
print("  smallest M (7/89) and smallest positive meas (563/105105).")
print()

# ============================================================
print("="*74)
print("PART 3 - EASY/HARD correspondence via covering families {C U {84m}}")
print("="*74)
print("User's claim: each HARD family corresponds (via section-shuffling) to a")
print("specific EASY section pattern of its 12-speed CORE C.")
print()
print("[3a] Covering family: S_m = {1..11,13} U {84m}.  84=lcm(12,14), so 84m==0 mod14")
print("     -> 84m sits in section 0 at EVERY unit (gridM=0).  Core C={1..11,13}.")
print("     Closed form (given): M(S_m)=7m/(84m+5).  Check exact + section picture:")
C = list(range(1, 12)) + [13]
print("     Core C residues mod14:", sorted(v % N for v in C),
      " (distinct & nonzero -> C is grid-SDR-able: an EASY core, by LRC(12))")
zC, usedC, emptyC, partC, HC = section_share_score(C)
print("     Core share-pattern: z=%d used=%d empty=%d part=%s (perfect SDR on its 12 residues)"
      % (zC, usedC, emptyC, partC))
print("     Core empty nonzero sections (besides 0):",
      sorted(set(range(1,14)) - set(v % N for v in C)),
      " <- section 12 is the lone hole the parked runner would want but can't use on-grid")
print()
for m in (1, 2, 3):
    S = C + [84*m]
    mm, at = M(S)
    pred = F(7*m, 84*m + 5)
    print("     m=%d: S=...{84*%d=%d}  M=%s (pred %s) match=%s  tau*=%s binding=%s"
          % (m, m, 84*m, mm, pred, mm == pred, at, binding(S, at)))
print()

print("[3b] Does (Z/14)* act compatibly: unit*S covering => covering, core->core?")
print("     Multiply the WHOLE set by a unit a (mod nothing - actual integer scaling")
print("     does NOT preserve covering; the RESIDUE action does).  Test the residue")
print("     action: 'covering' = contains a multiple of every q in 2..14.")
def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))
print("     covering(S_1)?", is_covering(C + [84]))
print("     Apply residue-unit a to mod-14 picture only (84m stays ==0):")
for a in units:
    R = tuple((v*a) % N for v in (C + [84]))
    # the section-0 occupant stays section 0 (0*a=0); core residues permute among 1..13
    z = R.count(0)
    print("       a=%2d: parked-in-sec0 count=%d  core sections=%s"
          % (a, z, sorted(set(R)-{0})))
print("     => the parked runner is FIXED in section 0 by every unit (0 is a fixed")
print("        point of mult), and the core's 12 sections are PERMUTED among {1..13}.")
print("        So (Z/14)* sends covering-core-config to covering-core-config: the")
print("        easy core stays an SDR (just relabeled sections), the hard runner stays")
print("        parked.  The action is COMPATIBLE with the easy/hard split.")
print()
print("[3c] Which EASY pattern does the hard family correspond to?  The core C is a")
print("     PERFECT SDR on 12 of the 13 nonzero sections, leaving EXACTLY ONE nonzero")
print("     section empty (the hole).  The hard family is obtained by adding a runner")
print("     that - instead of filling that empty section - is forced into section 0")
print("     (because 84m==0).  So: HARD = (easy 12-SDR with one hole) + (runner that")
print("     misses the hole and lands on the observer).  The 'shared section' that")
print("     indexes the family is section 0 itself (shared by observer + parked runner).")
print()

# Which nonzero section does the core leave empty, and what residue would fill it?
holes = sorted(set(range(1,14)) - set(v % N for v in C))
print("     Core C={1..11,13}: empty nonzero section(s)=%s." % holes)
print("     residue 12 is the missing one (12 not in {1..11,13} mod14).  A runner ==12")
print("     mod14 would COMPLETE the SDR (easy).  Instead 84m==0 -> parked.  The hard")
print("     family is the 'failed completion': the 13th runner aims at the hole(12)")
print("     but is dragged to 0.  DIFFERENT hard families <-> DIFFERENT holes/cores.")
print()
# enumerate: for each 12-subset-core that is an SDR with exactly one hole h, the
# matching hard family parks a runner in 0 instead of filling h.
print("[3d] Census: easy 12-cores that are perfect-SDR-with-one-hole, by hole h:")
# A 12-core SDR uses 12 distinct nonzero residues; the hole h is the unused one.
# (Z/14)* permutes the holes among {1..13}; count orbits of holes.
hole_orbit = {}
seen=set()
for h in range(1,14):
    if h in seen: continue
    orb = sorted(set((h*a)%N for a in units))
    for x in orb: seen.add(x)
    hole_orbit[h]=orb
print("     (Z/14)* orbits on possible holes {1..13}:")
for h,orb in hole_orbit.items():
    print("       orbit of hole %2d: %s (size %d)" % (h, orb, len(orb)))
print("     => holes fall into orbits; section 7 (the antipode) is its own orbit")
print("        (7*a==7 mod14 for all units), so a core whose hole is 7 is (Z/14)*-fixed.")
print()

# ============================================================
print("="*74)
print("PART 4 - Hall obstruction: trivial mod 14, restored on a finer grid?")
print("="*74)
print("[4a] mod 14: SDR-able <=> residues distinct & nonzero.  Bipartite graph")
print("     runners->sections is a single permutation (unit action), so Hall is")
print("     automatic or fails by a raw collision.  NO non-trivial deficiency. (S1)")
print()
print("[4b] FINER grid for a hard core.  S_1={1..11,13,84}, optimum off-grid at")
print("     tau*=37/89 (denominator 89=84+5).  Test whether on the mod-Q circle the")
print("     hard core admits a perfect SDR (distinct & nonzero residues mod Q):")
hard = list(range(1,12)) + [13, 84]
def sdr_units_modQ(S, Q):
    res=[]
    for a in range(1, Q):
        if gcd(a, Q) != 1: continue
        rr = [(v*a) % Q for v in S]
        if 0 not in rr and len(set(rr)) == len(rr):
            res.append(a)
    return res
for Q in (14, 28, 84, 89, 178):
    su = sdr_units_modQ(hard, Q)
    tot_units = sum(1 for a in range(1,Q) if gcd(a,Q)==1)
    print("     mod %3d: # units giving a perfect SDR = %d / %d units  (Q factors: %s)"
          % (Q, len(su), tot_units, Counter()))
print()
print("[4c] WHY mod 89 is the right finer grid (84 is invertible mod 89):")
print("     84 mod 89 =", 84 % 89, "; gcd(84,89)=", gcd(84,89), "(89 prime) -> 84 is a UNIT")
print("     so the parked runner is NO LONGER stuck in section 0 on the mod-89 circle.")
print("     A perfect SDR exists mod 89 -> the off-grid lonely time lives there.")
print()
print("[4d] Is there a genuine Hall/matching obstruction on the finer grid, or still")
print("     trivial?  On mod Q with all v invertible, mult-by-a is again a permutation")
print("     of Z/Q, so {v_i a} distinct <=> {v_i} distinct mod Q, and nonzero <=> v_i!=0.")
print("     => STILL trivial (distinct & nonzero mod Q).  The matching is forced; Hall")
print("     never bites.  BUT a real obstruction appears when we demand the SDR sections")
print("     avoid not just {0} but the whole FORBIDDEN BAND |section|<thr*Q around 0.")
print("     That is the true (off-grid) lonely condition.  Model it:")
print()

# The true lonely condition at modulus Q with threshold 1/14: tau=a/Q, runner v
# is >=1/14 from 0 iff its section s=(v a mod Q) satisfies thr*Q <= s <= (1-thr)*Q,
# i.e. s avoids the BAND B = {0,1,...,floor(Q/14)-?} and its mirror near Q.
# This is a REAL SDR-with-forbidden-zone (a "list" / Hall-with-deficiency) problem.
def band_forbidden(Q, thr=F(1,14)):
    # sections s in 0..Q-1 with ||s/Q|| < thr are forbidden (runner too close to 0)
    forb = set()
    for s in range(Q):
        if nrm(F(s, Q)) < thr:
            forb.add(s)
    return forb

print("[4e] BAND model on mod 89 for hard core S_1 (the FORBIDDEN ZONE, not just {0}):")
Q = 89
forb = band_forbidden(Q)
print("     forbidden sections (||s/89||<1/14):", sorted(forb))
print("     A lonely grid-89 time a needs EVERY runner's section v*a mod89 OUTSIDE")
print("     the band.  This is NOT a matching condition (runners may share sections")
print("     freely as long as none lands in the band) -> it's a COVERING-avoidance,")
print("     same shape as the mod-14 'avoid section 0' but with a thick band.")
goodA=[]
for a in range(1, Q):
    if gcd(a, Q) != 1: continue
    secs = [(v*a) % Q for v in hard]
    if all(s not in forb for s in secs):
        goodA.append(a)
print("     # units a (mod89) with ALL runners out of the band:", len(goodA),
      "e.g.", goodA[:8])
if goodA:
    a0 = goodA[0]
    print("     check a=%d: tau=%d/89, g(S,tau)=%s (>=1/14? %s)"
          % (a0, a0, g(hard, F(a0,89)), g(hard, F(a0,89)) >= F(1,14)))
m, at = M(hard)
print("     EXACT off-grid M(S_1)=%s=%.6f at tau*=%s ; is tau* an a/89? %s"
      % (m, float(m), at, at.denominator == 89))
print()
print("="*74)
print("PART 5 - VERDICT")
print("="*74)
print("""
CENSUS (Part 1):
 * The on-grid picture of a 13-runner config depends only on its residue
   multiset R mod 14, and the SECTION-SHARING PATTERN (z = #parked in section 0,
   plus the multiplicity partition on nonzero sections) is a (Z/14)*-INVARIANT.
 * Among ALL section-sharing patterns, PERFECT SDR = the UNIQUE loneable pattern
   with no sharing (z=0, partition 1^13).  {1..13} realizes it.  Every other
   loneable pattern shares >=1 section and leaves >=1 nonzero section empty.

MONOTONICITY (Part 2) - the naive proxy is REVERSED, be honest:
 * 'More section-sharing => harder' is FALSE.  Within z=0, the MAXIMALLY-SPREAD
   perfect SDR {1..13} (the AP) is the HARDEST z=0 config (M=1/14 exactly, the
   tight wall); sharing sections (collisions) RAISES M (1/13, 1/12, ...).  The
   AP/SDR sits at the BOTTOM of the well, not the top.  Section-sharing on the
   NONZERO sections is essentially orthogonal to / anti-correlated with hardness.
 * The real hardness coordinate is z (section-0 occupancy).  z=0 => M>=1/13
   (LRC(13)) so NOT a counterexample.  z>=1 (a runner == 0 mod 14) is the only
   way to a covering set; then gridM=0 and M is forced OFF-grid.  Among the
   covering family S_m={1..11,13,84m}, M(S_m)=7m/(84m+5) is INCREASING in m, so
   m=1 (M=7/89, the least-degraded) is the HARDEST and the closest approach to
   1/14 from above.  All of these are still > 1/14 (consistent with LRC(14)).
 * Tightness (M=1/14 exactly) is NOT unique to the SDR: Goddyn-Wong {1..11,13,24}
   ALSO hits M=1/14 (with a shared section, residue 10 doubled).  So 'tight' is a
   small finite locus {AP, Goddyn-Wong, ...}, not 'the SDR'.  meas(G_S)=0 exactly
   at these tight sets (the gap-1/14 lonely set is a measure-zero point set), and
   meas(G_S)>0 for every covering S_m (e.g. 563/105105 at m=1) - matching codex's
   meas(G_C) sharpening.

EASY/HARD CORRESPONDENCE (Part 3):
 * CONFIRMED in the precise sense: the hard family {C U {84m}} = (easy 12-core C
   that is a perfect SDR leaving exactly ONE nonzero section empty) + (a 13th
   runner forced into section 0 instead of filling that hole).  The 'shared
   section' indexing the family is section 0 (observer + parked runner).
 * (Z/14)* acts COMPATIBLY: 0 is a fixed point of multiplication, so the parked
   runner stays parked under every unit while the core's SDR is merely relabeled
   among sections {1..13}.  Covering->covering, core->core.  Holes fall into
   (Z/14)* orbits; the antipodal hole 7 gives a unit-FIXED core.

HALL (Part 4):
 * mod 14: trivial (distinct & nonzero), confirming HYP-2571b.  No Hall bite.
 * FINER grid mod 89 (84 is a unit mod 89): the parked runner is FREED, a perfect
   SDR exists, and the off-grid optimum lives there (tau*=37/89).  But the SDR
   condition is STILL trivial (distinct & nonzero mod 89) because mult-by-a is a
   permutation.  The genuine obstruction is NOT a matching/Hall condition: it is
   a BAND-AVOIDANCE (every section must miss the thick forbidden zone |s/Q|<1/14
   around 0), which is the same covering-avoidance shape as mod 14, just thicker.
   There is NO hidden Hall theorem; the hardness is covering, not matching.
""")
