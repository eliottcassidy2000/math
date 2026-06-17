#!/usr/bin/env python3
"""
ANGLE A: formalize and STRESS-TEST the hard<-easy reduction for LRC(14).

Claim (THM-523 + hand-in-hand mechanism): a covering 13-set S that contains a
"parked" runner w (a multiple of 14, the perfect-middle reversal fixed point) is
LRC(14) (M(S) >= 1/14) because:
  - removing w gives a core A = S\{w} that FAILS to cover some modulus q<=13
    (w was the only multiple of q in S),
  - so A is EASY: tau=1/q witnesses M(A) >= 1/q,
  - re-adding the perfect-middle runner w only DIPS the gap by a resonance amount,
  - if DIP = M(A)-M(S) <= SLACK = 1/q - 1/14, then M(S) >= 1/14.

We test DIP <= SLACK across many covering 13-sets. We report the worst
DIP/SLACK ratio and whether it ever reaches >= 1 (which would break the reduction).
We also confirm M(S)=M(A) for generic m (parked runner = 14m) and dips only at
resonances.

stdlib only.
"""
from fractions import Fraction as F
from itertools import combinations

# ---------- EXACT GAP TOOL ----------
def nrm(x):
    r = x - int(x)
    r = r + 1 if r < 0 else r
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
        if v > b:
            b = v; at = t
    return b, at

# ---------- helpers ----------
def covers(S, q):
    """Does S contain a multiple of q?"""
    return any(v % q == 0 for v in S)

def uncovered_moduli(S, qmax=14):
    return [q for q in range(2, qmax+1) if not covers(S, q)]

def smallest_uncovered(S, qmax=13):
    u = uncovered_moduli(S, qmax)
    return u[0] if u else None

ONE14 = F(1, 14)

# ---------- PART 1: generic vs resonant m for a fixed easy 12-core ----------
def part1():
    print("="*78)
    print("PART 1: M(A u {14m}) vs M(A) for an easy 12-core A, sweep m.")
    print("="*78)
    # A is an easy 12-core that does NOT cover q=12 (no multiple of 12).
    # Standard worst core from the prompt: A = {1,...,12} dropped... actually
    # {1..12} DOES contain 12. To make q uncovered we want a core that misses a
    # multiple of some modulus. Use A = {1,2,3,4,5,6,7,8,9,10,11,13} (drop 12):
    # then q=12 uncovered. Also test A={1..12} where parked runner is lcm.
    cores = {
        "drop-12 core {1..11,13}": [1,2,3,4,5,6,7,8,9,10,11,13],
        "{1..12}": list(range(1,13)),
    }
    for name, A in cores.items():
        MA, tA = M(A)
        q = smallest_uncovered(A, 13)
        print(f"\nCore {name}: M(A)={MA} ({float(MA):.5f}), smallest uncovered q={q}")
        if q is not None:
            print(f"   easy-witness 1/q = {F(1,q)} ({float(F(1,q)):.5f}); slack 1/q-1/14 = {F(1,q)-ONE14}")
        results = []
        for m in range(1, 40):
            w = 14*m
            S = A + [w]
            MS, tS = M(S)
            dip = MA - MS
            results.append((m, w, MS, dip))
        # show which m are resonant (dip>0)
        nonres = [r for r in results if r[3] == 0]
        res    = [r for r in results if r[3] != 0]
        print(f"   generic (dip=0): {len(nonres)}/{len(results)} of m in 1..39")
        print(f"   resonant (dip>0): m = {[r[0] for r in res]}")
        for (m,w,MS,dip) in res[:12]:
            sl = (F(1,q)-ONE14) if q else None
            ratio = (dip/sl) if (sl and sl>0) else None
            print(f"     m={m:>2} w={w:>4}: M(S)={MS} ({float(MS):.5f}) dip={dip} "
                  f"(={float(dip):.5f}) dip/slack={float(ratio):.4f}" if ratio is not None
                  else f"     m={m:>2} w={w:>4}: M(S)={MS} dip={dip}")
        # family min over m
        fammin = min(r[2] for r in results)
        print(f"   FAMILY MIN over m=1..39: M = {fammin} ({float(fammin):.6f}); "
              f">= 1/14 ? {fammin>=ONE14}")

# ---------- PART 2: stress test over MANY covering 13-sets ----------
def gen_covering_13sets():
    """
    Generate covering 13-sets S that contain a parked runner w (multiple of 14),
    built as A u {w} where A is a 12-set and w=14m. We want S to COVER all of
    2..14 while A FAILS to cover some q (w is the only multiple of q).
    We construct structured families plus random-ish small covering cores.
    """
    sets = []
    # Family 1: A = drop-one-residue cores {1..14}\{j} restricted to 12 elements,
    # plus parked runner. Build A from a 12-subset of {1..13}, parked = 14*m.
    base = list(range(1,14))  # 1..13
    # The 13 "drop one" 12-cores of {1..13}:
    for drop in base:
        A = [x for x in base if x != drop]   # 12 elements
        for m in range(1, 13):
            w = 14*m
            S = sorted(set(A + [w]))
            if len(S) != 13:   # parked runner could coincide if w in A (never, w>13)
                continue
            # require S covers all of 2..14
            if uncovered_moduli(S, 14):
                continue
            # require A fails some q (so A is easy)
            if smallest_uncovered(A, 13) is None:
                continue
            sets.append(("drop%d+14*%d"%(drop,m), tuple(A), tuple(S)))
    # Family 2: structured small covering sets. Take small covering 12-cores
    # that miss exactly one modulus q, parked = 14*m chosen to cover q and 14.
    # Cores missing q: choose A = small multiples covering 2..14 except q.
    # Use the {1..12} core (misses nothing under 13? it covers 2..12,
    # but misses 13 and 14). Parked w=14m must supply 13? no. So {1..12}+14m
    # never covers 13 -> not a covering set. Skip; rely on Family 1.
    # Family 3: prime-power covering cores. A built from {2,3,4,7,8,9,11,13, ...}
    extra_cores = [
        [2,3,4,5,6,7,8,9,10,11,12,13],   # covers 2..13, misses 14 alone
        [1,2,3,4,5,6,7,8,9,10,11,13],    # drop 12 -> misses 12
        [1,2,3,4,5,6,7,8,9,10,12,13],    # drop 11 -> misses 11
        [1,2,3,4,5,6,7,8,9,11,12,13],    # drop 10
        [2,4,6,8,10,12,3,9,5,11,7,13],   # mixed
    ]
    for A in extra_cores:
        A = sorted(set(A))
        if len(A) != 12: continue
        for m in range(1, 13):
            w = 14*m
            S = sorted(set(A + [w]))
            if len(S) != 13: continue
            if uncovered_moduli(S, 14): continue
            if smallest_uncovered(A, 13) is None: continue
            sets.append(("core%s+14*%d"%(A[0:1],m), tuple(A), tuple(S)))
    return sets

def part2():
    print("\n" + "="*78)
    print("PART 2: STRESS TEST DIP <= SLACK over covering 13-sets S = A u {14m}.")
    print("="*78)
    sets = gen_covering_13sets()
    # dedupe by S
    seen = {}
    for name, A, S in sets:
        if S not in seen:
            seen[S] = (name, A)
    print(f"distinct covering 13-sets generated: {len(seen)}")

    worst_ratio = F(0); worst_ex = None
    broke = []  # cases where M(S) < 1/14 (would actually break LRC14)
    ratio_ge_1 = []  # cases dip >= slack
    n_tested = 0
    minMS = None; minMS_ex = None
    for S,(name,A) in seen.items():
        S = list(S); A = list(A)
        MS, tS = M(S)
        MA, tA = M(A)
        dip = MA - MS
        q = smallest_uncovered(A, 13)        # largest slack via smallest q
        if q is None:
            continue
        slack = F(1,q) - ONE14
        n_tested += 1
        if minMS is None or MS < minMS:
            minMS = MS; minMS_ex = (name, S, MS, tS)
        # the reduction needs dip <= slack
        if slack > 0:
            ratio = dip / slack
        else:
            ratio = F(0) if dip <= 0 else F(10**9)  # slack==0 means q=14, no help
        if ratio > worst_ratio:
            worst_ratio = ratio
            worst_ex = (name, S, A, MS, MA, dip, q, slack, ratio)
        if dip >= slack and slack >= 0:
            ratio_ge_1.append((name, S, A, MS, MA, dip, q, slack))
        if MS < ONE14:
            broke.append((name, S, MS, tS))

    print(f"\ntested {n_tested} covering 13-sets with an easy core")
    print(f"GLOBAL MIN M(S) = {minMS} ({float(minMS):.6f}) at {minMS_ex[0]}, S={minMS_ex[1]}")
    print(f"  >= 1/14 ? {minMS >= ONE14}   (1/14 = {float(ONE14):.6f})")
    print(f"\nWORST DIP/SLACK ratio = {worst_ratio} ({float(worst_ratio):.6f})")
    if worst_ex:
        name,S,A,MS,MA,dip,q,slack,ratio = worst_ex
        print(f"  at {name}")
        print(f"  S={S}")
        print(f"  A={A}")
        print(f"  M(A)={MA} ({float(MA):.6f}), M(S)={MS} ({float(MS):.6f})")
        print(f"  smallest uncovered q={q}, 1/q={F(1,q)}, slack={slack} ({float(slack):.6f})")
        print(f"  DIP={dip} ({float(dip):.6f})  DIP/SLACK={float(ratio):.6f}")
    print(f"\ncases with DIP >= SLACK (reduction-as-stated fails to certify): {len(ratio_ge_1)}")
    for c in ratio_ge_1[:10]:
        name,S,A,MS,MA,dip,q,slack = c
        print(f"  {name}: M(S)={MS} M(A)={MA} dip={dip} q={q} slack={slack} "
              f"(but M(S)>=1/14? {MS>=ONE14})")
    print(f"\ncases that ACTUALLY break LRC14 (M(S) < 1/14): {len(broke)}")
    for c in broke[:10]:
        print(f"  {c[0]}: S={c[1]} M(S)={c[2]} ({float(c[2]):.6f}) at tau={c[3]}")

    return worst_ratio, worst_ex, broke, minMS, minMS_ex

# ---------- PART 3: honesty check -- does the reduction-as-stated always certify? ----------
def part3(worst_ratio, worst_ex):
    print("\n" + "="*78)
    print("PART 3: HONESTY -- when DIP > SLACK, is M(S) still >= 1/14 (by other means)?")
    print("="*78)
    # The reduction "DIP <= SLACK => M(S) >= 1/14" is a SUFFICIENT certificate.
    # If DIP > SLACK it does NOT prove M(S) < 1/14; M(S) could still be >= 1/14.
    # We already collected ratio_ge_1 cases in part2 and checked M(S)>=1/14 there.
    if worst_ratio >= 1:
        print("Worst ratio reached >= 1: the SIMPLE reduction does NOT certify all sets.")
    else:
        print(f"Worst ratio = {float(worst_ratio):.6f} < 1: DIP <= SLACK on EVERY tested set.")
        print("=> The hand-in-hand reduction certifies M(S) >= 1/14 for all tested covering 13-sets.")

if __name__ == "__main__":
    part1()
    wr, we, broke, minMS, minMS_ex = part2()
    part3(wr, we)
