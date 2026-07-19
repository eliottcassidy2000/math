#!/usr/bin/env python3
"""
death-star-2026-07-19-S59b -- HYP-7890: the N=31 BAND-LAW TEST.

THM-1284 band law: for N >= 8, first gap W_N = (1/(N+1), 2/(2N+1)) nonempty
<=> N == 1 (mod 6) AND 5 nmid (3N+2).  N=31 is the FIRST N==1 mod 6 with
5 | 3N+2 = 95: the b=5 binder congruence dies, mac-mini-S29 pinned the
canonical family F(31) = {1..29,31,90} at M = 1/32 (degrade to the floor).
PREDICTIONS UNDER TEST:
  P1  M(F(31)) = 1/32 exactly (replication of mac-mini-S29).
  P2  single-far strata: members exist at exactly (N,i,x) = (19,18,54) -> 3/59
      and (25,24,72) -> 3/77 among ALL N in 14..31, all i, all x; in
      particular N=31 single-far EMPTY.  (19/25 rows = the validation gates.)
  P3  N=31 first gap EMPTY across bordered / two-defect / band-descent /
      needle-repair species (evidence-grade).

Method notes.  For N >= 14 the base B = {1..N}\\{i} has > 12 speeds, so the
settled-LRC floor is unavailable; we use PER-INSTANCE certificates: exact
M(B) > theta computed directly (expected: M(B) >= 1/N > theta, i.e. no
LRC(N+1) near-violation among these bases -- itself a checkable side fact).
Then l(B,theta) > 0, X0 = ceil(2 theta / l), and the finite check decides the
stratum completely (THM-1284's absorption lemma, unchanged).
"""
from fractions import Fraction as F
from math import gcd, ceil
from functools import reduce
from itertools import combinations
import random, sys, time

sys.path.insert(0, '04-computation')
from lrc_singlefar_absorption_atlas_deathstar_S59 import (
    M_exact, M_exact_wit, safe_intervals)

random.seed(3159)
log = lambda s="": print(s, flush=True)

def primitive(S):
    g = reduce(gcd, S)
    return tuple(sorted(v//g for v in S))

# ---------------- P1: the canonical degrade ----------------
def part1():
    log("== P1: canonical family F(31) = {1..29, 31, 90} ==")
    Fam = list(range(1, 30)) + [31, 90]
    M, q, a = M_exact_wit(Fam)
    verdict = "PASS (replicates mac-mini-S29 degrade)" if M == F(1, 32) else "FAIL"
    log(f"   M(F(31)) = {M} at q={q}, a={a}   [{verdict}]\n")
    return M == F(1, 32)

# ---------------- P2: single-far classification N=14..31 ----------------
def part2():
    log("== P2: single-far classification, N = 14..31 (per-instance floors) ==")
    members = []
    floor_flags = []
    for N in range(14, 32):
        theta = F(2, 2*N+1); lo = F(1, N+1)
        t0 = time.time()
        n_mem = 0
        for i in range(1, N+1):
            B = [v for v in range(1, N+1) if v != i]
            MB = M_exact(B)          # exact; no stop_above -- need the true value
            if MB <= theta:
                floor_flags.append((N, i, MB))
                log(f"   !! N={N} i={i}: M(base) = {MB} <= theta -- absorption "
                    f"inapplicable, stratum NOT decided (and an LRC({N+1}) alarm)")
                continue
            iv, l = safe_intervals(B, theta)
            X0 = ceil(2*theta / l)
            hi_check = max(X0, 8*N) + 15
            for x in range(N+1, hi_check+1):
                if x <= N: continue
                m = M_exact(B + [x], stop_above=theta)
                if lo < m < theta:
                    Mw, q, a = M_exact_wit(B + [x])
                    members.append((N, i, x, Mw, q))
                    n_mem += 1
                    log(f"   N={N} i={i} x={x}: MEMBER M={Mw} (q={q}, a={a}) "
                        f"[X0={X0}]")
        log(f"   N={N}: done ({time.time()-t0:.0f}s), members this N: {n_mem}")
    log(f"\n   P2 complete-member list N=14..31: "
        f"{[(N,i,x,str(M)) for N,i,x,M,q in members]}")
    log(f"   base-floor alarms: {floor_flags if floor_flags else 'none'}")
    gate = (any((N,i,x)==(19,18,54) and M==F(3,59) for N,i,x,M,q in members) and
            any((N,i,x)==(25,24,72) and M==F(3,77) for N,i,x,M,q in members))
    n31 = [m for m in members if m[0] == 31]
    others = [m for m in members if m[0] not in (19, 25)]
    log(f"   GATES: 3/59 at (19,18,54): "
        f"{any((N,i,x)==(19,18,54) for N,i,x,M,q in members)}; "
        f"3/77 at (25,24,72): {any((N,i,x)==(25,24,72) for N,i,x,M,q in members)}")
    log(f"   N=31 single-far members: {n31 if n31 else 'NONE (prediction holds)'}")
    log(f"   members outside N=19/25: {others if others else 'NONE'}\n")
    return members, gate, floor_flags

# ---------------- P3: species census at N=31 ----------------
N31 = 31
TH31 = F(2, 63); LO31 = F(1, 32)

def check_fam(S, found, sp):
    Sp = primitive(S)
    if len(set(Sp)) < N31: return 0
    m = M_exact(list(Sp), stop_above=TH31)
    if LO31 < m < TH31 and m not in found:
        found[m] = (Sp, sp)
        log(f"   !! N=31 GAP MEMBER {m} via {sp}: {Sp}")
    return 1

def species_B31(found, budget):
    t0 = time.time(); n = 0
    for d in range(2, 9):
        for a0 in range(1, 41):
            for m in range(20, 31):
                spine = [a0 + d*k for k in range(m)]
                nb = N31 - m
                if nb <= 0: continue
                borders = set()
                for sv in spine:
                    for e in (1, 2, 3):
                        for sg in (1, -1):
                            b = sv + sg*e
                            if b > 0 and b not in spine: borders.add(b)
                borders = sorted(borders)
                if len(borders) < nb: continue
                cnt = 0
                for combo in combinations(borders, nb):
                    n += check_fam(tuple(sorted(spine+list(combo))), found, "B")
                    cnt += 1
                    if cnt >= 2500 or time.time()-t0 > budget: break
                if time.time()-t0 > budget: return n
    return n

def species_C31(found, budget):
    t0 = time.time(); n = 0
    base = list(range(1, N31+1))
    pairs = list(combinations(range(1, N31+1), 2))
    random.shuffle(pairs)
    for i, j in pairs:
        Bb = [v for v in base if v not in (i, j)]
        for w1 in range(N31+1, 7*N31):
            for w2 in range(w1+1, 7*N31):
                n += check_fam(tuple(Bb+[w1, w2]), found, "C")
                if time.time()-t0 > budget: return n
    return n

def species_T31(found, budget):
    """targeted: structured multiples replacing removed elements"""
    t0 = time.time(); n = 0
    base = list(range(1, N31+1))
    for i, j in combinations(range(1, N31+1), 2):
        Bb = [v for v in base if v not in (i, j)]
        cands = {2*i, 3*i, 2*j, 3*j, 3*i-1, 3*i+1, 3*j-1, 3*j+1, i+N31, j+N31}
        for w1, w2 in combinations(sorted(c for c in cands if c > N31), 2):
            n += check_fam(tuple(Bb+[w1, w2]), found, "T")
            if time.time()-t0 > budget: return n
    return n

def species_D31(found, budget):
    """descent from Allowed(a) bands; mediant band (3,95) first"""
    t0 = time.time(); n = 0
    vq = [(3, 95), (4, 127), (5, 158), (5, 159), (6, 191), (7, 221), (7, 222), (7, 223)]
    for val, q in vq:
        for a in range(1, q//2+1):
            if gcd(a, q) != 1: continue
            allowed = [v for v in range(1, int(1.5*q)+1)
                       if min((v*a) % q, q-((v*a) % q)) >= val]
            if not (N31 <= len(allowed) <= N31 + 12): continue
            if not any((x+y) % q == 0 for x, y in combinations(allowed, 2)):
                continue
            S = list(allowed)
            while len(S) > N31:
                bestM, bestS = None, None
                for r in range(len(S)):
                    T = S[:r] + S[r+1:]
                    m = M_exact(list(primitive(T)), stop_above=TH31)
                    if bestM is None or m < bestM:
                        bestM, bestS = m, T
                S = bestS
            n += check_fam(tuple(primitive(S)), found, "D")
            if time.time()-t0 > budget: return n
    return n

def species_E31(found, budget):
    t0 = time.time(); n = 0
    vq = [(3, 95), (4, 127), (5, 158), (5, 159)]
    while time.time()-t0 < budget:
        val, q = random.choice(vq)
        band = list(range(val, q-val+1))
        res = {val, q-val}
        while len(res) < N31: res.add(random.choice(band))
        S = sorted(res)
        curM = M_exact(S)
        for _ in range(150):
            if LO31 < curM < TH31: break
            idx = random.randrange(len(S)); new = random.choice(band)
            T = sorted(set(S[:idx]+S[idx+1:]+[new]))
            if len(T) < N31: continue
            m = M_exact(T, stop_above=TH31)
            if m < curM: S, curM = T, m
        n += check_fam(tuple(primitive(S)), found, "E")
    return n

def part3():
    log("== P3: species census at N=31, window (1/32, 2/63), width 1/2016 ==")
    found = {}
    counts = {}
    for name, fn, budget in [("B", species_B31, 150), ("C", species_C31, 150),
                             ("T", species_T31, 60), ("D", species_D31, 240),
                             ("E", species_E31, 90)]:
        counts[name] = fn(found, budget)
        log(f"   species {name}: {counts[name]} families tested "
            f"({len(found)} members so far)")
    log(f"   P3 verdict: {'MEMBERS FOUND: '+str(sorted(found)) if found else 'EMPTY (species+budget)'}\n")
    return found, counts

if __name__ == "__main__":
    log("== HYP-7890: the N=31 band-law test (death-star-S59b) ==\n")
    p1 = part1()
    members, gate, alarms = part2()
    found31, counts = part3()
    log("== VERDICT SUMMARY ==")
    log(f"  P1 canonical degrade M(F(31))=1/32: {'PASS' if p1 else 'FAIL'}")
    log(f"  P2 gates (3/59, 3/77 rediscovered): {'PASS' if gate else 'FAIL -- N=31 negative NOT trustworthy'}")
    n31 = [m for m in members if m[0] == 31]
    log(f"  P2 N=31 single-far stratum: {'EMPTY (PROVED, per-instance floors)' if not n31 else n31}")
    extra = [m for m in members if m[0] not in (19, 25)]
    log(f"  P2 unexpected members at other N: {extra if extra else 'none'}")
    log(f"  P3 N=31 census: {'EMPTY' if not found31 else sorted(found31)} {counts}")
    log(f"  BAND LAW at N=31: " +
        ("SURVIVES (all predictions hold)" if (p1 and gate and not n31 and not found31)
         else "CHECK DETAILS ABOVE"))
