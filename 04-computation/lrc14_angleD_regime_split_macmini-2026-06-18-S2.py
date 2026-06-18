#!/usr/bin/env python3
"""
ANGLE D — part 3: WHEN is the parked/cluster runner w ACTUALLY the binding big?
The regime split that decides whether private-q controls the crossing.

Part 2 found:
  - q in {7,8,9,10,11,12,13}: the M-binding SUM pair is exactly (flank, w); j grows
    linearly slope 14/gcd(q,14); j>=D/14 <=> 14>=q (NON-tautological derivation).
  - q in {2,3,4,5,6}: M is CONSTANT in m (e.g. 2/17, 2/19, 2/23); the binding pair is
    a SMALL-SMALL pair, w is NOT in it. The private-q of w is irrelevant.

This script:
  (1) confirms, per principal row, whether w is a member of the M-binding pair;
  (2) characterizes the split arithmetically (it is whether 7 | q, equiv g=gcd(q,14)=7,
      OR more precisely whether the *small* part {1..13}\{q} already has a small pair
      whose gap >= the (flank,w) gap);
  (3) tests the SHARP claim that decides Angle D:
        "If w (the SOLE multiple of its private q) is a MEMBER of the M-binding pair,
         then j >= ceil(D/14)."
      i.e. private-q controls the crossing ONLY in regime A, and there it WORKS.
  (4) the honest converse: in regime B, LRC holds for a DIFFERENT reason (the small
      part is itself lonely), so private-q is not needed — w is 'over-covered enough'.

We test (3) across the full S3 sample (not just principal): restrict to binding
records where the big member is the SOLE multiple of some q (private), and check
j>=ceil(D/14). If TRUE universally => clean SHARP LEMMA. We already know M=j/D so at
the *M*-crossing it's tautological; the content is whether 'big has private-q AND big
is in the binding pair' can be detected and used. We separate:
   (3a) does the M-binding pair ever have its big = a private-q owner with j<ceil(D/14)?
        (must be 0, else LRC broken)
   (3b) NON-M binding crossings where big is a private-q owner and (flank,big) is the
        UNIQUE smallest-gap pair (i.e. would-be M if others cleared): does j>=ceil(D/14)?
        This is the genuinely arithmetic statement.
stdlib only, exact.
"""
from __future__ import annotations
from fractions import Fraction as F
from functools import lru_cache, reduce
from math import gcd
import itertools, random

N = 14
LEVEL = F(1, N)
AP13 = tuple(range(1, 14))
RNG = random.Random(424242)

def lcm(a, b): return a*b//gcd(a, b)
def gcd_all(v): return reduce(gcd, v, 0)
def norm1(x: F) -> F:
    r = x - int(x)
    if r < 0: r += 1
    return r if r <= F(1, 2) else 1 - r
def g_value(S, t): return min(norm1(v*t) for v in S)
def fold(r, D):
    r %= D; return min(r, D-r)

@lru_cache(maxsize=None)
def candidate_taus(S):
    S = tuple(sorted(set(S))); out = {F(1, 2)}
    for v in S:
        k = 0
        while True:
            t = F(2*k+1, 2*v)
            if t > F(1, 2): break
            out.add(t); k += 1
    for a, b in itertools.combinations(S, 2):
        for d in (a+b, abs(b-a)):
            if d <= 0: continue
            k = 1
            while True:
                t = F(k, d)
                if t > F(1, 2): break
                out.add(t); k += 1
    return frozenset(out)

@lru_cache(maxsize=None)
def exact_M(S):
    best = F(0); ats = []
    for t in candidate_taus(S):
        val = g_value(S, t)
        if val > best: best = val; ats = [t]
        elif val == best: ats.append(t)
    return best, tuple(sorted(ats))

def cover_counts(S): return {q: sum(1 for v in S if v % q == 0) for q in range(2, 15)}
def private_owner(S):
    """v -> tuple of q it solely covers."""
    cc = cover_counts(S); out = {}
    for v in S:
        p = tuple(q for q in range(2, 15) if v % q == 0 and cc[q] == 1)
        if p: out[v] = p
    return out

def binding_sum_records(S, tau):
    val = g_value(S, tau)
    binders = [v for v in S if norm1(v*tau) == val]
    out = []
    for a, b in itertools.combinations(sorted(binders), 2):
        d = a+b
        if norm1(d*tau) == 0:
            num = tau.numerator*(d//tau.denominator)
            out.append({"pair": (a, b), "D": d, "num": num, "j": int(val*d)})
    return out

# ---------- Part A: principal tower membership ----------
def principal_membership():
    print("="*82)
    print("Part A: is w a member of the M-binding pair? (principal single-drop towers)")
    print("="*82)
    print(f"{'q':>2} {'m':>2} {'w':>5} {'M':>11} {'binding pair':>14} {'w in pair?':>10} "
          f"{'regime':>8} {'7|q?':>5}")
    for q in range(2, 14):
        Lq = lcm(q, 14); core = tuple(v for v in AP13 if v != q)
        for m in range(1, 5):
            w = Lq*m
            S = tuple(sorted(set(core+(w,))))
            if len(S) != 13 or gcd_all(S) != 1: continue
            M, taus = exact_M(S)
            recs = []
            for tau in taus: recs.extend(binding_sum_records(S, tau))
            if not recs: continue
            r = min(recs, key=lambda r: (r["D"], r["pair"]))
            pair = r["pair"]; win = w in pair
            regime = "A(w-big)" if win else "B(small)"
            print(f"{q:2d} {m:2d} {w:5d} {str(M):>11} {str(pair):>14} {str(win):>10} "
                  f"{regime:>8} {str(q%7==0):>5}")

# ---------- Part B: the sharp lemma over the full S3 sample ----------
def enumerate_s3(vmax_cap=110, target=500):
    rows = set()
    for q in range(2, 14):
        Lq = lcm(q, 14); core = tuple(v for v in AP13 if v != q)
        for m in range(1, 10):
            w = Lq*m
            if w > 1600: break
            S = tuple(sorted(set(core+(w,))))
            if len(S) == 13 and gcd_all(S) == 1 and all(any(v % qq == 0 for v in S) for qq in range(2, 15)):
                rows.add(S)
    for q in range(2, 14):
        core = tuple(v for v in AP13 if v != q)
        cnt = 0
        for w1 in range(14, vmax_cap+1):
            for w2 in range(w1+1, vmax_cap+1):
                S = tuple(sorted(set(core+(w1, w2))))
                if len(S) != 13 or gcd_all(S) != 1: continue
                if all(any(v % qq == 0 for v in S) for qq in range(2, 15)):
                    rows.add(S); cnt += 1
            if cnt > 200: break
    # seeded random
    att = 0
    while len([r for r in rows if max(r) > 13]) < target and att < 80000:
        att += 1
        vals = set()
        for q in range(2, 15): vals.add(q*RNG.randint(1, max(1, vmax_cap//q)))
        vals.add(1)
        while len(vals) < 13: vals.add(RNG.randint(1, vmax_cap))
        if len(vals) > 13: vals = set(RNG.sample(sorted(vals), 13))
        S = tuple(sorted(vals))
        if len(S) == 13 and gcd_all(S) == 1 and max(S) > 13 and all(any(v % qq == 0 for v in S) for qq in range(2, 15)):
            rows.add(S)
    return sorted(r for r in rows if max(r) > 13 and any(v <= 13 for v in r))

def sharp_lemma_test():
    print("\n" + "="*82)
    print("Part B: SHARP LEMMA — if a private-q owner w is a MEMBER of the M-binding pair,")
    print("        is j >= ceil(D/14)? (this is M>=1/14 read at the crossing; verify 0 breaks)")
    print("="*82)
    rows = enumerate_s3()
    print(f"S3 rows: {len(rows)}")
    member_priv = 0          # M-binding pairs where a member privately owns q
    member_priv_pass = 0
    nonmember_lrc = 0        # rows where NO M-binding member is a private owner (LRC by other means)
    breaks = 0
    examples_nonmember = []
    for S in rows:
        M, taus = exact_M(S)
        if M < LEVEL: breaks += 1; continue
        priv = private_owner(S)
        recs = []
        for tau in taus: recs.extend(binding_sum_records(S, tau))
        if not recs:
            continue
        r = min(recs, key=lambda r: (r["D"], r["pair"]))
        a, b = r["pair"]; D = r["D"]; j = r["j"]; cD = -(-D//14)
        a_priv = priv.get(a, ()); b_priv = priv.get(b, ())
        if a_priv or b_priv:
            member_priv += 1
            if j >= cD: member_priv_pass += 1
        else:
            nonmember_lrc += 1
            if len(examples_nonmember) < 8:
                examples_nonmember.append((S, r["pair"], D, j, M))
    print(f"  M-binding pair has a private-q member: {member_priv}/{len(rows)}")
    print(f"     ...and j>=ceil(D/14): {member_priv_pass}/{member_priv}  (must equal: tautological at M)")
    print(f"  M-binding pair has NO private-q member (LRC via OTHER structure): {nonmember_lrc}/{len(rows)}")
    print(f"  LRC breaks: {breaks}")
    print("  Examples where M-binding pair carries NO private-q (so private-q is NOT the lever):")
    for s in examples_nonmember:
        print("    S=", s[0], " pair=", s[1], " D=", s[2], " j=", s[3], " M=", s[4])

    # ---- Part C: the NON-tautological arithmetic statement ----
    # For each S, consider the (flank, w) SUM pair where w privately owns some q and w>13.
    # Define the 'pair-only M' = best over crossings num of (gap of THIS pair) subject to
    # being the running min — but to make it arithmetic, ask: is there ANY crossing num
    # where (flank,w) is the WIDEST pair gap AND that gap < 1/14 while ALL OTHER speeds
    # also clear (>= that gap)? That would be a counterexample to LRC, so must be empty.
    # We already know empty. The real arithmetic q: does private-q give an a-priori bound?
    # ANSWER from part-2: only when (flank,w) IS the M pair, i.e. regime A. Quantify the
    # split by 'small-part lonely value' L_small = M(P) of the small part alone.
    print("\n" + "="*82)
    print("Part C: the regime split is governed by the SMALL part's own loneliness M(P)")
    print("="*82)
    print("  In regime B, the small part P=S∩{1..13} is ALREADY lonely (M(P)>=1/14) at a")
    print("  small-small crossing that the cluster does not spoil -> private-q irrelevant.")
    print("  In regime A, M(P) is spoiled by the cluster and the (flank,w) crossing decides.")
    regA = regB = 0
    sampleA = []; sampleB = []
    for S in rows:
        M, taus = exact_M(S)
        if M < LEVEL: continue
        P = tuple(v for v in S if v <= 13)
        MP, _ = exact_M(P)
        priv = private_owner(S)
        recs = []
        for tau in taus: recs.extend(binding_sum_records(S, tau))
        if not recs: continue
        r = min(recs, key=lambda r: (r["D"], r["pair"]))
        a, b = r["pair"]
        big = max(a, b)
        w_in_pair_priv = (big > 13 and big in priv)
        if w_in_pair_priv:
            regA += 1
            if len(sampleA) < 4: sampleA.append((S, r["pair"], M, MP))
        else:
            regB += 1
            if len(sampleB) < 4: sampleB.append((S, r["pair"], M, MP))
    print(f"  regime A (M-binding big > 13 and privately owns q): {regA}")
    print(f"  regime B (otherwise; small or non-private big):      {regB}")
    print("  regime A examples (S, pair, M, M(P)):")
    for s in sampleA: print("    ", s)
    print("  regime B examples (S, pair, M, M(P)):")
    for s in sampleB: print("    ", s)

if __name__ == "__main__":
    principal_membership()
    sharp_lemma_test()
