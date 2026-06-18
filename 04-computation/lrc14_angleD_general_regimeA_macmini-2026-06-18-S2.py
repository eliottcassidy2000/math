#!/usr/bin/env python3
"""
ANGLE D — part 4: does the regime-A arithmetic derivation EXTEND beyond principal towers?

Regime A = the M-binding SUM pair is (flank, w) with w>13 and w privately owning some q.
On principal towers, j>=ceil(D/14) follows NON-tautologically from:
    D = flank + w,  w = lcm(q,14)*m,  j = (14/g)*m   (g=gcd(q,14)),
    14*j - D = 14*(14/g)*m - flank - (14q/g)*m = (196/g - 14q/g)*m - flank
             = (14/g)*(14-q)*m - flank  >= 0  iff  (14/g)*(14-q)*m >= flank.
Since q<=13 and m>=1: (14/g)*(14-q) >= 14/g >= flank-ceiling? We must check the flank.
The principal flanks are small (1,3,5,7,8), and (14-q)*(14/g) dominates -> holds with margin.

GENERAL regime A: w is ANY private-q owner > 13 (a cluster runner), flank ANY smaller
speed. We test whether the SAME inequality structure 14*j >= D holds and WHETHER it is
forced by the private-q residue independent of M.

The honest question: in general regime A is there a UNIFORM arithmetic reason, or does it
still rely on 'others clear' (i.e. is regime A itself only certified a-posteriori)?

We:
 (1) collect all M-binding-pair regime-A records (flank, w, D=flank+w, num, j, priv q of w);
 (2) check 14*j >= D (=M>=1/14) — must hold (0 breaks);
 (3) test the PREDICTIVE form: j = fold(flank*num, D). The crossing index num is forced by
     'others clear'. We test if num is PINNED by a private-q congruence: at tau=num/D the
     private q | w gives w*num/D = (D-flank)*num/D = num - flank*num/D, so
     ||w*num/D|| = ||flank*num/D|| (the pair binds, trivially). The real pin is on the SMALL
     speeds. So the arithmetic of j requires the WHOLE small part, NOT just private-q.
 (4) CONCLUSION metric: fraction of regime-A rows where j (hence M) is determined by the
     private-q owner's lcm structure ALONE (principal-like) vs requiring the full small part.
stdlib only, exact.
"""
from __future__ import annotations
from fractions import Fraction as F
from functools import lru_cache, reduce
from math import gcd
import itertools, random

N = 14; LEVEL = F(1, N); AP13 = tuple(range(1, 14)); RNG = random.Random(99)
def lcm(a, b): return a*b//gcd(a, b)
def gcd_all(v): return reduce(gcd, v, 0)
def norm1(x):
    r = x-int(x)
    if r < 0: r += 1
    return r if r <= F(1, 2) else 1-r
def gv(S, t): return min(norm1(v*t) for v in S)
def fold(r, D):
    r %= D; return min(r, D-r)
@lru_cache(None)
def ctaus(S):
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
@lru_cache(None)
def exact_M(S):
    best = F(0); ats = []
    for t in ctaus(S):
        v = gv(S, t)
        if v > best: best = v; ats = [t]
        elif v == best: ats.append(t)
    return best, tuple(sorted(ats))
def cover_counts(S): return {q: sum(1 for v in S if v % q == 0) for q in range(2, 15)}
def private_owner(S):
    cc = cover_counts(S); out = {}
    for v in S:
        p = tuple(q for q in range(2, 15) if v % q == 0 and cc[q] == 1)
        if p: out[v] = p
    return out
def bsum(S, tau):
    val = gv(S, tau); binders = [v for v in S if norm1(v*tau) == val]; out = []
    for a, b in itertools.combinations(sorted(binders), 2):
        d = a+b
        if norm1(d*tau) == 0:
            num = tau.numerator*(d//tau.denominator)
            out.append({"pair": (a, b), "D": d, "num": num, "j": int(val*d)})
    return out

def enum_s3(vmax=120, target=700):
    rows = set()
    for q in range(2, 14):
        Lq = lcm(q, 14); core = tuple(v for v in AP13 if v != q)
        for m in range(1, 12):
            w = Lq*m
            if w > 2000: break
            S = tuple(sorted(set(core+(w,))))
            if len(S) == 13 and gcd_all(S) == 1 and all(any(v % qq == 0 for v in S) for qq in range(2, 15)):
                rows.add(S)
    for q in range(2, 14):
        core = tuple(v for v in AP13 if v != q); cnt = 0
        for w1 in range(14, vmax+1):
            for w2 in range(w1+1, vmax+1):
                S = tuple(sorted(set(core+(w1, w2))))
                if len(S) != 13 or gcd_all(S) != 1: continue
                if all(any(v % qq == 0 for v in S) for qq in range(2, 15)):
                    rows.add(S); cnt += 1
            if cnt > 250: break
    att = 0
    while len(rows) < target and att < 120000:
        att += 1
        vals = set()
        for q in range(2, 15): vals.add(q*RNG.randint(1, max(1, vmax//q)))
        vals.add(1)
        while len(vals) < 13: vals.add(RNG.randint(1, vmax))
        if len(vals) > 13: vals = set(RNG.sample(sorted(vals), 13))
        S = tuple(sorted(vals))
        if len(S) == 13 and gcd_all(S) == 1 and max(S) > 13 and all(any(v % qq == 0 for v in S) for qq in range(2, 15)):
            rows.add(S)
    return sorted(r for r in rows if max(r) > 13 and any(v <= 13 for v in r))

def main():
    print("="*82)
    print("ANGLE D part 4: general regime-A — is the j>=D/14 derivation arithmetic or a-posteriori?")
    print("="*82)
    rows = enum_s3()
    print(f"S3 rows: {len(rows)}")
    regA = []  # (S, flank, w, D, num, j, privq_of_w)
    breaks = 0
    for S in rows:
        M, taus = exact_M(S)
        if M < LEVEL: breaks += 1; continue
        priv = private_owner(S)
        recs = []
        for tau in taus: recs.extend(bsum(S, tau))
        if not recs: continue
        r = min(recs, key=lambda r: (r["D"], r["pair"]))
        a, b = r["pair"]; big = max(a, b); flank = min(a, b)
        if big > 13 and big in priv:
            regA.append((S, flank, big, r["D"], r["num"], r["j"], priv[big]))
    print(f"regime-A rows (M-binding big>13 privately owns q): {len(regA)}  (LRC breaks={breaks})")

    # (2) check 14*j >= D always
    bad = [t for t in regA if 14*t[5] < t[3]]
    print(f"14*j >= D fails: {len(bad)} (must be 0)")

    # (3) test the PRINCIPAL-LIKE arithmetic: is w = lcm(q,14)*m for the private q, and
    #     does j = (14/g)*m? i.e. is the big a 'clean lcm multiple' giving the slope law?
    principal_like = 0
    nonprincipal = 0
    npe = []
    for (S, flank, w, D, num, j, pq) in regA:
        # the private q with smallest lcm
        ok = False
        for q in pq:
            g = gcd(q, 14); Lq = lcm(q, 14)
            if w % Lq == 0:
                m = w // Lq
                if j == (14 // g) * m and D == flank + w:
                    ok = True; break
        if ok: principal_like += 1
        else:
            nonprincipal += 1
            if len(npe) < 10: npe.append((S, flank, w, D, num, j, pq))
    print(f"\n--- regime-A decomposition ---")
    print(f"  principal-like (w=lcm(q,14)*m AND j=(14/g)*m): {principal_like}/{len(regA)}")
    print(f"  non-principal (j NOT given by the lcm slope law): {nonprincipal}/{len(regA)}")
    print("  non-principal regime-A examples (S, flank, w, D, num, j, privq):")
    for e in npe: print("    ", e)

    # (4) For the principal-like, RESTATE the clean inequality and verify margin > 0.
    print(f"\n--- the clean arithmetic inequality on principal-like regime-A ---")
    print("  14*j - D = (14/g)*(14-q)*m - flank.  margin must be >=0; report min over sample:")
    margins = []
    for (S, flank, w, D, num, j, pq) in regA:
        for q in pq:
            g = gcd(q, 14); Lq = lcm(q, 14)
            if w % Lq == 0:
                m = w // Lq
                if j == (14//g)*m and D == flank + w:
                    margin = (14//g)*(14-q)*m - flank   # = 14*j - D
                    assert margin == 14*j - D
                    margins.append((margin, q, m, flank, S))
                    break
    if margins:
        margins.sort()
        print(f"  min margin (14*j - D) = {margins[0][0]} at q={margins[0][1]} m={margins[0][2]} flank={margins[0][3]}")
        print(f"     row: {margins[0][4]}")
        print(f"  these are all >= 0; the inequality (14/g)(14-q)m >= flank holds because")
        print(f"     (14/g)(14-q) >= (14/g) >= 1 and flank <= 13 < (14/g)(14-q)m in every principal case.")
        # tightest principal cases: q=13 (g=1): (14)(1)*m=14m >= flank=1 ok; q=12(g=2):7*2*m=14m>=5 ok
        tight = [mr for mr in margins if mr[0] < 14]
        print(f"  tightest principal margins (<14): {[(mr[0],mr[1],mr[2],mr[3]) for mr in tight[:10]]}")

    print("\n" + "="*82)
    print("CONCLUSION")
    print("="*82)
    print("  Regime A splits into PRINCIPAL-LIKE (clean lcm slope law, j>=D/14 PROVABLE from")
    print("  private-q arithmetic) and NON-PRINCIPAL (j set by the full small-part alignment,")
    print("  private-q insufficient). The literal HYP-2579 conjecture 'private-q => j>=D/14' is")
    print("  TRUE-but-tautological at the M-crossing, PROVABLE-non-tautologically ONLY on the")
    print("  principal-like sub-family, and FALSE as a per-record arithmetic statement (part 1).")

if __name__ == "__main__":
    main()
