#!/usr/bin/env python3
"""
ANGLE D — the arithmetic dual (codex HYP-2579).

Question: does the parked runner's PRIVATE q-obligation force the THM-524 binding
crossing index j >= D/14 (=> M(S) >= 1/14)?

THM-524: M(S) = max_tau min_v ||v tau|| is attained at a binding pair crossing
   tau* = num/D  with D = v_a +/- v_b  (a binding pair), and
   M(S) = j/D  where  j = (flank * num) mod D  in the sum case
   (more precisely: ||flank * tau*|| = M, and at tau*=num/D this equals
    dist(flank*num mod D, 0)/D = min(r, D-r)/D where r = flank*num mod D).

The covering S3 residual: cluster contains a runner w that is the SOLE multiple
of some q in 2..14 ("private q-obligation").

KEY HONESTY POINT (codex's own caveat, classifier line ~455):
   "j>=D/14 is tautological unless j is derived from non-M arithmetic."
   Because M = j/D, proving "j >= ceil(D/14)" IS proving M >= 1/14. So the whole
   content is whether the private-q structure ARITHMETICALLY predicts / lower-bounds j
   WITHOUT first computing M. We test exactly that.

We:
  1. Enumerate many primitive q-covering case-S3 13-sets.
  2. Compute exact M, the binding crossing (pair, kind, D, num, j), and the flank.
  3. Identify private-q obligations of every parked runner (w%14==0) AND of every
     cluster runner (the THM-527 cluster L = speeds > 13), and the SOLE-multiple-of-q
     structure generally.
  4. Test the literal conjecture j >= ceil(D/14) on all binding crossings.
  5. THE ARITHMETIC DUAL TEST: try to predict j (or a lower bound on j) from the
     private-q residues of the flank/D alone, independent of M. Specifically:
       - flank, D, num determine j = flank*num mod D (folded). num is coprime-ish to D.
       - if a private q | w forces a residue relation on D = flank +/- w, can we bound
         (flank*num mod D)?
  6. Correlate j/D with private-q structure; look for the clean lemma or obstruction.

stdlib only; EXACT Fractions / ints.
"""
from __future__ import annotations
from fractions import Fraction as F
from functools import lru_cache, reduce
from math import gcd
import itertools, random

N = 14
LEVEL = F(1, N)
QS = tuple(range(2, N + 1))
AP13 = tuple(range(1, 14))
RNG = random.Random(20260618)


def lcm(a, b): return a * b // gcd(a, b)
def gcd_all(vals): return reduce(gcd, vals, 0)

def norm1(x: F) -> F:
    r = x - int(x)
    if r < 0: r += 1
    return r if r <= F(1, 2) else 1 - r

def g_value(S, tau): return min(norm1(v * tau) for v in S)

@lru_cache(maxsize=None)
def candidate_taus(S):
    S = tuple(sorted(set(S)))
    out = {F(1, 2)}
    for v in S:
        k = 0
        while True:
            t = F(2 * k + 1, 2 * v)
            if t > F(1, 2): break
            out.add(t); k += 1
    for a, b in itertools.combinations(S, 2):
        for d in (a + b, abs(b - a)):
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
        if val > best:
            best = val; ats = [t]
        elif val == best:
            ats.append(t)
    return best, tuple(sorted(ats))

def is_q_covering(S):
    return all(any(v % q == 0 for v in S) for q in QS)

def cover_counts(S):
    return {q: sum(1 for v in S if v % q == 0) for q in QS}

def private_q_of(S, w):
    """q in 2..14 such that q | w and w is the SOLE multiple of q in S."""
    cc = cover_counts(S)
    return tuple(q for q in QS if w % q == 0 and cc[q] == 1)

def all_private(S):
    """map speed -> tuple of q it privately owns (sole multiple)."""
    cc = cover_counts(S)
    out = {}
    for v in S:
        priv = tuple(q for q in QS if v % q == 0 and cc[q] == 1)
        if priv:
            out[v] = priv
    return out

def case_of(S):
    """S3 = mixed: has small part P=S∩{1..13} nonempty AND cluster L={v>13} with |L|>=1,
       and is 'genuinely mixed' covering. We classify by k=|cluster| and spread."""
    P = tuple(v for v in S if v <= 13)
    L = tuple(v for v in S if v > 13)
    return P, L

def binding_at(S, tau):
    """Return all binding-pair records at tau where M is attained.
       Each record: pair (a,b), kind sum/diff, D, num (tau=num/D reduced), j (=M*D folded)."""
    val = g_value(S, tau)
    binders = [v for v in S if norm1(v * tau) == val]
    recs = []
    for a, b in itertools.combinations(sorted(binders), 2):
        for kind, d in (("sum", a + b), ("diff", abs(b - a))):
            if d <= 0: continue
            if norm1(d * tau) == 0:
                # tau = num/d reduced
                num = tau.numerator
                den = tau.denominator
                # tau = num/den; d*tau integer means den | d
                # express tau as p/d: p = num * d // den
                assert d % den == 0
                p = num * (d // den)
                # j = (val * d) but we also want the *crossing index* p (num) and the folded residue
                jD_frac = val * d  # = M*D, mathematically an integer (folded gap index)
                assert jD_frac.denominator == 1
                jD = int(jD_frac)
                recs.append({
                    "pair": (a, b), "kind": kind, "D": d, "num": p,
                    "tau": tau, "j": jD, "M_times_D": jD,
                })
    return recs

def flank_of_pair(a, b, cluster):
    """The 'flank' = the small/non-cluster member of the binding pair (the one that is
       NOT the big parked/cluster runner). Returns (flank, big)."""
    # big = the larger; flank = smaller, per THM-524 'small binding flank' language.
    big = max(a, b); flank = min(a, b)
    return flank, big


def fold(r, D):
    """folded residue: min(r mod D, D - r mod D)."""
    r = r % D
    return min(r, D - r)


def enumerate_covering_s3(max_extra=2, vmax_cap=120, target=400):
    """Enumerate primitive q-covering 13-sets that are case S3 (mixed: small part +
       a large cluster k>=1 with Vmax>13). We use:
        - principal single-drop towers {1..13}\{q} ∪ {lcm(q,14)*k}
        - drop-one + cluster (two big runners)
        - drop-two cores + two/three big runners
        - seeded random covering rows
    Returns set of tuples."""
    rows = set()
    # principal single-drop
    for q in range(2, 14):
        Lc = lcm(q, 14)
        core = tuple(v for v in AP13 if v != q)
        for k in range(1, 12):
            w = Lc * k
            if w > vmax_cap * 6: break
            S = tuple(sorted(set(core + (w,))))
            if len(S) == 13 and gcd_all(S) == 1 and is_q_covering(S):
                rows.add(S)
    # drop-one + 2 big runners (mixed clusters k>=2)
    for q in range(2, 14):
        core = tuple(v for v in AP13 if v != q)  # 12 small
        for w1 in range(14, vmax_cap + 1):
            for w2 in range(w1 + 1, vmax_cap + 1):
                S = tuple(sorted(set(core + (w1, w2))))
                if len(S) != 13: continue
                if gcd_all(S) != 1: continue
                if is_q_covering(S):
                    rows.add(S)
            if len(rows) > target * 3: break
        if len(rows) > target * 3: break
    # drop-two cores + cluster
    drops = list(itertools.combinations(range(2, 14), 2))
    RNG.shuffle(drops)
    for (q1, q2) in drops[:40]:
        core = tuple(v for v in AP13 if v not in (q1, q2))  # 11 small
        for _ in range(300):
            ws = sorted(RNG.sample(range(14, vmax_cap + 1), 2))
            S = tuple(sorted(set(core + tuple(ws))))
            if len(S) != 13: continue
            if gcd_all(S) != 1: continue
            if is_q_covering(S):
                rows.add(S)
        if len(rows) > target * 3: break
    # seeded random covering
    attempts = 0
    while len([r for r in rows if max(r) > 13]) < target and attempts < 60000:
        attempts += 1
        vals = set()
        for q in QS:
            vals.add(q * RNG.randint(1, max(1, vmax_cap // q)))
        vals.add(1)
        while len(vals) < 13:
            vals.add(RNG.randint(1, vmax_cap))
        if len(vals) > 13:
            vals = set(RNG.sample(sorted(vals), 13))
        S = tuple(sorted(vals))
        if len(S) == 13 and gcd_all(S) == 1 and is_q_covering(S) and max(S) > 13:
            rows.add(S)
    # keep only S3 (mixed: small part nonempty AND cluster nonempty), and Vmax>13
    out = []
    for S in rows:
        P, L = case_of(S)
        if len(P) >= 1 and len(L) >= 1 and max(S) > 13:
            out.append(S)
    return sorted(out, key=lambda S: (max(S), S))


def main():
    print("=" * 80)
    print("ANGLE D — arithmetic dual of HYP-2579: does private-q force j >= ceil(D/14)?")
    print("=" * 80)
    rows = enumerate_covering_s3(target=500)
    print(f"enumerated primitive q-covering case-S3 13-sets: {len(rows)}")

    # --- Part 1: the literal conjecture on the M-binding crossing ---
    # For each S, at every argmax tau, list binding records; the M crossing has
    # M = j/D for the binding pair. Conjecture: there is a binding pair whose big
    # member has a private-q obligation and j >= ceil(D/14).
    n_break = 0
    n_total = 0
    literal_pass = 0   # rows where SOME sum-binding pair (flank,big) has big with private-q AND j>=ceil(D/14)
    literal_have_private_pair = 0
    rows_priv_big_in_binding = 0
    examples = []
    j_over_D_all = []
    big_private_records = []  # (M, D, j, flank, big, private_qs_of_big, jD>=ceil(D/14))
    for S in rows:
        n_total += 1
        M, taus = exact_M(S)
        if M < LEVEL:
            n_break += 1
            examples.append(("BREAK", S, M))
            continue
        priv = all_private(S)
        # collect all binding records over all argmax taus
        recs = []
        for tau in taus:
            recs.extend(binding_at(S, tau))
        # for each binding record compute flank/big and j/D
        any_priv_pair = False
        any_priv_pass = False
        for r in recs:
            a, b = r["pair"]
            D = r["D"]; jD = r["j"]
            flank, big = flank_of_pair(a, b, None)
            jr = jD  # M*D
            jod = F(jr, D)
            j_over_D_all.append(jod)
            ceilD14 = -(-D // 14)  # ceil(D/14)
            big_priv = priv.get(big, ())
            flank_priv = priv.get(flank, ())
            if big_priv:
                any_priv_pair = True
                big_private_records.append((M, D, jr, flank, big, big_priv, jr >= ceilD14, r["kind"]))
                if jr >= ceilD14 and r["kind"] == "sum":
                    any_priv_pass = True
        if any_priv_pair:
            rows_priv_big_in_binding += 1
        if any_priv_pass:
            literal_pass += 1

    print(f"\n--- Part 1: literal conjecture status ---")
    print(f"LRC breaks (M<1/14): {n_break}/{n_total}")
    print(f"rows where a binding pair's BIG member has a private-q: {rows_priv_big_in_binding}/{n_total}")
    print(f"rows with a sum-binding (flank,big-with-private-q) having j>=ceil(D/14): {literal_pass}/{n_total}")
    # Since M = j/D always, j/D >= 1/14 iff M-binding. So check that:
    nfold_ge = sum(1 for jod in j_over_D_all if jod >= LEVEL)
    print(f"binding records with j/D >= 1/14: {nfold_ge}/{len(j_over_D_all)}  (these are the M-realizing ones)")

    # --- Part 2: THE TAUTOLOGY TEST ---
    # M = j/D exactly at the binding crossing. So "j>=ceil(D/14)" <=> "M>=1/14".
    # The arithmetic content can ONLY be: predict j (or a lower bound) from
    # private-q residues of (flank, D, num) WITHOUT knowing M.
    # Here j_realized = fold(flank*num, D)  [= M*D].  We verify this identity, then
    # test whether private-q gives an independent lower bound on fold(flank*num,D).
    print(f"\n--- Part 2: verify j = fold(flank*num, D) at the M crossing (the identity) ---")
    id_ok = 0; id_bad = 0
    sample_idrows = []
    for S in rows:
        M, taus = exact_M(S)
        if M < LEVEL: continue
        for tau in taus:
            for r in binding_at(S, tau):
                a, b = r["pair"]; D = r["D"]; num = r["num"]; jD = r["j"]
                flank, big = flank_of_pair(a, b, None)
                # at tau=num/D, ||flank*tau|| = fold(flank*num, D)/D
                fr = fold(flank * num, D)
                if fr == jD:
                    id_ok += 1
                else:
                    # the folded residue of the OTHER member equals jD too (pair binds)
                    fb = fold(big * num, D)
                    if fb == jD:
                        id_ok += 1
                    else:
                        id_bad += 1
                        if len(sample_idrows) < 5:
                            sample_idrows.append((S, r, flank, big, fr, fb, jD))
    print(f"  identity j == fold(member*num, D) holds: {id_ok} ok, {id_bad} bad")
    for s in sample_idrows:
        print("   MISMATCH:", s)

    # --- Part 3: arithmetic dual — can private-q lower-bound fold(flank*num,D)? ---
    # The private q | big means big ≡ 0 mod q and NO other speed ≡ 0 mod q.
    # D = flank ± big. For the SUM case D = flank+big. num/D = tau* with gcd(num,D) small.
    # The crossing index num is the THM-524 'k' s.t. tau*=k/D. We probe what private-q
    # tells us about num and hence j=fold(flank*num,D).
    print(f"\n--- Part 3: arithmetic-dual probe ---")
    # Tabulate, for sum-binding records where big has private-q: relation among
    # (private q of big), D mod q, num mod q, j, ceil(D/14).
    print("  sum-binding records (big has private-q), sorted by j/D ascending — hardest:")
    sumrecs = [t for t in big_private_records if t[7] == "sum"]
    sumrecs.sort(key=lambda t: F(t[2], t[1]))
    seen = 0
    for (M, D, jr, flank, big, bp, ge, kind) in sumrecs:
        if seen >= 18: break
        seen += 1
        # private q vs D
        qrel = [(q, D % q, big % q, flank % q) for q in bp]
        print(f"   M={str(M):>9} jod={str(F(jr,D)):>8} D={D:4d} j={jr:3d} ceil(D/14)={-(-D//14):3d} "
              f"j>=ceilD14={ge} flank={flank:3d} big={big:4d} privq(big)={bp} "
              f"[D%q,big%q,flank%q]={qrel}")

    # --- Part 4: does EVERY hard residual have the binding-big carrying private-q? ---
    # The crux of HYP-2579: arc-width residuals have a parked runner with private q.
    # Test: among rows where M is small (slack small), is the binding-pair big member
    # ALWAYS privately-owning some q? And is the flank in {2,4,5,13} (THM-524 part C)?
    print(f"\n--- Part 4: the binding flank and private-big on the hardest rows ---")
    hard = []
    for S in rows:
        M, taus = exact_M(S)
        if M < LEVEL: continue
        priv = all_private(S)
        # the canonical M crossing: take the sum-binding with smallest D (THM-524 part C language)
        recs = []
        for tau in taus:
            recs.extend([r for r in binding_at(S, tau) if r["kind"] == "sum"])
        if not recs: continue
        r = min(recs, key=lambda r: (r["D"], r["pair"]))
        a, b = r["pair"]; D = r["D"]
        flank, big = flank_of_pair(a, b, None)
        hard.append((M, S, flank, big, D, r["j"], priv.get(big, ()), priv.get(flank, ())))
    hard.sort(key=lambda t: (t[0], t[1]))
    print(f"  {'M':>10} {'flank':>5} {'big':>5} {'D':>5} {'j':>4} {'priv(big)':>14} {'priv(flank)':>12}")
    flank_palette = set()
    big_has_priv = 0
    for (M, S, flank, big, D, j, bp, fp) in hard[:25]:
        flank_palette.add(flank)
        if bp: big_has_priv += 1
        print(f"  {str(M):>10} {flank:5d} {big:5d} {D:5d} {j:4d} {str(bp):>14} {str(fp):>12}  S={S}")
    print(f"  flank palette (hardest 25): {sorted(flank_palette)}")
    print(f"  big member carries a private-q (hardest 25): {big_has_priv}/25")

    # --- Part 5: counterexample hunt for the literal claim across ALL binding records ---
    # The literal claim 'private q-debt on big => j>=ceil(D/14)' as a PER-RECORD statement
    # (not just at the M-optimum). Find a sum-binding record (anywhere, not nec. M-optimal)
    # where big has private-q but fold(flank*num,D) < ceil(D/14). This is the real test of
    # the *arithmetic* (non-tautological) claim — because non-M-optimal crossings also exist.
    print(f"\n--- Part 5: per-record arithmetic test (NON-tautological) ---")
    print("  Over ALL sum-crossings tau=num/D (num=1..floor(D/2), gcd(num,D)=1) with big|q private,")
    print("  is fold(flank*num,D) >= ceil(D/14) ALWAYS? (this is the arithmetic statement, not M)")
    per_record_fail = []
    checked = 0
    for S in rows:
        priv = all_private(S)
        Sset = set(S)
        for big, bp in priv.items():
            if big <= 13:  # we care about parked/cluster big runners
                continue
            for flank in S:
                if flank >= big: continue
                D = flank + big
                ceilD14 = -(-D // 14)
                # all crossings tau=num/D, num=1..D//2
                for num in range(1, D // 2 + 1):
                    if gcd(num, D) != 1:
                        continue
                    checked += 1
                    fr = fold(flank * num, D)
                    if fr < ceilD14:
                        # this is a crossing where the (flank,big) sum-pair gap < 1/14
                        # does it correspond to a true M? Not necessarily — others may clear.
                        per_record_fail.append((S, flank, big, D, num, fr, ceilD14, bp))
    print(f"  total sum-crossings checked (flank<big, big privately owns q): {checked}")
    print(f"  crossings where fold(flank*num,D) < ceil(D/14): {len(per_record_fail)}")
    print("  => the PER-RECORD arithmetic claim is FALSE if this is >0 (expected): the pair gap")
    print("     can be tiny; LRC is saved only because OTHER runners are NOT clear there.")
    for s in per_record_fail[:6]:
        print("    e.g.", s)

    # --- Part 6: the RIGHT arithmetic question ---
    # At a sum-crossing tau=num/D with D=flank+big, ALL of S must be 'binding-or-worse'
    # for it to be the M optimum, i.e. min_v ||v num/D|| = fold(flank*num,D).
    # The private-q of big controls big mod stuff; but j is determined by flank*num mod D.
    # Test the CONTRAPOSITIVE structure: at the TRUE M crossing, is num always ≡ ±1-ish
    # giving j = fold(flank*num,D) and does private-q pin num? Report num distribution.
    print(f"\n--- Part 6: at the TRUE M-crossing, distribution of num and num mod (private q) ---")
    from collections import Counter
    num_mod_q = Counter()
    num_vals = Counter()
    for S in rows:
        M, taus = exact_M(S)
        if M < LEVEL: continue
        priv = all_private(S)
        recs = []
        for tau in taus:
            recs.extend([r for r in binding_at(S, tau) if r["kind"] == "sum"])
        if not recs: continue
        r = min(recs, key=lambda r: (r["D"], r["pair"]))
        a, b = r["pair"]; D = r["D"]; num = r["num"]
        flank, big = flank_of_pair(a, b, None)
        num_vals[num] += 1
        for q in priv.get(big, ()):
            num_mod_q[(q, num % q)] += 1
    print(f"  num value histogram at M-crossing (top): {dict(num_vals.most_common(12))}")
    print(f"  (q, num mod q) histogram for private-q of big (top): {dict(num_mod_q.most_common(15))}")

    # --- conclusion lines ---
    print("\n" + "=" * 80)
    print("VERDICT")
    print("=" * 80)
    print(f"LRC breaks: {n_break} (0 expected — LRC(14) holds on all sampled S3 rows)")
    print(f"M-crossing binding-big carries private-q on hardest 25: {big_has_priv}/25")
    print("Literal 'j>=ceil(D/14)' at M-crossing is TAUTOLOGICAL (M=j/D).")
    print(f"Per-record arithmetic claim FALSE: {len(per_record_fail)} sum-crossings have pair-gap<1/14.")


if __name__ == "__main__":
    main()
