#!/usr/bin/env python3
"""
lrc14_midband_pairsum_anatomy_klein_S209.py

HYP-5731: MID-BAND x PAIR-SUM witness anatomy + certificate stress on the
ACTUAL open family (k>=8 covering sets with mid-band members), complementing
mac-mini HYP-5730 (whose census used covering [1,18] = k<=5 and unstratified
random sets).

Per THM-668 (pair-sum ruler theorem): M(S) is attained at t = p/q with
q = v_i + v_j (i <= j), and LRC(14) for S <=> some pair-sum ruler q has a
LIVE MULTIPLIER p: all v_l * p mod q in the middle band [ceil(q/14), q - ceil(q/14)].

QUESTIONS:
 (Q1) GROUND TRUTH + ANATOMY. For mid-band-heavy instances: where does the
      witness (argmax) live? Does the witness pair ENTANGLE a mid-band member
      (the "self-hosting" mechanism: on ruler q = u + w the mid-band runner u
      merges with w, since v ≡ -u mod q)? Compare against the base rate.
 (Q2) WITNESS SUPPLY. Exact live-multiplier counts LM(q) per ruler and
      #live rulers per instance, stratified by mid-band population s_mid.
      Is the supply ever thin (min over instances of max_q LM)?
 (Q3) CERTIFICATE STRESS. mac-mini's C1 gcd-exact union ledger
      (live if sum_classes(|B|-1) < q-1, classes merged iff v_l ≡ ±v_k mod q):
      does it certify on this stratum? Where it fails, does C4 = HUNTER LEDGER
      (second-order: union <= sum|B| - max-spanning-tree of pairwise overlaps,
      overlaps computed EXACTLY) certify? Report the certificate gap.

All computations exact (integers only).
"""

import random
from math import gcd
from fractions import Fraction

random.seed(20260709)

QS = list(range(2, 15))


def is_covering(S):
    return all(any(s % q == 0 for s in S) for q in QS)


# ---------------------------------------------------------------- generation
def gen_instance(P, V, style, k):
    lo_map = {'confined': max(14, (66 * V) // 100), 'uniform': 14, 'slowheavy': 14}
    lo = lo_map[style]
    L = {V}
    missed = [q for q in QS if not any(p % q == 0 for p in P)]
    for q in missed:
        if any(u % q == 0 for u in L):
            continue
        mlo = (lo + q - 1) // q
        mhi = V // q
        if mlo > mhi:
            return None
        L.add(q * random.randint(mlo, mhi))
    if style == 'slowheavy':
        for _ in range(3):
            if len(L) < k:
                L.add(random.randint(max(14, V // 14 + 1), max(16, (9 * V) // 14 - 1)))
    tries = 0
    while len(L) < k and tries < 500:
        L.add(random.randint(lo, V))
        tries += 1
    if len(L) != k:
        return None
    S = sorted(set(P) | L)
    if len(S) != len(P) + k or not is_covering(S):
        return None
    return S


# ---------------------------------------------------------------- exact ruler analysis
def ruler_analysis(S, q):
    """Exact per-ruler data: kill-union (bytearray), LM(q), C1 ledger bound,
    C4 Hunter bound. Kill window K = {r : r < c or r > q - c}, c = ceil(q/14)
    (residue r is SAFE iff c <= r <= q - c, i.e. dist(r,{0,q}) >= q/14...
    careful: safe means 14*r >= q AND 14*(q-r) >= q, i.e. r in
    [ceil(q/14), floor(13q/14)]; kill = complement)."""
    c = -(-q // 14)  # ceil(q/14)
    # kill residues: r in [0, c-1] or 14*(q-r) < q  <=> r > q - c ... handle exact:
    # safe iff 14*r >= q and 14*(q - r) >= q
    def r_safe(r):
        return 14 * r >= q and 14 * (q - r) >= q

    # classes: merge l,k iff v_l ≡ ±v_k (mod q)
    reps = []          # class representatives (residue up to sign)
    class_g = []       # gcd of the class rep with q
    seen = {}
    for v in S:
        r = v % q
        key = min(r, (q - r) % q)
        if key in seen:
            continue
        seen[key] = True
        reps.append(r)
        class_g.append(gcd(r, q) if r else q)
    nclasses = len(reps)

    # exact kill union over nonzero p
    bad = bytearray(q)
    bad[0] = 1
    Bsets = []
    for r in reps:
        Bl = []
        if r == 0:
            # runner ≡ 0 mod q: kills EVERY p (dead ruler)
            for p in range(q):
                bad[p] = 1
            Bsets.append(None)
            continue
        g = gcd(r, q)
        # v*p mod q takes value m*g; p solutions of r*p ≡ s (mod q): g sols if g|s
        # enumerate kill residues s (multiples of g failing safety), mark their p-fibers
        rr = r // g
        qq = q // g
        inv = pow(rr, -1, qq)
        for m in range(qq):
            s = m * g
            if not r_safe(s):
                p0 = (m * inv) % qq
                for t in range(g):
                    p = p0 + t * qq
                    bad[p] = 1
                    Bl.append(p)
        Bsets.append(set(Bl))
    LM = q - sum(bad)  # count of live p (bad[0]=1 accounts p=0)

    # C1 ledger: sum_classes (|B_class| - 1) < q - 1  => live certified
    c1_ok = None
    if all(b is not None for b in Bsets):
        tot = sum(len(b) for b in Bsets)
        c1_ok = (tot - nclasses) < (q - 1)
        c1_slack = (q - 1) - (tot - nclasses)
    else:
        c1_ok, c1_slack = False, None

    # C4 Hunter: bad_nonzero <= sum|B\{0}| - max spanning tree of |B∩B'| (nonzero p)
    c4_ok = None
    c4_slack = None
    if all(b is not None for b in Bsets) and nclasses >= 2:
        Bs = [b - {0} for b in Bsets]
        n = len(Bs)
        # max spanning tree on overlap weights (Prim)
        w = [[0] * n for _ in range(n)]
        for a in range(n):
            for bidx in range(a + 1, n):
                w[a][bidx] = w[bidx][a] = len(Bs[a] & Bs[bidx])
        in_tree = [False] * n
        in_tree[0] = True
        best = [w[0][i] for i in range(n)]
        tree_w = 0
        for _ in range(n - 1):
            j = max((i for i in range(n) if not in_tree[i]), key=lambda i: best[i])
            tree_w += best[j]
            in_tree[j] = True
            for i in range(n):
                if not in_tree[i] and w[j][i] > best[i]:
                    best[i] = w[j][i]
        hunter_bound = sum(len(b) for b in Bs) - tree_w
        c4_ok = hunter_bound < (q - 1)
        c4_slack = (q - 1) - hunter_bound
    return LM, c1_ok, c1_slack, c4_ok, c4_slack


def full_anatomy(S):
    """Exact M(S) + argmax over all pair-sum points, all live rulers, and
    certificate results per ruler."""
    S = sorted(S)
    V = S[-1]
    n = len(S)
    rulers = {}
    for i in range(n):
        for j in range(i, n):
            q = S[i] + S[j]
            rulers.setdefault(q, []).append((S[i], S[j]))
    bestM = Fraction(0)
    bestat = None
    per_ruler = {}
    for q, pairs in sorted(rulers.items()):
        # exact max over p of min_l dist(v_l p mod q, {0,q})/q, via integer scan
        best_num = -1
        best_p = None
        for p in range(1, q):
            mn = q  # min over runners of dist*14-ish; store raw dist numerator
            for v in S:
                r = (v * p) % q
                d = r if r <= q - r else q - r
                if d < mn:
                    mn = d
                    if 14 * mn < q and mn < best_num:
                        break
            if mn > best_num:
                best_num = mn
                best_p = p
        Mq = Fraction(best_num, q)
        if Mq > bestM:
            bestM = Mq
            bestat = (q, best_p)
        LM, c1, c1s, c4, c4s = ruler_analysis(S, q)
        per_ruler[q] = dict(pairs=pairs, M=Mq, LM=LM, c1=c1, c4=c4)
    return bestM, bestat, per_ruler


def midband(S):
    V = max(S)
    return [v for v in S if 14 * v > V and 14 * v < 9 * V]


def main():
    print("=" * 78)
    print("HYP-5731: mid-band x pair-sum anatomy + certificate stress (klein-S209)")
    print("=" * 78)

    Ps = [(8, 9, 10, 12), (7, 9, 10, 11, 12), (11, 12, 13), (10, 11, 12, 13), (9, 11, 13)]
    instances = []
    for P in Ps:
        for style in ('slowheavy', 'uniform', 'confined'):
            got, att = 0, 0
            while got < 5 and att < 80:
                att += 1
                V = random.choice([100, 140, 180, 240, 300])
                k = 13 - len(P)
                if k < 8:
                    continue
                S = gen_instance(list(P), V, style, k)
                if S is None:
                    continue
                instances.append(S)
                got += 1
    # the 7-structured @91 cluster as a speed set (mac-mini/kps hard instance)
    hard91 = sorted(91 - e for e in [0, 7, 14, 21, 26, 29, 37, 44, 51, 58, 67, 75, 82])
    if len(hard91) == 13 and is_covering(hard91):
        instances.append(hard91)
        print(f"  [included 7-structured @91 speed set: {hard91}, covering]")
    else:
        print(f"  [7-structured @91 speed set covering={is_covering(hard91)}; included anyway for anatomy]")
        instances.append(hard91)

    print(f"\ninstances: {len(instances)}")
    strata = {}
    cert_gap_instances = []
    entangle_num = entangle_den = 0
    base_rates = []
    supply_by_smid = {}
    c1_fail_rulers_live = 0
    c1_fail_c4_ok = 0
    total_rulers = 0

    for S in instances:
        V = max(S)
        mb = midband(S)
        smid = len(mb)
        M, at, per = full_anatomy(S)
        q_w, p_w = at
        pairs_w = per[q_w]['pairs']
        # entanglement: any pair representation of the witness q containing a mid-band member
        ent = any((a in mb) or (b in mb) for (a, b) in pairs_w)
        if smid > 0:
            entangle_den += 1
            entangle_num += 1 if ent else 0
            # base rate: fraction of the 91 pairs with a mid-band member
            n_mb = smid
            base = 1 - (13 - n_mb) * (12 - n_mb) / (13 * 12)
            base_rates.append(base)
        live_rulers = [q for q, d in per.items() if d['LM'] > 0]
        maxLM = max(d['LM'] for d in per.values())
        c1_live = [q for q, d in per.items() if d['c1']]
        c4_live = [q for q, d in per.items() if d['c4']]
        for q, d in per.items():
            total_rulers += 1
            if d['LM'] > 0 and not d['c1']:
                c1_fail_rulers_live += 1
                if d['c4']:
                    c1_fail_c4_ok += 1
        st = strata.setdefault(min(smid, 5), dict(n=0, Mmin=None, ent=0, enttot=0,
                                                  supply=[], c1cert=0, c4cert=0, anycert=0))
        st['n'] += 1
        st['Mmin'] = M if st['Mmin'] is None or M < st['Mmin'] else st['Mmin']
        if smid > 0:
            st['enttot'] += 1
            st['ent'] += 1 if ent else 0
        st['supply'].append((len(live_rulers), maxLM))
        st['c1cert'] += 1 if c1_live else 0
        st['c4cert'] += 1 if (c1_live or c4_live) else 0
        if not c1_live and not c4_live:
            cert_gap_instances.append((S, M, at, maxLM, len(live_rulers)))
        supply_by_smid.setdefault(min(smid, 5), []).append(maxLM)

    print("\n[Q1/Q2] strata by mid-band population s_mid (5 = '>=5'):")
    print(f"{'s_mid':>6} {'n':>4} {'min M':>10} {'M>=1/14':>8} {'witness-entangles-mb':>22} "
          f"{'base':>6} {'avg #live rulers':>17} {'min maxLM':>10}")
    for smid in sorted(strata):
        st = strata[smid]
        nlive = [a for a, b in st['supply']]
        maxlms = [b for a, b in st['supply']]
        entfrac = (st['ent'] / st['enttot']) if st['enttot'] else float('nan')
        print(f"{smid:>6} {st['n']:>4} {str(st['Mmin']):>10} "
              f"{'YES' if st['Mmin'] >= Fraction(1,14) else 'NO':>8} "
              f"{entfrac if st['enttot'] else '-':>22.2} "
              f"{(sum(base_rates)/len(base_rates) if base_rates else 0):>6.2f} "
              f"{sum(nlive)/len(nlive):>17.1f} {min(maxlms):>10}")

    print(f"\n  witness-pair entanglement overall: {entangle_num}/{entangle_den} "
          f"= {entangle_num/max(1,entangle_den):.2f} vs mean base rate "
          f"{sum(base_rates)/max(1,len(base_rates)):.2f}")

    print(f"\n[Q3] certificate stress over all {total_rulers} rulers:")
    print(f"  rulers exactly-live but NOT C1-certified: {c1_fail_rulers_live}")
    print(f"     ... of which C4 (Hunter ledger) certifies: {c1_fail_c4_ok}")
    per_inst_c1 = sum(st['c1cert'] for st in strata.values())
    per_inst_any = sum(st['c4cert'] for st in strata.values())
    print(f"  instances with >=1 C1-certified ruler: {per_inst_c1}/{len(instances)}")
    print(f"  instances with >=1 C1-or-C4-certified ruler: {per_inst_any}/{len(instances)}")
    print(f"  instances with NO certificate (but exact truth lonely): {len(cert_gap_instances)}")
    for (S, M, at, maxLM, nl) in cert_gap_instances[:6]:
        print(f"    GAP: S={S} M={M} witness={at} maxLM={maxLM} #live={nl}")

    # anatomy detail on the @91 hard set + 3 slowheavy examples
    print("\n[detail] witness anatomy on named/hard instances:")
    for S in [hard91] + [x for x in instances if len(midband(x)) >= 3][:3]:
        V = max(S)
        mb = midband(S)
        M, at, per = full_anatomy(S)
        q_w, p_w = at
        print(f"  S={S}")
        print(f"    V={V} mid-band={mb} (s={len(mb)})")
        print(f"    M={M}={float(M):.5f} witness q={q_w} (q/V={q_w/V:.2f}) p={p_w} "
              f"pairs={per[q_w]['pairs']} r(q)={len(per[q_w]['pairs'])}")
        top = sorted(per.items(), key=lambda kv: -kv[1]['LM'])[:3]
        print(f"    top-supply rulers: " + ", ".join(
            f"q={q} LM={d['LM']} c1={'Y' if d['c1'] else 'n'}{'/c4Y' if d['c4'] else ''}"
            for q, d in top))

    print("\nDONE.")


if __name__ == '__main__':
    main()
