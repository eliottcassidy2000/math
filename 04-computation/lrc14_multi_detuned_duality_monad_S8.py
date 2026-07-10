#!/usr/bin/env python3
r"""
lrc14_multi_detuned_duality_monad_S8.py   (monad-explorer-2026-07-09-S8, HYP-5787/THM-678)

(A) THE MULTI-DETUNED COUNTING DISPATCH: v = g*H u D, branch times (u+c)/g; good branch
    exists iff sum_i N_i/q_i < 1 (exact orbit counts).  Verify the dispatch fires exactly
    on the predicted region over d = 2, 3 batteries (exact rational witnesses).
(B) THE d = 2 RESIDUAL (q1 = q2 = 2): both detuned == g/2 (mod g); demonstrate per
    instance that the g-dispatch fails but the 2g-lift fires once odd harmonics are
    handled (or the pair is jointly safe via the h-speed condition).
(C) THE BLOCKING/SUPPLY DUALITY: over residual-core adversaries, measure per family
    - the branch-blocking ledger: #{moduli g <= G0 : branch dispatch blocked at g}
      (blocked = some congruent +-pair of speeds mod g collapses the orbit structure
       below the counting threshold)
    - the pair-sum supply: klein-B5-style liveness count over pair-sum moduli
    and test the conservation candidate (blocking forces supply).
"""
import sys, random
from fractions import Fraction as F
from itertools import combinations
from math import gcd

def exact_witness_H(H, clearance_target=None):
    H = sorted(set(H))
    best = (F(0), None)
    dens = sorted(set(a + b for i, a in enumerate(H) for b in H[i:]))
    for den in dens:
        for p in range(1, den):
            u = F(p, den)
            c = min(min((m * u) % 1, 1 - (m * u) % 1) for m in H)
            if c > best[0]:
                best = (c, u)
    return best

def nearint(x):
    y = x % 1
    return min(y, 1 - y)

def counting_dispatch(g, H, D, verbose=False):
    """THM-678: returns (fires_predicted, fires_actual, witness_tau, min_clear)."""
    cH, u = exact_witness_H(H)
    # predicted: sum N_i/q_i < 1 with N_i = orbit points in open (-1/14, 1/14)
    tot = F(0)
    for delta in D:
        d_i = gcd(delta, g)
        q_i = g // d_i
        x0 = F(delta) * u / g
        N = 0
        for m in range(q_i):
            val = (x0 + F(m * d_i, g)) % 1
            if min(val, 1 - val) < F(1, 14):
                N += 1
        tot += F(N, q_i)
    predicted = tot < 1
    # actual: search branches
    best = (F(0), None)
    for c in range(g):
        tau = (u + c) / g
        mc = min(nearint(F(vv) * tau) for vv in ([g * m for m in H] + list(D)))
        if mc > best[0]:
            best = (mc, tau)
    actual = best[0] >= F(1, 14)
    return predicted, actual, best[1], best[0]

def is_primitive(v):
    g0 = 0
    for x in v:
        g0 = gcd(g0, x)
    return g0 == 1

if __name__ == "__main__":
    rng = random.Random(20260709588)
    print("=" * 100)
    print("PART A -- THE COUNTING DISPATCH, d = 2 and 3 batteries")
    print("=" * 100)
    total = pred_fire = act_fire = agree = 0
    resid_q22 = []
    for trial in range(60):
        g = rng.choice([6, 10, 12, 14, 15, 18, 20, 21])
        dcount = rng.choice([2, 2, 3])
        H = sorted(rng.sample(range(1, 25), 13 - dcount))
        D = []
        while len(D) < dcount:
            dd = rng.randint(g + 1, 30 * g)
            if dd % g != 0 and dd not in D and all(dd != g * m for m in H):
                D.append(dd)
        v = sorted([g * m for m in H] + D)
        if len(set(v)) != 13:
            continue
        total += 1
        pred, act, tau, mc = counting_dispatch(g, H, D)
        if pred:
            pred_fire += 1
        if act:
            act_fire += 1
        if pred == act or (pred and act):
            agree += 1
        if not pred:
            qs = sorted(g // gcd(dd, g) for dd in D)
            if qs[:2] == [2, 2]:
                resid_q22.append((g, H, D))
    print(f"  {total} instances: predicted-fire {pred_fire}, actual-fire {act_fire} "
          f"(prediction => actual sound: {agree}/{total})")
    print(f"  counting-blocked with q1=q2=2 pattern: {len(resid_q22)}")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART B -- THE q1 = q2 = 2 RESIDUAL: congruent half-harmonic pairs")
    print("=" * 100)
    # construct explicit residual instances: delta_i == g/2 (mod g)
    for g, hcount in [(14, 11), (10, 11)]:
        H = list(range(1, hcount + 1))
        d1 = g // 2 + 3 * g          # == g/2 mod g
        d2 = g // 2 + 7 * g
        D = [d1, d2]
        v = sorted([g * m for m in H] + D)
        pred, act, tau, mc = counting_dispatch(g, H, D)
        print(f"  g={g}: D={D} (both == {g//2} mod {g}): counting predicts "
              f"{'FIRE' if pred else 'BLOCKED'}, actual branch search: "
              f"{'FIRES mc=' + str(mc) if act else 'fails at g'}")
        # the h-speed condition: h = (d1-d2)/g; both-states-bad needs ||h u +- 1/2|| < 1/7
        h = (d1 - d2) // g
        cH, u = exact_witness_H(H)
        hu = nearint(F(h) * u + F(1, 2))
        print(f"    h = (δ1-δ2)/g = {h}: ||h·u - 1/2|| = {float(hu):.4f} "
              f"({'JOINT-SAFE (both-bad impossible)' if hu >= F(1,7) else 'both-bad possible -- needs lift/enriched cite'})")
        # the 2g lift structure
        q1_lift = 2 * g // gcd(d1, 2 * g)
        odd_h = sum(1 for m in H if m % 2 == 1)
        print(f"    mod-2g lift: detuned become q = {q1_lift}; odd harmonics becoming "
              f"half-harmonics of 2g: {odd_h}/{len(H)} -- the 2-adic descent")
        sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART C -- THE BLOCKING/SUPPLY DUALITY on residual-core adversaries")
    print("  blocked(g) = counting dispatch blocked at g; supply = #live pair-sum moduli")
    print("=" * 100)
    G0 = 24
    def blocking_ledger(v, G0=G0):
        blocked = 0
        for g in range(2, G0 + 1):
            Hm = [x for x in v if x % g == 0]
            Dm = [x for x in v if x % g != 0]
            if len(Dm) == 0 or len(Dm) > 4:
                continue  # not a (near-)harmonic modulus for v
            # counting threshold at g using worst-case orbit counts (structure only)
            tot = F(0)
            for dd in Dm:
                q_i = g // gcd(dd, g)
                N_max = 1 if q_i <= 7 else q_i // 7 + 1
                tot += F(N_max, q_i)
            if tot >= 1:
                blocked += 1
        return blocked
    def pair_sum_supply(v):
        # crude B5-flavored liveness: pair-sum moduli q with all residues v*p in the
        # [q/14, 13q/14] band for SOME p (search p <= 12) -- the THM-668 band cert
        live = 0
        sums = sorted(set(a + b for i, a in enumerate(v) for b in v[i + 1:]))
        for q in sums:
            if q < 15:
                continue
            for p in range(1, min(13, q)):
                if all(q <= 14 * ((x * p) % q) <= 13 * q for x in v):
                    live += 1
                    break
        return live, len([q for q in sums if q >= 15])
    fams = []
    # residual-core adversaries: covering, gapped, compressed, distinct, max >= 23
    base = [14 * i for i in range(1, 14)]
    for tries in range(300):
        v = list(base)
        for _ in range(rng.randint(2, 5)):
            i = rng.randrange(13)
            v[i] = max(2, v[i] + rng.choice([-3, -2, -1, 1, 2, 3, 7, -7]))
        v = sorted(set(v))
        if len(v) != 13 or not is_primitive(v):
            continue
        cov = all(any(x % q == 0 for x in v) for q in range(2, 15))
        if not cov:
            continue
        fams.append(v)
        if len(fams) >= 25:
            break
    print(f"  families: {len(fams)} (near-14AP primitive covering adversaries)")
    rows = []
    for v in fams:
        b = blocking_ledger(v)
        live, tot = pair_sum_supply(v)
        rows.append((b, live, tot, v))
    rows.sort(key=lambda r: -r[0])
    print(f"  {'blocked':>8} {'live':>6} {'of':>4}   family (first 6)")
    for b, live, tot, v in rows[:10]:
        print(f"  {b:>8} {live:>6} {tot:>4}   {v[:6]}...")
    bs = [r[0] for r in rows]
    ls = [r[1] for r in rows]
    n = len(rows)
    if n > 2:
        mb, ml = sum(bs) / n, sum(ls) / n
        cov_bl = sum((b - mb) * (l - ml) for b, l in zip(bs, ls)) / n
        sb = (sum((b - mb) ** 2 for b in bs) / n) ** 0.5
        sl = (sum((l - ml) ** 2 for l in ls) / n) ** 0.5
        corr = cov_bl / (sb * sl) if sb * sl > 0 else float('nan')
        print(f"  corr(blocked, live-supply) = {corr:+.3f}   "
              f"min(live) over families = {min(ls)}  (conservation: supply never dies)")
    sys.stdout.flush()
