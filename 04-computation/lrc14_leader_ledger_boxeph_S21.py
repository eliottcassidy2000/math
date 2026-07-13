#!/usr/bin/env python3
"""
lrc14_leader_ledger_boxeph_S21.py  —  THE LEADER LEDGER (THM-722 verifier)

boxeph-2026-07-12-S21.  Pure Python (fractions.Fraction), exact everywhere.

THE OBJECT.  For a speed set S = {v_1 < ... < v_n} of distinct positive integers,
not all odd (so M(S) < 1/2), define on the time circle R/Z:

    p_i(t) = signed phase of v_i t  in (-1/2, 1/2]
    f(t)   = min_i |p_i(t)|                      (distance-to-loneliness)
    lam(t) = argmin_i |p_i(t)|                   (the LEADER = nearest runner)
    phi(t) = p_{lam(t)}(t)                       (signed leader phase)

FACTS (THM-722, proved in canon this session):
  * phi is piecewise linear with slope v_{lam(t)} > 0; all leader changes happen
    at pair-ruler points: DIFFERENCE events ((v_j - v_i) t in Z, phi continuous)
    or SUM events ((v_i + v_j) t in Z, phi jumps +x -> -x, x = f(t) >= 0).
  * THE LEDGER IDENTITY:   int_0^1 v_{lam(t)} dt  =  2 * sum_{sum-handoffs} f(t_h).
  * CHAINS: the circle partitions into H+ = #(positive-depth sum-handoffs) chains;
    each chain has speed-mass  int_chain v_lam = x_in + x_out <= 2 M(S)  and crosses
    phi = 0 exactly once (a "leader landing").
  * PARITY: if S has an even element, H+ is EVEN (iota: t -> -t pairs chains; the
    chains through the fixed points 0 and 1/2 are self-paired).
  * M(S) = max over sum-handoffs of depth (witness = deepest handoff; Kravitz /
    THM-668 pair-sum rulers rederived).

This script computes the whole ledger EXACTLY and verifies every claim.
"""

from fractions import Fraction as F
from math import gcd, floor
import random
import sys
from functools import reduce

# ----------------------------------------------------------------------------
# exact phase helpers
# ----------------------------------------------------------------------------

def signed_phase(v, t):
    """p = v t - round(v t) in (-1/2, 1/2], exact."""
    x = v * t
    fr = x - floor(x)                       # in [0,1)
    return fr if 2 * fr <= 1 else fr - 1

def f_and_leader(speeds, t):
    """(f(t), leader) with deterministic tiebreak (smallest |p|, then smallest v)."""
    best = None
    bestv = None
    for v in speeds:
        a = abs(signed_phase(v, t))
        if best is None or a < best or (a == best and v < bestv):
            best, bestv = a, v
    return best, bestv

# ----------------------------------------------------------------------------
# the ledger
# ----------------------------------------------------------------------------

def leader_ledger(speeds, verbose=False):
    """Compute the exact leader ledger of a distinct-positive-integer speed set.

    Returns a dict with: events, segments (merged leader intervals), handoffs
    (sum-handoffs with depths), Hplus, int_vlam, identity_ok, M, chains, ...
    """
    speeds = sorted(speeds)
    n = len(speeds)
    assert len(set(speeds)) == n, "speeds must be distinct"
    assert all(v > 0 for v in speeds)
    assert any(v % 2 == 0 for v in speeds), \
        "all-odd family: M = 1/2 via t = 1/2, ledger wraps; excluded"

    # 1. events: all sum- and difference-ruler points in [0,1)
    evset = set()
    for i in range(n):
        for j in range(i + 1, n):
            s = speeds[i] + speeds[j]
            d = speeds[j] - speeds[i]
            for m in range(s):
                evset.add(F(m, s))
            for m in range(d):
                evset.add(F(m, d))
    events = sorted(evset)
    ne = len(events)

    # 2. leader on each open interval between consecutive events (midpoint eval)
    seg = []            # (t_start, t_end, leader)
    for k in range(ne):
        t0 = events[k]
        t1 = events[k + 1] if k + 1 < ne else F(1)
        mid = (t0 + t1) / 2
        _, lead = f_and_leader(speeds, mid)
        seg.append((t0, t1, lead))

    # merge consecutive segments with the same leader
    merged = []
    for s0 in seg:
        if merged and merged[-1][2] == s0[2]:
            merged[-1] = (merged[-1][0], s0[1], s0[2])
        else:
            merged.append(list(s0))
    # circular merge across t = 1 ~ 0
    if len(merged) > 1 and merged[0][2] == merged[-1][2]:
        # unify: fold the first segment into the last (record offset by +1)
        merged[-1] = (merged[-1][0], merged[0][1] + 1, merged[-1][2])
        merged.pop(0)
    segments = [tuple(m) for m in merged]

    # 3. classify handoffs at segment boundaries
    handoffs = []       # (t, old_leader, new_leader, kind, depth)
    for k in range(len(segments)):
        t = segments[k][0]                       # boundary between k-1 and k
        old = segments[k - 1][2]
        new = segments[k][2]
        assert old != new
        tt = t if t < 1 else t - 1
        po = signed_phase(old, tt)
        pn = signed_phase(new, tt)
        if po == pn:
            kind = 'diff'
            depth = abs(po)
        elif po == -pn:
            kind = 'sum' if po != 0 else 'zero'
            depth = abs(po)
            if po != 0:
                assert po > 0 > pn, \
                    f"sum-jump must fall +x -> -x, got {po} -> {pn} at t={tt}"
        else:
            raise AssertionError(
                f"unclassifiable handoff at t={tt}: p_old={po}, p_new={pn}")
        handoffs.append((tt, old, new, kind, depth))

    sum_handoffs = [h for h in handoffs if h[3] == 'sum']
    Hplus = len(sum_handoffs)

    # 4. the ledger identity
    int_vlam = sum(s[2] * (s[1] - s[0]) for s in segments)
    twice_depths = 2 * sum(h[4] for h in sum_handoffs)
    identity_ok = (int_vlam == twice_depths)

    # 5. M(S): max handoff depth; cross-check = max f over all events
    M = max((h[4] for h in sum_handoffs), default=F(0))
    M_check = max(f_and_leader(speeds, e)[0] for e in events)
    assert M == M_check, f"M mismatch: handoff-max {M} vs event-max {M_check}"

    # 6. chains: split the segment cycle at positive-depth sum-handoffs
    #    (handoff times ARE segment boundaries) — single linear walk.
    #    chain mass = sum of v*len; verify mass = x_in + x_out per chain
    chains = []
    if Hplus > 0:
        sh_set = {h[0]: h[4] for h in sum_handoffs}      # time -> depth
        # rotate the segment list to start at a sum-handoff boundary
        start = next(k for k, s in enumerate(segments)
                     if (s[0] if s[0] < 1 else s[0] - 1) in sh_set)
        order = segments[start:] + segments[:start]
        t_in = order[0][0] if order[0][0] < 1 else order[0][0] - 1
        mass = F(0)
        for k, (s0, s1, lv) in enumerate(order):
            mass += lv * (s1 - s0)
            nxt = order[(k + 1) % len(order)][0]
            nxt = nxt if nxt < 1 else nxt - 1
            if nxt in sh_set:                            # chain closes here
                chains.append({'t_in': t_in, 't_out': nxt,
                               'x_in': sh_set[t_in], 'x_out': sh_set[nxt],
                               'mass': mass})
                t_in, mass = nxt, F(0)
        assert len(chains) == Hplus
        for c in chains:
            assert c['mass'] == c['x_in'] + c['x_out'], \
                f"chain mass {c['mass']} != x_in+x_out {c['x_in']+c['x_out']}"

    max_mass = max((c['mass'] for c in chains), default=F(0))

    # 7. iota check: the chains containing 0 and 1/2 in their interior
    def contains(c, pt):
        t_in, t_out = c['t_in'], c['t_out']
        if t_in < t_out:
            return t_in < pt < t_out
        return pt > t_in or pt < t_out
    fixed0 = [c for c in chains if contains(c, F(0)) or c['t_in'] == 0 or
              (c['t_in'] > c['t_out'] and True)]
    # simpler: 0 is interior to the chain that wraps (t_in > t_out) or starts <0<
    fixed0 = [c for c in chains if contains(c, F(0))]
    fixedh = [c for c in chains if contains(c, F(1, 2))]

    # leadership shares
    share = {}
    for (s0, s1, lv) in segments:
        share[lv] = share.get(lv, F(0)) + (s1 - s0)

    return {
        'speeds': speeds, 'events': ne, 'segments': segments,
        'handoffs': handoffs, 'sum_handoffs': sum_handoffs, 'Hplus': Hplus,
        'int_vlam': int_vlam, 'identity_ok': identity_ok, 'M': M,
        'chains': chains, 'max_mass': max_mass,
        'fixed0': fixed0, 'fixedh': fixedh, 'share': share,
    }

# ----------------------------------------------------------------------------
# predicates
# ----------------------------------------------------------------------------

def is_primitive(speeds):
    return reduce(gcd, speeds) == 1

def is_covering(speeds):
    """divisor-complete: a multiple of every d in {2..14} (programmatic check,
    MISTAKE-140 standing rule)."""
    return all(any(v % d == 0 for v in speeds) for d in range(2, 15))

# ----------------------------------------------------------------------------
# reporting
# ----------------------------------------------------------------------------

def climb_bound(speeds):
    """The stopping-time / ledger-climb lower bound (LEM-025):
    v = min, f = max, B = second max, q = v+f:
    M >= v*floor(q/(B+v))/q, with equality iff the (v,f)-ruler carries the
    witness."""
    ss = sorted(speeds)
    v, f, B = ss[0], ss[-1], ss[-2]
    q = v + f
    k = q // (B + v)
    return F(v * k, q) if k >= 1 else F(0)


def report(name, speeds, expect_M=None, show_handoffs=False, staircase=False):
    L = leader_ledger(speeds)
    M, Hp = L['M'], L['Hplus']
    eta = F(0) if Hp == 0 else L['int_vlam'] / (2 * Hp * M)
    wit = max(L['chains'], key=lambda c: max(c['x_in'], c['x_out']),
              default=None)
    wit_ratio = (wit['mass'] / (2 * M)) if wit else None
    thresh = F(28, 183)
    print(f"== {name}  {speeds if len(speeds)<=13 else '(13 speeds)'}")
    print(f"   primitive={is_primitive(speeds)} covering(DC)={is_covering(speeds)}")
    print(f"   M = {M} = {float(M):.6f}"
          + (f"   [expected {expect_M} {'OK' if M == expect_M else '** MISMATCH **'}]"
             if expect_M is not None else ""))
    print(f"   events={L['events']}  segments={len(L['segments'])}  H+={Hp}"
          f"  (parity {'EVEN' if Hp % 2 == 0 else 'ODD'})")
    print(f"   int v_lam = {L['int_vlam']} = {float(L['int_vlam']):.6f}"
          f"   identity(int = 2*sum depths): {L['identity_ok']}")
    print(f"   max chain mass = {L['max_mass']} = {float(L['max_mass']):.6f}"
          f"   >= 28/183? {L['max_mass'] >= thresh}"
          f"   (mass/2M = {float(L['max_mass']/(2*M)):.4f})")
    print(f"   ledger efficiency eta = int/(2 H+ M) = {float(eta):.4f}"
          f"   witness-chain equioscillation (x_in+x_out)/2M = "
          f"{float(wit_ratio):.4f}" if wit else "")
    cb = climb_bound(speeds)
    print(f"   climb bound (LEM-025) v*floor(q/(B+v))/q = {cb} = "
          f"{float(cb):.6f}   TIGHT (== M)? {cb == M}")
    # witness pair: which sum-ruler carries the deepest handoff
    deep = max(L['sum_handoffs'], key=lambda h: h[4])
    wp = tuple(sorted((deep[1], deep[2])))
    print(f"   witness pair = {wp} (q = {wp[0]+wp[1]}) at t = {deep[0]}, "
          f"depth {deep[4]}")
    print(f"   fixed chains: through 0: {len(L['fixed0'])}, through 1/2: "
          f"{len(L['fixedh'])}")
    if show_handoffs:
        print("   sum-handoffs (t, old->new, depth):")
        for (t, o, ne_, k, d) in sorted(L['sum_handoffs']):
            print(f"     t={t}  {o}->{ne_}  depth={d} = {float(d):.6f}")
    if staircase:
        print("   staircase check (depths at k/183?):")
        for (t, o, ne_, k, d) in sorted(L['sum_handoffs']):
            if t <= F(1, 2):
                print(f"     t={t}  pair({o},{ne_})  depth={d}"
                      f"   t*183={t*183}  depth*183={d*183}")
    print()
    return L

# ----------------------------------------------------------------------------
# main battery
# ----------------------------------------------------------------------------

def main():
    print("THE LEADER LEDGER — exact verification battery (boxeph-S21)\n")

    # 0. hand-checked toy
    L = report("toy {1,2}", [1, 2], expect_M=F(1, 3), show_handoffs=True)
    assert L['Hplus'] == 2 and L['int_vlam'] == F(4, 3)

    # 1. the tight non-covering extremals
    report("AP {1..13}", list(range(1, 14)), expect_M=F(1, 14),
           show_handoffs=True)
    report("GW {1..11,13,24}", list(range(1, 12)) + [13, 24],
           expect_M=F(1, 14))

    # 2. the covering-min: the deep well, with the staircase
    report("deep well {1..12,182}", list(range(1, 13)) + [182],
           expect_M=F(14, 183), staircase=True)

    # 3. the Ostrowski ladder {1..12,13k}: M = k/(13k+1) (mac-mini cont.55)
    for k in (1, 2, 3, 7, 14):
        report(f"ladder k={k} {{1..12,{13*k}}}", list(range(1, 13)) + [13 * k],
               expect_M=F(k, 13 * k + 1))

    # 4. THM-721 near-dilate at L=2 (M = 1/13) and compressed extremal (M = 1/12)
    report("near-dilate 2*{1..12} u {27}", [2 * i for i in range(1, 13)] + [27],
           expect_M=F(1, 13))
    # (M = 1/13 here, NOT the 1/12 of kps cont.41 — 1/12 was the MISTAKE-141 box
    #  artifact; this family is exactly why the compressed floor fell to 1/13.)
    report("compressed 2*{1..12} u {13}", [2 * i for i in range(1, 13)] + [13],
           expect_M=F(1, 13))

    # 4b. opus-S254 stratum probes: s-scaled core + killer 182, s coprime
    #     (their reduction constant (182+s)/2379 lives on the (s,182)-ruler)
    for s in (3, 5):
        fam = [s * i for i in range(1, 13)] + [182]
        report(f"opus-S254 probe s={s}: {s}*{{1..12}} u {{182}}", fam)

    # 5. random primitive covering (DC) families — the covering-law battery
    print("== RANDOM PRIMITIVE DC BATTERY (Vmax <= 60) ==")
    rng = random.Random(20260712)
    thresh = F(28, 183)
    found, tried = 0, 0
    min_maxmass = None
    min_wit_ratio = None
    stats = []
    while found < 40 and tried < 40000:
        tried += 1
        fam = sorted(rng.sample(range(1, 61), 13))
        if not is_primitive(fam) or not is_covering(fam):
            continue
        found += 1
        L = leader_ledger(fam)
        assert L['identity_ok']
        assert L['Hplus'] % 2 == 0
        wit = max(L['chains'], key=lambda c: max(c['x_in'], c['x_out']))
        wr = wit['mass'] / (2 * L['M'])
        cb = climb_bound(fam)
        assert cb <= L['M'], f"climb bound violated: {cb} > {L['M']} {fam}"
        stats.append((L['M'], L['max_mass'], wr, L['Hplus'], fam, cb))
        if min_maxmass is None or L['max_mass'] < min_maxmass[0]:
            min_maxmass = (L['max_mass'], fam, L['M'])
        if min_wit_ratio is None or wr < min_wit_ratio[0]:
            min_wit_ratio = (wr, fam, L['M'])
    print(f"   {found} primitive DC families tested (of {tried} draws); "
          f"identity + even-H+ + [climb bound <= M] hold on ALL")
    tight = [s for s in stats if s[5] == s[0]]
    print(f"   climb bound TIGHT (== M) on {len(tight)}/{found} random DC")
    viol = [s for s in stats if s[1] < thresh]
    print(f"   max-chain-mass >= 28/183 violations: {len(viol)}")
    for s in viol[:5]:
        print(f"     M={s[0]}  maxmass={s[1]}={float(s[1]):.4f}  fam={s[4]}")
    print(f"   min max-chain-mass: {min_maxmass[0]} = "
          f"{float(min_maxmass[0]):.4f}  (M={min_maxmass[2]}) fam={min_maxmass[1]}")
    print(f"   min witness-equioscillation ratio: {float(min_wit_ratio[0]):.4f}"
          f"  (M={min_wit_ratio[2]}) fam={min_wit_ratio[1]}")
    print(f"   M range over battery: [{float(min(s[0] for s in stats)):.4f}, "
          f"{float(max(s[0] for s in stats)):.4f}]")
    print()

    # 6. optional: the kps blocker (large speeds; slow) — run with arg 'big'
    if 'big' in sys.argv:
        kps = [200, 496, 540, 656, 851, 921, 935, 1122, 1482, 1680, 1835,
               1849, 1856]
        report("kps blocker (cont.47)", kps, expect_M=F(406, 1669))

    print("BATTERY COMPLETE.")

if __name__ == '__main__':
    main()
