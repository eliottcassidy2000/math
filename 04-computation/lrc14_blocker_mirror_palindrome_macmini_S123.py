#!/usr/bin/env python3
"""HYP-7965 referee — the blocker mirror palindrome (mac-mini-2026-07-19-S123).

Frame = THM-1240 (codex): carrier speed c with complete safe gaps
    G_k(c) = [(14k+1)/(14c), (14k+13)/(14c)],  k = 0..c-1,  center t0 = (k+1/2)/c.
For each fast speed d in the six-pack, the centered spoke phase is
    q = c+d,  p = a nearest integer to q*t0 (BOTH on an exact half-integer
    tie — THM-1240 allows either),  t(d,k) = p/q,
and the blocker set of a choice is
    B(p/q) = { d' in pack, d' != d : ||d' p/q|| < 1/14 }.

CLAIMS (proved by the mirror symmetry, verified exactly here):

  (P1') SET-VALUED PALINDROME LAW: p |-> q-p is a bijection of the
        nearest-integer choices at gap k onto those at gap c-1-k, and
        B(p/q) = B((q-p)/q).  Proof: t0(c-1-k) = 1 - t0(k), so
        |q t0(c-1-k) - (q-p)| = |q t0(k) - p|; and ||w(1-t)|| = ||wt||.
        The gap-indexed word of choice-sets is a palindrome.
        (A FIRST DRAFT claimed the two tie-choices AT THE SAME GAP have
        equal blocker sets — FALSE off-center: c=3, pack {8..13}, k=0,
        d=12 gives q=15, q*t0=5/2, choices 2/15 vs 3/15 with blocker sets
        {8} vs {10}.  The referee below checks the TRUE mirror pairing.)

  (P2)  CENTRAL COLUMN LAW (c odd, d odd): q = c+d even, q*t0 = q/2 is an
        EXACT integer — no tie — the spoke phase is exactly tau = 1/2 and
        B = {even speeds of the pack} \ {d} exactly (||d'/2|| = 0 for even
        d', 1/2 for odd d').  The mirror-fixed column of the palindrome IS
        LEM-020's depth-1 (2-adic) layer.

  (P3') CENTRAL TIES (c odd, d even): q odd, q*t0 = q/2 is a half-integer
        tie, and the two choices (q+-1)/(2q) = 1/2 -+ 1/(2q) are mirror
        images OF EACH OTHER, hence carry EQUAL blocker sets — central
        ties are harmless (this is P1' at the self-mirrored gap).

  (P4') PARITY LOCALIZATION: fix any MIRROR-EQUIVARIANT selection of one
        choice per (d,k) (here: the choice whose blocker set is
        lexicographically smallest; equivariant because mirrored choice
        sets carry equal blocker-set families).  For any gap predicate Q
        on the selected configuration,
            #{k : Q(k)} == Q(central)  (mod 2)   if c odd,
            #{k : Q(k)} == 0           (mod 2)   if c even.
        Any ODD total of gap events forces the event at the central gap —
        onto tau = 1/2, i.e. onto the even speeds.

Also recorded: THM-1240 forced-cycle statistics (gauge b(i) = min blocker on
the selected configuration) and the cycle-length parity profile.
"""

from fractions import Fraction as F
from math import floor

ONE_14 = F(1, 14)


def dist(x: F) -> F:
    fr = x - floor(x)
    return min(fr, 1 - fr)


def round_half_ties(x: F):
    fl = floor(x)
    fr = x - fl
    if fr < F(1, 2):
        return [fl]
    if fr > F(1, 2):
        return [fl + 1]
    return [fl, fl + 1]


def choices(c: int, pack, k: int):
    """dict d -> list of (p, q, frozenset blockers), one per nearest-int."""
    t0 = F(2 * k + 1, 2 * c)
    out = {}
    for d in pack:
        q = c + d
        entries = []
        for p in round_half_ties(q * t0):
            t = F(p, q)
            B = frozenset(d2 for d2 in pack
                          if d2 != d and dist(d2 * t) < ONE_14)
            entries.append((p, q, B))
        out[d] = entries
    return out


def select(entries):
    """Mirror-equivariant selection: lex-min blocker set (then smallest p —
    only reached when both tie-choices have EQUAL blocker sets, so the
    selected blocker set is still mirror-equivariant)."""
    return min(entries, key=lambda e: (sorted(e[2]), e[0]))


def forced_cycle_length(sel):
    labels = sorted(sel)
    nxt = {}
    for d in labels:
        B = sel[d][2]
        if not B:
            return None
        nxt[d] = min(B)
    seen = {}
    v = labels[0]
    i = 0
    while v not in seen:
        seen[v] = i
        v = nxt[v]
        i += 1
    return i - seen[v]


def analyze(name, W, n_fast=6):
    W = sorted(W)
    pack = W[-n_fast:]
    slows = W[:-n_fast]
    print(f"\n-- {name}: pack={pack}, carriers={slows}")
    for c in slows:
        per_gap = [choices(c, pack, k) for k in range(c)]

        # (P1') set-valued palindrome: p <-> q-p matches choices across
        # mirrored gaps with EQUAL blocker sets, choice by choice.
        pal = True
        for k in range(c):
            km = c - 1 - k
            for d in pack:
                A = {(p, q, B) for (p, q, B) in per_gap[k][d]}
                Bm = {(q - p, q, B) for (p, q, B) in per_gap[km][d]}
                if A != Bm:
                    pal = False
        # (P2)+(P3') central column
        central_note = ""
        if c % 2 == 1:
            k0 = (c - 1) // 2
            evens = frozenset(d2 for d2 in pack if d2 % 2 == 0)
            ok2 = all(
                per_gap[k0][d][0][2] == (evens - {d})
                for d in pack if d % 2 == 1
            ) and all(
                len(per_gap[k0][d]) == 1 and per_gap[k0][d][0][0] * 2 ==
                per_gap[k0][d][0][1]
                for d in pack if d % 2 == 1
            )
            ok3 = all(
                len(per_gap[k0][d]) == 2 and
                per_gap[k0][d][0][2] == per_gap[k0][d][1][2]
                for d in pack if d % 2 == 0
            )
            central_note = (f" central k={k0}: odd-spokes at tau=1/2 with "
                            f"evens-block: {'YES' if ok2 else 'NO'}; even-spoke "
                            f"tie-pair equal-B: {'YES' if ok3 else 'NO'}")
            if not (ok2 and ok3):
                raise SystemExit(f"CENTRAL CLAIM FAILED {name} c={c}")

        # (P4') parity localization on the selected configuration
        cyc = []
        tie_count = 0
        for k in range(c):
            sel = {d: select(per_gap[k][d]) for d in pack}
            tie_count += sum(1 for d in pack if len(per_gap[k][d]) == 2)
            cyc.append(forced_cycle_length(sel))
        covered = [k for k, L in enumerate(cyc) if L is not None]
        odd_cyc = [k for k, L in enumerate(cyc) if L is not None and L % 2 == 1]
        if c % 2 == 0:
            par_ok = len(covered) % 2 == 0 and len(odd_cyc) % 2 == 0
        else:
            k0 = (c - 1) // 2
            par_ok = (len(covered) % 2 == (1 if k0 in covered else 0) and
                      len(odd_cyc) % 2 == (1 if k0 in odd_cyc else 0))

        print(f"   c={c:>3}: palindrome={'OK' if pal else 'FAIL'}"
              f"  ties={tie_count:>3}"
              f"  full-blocker gaps={len(covered)}/{c}"
              f"  odd-cycle gaps={len(odd_cyc)}"
              f"  parity-loc={'OK' if par_ok else 'FAIL'}{central_note}")
        if not pal or not par_ok:
            raise SystemExit(f"CLAIM FAILED at {name}, c={c}")


def main():
    fams = [
        ("AP (tight)",              list(range(1, 14))),
        ("GW doubling (tight)",     list(range(1, 12)) + [13, 24]),
        ("deep well",               list(range(1, 13)) + [182]),
        ("slack-1 D=3",             list(range(1, 12)) + [13, 36]),
        ("slack-1 D=2",             list(range(1, 13)) + [26]),
        ("k=10 multi-killer",       list(range(1, 11)) + [13, 22, 84]),
        ("BF38 band family",        [3, 5, 7, 8, 9, 10, 11, 12, 13, 15, 21, 35]),
        ("translate {2..14}",       list(range(2, 15))),
    ]
    for name, W in fams:
        analyze(name, W)
    print("\nAll set-valued palindrome / central-column / parity-localization")
    print("claims hold.  LEM-020's tau <-> 1-tau involution acts on THM-1240's")
    print("blocker data as an exact palindrome of choice-sets; the mirror-fixed")
    print("column (odd carriers) sits at tau = 1/2 where odd spokes are blocked")
    print("by exactly the EVEN pack speeds (depth-1 2-adic layer), and any odd")
    print("count of odd-cycle gaps localizes to that column.")


if __name__ == "__main__":
    main()
