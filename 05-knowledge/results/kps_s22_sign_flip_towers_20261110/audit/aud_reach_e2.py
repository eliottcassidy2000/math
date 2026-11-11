"""aud_reach_e2.py - independent recomputation of the pre-registered run E2.
(1) the start states: all pairs (B, O), B LABELLED of order <= 7, every O, offset +1 or -1, computed by
    aud_unit_states.c (file aud_unit_states_le7.txt, sorted unique) are compared with runs/e2_states_order_le7.txt;
(2) at most two tower levels (box delta in {1,2}, 0 <= x <= 10, 0 <= m <= 4) with the auditor's level map;
    counts of odd c <= 4097 by sign and depth; the 33 exceptional values of HYP-9658 (B)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
EXC33 = [245, 853, 909, 1049, 1407, 1467, 1515, 1577, 1629, 1829, 1887, 2023, 2037, 2085, 2343, 2393, 2703, 3001, 3033,
         3123, 3313, 3411, 3427, 3545, 3713, 3733, 3753, 3773, 3907, 3909, 3985, 3991, 4023]
mine = sorted(set(tuple(map(int, l.split())) for l in open(os.path.join(HERE, "aud_unit_states_le7.txt")) if l.strip()))
theirs = sorted(set(tuple(map(int, l.split())) for l in open(os.path.join(HERE, "..", "runs", "e2_states_order_le7.txt")) if l.strip()))
print("distinct unit states of order <= 7: auditor %d, packet %d, identical sets: %s" % (len(mine), len(theirs), mine == theirs))
if mine != theirs:
    print("   only auditor:", [s for s in mine if s not in theirs][:10])
    print("   only packet:", [s for s in theirs if s not in mine][:10])


def odd(M):
    v = (M & -M).bit_length() - 1
    return M >> v, v


def level(s, delta, x, m):
    r, g, M, E = s                       # g = h + r
    h = g - r
    al, ga = 2 ** delta + 1, 3 ** delta + 2 ** delta
    Mp = al * 2 ** m * (h + r) + 2 * ga * 2 ** x * 3 ** m * M - (E if x == 0 else 0)
    hp = r + 2 ** (x + 1 + delta + m) * M
    Ep = al * 2 ** x * M + ((r - h) if m == 0 else 0)
    return (-r, hp - r, Mp, Ep)


levels = [(d, x, m) for d in (1, 2) for x in range(11) for m in range(5)]
best = {1: {}, -1: {}}                  # c -> least depth


def rec(s, depth):
    r, g, M, E = s
    if M > 0:
        c, v = odd(M)
        if c <= 4097 and depth < best[r].get(c, 99):
            best[r][c] = depth


frontier = mine
for s in frontier:
    rec(s, 0)
for depth in (1, 2):
    nxt = []
    for s in frontier:
        for (d, x, m) in levels:
            t = level(s, d, x, m)
            rec(t, depth)
            if depth < 2:
                nxt.append(t)
    frontier = nxt
for d in (0, 1, 2):
    print("depth <= %d: offset +1: %d odd c <= 4097 ; offset -1: %d" % (d, sum(1 for c, dd in best[1].items() if dd <= d),
                                                                       sum(1 for c, dd in best[-1].items() if dd <= d)))
print("packet (runs/e2_d2.err): depth <= 0: 89 / 70 ; depth <= 1: 1028 / 1124 ; depth <= 2: 1872 / 1782")
print("exceptional values of HYP-9658 (B) with a seed:", sum(1 for c in EXC33 if c in best[1]), "of 33 ; missing:",
      [c for c in EXC33 if c not in best[1]])
