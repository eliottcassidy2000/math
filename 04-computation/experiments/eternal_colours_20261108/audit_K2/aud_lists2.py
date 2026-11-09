#!/usr/bin/env python3
"""audit K2 of THM-4930 statement 5 / HYP-9520 (2026-11-08). Independent recomputation of the eternal-colour lists
below 2001 from the auditor's own complete inventories (written from scratch; extends audit K1's aud_lists.py).

Inputs:
  K1 dir: census_vals_n.txt (n <= 9, re-verified census), anti_10_red.txt (order 10, complete), o11_full.txt (order 11,
          complete, all H)
  K2 dir: o1213_lt1125.txt  aud_sc13.c 13 1125, merged: "E h count strong" (order 12) and "O h count strong" (order 13),
          complete for H <= 1124
  packet dir (comparison and witnesses only): sc_odd13_1124.out, sc_witness.out, sc_walk.out (the witness files were
          re-verified line by line by aud_witness.c / aud_witness16.c)
Decision: h is blue iff h = d^2 e (d odd, d not 7 or 21) with e = 1 or e a strong self-converse value; e is decided when
every order m <= B(e) is complete at e, except even m when e = 3 (mod 4) (THM-791(a)).
usage: aud_lists2.py K1DIR K2DIR PKTDIR
"""
import sys, os, collections
K1, K2, PKT = sys.argv[1], sys.argv[2], sys.argv[3]
INF = float("inf")
def f(m):
    q, r = divmod(m, 3)
    return {0: 3 * 5 ** (q - 1), 1: 5 ** q, 2: 9 * 5 ** (q - 1)}[r]
def B(e):
    m = 3
    while f(m + 1) <= e: m += 1
    return m
SC = collections.defaultdict(dict)        # SC[m][h] = 1 if some realisation is strong
spec = collections.defaultdict(set)
comp = collections.defaultdict(int)
for n in range(1, 10):
    for line in open(os.path.join(K1, f"census_vals_{n}.txt")):
        _, h, s, ss, a = map(int, line.split()); spec[n].add(h)
        if s: SC[n][h] = int(ss > 0)
    comp[n] = INF
for line in open(os.path.join(K1, "anti_10_red.txt")):
    h, c, s = map(int, line.split()); SC[10][h] = int(s > 0)
comp[10] = INF
for line in open(os.path.join(K1, "o11_full.txt")):
    h, c, s = map(int, line.split()); SC[11][h] = int(s > 0)
comp[11] = INF
for line in open(os.path.join(K2, "o1213_lt1125.txt")):
    t, h, c, s = line.split(); m = 12 if t == "E" else 13
    SC[m][int(h)] = int(int(s) > 0)
comp[12] = comp[13] = 1124
print("inventories complete to:", {m: comp[m] for m in sorted(comp)})
# --- comparison with the packet's sc_odd13_1124.out ---
pk = {12: {}, 13: {}}
for line in open(os.path.join(PKT, "sc_odd13_1124.out")):
    t, h, s = line.split(); pk[12 if t == "E" else 13][int(h)] = int(s)
for m in (12, 13):
    mine = {h: s for h, s in SC[m].items() if h <= 1124}
    print(f"order {m}, H <= 1124: auditor {len(mine)} values ({sum(mine.values())} strong); packet {len(pk[m])} "
          f"({sum(pk[m].values())} strong); identical incl. strong flags: {mine == pk[m]}")
    if mine != pk[m]:
        print("   only auditor:", sorted(set(mine) - set(pk[m]))[:20], " only packet:", sorted(set(pk[m]) - set(mine))[:20],
              " flag diffs:", [h for h in mine if h in pk[m] and mine[h] != pk[m][h]][:20])
print("THM-791 at order 12 (all values = 1 mod 4):", all(h % 4 == 1 for h in SC[12]))
strong = {m: sorted(h for h, s in SC[m].items() if s) for m in SC}
print("order 12 strong SC values <= 1124:", len([h for h in strong[12] if h <= 1124]), "least", min(strong[12]))
print("order 13 strong SC values <= 1124:", len([h for h in strong[13] if h <= 1124]), [h for h in strong[13] if h <= 1124])
print("least strong SC value per order m = 3..13 (f_sc) vs Busch f:")
for m in range(3, 14):
    fs = min(strong[m]) if strong.get(m) else None
    print(f"   m = {m:2d}: f_sc = {fs}, f = {f(m)}, equal: {fs == f(m)}, m = 0 mod 6: {m % 6 == 0}")
Q = set().union(*[set(v) for v in strong.values()])
birth = {}
for n in sorted(spec):
    for h in spec[n]: birth.setdefault(h, n)
attain = lambda d: d % 2 == 1 and d not in (7, 21)
def e_decided(e):
    if e == 1: return True
    return all(comp[m] >= e or (m % 2 == 0 and e % 4 == 3) for m in range(3, B(e) + 1))
def colour(h):
    und = []
    d = 1
    while d * d <= h:
        if h % (d * d) == 0 and attain(d):
            e = h // (d * d)
            if e == 1 or e in Q: return "blue"
            if not e_decided(e): und.append(e)
        d += 2
    return "open" if und else "black"
HTOP = 2001
col = {h: colour(h) for h in range(1, HTOP, 2)}
LIMX = next(h for h in range(1, HTOP + 2, 2) if h >= HTOP or col[h] == "open")
print(f"\nEXACT RANGE: every odd h < {LIMX} is decided")
blk = [h for h in range(1, LIMX, 2) if col[h] == "black"]
print(f"eternally black below {LIMX}: {len(blk)}")
b3 = [h for h in blk if h >= 625]
print(f"   in [625, {LIMX}): {len(b3)}: {b3}")
print("   mod 4:", dict(sorted(collections.Counter(h % 4 for h in b3).items())), "; birth orders:",
      dict(sorted(collections.Counter(birth.get(h) for h in b3).items())))
print("   all black below 625 equal K1's 62:", [h for h in blk if h < 625] == [7, 19, 21, 23, 35, 39, 43, 47, 53, 55, 59, 63,
      101, 103, 107, 109, 113, 119, 131, 163, 167, 181, 183, 185, 203, 205, 207, 209, 215, 231, 235, 239, 251, 267, 271, 291,
      299, 301, 303, 311, 317, 321, 327, 331, 335, 347, 359, 377, 415, 427, 431, 439, 451, 463, 483, 487, 503, 515, 519, 523,
      527, 563])
print("   mod 4 of all", len(blk), ":", dict(sorted(collections.Counter(h % 4 for h in blk).items())))
def first_sc(h):
    for m in range(1, 40):
        if h in SC[m]: return m
        if comp[m] < h and not (m % 2 == 0 and h % 4 == 3): return None
    return None
for lo, hi in ((3, 625), (625, LIMX)):
    rows = [h for h in range(lo, hi, 2) if h in birth]
    bb = [h for h in rows if h in SC[birth[h]]]
    late = [h for h in rows if h not in SC[birth[h]] and col[h] == "blue"]
    eb = [h for h in rows if col[h] == "black"]
    dl = collections.Counter((first_sc(h) - birth[h]) if first_sc(h) is not None else None for h in rows if col[h] == "blue")
    print(f"\n[{lo}, {hi}): attainable {len(rows)}; blue at birth {len(bb)}; black at birth {len(rows) - len(bb)} = "
          f"eternally black {len(eb)} + later blue {len(late)}")
    print("   delays:", dict(sorted(dl.items(), key=lambda t: (t[0] is None, t[0] or 0))))
    print("   delay >= 4:", [(h, birth[h], first_sc(h)) for h in rows if col[h] == "blue" and (first_sc(h) is None or first_sc(h) - birth[h] >= 4)])
comp_only = [h for h in range(3, LIMX, 2) if col[h] == "blue" and h not in Q]
print(f"\ncomposite-only blue below {LIMX}:", comp_only)
for h in comp_only:
    ds = [(d, h // (d * d)) for d in range(1, int(h ** 0.5) + 2, 2) if d * d <= h and h % (d * d) == 0]
    info = [(d, e, sorted(m for m in strong if e in strong[m])) for d, e in ds]
    print(f"   {h}: decompositions (d, e, orders with e strong SC) {info}; first SC order {first_sc(h)}; e = h decided: {e_decided(h)}")
print("\n79 and 95 strong SC at orders:", [m for m in strong if 79 in strong[m]], [m for m in strong if 95 in strong[m]],
      "; 1011:", [m for m in SC if 1011 in SC[m]], "strong at", [m for m in strong if 1011 in strong[m]])
print("1125, 1127, 1129 strong SC at order 9:", [h in strong[9] for h in (1125, 1127, 1129)])
# --- [LIMX, 2001): witnesses (re-verified) + complete data + d^2 closure ---
known = {h for m in SC for h in SC[m]} | {1}
wit = collections.Counter()
for fn in ("sc_witness.out", "sc_walk.out"):
    p = os.path.join(PKT, fn)
    if not os.path.exists(p): continue
    for line in open(p):
        if line.startswith("W "):
            h = int(line.split()[1])
            if h not in known: wit[fn] += 1
            known.add(h)
ch = True
while ch:
    ch = False
    for e in sorted(known):
        d = 3
        while d * d * e <= HTOP:
            if attain(d) and d * d * e not in known: known.add(d * d * e); ch = True
            d += 2
unres = [h for h in range(LIMX, HTOP, 2) if h not in known]
print(f"\n[{LIMX}, {HTOP}): odd {len(range(LIMX, HTOP, 2))}; blue {len(range(LIMX, HTOP, 2)) - len(unres)}; unresolved {len(unres)}: {unres}")
print("   unresolved mod 4:", dict(sorted(collections.Counter(h % 4 for h in unres).items())), "; witness files added:", dict(wit))
# --- the walk claims ---
p = os.path.join(PKT, "sc_walk.out")
if os.path.exists(p):
    ws = collections.defaultdict(set)
    for line in open(p):
        if line.startswith("W "):
            t = line.split(); h, n, st = int(t[1]), int(t[2]), int(t[3])
            if st: ws[n].add(h)
    for m in (12, 13):
        ex = {h for h in strong[m] if h <= 1124}; wl = {h for h in ws[m] if h <= 1124}
        print(f"walk vs exact strong inventory, order {m}, H <= 1124: walk {len(wl)}, exact {len(ex)}, missing {sorted(ex - wl)}, extra {sorted(wl - ex)}")
    print("walk strong values by order:", {m: len(ws[m]) for m in sorted(ws)}, "; above 1124 at order 13:", len([h for h in ws[13] if h > 1124]),
          "; order 14 values <= 1300:", sorted(h for h in ws[14] if h <= 1300))
    print("walk hits among the unresolved:", sorted(set(unres) & set().union(*ws.values())))
print("\nshares: [625,1131)", f"{len(b3)}/{len(range(625, 1131, 2))} = {len(b3) / len(range(625, 1131, 2)):.4f};",
      f"[1131,2001) {len(unres)}/{len(range(1131, 2001, 2))} = {len(unres) / len(range(1131, 2001, 2)):.4f}")
