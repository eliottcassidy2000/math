#!/usr/bin/env python3
# mac-mini-2026-11-08-firstrealizers: the eternal colour of every odd h, i.e. whether h is the Hamiltonian-path count of
# SOME self-converse tournament (of any order).
# Criterion (palindromic condensation; THM-4866(E), THM-4743 (1b), Session24 eq. (2), Session224 sec. 5):
#   h is self-converse-attainable  <=>  h = d^2 e with d attainable (d = 1 allowed) and e = 1 or e the H of a STRONG
#   self-converse tournament (the centre of the palindrome).
# Decidability (KNOWN, Session24 PROOF sec. 1; Busch's minimum, THM-4744, CITED: f(3q)=3*5^(q-1), f(3q+1)=5^q,
#   f(3q+2)=9*5^(q-1)): a strong tournament with H = e has order <= B(e) = max{n : f(n) <= e}. Even-order self-converse
#   tournaments have H = 1 (mod 4) (THM-791(a) via the anti-involution). So e is decided once the strong self-converse
#   inventory of every order m <= B(e) (m odd, or e = 1 mod 4) is complete at e.
# Inputs (directory argv[1]); "complete to X" = every self-converse tournament of that order with H <= X is listed:
#   cls_n.txt, n <= 9           per-class census (THM-4922's classgen.c: 'code H c3 strong SC aut ccode | ...'); complete
#   sc_even10.out               order 10 ('S h strong'; sc_even10.c, all 2^25 anti-invariant labellings); complete
#   sc_odd11.out                order 11 ('E h strong'; sc_odd11.c / sc_odd11_fast.c); complete to 624
#   sc_full1011_le2001.out      (optional) orders 10 ('E') and 11 ('O') from sc_odd13.c -DM=4 with no pruning, lines
#                               with H <= 2001 kept; complete to 2001
#   sc_even12_377.out           (optional) order 12 ('P h strong'; sc_even_pair.c); complete to 377
#   sc_odd13_1124.out           (optional) orders 12 ('E') and 13 ('O') from sc_odd13.c -DM=5; complete to 1124
#   sc_witness.out, sc_walk.out (optional) explicit self-converse witnesses ('W h n strong ...'); witnesses only
#   audit_K1/out/o11_full.txt, audit_K1/out/o12_lt625.txt  (optional) audit K1's complete order 11 and order 12 below
#                               625 ('h labelled strong'); used only for consistency checks
# Output: sc_holes.out
import sys, os, collections
D = sys.argv[1] if len(sys.argv) > 1 else "."
INF = float("inf")
def f(n):
    q, r = divmod(n, 3)
    return {0: 3 * 5 ** (q - 1), 1: 5 ** q, 2: 9 * 5 ** (q - 1)}[r]
def B(e): return max(n for n in range(3, 40) if f(n) <= e)
spec = collections.defaultdict(set)                           # all values per order (census only)
SC = collections.defaultdict(dict)                            # SC[m][h] = 1 if some realisation is strong, else 0
comp = collections.defaultdict(int)                           # order m is complete for H <= comp[m]
def add(m, h, s):
    SC[m][h] = SC[m].get(h, 0) | int(s)
def readfile(fn, tagmap, cap=INF):                           # tagmap: line tag -> order
    got = collections.defaultdict(dict)
    for line in open(os.path.join(D, fn)):
        t = line.split()
        if t and t[0] in tagmap and int(t[1]) <= cap: got[tagmap[t[0]]][int(t[1])] = int(t[2])
    return got
for n in range(1, 10):
    for line in open(os.path.join(D, f"cls_{n}.txt")):
        hd = line.split("|")[0].split(); H, strong, isSC = int(hd[1]), int(hd[3]), int(hd[4])
        spec[n].add(H)
        if isSC: add(n, H, strong)
    comp[n] = INF
data = {}
data["sc_even10"] = readfile("sc_even10.out", {"S": 10}); comp[10] = INF
data["sc_odd11"] = readfile("sc_odd11.out", {"E": 11}); comp[11] = 624
opt = [("sc_full1011_le2001.out", {"E": 10, "O": 11}, {10: INF, 11: 2001}),
       ("sc_even12_377.out", {"P": 12}, {12: 377}),
       ("sc_odd13_1124.out", {"E": 12, "O": 13}, {12: 1124, 13: 1124})]
for fn, tags, cmp_ in opt:
    if os.path.exists(os.path.join(D, fn)):
        data[fn] = readfile(fn, tags)
        for m, c in cmp_.items(): comp[m] = max(comp[m], c)
for src in data.values():
    for m, dct in src.items():
        for h, s in dct.items(): add(m, h, s)
# consistency of overlapping complete inventories
print("CONSISTENCY of overlapping inventories (value sets and strong flags on the common complete range):")
def cmp(name, A, Bd, cap):
    a = {h: s for h, s in A.items() if h <= cap}; b = {h: s for h, s in Bd.items() if h <= cap}
    print(f"   {name}: {len(a)} / {len(b)} values, identical: {a == b}")
if "sc_full1011_le2001.out" in data:
    cmp("order 10, sc_even10 vs sc_odd13 -DM=4 (H <= 2001)", data["sc_even10"][10], data["sc_full1011_le2001.out"][10], 2001)
    cmp("order 11, sc_odd11_fast vs sc_odd13 -DM=4 (H < 625)", data["sc_odd11"][11], data["sc_full1011_le2001.out"][11], 624)
if "sc_even12_377.out" in data and "sc_odd13_1124.out" in data:
    cmp("order 12, sc_even_pair vs sc_odd13 (H <= 377)", data["sc_even12_377.out"][12], data["sc_odd13_1124.out"][12], 377)
for fn, m, cap, other in (("audit_K1/out/o11_full.txt", 11, 2001, "sc_full1011_le2001.out"),
                          ("audit_K1/out/o12_lt625.txt", 12, 624, "sc_odd13_1124.out")):
    if os.path.exists(os.path.join(D, fn)) and other in data:
        aud = {}
        for line in open(os.path.join(D, fn)):
            h, c, s = map(int, line.split())
            if h <= cap: aud[h] = int(s > 0)
        cmp(f"order {m}, audit K1 vs {other} (H <= {cap})", aud, data[other][m], cap)
strongSC = {m: {h for h, s in SC[m].items() if s} for m in SC}
Q = set().union(*strongSC.values())                          # strong self-converse values found
birth = {}
for n in sorted(spec):
    for h in spec[n]: birth.setdefault(h, n)
attain = lambda d: d % 2 == 1 and d not in (7, 21)          # all odd d <= 2^32 except 7, 21 (THM-4875, CITED)
def exact_e(e):                                               # is strong-SC membership of e decided by the data?
    return all(comp[m] >= e for m in range(3, B(e) + 1) if not (m % 2 == 0 and e % 4 == 3))
def classify(h):
    """('blue', witness) | ('black', None) exact | ('open', undecided centres)"""
    undecided = []
    for d in range(1, int(h ** 0.5) + 1, 2):
        if h % (d * d) or not attain(d): continue
        e = h // (d * d)
        if e == 1: return "blue", f"{d}^2 sandwich B => B^op, H(B) = {d}"
        if e in Q:
            n0 = min(n for n in strongSC if e in strongSC[n])
            return "blue", (f"{d}^2 * " if d > 1 else "") + f"strong self-converse centre H = {e} (order {n0})"
        if not exact_e(e): undecided.append(e)
    return ("open", undecided) if undecided else ("black", None)
HTOP = 2001
res = {h: classify(h) for h in range(1, HTOP, 2)}
LIMX = next(h for h in range(1, HTOP + 2, 2) if h >= HTOP or res[h][0] == "open")   # every odd h < LIMX is decided
print("\nBusch f(n), n = 3..15:", [f(n) for n in range(3, 16)])
print("inventories complete to:", {m: comp[m] for m in sorted(comp)})
print(f"strong self-converse values known: {len(Q)}; EXACT RANGE: every odd h < {LIMX} is decided")
def report(lo, hi):
    blk = [h for h in range(lo, hi, 2) if res[h][0] == "black"]; op = [h for h in range(lo, hi, 2) if res[h][0] == "open"]
    print(f"\nETERNALLY BLACK odd h in [{lo}, {hi}) (never the H of a self-converse tournament): {len(blk)} values")
    print("  ", blk)
    print("   mod 4:", dict(sorted(collections.Counter(h % 4 for h in blk).items())),
          "; by birth order:", dict(sorted(collections.Counter(birth.get(h, 0) for h in blk).items())))
    if op: print(f"   OPEN: {len(op)}: {op}")
    return blk
b1 = report(1, 625)
print("   below 225:", [h for h in b1 if h < 225], f"({sum(1 for h in b1 if h < 225)});",
      f"[225, 625): {sum(1 for h in b1 if h >= 225)}")
b2 = report(625, LIMX) if LIMX > 625 else []
# birth colour against eternal colour; delays (first self-converse order minus birth order)
def first_sc(h):                                              # least self-converse order, if the inventories decide it
    for m in range(1, 40):
        if h in SC[m]: return m
        if comp[m] < h: return None
def last_complete(h):                                         # largest m with every order <= m complete at h
    m = 0
    while comp[m + 1] >= h: m += 1
    return m
for lo, hi in ((3, 225), (3, 625), (625, LIMX)):
    if hi <= lo: continue
    rows = [h for h in range(lo, hi, 2) if h in birth and birth[h] <= 9]
    bb = [h for h in rows if h in SC[birth[h]]]
    late = [h for h in rows if h not in SC[birth[h]] and res[h][0] == "blue"]
    eb = [h for h in rows if res[h][0] == "black"]
    delay = collections.Counter()
    for h in rows:
        if res[h][0] != "blue": continue
        fs = first_sc(h)
        delay[(fs - birth[h]) if fs is not None else ">" + str(last_complete(h) - birth[h])] += 1
    print(f"\nBIRTH vs ETERNAL on [{lo}, {hi}) (values born by order 9: {len(rows)}; born later: "
          f"{sum(1 for h in range(lo, hi, 2) if attain(h) and (h not in birth or birth[h] > 9))}): blue at birth {len(bb)}, "
          f"black at birth {len(rows) - len(bb)} = eternally black {len(eb)} + later blue {len(late)}")
    print("   delays (first self-converse order - birth order):", dict(sorted(delay.items(), key=lambda kv: str(kv[0]).zfill(4))))
    big = [(h, birth[h], first_sc(h), res[h][1]) for h in late if first_sc(h) is None or first_sc(h) - birth[h] >= 4]
    if big: print("   delay >= 4 or beyond the inventories:", big)
comp_only = [h for h in range(3, LIMX, 2) if res[h][0] == "blue" and h not in Q]
print(f"\nCOMPOSITE-ONLY blue values below {LIMX} (blue, but not the H of any strong self-converse tournament):", comp_only)
# witnesses above the exact range
known = {h for m in SC for h in SC[m]} | {1}
srcs = collections.Counter()
for fn in ("sc_witness.out", "sc_walk.out"):
    if not os.path.exists(os.path.join(D, fn)): continue
    for line in open(os.path.join(D, fn)):
        if line.startswith("W "):
            h = int(line.split()[1])
            if h not in known: srcs[fn] += 1
            known.add(h)
changed = True
while changed:
    changed = False
    for e in sorted(known):
        for d in range(3, int((HTOP / e) ** 0.5) + 1, 2):
            if attain(d) and d * d * e <= HTOP and d * d * e not in known: known.add(d * d * e); changed = True
lo = max(LIMX, 625)
cand = [h for h in range(lo, HTOP, 2) if h not in known]
print(f"\n[{lo}, {HTOP}): blue by explicit witness {sum(1 for h in range(lo, HTOP, 2) if h in known)} of {len(range(lo, HTOP, 2))}; "
      f"unresolved {len(cand)} (witness files added {dict(srcs)} beyond the inventories)")
print("   unresolved:", cand)
print("   unresolved mod 4:", dict(sorted(collections.Counter(h % 4 for h in cand).items())))
nb = len(range(625, HTOP, 2)); ub = len(b2) + len(cand)
print(f"   on [625, {HTOP}): eternally black (exact, below {LIMX}) {len(b2)} + unresolved {len(cand)} = at most {ub} of {nb} "
      f"= {ub / nb:.4f}")
