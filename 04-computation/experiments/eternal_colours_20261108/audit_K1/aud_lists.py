#!/usr/bin/env python3
"""audit K1 of THM-4930 / HYP-9520 (2026-11-08). Independent recomputation of the eternal-colour lists.

Inputs (directory W), all produced by the auditor's own programs:
  census_vals_n.txt (n <= 9)  aud_census.c on THM-4922's class files (every field re-verified class by class):
                              'n h SCclasses strongSCclasses allclasses'
  anti_10_red.txt             aud_sc.c order 10 (all self-converse classes): 'h count strongcount'
  o11_full.txt                aud_sc.c order 11, ALL self-converse classes, no H bound and no monotonicity pruning
  o12_lt625.txt (optional)    aud_sc12.c: all self-converse 12-tournaments with H < 625: 'h count strongcount'
Packet files (for comparison only): sc_odd11.out, sc_even10.out, sc_witness.out (re-checked by aud_witness.c).

The decision uses only:
  (i)  the palindromic criterion: h in Sigma_sc iff h = d^2 e, d odd path count (all odd d except 7, 21; for the d that
       occur here, d <= 45, this is read off the census), e = 1 or e the H of a strong self-converse tournament;
  (ii) Busch: a strong m-tournament has H >= f(m) (THM-4744, CITED); the auditor checks f against the census minima;
  (iii) THM-791(a) (even-order self-converse => H = 1 mod 4), used only for orders that are NOT enumerated.
"""
import sys, os, collections
args = [a for a in sys.argv[1:] if not a.startswith("--")]
NO12 = "--no12" in sys.argv
W = args[0]
PKT = args[1] if len(args) > 1 else None

def f(m):
    q, r = divmod(m, 3)
    return {0: 3 * 5 ** (q - 1), 1: 5 ** q, 2: 9 * 5 ** (q - 1)}[r]

spec = collections.defaultdict(set)      # all H at order n (n <= 9)
SC = collections.defaultdict(set)        # self-converse H at order n
SSC = collections.defaultdict(set)       # strong self-converse H at order n
minstrong = {}
for n in range(1, 10):
    for line in open(os.path.join(W, f"census_vals_{n}.txt")):
        _, h, s, ss, a = map(int, line.split())
        spec[n].add(h)
        if s: SC[n].add(h)
        if ss: SSC[n].add(h)
for line in open(os.path.join(W, "anti_10_red.txt")):
    h, c, s = map(int, line.split()); SC[10].add(h)
    if s: SSC[10].add(h)
for line in open(os.path.join(W, "o11_full.txt")):
    h, c, s = map(int, line.split()); SC[11].add(h)
    if s: SSC[11].add(h)
have12 = os.path.exists(os.path.join(W, "o12_lt625.txt")) and not NO12
print("order-12 data (H < 625) used:", have12)
if have12:
    for line in open(os.path.join(W, "o12_lt625.txt")):
        h, c, s = map(int, line.split()); SC[12].add(h)
        if s: SSC[12].add(h)
ENUM_FULL = set(range(1, 12))            # orders whose self-converse inventory is complete for every H
ENUM_625 = {12} if have12 else set()     # orders complete below 625 only

# Busch check against the data: no strong self-converse value below f(m) at any enumerated order
print("Busch f(m), m = 3..15:", [f(m) for m in range(3, 16)])
for n in range(3, 12):
    if SSC[n]:
        print(f"  order {n}: least strong self-converse H = {min(SSC[n])}, f({n}) = {f(n)}, consistent: {min(SSC[n]) >= f(n)}")

def B(e):
    m = 3
    while f(m + 1) <= e: m += 1
    return m

birth = {}
for n in sorted(spec):
    for h in spec[n]: birth.setdefault(h, n)
attain = lambda d: d % 2 == 1 and d not in (7, 21)
# check the attainability of the d actually used (d^2 < 2001 -> d <= 43) directly from the census
for d in range(1, 46, 2):
    assert attain(d) == (d in birth), d

SSC_all = set().union(*SSC.values())
def e_decided(e, lim_orders=True):
    """membership of e in Sigma_sc^str is decided by the enumerated orders (+ THM-791 for even orders not enumerated)"""
    if e == 1: return True
    for m in range(3, B(e) + 1):
        if m in ENUM_FULL: continue
        if m in ENUM_625 and e < 625: continue
        if m % 2 == 0 and e % 4 == 3: continue           # THM-791(a)
        return False
    return True

def colour(h):
    und = []
    for d in range(1, int(h ** 0.5) + 2, 2):
        if d * d > h or h % (d * d) or not attain(d): continue
        e = h // (d * d)
        if e == 1 or e in SSC_all: return "blue"
        if not e_decided(e): und.append(e)
    return "open" if und else "black"

LIM = 625
col = {h: colour(h) for h in range(1, LIM, 2)}
black = [h for h in col if col[h] == "black"]
openv = [h for h in col if col[h] == "open"]
print(f"\neternally black odd h < {LIM}: {len(black)}")
print("  ", black)
print("  below 225:", len([h for h in black if h < 225]), [h for h in black if h < 225])
print("  in [225, 625):", len([h for h in black if h >= 225]), [h for h in black if h >= 225])
print("  mod 4:", dict(collections.Counter(h % 4 for h in black)))
print("  birth orders:", dict(sorted(collections.Counter(birth.get(h, 0) for h in black).items())))
print("open:", openv)

# exact-range claims: every odd h < 375, and h = 3 mod 4 below 625, decided WITHOUT order 12
saveF = set(ENUM_625); ENUM_625.clear()
und_without12 = [h for h in range(1, LIM, 2) if colour(h) == "open"]
ENUM_625 |= saveF
print("undecided without any order-12 data:", und_without12)
print("  all = 1 mod 4 and >= 375:", all(h % 4 == 1 and h >= 375 for h in und_without12))

# birth vs eternal colours
def first_sc(h):
    return min((n for n in SC if h in SC[n]), default=None)
rows = []
for h in range(3, LIM, 2):
    if h not in birth: continue
    b = birth[h]; bb = h in SC[b]
    rows.append((h, b, bb, col[h], first_sc(h)))
for top in (225, 625):
    R = [r for r in rows if r[0] < top]
    nb = sum(1 for r in R if r[2]); nk = len(R) - nb
    eb = [r for r in R if not r[2] and r[3] == "black"]
    lb = [r for r in R if not r[2] and r[3] == "blue"]
    op = [r for r in R if not r[2] and r[3] == "open"]
    print(f"\nattainable odd 3 <= h < {top}: {len(R)}; blue at birth {nb}; black at birth {nk} = eternally black {len(eb)} "
          f"+ later blue {len(lb)} + open {len(op)}")
    par = [r for r in R if not r[2] and r[1] % 2 == 0 and r[0] % 4 == 3]
    print(f"  parity-forced black at birth (even birth order, h = 3 mod 4): {len(par)}; eternally black "
          f"{sum(1 for r in par if r[3] == 'black')}, later blue {sum(1 for r in par if r[3] == 'blue')}")
    if top == 225:
        print("    parity-forced, eternally black:", [r[0] for r in par if r[3] == "black"])
        print("    parity-forced, later blue:", [(r[0], r[4]) for r in par if r[3] == "blue"])
        npar = [r for r in R if not r[2] and not (r[1] % 2 == 0 and r[0] % 4 == 3)]
        print("    not forced, eternally black:", [r[0] for r in npar if r[3] == "black"])
        print("    not forced, later blue:", [(r[0], r[4]) for r in npar if r[3] == "blue"])
        print("    later blue (h, birth, first SC order):", [(r[0], r[1], r[4]) for r in lb])
    delays = collections.Counter((r[4] - r[1]) if r[4] is not None else None for r in R if r[3] == "blue")
    print("  delays (first SC order - birth order) over eternally blue values:", dict(sorted(delays.items(), key=lambda t: (t[0] is None, t[0] or 0))))
    print("  eternally blue values with no SC realisation at enumerated orders:", [r[0] for r in R if r[3] == "blue" and r[4] is None])

# composite-only blue: eternally blue but not the H of any strong self-converse tournament
comp = []
for h in range(3, LIM, 2):
    if col[h] != "blue" or h in SSC_all: continue
    comp.append((h, "decided" if e_decided(h) else "order-12 strong SC unknown"))
print("\ncomposite-only blue candidates below 625 (blue, no strong SC found at enumerated orders):", comp)

# 275: realisations
print("\n275: decompositions d^2 e:", [(d, 275 // (d * d)) for d in range(1, 17, 2) if 275 % (d * d) == 0])
print("  275 SC at orders:", [n for n in sorted(SC) if 275 in SC[n]], "; 11 strong SC at orders:", [n for n in sorted(SSC) if 11 in SSC[n]],
      "; 5 strong at order 4:", 5 in spec[4])

# order 9 share quoted in HYP-9520
mx9 = max(spec[9]); odd9 = (mx9 + 1) // 2
print(f"\norder 9: max H = {mx9}; odd values 1..{mx9}: {odd9}; self-converse values at order 9: {len(SC[9])}")

# values in [625, 2001): blue by data (orders <= 11 complete, witnesses, closure under d^2) vs unresolved
HTOP = 2001
known = set().union(*SC.values())
if PKT:
    for line in open(os.path.join(PKT, "sc_witness.out")):
        if line.startswith("W "): known.add(int(line.split()[1]))
known.add(1)
ch = True
while ch:
    ch = False
    for e in sorted(known):
        d = 3
        while d * d * e <= HTOP:
            if attain(d) and d * d * e not in known: known.add(d * d * e); ch = True
            d += 2
unres = [h for h in range(LIM, HTOP, 2) if h not in known]
print(f"[625, {HTOP}): odd values {len(range(LIM, HTOP, 2))}; blue by data {len(range(LIM, HTOP, 2)) - len(unres)}; unresolved {len(unres)}")
print("  unresolved:", unres)
print("  unresolved mod 4:", dict(collections.Counter(h % 4 for h in unres)))
if PKT:
    pk = [671, 683, 707, 719, 723, 755, 785, 795, 815, 851, 869, 887, 907, 931, 949, 983, 1003, 1009, 1067, 1087, 1115,
          1131, 1145, 1231, 1245, 1263, 1391, 1459, 1475, 1535, 1551, 1567, 1583, 1627, 1699, 1891]
    newly = sorted(set(pk) - set(unres))
    print("  of the packet's 36 unresolved, now blue:", newly, [ (h, [n for n in sorted(SC) if h in SC[n]]) for h in newly])
    print("  unresolved not in packet's list:", sorted(set(unres) - set(pk)))

# compare with the packet's order-11 list below 625
if PKT:
    pk11 = {}
    for line in open(os.path.join(PKT, "sc_odd11.out")):
        _, h, s = line.split(); pk11[int(h)] = int(s)
    mine = {h: int(h in SSC[11]) for h in SC[11] if h < 625}
    print(f"\norder 11 below 625: auditor {len(mine)} values ({sum(mine.values())} strong); packet {len(pk11)} ({sum(pk11.values())} strong); identical incl. strong flags: {mine == pk11}")
    if mine != pk11:
        print("  only auditor:", sorted(set(mine) - set(pk11)), " only packet:", sorted(set(pk11) - set(mine)),
              " flag diffs:", [h for h in mine if h in pk11 and mine[h] != pk11[h]])
