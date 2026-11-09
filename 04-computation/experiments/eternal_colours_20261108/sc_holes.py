#!/usr/bin/env python3
# mac-mini-2026-11-08-firstrealizers: the eternal colour of every odd h, i.e. whether h is the Hamiltonian-path count of
# SOME self-converse tournament (of any order).
# Criterion (palindromic condensation; THM-4866(E), Session24 eq. (2), Session224 sec. 5):
#   h is self-converse-attainable  <=>  h = d^2 e with d attainable (d = 1 allowed) and e = 1 or e the H of a STRONG
#   self-converse tournament (the centre of the palindrome).
# Decidability (Busch's minimum, THM-4744, CITED: f(3q)=3*5^(q-1), f(3q+1)=5^q, f(3q+2)=9*5^(q-1)): a strong tournament
#   with H = e has order <= B(e) = max{n : f(n) <= e}. Even-order self-converse tournaments have H = 1 (mod 4) (THM-791).
# Inputs (directory argv[1]):
#   cls_n.txt, n <= 9   per-class census (classgen.c of THM-4922's packet: 'code H c3 strong SC aut ccode | ...')
#   sc_even10.out       all self-converse tournaments of order 10 ('S h strong'; sc_even10.c)
#   sc_odd11.out        all self-converse tournaments of order 11 with H < 625 ('E h strong'; sc_odd11.c)
#   sc_witness.out      near-transitive mirror-symmetric self-converse witnesses, orders <= 13 ('W h n strong ...')
# Exactness: strong self-converse inventories are complete for e < 375 (orders <= 11 cover B(e) <= 11) and, for
#   e = 3 (mod 4), below 625 (order 12 is even, order >= 13 needs e >= 625).
# Output: sc_holes.out
import sys, os, collections
D = sys.argv[1] if len(sys.argv) > 1 else "."
def f(n):
    q, r = divmod(n, 3)
    return {0: 3 * 5 ** (q - 1), 1: 5 ** q, 2: 9 * 5 ** (q - 1)}[r]
spec = collections.defaultdict(set); SCany = collections.defaultdict(set); strongSC = collections.defaultdict(set)
for n in range(1, 10):
    for line in open(os.path.join(D, f"cls_{n}.txt")):
        hd = line.split("|")[0].split(); H, strong, SC = int(hd[1]), int(hd[3]), int(hd[4])
        spec[n].add(H)
        if SC:
            SCany[n].add(H)
            if strong: strongSC[n].add(H)
for line in open(os.path.join(D, "sc_even10.out")):
    if line.startswith("S "):
        _, h, st = line.split(); SCany[10].add(int(h))
        if int(st): strongSC[10].add(int(h))
for line in open(os.path.join(D, "sc_odd11.out")):
    if line.startswith("E "):
        _, h, st = line.split(); SCany[11].add(int(h))
        if int(st): strongSC[11].add(int(h))
birth = {}
for n in sorted(spec):
    for h in spec[n]: birth.setdefault(h, n)
attain = lambda d: d % 2 == 1 and d not in (7, 21)          # all odd d <= 2^32 except 7, 21 (THM-4875, CITED)
Q = set().union(*strongSC.values())                          # strong self-converse values found (complete as stated)
def exact_e(e):                                               # is strong-SC membership of e decided by the data?
    return e < 375 or (e < 625 and e % 4 == 3)
def classify(h):
    """('blue', witness) | ('black', None) exact | ('open', reason)"""
    undecided = []
    for d in range(1, int(h ** 0.5) + 1, 2):
        if h % (d * d) or not attain(d): continue
        e = h // (d * d)
        if e == 1: return "blue", f"{d}^2 sandwich B => B^op, H(B) = {d}"
        if e in Q:
            n0 = min(n for n in strongSC if e in strongSC[n])
            return "blue", (f"{d}^2 * " if d > 1 else "") + f"strong self-converse centre H = {e} (order {n0})"
        if not exact_e(e): undecided.append(e)
    return ("open", f"strong self-converse centres with H in {undecided} need order 12") if undecided else ("black", None)
LIM = 625
res = {h: classify(h) for h in range(1, LIM, 2)}
black = [h for h in res if res[h][0] == "black"]; openv = [h for h in res if res[h][0] == "open"]
print("Busch f(n), n = 3..13:", [f(n) for n in range(3, 14)])
print(f"self-converse values found: orders 10 / 11 (H < 625): {len(SCany[10])} / {len(SCany[11])}; strong centres known: {len(Q)}")
print(f"\nETERNALLY BLACK odd h < {LIM} (exact; never the H of a self-converse tournament): {len(black)} values")
print("  ", black)
print("   of which below 225:", [h for h in black if h < 225])
print("   mod 4:", dict(collections.Counter(h % 4 for h in black)))
print("   by birth order:", dict(sorted(collections.Counter(birth.get(h, 0) for h in black).items())))
print(f"OPEN (h = 1 mod 4 in [375, 625) needing order-12 data): {len(openv)}: {openv}")
# birth colour against eternal colour
late = []
for h in range(3, LIM, 2):
    if h not in birth or birth[h] > 9: continue
    b = birth[h]
    if h not in SCany[b] and res[h][0] == "blue":
        nsc = min((n for n in SCany if h in SCany[n]), default=None)
        late.append((h, b, nsc, res[h][1]))
print(f"\nBLACK AT BIRTH, BLUE LATER (h < {LIM}, birth order <= 9): {len(late)}")
for row in late: print("  ", row)
tot = sum(1 for h in range(3, LIM, 2) if h in birth and birth[h] <= 9)
bb = sum(1 for h in range(3, LIM, 2) if h in birth and birth[h] <= 9 and h in SCany[birth[h]])
print(f"summary (odd 3 <= h < {LIM} born by order 9): {tot}; blue at birth {bb}; black at birth {tot - bb} = eternally black "
      f"{sum(1 for h in black if h in birth)} + late blue {len(late)} + open {sum(1 for h in openv if h in birth and h not in SCany[birth[h]])}")
# witnesses above the exact range
known = set().union(*SCany.values()) | {1}
for line in open(os.path.join(D, "sc_witness.out")):
    if line.startswith("W "): known.add(int(line.split()[1]))
HTOP = 2001
changed = True
while changed:
    changed = False
    for e in sorted(known):
        for d in range(3, int((HTOP / e) ** 0.5) + 1, 2):
            if attain(d) and d * d * e <= HTOP and d * d * e not in known: known.add(d * d * e); changed = True
cand = [h for h in range(LIM, HTOP, 2) if h not in known]
print(f"\n{LIM} <= h < {HTOP}: blue by explicit witness {sum(1 for h in range(LIM, HTOP, 2) if h in known)}; unresolved {len(cand)}")
print("   unresolved:", cand)
print("   unresolved mod 4:", dict(collections.Counter(h % 4 for h in cand)))
