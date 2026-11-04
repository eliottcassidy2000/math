#!/usr/bin/env python3
"""Order <= 10, H <= 1125 strong-class census (strongB1125_classes.txt from ./sp.exe strongBpar 10 1125):
count, |Aut| values, values attained only by symmetric strong classes (cumulative by order), |Aut| | H,
and the product closure over asymmetric strong classes on [225, 1125]."""
from collections import Counter, defaultdict
rows = []
for l in open('strongB1125_classes.txt'):
    if l.startswith('#'):
        continue
    n, H, aut, cyc, code = l.split()
    rows.append((int(n), int(H), int(aut), int(cyc)))
c = Counter(n for n, *_ in rows)
print("strong classes with H <= 1125 by order:", sorted(c.items()))
print("total incl. order 1:", len(rows), "  excl. order 1:", len(rows) - 1, "  (author's 'strong classes read': 107773)")
print("|Aut| values:", sorted({a for _, _, a, _ in rows}), "  |Aut| | H violations:", sum(1 for _, H, a, _ in rows if H % a))
print("non-cyclic Aut by order:", sorted(Counter(n for n, H, a, cy in rows if not cy).items()))
for NM in range(7, 11):
    sv = defaultdict(set)
    for n, H, a, cy in rows:
        if n <= NM:
            sv[H].add(a == 1)
    need = sorted(H for H, s in sv.items() if True not in s)
    print("orders <= %2d: strong H values <= 1125 attained only with nontrivial Aut: %s" % (NM, need))
asym = sorted({H for n, H, a, cy in rows if a == 1 and H > 1})
reach = {1}
ch = True
while ch:
    ch = False
    for x in list(reach):
        for v in asym:
            if x * v <= 1125 and x * v not in reach:
                reach.add(x * v); ch = True
print("odd H <= 1125 not in the product closure of asymmetric strong classes of order <= 10:",
      sorted(set(range(1, 1126, 2)) - reach))
