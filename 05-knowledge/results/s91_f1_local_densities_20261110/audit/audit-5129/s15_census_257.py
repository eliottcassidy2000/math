# audit-5129: the census of order 9 against the product-law lists for EVERY offset D <= 257 (both directions)
import pickle
import numpy as np
import supply as S
rows = pickle.load(open("strong_le257.pkl", "rb"))
Root, PI, Aord, Suf, S0, minord = S.build(rows, 257)
p9 = np.load("pairs9.npy")
D = (p9 >> 32).astype(np.int64); U = (p9 & 0xFFFFFFFF).astype(np.int64)
sel = D <= 257
cen = set(zip(D[sel].tolist(), U[sel].tolist()))
bad = sum(1 for (d, u) in cen if d not in S0 or S.minorder(S0[d], u) is None or S.minorder(S0[d], u) > 9)
pred = set()
for d in S0:
    for s0, o in S0[d].items():
        j = 0
        while o + j <= 9:
            pred.add((d, s0 << j))
            j += 1
print("order-9 census, D <= 257: %d pairs ; not predicted (or predicted only at a larger order): %d ; "
      "predicted with least order <= 9: %d ; predicted but absent: %d" % (len(cen), bad, len(pred), len(pred - cen)))
for d in (129, 257):
    print("offset %d: %d slopes in the census of order <= 9" % (d, sum(1 for (a, b) in cen if a == d)))
# number of strong classes by H for the correlation statement at larger range (bonus)
