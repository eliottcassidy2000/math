# Complete one-vertex parent sets P(h) = { H(T - x) < h : H(T) = h } for odd h <= 223.
# Census part (lift_par_order9.txt): every tournament T of order <= 10 and every vertex x.
# A realiser of h <= 223 has strong components of order <= 10 (Busch: a strong tournament of order 11 has H >= 225),
# but its total order may exceed 10 when it has several nontrivial components.  For T = Q + T'' (ordinal sum, either
# order) and x in T'':  H(T - x) = H(Q) H(T'' - x).  Hence the closure
#     P(h) = P_census(h)  u  { q d' : h = q h', q > 1 realisable, h' > 1, d' in P(h') } .
import statistics
par = {}
for line in open('lift_par_order9.txt'):
    if line.startswith('PAR'):
        h, rest = line[4:].split(':'); par[int(h)] = set(int(t) for t in rest.split())
real = set(h for h in par if h <= 223) | {1}
P = {}
for h in sorted(x for x in par if x <= 223):
    S = set(par[h])
    for q in range(3, h, 2):
        if h % q == 0 and q in real and (h // q) in P: S |= {q * d for d in P[h // q]}
    P[h] = S
added = {h: sorted(P[h] - par[h]) for h in P if P[h] - par[h]}
print("counts <= 223 with a realiser:", len(P), " (odd numbers 3..223 missing:", [h for h in range(3, 224, 2) if h not in P], ")")
print("product parents not seen at order <= 10:", added)
hs = sorted(P)
npar = [len(P[h]) for h in hs]
print("number of parents: min %d (h=%d), median %s, max %d (h=%d)" % (min(npar), hs[npar.index(min(npar))], statistics.median(npar), max(npar), hs[npar.index(max(npar))]))
print("fewest parents:", sorted(((len(P[h]), h) for h in hs))[:12])
dist = {}
for h in hs:
    if h == 3: continue
    m = min(abs(h - 2 * d) for d in P[h]); dist.setdefault(m, []).append(h)
print("least |h - 2d| over parents (h >= 5):", {k: len(v) for k, v in sorted(dist.items())})
for k, v in sorted(dist.items()):
    if k > 1: print("  defect %d:" % k, v)
def D(h): return (h + 1) // 2 if h % 4 == 1 else ((3 * h + 3) // 4 if h % 8 == 3 else (3 * h - 1) // 4)
halv = [h for h in hs if h % 4 == 1]
print("h = 1 (mod 4): %d counts; D(h) = (h+1)/2 is a parent for %d; not for %s" % (len(halv), sum(D(h) in P[h] for h in halv), [h for h in halv if D(h) not in P[h]]))
ris = [h for h in hs if h % 4 == 3 and h > 3]
print("h = 3 (mod 4), h > 3: %d counts; D(h) is a parent for %d (forbidden by the floor)" % (len(ris), sum(D(h) in P[h] for h in ris)))
two = [h for h in hs if h > 3 and not ({(h - 1) // 2, (h + 1) // 2} & P[h])]
print("neither (h-1)/2 nor (h+1)/2 is a parent: %d counts: %s" % (len(two), two))
rat = sorted((max(P[h]) / h, h) for h in hs if h > 3)
print("largest parent / h: max %.4f at h=%d (parent %d); min %.4f at h=%d" % (rat[-1][0], rat[-1][1], max(P[rat[-1][1]]), rat[0][0], rat[0][1]))
viol = [(h, d) for h in hs for d in P[h] if 2 * h < 3 * d + 1]
print("violations of the floor h >= (3d+1)/2:", viol)
tight = sorted((2 * h - 3 * d, d, h) for h in hs for d in P[h] if 2 * h - 3 * d <= 3)
print("pairs (2h-3d, d, h) with 2h - 3d <= 3:", tight)
with open('parents_closed_223.txt', 'w') as f:
    for h in hs: f.write("%d: %s\n" % (h, ' '.join(map(str, sorted(P[h])))))
