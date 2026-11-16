# Port spectra.  Pi(d) = { delta odd, |delta| <= 15 : some tournament T with H(T) = d has a one-vertex extension
# with H = 2d + delta }.  Source: lift_near_order9.txt (every tournament of order <= 9 as base; complete for
# 2d + 15 <= 223, i.e. d <= 104, because realisers of counts <= 223 have strong cores of order <= 10 and a product
# parent q d' of q h' has defect q (h' - 2d'), which the closure of the parent sets shows is never within 15 here).
import math
near = {}
for line in open('lift_near_order9.txt'):
    if line.startswith('NEAR'):
        h, rest = line[5:].split(':'); near[int(h)] = set(int(t) for t in rest.split())
counts = [d for d in range(1, 105, 2) if d not in (7, 21)]
Pi = {d: sorted(dl for dl in range(-15, 16, 2) if (2 * d + dl) in near and dl in near[2 * d + dl]) for d in counts}
with open('port_spectra_104.txt', 'w') as f:
    for d in counts: f.write("%d: %s\n" % (d, ' '.join(map(str, Pi[d]))))
both = [d for d in counts if -1 in Pi[d] and 1 in Pi[d]]
nom = [d for d in counts if -1 not in Pi[d]]; nop = [d for d in counts if 1 not in Pi[d]]
print("counts d <= 104:", len(counts)); print("both digits -1 and +1 available:", len(both))
print("no -1 (2d-1 is not a child):", nom); print("no +1 (2d+1 is not a child):", nop)
print("neither:", [d for d in counts if d in nom and d in nop])
print("Pi for the first counts:", {d: Pi[d] for d in counts[:14]})
print("sizes of Pi: min %d max %d" % (min(len(Pi[d]) for d in counts), max(len(Pi[d]) for d in counts)))
# commutator of the dyadic step D (h -> 2h-1 at a (-1)-port) and the triadic step T (h -> 3h, always physical)
ok = []
for m in counts:
    if 3 * m > 104: break
    td = (-1 in Pi[m]); dt = (-1 in Pi[3 * m]); ok.append((m, 6 * m - 3 if td else None, 6 * m - 1 if dt else None))
print("m : T(D(m)) = 6m-3 , D(T(m)) = 6m-1 (None = the needed (-1)-port does not exist):", ok)
# the physical D/T closure from the count 1 inside counts <= 209, against the arithmetic semigroup <2n, 3n+1> on half-counts
arith = {0}; frontier = [0]
while frontier:
    n = frontier.pop()
    for v in (2 * n, 3 * n + 1):
        if v <= 104 and v not in arith: arith.add(v); frontier.append(v)
phys = {1}; frontier = [1]
while frontier:
    h = frontier.pop()
    nxt = [3 * h]
    if h <= 104 and -1 in Pi.get(h, []): nxt.append(2 * h - 1)
    for v in nxt:
        if v <= 209 and v not in phys: phys.add(v); frontier.append(v)
ar_counts = sorted(2 * n + 1 for n in arith)
print("arithmetic semigroup <2n, 3n+1> from 0, as counts <= 209:", ar_counts)
print("physically reached by T and (-1)-tail steps:", sorted(phys))
print("in the semigroup but not reached:", sorted(set(ar_counts) - phys))
# Erdos-Klarner-Rado exponent
lo, hi = 0.5, 1.0
for _ in range(60):
    mid = (lo + hi) / 2
    if 2 ** -mid + 3 ** -mid > 1: lo = mid
    else: hi = mid
print("sigma with 2^-s + 3^-s = 1: %.6f" % lo)
S = {0}; fr = [0]; X = 10 ** 7
while fr:
    n = fr.pop()
    for v in (2 * n, 3 * n + 1):
        if v <= X and v not in S: S.add(v); fr.append(v)
for e in (4, 5, 6, 7):
    c = sum(1 for v in S if v <= 10 ** e); print("  semigroup elements <= 10^%d: %d  (log count / log X = %.4f)" % (e, c, math.log(c) / math.log(10 ** e)))
