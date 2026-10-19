import sys, math, json
from collections import defaultdict, Counter
rows = {}
for line in open("mersenne_horizon_50000.txt"):
    K, odd, ter, lj, peak, pay = map(int, line.split()); rows[K] = (odd, ter, lj, peak, pay)
KMAX = max(rows); sig = {K: (rows[K][1]-K, rows[K][0]) for K in rows}
clusters = defaultdict(list)
for K in sorted(rows): clusters[sig[K]].append(K)
cl_of = {K: clusters[sig[K]] for K in rows}
# terminal channel per cluster
chan = Counter(); mixed = 0
for C in clusters.values():
    js = {rows[k][2] for k in C}
    if len(js) > 1: mixed += 1
    chan[min(js)] += 1
print("terminal channel per cluster (j: 4->5, 8->85, 10->341):", sorted(chan.items()), "clusters with mixed channels:", mixed)
# holes
tot_span = sum(max(C)-min(C)+1 for C in clusters.values()); tot_size = sum(len(C) for C in clusters.values())
print(f"total span {tot_span}, total size {tot_size}, hole fraction {1-tot_size/tot_span:.4f}")
# cluster size vs sqrt(K): per cluster (min K), size/sqrt(min)
ratios = [(len(C)/math.sqrt(min(C)), min(C), len(C)) for C in clusters.values() if min(C) >= 1000]
ratios.sort()
print("size/sqrt(Kmin) for clusters with Kmin>=1000: median %.3f, mean %.3f, max %.3f (at Kmin=%d, size %d)" % (ratios[len(ratios)//2][0], sum(r[0] for r in ratios)/len(ratios), ratios[-1][0], ratios[-1][1], ratios[-1][2]))
# least-shift law
leastD = Counter()
for K in range(3, KMAX+1, 2):
    sm = [k for k in cl_of[K] if k < K]
    if sm: leastD[K-max(sm)] += 1
tot = sum(leastD.values())
print("least shift D distribution (odd K with a smaller mate, K<=50000): total", tot)
xs=[];ys=[]
for D in sorted(leastD):
    if D <= 41: print(f"  D={D}: {leastD[D]} ({leastD[D]/tot:.5f})", end=";")
    if 3 <= D <= 41 and leastD[D] >= 5: xs.append(math.log(D)); ys.append(math.log(leastD[D]))
print()
n=len(xs); mx=sum(xs)/n; my=sum(ys)/n; sl=sum((x-mx)*(y-my) for x,y in zip(xs,ys))/sum((x-mx)**2 for x in xs)
print(f"power-law fit of least-shift counts for 3<=D<=41: exponent {sl:.3f}")
print("max least shift:", max(leastD), " shifts > 100:", sum(v for D,v in leastD.items() if D>100))
# orphans beyond 12800
orph = [K for K in range(12801, KMAX+1, 2) if min(cl_of[K]) == K]
print(f"orphans in (12800, 50000]: {len(orph)}; list: {orph}")
# orphan rate pooled by two-run length J = 1 + floor(v2(K-1)/2)
def v2(x): return (x & -x).bit_length() - 1
byJ = defaultdict(lambda: [0,0])
for K in range(1001, KMAX+1, 2):
    J = 1 + v2(K-1)//2; J = min(J, 4); byJ[J][0] += 1
    if min(cl_of[K]) == K: byJ[J][1] += 1
print("orphan rate by two-run length J (K in [1001,50000]):", {J: (n_, o_, round(o_/n_,4)) for J,(n_,o_) in sorted(byJ.items())})
# orphan fraction * sqrt(K) in sliding windows of 5000
print("orphan fraction * sqrt(Kmid) in consecutive windows of 5000 (odd K):")
for lo in range(5001, KMAX, 5000):
    hi = min(lo+4999, KMAX); odds = list(range(lo, hi+1, 2)); o = sum(1 for K in odds if min(cl_of[K]) == K)
    print(f"  {lo}-{hi}: {o}/{len(odds)} = {o/len(odds):.4f}, *sqrtK = {o/len(odds)*math.sqrt((lo+hi)/2):.3f}")
