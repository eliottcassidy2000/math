"""Exact horizon analysis of the Mersenne line M_K = 2^K - 1, K = 2..KMAX, from mersenne_horizon_50000.txt
(columns: K odd terras last_j peak pay).  Clusters = level sets of the deletion invariant (terras - K, odd)
(barriers note Prop. 4.1 / HYP-9242 criterion).  Orphan = odd K that is the minimum of its cluster."""
import sys, math, json
from collections import defaultdict, Counter
fn = sys.argv[1] if len(sys.argv) > 1 else "mersenne_horizon_50000.txt"
rows = {}
for line in open(fn):
    K, odd, ter, lj, peak, pay = map(int, line.split())
    rows[K] = (odd, ter, lj, peak, pay)
KMAX = max(rows)
sig = {K: (rows[K][1] - K, rows[K][0]) for K in rows}
clusters = defaultdict(list)
for K in sorted(rows): clusters[sig[K]].append(K)
cl_of = {K: clusters[sig[K]] for K in rows}
print(f"K range 2..{KMAX}; {len(clusters)} clusters (level sets of (sigma_T - K, odd count))")

# sanity: reset pairs (THM-4556 ii): even K>=4 shares its cluster with K-1
bad = [K for K in range(4, KMAX+1, 2) if sig[K] != sig[K-1]]
print("reset pairs (even K with K-1 in the same cluster): violations =", bad[:5], "count", len(bad))
# least shift to a smaller cluster-mate is odd (THM-4556 iii), over odd K
viol = 0; leastD = Counter()
for K in range(3, KMAX+1, 2):
    C = cl_of[K]; sm = [k for k in C if k < K]
    if sm:
        D = K - max(sm); leastD[D] += 1
        if D % 2 == 0: viol += 1
print("odd K with a smaller cluster-mate: least shift D even in", viol, "cases; least-D histogram (first 12):", sorted(leastD.items())[:12])

def window_stats(lo, hi):
    odds = [K for K in range(lo | 1, hi+1, 2) if K >= 3]
    orph = [K for K in odds if min(cl_of[K]) == K]
    # clusters with minimum in the window
    cmins = [C for C in clusters.values() if lo <= min(C) <= hi]
    mean_size = sum(len(C) for C in cmins) / max(1, len(cmins))
    size_biased = sum(len(cl_of[K]) for K in odds) / len(odds)
    below = sum(sum(1 for k in cl_of[K] if k < K) for K in odds) / len(odds)
    T = sum(rows[K][1] - (K - 1) for K in odds) / len(odds)   # post-run Terras time sigma_T(x_K)
    Kmid = (lo + hi) / 2
    return dict(lo=lo, hi=hi, n_odd=len(odds), orphans=len(orph), frac=len(orph)/len(odds),
                frac_sqrtK=len(orph)/len(odds)*math.sqrt(Kmid), mean_cluster=mean_size,
                two_over_mean_sqrtK=2/mean_size*math.sqrt(Kmid), size_biased=size_biased,
                mean_below=below, mean_below_over_sqrtT=below/math.sqrt(T), T=T)

print("\n== windows of odd K: orphan fraction, mean cluster size (clusters with min in window), size-biased, members below K ==")
print("window         n_odd orphans  frac    frac*sqrtK  meanCl  2/meanCl*sqrtK  sizeBiased  below  below/sqrtT   T(post-run Terras)")
wins = [(3,400),(401,800),(801,1600),(1601,3200),(3201,6400),(6401,12800),(12801,19200),(19201,25600),(25601,38400),(38401,50000),
        (1000,2000),(2000,4000),(4000,8000),(8000,16000),(16000,32000),(32000,50000)]
W = []
for lo, hi in wins:
    s = window_stats(lo, hi); W.append(s)
    print(f"{lo:6d}-{hi:6d} {s['n_odd']:6d} {s['orphans']:6d} {s['frac']:.4f} {s['frac_sqrtK']:8.3f} {s['mean_cluster']:8.2f} {s['two_over_mean_sqrtK']:9.3f} {s['size_biased']:10.1f} {s['mean_below']:7.1f} {s['mean_below_over_sqrtT']:9.4f} {s['T']:10.0f}")

# exponent fit of the orphan fraction: log frac vs log K over dyadic windows K >= 400
import statistics
xs = []; ys = []
for lo, hi in [(401,800),(801,1600),(1601,3200),(3201,6400),(6401,12800),(12801,25600),(25601,50000)]:
    s = window_stats(lo, hi); xs.append(math.log(math.sqrt(lo*hi))); ys.append(math.log(s['frac']))
n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
slope = sum((x-mx)*(y-my) for x, y in zip(xs, ys)) / sum((x-mx)**2 for x in xs)
print(f"\norphan-fraction exponent (least squares over dyadic windows 400..50000): alpha = {-slope:.3f}")
# per-pair local exponents
for i in range(n-1):
    print(f"  local exponent between windows {i},{i+1}: {-(ys[i+1]-ys[i])/(xs[i+1]-xs[i]):.3f}")

# interval structure
holes = Counter(); interval = 0; spans = []
depth = [0]*(KMAX+2)
for C in clusters.values():
    sp = max(C) - min(C) + 1; h = sp - len(C); holes[h] += 1; spans.append(sp)
    if h == 0: interval += 1
    for k in range(min(C), max(C)+1): depth[k] += 1
print(f"\n== interval structure == clusters: {len(clusters)}; intervals (no holes): {interval} ({interval/len(clusters):.3f}); holes histogram (first 10): {sorted(holes.items())[:10]}")
print("largest clusters (size, min, max, holes):", sorted(((len(C), min(C), max(C), max(C)-min(C)+1-len(C)) for C in clusters.values()), reverse=True)[:8])
print("max interleaving depth (number of clusters whose span covers K) by window:")
for lo, hi in [(3,1000),(1000,5000),(5000,12800),(12800,25000),(25000,50000)]:
    d = depth[lo:hi+1]; print(f"  {lo}-{hi}: max {max(d)}, mean {sum(d)/len(d):.2f}")

# distinct odd counts
print("\n== distinct odd-step counts sigma(M_K) for K <= N (THM-4556 vi: 23,37,58,104 at 100,400,1200,6000) ==")
for N in (100, 400, 1200, 6000, 12800, 25000, 50000):
    print(f"  N={N}: {len({rows[K][0] for K in range(2, N+1)})}")

# last odd value before 1
lj = Counter(rows[K][2] for K in rows if K >= 3)
print("\n== last odd value before 1 is (2^j-1)/3: histogram of j (K>=3) ==", sorted(lj.items()))
print("   j=4 <-> 5, j=8 <-> 85, j=10 <-> 341, j=14 <-> 5461, j=16 <-> 21845 ... (j = 6, 12, 18: (4^m-1)/3 divisible by 3, unreachable)")

# payment law
r = [(rows[K][4] - K) / K for K in range(100, KMAX+1)]
print(f"\n== payment (first Terras time below the source): (pay-K)/K over K>=100: mean {statistics.mean(r):.4f}, sd {statistics.pstdev(r):.4f}, min {min(r):.3f}, max {max(r):.3f}; prediction 2/(1-log2(3/2)) - 1 ... drift: 1/(1-log2(3)/2) = {1/(1-math.log2(3)/2):.4f}; expected (pay-K)/K = (log2 3 - 1)/(1 - log2(3)/2) = {(math.log2(3)-1)/(1-math.log2(3)/2):.4f}")
pk = [(rows[K][3]) / K for K in range(100, KMAX+1)]
print(f"   peak bits / K: mean {statistics.mean(pk):.4f}, max {max(pk):.4f} (run end = log2(3) = {math.log2(3):.4f})")

# orphans vs v2(K-1)
def v2(x):
    return (x & -x).bit_length() - 1
print("\n== orphan rate among odd K in [1000, 50000] by v2(K-1) (two-run length J = 1 + floor(v2/2)) ==")
byv = defaultdict(lambda: [0, 0])
for K in range(1001, KMAX+1, 2):
    v = min(v2(K-1), 7); byv[v][0] += 1
    if min(cl_of[K]) == K: byv[v][1] += 1
for v in sorted(byv):
    n_, o_ = byv[v]; print(f"  v2(K-1)={v}{'+' if v==7 else ''}: n={n_}, orphans={o_}, rate={o_/n_:.4f}")

json.dump(dict(KMAX=KMAX, clusters=len(clusters), windows=W, alpha=-slope, intervals=interval,
               distinct_odd={N: len({rows[K][0] for K in range(2, N+1)}) for N in (100,400,1200,6000,12800,25000,50000)},
               last_j=sorted(lj.items())), open("analyze_horizon_summary.json", "w"), indent=1)
