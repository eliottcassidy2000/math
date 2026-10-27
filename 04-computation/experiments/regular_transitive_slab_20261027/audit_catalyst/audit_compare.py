"""Compare the audit's own computations with THM-4820 and the producer's output files.

Inputs (audit, produced by audit_lattice / audit_dp in this directory):
  audit_lat_W2500.txt (mode 6, closure to 2600), audit_lat_W2500_m3.txt (three-unit seeds), audit_lat_W3500.txt,
  audit_lat_W7500.txt (closure to 7500), audit_lat_W400.txt, audit_lat_W900.txt, audit_dp_76_62_62.out
Inputs (producer, read only for comparison): ../catalyst/lattice_W2500.txt, ../catalyst/lattice_W7500.txt,
  ../catalyst/merge_debt_11.out
"""
import re
import sys
from collections import defaultdict

FAIL = []


def check(c, m):
    if not c:
        FAIL.append(m)
        print("FAIL:", m)


def Q(a, b):
    return a * a - a * b + b * b


PROD = "../catalyst/"


def load_audit(fn):
    pts, rat = {}, {}
    for line in open(fn):
        t = line.split()
        if t[0] == "P":
            a, b, q, g, ing = map(int, t[1:])
            pts[(a, b)] = (q, g, ing)
        elif t[0] == "N":
            rat[int(t[1])] = int(t[2])
    return pts, rat


def load_prod_lattice(fn):
    pts, rat = {}, {}
    for line in open(fn):
        if line.startswith("P "):
            a, b, q, g = map(int, line.split()[1:])
            pts[(a, b)] = (q, g)
        m = re.match(r"g\((\d+)\) = (\d+)", line)
        if m:
            rat[int(m.group(1))] = int(m.group(2))
    return pts, rat


A25, R25 = load_audit("audit_lat_W2500.txt")
A25m3, _ = load_audit("audit_lat_W2500_m3.txt")
A35, R35 = load_audit("audit_lat_W3500.txt")
A75, R75 = load_audit("audit_lat_W7500.txt")
_, R4 = load_audit("audit_lat_W400.txt")
_, R9 = load_audit("audit_lat_W900.txt")
P25, PR25 = load_prod_lattice(PROD + "lattice_W2500.txt")
_, PR75 = load_prod_lattice(PROD + "lattice_W7500.txt")

print("=== (D) lattice table ===")
# 1. point-by-point agreement with the producer's W = 2500 table
w25 = {x: v for x, v in A25.items() if v[0] <= 2500}
check(set(w25) == set(P25), "same point set (norm <= 2500)")
mism = [x for x in P25 if P25[x][1] != w25[x][1]]
mism400 = [x for x in mism if P25[x][0] <= 400]
n400 = sum(1 for x in P25 if P25[x][0] <= 400)
print(f"producer lattice_W2500.txt vs audit Knuth-Dijkstra W=2500: {len(P25)} points, {len(mism)} mismatches "
      f"({n400} points of norm <= 400, {len(mism400)} mismatches there)")
check(not mism, "lattice mismatch")
# 2. three-unit seeds give the same table (A4 at the lattice level)
m3 = [x for x in w25 if A25m3[x][1] != w25[x][1]]
print(f"three-unit seeds (1, w, w^2 only) vs six-unit seeds at W=2500: {len(m3)} differences")
check(not m3, "A4 lattice")
# 3. window stability W=2500 vs W=3500 on norms <= 2500
st = [x for x in w25 if A35[x][1] != w25[x][1]]
print(f"W=2500 vs W=3500 on norms <= 2500: {len(st)} differences")
check(not st, "window stability")

# 4. two values only on norms <= 2000; exceptional norms
exc_pts = defaultdict(list)
all_pts_by_norm = defaultdict(list)
excess_hist = defaultdict(int)
for x, (q, g, ing) in w25.items():
    if q > 2000:
        continue
    e = g - q // 3
    excess_hist[e] += 1
    all_pts_by_norm[q].append(x)
    if e >= 1:
        exc_pts[q].append(x)
print(f"norms <= 2000: excess g - floor(Q/3) histogram {dict(excess_hist)}")
check(set(excess_hist) <= {0, 1}, "two values only")
exc_norms = sorted(exc_pts)
claimed = [121, 229, 313, 421, 484, 529, 625, 733, 841, 853, 1033, 1093, 1213, 1429, 1597, 1609, 1681, 1693, 1825, 1936]
print(f"exceptional norms <= 2000 ({len(exc_norms)}): {exc_norms}")
check(exc_norms == claimed, "exceptional norm list")
print(f"    all exceptional at Q = 1 mod 3: {all(q % 3 == 1 for q in exc_norms)}")
n1 = sorted(q for q in all_pts_by_norm if q % 3 == 1)
print(f"    represented norms = 1 mod 3 up to 2000: {len(n1)} including the unit norm 1, {len(n1) - 1} excluding it "
      f"(theorem: 328, i.e. norm 1 excluded as in tight_splits.py)")
check(len(n1) - 1 == 328 and n1[0] == 1, "328 norms (excluding norm 1)")
mixed = [q for q in exc_norms if len(exc_pts[q]) != len(all_pts_by_norm[q])]
print(f"    exceptional norms where NOT every point of that norm is exceptional: {mixed}")
for q in exc_norms:
    print(f"      norm {q}: {len(exc_pts[q])}/{len(all_pts_by_norm[q])} points exceptional")

# 5. window-free closure vs window table
dis = [(x, v) for x, v in A25.items() if v[2] != -1 and v[0] <= 2500 and ((v[1] - v[0] // 3 >= 1) != (v[2] == 0))]
nexact = sum(1 for x, v in A25.items() if v[2] != -1 and v[0] <= 2500)
print(f"window-free tight closure (all points to norm 2600) vs W=2500 exceptional set: {nexact} points compared, "
      f"{len(dis)} disagreements")
check(not dis, "closure vs window")
dis75 = [(x, v) for x, v in A75.items() if v[2] != -1 and ((v[1] - v[0] // 3 >= 1) != (v[2] == 0))]
n75 = sum(1 for x, v in A75.items() if v[2] != -1)
print(f"window-free tight closure (to norm 7500) vs W=7500 exceptional set: {n75} points compared, {len(dis75)} disagreements")
check(not dis75, "closure vs window 7500")

# 6. W = 7500: excess histogram on all norms <= 7499, first norm = 0 mod 3 exception
h75 = defaultdict(int)
exc75 = defaultdict(list)
for x, (q, g, ing) in A75.items():
    if q >= 7500:
        continue
    e = g - q // 3
    h75[e] += 1
    if e >= 1:
        exc75[q].append(x)
print(f"W=7500, all points of norm <= 7499: excess histogram {dict(h75)}")
z0 = sorted(q for q in exc75 if q % 3 == 0)
print(f"    exceptional norms = 0 mod 3 below 7500: {z0[:12]}{' ...' if len(z0) > 12 else ''} (count {len(z0)})")
print(f"    number of exceptional norms below 7500: {len(exc75)}; total exceptional points {sum(len(v) for v in exc75.values())}")
ndist = len({q for (x, (q, g, ing)) in A75.items() if q < 7500})
print(f"    represented norms below 7500: {ndist}")

# 7. rational integers
print("=== rational n ===")
for n in range(1, 17):
    check(R4[n] == R9[n] == R25[n] == PR25[n], f"W=400/900/2500 agreement n={n}")
print("g(n), n=1..16, audit W=400, 900, 2500 and producer W=2500 agree:", [R25[n] for n in range(1, 17)])
check([R25[n] for n in range(1, 17)] == [0, 1, 3, 5, 8, 12, 16, 21, 27, 33, 41, 48, 56, 65, 75, 85], "g sequence")
d75 = [n for n in PR75 if PR75[n] != R75[n]]
print(f"producer lattice_W7500.txt vs audit W=7500, n <= 86: {len(d75)} mismatches")
check(not d75, "W7500 rational")
check(all(R75[n] == R25[n] for n in range(1, 51)), "W7500 vs W2500 n<=50")
exc_rat = [n for n in range(1, 87) if R75[n] - (n * n) // 3 >= 1]
print(f"exceptional rational n <= 86 (W=7500): {exc_rat}")
check([n for n in exc_rat if n <= 50] == [11, 22, 23, 25, 29, 41, 44, 46, 47, 50], "rational n<=50 list")
check([n for n in exc_rat if 50 < n] == [51, 58, 59, 82, 83, 85], "rational 51..86 list")
print(f"    max excess on n <= 86: {max(R75[n] - n * n // 3 for n in range(1, 87))}")
inG = {n: A75[(n, 0)][2] for n in range(1, 87)}
print(f"    window-free closure status available (B(n) <= 7500) for all n <= 86: {all(v != -1 for v in inG.values())}; "
      f"closure-exceptional n: {[n for n in range(1, 87) if inG[n] == 0]}")

# 8. tight splits of the exceptional points: every tight split has a non-closure part (propagation)
Gset = {x for x, v in A75.items() if v[2] == 1}
Gknown = {x for x, v in A75.items() if v[2] != -1}


def tight_splits(x):
    qx = Q(*x)
    L = int((4 * (qx + 2) / 3) ** 0.5) + 2
    out = []
    for a in range(-L, L + 1):
        for b in range(-L, L + 1):
            y = (a, b)
            z = (x[0] - a, x[1] - b)
            if y == (0, 0) or z == (0, 0):
                continue
            qy, qz = Q(*y), Q(*z)
            bt = qx - qy - qz
            if max(bt, 0) + qy // 3 + qz // 3 - qx // 3 == 0:
                out.append((y, z, bt))
    return out


print("=== propagation structure of the exceptional points (norms <= 2000, one representative per norm) ===")
for q in exc_norms:
    x = exc_pts[q][0]
    ts = tight_splits(x)
    both_good = [s for s in ts if s[0] in Gset and s[1] in Gset]
    check(not both_good, f"exceptional {x} has a tight split with both parts in G")
    print(f"  norm {q}, x={x}: {len(ts)} ordered tight splits, none with both parts in G")
x229 = None
for x in exc_pts[229]:
    for (y, z, bt) in tight_splits(x):
        if y == (-11, 0) and z == (-6, -12):
            x229 = (x, bt)
print("  theorem's example 229 = (-11) + (-6 - 12w):", x229)
for n in (51, 59):
    ts = tight_splits((n, 0))
    print(f"  n={n}: {len(ts)} ordered tight splits, e.g. {ts[:2]}; any with both parts in G: "
          f"{any(s[0] in Gset and s[1] in Gset for s in ts)}")

# 9. rows of f and thresholds from the audit DP
print("=== (C)/(E) f(n+k,k,k) from audit_dp 76 62 62 ===")
F = {}
for line in open("audit_dp_76_62_62.out"):
    if line.startswith("F "):
        n, k, v = map(int, line.split()[1:])
        F[(n, k)] = v
row = [F[(11, k)] for k in range(0, 66) if (11, k) in F]
print("f(11+k,k,k), k = 0..%d:" % (len(row) - 1), row)
thm_row = list(map(int, "110 109 106 103 100 97 94 92 89 87 84 82 80 78 75 73 71 69 67 65 63 61 59 57 55 53 52 50 49 48 47 46 45 44 43 42 41".split()))
check(row[:37] == thm_row, "theorem row k <= 36")
check(all(v == 41 for v in row[36:]), "41 for k >= 36 in the box")
prod = {}
for line in open(PROD + "merge_debt_11.out"):
    m = re.match(r"k=(\d+) f\(\d+,\d+,\d+\)=(\d+)", line)
    if m:
        prod[int(m.group(1))] = int(m.group(2))
check(all(prod[k] == row[k] for k in prod), "producer merge_debt_11.out row")
print(f"   agrees with the theorem's printed row (k <= 36) and with merge_debt_11.out (k <= 60); "
      f"constant 41 for 36 <= k <= {len(row) - 1}; first k with 41: {row.index(41)}")
kstar = {}
for n in range(1, 15):
    g = R25[n]
    ks = [k for k in range(0, 63) if (n, k) in F and F[(n, k)] == g]
    kmax = max(k for k in range(0, 63) if (n, k) in F)
    check(all(F[(n, k)] >= g for k in range(0, kmax + 1)), f"f >= g for n={n}")
    check(all(F[(n, k)] >= F[(n, k + 1)] for k in range(0, kmax)), f"non-increasing n={n}")
    kstar[n] = ks[0] if ks else None
    lb = (n * n - n - g + 2) // 3  # ceil((n^2 - n - g)/3)
    check(kstar[n] is not None and kstar[n] >= lb, f"lower bound n={n}")
print("k*(n), n=1..14:", [kstar[n] for n in range(1, 15)])
check([kstar[n] for n in range(1, 15)] == [0, 1, 2, 4, 7, 10, 14, 19, 24, 31, 36, 44, 52, 61], "thresholds")
print("ceil lower bounds (n^2-n-g(n))/3:", [-(-(n * n - n - R25[n]) // 3) for n in range(1, 15)])
print("f(n+K,K,K) at K = 4n+8 for n = 15, 16 (theorem: K insufficient):",
      {n: (F.get((n, 4 * n + 8)), R25[n]) for n in (15, 16)})
print("f(15+k,k,k), k = 55..61:", [F[(15, k)] for k in range(55, 62)], " g(15) =", R25[15])

print()
print("TOTAL FAILURES:", len(FAIL))
sys.exit(1 if FAIL else 0)
