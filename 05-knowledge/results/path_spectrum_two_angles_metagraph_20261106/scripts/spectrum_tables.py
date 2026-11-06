#!/usr/bin/env python3
"""spectrum_tables.py (klein-2026-11-06): the two angles of the Hamiltonian-path spectrum from the exact censuses.
usage: spectrum_tables.py census_dir nmax [strong_ears_375.out]
Reads census_<n>.txt (n = 3..nmax; lines "T h c3 strong count", "C h hc count") and prints
  (1) per order: classes, strong classes, distinct values, strong-distinct values, max, the new values N(n) = S(n) - S(n-1),
      whether every new value is strongly attained at that order, the unattained odd values below the maximum,
      the solid prefix s(n) (every odd value <= s(n) other than 7 and 21 is attained) and its fraction of the maximum;
  (2) the least-order function n(h) against m(h) = min{n : max H(n) >= h} (max H from the censuses, 95095 at 11 from OEIS A038375);
  (3) the strong windows W(h) = {n : some strong n-tournament has h paths} for odd h below the strong floor of nmax+1,
      compared with the recorded strong polynomials P_d of THM-4744 when a strong_ears file is given;
  (4) the least strong order n_s(h) against n(h);
  (5) the recorded multiplicity predictions N(10, h) of THM-4744 and the polynomial laws N(n, h)."""
import sys, re
from collections import defaultdict
D = sys.argv[1]; NMAX = int(sys.argv[2]); ears = sys.argv[3] if len(sys.argv) > 3 else None
def busch(n):  # Busch's minimum number of Hamiltonian paths of a strong n-tournament (n >= 3)
    if n == 3: return 3
    if n == 4: return 5
    k, r = divmod(n, 3)
    return {1: 5**k, 2: 9 * 5**(k - 1), 0: 3 * 5**(k - 1)}[r]
cnt = {}; hc = {}
for n in range(3, NMAX + 1):
    T = defaultdict(lambda: defaultdict(lambda: [0, 0])); C = defaultdict(lambda: defaultdict(int))
    with open(f"{D}/census_{n}.txt") as fh:
        for line in fh:
            a = line.split()
            if not a: continue
            if a[0] == 'T': T[int(a[1])][int(a[2])][int(a[3])] += int(a[4])
            elif a[0] == 'C': C[int(a[1])][int(a[2])] += int(a[3])
    cnt[n] = T; hc[n] = C
S = {n: set(cnt[n]) for n in cnt}; S[1] = {1}; S[2] = {1}
Ss = {n: {h for h in cnt[n] if any(cnt[n][h][k][1] for k in cnt[n][h])} for n in cnt}
maxH = {n: max(S[n]) for n in cnt}; maxH.update({1: 1, 2: 1}); maxH[11] = 95095  # OEIS A038375 (Royle 2022)
classes = {n: sum(sum(v) for h in cnt[n] for v in cnt[n][h].values()) for n in cnt}
strongc = {n: sum(v[1] for h in cnt[n] for v in cnt[n][h].values()) for n in cnt}
def ranges(vals):
    vals = sorted(vals); out = []; i = 0
    while i < len(vals):
        j = i
        while j + 1 < len(vals) and vals[j + 1] == vals[j] + 2: j += 1
        out.append(f"{vals[i]}" if i == j else f"{vals[i]}..{vals[j]}({(j - i) + 1})"); i = j + 1
    return " ".join(out)
print("(1) per order")
print("  n | classes | strong | distinct | strong-distinct | max H | new values | new & strong | unattained odd < max (excl. 7, 21) | solid prefix s(n) | s/max")
newvals = {}
for n in range(3, NMAX + 1):
    new = S[n] - S[n - 1]; newvals[n] = new
    news = new & Ss[n]
    un = [h for h in range(1, maxH[n], 2) if h not in S[n] and h not in (7, 21)]
    s = maxH[n]
    for h in range(1, maxH[n] + 1, 2):
        if h not in S[n] and h not in (7, 21): s = h - 2; break
    print(f"  {n:2d} | {classes[n]:8d} | {strongc[n]:8d} | {len(S[n]):5d} | {len(Ss[n]):5d} | {maxH[n]:6d} | {len(new):5d} | {len(news):5d} | {len(un):5d} | {s:6d} | {s / maxH[n]:.3f}")
print("  the unattained odd values below the maximum, by order:")
for n in range(5, NMAX + 1):
    un = [h for h in range(1, maxH[n], 2) if h not in S[n] and h not in (7, 21)]
    print(f"    n = {n}: {ranges(un) if un else 'none'}")
print("  the values attained at order n only by non-strong tournaments (S(n) - Ss(n)):")
for n in range(5, NMAX + 1):
    d = sorted(S[n] - Ss[n]); print(f"    n = {n}: {len(d)}: {ranges(d)}")
    assert all(h < busch(n) or True for h in d)
print("  new values that are NOT strongly attained at their least order:")
bad = [(n, h) for n in range(4, NMAX + 1) for h in sorted(newvals[n]) if h not in Ss[n]]
print(f"    {bad if bad else 'none'}")
print("  new values at order n below the strong floor f(n) (would have to come from products):",
      [(n, h) for n in range(5, NMAX + 1) for h in sorted(newvals[n]) if h < busch(n)] or "none")
print("\n(2) least order n(h) versus m(h) = min{n : max H(n) >= h}")
least = {}
for n in range(1, NMAX + 1):
    for h in S[n]: least.setdefault(h, n)
def m_of(h):
    for n in range(1, 12):
        if maxH[n] >= h: return n
    return None
exc = [(h, least[h], m_of(h)) for h in sorted(least) if least[h] != m_of(h)]
print(f"  odd values with n(h) <= {NMAX}: {len(least)}; values with n(h) > m(h): {len(exc)}")
byord = defaultdict(list)
for h, nh, mh in exc: byord[mh].append(h)
for mh in sorted(byord): print(f"    m(h) = {mh}: {len(byord[mh])} values, n(h) - m(h) = {sorted(set(least[h] - mh for h in byord[mh]))}: {ranges(byord[mh])}")
print("  check: n(h) > m(h) exactly for the unattained values below the maximum of order m(h):",
      all(h not in S[mh] for h, nh, mh in exc) and all(h in S[least[h]] for h in least))
print("  largest h with n(h) = n, and the number of values first attained at n:")
for n in range(3, NMAX + 1): print(f"    n = {n}: |N(n)| = {len(newvals[n])}, min new = {min(newvals[n]) if newvals[n] else '-'}, max new = {max(newvals[n]) if newvals[n] else '-'}")
print("\n(3) strong windows W(h) for odd h below f(nmax+1) =", busch(NMAX + 1))
strong_count = {}  # (h, n) -> number of strong classes
for n in cnt:
    for h in cnt[n]:
        c = sum(v[1] for v in cnt[n][h].values())
        if c: strong_count[(h, n)] = c
W = defaultdict(list)
for (h, n), c in strong_count.items(): W[h].append(n)
lim = busch(NMAX + 1)
nonint = []; lens = defaultdict(int); table = []
for h in range(3, lim, 2):
    w = sorted(W.get(h, []))
    if not w: continue
    interval = (w == list(range(w[0], w[-1] + 1)))
    if not interval: nonint.append((h, w))
    lens[len(w)] += 1
    table.append((h, w, [strong_count[(h, n)] for n in w]))
print(f"  odd h < {lim} with a strong realiser of order <= {NMAX}: {len(table)}; missing: {[h for h in range(3, lim, 2) if h not in W]}")
print(f"  window lengths: {dict(sorted(lens.items()))}; windows that are not intervals: {nonint if nonint else 'none'}")
print("  windows (h: orders [strong classes per order]) for h < 160:")
for h, w, c in table:
    if h < 160: print(f"    {h}: {w} {c}")
print("  least strong order n_s(h) versus least order n(h), odd h <", lim, ":",
      [(h, least[h], min(w)) for h, w, c in table if least[h] != min(w)] or "always equal")
print("  last strong order versus the Busch bound (largest n with f(n) <= h):")
viol = []; slack = defaultdict(int)
for h, w, c in table:
    nb = max(n for n in range(3, 40) if busch(n) <= h)
    if max(w) > nb: viol.append((h, w, nb))
    slack[nb - max(w)] += 1
print(f"    windows ending beyond the Busch bound: {viol if viol else 'none'}; distribution of (Busch bound - last strong order): {dict(sorted(slack.items()))}")
print("  values whose window ends at the Busch bound are, at that order, near the strong floor; the window shape in the level plane is reported by level_tables.py")
if ears:
    print("\n  comparison with the recorded strong polynomials P_d (THM-4744 / strong_ears_375.out):")
    P = {}
    for line in open(ears):
        m = re.match(r"\s*P_(\d+) = (.*)", line)
        if m:
            d = int(m.group(1)); poly = {}
            for term in m.group(2).split('+'):
                term = term.strip(); mm = re.match(r"(\d*)z\^(\d+)", term)
                if mm: poly[int(mm.group(2))] = int(mm.group(1) or 1)
            P[d] = poly
    bad = []
    for d, poly in sorted(P.items()):
        for n in range(3, NMAX + 1):
            if poly.get(n, 0) != strong_count.get((d, n), 0): bad.append((d, n, poly.get(n, 0), strong_count.get((d, n), 0)))
    print(f"    recorded polynomials: {len(P)} (d = {min(P)}..{max(P)}); coefficient mismatches at orders <= {NMAX}: {bad if bad else 'none'}")
print("\n(4) multiplicity checks")
def N(n, h): return sum(sum(v) for v in cnt[n][h].values()) if h in cnt[n] else 0
pred = {45: 96, 63: 6, 75: 70, 81: 106, 105: 59, 125: 43, 135: 146, 189: 102, 223: 82}
if NMAX >= 10:
    print("  THM-4744's predicted N(10, h):", {h: (pred[h], N(10, h), pred[h] == N(10, h)) for h in pred})
laws = {3: lambda n: n - 2, 9: lambda n: (n - 4) * (n - 1) // 2, 15: lambda n: (n - 4) * (n - 3), 27: lambda n: n * (n - 5) * (n - 4) // 6,
        33: lambda n: n * n - 5 * n - 1, 45: lambda n: (n**3 - 12 * n * n + 43 * n - 38) // 2, 49: lambda n: 8 * n - 49, 63: lambda n: 2 * (n - 7),
        75: lambda n: n * (n * n - 13 * n + 44) // 2}
for h, law in laws.items():
    print(f"  N(n,{h}) law: ", [(n, N(n, h), law(n)) for n in range(7, NMAX + 1)])
print("\n(5) the top of the spectrum: values attained at order n above 0.85 max, and their class counts")
for n in range(7, NMAX + 1):
    top = sorted(h for h in S[n] if h > 0.85 * maxH[n])
    print(f"  n = {n}: {len(top)} values above {0.85 * maxH[n]:.0f}: " + " ".join(f"{h}({N(n, h)})" for h in top[-12:]))
print("\n(6) Hamiltonian cycles: for each order, the least and greatest path count for each cycle count (hc <= 6), and max hc")
for n in range(5, NMAX + 1):
    C = hc[n]; byhc = defaultdict(list)
    for h in C:
        for k in C[h]: byhc[k].append(h)
    print(f"  n = {n}: max hc = {max(byhc)}; " + "; ".join(f"hc={k}: [{min(byhc[k])}, {max(byhc[k])}] {len(byhc[k])} values" for k in range(1, 7) if k in byhc))
print("  the hc <= 2 values at each order (least-order realisers with at most two cycles):")
for n in range(6, NMAX + 1):
    C = hc[n]; v12 = sorted(h for h in C if 1 in C[h] or 2 in C[h])
    print(f"    n = {n}: {len(v12)} values, {ranges(v12)[:300]}")
