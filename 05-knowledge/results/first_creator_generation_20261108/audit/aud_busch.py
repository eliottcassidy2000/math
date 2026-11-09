"""aud_busch.py - auditor's recomputation of PROOF.md section 4.2 (values born at their Busch bound) and of the
S(n) agreement with klein's least_order_11.txt.  Reads klein's table and the auditor's own aud_S<n>.txt (aud_gen.c).
usage: python aud_busch.py KLEIN_least_order_11.txt"""
import sys
def f(n):                        # Moon-Busch floor, n >= 3: min 3^a 5^b with 2a + 3b = n - 1
    best = None
    for b in range(0, n):
        r = n - 1 - 3 * b
        if r < 0: break
        if r % 2 == 0:
            v = 3 ** (r // 2) * 5 ** b
            best = v if best is None or v < best else best
    return best
def B(h):
    n = 3
    while f(n + 1) <= h: n += 1
    return n if f(3) <= h else None
assert [f(n) for n in range(3, 14)] == [3, 5, 9, 15, 25, 45, 75, 125, 225, 375, 625], [f(n) for n in range(3, 14)]
print("f(3..26):", [f(n) for n in range(3, 27)])
nh, unatt = {}, []
for line in open(sys.argv[1]):
    if line.startswith('#'): continue
    p = line.split()
    h = int(p[0])
    if p[1].startswith('>'): unatt.append(h)
    else: nh[h] = int(p[1])
odd_all = list(range(1, 95096, 2))
print("table rows:", len(nh) + len(unatt), " odd values <= 95095:", len(odd_all), " numeric n(h):", len(nh), " '>11':", len(unatt))
born = [h for h in sorted(nh) if h >= 3 and nh[h] == B(h)]
print("n(h) = B(h) among h with n(h) <= 11:", born)
un = [h for h in unatt if h not in (7, 21)]
print("unattained at order 11 (besides 7, 21):", len(un), " min", min(un), " max", max(un), " B(h) range:", sorted(set(B(h) for h in un)))
print("largest h with every odd value in [1, h] \ {7, 21} attained at order <= 11:", min(un) - 2)
for n in range(1, 12):
    S = {h for h in nh if nh[h] <= n}
    try:
        mine = {int(x) for x in open('aud_S%d.txt' % n).read().split()}
        print("S(%d): klein %d values (max %d); auditor %s" % (n, len(S), max(S), "IDENTICAL" if mine == S else "DIFFERENT (%d)" % len(mine ^ S)))
    except FileNotFoundError:
        print("S(%d): klein %d values (max %d)" % (n, len(S), max(S)))
for h in (11475, 14787, 15213, 15273, 15561, 15637, 15703, 385, 2901, 63, 225):
    print("  h=%d: n(h)=%s B(h)=%s" % (h, nh.get(h, '>11'), B(h)))
