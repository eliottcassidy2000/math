# Two-dimensional deletion-pair census: exact window statements for bases of order <= 10.
import sys
HM, KM = 1024, 4096
def tab(fn):
    b = open(fn,'rb').read(); t = {}
    for h in range(1, HM, 2):
        row = b[(h>>1)*(KM//8):(h>>1)*(KM//8)+KM//8]
        if any(row): t[h] = int.from_bytes(row,'little')
    return t
A = tab(sys.argv[1]); Z = tab(sys.argv[2])
allodd = [h for h in range(1, HM, 2)]
print("realised h < 1024:", len(A), " unrealised odd h:", [h for h in allodd if h not in A])
def missing(m, h, lo, hi):   # odd n in [lo*h, hi*h) not of the form h + 2k with bit k set
    out = []
    n = lo*h if (lo*h) % 2 == 1 else lo*h + 1
    while n < hi*h:
        k = (n - h)//2
        if k < KM and not (m >> k & 1): out.append(n)
        n += 2
    return out
exc35 = {h: missing(m, h, 3, 5) for h, m in A.items()}
exc35 = {h: v for h, v in exc35.items() if v}
print("\n(W1) window [3h,5h): complete for every realised h except", len(exc35), "values; largest exception h =", max(exc35))
for h in sorted(exc35): print("   h=%4d missing n: %s" % (h, exc35[h]))
exc39 = {h: missing(m, h, 3, 9) for h, m in A.items() if 9*h - h < 2*KM}
exc39 = {h: v for h, v in exc39.items() if v}
print("\n(W2) window [3h,9h): complete for every realised h (with 9h in table range, h <= %d) except %d values; largest exception h = %d" % (max(h for h in A if 9*h-h < 2*KM), len(exc39), max(exc39)))
# where does fullness start?  n0(h) = least odd n0 such that every odd n in [n0, 5h) is an extension value
import statistics
ratios = []
for h, m in A.items():
    if h < 237: continue
    n = 5*h - 1 if (5*h - 1) % 2 == 1 else 5*h - 2
    while n > h and (m >> ((n-h)//2) & 1): n -= 2
    ratios.append(((n+2)/h, h))
print("\nfull-from ratio n0(h)/h for 237 <= h <= 1023: min %.3f  median %.3f  max %.3f (at h=%d)" % (min(r for r,_ in ratios), statistics.median(r for r,_ in ratios), max(ratios)[0], max(ratios)[1]))
# least positive increment
print("least extension ratio n_min(h)/h (over h>=237): min %.4f" % min((h + 2*next(k for k in range(1, KM) if A[h] >> k & 1))/h for h in A if h >= 237))
# out-degree-one extensions only
def dens(T, lo, hi, hmin):
    tot = pres = 0
    for h, m in T.items():
        if h >= hmin:
            n = -int(-lo*h // 1); n += (n % 2 == 0)      # least odd integer >= lo*h
            while n < hi*h:
                k = (n-h)//2
                if k < KM: tot += 1; pres += m >> k & 1
                n += 2
    return pres/tot
print("\nout-degree-one extensions only (h >= 237): density in [3h,5h): %.4f   in [2h,3h): %.4f" % (dens(Z, 3, 5, 237), dens(Z, 2, 3, 237)))
print("all extensions (h >= 237): density in [2h,3h): %.5f   in [3h/2,2h): %.4f" % (dens(A, 2, 3, 237), dens(A, 1.5, 2, 237)))
