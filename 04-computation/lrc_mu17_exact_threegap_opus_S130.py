from fractions import Fraction as F

def mu17_exact(k):
    """EXACT meas{ x in (0,1) : circular max-gap of {frac(jx): j=1..k} > 1/7 }.
    max-gap(x) is piecewise-linear; breakpoints at Farey fractions m/d, d<=k
    (phase coincidences |i-j|<=k-1 and wraps i<=k). On each order-interval every
    gap is linear, so {some gap > 1/7} is a union of rational sub-intervals."""
    thr = F(1,7)
    # collect breakpoints
    bps = {F(0), F(1)}
    for d in range(1, k+1):
        for m in range(1, d):
            bps.add(F(m, d))
    bps = sorted(bps)
    total = F(0)
    for a, b in zip(bps, bps[1:]):
        mid = (a + b) / 2
        # phases at mid: p_j = j*mid - floor(j*mid)
        fl = {}
        ph = []
        for j in range(1, k+1):
            fj = (j*mid).__floor__()
            fl[j] = fj
            ph.append((j*mid - fj, j))
        ph.sort()  # sorted by phase value at mid
        order = [j for _, j in ph]
        # linear phase p_j(x) = j*x - fl[j]; gaps between consecutive sorted + wrap
        seglen = F(0)
        # each gap g(x) = c*x + b0 ; find sub-interval of (a,b) where g > 1/7, union them
        subs = []  # list of (lo,hi) where some gap>thr
        gaps = []
        for s in range(k):
            j1 = order[s]; j2 = order[(s+1) % k]
            if s < k-1:
                c = F(j2 - j1); b0 = F(-(fl[j2] - fl[j1]))
            else:  # wrap: p_first + 1 - p_last
                c = F(j1_first := order[0]) - F(order[-1]); b0 = F(-(fl[order[0]] - fl[order[-1]])) + 1
                c = F(order[0] - order[-1]); b0 = F(-(fl[order[0]] - fl[order[-1]]) + 1)
            gaps.append((c, b0))
        # region where max gap > thr = union over gaps of {c*x+b0>thr} ∩ (a,b)
        ivs = []
        for c, b0 in gaps:
            # c*x + b0 > thr
            if c == 0:
                if b0 > thr: ivs.append((a, b))
            elif c > 0:
                x0 = (thr - b0)/c
                lo = max(a, x0)
                if lo < b: ivs.append((lo, b))
            else:
                x0 = (thr - b0)/c
                hi = min(b, x0)
                if a < hi: ivs.append((a, hi))
        # union length
        ivs.sort()
        cur = None
        for lo,hi in ivs:
            if lo>=hi: continue
            if cur is None: cur=[lo,hi]
            elif lo<=cur[1]: cur[1]=max(cur[1],hi)
            else: seglen += cur[1]-cur[0]; cur=[lo,hi]
        if cur is not None: seglen += cur[1]-cur[0]
        total += seglen
    return total

print("=== EXACT mu_17(AP {1..k}) via three-gap piecewise-linear breakpoints (opus-S130) ===")
claimed = {13: F(477,1078)}
for k in range(3,14):
    m = mu17_exact(k)
    tag = ""
    if k in claimed: tag = f"   claimed rhoGlobFloorRat={claimed[k]}={float(claimed[k]):.5f}  MATCH={m==claimed[k]}"
    note = "  (>=1/7 always: k<=7 pigeonhole => mu=1)" if k<=7 else ""
    print(f"  k={k:>2}: mu_17 = {m} = {float(m):.6f}{note}{tag}")
print("\n  => exact floor values for k=8..13 (the AP minimizer). These are the rigorous witness-floor")
print("     constants; the density floor mu_17(E) >= mu_17(AP) > 0 reduces to (i) AP-minimality + (ii) these.")
