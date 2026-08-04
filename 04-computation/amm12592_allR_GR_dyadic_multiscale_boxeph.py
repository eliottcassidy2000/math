#!/usr/bin/env python3
"""ANGLE B3 step 5 -- dyadic multiscale (product-telescope) expansion of G_R.

Substitute u = -x (so q = 1+u).  With R = 2^t, R-1 = 2^t - 1 = sum 2^j:
  (1+u)^{R-1} = prod_{j<t} (1+u)^{2^j},    V := prod_{j<t} (1+u^{2^j})
             = 1 + u + ... + u^{R-1}   (binary digits),
and the telescope  prod a - prod b = sum_l (prod_{j<l} a)(a_l - b_l)(prod_{j>l} b)
with a_l - b_l = (1+u)^{2^l} - 1 - u^{2^l} = 2 W_l,
  W_l := sum_{0<n<2^l} (C(2^l,n)/2) u^n   (integer: interior dyadic Pascal rows
  are even -- this IS the dyadic parity miracle in product form) gives

  (IDENTITY M)   Ghat_R = 1 + sum_{j odd, j<R} u^j
        + sum_{l=0}^{t-1} (1+u)^{2^l - 1} W_l  sum_{m < 2^{t-l-1}} u^{m 2^{l+1}}

where Ghat_R(u) = G_R(-u), i.e. [u^n] Ghat = (C(R-1,n) - (-1)^n)/2 >= 0.
ALL COEFFICIENTS on the right are nonnegative: this is a POSITIVE multiscale
decomposition: scale-l pieces P_l := (1+u)^{2^l-1} W_l (support [1, 2^{l+1}-2],
magnitude ~ 2^{2^{l+1}}) planted at every 2^{l+1}-aligned offset.  Each
coefficient receives O(log R) pieces (one per scale).  The top scale l = t-1 is
(1+u)^{R/2 - 1} W_{t-1} -- an epoch-R/2 binomial times the half-Pascal row:
the expansion is a doubling recursion in disguise (route (c) skeleton).

This script verifies (IDENTITY M) exactly for R = 4..1024 and tabulates piece
magnitudes vs the gamma* row-capacity envelope.
"""
from math import comb

WT = "/tmp/math-wt-boxeph-multifront"
OUT = WT + "/05-knowledge/results/amm12592_allR_GR_dyadic_multiscale_boxeph.out"
log_lines = []
def log(s=""):
    print(s); log_lines.append(s)
def flush():
    with open(OUT, "w") as f: f.write("\n".join(log_lines) + "\n")
FAIL = 0
def check(name, ok, detail=""):
    global FAIL
    if not ok: FAIL += 1
    log("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- "+detail) if detail else ""))

def pmul(a, b):
    if not a or not b: return []
    r = [0]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        if u:
            for j,v in enumerate(b): r[i+j] += u*v
    while r and r[-1] == 0: r.pop()
    return r
def padd(a, b):
    r = [0]*max(len(a), len(b))
    for i,v in enumerate(a): r[i] += v
    for i,v in enumerate(b): r[i] += v
    while r and r[-1] == 0: r.pop()
    return r
def onepu(n):  # (1+u)^n
    return [comb(n,k) for k in range(n+1)]

ok_all = True
for t in range(2, 11):
    R = 2**t
    # target Ghat
    # [u^0] Ghat = 1 (the j = 0 special case of S3: [x^0] 2G_R = 2)
    Ghat = [1] + [(comb(R-1, n) - (-1)**n)//2 for n in range(1, R)]
    assert all((comb(R-1,n) - (-1)**n) % 2 == 0 for n in range(1, R))
    # rhs
    rhs = [0]*R
    rhs[0] += 1
    for j in range(1, R, 2): rhs[j] += 1
    for l in range(t):
        Wl = [0]*(2**l)
        for n in range(1, 2**l):
            c = comb(2**l, n)
            assert c % 2 == 0
            Wl[n] = c//2
        while Wl and Wl[-1] == 0: Wl.pop()
        if not Wl: continue
        P = pmul(onepu(2**l - 1), Wl)
        for m in range(2**(t-l-1)):
            off = m * 2**(l+1)
            for s, v in enumerate(P):
                if off + s < R: rhs[off+s] += v
                elif v: raise AssertionError("overflow")
    if rhs != Ghat: ok_all = False; break
check("(IDENTITY M) dyadic multiscale expansion of Ghat_R, R = 4..1024", ok_all)
# positivity + overlap count
t = 9; R = 512
overl = [0]*R
for l in range(t):
    for m in range(2**(t-l-1)):
        off = m*2**(l+1)
        lo, hi = off+1, off + 2**(l+1) - 2
        for s in range(lo, min(hi, R-1)+1): overl[s] += 1
check("every coefficient receives <= t = log2 R pieces (R = 512)",
      max(overl) <= t, "max overlap = %d" % max(overl))
log("")
log("piece magnitude vs gamma* envelope (R = 512): scale l, piece sup bits, ")
log("  window width 2^(l+1), first-window row cap bits ~ 0.598*(R+start)")
for l in range(9):
    Wl = [comb(2**l, n)//2 for n in range(1, 2**l)]
    if not Wl:
        log("  l = 0: W empty (scale trivial)"); continue
    P = pmul(onepu(2**l - 1), [0]+Wl)
    sup = max(P)
    log("  l = %d: sup(P_l) bits = %4d, width = %4d, cap bits at row 1 ~ %d"
        % (l, sup.bit_length(), 2**(l+1), int(0.598*(512+1))))
log("")
log("top scale l = t-1 is (1+u)^{R/2-1} W_{t-1}: the R/2-epoch binomial times")
log("the half-Pascal row => doubling recursion skeleton for route (c).")
log("SUMMARY: %s" % ("ALL CHECKS PASSED" if FAIL == 0 else "FAILURES"))
flush()
