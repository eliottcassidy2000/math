#!/usr/bin/env python3
"""ANGLE B3 step 3A -- exact capacity ledger for the diagonal-atom transportation.

Reduced problem (proved in amm12592_allR_GR_ybox_reduction_boxeph.py):
  cells y_{i,k} in [-C(d_i-1,k), C(d_i-1,k-1)],  sum_i x^i gamma_i = G_R.
Lucas diagonal expansion (amm12592_allR_GR_lucas_diagonal_boxeph.py):
  G_R = 1 + m_1 p + sum_{k=2}^{R/2} m_k p^k q^{k-1},  m_k = a_k/2,
  a_k = (-1)^k [C(R-k,k)+C(R-k-1,k-1)],  m_1 = -R/2.

Atom k placed (as far as possible) in row i <= k has row content
mu * p^{k-i} q^{k-1}; its cells at degree d = d_i are mu * s_j,
  s_j = C(M, j-(k-1)),  M := d - (2k-1-i)  (need M >= 0),
support j in [k-1, d-(k-i)].  The EXACT per-row coefficient capacity is
  cap+(i,k) = min_j floor( C(d-1, j-1) / s_j )   (mu > 0)
  cap-(i,k) = min_j floor( C(d-1, j)   / s_j )   (mu < 0)
(o the empty-box convention C(d-1,-1) = C(d-1,d) = 0 kills nothing here since
support starts at j = k-1 >= 1 and ends at j = d-(k-i) <= d - ... only if
k > i; for i = k the support reaches j = d: top cell forces cap+ <= 1 there).

LEDGER per dyadic R: for each k: |m_k|, sign, best single-row cap (arg i),
window sum of caps over ALL legal rows, ratio window/|m_k|.  Also the head
ledger: pure-p demand c = R/2 - 1 (one -1 marker per row, value ledger at x=1)
and the correction demand N_u = coefficient of -p^u q after the marker
telescope with markers at odd rows (exponents i+1).

VERDICTS this script decides EXACTLY:
  (V1) is single-row placement feasible for every k?  (expected: NO near peak)
  (V2) is the windowed (all-rows, private-capacity) sum >= |m_k| for every k?
       -- necessary for any per-atom split WITHOUT cross-row cancellation.
  (V3) what fraction of each row's box does the naive proportional split use
       (upper bound: sum over k of |share| * s_j vs box, worst cell) --
       NOT computed here in full; step 3B does the actual explicit split.

Search negatives here are NOT infeasibility of (*) -- cross-row cancellation
(carries) is excluded by design in this ledger; the point is to map EXACTLY
where the diagonal route needs carries.
"""
import sys
from math import comb

WT = "/tmp/math-wt-boxeph-multifront"
OUT = WT + "/05-knowledge/results/amm12592_allR_GR_capacity_ledger_boxeph.out"

log_lines = []
def log(s=""):
    print(s); log_lines.append(s)
def flush():
    with open(OUT, "w") as f: f.write("\n".join(log_lines) + "\n")

def _fp(n):
    if n == 0: return (0,1)
    a,b = _fp(n>>1); c = a*(2*b-a); d2 = a*a+b*b
    return (d2, c+d2) if n & 1 else (c, d2)
def _le(d, m):
    if d < 0: return True
    F, F1 = _fp(2*m); L = 2*F1 - F
    A = 2*5**d - L
    return A <= 0 or A*A < 5*F*F
def floorgs(m):
    d = int(m*0.5979874356654402)
    while _le(d+1, m): d += 1
    while not _le(d, m): d -= 1
    return d

def m_coef(R, k):
    if k == 1: return -R//2
    a = (-1)**k * (comb(R-k, k) + comb(R-k-1, k-1))
    assert a % 2 == 0
    return a//2

def caps_for(R, k, D0=0):
    """returns list of (i, cap) for the correct sign of m_k, plus best/(sum)."""
    mk = m_coef(R, k)
    sgn = 1 if mk > 0 else -1
    out = []
    for i in range(0, min(k, R-2) + 1):
        d = floorgs(R + i) + D0
        M = d - (2*k - 1 - i)
        if M < 0: continue
        lo = k - 1; hi = d - (k - i)
        cap = None
        for j in range(lo, hi + 1):
            s = comb(M, j - (k-1))
            box = comb(d-1, j-1) if sgn > 0 else comb(d-1, j)
            # box entries: upper C(d-1,j-1) (j>=1), lower C(d-1,j) (j<=d-1)
            if sgn > 0 and j == 0: box = 0
            if sgn < 0 and j == d: box = 0
            c = box // s
            cap = c if cap is None else min(cap, c)
            if cap == 0: break
        out.append((i, cap))
    return mk, out

def ledger(R, D0=0, kstep=1):
    log("="*90)
    log("R = %d, D0 = %d   (profile d_i = floor(gamma*(R+i)) + D0)" % (R, D0))
    log("="*90)
    log("  k | sign |m_k| bits | best single cap bits @row | winsum bits | single-OK | winsum-OK")
    v1 = True; v2 = True; worst = None
    for k in range(2, R//2 + 1, kstep):
        mk, caps = caps_for(R, k, D0)
        need = abs(mk)
        if not caps:
            log("  %4d : NO LEGAL ROW" % k); v1 = v2 = False; continue
        best_i, best = max(((i, c) for i, c in caps), key=lambda t: t[1])
        tot = sum(c for _, c in caps)
        ok1 = best >= need; ok2 = tot >= need
        v1 = v1 and ok1; v2 = v2 and ok2
        ratio_bits = (tot.bit_length() - need.bit_length())
        if worst is None or ratio_bits < worst[1]:
            worst = (k, ratio_bits)
        if (k <= 8 or k % max(1, R//32) == 0 or not ok1 or not ok2):
            log("  %4d | %+d | %6d | %6d @ i=%-4d | %8d | %5s | %5s" %
                (k, 1 if mk > 0 else -1, need.bit_length(), best.bit_length(),
                 best_i, tot.bit_length(), ok1, ok2))
    log("  VERDICT V1 (every k single-row placeable): %s" % v1)
    log("  VERDICT V2 (every k window-sum placeable, private caps): %s" % v2)
    if worst: log("  tightest window margin: k = %d, winsum/|m_k| ~ 2^%d" % worst)
    # head ledger
    c = R//2 - 1
    log("  head: pure-p demand c = R/2-1 = %d markers (rows with gamma_i(1) = -1)" % c)
    Ns = {}
    for i in range(1, R-2, 2):           # markers at odd rows, exponent E = i+1
        for u in range(1, i+1):
            Ns[u] = Ns.get(u, 0) + 1
    if Ns:
        us = sorted(Ns)
        log("  correction demand -N_u p^u q: u=1..%d, N_1 = %d, N_2 = %d, max N = %d"
            % (max(us), Ns.get(1,0), Ns.get(2,0), max(Ns.values())))
        # exact caps for shape p^a q (a = u - i) rows i<u -- the u=1 obstruction
        for u in (1, 2, 3):
            if u not in Ns: continue
            tot = 0
            for i in range(0, u):
                d = floorgs(R+i) + D0
                a = u - i; M = d - a - 1
                if M < 0: continue
                capn = None
                for j in range(1, d - a + 1):
                    s = comb(M, j-1)
                    box = comb(d-1, j) if j <= d-1 else 0
                    cc = box // s
                    capn = cc if capn is None else min(capn, cc)
                    if capn == 0: break
                tot += capn or 0
            log("    u = %d: demand N_u = %d, total private cap over rows i<u = %d  %s"
                % (u, Ns[u], tot, "OK" if tot >= Ns[u] else "NEEDS CARRIES"))
    flush()

for R in (8, 16, 32, 64, 128):
    ledger(R)
ledger(256, kstep=1)
log("")
log("NOTE: V1/V2 are about PRIVATE per-atom capacity without cross-row")
log("cancellation.  V2 = True for all k would make an explicit split-only")
log("construction plausible; False marks exactly where carries are mandatory.")
flush()
