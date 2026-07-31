#!/usr/bin/env python3
"""
laneA_referee.py -- Lane A (mechanism extraction) exact-rational referee.

Verifies the "spine reformulation" of AMM 12592 (THM-2160/THM-2225 canon) on
THE known scheme: THM-2225's cyclic-checksum extractor with pathwise deadline
tau <= max(2, 2n-1)  (C = 2, D = -1).

All arithmetic is exact (fractions.Fraction).  No floats anywhere.

Sections
  S0  polynomial utilities (dense Fraction coefficient lists in p, q = 1-p)
  S1  the extractor (THM-2225 rule), implemented from the theorem text
  S2  spine polynomials W_m (node 0^m 1) and V_m (node 1^m 0), m = 1..31,
      via full-shell enumeration AND independent lazy-stopping recursion
  S3  depth accounting: nominal / needed / budget d_m = (C-1)m + D - 1
  S4  deficit tables Delta_m = W_m - 1/2 (exact coefficients, forced cells)
  S5  exact polynomial identities: per-spine block sums, block closure
      D_M = 0, cross-spine transfer +-p^{m'}q^{m'}/2, R_{2^K-1} closed form
  S6  truncation identity (S): R_N in [0, p^{N+1}+q^{N+1}] at 60 rational p,
      and structural enclosure R_N in p^{N+1} G0 + q^{N+1} G1 via exact
      interval chains for the spine values G_{0^t}, G_{1^t}
  S7  cancellation ledger per dyadic block: Hamming-column deficit matrix
  S8  auxiliary diagnostics (bit-flip equivariance, single-defect words)

Exit code 0 and final line "LANE-A REFEREE: ALL CHECKS PASSED" iff everything
verifies.
"""

from fractions import Fraction as F
import math
import sys

CHECKS = []  # (name, ok, detail)


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    tag = "PASS" if ok else "FAIL"
    print(f"[{tag}] {name}" + (f"  -- {detail}" if detail else ""))
    return ok


# ----------------------------------------------------------------------
# S0. Exact polynomial utilities.  Poly = list of Fraction, index = power
#     of p.  q means (1-p) after expansion.
# ----------------------------------------------------------------------

def ptrim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def padd(a, b):
    n = max(len(a), len(b))
    return ptrim([ (a[i] if i < len(a) else F(0)) + (b[i] if i < len(b) else F(0))
                   for i in range(n) ])


def psub(a, b):
    return padd(a, [-c for c in b])


def pmul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    out[i + j] += ai * bj
    return ptrim(out)


def pscal(c, a):
    return ptrim([F(c) * x for x in a])


def ppow(a, n):
    r = [F(1)]
    base = list(a)
    while n:
        if n & 1:
            r = pmul(r, base)
        base = pmul(base, base)
        n >>= 1
    return r


PP = [F(0), F(1)]     # p
QQ = [F(1), F(-1)]    # q = 1 - p


def mon(z, o):
    """p^z q^o expanded in p."""
    return pmul(ppow(PP, z), ppow(QQ, o))


def peval(a, x):
    x = F(x)
    r = F(0)
    for c in reversed(a):
        r = r * x + c
    return r


def pzero(a):
    return all(c == 0 for c in a)


def pstr(a):
    if pzero(a):
        return "0"
    terms = []
    for i, c in enumerate(a):
        if c == 0:
            continue
        if i == 0:
            terms.append(str(c))
        elif i == 1:
            terms.append(f"{c}*p" if c != 1 else "p")
        else:
            terms.append(f"{c}*p^{i}" if c != 1 else f"p^{i}")
    return " + ".join(terms).replace("+ -", "- ")


def pcompose(a, b):
    """a(b(p)) -- Horner in poly arithmetic."""
    r = [F(0)]
    for c in reversed(a):
        r = padd(pmul(r, b), [F(c)])
    return r


HALF = [F(1, 2)]

# ----------------------------------------------------------------------
# S1. THM-2225 cyclic-checksum extractor, straight from the theorem text.
#     Input words use bit values {0,1}; P(0)=p, P(1)=q=1-p.
# ----------------------------------------------------------------------

def run_len(w):
    c = w[0]
    n = 1
    while n < len(w) and w[n] == c:
        n += 1
    return n


def blockM(n):
    """M = 2^ceil(log2(n+1)): the unique power of two with M/2 <= n <= M-1."""
    M = 1
    while M < n + 1:
        M <<= 1
    return M


def output_full(w):
    """THM-2225 output on a full M-bit word of S_M.  Returns 'H' or 'T'."""
    n = run_len(w)
    assert n < len(w), "constant word has no output"
    M = blockM(n)
    assert len(w) == M, (len(w), M)
    if M == 2:
        return 'H' if w == (0, 1) else 'T'
    m2 = M // 2
    # second half y_i = w_{m2+i}, checksum s = sum i*y_i mod m2  (thm (7))
    s = sum(i * w[m2 + i - 1] for i in range(1, m2 + 1)) % m2
    return 'H' if s < m2 // 2 else 'T'          # thm (8)


def popcount(x):
    return bin(x).count("1")


# ----------------------------------------------------------------------
# S2. Spine polynomials.
# ----------------------------------------------------------------------

MMAX = 31   # complete blocks M = 2,4,8,16,32

class SpineNode:
    pass


def spine_data(c, m):
    """Node prefix = c^m (1-c)  (critical value exactly m on the c-spine).
    Returns SpineNode with:
      .M .f (nominal free depth) .table (output per free-bit word, temporal
      order, MSB = first free bit) .wvec .poly .needed .lazy_depth .lazypoly
    """
    d = SpineNode()
    M = blockM(m)
    m2 = M // 2
    prefix = (c,) * m + (1 - c,)
    d.c, d.m, d.M = c, m, M
    if m == M - 1:
        d.f = 0
        d.table = [output_full(prefix)]      # prefix itself has length M
    else:
        f = M - m - 2                        # free bits at positions m+2..M-1
        d.f = f
        table = []
        for idx in range(1 << f):
            bits = tuple((idx >> (f - 1 - i)) & 1 for i in range(f))
            word = prefix + bits
            o0 = output_full(word + (0,))
            o1 = output_full(word + (1,))
            assert o0 == o1, ("bit-M relevance!", c, m, idx)
            table.append(o0)
        d.table = table
    f = d.f
    # w-vector: k = number of ONES among free bits -> weight p^{f-k} q^k
    wvec = [0] * (f + 1)
    for idx, out in enumerate(d.table):
        if out == 'H':
            wvec[popcount(idx)] += 1
    d.wvec = wvec
    poly = [F(0)]
    for k, wk in enumerate(wvec):
        if wk:
            poly = padd(poly, pscal(wk, mon(f - k, k)))
    d.poly = poly
    # needed depth = max temporal coordinate the decision depends on
    needed = 0
    for i in range(1, f + 1):
        maskbit = 1 << (f - i)
        dep = any(d.table[x] != d.table[x | maskbit]
                  for x in range(1 << f) if not (x & maskbit))
        if dep:
            needed = i
    d.needed = needed
    # independent lazy-stopping recursion (minimal stopping rule refinement)
    maxdep = [0]

    def rec(lo, hi, depth):
        first = d.table[lo]
        if all(d.table[i] == first for i in range(lo, hi)):
            maxdep[0] = max(maxdep[0], depth)
            return [F(1)] if first == 'H' else [F(0)]
        mid = (lo + hi) // 2
        return padd(pmul(PP, rec(lo, mid, depth + 1)),
                    pmul(QQ, rec(mid, hi, depth + 1)))

    d.lazypoly = rec(0, 1 << f, 0)
    d.lazy_depth = maxdep[0]
    return d


print("=" * 78)
print("S2. Spine polynomials W_m (node 0^m 1) and V_m (node 1^m 0), m=1..%d" % MMAX)
print("=" * 78)

Wn = {}
Vn = {}
for m in range(1, MMAX + 1):
    Wn[m] = spine_data(0, m)
    Vn[m] = spine_data(1, m)

ok_lazy = all(Wn[m].poly == Wn[m].lazypoly and Vn[m].poly == Vn[m].lazypoly
              for m in range(1, MMAX + 1))
check("S2.route-equality: full-shell enumeration == lazy-stopping recursion, "
      "all m<=%d, both spines" % MMAX, ok_lazy)
ok_lazyd = all(Wn[m].lazy_depth == Wn[m].needed and Vn[m].lazy_depth == Vn[m].needed
               for m in range(1, MMAX + 1))
check("S2.lazy worst-path depth == max dependent coordinate (needed depth)",
      ok_lazyd)

ok_fact = True
for m in range(1, MMAX + 1):
    for d in (Wn[m], Vn[m]):
        for k, wk in enumerate(d.wvec):
            if not (0 <= wk <= math.comb(d.f, k)):
                ok_fact = False
check("S2.FACT containment: 0 <= w_k <= binom(f,k) for every cell (m,k), "
      "both spines", ok_fact)

for m in range(1, 16):
    print(f"  W_{m:<2} = {pstr(Wn[m].poly):<40}  V_{m:<2} = {pstr(Vn[m].poly)}")

# ----------------------------------------------------------------------
# S3. Depth accounting.
# ----------------------------------------------------------------------

print()
print("=" * 78)
print("S3. Depth accounting.  budget d_m = T(m)-(m+1), T(n)=max(2,2n-1)")
print("    nominal = theorem stop depth = M-m-2 (m<=M-2) or 0 (m=M-1)")
print("=" * 78)
print(f"{'m':>3} {'M':>3} {'budget':>6} {'nominalW':>8} {'neededW':>7} "
      f"{'degW':>4} {'nominalV':>8} {'neededV':>7} {'degV':>4}  binds")
ok_depth = True
bind_rows = []
for m in range(1, MMAX + 1):
    M = blockM(m)
    budget = 0 if m == 1 else m - 2
    nomW = Wn[m].f
    nomV = Vn[m].f
    degW = len(Wn[m].poly) - 1 if not pzero(Wn[m].poly) else 0
    degV = len(Vn[m].poly) - 1 if not pzero(Vn[m].poly) else 0
    if not (Wn[m].needed <= nomW <= budget and Vn[m].needed <= nomV <= budget):
        ok_depth = False
    binds = []
    if nomW == budget and m >= 2:
        binds.append("nomW")
    if Wn[m].needed == budget and m >= 2:
        binds.append("needW")
    if Vn[m].needed == budget and m >= 2:
        binds.append("needV")
    bind_rows.append((m, binds))
    print(f"{m:>3} {M:>3} {budget:>6} {nomW:>8} {Wn[m].needed:>7} {degW:>4} "
          f"{nomV:>8} {Vn[m].needed:>7} {degV:>4}  {','.join(binds) or '-'}")
check("S3.depths: needed <= nominal <= budget=(C-1)m+D-1 for C=2,D=-1 "
      "(m>=2; m=1 via max(2,.) clause), all m, both spines", ok_depth)

# theorem stop bound pathwise
ok_stop = True
for m in range(1, MMAX + 1):
    M = blockM(m)
    absstop_nominal = M if m == M - 1 else M - 1
    if absstop_nominal > max(2, 2 * m - 1):
        ok_stop = False
check("S3.pathwise deadline: nominal absolute stop <= max(2,2m-1) for all "
      "critical values m<=%d" % MMAX, ok_stop)

# needed-depth binding pattern at block openers
openers = [2, 4, 8, 16]
print("\n  block openers m'=2^k: needed depths  "
      + "  ".join(f"m={m}: W {Wn[m].needed}/V {Vn[m].needed} (budget {m-2})"
                  for m in openers))
shellLB = {m: m // 2 - 1 for m in openers}   # THM-2160 (27): d>=h/2 => needed>=h/2-1
ok_lb = all(Wn[m].needed >= shellLB[m] and Vn[m].needed >= shellLB[m]
            for m in openers if m >= 2)
check("S3.shell-model floor: needed >= m/2-1 at openers (THM-2160 sec 6.2, "
      "d>=h/2 with d=needed+1)", ok_lb,
      detail=" ".join(f"m={m}:W{Wn[m].needed},V{Vn[m].needed}>=" +
                      str(shellLB[m]) for m in openers))

# ----------------------------------------------------------------------
# S4. Deficits Delta_m = W_m - 1/2.
# ----------------------------------------------------------------------

print()
print("=" * 78)
print("S4. Deficit vectors delta_{m,k} = w_k - binom(f,k)/2 (natural depth f)")
print("=" * 78)

GRID = [F(i, 240) for i in range(1, 240)]


def supgrid(poly):
    best, arg = F(0), F(0)
    for x in GRID:
        v = abs(peval(poly, x))
        if v > best:
            best, arg = v, x
    return best, arg


for m in range(1, 16):
    for lab, d in (("W", Wn[m]), ("V", Vn[m])):
        f = d.f
        binoms = [math.comb(f, k) for k in range(f + 1)]
        delta = [F(w) - F(b, 2) for w, b in zip(d.wvec, binoms)]
        forced = [k for k, b in enumerate(binoms) if b % 2 == 1]
        dpoly = psub(d.poly, HALF)
        s, arg = supgrid(dpoly)
        print(f"  m={m:<2} {lab}: f={f}  w={d.wvec}  binom={binoms}")
        print(f"        delta={[str(x) for x in delta]}  forced(odd binom) k={forced}")
        print(f"        Delta(p) = {pstr(dpoly)}   sup_grid|Delta|={s} at p={arg}")

# half-integer forcing is real: every cell with odd binom(f,k) has
# delta in (1/2)Z \ Z
ok_forced = True
for m in range(1, MMAX + 1):
    for d in (Wn[m], Vn[m]):
        for k in range(d.f + 1):
            b = math.comb(d.f, k)
            delta2 = 2 * d.wvec[k] - b   # = 2*delta, integer
            if b % 2 == 1 and delta2 % 2 == 0:
                ok_forced = False
check("S4.parity: cells with odd binom(f,k) carry half-integer deficits "
      "(2*delta odd), all m<=%d" % MMAX, ok_forced)

# ----------------------------------------------------------------------
# S5. Exact polynomial identities (proved by exact coefficient arithmetic).
# ----------------------------------------------------------------------

print()
print("=" * 78)
print("S5. Exact identities (coefficient-level, hence polynomial identities)")
print("=" * 78)


def blockW(mp):
    s = [F(0)]
    for m in range(mp, 2 * mp):
        s = padd(s, pmul(mon(m, 1), Wn[m].poly))
    return s


def blockV(mp):
    s = [F(0)]
    for m in range(mp, 2 * mp):
        s = padd(s, pmul(mon(1, m), Vn[m].poly))
    return s


ok = True
for mp in (2, 4, 8, 16):
    tgtW = pscal(F(1, 2), pmul(mon(mp, 0),
                               psub(psub([F(1)], mon(mp, 0)), mon(0, mp))))
    tgtV = pmul(mon(0, mp),
                padd(mon(mp, 0),
                     pscal(F(1, 2), psub(psub([F(1)], mon(mp, 0)), mon(0, mp)))))
    okW = blockW(mp) == tgtW
    okV = blockV(mp) == tgtV
    ok = ok and okW and okV
    print(f"  m'={mp:<2}: sum p^m q W_m = p^m'(1-p^m'-q^m')/2 : {okW}   "
          f"sum q^m p V_m = q^m'(p^m'+(1-p^m'-q^m')/2) : {okV}")
check("S5.per-spine block sums (closed forms), m'=2,4,8,16", ok)

# M=2 anomaly
okA = blockW(1) == mon(1, 1) and pzero(blockV(1))
check("S5.M=2 anomaly: blockW(1)=pq (0-spine SURPLUS +q/2 rel.), blockV(1)=0 "
      "(1-spine deficit) -- sign REVERSED vs all later blocks", okA)

# block closure and cross-spine transfer
ok1 = ok2 = ok3 = True
for mp in (1, 2, 4, 8, 16):
    tot = padd(blockW(mp), blockV(mp))
    tgt = pscal(F(1, 2), psub(padd(mon(mp, 0), mon(0, mp)),
                              padd(mon(2 * mp, 0), mon(0, 2 * mp))))
    if tot != tgt:
        ok1 = False
    # deficit form: D_M = blockW+blockV - (1/2) sum (p^m q + q^m p) = 0
    fair = [F(0)]
    for m in range(mp, 2 * mp):
        fair = padd(fair, pscal(F(1, 2), padd(mon(m, 1), mon(1, m))))
    if not pzero(psub(tot, fair)):
        ok2 = False
    # per-spine deficit flow
    fw = psub(blockW(mp), pscal(F(1, 2), psub(mon(mp, 0), mon(2 * mp, 0))))
    fv = psub(blockV(mp), pscal(F(1, 2), psub(mon(0, mp), mon(0, 2 * mp))))
    trans = pscal(F(1, 2), mon(mp, mp))
    if mp == 1:
        if not (fw == trans and fv == pscal(-1, trans)):
            ok3 = False
    else:
        if not (fw == pscal(-1, trans) and fv == trans):
            ok3 = False
check("S5.block closure: blockW+blockV = (p^m'+q^m'-p^2m'-q^2m')/2, "
      "m'=1,2,4,8,16", ok1)
check("S5.block deficit identity D_M = 0 exactly (all deficit cancellation "
      "is INTRA-block)", ok2)
check("S5.cross-spine transfer: sum p^m q Delta_m = -p^m' q^m'/2 and "
      "sum q^m p Delta'_m = +p^m' q^m'/2 per block (sign flipped at M=2)", ok3)

# R_N closed form at dyadic N
Tm = {m: pmul(mon(m, 1), Wn[m].poly) for m in Wn}
Um = {m: pmul(mon(1, m), Vn[m].poly) for m in Vn}
RN = {}
acc = [F(0)]
for m in range(1, MMAX + 1):
    acc = padd(acc, padd(Tm[m], Um[m]))
    RN[m] = psub(HALF, acc)
RN[0] = HALF

okR = True
for K in (1, 2, 3, 4, 5):
    N = (1 << K) - 1
    tgt = pscal(F(1, 2), padd(mon(1 << K, 0), mon(0, 1 << K)))
    good = RN[N] == tgt
    okR = okR and good
    print(f"  K={K}: R_(2^{K}-1) == (p^(2^{K})+q^(2^{K}))/2 : {good}")
check("S5.dyadic remainder closed form R_{2^K-1} = (p^{2^K}+q^{2^K})/2, "
      "K=1..5", okR)

# ----------------------------------------------------------------------
# S6. Truncation identity (S) at rational points + structural enclosures.
# ----------------------------------------------------------------------

print()
print("=" * 78)
print("S6. R_N bounds and structural spine-value enclosures at 60 rational p")
print("=" * 78)

PTS = [F(k, 61) for k in range(1, 61)]
JTOP = 13     # tail cut: G_{0^{2^JTOP}} in [0,1]

ok_bounds = True
for N in range(0, MMAX):
    for x in PTS:
        r = peval(RN[N], x)
        hi = x ** (N + 1) + (1 - x) ** (N + 1)
        if not (0 <= r <= hi):
            ok_bounds = False
check("S6.demanded bound: 0 <= R_N <= p^{N+1}+q^{N+1}  (so R_N/(p^{N+1}+"
      "q^{N+1}) in [0,1]) for N=0..%d at 60 rational p" % (MMAX - 1), ok_bounds)


def chain_u(x, spine):
    """Exact enclosures [lo,hi] for u_j = G_{0^{2^j}} (spine=0) or
    v_j = G_{1^{2^j}} (spine=1), j = 1..5, via the proven per-spine block
    recursion; tail value in [0,1] at level JTOP."""
    q = 1 - x
    lo, hi = F(0), F(1)
    vals = {}
    for j in range(JTOP - 1, 0, -1):
        e = 1 << j
        pe, qe = x ** e, q ** e
        if spine == 0:
            base = (1 - pe - qe) / 2
            lo, hi = base + pe * lo, base + pe * hi
        else:
            base = (1 + pe - qe) / 2
            lo, hi = base + qe * lo, base + qe * hi
        if j <= 5:
            vals[j] = (lo, hi)
    return vals


ok_struct = True
worst_gap = F(0)
for x in PTS:
    q = 1 - x
    U = chain_u(x, 0)
    V = chain_u(x, 1)
    for N in range(0, MMAX):
        t = N + 1
        Mb = blockM(t)
        j = Mb.bit_length() - 1
        # G_{0^t} = sum_{m=t}^{Mb-1} p^{m-t} q W_m + p^{Mb-t} u_j
        f0 = sum((x ** (m - t)) * q * peval(Wn[m].poly, x)
                 for m in range(t, Mb))
        f1 = sum((q ** (m - t)) * x * peval(Vn[m].poly, x)
                 for m in range(t, Mb))
        g0lo, g0hi = f0 + x ** (Mb - t) * U[j][0], f0 + x ** (Mb - t) * U[j][1]
        g1lo, g1hi = f1 + q ** (Mb - t) * V[j][0], f1 + q ** (Mb - t) * V[j][1]
        lo = x ** t * g0lo + q ** t * g1lo
        hi = x ** t * g0hi + q ** t * g1hi
        r = peval(RN[N], x)
        if not (lo <= r <= hi):
            ok_struct = False
        worst_gap = max(worst_gap, hi - lo)
        # sanity: G enclosures inside [0,1]
        if not (0 <= g0lo <= g0hi <= 1 and 0 <= g1lo <= g1hi <= 1):
            ok_struct = False
check("S6.structural: R_N inside p^{N+1}[G_{0^{N+1}}] + q^{N+1}[G_{1^{N+1}}] "
      "exact interval enclosures, N=0..%d, 60 points" % (MMAX - 1), ok_struct,
      detail=f"max enclosure width {float(worst_gap):.3e}")

# root fairness through the spine decomposition (N=0 row is exactly that)
check("S6.root fairness: 1/2 = p G_0 + q G_1 within enclosures (the N=0 row)",
      ok_struct)

# ----------------------------------------------------------------------
# S7. Cancellation ledger per dyadic block.
# ----------------------------------------------------------------------

print()
print("=" * 78)
print("S7. Hamming-column deficit ledger per block (entry = #H - |cell|/2)")
print("=" * 78)


def ledger(M):
    """columns j=1..M-1 (total ones of the M-word), rows (spine c, run m)."""
    m2 = M // 2
    rows = {}
    for c in (0, 1):
        for y in range(1 << m2):
            bits = tuple((y >> (m2 - 1 - i)) & 1 for i in range(m2))
            if all(b == c for b in bits):
                continue                      # constant word excluded
            w = (c,) * m2 + bits
            m = run_len(w)
            j = sum(w)
            out = output_full(w)
            key = (c, m)
            rows.setdefault(key, {})
            h, n = rows[key].get(j, (0, 0))
            rows[key][j] = (h + (1 if out == 'H' else 0), n + 1)
    return rows


ok_led = True
for M in (4, 8, 16, 32):
    rows = ledger(M)
    m2 = M // 2
    colsum = {}
    for (c, m), cols in rows.items():
        for j, (h, n) in cols.items():
            colsum[j] = colsum.get(j, F(0)) + (F(h) - F(n, 2))
    if any(v != 0 for v in colsum.values()):
        ok_led = False
    # column purity
    for (c, m), cols in rows.items():
        for j in cols:
            if j < m2 and c != 0:
                ok_led = False
            if j > m2 and c != 1:
                ok_led = False
    if M <= 16:
        print(f"\n  block M={M} (m'={m2}); rows (spine,m), cols j; "
              f"entries 2*(#H-|cell|/2):")
        hdr = "   ".join(f"j={j}" for j in range(1, M))
        print(f"      {'row':>10}  {hdr}")
        for c in (0, 1):
            for m in range(m2, M):
                key = (c, m)
                if key not in rows:
                    continue
                ent = []
                for j in range(1, M):
                    if j in rows[key]:
                        h, n = rows[key][j]
                        ent.append(f"{2*h-n:>3}")
                    else:
                        ent.append("  .")
                print(f"      ({c},{m:>2})    " + "   ".join(ent))
    else:
        tot = sum(abs(F(h) - F(n, 2)) for cols in rows.values()
                  for (h, n) in cols.values())
        print(f"\n  block M={M}: {len(rows)} rows, total |deficit| mass "
              f"sum = {tot}; all {M-1} columns sum to 0: "
              f"{all(v == 0 for v in colsum.values())}")
check("S7.ledger: every Hamming column of every block M=4,8,16,32 sums to 0; "
      "columns j<m' pure 0-spine, j>m' pure 1-spine, j=m' = boundary pair",
      ok_led)

# middle-column boundary pair carries the whole cross-spine transfer
ok_pair = True
for M in (4, 8, 16, 32):
    m2 = M // 2
    w0 = (0,) * m2 + (1,) * m2
    w1 = (1,) * m2 + (0,) * m2
    if not (output_full(w0) == 'T' and output_full(w1) == 'H'):
        ok_pair = False
check("S7.boundary pair: 0^m'1^m' -> T and 1^m'0^m' -> H for M=4,8,16,32 "
      "(the transfer words; note M=2 is 01->H, 10->T, reversed)", ok_pair)

# ----------------------------------------------------------------------
# S8. Auxiliary diagnostics.
# ----------------------------------------------------------------------

print()
print("=" * 78)
print("S8. Diagnostics")
print("=" * 78)

# bit-flip equivariance V_m(p) ?= 1 - W_m(1-p)
eqv = []
for m in range(1, MMAX + 1):
    flip = psub([F(1)], pcompose(Wn[m].poly, QQ))
    eqv.append((m, Vn[m].poly == flip))
eq_true = [m for m, e in eqv if e]
eq_false = [m for m, e in eqv if not e]
print(f"  bit-flip equivariance V_m(p) == 1 - W_m(1-p): holds for m={eq_true}")
print(f"    fails for m={eq_false}")
check("S8.equivariance holds EXACTLY for blocks M<=4 plus the trivial "
      "block-tail constants m=M-1; broken everywhere else from M=8 on",
      set(eq_true) == {1, 2, 3, 7, 15, 31})

# single-defect words (the p->1 / p->0 boundary constraint carriers)
seqW = "".join(
    output_full((0,) * m + (1,) + (0,) * (blockM(m) - m - 1))
    for m in range(1, MMAX + 1))
seqV = "".join(
    output_full((1,) * m + (0,) + (1,) * (blockM(m) - m - 1))
    for m in range(1, MMAX + 1))
print(f"  single-defect 0-spine words 0^m 1 0^*: {seqW}")
print(f"  single-defect 1-spine words 1^m 0 1^*: {seqV}")
ok_sd = True
pos = 0
for K in (1, 2, 3, 4, 5):
    m2 = 1 << (K - 1)
    seg = seqW[m2 - 1: 2 * m2 - 1]
    segV = seqV[m2 - 1: 2 * m2 - 1]
    hW, hV = seg.count('H'), segV.count('H')
    print(f"    block m'={m2:>2}: 0-spine {seg} (#H={hW}/{len(seg)})   "
          f"1-spine {segV} (#H={hV}/{len(segV)})")
    if m2 >= 2 and (2 * hW != len(seg) or 2 * hV != len(segV)):
        ok_sd = False
check("S8.single-defect balance: within every block M>=4 the class-j=1 "
      "(resp. j=M-1) words are exactly half H (Abel constraint at p->1, "
      "p->0 satisfied block-locally)", ok_sd)

# needed-depth profile summary for findings
print("\n  needed-depth profiles (W;V) by block:")
for K in (1, 2, 3, 4, 5):
    m2 = 1 << (K - 1)
    prof = "  ".join(f"m={m}:{Wn[m].needed};{Vn[m].needed}"
                     for m in range(m2, 2 * m2))
    print(f"    M={2*m2:>2}: {prof}")

# ----------------------------------------------------------------------
# verdict
# ----------------------------------------------------------------------

print()
bad = [name for name, ok, _ in CHECKS if not ok]
if bad:
    print("LANE-A REFEREE: FAILURES:")
    for b in bad:
        print("  -", b)
    sys.exit(1)
print(f"LANE-A REFEREE: ALL CHECKS PASSED ({len(CHECKS)} checks)")
