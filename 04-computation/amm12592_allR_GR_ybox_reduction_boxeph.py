#!/usr/bin/env python3
"""ANGLE B3 step 1 -- exact bookkeeping for the direct G_R representation.

Conventions (THM-3302 / THM-3329 line): p = x, q = 1-x, epoch identity
  (*)  sum_{i=0}^{R-1} x^i Delta_i = q^{R-1},
Delta_i admissible at d_i (Lucas box: |delta_k| <= C(d,k), delta_k == C(d,k) mod 2)
in Bernstein basis B_{d,k} = x^{d-k} q^k.  Ballot normal form (S2):
  Delta_i = (p-q) + c_i (i <= R-2),  Delta_{R-1} = -1,
  (*)  <=>  sum_{i<=R-2} x^i c_i = q^{R-1} - E_{R-1} =: 2 G_R,  E_m = -1+x+...+x^m.

THIS SCRIPT proves/verifies exactly:

(A) [S3 recompute] [x^0] 2G_R = 2 and [x^j] 2G_R = (-1)^j C(R-1,j) - 1 (j>=1);
    2G_R even coefficientwise iff R = 2^t.

(B) [Y-BOX LEMMA] Write c_i = 2 gamma_i.  Let y_k be the (unique) Bernstein
    cells of gamma_i at degree d.  Then
       Delta_i = (p-q) + 2 gamma_i  is Lucas-box admissible at d
       <=>   -C(d-1,k) <= y_k <= C(d-1,k-1)   for all 0 <= k <= d
    (C(d-1,-1) = C(d-1,d) = 0).  Parity disappears entirely.
    Proof: ballot cells b_k = C(d-1,k) - C(d-1,k-1) satisfy b_k == C(d,k) mod 2
    (Pascal), so admissibility of b + 2y <=> |b_k + 2 y_k| <= C(d,k)
    <=> -C(d,k) - b_k <= 2 y_k <= C(d,k) - b_k, and C(d,k) = C(d-1,k)+C(d-1,k-1)
    gives exactly -2C(d-1,k) <= 2y_k <= 2C(d-1,k-1).  Verified exhaustively for
    small d and on random blocks for larger d.

(C) [CORNERS] lower corner y_k = -C(d-1,k) is gamma = -p (Delta = -1);
    upper corner y_k = +C(d-1,k-1) is gamma = +q (Delta = +1).

(D) [BALLOT CUT in y] y_0 = gamma(1) in {-1,0} and y_d = gamma(0) in {0,1}.

(E) [ATOM FEASIBILITY TABLE] for atoms epsilon * p^a q^b (epsilon=+-1, a+b <= d):
    cells epsilon*C(d-a-b, k-b); feasible alone iff (epsilon=+1 and b>=1) or
    (epsilon=-1 and a>=1).  Verified by direct box check.

(F) [WITNESS MAP] the R = 8 rule witness (and R = 16, 32, 64, 128 floor
    witnesses) map to y-cell arrays: gamma_i = (Delta_i - (p-q))/2, all rows in
    the y-box, and sum_i x^i gamma_i = G_R exactly (rows 0..R-2; row R-1 is the
    constant -1 = lower corner at gamma = -p ... actually Delta_{R-1} may be a
    general admissible block in a witness; we verify the normal-form reduction
    on witnesses that carry the backbone, and verify the master equation
    directly in all cases).

Exact int arithmetic only.  No numpy.
"""
import json, os, random, sys
from math import comb

WT = "/tmp/math-wt-boxeph-multifront"
CP = WT + "/04-computation"
OUT = WT + "/05-knowledge/results/amm12592_allR_GR_ybox_reduction_boxeph.out"

log_lines = []
def log(s=""):
    print(s)
    log_lines.append(s)
def flush():
    with open(OUT, "w") as f:
        f.write("\n".join(log_lines) + "\n")

FAIL = 0
def check(name, ok, detail=""):
    global FAIL
    tag = "PASS" if ok else "FAIL"
    if not ok: FAIL += 1
    log("  [%s] %s%s" % (tag, name, ("  -- " + detail) if detail else ""))
    return ok

# ---------- exact poly ops (independent re-implementation) ----------
def ptrim(a):
    a = list(a)
    while a and a[-1] == 0: a.pop()
    return a
def padd(a, b):
    r = [0]*max(len(a), len(b))
    for i,v in enumerate(a): r[i] += v
    for i,v in enumerate(b): r[i] += v
    return ptrim(r)
def pneg(a): return [-v for v in a]
def psub(a, b): return padd(a, pneg(b))
def pmul(a, b):
    if not a or not b: return []
    r = [0]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        if u:
            for j,v in enumerate(b): r[i+j] += u*v
    return ptrim(r)
def pshift(a, s): return ([0]*s + list(a)) if a else []
def pscale(a, c): return [c*v for v in a] if c else []
def qpow(m): return [((-1)**k)*comb(m,k) for k in range(m+1)]

def bern_to_poly_slow(cells, d):
    P = []
    for k, c in enumerate(cells):
        if c: P = padd(P, pscale(pshift(qpow(k), d-k), c))
    return P

def bern_to_poly(cells, d):
    """Horner in q: U_d = delta_d; U_k = delta_k x^{d-k} + q U_{k+1}; Delta = U_0.
    q-multiply is shift-subtract (big-int adds only)."""
    if len(cells) != d+1:
        return bern_to_poly_slow(cells, d)
    U = [0]*(d+1)
    U[0] = cells[d]
    for k in range(d-1, -1, -1):
        # U <- q*U  (q = 1 - x): new[i] = U[i] - U[i-1]
        prev = 0
        for i in range(d+1):
            cur = U[i]
            U[i] = cur - prev
            prev = cur
        U[d-k] += cells[k]
    return ptrim(U)
def poly_to_bern(P, d):
    assert len(P) <= d+1, (len(P), d)
    return [sum(P[j]*comb(d-j, t) for j in range(len(P))) for t in range(d+1)]
def ballot(d):
    return [comb(d-1,k) - (comb(d-1,k-1) if k >= 1 else 0) for k in range(d+1)]
def admissible(delta, d):
    return len(delta) == d+1 and all(
        abs(c) <= comb(d,k) and (c - comb(d,k)) % 2 == 0 for k,c in enumerate(delta))
def ybox_ok(y, d):
    """the reduced box: -C(d-1,k) <= y_k <= C(d-1,k-1)"""
    if len(y) != d+1: return False
    for k, v in enumerate(y):
        lo = -comb(d-1, k) if k <= d-1 else 0
        hi = comb(d-1, k-1) if k >= 1 else 0
        if not (lo <= v <= hi): return False
    return True

def Em(m): return ptrim([-1] + [1]*m) if m >= 0 else []
def twoG(R): return psub(qpow(R-1), Em(R-1))

# exact floor(m*gamma*), gamma* = log_5(phi^2), via Fibonacci/Lucas (independent)
def _fib_pair(n):
    if n == 0: return (0, 1)
    a, b = _fib_pair(n >> 1)
    c = a*(2*b - a); d2 = a*a + b*b
    return (d2, c + d2) if n & 1 else (c, d2)
def five_pow_le_phi2m(d, m):
    if d < 0: return True
    F, F1 = _fib_pair(2*m); L = 2*F1 - F
    A = 2*5**d - L
    if A <= 0: return True
    return A*A < 5*F*F
def floor_gamma_star(m):
    d = int(m*0.5979874356654402)
    while five_pow_le_phi2m(d+1, m): d += 1
    while not five_pow_le_phi2m(d, m): d -= 1
    return d

log("="*78)
log("(A) S3 recompute: coefficients of 2G_R = q^{R-1} - E_{R-1}")
log("="*78)
ok = True
for R in range(2, 1025):
    K = twoG(R)
    Kfull = K + [0]*(R - len(K))
    if Kfull[0] != 2: ok = False; break
    for j in range(1, R):
        want = (-1)**j * comb(R-1, j) - 1
        if (Kfull[j] if j < len(Kfull) else 0) != want: ok = False; break
    if not ok: break
check("[x^0]2G_R = 2 and [x^j]2G_R = (-1)^j C(R-1,j) - 1 (j>=1), R=2..1024", ok)
dy = [R for R in range(2, 1025) if all(v % 2 == 0 for v in twoG(R))]
check("2G_R even coefficientwise iff R = 2^t (R=2..1024)",
      dy == [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024], str(dy[:12]))
flush()

log("")
log("="*78)
log("(B) Y-BOX LEMMA: Delta=(p-q)+2*gamma admissible at d  <=>  cells y of gamma")
log("    at degree d satisfy -C(d-1,k) <= y_k <= C(d-1,k-1).  Parity vanishes.")
log("="*78)
# proof-side identity: ballot parity b_k == C(d,k) mod 2 and capacity split
ok = True
for d in range(1, 200):
    b = ballot(d)
    for k in range(d+1):
        if (b[k] - comb(d,k)) % 2 != 0: ok = False
        if comb(d,k) != comb(d-1,k) + (comb(d-1,k-1) if k >= 1 else 0):
            if not (k == 0 and comb(d,0) == comb(d-1,0)): ok = False
check("ballot parity b_k == C(d,k) mod 2 and Pascal split, d=1..199", ok)

# exhaustive equivalence for small d: enumerate ALL admissible Delta
ok = True; count = 0
for d in range(1, 6):
    b = ballot(d)
    # enumerate all admissible delta vectors
    def rec(k, delta):
        global ok, count
        if k == d+1:
            count += 1
            y = [(delta[t] - b[t])//2 for t in range(d+1)]
            # y must be the Bernstein cells of gamma=(Delta-(p-q))/2 -- check
            gam = [(v) for v in y]
            P = bern_to_poly(y, d)
            y2 = poly_to_bern(P, d)
            if y2 != y: ok = False
            if not ybox_ok(y, d): ok = False
            return
        c0 = comb(d,k)
        for v in range(-c0, c0+1):
            if (v - c0) % 2 == 0:
                rec(k+1, delta + [v])
    rec(0, [])
check("all admissible Delta at d=1..5 map to in-box y (exhaustive, %d blocks)" % count, ok)

# converse: every in-box y gives admissible Delta (exhaustive small d)
ok = True; count = 0
for d in range(1, 6):
    b = ballot(d)
    def rec2(k, y):
        global ok, count
        if k == d+1:
            count += 1
            delta = [b[t] + 2*y[t] for t in range(d+1)]
            if not admissible(delta, d): ok = False
            return
        lo = -comb(d-1,k) if k <= d-1 else 0
        hi = comb(d-1,k-1) if k >= 1 else 0
        for v in range(lo, hi+1):
            rec2(k+1, y + [v])
    rec2(0, [])
check("all in-box y at d=1..5 map to admissible Delta (exhaustive, %d cells arrays)" % count, ok)

# random large-d cross-check
random.seed(20260803)
ok = True
for trial in range(4000):
    d = random.randint(6, 40)
    y = []
    for k in range(d+1):
        lo = -comb(d-1,k) if k <= d-1 else 0
        hi = comb(d-1,k-1) if k >= 1 else 0
        y.append(random.randint(lo, hi))
    delta = [ballot(d)[t] + 2*y[t] for t in range(d+1)]
    if not admissible(delta, d): ok = False; break
    # perturb out of box -> must break admissibility
    k = random.randint(0, d)
    y2 = list(y)
    lo = -comb(d-1,k) if k <= d-1 else 0
    hi = comb(d-1,k-1) if k >= 1 else 0
    y2[k] = hi + 1 if random.random() < 0.5 else lo - 1
    delta2 = [ballot(d)[t] + 2*y2[t] for t in range(d+1)]
    if admissible(delta2, d): ok = False; break
check("random d<=40: in-box y admissible AND out-of-box y inadmissible (4000+4000)", ok)
flush()

log("")
log("="*78)
log("(C) corners: y = -C(d-1,k) is gamma=-p (Delta=-1); y = +C(d-1,k-1) is gamma=q (Delta=+1)")
log("="*78)
ok = True
for d in range(1, 60):
    ylo = [-comb(d-1,k) if k <= d-1 else 0 for k in range(d+1)]
    yhi = [comb(d-1,k-1) if k >= 1 else 0 for k in range(d+1)]
    if bern_to_poly(ylo, d) != [0,-1]: ok = False           # -p = -x
    if bern_to_poly(yhi, d) != [1,-1]: ok = False           # q = 1-x
    Dlo = [ballot(d)[t] + 2*ylo[t] for t in range(d+1)]
    Dhi = [ballot(d)[t] + 2*yhi[t] for t in range(d+1)]
    if bern_to_poly(Dlo, d) != [-1]: ok = False
    if bern_to_poly(Dhi, d) != [1]: ok = False
check("corners are -p / +q (Delta = -1 / +1), d = 1..59", ok)

log("")
log("="*78)
log("(D) boundary cells: y_0 = gamma(1) in {-1,0}, y_d = gamma(0) in {0,1}")
log("="*78)
ok = True
for d in range(1, 40):
    for _ in range(50):
        y = []
        for k in range(d+1):
            lo = -comb(d-1,k) if k <= d-1 else 0
            hi = comb(d-1,k-1) if k >= 1 else 0
            y.append(random.randint(lo, hi))
        P = bern_to_poly(y, d)
        v1 = sum(P)                      # gamma(1)
        v0 = P[0] if P else 0            # gamma(0)
        if not (v1 == y[0] and -1 <= y[0] <= 0): ok = False
        if not (v0 == y[d] and 0 <= y[d] <= 1): ok = False
check("y_0 = gamma(1) in {-1,0} and y_d = gamma(0) in {0,1} (random)", ok)
flush()

log("")
log("="*78)
log("(E) atom feasibility: eps*p^a*q^b at degree d, cells eps*C(d-a-b,k-b)")
log("="*78)
ok = True; tab = []
for d in range(2, 16):
    for a in range(0, d+1):
        for b in range(0, d+1-a):
            atom = pscale(pshift(qpow(b), a), 1)  # p^a q^b
            cells = poly_to_bern(atom, d)
            want = [comb(d-a-b, k-b) if b <= k <= d-a else 0 for k in range(d+1)]
            if cells != want: ok = False
            fplus = ybox_ok(cells, d)
            fminus = ybox_ok([-v for v in cells], d)
            wplus = (b >= 1)
            wminus = (a >= 1)
            if (fplus != wplus) or (fminus != wminus):
                ok = False
                tab.append((d,a,b,fplus,fminus))
check("cells(p^a q^b) = C(d-a-b,k-b); +feasible iff b>=1, -feasible iff a>=1 (d<=15)",
      ok, str(tab[:5]))
flush()

log("")
log("="*78)
log("(F) witness map: known witnesses -> y-cells, master equation sum x^i gamma_i = G_R")
log("="*78)

def load_blocks(fn):
    d = json.load(open(os.path.join(CP, fn)))
    return d["R"], d["profile"], d["blocks"]

WITNESS_FILES = [
    ("amm12592_floor_witness_R8_ruleB_boxeph.json", None),
    ("amm12592_floor_witness_R16_ruleB_boxeph.json", None),
    ("amm12592_witness_R128_ruleB_D0_0_boxeph.json", None),
    ("amm12592_witness_R256_ruleA_D0_1_boxeph.json", None),
    ("amm12592_witness_R512_ruleA_D0_8_boxeph.json", None),
]
# plus the combined R8/R16/R32 file
combined = json.load(open(os.path.join(CP, "amm12592_floor_witnesses_R8_R16_R32.json")))

def analyze(R, prof, blocks, label):
    # row polynomials once (Horner), epoch identity
    rowP = [bern_to_poly(blocks[i], prof[i]) for i in range(R)]
    S = []
    for i in range(R):
        S = padd(S, pshift(rowP[i], i))
    ok1 = (S == qpow(R-1))
    # admissibility every row
    ok2 = all(admissible(blocks[i], prof[i]) for i in range(R))
    # profile flag: is prof within floor+D0 for some small D0?
    D0 = max(prof[i] - floor_gamma_star(R+i) for i in range(R))
    # ballot normal form? Delta_{R-1} == -1 and rows 0..R-2 -> gamma in box
    nf = (rowP[R-1] == [-1])
    okY = True; okG = None
    if nf:
        Gsum = []
        for i in range(R-1):
            d = prof[i]; b = ballot(d)
            yy = []
            bad = False
            for t in range(d+1):
                diff = blocks[i][t] - b[t]
                if diff % 2 != 0: bad = True; break
                yy.append(diff//2)
            if bad or not ybox_ok(yy, d): okY = False; break
            # gamma_i = (Delta_i - (p-q))/2 as polynomial (cheap, from rowP)
            gpoly = pscale(psub(rowP[i], [-1, 2]), 1)
            assert all(v % 2 == 0 for v in gpoly)
            gpoly = [v//2 for v in gpoly]
            Gsum = padd(Gsum, pshift(gpoly, i))
        if okY:
            okG = (pscale(Gsum, 2) == twoG(R))
    check("%s: identity/admissible/D0=%d/normal-form=%s/y-box=%s/sum=G_R:%s" %
          (label, D0, nf, okY if nf else "-", okG),
          ok1 and ok2 and ((not nf) or (okY and okG)))
    return nf

for w in combined:
    analyze(w["R"], w["profile"], w["blocks"], "combined R=%d" % w["R"])
for fn, _ in WITNESS_FILES:
    try:
        R, prof, blocks = load_blocks(fn)
        analyze(R, prof, blocks, fn.replace("amm12592_", "").replace("_boxeph.json",""))
    except FileNotFoundError:
        log("  [SKIP] %s not found" % fn)
flush()

log("")
log("="*78)
log("SUMMARY: %s" % ("ALL CHECKS PASSED" if FAIL == 0 else "%d FAILURES" % FAIL))
log("="*78)
log("""
REDUCED PROBLEM (proved above, no parity left):
  find integer cell arrays y_i (i = 0..R-2), y_{i,k} in
      [-C(d_i-1,k), C(d_i-1,k-1)],
  with  sum_{i=0}^{R-2} x^i * (sum_k y_{i,k} B_{d_i,k}) = G_R,
  where [x^0]G_R = 1, [x^j]G_R = ((-1)^j C(R-1,j) - 1)/2.
Equivalently gamma_i = q*beta_i - p*alpha_i with alpha_i,beta_i subunital at
degree d_i - 1 (cells in [0, C(d_i-1,k)]).  Corners: gamma_i = -p ... +q.
""")
flush()
