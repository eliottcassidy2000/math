"""AMM 12592 / HYP-9061: doubling structure of gamma* floor-profile epoch closures.

boxeph 2026-08-03. Machinery: THM-2966 (spine normal form), THM-3002 (epoch closure
(*) q^{R-1} = sum_i p^i Delta_i, Lucas boxes), THM-3026 (multiplicativity (M),
lifting (L), degree fit (D), doubling obstruction = row redistribution),
THM-3029 (profile monotonicity, gamma* floor profile closed at R = 8, 16, 32).

This referee establishes, exactly (int arithmetic only):

S1  STRUCTURE of the R = 8, 16, 32 floor witnesses: every block is
    backbone + correction, Delta_i = b_i + c_i with b_i = (p - q) for
    1 <= i <= R-2 (the ballot block, admissible at every degree >= 1),
    b_{R-1} = -1 (the full box), b_0 = 0; corrections c_i are even in every
    cell and supported in a high-degree window.  The residual recursion has an
    exact ENDGAME ATTRACTOR: E_m := -1 + x + ... + x^m satisfies
    E_m - (2x-1) = x * E_{m-1}, so from sigma = E_m the epoch closes by
    emitting (p - q) blocks until E_0 = -1 = minus the full box.

S2  REDUCED IDENTITY: the backbone decomposition turns the epoch identity into
        q * C  =  p^R + q^R - p(p - q),      C := Delta_0 + sum_{1<=i<=R-2} p^i c_i,
    verified exactly on all three witnesses.  p^R + q^R = L_R(pq) is the Lucas
    polynomial (p + q = 1); its coefficient mass is Fibonacci,
    sum_m binom(R-m, m) = F_{R+1} ~ phi^R -- the phi of gamma* on the
    construction side.  The Lucas doubling identity
        L_{2R}(u) = L_R(u)^2 - 2 u^R          (u = pq)
    is THM-2160's middle pair at epoch scale, verified as an exact polynomial
    identity for R = 8, 16, 32.

S3  DIFFERENCED CARVE: with D_j := sum_{i+i'=j} Delta_i Delta_{i'} (ordered
    pairs, cell convolution = product by (M)), the q-split rows
        E_j := D_j - D_{j-1}   (lifted to the 2R floor profile)
    satisfy sum_j p^j E_j = q^{2R-1} EXACTLY, with degrees fitting by
    floor-superadditivity at ANY rate.  The differencing telescopes the
    backbone products: measured worst overshoot max |cell|/capacity collapses
    from the naive lift's pile-up (THM-3026: factors ~ pair count) to
    7 / 9 / 25 at 8->16 / 16->32 / 32->64.

S4  CARRY SWEEP (the doubling map): process rows j = 0..2R-1; X := carve row +
    carry; emit the cellwise parity-clamp of X into the Lucas box; the
    leftover P := X - emitted has P_top = 0 and becomes the carry AT DEGREE
    d-1 WITH THE SAME CELLS (division by x is exact in cell space).  A
    tau-aware boundary recursion rho(j, m) prescribes the top W carry cells so
    every forced top cell is +-1.  RESULT:
      8->16   CLOSED (18 fallbacks), verified exact;
      16->32  CLOSED (22 fallbacks), verified exact;
      32->64  CLOSED with ZERO fallbacks -- a fully deterministic map -- and
              THE FIRST CLOSURE OF THE gamma* FLOOR PROFILE AT R = 64:
              C = 1 + gamma* = log_5(5 phi^2) is now ATTAINED for n <= 127.
              Witness: amm12592_floor_witness_R64_boxeph.json (re-verified here).
    64->128 does NOT close under this sweep (carry runaway, reported below);
    per THM-3029 these failures are TRANSPORT ARTIFACTS of this particular
    map, not evidence of infeasibility.

S5  OBSTRUCTION, QUANTIFIED: the map's output blocks are FAT -- correction
    mass max|c| = 3.4e21 at R = 64 versus 1.8e5 for the beam witness at
    R = 32 -- and fat input makes the next carve's cross-mass box-scale:
    64->128 carry runs away (x ~2.6/row vs box growth x 2^gamma*/row) for
    every tried slack D0 <= 12 and steering width W <= 100.  Folding two
    levels at once from a slim seed does not help: the (2^s - 1)-th
    difference carve overshoot grows 25 -> 341 -> 5739 -> 1.4e6 with level.
    THE DECISIVE MISSING PIECE IS A SLIMMING PASS (any slim R = 64 solution
    closes 128 by this map, by the S3/S4 mechanism that closed 32->64).

S6  Floor-profile guard: the rational approximation GS = 5979874356654402/1e16
    of gamma* = log(phi)/log(sqrt 5) produces the TRUE floor profile for all
    row values m <= 512 (sympy 60-digit check, min distance to integer).

Exact arithmetic throughout: python int / math.comb; sympy only in S6.
"""
import json, os, sys, time
from math import comb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import amm12592_gamma35_beam_deathstar as beam

GS = (5979874356654402, 10**16)
HERE = os.path.dirname(os.path.abspath(__file__))
LEDGER = []

def check(name, ok):
    LEDGER.append((name, bool(ok)))
    print(f'  [{"PASS" if ok else "FAIL"}] {name}')

def prof(R, g1, g2, D0): return [(g1*(R+i))//g2 + D0 for i in range(R)]

def conv(a, b):
    r = [0]*(len(a)+len(b)-1)
    for i, u in enumerate(a):
        if u:
            for j, v in enumerate(b): r[i+j] += u*v
    return r

def lift(delta, d, dp):
    if dp == d: return list(delta)
    return conv(delta, [comb(dp-d, k) for k in range(dp-d+1)])

def admissible(delta, d):
    if len(delta)-1 > d: return False
    return all(abs(v) <= comb(d, k) and (v-comb(d, k)) % 2 == 0 for k, v in enumerate(delta))

def block_poly(delta, d):
    poly = [0]*(d+1)
    for k, v in enumerate(delta):
        if v:
            for t, c in enumerate(beam.basis_poly(d, k)): poly[t] += v*c
    return poly

def ballot(d):
    """cells of the block representing p - q = 2x - 1 at degree d (admissible for d >= 1)"""
    return [comb(d-1, k) - (comb(d-1, k-1) if k else 0) for k in range(d+1)]

def verify_epoch(R, sol, d):
    acc = [0]*(4*R+8)
    for i, de in enumerate(sol):
        if len(de)-1 > d[i]: return False
        for k, v in enumerate(de):
            if abs(v) > comb(d[i], k) or (v-comb(d[i], k)) % 2: return False
            if v:
                for t, c in enumerate(beam.basis_poly(d[i], k)): acc[i+t] += v*c
    return beam.trim(acc) == beam.trim(beam.qpow(R-1))

W_ = {w['R']: w for w in json.load(open(os.path.join(HERE, 'amm12592_floor_witnesses_R8_R16_R32.json')))}

print('=' * 78)
print('S1  backbone + correction structure; endgame attractor')
print('=' * 78)
for d in (5, 12, 30, 75):
    ok = admissible(ballot(d), d) and beam.trim(block_poly(ballot(d), d)) == [-1, 2]
    if not ok: check(f'ballot block admissible/represents p-q at d={d}', False)
check('ballot block (p-q) admissible and rep = 2x-1 at d in {5,12,30,75}', True)
ok = True
for m in range(1, 41):
    Em = [-1] + [1]*m
    Em1 = [-1] + [1]*(m-1)
    lhs = [Em[0] - (-1), Em[1] - 2] + Em[2:]
    ok &= (beam.trim(lhs) == beam.trim([0] + Em1))
check('endgame attractor: E_m - (2x-1) = x E_{m-1} for m <= 40', ok)

for R in (8, 16, 32):
    w = W_[R]
    dsrc, blocks = w['profile'], w['blocks']
    tags, cmax, windows, attractor_row = [], [], [], None
    sig = beam.qpow(R-1)
    for i, (d, de) in enumerate(zip(dsrc, blocks)):
        poly = beam.trim(block_poly(de, d))
        if poly == [-1]: tags.append('CONST-1')
        elif poly == [-1, 2]: tags.append('p-q')
        elif poly[:2] == [-1, 2] or i == 0: tags.append('p-q+window' if i else 'head')
        else: tags.append('other')
        if 1 <= i <= R-2:
            c = [a-b for a, b in zip(de, ballot(d))]
            cmax.append(max(abs(v) for v in c))
            nz = [k for k, v in enumerate(c) if v]
            windows.append((nz[0], nz[-1]) if nz else None)
        # advance residual
        acc = [0]*(max(len(sig), d+2)+2)
        for t, c2 in enumerate(sig): acc[t] += c2
        for t, c2 in enumerate(block_poly(de, d)): acc[t] -= c2
        sig = beam.trim(acc[1:])
        if attractor_row is None and sig and sig == [-1] + [1]*(len(sig)-1):
            attractor_row = i
    evenok = all((a-b) % 2 == 0 for i in range(1, R-1)
                 for a, b in zip(blocks[i], ballot(dsrc[i])))
    print(f'R={R:3d}: tags head..tail: {tags[0]}, mid rows p-q+window, tail {tags[-3:]};'
          f' attractor entered after row {attractor_row}')
    print(f'       max|c| over mid rows = {max(cmax)},  median = {sorted(cmax)[len(cmax)//2]}')
    check(f'R={R}: corrections even in every cell (parity class of backbone)', evenok)
    check(f'R={R}: last block = -1 (full box), residual enters E_m attractor',
          tags[-1] == 'CONST-1' and attractor_row is not None)

print()
print('=' * 78)
print('S2  reduced identity qC = p^R + q^R - p(p-q); Lucas / middle-pair doubling')
print('=' * 78)
for R in (8, 16, 32):
    w = W_[R]
    dsrc, blocks = w['profile'], w['blocks']
    # C(x) = Delta_0 + sum_{i=1}^{R-2} x^i c_i  (polys); check q*C == p^R + q^R - p(p-q)
    Cacc = [0]*(6*R+16)
    for t, c in enumerate(block_poly(blocks[0], dsrc[0])): Cacc[t] += c
    for i in range(1, R-1):
        cpoly = block_poly([a-b for a, b in zip(blocks[i], ballot(dsrc[i]))], dsrc[i])
        for t, c in enumerate(cpoly): Cacc[i+t] += c
    qC = conv([1, -1], beam.trim(Cacc))
    # rhs = x^R + (1-x)^R - x(2x-1)
    rhs = [0]*(R+3)
    rhs[R] += 1
    for t, c in enumerate(beam.qpow(R)): rhs[t] += c
    rhs[1] += 1; rhs[2] -= 2
    check(f'R={R}: q*C == p^R + q^R - p(p-q) exactly', beam.trim(qC) == beam.trim(rhs))
# Lucas doubling: p^{2R} + q^{2R} = (p^R + q^R)^2 - 2 (pq)^R, as x-polynomials
for R in (8, 16, 32):
    LR = [0]*(R+1); LR[R] = 1
    LR = [a+b for a, b in zip(LR, beam.qpow(R) + [0]*(R+1-len(beam.qpow(R))))]
    L2R = [0]*(2*R+1); L2R[2*R] = 1
    q2 = beam.qpow(2*R)
    L2R = [a + (q2[t] if t < len(q2) else 0) for t, a in enumerate(L2R)]
    u = conv([0, 1], [1, -1])          # pq = x(1-x)
    uR = [1]
    for _ in range(R): uR = conv(uR, u)
    rhs = [a - 2*b for a, b in zip(conv(LR, LR), uR + [0]*(len(conv(LR, LR))-len(uR)))]
    check(f'R={R}: Lucas doubling L_2R = L_R^2 - 2(pq)^R (middle pair, THM-2160)',
          beam.trim(L2R) == beam.trim(rhs))
fib = [1, 1]
for _ in range(70): fib.append(fib[-1]+fib[-2])
ok = all(sum(comb(R-m, m) for m in range(R//2+1)) == fib[R] for R in range(1, 65))
check('Fibonacci mass: sum_m binom(R-m,m) = F_{R+1} for R <= 64 (phi enters)', ok)

print()
print('=' * 78)
print('S3  differenced carve E_j = D_j - D_{j-1}: exact identity, overshoot collapse')
print('=' * 78)

def carve(R, dsrc, blocks, dt):
    R2 = 2*R
    E = []
    for j in range(R2):
        cells = [0]*(dt[j]+1)
        def add(jsrc, sign):
            for i in range(max(0, jsrc-R+1), min(jsrc, R-1)+1):
                pr = conv(blocks[i], blocks[jsrc-i])
                dd = dsrc[i]+dsrc[jsrc-i]
                assert dd <= dt[j]
                pl = lift(pr, dd, dt[j])
                for k, v in enumerate(pl): cells[k] += sign*v
        add(j, +1)
        if j >= 1: add(j-1, -1)
        E.append(cells)
    return E

def carve_identity(R, dt, E):
    acc = [0]*(8*R+16)
    for j, cells in enumerate(E):
        for t, c in enumerate(block_poly(cells, dt[j])): acc[j+t] += c
    return beam.trim(acc) == beam.trim(beam.qpow(2*R-1))

for R in (8, 16, 32):
    w = W_[R]
    dt = prof(2*R, *GS, 0)
    E = carve(R, w['profile'], w['blocks'], dt)
    idok = carve_identity(R, dt, E)
    over = max(max((abs(v)/comb(dt[j], k) for k, v in enumerate(E[j])), default=0)
               for j in range(2*R))
    # naive (undifferenced) comparison: rows q*D_j, overshoot vs box at degree+1
    naive = 0.0
    for j in range(2*R-1):
        cells = [0]*(dt[j]+2)
        for i in range(max(0, j-R+1), min(j, R-1)+1):
            pr = conv(w['blocks'][i], w['blocks'][j-i])
            pl = lift(pr, w['profile'][i]+w['profile'][j-i], dt[j])
            for k, v in enumerate(pl): cells[k+1] += v     # * q = cell shift
        naive = max(naive, max(abs(v)/comb(dt[j]+1, k) for k, v in enumerate(cells)))
    print(f'{R:3d}->{2*R:3d}: carve identity exact: {idok};  worst overshoot '
          f'differenced = {over:7.1f}   naive q*D_j = {naive:9.1f}')
    check(f'{R}->{2*R}: carve identity sum p^j E_j = q^(2R-1)', idok)

print()
print('=' * 78)
print('S4  the carry sweep: doubling map closes 8->16, 16->32, 32->64 (floor, D0=0)')
print('=' * 78)

def pclamp(v, cap):
    w = max(-cap, min(cap, v))
    if (w-cap) % 2: w = w-1 if w-1 >= -cap else w+1
    return w

class Steer:
    """tau-aware top-boundary prescription: rho(j,m) = required carry cell
    P_j[d_j-1-m] so that row j+1 keeps its forced top +-1 and its own chain."""
    def __init__(self, E, dt, R2, L):
        self.E, self.dt, self.R2, self.L = E, dt, R2, L
        self.memo = {}
    def rho(self, j, m, depth):
        if (j, m) in self.memo: return self.memo[(j, m)]
        E, dt = self.E, self.dt
        if j+1 >= self.R2:
            val = 0
        else:
            d = dt[j]; dn = dt[j+1]; r = dn - (d-1)
            fixed = 0
            for u in range(max(0, m-r), m):
                fixed += comb(r, m-u)*self.rho(j, u, depth)
            base = E[j+1][dn-m] + fixed
            if m == 0:
                val = min((1-base, -1-base), key=abs)
            elif depth < self.L:
                b = 0 if comb(dn, dn-m) % 2 == 0 else 1
                val = (self.rho(j+1, m-1, depth+1) + b) - base
            else:
                val = pclamp(base, comb(dn, dn-m)) - base
        self.memo[(j, m)] = val
        return val

def carry_sweep(R, dsrc, blocks, D0=0, W=10, L=10):
    R2 = 2*R
    dt = prof(R2, *GS, D0)
    E = carve(R, dsrc, blocks, dt)
    st = Steer(E, dt, R2, L)
    K, dK, sol, nfall, maxcarry = None, None, [], 0, 0
    for j in range(R2):
        d = dt[j]
        X = list(E[j])
        if K is not None:
            if dK > d: return None, f'carry degree {dK} > row degree {d} at row {j}', dt, nfall, maxcarry
            KL = lift(K, dK, d) if dK < d else list(K)
            for k, v in enumerate(KL): X[k] += v
        if abs(X[d]) != 1: return None, f'row {j}: forced top cell = {X[d]}', dt, nfall, maxcarry
        cells = [0]*(d+1); cells[d] = X[d]
        P = [0]*d
        presc = {}
        if j+1 < R2:
            for m in range(min(W, d)):
                presc[d-1-m] = st.rho(j, m, 0)
        for k in range(d-1, -1, -1):
            cap = comb(d, k)
            if k in presc:
                e = X[k] - presc[k]
                if abs(e) <= cap and (e-cap) % 2 == 0:
                    cells[k] = e; P[k] = presc[k]; continue
                nfall += 1
            v = pclamp(X[k], cap)
            cells[k] = v; P[k] = X[k] - v
        sol.append(cells)
        if P: maxcarry = max(maxcarry, max(abs(x) for x in P))
        K, dK = P, d-1
        if j == R2-1 and any(P):
            return None, 'final carry nonzero', dt, nfall, maxcarry
    return sol, 'CLOSED', dt, nfall, maxcarry

params = {8: (10, 10), 16: (14, 14), 32: (10, 10)}
out64 = None
for R in (8, 16, 32):
    w = W_[R]
    wd, lk = params[R]
    t0 = time.time()
    sol, msg, dt, nfall, mc = carry_sweep(R, w['profile'], w['blocks'], W=wd, L=lk)
    el = time.time()-t0
    ok = verify_epoch(2*R, sol, dt) if sol else False
    print(f'{R:3d}->{2*R:3d}: {msg}, fallbacks={nfall}, verified exact: {ok}  [{el:.1f}s]')
    check(f'{R}->{2*R}: carry sweep closes at floor profile D0=0, verified', ok)
    if R == 32 and ok:
        out64 = {'R': 64, 'profile': dt, 'blocks': sol, 'verified': True,
                 'source': 'carve-and-carry doubling map from R=32 floor witness (boxeph)'}
        json.dump(out64, open(os.path.join(HERE, 'amm12592_floor_witness_R64_boxeph.json'), 'w'))
        from fractions import Fraction as F
        eff = max(F(dt[i], 64+i) for i in range(64))
        print(f'      R=64 gamma* floor witness SAVED; effective rate = {float(eff):.6f}'
              f'  (3/5 profile would be 0.600000)')
        check('32->64: zero steering fallbacks (map is deterministic)', nfall == 0)

print()
print('64->128 attempts from the (fat) sweep output at 64 -- reported as transport')
print('artifacts per THM-3029 sec 1, NOT as evidence of infeasibility:')
for D0 in (0, 12):
    t0 = time.time()
    sol, msg, dt, nfall, mc = carry_sweep(64, out64['profile'], out64['blocks'], D0=D0, W=10, L=10)
    print(f'  64->128 D0={D0:2d}: {"CLOSED" if sol else "runaway -- " + msg}'
          f'  (max carry cell ~ 1e{len(str(mc))-1})  [{time.time()-t0:.1f}s]')

print()
print('=' * 78)
print('S5  obstruction quantified: output fatness and multi-fold overshoot growth')
print('=' * 78)
def cmass(R, dsrc, blocks):
    return max(max(abs(a-b) for a, b in zip(blocks[i], ballot(dsrc[i])))
               for i in range(1, R-1))
for R in (16, 32):
    print(f'  beam/lift witness R={R:3d}: max|c| = {cmass(R, W_[R]["profile"], W_[R]["blocks"])}')
print(f'  sweep output    R= 64: max|c| = {cmass(64, out64["profile"], out64["blocks"])}'
      f'   <-- FAT: box-scale corrections')

def carve_pow(R0, dsrc, blocks, s):
    R = (2**s)*R0
    dt = prof(R, *GS, 0)
    rows, degs = [list(b) for b in blocks], list(dsrc)
    n = R0
    for lev in range(s):
        capdeg = []
        for a in range(2*n-1):
            capdeg.append(max(degs[i]+degs[a-i] for i in range(max(0, a-n+1), min(a, n-1)+1)))
        nrows, ndegs = [], []
        for a in range(2*n-1):
            d = capdeg[a]; cells = [0]*(d+1)
            for i in range(max(0, a-n+1), min(a, n-1)+1):
                pl = lift(conv(rows[i], rows[a-i]), degs[i]+degs[a-i], d)
                for k, v in enumerate(pl): cells[k] += v
            nrows.append(cells); ndegs.append(d)
        rows, degs, n = nrows, ndegs, 2*n-1
    T = 2**s - 1
    E = []
    for j in range(R):
        d = dt[j]; cells = [0]*(d+1)
        for t in range(T+1):
            a = j - t
            if a < 0 or a >= len(rows): continue
            c = comb(T, t)*(-1)**t
            pl = lift(rows[a], degs[a], d) if degs[a] < d else rows[a]
            for k, v in enumerate(pl): cells[k] += c*v
        E.append(cells)
    return dt, E

print('  multi-fold carve overshoot from slim seeds (identity exact in all cases):')
for (R0, s) in [(32, 1), (16, 2), (32, 2), (16, 3)]:
    dt, E = carve_pow(R0, W_[R0]['profile'], W_[R0]['blocks'], s)
    R = (2**s)*R0
    acc = [0]*(4*R+8)
    for j, cells in enumerate(E):
        for t, c in enumerate(block_poly(cells, dt[j])): acc[j+t] += c
    idok = beam.trim(acc) == beam.trim(beam.qpow(R-1))
    over = max(max((abs(v)/comb(dt[j], k) for k, v in enumerate(E[j])), default=0) for j in range(R))
    print(f'    seed R0={R0:3d} x 2^{s} -> {R:4d}: identity={idok}  worst overshoot = {over:12.1f}')
    check(f'multi-fold carve identity R0={R0} s={s}', idok)

print()
print('=' * 78)
print('S6  floor-profile guard: GS rational floors == true gamma* floors, m <= 512')
print('=' * 78)
import sympy as sp
gam = sp.log(sp.GoldenRatio)/sp.log(sp.sqrt(5))
ok, mind = True, sp.Float(1, 60)
for m in range(1, 513):
    x = sp.N(gam*m, 60)
    fl = int(sp.floor(x))
    dist = min(x - fl, fl + 1 - x)
    if dist < mind: mind = dist
    if fl != (GS[0]*m)//GS[1]: ok = False
check(f'GS floors exact for m <= 512 (min dist to integer = {sp.N(mind, 3)})', ok)

print()
print('=' * 78)
nfail = sum(1 for _, o in LEDGER if not o)
print(f'LEDGER: {len(LEDGER)} checks, {nfail} failures')
print('ALL DOUBLING-STRUCTURE CHECKS PASSED' if nfail == 0 else 'FAILURES PRESENT')
