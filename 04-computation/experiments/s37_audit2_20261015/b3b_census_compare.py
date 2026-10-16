#!/usr/bin/env python3
"""S37 audit 2, script b3b: assemble the independent census of b3 (runs over [1, 2^36) and [2^36, 2^37)) with exact
Python integers and compare it with the session's census extension and block moments.

Checks:
 * class sums modulo 243, H(N_j) and the exact prefix sums sum_(n<=N_j) H(n), sum_(n<=N_j) H(n)^2 at every grid point
   N_j = round(2^(10 + j/64)) < 2^37 against 05-knowledge/results/collatz_census_extension_2p37_20261015.npz
   (exact decimal strings preH1_str, preH2_str);
 * exact octave sums of H and H^2 for k = 20..36 against collatz_census_extension_2p37_20261015.json;
 * all 2^20 block moments (sum t, sum (n - jB) t), B = 2^17, against collatz_boundary_moments_2p37_20261015.npz,
   and S(2^37), H(2^37) (inclusive of n = 2^37, t = -1) against collatz_boundary_moments_2p37_20261015.json;
 * exact dyadic energies e_X = E_X/X^2 and octave means m_X = sum_(X<=n<2X) H(n)/(X sqrt X) for every k = 0..36
   (k < 20 from a direct exact census), the note's section-4 table, the window energies for starts in (2^33, 2^36].
Usage: python b3b_census_compare.py A.npz B.npz census_ext.npz census_ext.json moments37.npz moments37.json
"""
import sys, json, math
import numpy as np

A = np.load(sys.argv[1]); Bz = np.load(sys.argv[2])
CEN = np.load(sys.argv[3]); CJ = json.load(open(sys.argv[4]))
MOMZ = np.load(sys.argv[5]); MJ = json.load(open(sys.argv[6]))
BLK = 1 << 20
NB = (1 << 37) // BLK
LOC = np.concatenate([A['LOC'], Bz['LOC']]).astype(np.int64)
HEND = np.concatenate([A['HEND'], Bz['HEND']])
SH1 = np.concatenate([A['SH1'], Bz['SH1']])
SH2 = np.concatenate([A['SH2'], Bz['SH2']])
MOM = np.concatenate([A['MOM'], Bz['MOM']])
assert LOC.shape[0] == NB and int(A['b1']) == int(Bz['b0'])
grid = A['grid']
assert np.array_equal(grid, Bz['grid']) and np.array_equal(grid, CEN['grid'])
# my own grid construction must match the session's
raw = sorted(set(int(round(2 ** (10 + j / 64))) for j in range(0, 27 * 64 + 1)))
mygrid = np.array([g for g in raw if g < (1 << 37)], dtype=np.int64)
print(f'grid: {len(grid)} points; identical to an independent construction: {np.array_equal(mygrid, grid)}')
gselA = A['gsel']; gselB = Bz['gsel']
GLOC = np.where(gselA[:, None], A['GLOC'], Bz['GLOC'])
GH = np.where(gselA, A['GH'], Bz['GH'])
GSH1 = np.where(gselA, A['GSH1'], Bz['GSH1'])
GSH2 = np.where(gselA, A['GSH2'], Bz['GSH2'])
assert (gselA ^ gselB).all()

# exact prefix over blocks
pref = np.zeros((NB + 1, 243), dtype=np.int64)
pref[1:] = np.cumsum(LOC, axis=0)
H0 = [0] * (NB + 1)
cum1 = [0] * (NB + 1)
cum2 = [0] * (NB + 1)
h = 0
for b in range(NB):
    H0[b] = h
    size = BLK if b > 0 else BLK - 1
    s1 = int(SH1[b]); s2 = int(SH2[b])
    cum1[b + 1] = cum1[b] + size * h + s1
    cum2[b + 1] = cum2[b] + size * h * h + 2 * h * s1 + s2
    h += int(HEND[b])
H0[NB] = h
# H at block starts is also 3 C_2 - S from the class prefix
chk = all(H0[b] == 3 * int(pref[b, 2::3].sum()) - int(pref[b].sum()) for b in range(0, NB, 997))
print(f'H at block starts consistent with class prefix sums: {chk}')

bidx = grid // BLK
mycls = np.zeros((len(grid), 243), dtype=np.int64)
myH = np.zeros(len(grid), dtype=object)
myP1 = []; myP2 = []
for i, g in enumerate(grid.tolist()):
    b = g // BLK
    mycls[i] = pref[b] + GLOC[i]
    hb = H0[b]
    myH[i] = hb + int(GH[i])
    cnt = g - max(b * BLK, 1) + 1
    myP1.append(cum1[b] + cnt * hb + int(GSH1[i]))
    myP2.append(cum2[b] + cnt * hb * hb + 2 * hb * int(GSH1[i]) + int(GSH2[i]))
same_cls = np.array_equal(mycls, CEN['classsums'])
same_H = all(int(myH[i]) == int(CEN['H'][i]) for i in range(len(grid)))
p1s = [int(x) for x in CEN['preH1_str']]; p2s = [int(x) for x in CEN['preH2_str']]
same_p1 = all(a == b for a, b in zip(myP1, p1s)); same_p2 = all(a == b for a, b in zip(myP2, p2s))
nbad_cls = int((mycls != CEN['classsums']).any(axis=1).sum())
print(f'grid comparison with the session census (all {len(grid)} points): class sums mod 243 identical {same_cls} '
      f'(rows differing: {nbad_cls}); H identical {same_H}; prefix sums of H identical {same_p1}; of H^2 identical {same_p2}')
S_at = mycls.sum(axis=1)
i34 = int(np.nonzero(grid == 1 << 34)[0][0])
print(f'S(2^34) = {int(S_at[i34])}, H(2^34) = {int(myH[i34])}')

# octave sums k = 20..36
mismatch = []
oct_my = {}
for k in range(20, 37):
    b0, b1 = (1 << k) // BLK, (1 << (k + 1)) // BLK
    s1 = cum1[b1] - cum1[b0]; s2 = cum2[b1] - cum2[b0]
    oct_my[k] = (s1, s2)
    ref = CJ['octave_sums_H_H2'][str(k)]
    if (str(s1), str(s2)) != (ref[0], ref[1]):
        mismatch.append(k)
print(f'exact octave sums of H and H^2, k = 20..36: mismatches with the session JSON: {mismatch}')

# block moments
mom_my = MOM.reshape(-1, 2)
mom_ref = MOMZ['moments']
print(f'block moments: {mom_my.shape[0]} rows; identical to the session moments file: {np.array_equal(mom_my, mom_ref)}'
      f' (rows differing: {int((mom_my != mom_ref).any(axis=1).sum())})')
N37 = 1 << 37
S37 = int(pref[NB].sum()) - 1        # n = 2^37: tau = 37, t = -1
C2 = int(pref[NB, 2::3].sum()) + (-1 if N37 % 3 == 2 else 0)
H37 = 3 * C2 - S37
print(f'S(2^37) = {S37} (session {MJ["S"]}), H(2^37) = {H37} (session {MJ["H"]}); 2^37 mod 3 = {N37 % 3}')

# small octaves from a direct census
def tau_par(n):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        s += 1
    return s & 1
LIM = 1 << 21
tt = np.zeros(LIM + 1, dtype=np.int64)
tt[1] = 1
for n in range(2, LIM + 1):
    if n % 2 == 0:
        tt[n] = -tt[n // 2]
    else:
        x, s = n, 0
        while x >= n:
            x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
            s += 1
        tt[n] = tt[x] * (1 - 2 * (s & 1))
w = np.where(np.arange(LIM + 1) % 3 == 2, 2, -1); w[0] = 0
Hs = np.cumsum(w * tt)
print('exact dyadic energies e_X = E_X/X^2 and octave means m_X (k = 0..36):')
emax = (0, None)
rows = []
for k in range(0, 37):
    X = 1 << k
    if k < 20:
        seg = Hs[X:2 * X].astype(object)
        s1 = int(sum(int(v) for v in seg)); s2 = int(sum(int(v) ** 2 for v in seg))
        if k >= 20:
            pass
    else:
        s1, s2 = oct_my[k]
    e = s2 / X ** 2
    mX = s1 / (X * math.sqrt(X))
    rows.append((k, e, mX))
    if k >= 10 and e > emax[0]:
        emax = (e, k)
print('  ' + '  '.join(f'{k}:{e:.4f}/{m:+.3f}' for k, e, m in rows))
big = [(k, round(e, 3)) for k, e, m in rows if e > 2]
print(f'  windows with e_X > 2: {big}; max over 10 <= k <= 36: {emax[0]:.6f} at k = {emax[1]}; '
      f'max over k = 33..36: {max(e for k, e, m in rows if k >= 33):.4f}')
# cross-check k = 20 from grid prefix sums
g = grid.tolist()
def E_grid(k):
    i0 = g.index(1 << k); i1 = g.index(1 << (k + 1))
    return (myP2[i1] - int(myH[i1]) ** 2) - (myP2[i0] - int(myH[i0]) ** 2)
print(f'  octave sums from the grid prefix agree with block sums for k = 20..35: '
      f'{all(E_grid(k) == oct_my[k][1] for k in range(20, 36))}')
# window energies, starts in (2^33, 2^36]
N = grid.astype(np.float64)
p2 = np.array([float(x) for x in myP2]); Hf = np.array([float(x) for x in myH])
jj = np.arange(len(grid) - 64)
E = (p2[jj + 64] - Hf[jj + 64] ** 2) - (p2[jj] - Hf[jj] ** 2)
e_w = E / N[jj] ** 2
uY = np.log2(N[jj])
sel = uY > 33 + 1e-9
print(f'window energies for starts in (2^33, 2^36]: n = {int(sel.sum())}, mean {e_w[sel].mean():.4f}, max {e_w[sel].max():.4f}; '
      f'in sample (2^14 <= start, end <= 2^34): mean {e_w[(uY >= 14) & (uY <= 33 + 1e-9)].mean():.4f}, '
      f'max {e_w[(uY >= 14) & (uY <= 33 + 1e-9)].max():.4f}')
