# Audit: THM-4780 (5) at prime powers WITHOUT modular input.  Over the audit's own GF(q) (audit_gf.py), brute-force
# S_W4(q) and S_(K6-P6)(q) with exact injectivity, then check
#   2^10 TT5(P_q) = q(q-1)[(q-3)(q^2-21q+105) + S_W4]  and  2^15 TT6(P_q) = q(q-1)[(q-3)(q^3-46q^2+682q-3239)
#   + 7(q-15) S_W4 - 8 S_(K6-P6)]
# against the audit's census values (audit_tt_values.out); and the Frobenius-trace readings
#   S_W4 = -4 - (alpha^r + beta^r),  S_(K6-P6) = h_r^2 - q^2 - 7   (PARI a_p, b(p)).   Audit 2026-10-25.
import re, sys, time
import numpy as np
sys.path.insert(0, '.')
from audit_gf import GF

def S_W4_gf(F):
    q = F.q; chi = F.chi
    D = F.sub(np.arange(q)[None, :], np.arange(q)[:, None])     # D[x, y] = y - x
    one = 1
    Y, Z = np.meshgrid(np.arange(q), np.arange(q), indexing='ij')
    tot = 0
    for a in range(q):
        if a in (0, one): continue
        mask = (Y != Z)
        for z0 in (0, one, a):
            mask &= (Y != z0) & (Z != z0)
        fy = chi[Y] * chi[D[one][Y]] * chi[D[a][Y]]
        fz = chi[Z] * chi[D[one][Z]] * chi[D[a][Z]]
        tot += int(chi[a]) * int((fy * fz)[mask].sum())
    return tot

def S_KP_gf(F):
    """K6 - P6, vertices 0..5, path 0-1-2-3-4-5 removed, order orientation; x0 = 0, x1 = 1"""
    q = F.q; chi = F.chi
    D = F.sub(np.arange(q)[None, :], np.arange(q)[:, None])     # D[x, y] = y - x
    one = 1
    X4, X5 = np.meshgrid(np.arange(q), np.arange(q), indexing='ij')
    base = (X4 != X5) & (X4 != 0) & (X4 != one) & (X5 != 0) & (X5 != one)
    # edges at 4: 04, 14, 24 ; at 5: 05, 15, 25, 35 ; 45 is a path edge (absent)
    c4_0 = chi[X4]; c4_1 = chi[D[one][X4]]
    c5_0 = chi[X5]; c5_1 = chi[D[one][X5]]
    tot = 0
    for x2 in range(q):
        if x2 in (0, one): continue
        s2 = int(chi[x2])                                         # edge 02 ; 12 absent
        for x3 in range(q):
            if x3 in (0, one, x2): continue
            s3 = int(chi[x3] * chi[D[one][x3]])                   # edges 03, 13 ; 23 absent
            mask = base & (X4 != x2) & (X4 != x3) & (X5 != x2) & (X5 != x3)
            prod = c4_0 * c4_1 * chi[D[x2][X4]] * c5_0 * c5_1 * chi[D[x2][X5]] * chi[D[x3][X5]]   # 04,14,24 ; 05,15,25,35
            tot += s2 * s3 * int(prod[mask].sum())
    return tot

vals = {}
for l in open('audit_tt_values.out'):
    m = re.match(r'q=(\d+) TT5=(\d+) TT6=(\d+)', l)
    if m: vals[int(m.group(1))] = (int(m.group(2)), int(m.group(3)))
gp = {}
for l in open('modular_coeffs_gp.out'):
    m = re.match(r'p=(\d+) a8=(-?\d+) bre=(-?\d+) bim=(-?\d+)', l)
    if m: gp[int(m.group(1))] = (int(m.group(2)), int(m.group(3)), int(m.group(4)))
def power_sum(s1, prod, r):
    s = [2, s1]
    for _ in range(2, r + 1): s.append(s1 * s[-1] - prod * s[-2])
    return s[r]
out = []
t0 = time.time()
for (p, r, do6) in [(3, 3, True), (3, 5, True), (7, 3, False)]:
    F = GF(p, r); q = F.q
    sw = S_W4_gf(F)
    tt5, tt6 = vals[q]
    ok5 = q * (q - 1) * ((q - 3) * (q * q - 21 * q + 105) + sw) == 1024 * tt5
    Aq = power_sum(gp[p][0], p ** 3, r)
    line = f'q = {p}^{r} = {q}: S_W4 (brute force) = {sw}; THM-4780 (5) TT5: {ok5}; S_W4 = -4 - A_q (Frobenius, A_q = {Aq}): {sw == -4 - Aq}'
    if do6:
        skp = S_KP_gf(F)
        ok6 = q * (q - 1) * ((q - 3) * (q ** 3 - 46 * q * q + 682 * q - 3239) + 7 * (q - 15) * sw - 8 * skp) == 32768 * tt6
        hr = power_sum(gp[p][2], p * p, r)
        line += f'; S_(K6-P6) (brute force) = {skp}; THM-4780 (5) TT6: {ok6}; S_(K6-P6) = h_r^2 - q^2 - 7: {skp == hr * hr - q * q - 7}'
    out.append(line + f'  [{time.time() - t0:.0f}s]')
    print(out[-1], flush=True)
open('audit_primepower_sums.out', 'w').write('\n'.join(out) + '\n')
