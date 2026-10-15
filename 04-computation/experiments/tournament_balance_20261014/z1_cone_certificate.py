"""Certificate for HYP-9279 (Z1) when every strong component of T - z is small (THM-4655 (5) chain).

Normalised block parameters (divide the M_j entries by h(D)): e, beta, gamma, s, sigma, alpha; for v in the block the
'v in A' parts beta(v), gamma(v), alpha(v) (and sigma/2).  Left vector after >= 1 block: (pre, gap, mid, post) = (1, g, P, Lp);
right vector from the next block on: (Rp, s', Q, 1).  For v in an interior block i:
   N = Rp + e s' + beta(v) Q + gamma(v) + P (sigma Q / 2 + alpha(v))
   H = Rp + e s' + beta Q + gamma + g s + P (sigma Q + alpha) + Lp
and (Z1) says 3N - H >= 0 and 2H - 3N >= 0.

Invariant region (all preserved by every block, proved in the theorem text):
   0 <= g <= 1,  P >= Pmin,  cL <= Lp <= P;     0 <= s' <= 1,  Q >= Qmin,  cR <= Rp <= Q,
with Pmin = min beta over first blocks (z has an out-neighbour in D_1), cL = min (gamma + s) over first blocks,
Qmin = min alpha over last blocks (D_k has an in-neighbour of z), cR = min (e + gamma) over last blocks.
The theorem uses Pmin = cL = Qmin = cR = 1 for end blocks of ANY size (THM-4655 (7)(ii): cut a Hamiltonian path of the block
at its first out-neighbour / last in-neighbour of z); the values computed here over the enumerated small blocks are also 1,
so the check is exactly the one the proof needs (audit E part 2).
Both expressions are affine in each of g, Lp, Rp, s' (for fixed P, Q), so the minimum over the region is attained at
extreme values of these; for fixed extremes they are bilinear a PQ + b P + c Q + d in (P, Q) (Lp = P and Rp = Q are allowed
extremes, giving extra P, Q terms), and a PQ + bP + cQ + d >= 0 on [Pmin, inf) x [Qmin, inf) iff a >= 0, a Qmin + b >= 0,
a Pmin + c >= 0 and the value at (Pmin, Qmin) is >= 0.
Usage: python3 z1_cone_certificate.py MAXSIZE   (blocks = all strong tournaments with <= MAXSIZE vertices, all z-patterns)
"""
import sys, subprocess
from fractions import Fraction as Fr
from itertools import product

def strong_blocks(maxsize):
    out = [((), 1)]                                    # singleton
    for s in range(3, maxsize + 1):
        for line in subprocess.run(['gentourng', '-q', '-c', str(s)], capture_output=True, text=True).stdout.split():
            A = [[False] * s for _ in range(s)]; k = 0
            for i in range(s):
                for j in range(i + 1, s):
                    if line[k] == '1': A[i][j] = True
                    else: A[j][i] = True
                    k += 1
            out.append((tuple(map(tuple, A)), s))
    return out

def paths_counts(A, s, zout):
    N = 1 << s
    f = [[0] * s for _ in range(N)]; g = [[0] * s for _ in range(N)]
    for v in range(s): f[1 << v][v] = g[1 << v][v] = 1
    for S in range(1, N):
        if S & (S - 1) == 0: continue
        for v in range(s):
            if S >> v & 1:
                R = S ^ (1 << v)
                f[S][v] = sum(f[R][u] for u in range(s) if R >> u & 1 and A[u][v])
                g[S][v] = sum(g[R][u] for u in range(s) if R >> u & 1 and A[v][u])
    h = [sum(f[S]) if S else 1 for S in range(N)]
    hp = [sum(g[S][v] for v in range(s) if zout[v]) if S else 1 for S in range(N)]
    hm = [sum(f[S][v] for v in range(s) if not zout[v]) if S else 1 for S in range(N)]
    return h, hp, hm

def params(A, s, zout):
    if s == 1:
        A = ((False,),)
    h, hp, hm = paths_counts(A, s, zout); full = (1 << s) - 1; H = h[full]
    P = {}
    P['e'] = Fr(hm[full], H); P['s'] = Fr(hp[full], H)
    P['beta'] = Fr(sum(h[X] * hp[full ^ X] for X in range(full)), H)
    P['gamma'] = Fr(sum(hm[X] * hp[full ^ X] for X in range(1, full)), H)
    P['sigma'] = Fr(sum(h[X] * h[full ^ X] for X in range(full + 1)), H)
    P['alpha'] = Fr(sum(hm[X] * h[full ^ X] for X in range(1, full + 1)), H)
    vparts = []
    for v in range(s):
        vparts.append({'beta': Fr(sum(h[X] * hp[full ^ X] for X in range(full) if X >> v & 1), H),
                       'gamma': Fr(sum(hm[X] * hp[full ^ X] for X in range(1, full) if X >> v & 1), H),
                       'alpha': Fr(sum(hm[X] * h[full ^ X] for X in range(1, full + 1) if X >> v & 1), H)})
    return P, vparts

def bilinear_nonneg(a, b, c, d, Pm, Qm):
    return a >= 0 and a * Qm + b >= 0 and a * Pm + c >= 0 and a * Pm * Qm + b * Pm + c * Qm + d >= 0

def check(maxsize):
    blocks = []
    for A, s in strong_blocks(maxsize):
        for zp in product((False, True), repeat=s):
            P, vp = params(A, s, zp)
            blocks.append((s, A, zp, P, vp))
    first = [b for b in blocks if any(b[2])]           # z beats some vertex of D_1
    last = [b for b in blocks if not all(b[2])]        # some vertex of D_k beats z
    Pmin = min(b[3]['beta'] for b in first); cL = min(b[3]['gamma'] + b[3]['s'] for b in first)
    Qmin = min(b[3]['alpha'] for b in last); cR = min(b[3]['e'] + b[3]['gamma'] for b in last)
    print(f"blocks up to size {maxsize}: {len(blocks)} (type, z-pattern) pairs; Pmin = {Pmin}, cL = {cL}, Qmin = {Qmin}, cR = {cR}")
    fails = []
    for (s, A, zp, P, vp) in blocks:
        e, be, ga, ss, si, al = P['e'], P['beta'], P['gamma'], P['s'], P['sigma'], P['alpha']
        for v, V in enumerate(vp):
            bv, gv, av = V['beta'], V['gamma'], V['alpha']
            # 3N - H = 2Rp + 2 e s' + (3bv - be) Q + 3gv - ga - g s + P (si Q / 2 + 3 av - al) - Lp
            # 2H - 3N = -Rp - e s' + (2 be - 3 bv) Q + 2 ga - 3 gv + 2 g s + P (si Q / 2 + 2 al - 3 av) + 2 Lp
            for g in (Fr(0), Fr(1)):
                for sp in (Fr(0), Fr(1)):
                    for Lp_is_P in (False, True):
                        for Rp_is_Q in (False, True):
                            # lower: coefficients a PQ + b P + c Q + d ; worst Lp = P (coefficient -1) unless Lp = cL better? Lp enters with -1: worst is Lp = P
                            # we enumerate both extremes anyway
                            Lp_b = -1 if Lp_is_P else 0; Lp_d = 0 if Lp_is_P else -cL
                            Rp_c = 2 if Rp_is_Q else 0; Rp_d = 0 if Rp_is_Q else 2 * cR
                            a = si / 2; b = 3 * av - al + Lp_b; c = 3 * bv - be + Rp_c; d = 2 * e * sp + 3 * gv - ga - g * ss + Lp_d + Rp_d
                            if not bilinear_nonneg(a, b, c, d, Pmin, Qmin): fails.append(('lower', s, A, zp, v, g, sp, Lp_is_P, Rp_is_Q))
                            # upper: Lp enters +2 (worst Lp = cL), Rp enters -1 (worst Rp = Q)
                            Lp_b = 2 if Lp_is_P else 0; Lp_d = 0 if Lp_is_P else 2 * cL
                            Rp_c = -1 if Rp_is_Q else 0; Rp_d = 0 if Rp_is_Q else -cR
                            a = si / 2; b = 2 * al - 3 * av + Lp_b; c = 2 * be - 3 * bv + Rp_c; d = -e * sp + 2 * ga - 3 * gv + 2 * g * ss + Lp_d + Rp_d
                            if not bilinear_nonneg(a, b, c, d, Pmin, Qmin): fails.append(('upper', s, A, zp, v, g, sp, Lp_is_P, Rp_is_Q))
    print(f"  failures of the certificate: {len(fails)}")
    for f in fails[:8]: print("   ", f)
    return fails

if __name__ == '__main__':
    check(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
