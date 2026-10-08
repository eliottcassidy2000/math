#!/usr/bin/env python3
"""audit_F (C1): congruence obstructions -- contracting, 'non-degenerate' maps in which y and y+1 NEVER merge.

Lemma (elementary).  Let T(x) = (m_i x + r_i)/d on x = i mod d, and let p be a prime with p not dividing d * prod m_i.
Suppose there is c in Z_p with (m_i - d) c + r_i = 0 mod p for every i (for instance c = 0 when p | r_i for all i).
Then along the pair chain u_n = M_n v_n + e_n started at (M, e) = (1, e_0),
        e_n - (1 - M_n) c  =  P_n (e_0)   (mod p),   P_n = prod_k m_{i_k}/d^n  (a p-adic unit),
so at every time with M_n = 1 we have e_n = P_n e_0 != 0 mod p whenever p does not divide e_0: the state (1, 0) is never
reached, and a meeting with M_n != 1 is a Haar-null event.  Hence P(y and y + e_0 ever meet at equal times) = 0.
(Proof: e' - (1 - M')c = [m_i e + r_i - (m_i/m_j) M r_j]/d - (1 - M m_i/m_j) c; substitute r_i = (d - m_i) c,
 r_j = (d - m_j) c mod p and simplify to (m_i/d)(e - (1 - M) c).)

This script
 (1) verifies the identity e_n - (1 - M_n) c = P_n e_0 mod p EXACTLY along random pair-chain paths (exact rationals),
     for the maps below;
 (2) confirms by direct big-integer orbits (bigint_pairs.run_sample) that no merge occurs;
 (3) runs the same multipliers with an unobstructed sheet, which merges like the HYP predicts.
Maps:
  3x+5 on Z_2             m = (1, 3),    r = (0, 5)      rank 1, Lambda < 0, p = 5, c = 0   (conjugate to 3x+1 with offset 1/5)
  Z_3 (1,1,5)  r=(0,2,2)  rank 1, Lambda < 0, p = 2, c = 0   (T preserves parity in Z_(2))
  Z_3 (1,2,5)  r=(0,7,14) rank 2, Lambda < 0, p = 7, c = 0   (same multipliers as the HYP's rank-2 example)
  Z_3 (1,2,5)  r=(9,1,5)  rank 2, Lambda < 0, p = 7, c = 1   (affine invariant, r_0 != 0)
  Z_5 (1,2,3,7,1) r=(0,3,4,4,1) + 5*(...) chosen = 0 mod 11:  rank 3 (control: obstruction is rank-independent)
Unobstructed controls: 3x+1; Z_3 (1,1,5) r = (0,-1,-1); Z_3 (1,2,5) r = (0,1,2).
"""
import sys, random, math
from fractions import Fraction as Fr
sys.path.insert(0, '.')
from bigint_pairs import Map, run_sample


def crt_r(d, m, i, p, c=0):
    """smallest nonnegative r with r = -m i mod d and r = -(m - d) c mod p"""
    for r in range(d * p):
        if (m * i + r) % d == 0 and ((m - d) * c + r) % p == 0:
            return r
    raise ValueError


def chain_identity_check(d, br, p, c, e0, nsteps, npaths, rnd):
    """exact pair chain from (1, e0); check e_n - (1 - M_n) c == P_n e0 (mod p) at every step"""
    def modp(q):
        return (q.numerator * pow(q.denominator, -1, p)) % p
    def modd(q):
        return (q.numerator * pow(q.denominator, -1, d)) % d
    bad = 0; hits10 = 0; tot = 0
    for _ in range(npaths):
        M, e, P = Fr(1), Fr(e0), Fr(1)
        for n in range(nsteps):
            j = rnd.randrange(d)
            i = (modd(M) * j + modd(e)) % d
            mi, ri = br[i]; mj, rj = br[j]
            e = (mi * e + ri - Fr(mi, mj) * M * rj) / d
            M = M * Fr(mi, mj)
            P = P * Fr(mi, d)
            tot += 1
            if (modp(e - (1 - M) * c) - modp(P * e0)) % p != 0:
                bad += 1
            if M == 1 and e == 0:
                hits10 += 1
                break
    return bad, hits10, tot


MAPS = [
    ('3x+5 (p=5,c=0)', 2, [(1, 0), (3, 5)], 5, 0),
    ('Z3 (1,1,5) r=(0,2,2) (p=2,c=0)', 3, [(1, 0), (1, 2), (5, 2)], 2, 0),
    ('Z3 (1,2,5) r=(0,7,14) (p=7,c=0)', 3, [(1, 0), (2, 7), (5, 14)], 7, 0),
    ('Z3 (1,2,5) r=(9,1,5) (p=7,c=1)', 3, [(1, 9), (2, 1), (5, 5)], 7, 1),
]
# a rank-3 obstructed map: multipliers (1,2,3,7,1) on Z_5, r_i = -m_i i mod 5 and = 0 mod 11
_m = [1, 2, 3, 7, 1]
MAPS.append(('Z5 (1,2,3,7,1) r=0 mod 11 (p=11,c=0)', 5, [(m, crt_r(5, m, i, 11)) for i, m in enumerate(_m)], 11, 0))
CONTROLS = [
    ('3x+1 (control)', 2, [(1, 0), (3, 1)]),
    ('Z3 (1,1,5) r=(0,-1,-1) (control)', 3, [(1, 0), (1, -1), (5, -1)]),
    ('Z3 (1,2,5) r=(0,1,2) (control = HYP map)', 3, [(1, 0), (2, 1), (5, 2)]),
    ('Z5 (1,2,3,7,1) r=(0,3,4,4,1) (control = HYP rank-3 map)', 5, [(1, 0), (2, 3), (3, 4), (7, 4), (1, 1)]),
]

if __name__ == '__main__':
    NS = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    TM = int(sys.argv[2]) if len(sys.argv) > 2 else 4096
    rnd = random.Random(2026)
    print(f"(1) exact identity e_n - (1-M_n)c = P_n e_0 (mod p) along 200 random chain paths of 300 steps, start (1, 1):")
    for name, d, br, p, c in MAPS:
        mp = Map(d, br)
        assert all(((m - d) * c + r) % p == 0 for m, r in br), name
        bad, hits, tot = chain_identity_check(d, br, p, c, 1, 300, 200, rnd)
        print(f"   {name:42s} m={mp.ms} r={mp.rs} Lambda={mp.lam:+.3f}: {tot} steps, identity violations {bad}, absorptions {hits}")
    print(f"(2)/(3) direct big-integer orbits of y, y+1 (y uniform mod d^(T+64)), {NS} samples, T = {TM}:")
    for item in MAPS + CONTROLS:
        name, d, br = item[0], item[1], item[2]
        mp = Map(d, br)
        res = [run_sample(mp, TM, TM + 64, rnd, 128, False, 1) for _ in range(NS)]
        met = [r for r in res if r[0] is not None]
        nonzero = sum(1 for r in met if any(r[1]))
        q = 1 - len(met) / NS
        qs = {T: sum(1 for r in res if r[0] is None or r[0] > T) / NS for T in (16, 256, TM)}
        print(f"   {name:58s} met {len(met):4d}/{NS} (nonzero-debt meetings {nonzero}); q(16)={qs[16]:.3f} q(256)={qs[256]:.3f} q({TM})={qs[TM]:.3f}", flush=True)
