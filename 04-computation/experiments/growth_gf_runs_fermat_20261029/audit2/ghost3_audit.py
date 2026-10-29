"""Audit of THM-4846 (D') (mod-3 ghost lemma), independent of the producer's runs.py / ghost3_lemma.py.

Conventions (derived here, then CHECKED against the data):
  ordered triple (b, g) = (b[g(0)], b[g(1)], b[g(2)]), b a sorted 3-set, g in S_3 as an image tuple;
  permuting run positions by s sends (b, g) to (b, g o s); M3 commutes with it (checked).
  rho(g)_ij = [g(j) = i] - [g(2) = i]  (action on the augmentation basis e0-e2, e1-e2);
  chart (per orbit) c -> (u, v, W) = (sum c_g, sum sgn(g) c_g, sum c_g rho(g)), W = [[a,b],[c,d]].
  With k_ij(g) = M3[(i,id),(j,g)]:  P_ij = sum k, S_ij = sum sgn k, L_ij = sum k_ij(g) rho(g^-1) and
  N := F M F^-1 must act by u -> P u, v -> S v, W_i -> sum_j L_ij W_j (LEFT multiplication, columns separately).
Checks:
  (a) F from the definition equals THREE_RUN_GLUE (R2); det 54, Smith 1,1,1,3,3,6; H = 6F^-1 = (R4);
      image = (R3) lattice; on (Z/3)^6: (R3) mod 3 <=> Q^-1 W Q = [[u,*],[0,v]] mod 3.
  (b) for every class: N has exactly the claimed block form (verifies left multiplication / coordinates);
      L' = Q^-1 L Q has L21 = 0, L11 = P, L22 = S mod 3; the matrix of N on the claimed Z_(3)-basis of the image,
      ordered (u, y21, y11 | v, w12, y22), reduces mod 3 to diag(I, II) EXACTLY (entrywise).
  (c) rank_F3 M = rank I + rank II; rank_F3 N = rk P + rk S + 2 rk[[P,L12],[0,S]]; spur_3 formula; and
      spur_3 (rank version) = excess count of invariant factors divisible by 3 (Smith forms, FLINT) for n <= 6.
  (d) when S is invertible mod 3: spur_3 = rank(E : ker P -> coker P), E = D1 - L12 S^-1 D2; T5 specifics.
"""
import sys
import itertools
from fractions import Fraction
import numpy as np
from alib import (gentourng, adj_from_ascii, circulant, path_table, m3_matrix, rank_mod_p, nullspace_mod_p,
                  inv_mod_p, snf_diag, snf_summary, fmt_snf, S3)

p = 3


def compose(g, s):          # (g o s)(k) = g(s(k))
    return tuple(g[s[k]] for k in range(3))


def inverse(g):
    out = [0, 0, 0]
    for k in range(3):
        out[g[k]] = k
    return tuple(out)


def sign(g):
    inv = sum(1 for i in range(3) for j in range(i + 1, 3) if g[i] > g[j])
    return -1 if inv % 2 else 1


def rho(g):
    return np.array([[(1 if g[j] == i else 0) - (1 if g[2] == i else 0) for j in range(2)] for i in range(2)],
                    dtype=np.int64)


# ---------------- (a) the chart ----------------
F = np.zeros((6, 6), dtype=np.int64)
for col, g in enumerate(S3):
    R = rho(g)
    F[:, col] = [1, sign(g), R[0, 0], R[0, 1], R[1, 0], R[1, 1]]
F_R2 = np.array([[1, 1, 1, 1, 1, 1], [1, -1, -1, 1, 1, -1], [1, 1, 0, -1, 0, -1], [0, 0, 1, -1, 1, -1],
                 [0, -1, 1, 1, -1, 0], [1, -1, 0, 0, -1, 1]], dtype=np.int64)
H_R4 = np.array([[1, 1, 2, 0, 0, 2], [1, -1, 2, -2, 0, -2], [1, -1, 0, 2, 2, 0], [1, 1, 0, -2, 2, -2],
                 [1, 1, -2, 2, -2, 0], [1, -1, -2, 0, -2, 2]], dtype=np.int64)
print('(a) F from definition == (R2):', bool((F == F_R2).all()))
for g in S3:
    for h in S3:
        assert (rho(compose(g, h)) == rho(g) @ rho(h)).all()
print('    rho is a homomorphism for composition (g o h)(k) = g(h(k)): True')
import flint
fF = flint.fmpz_mat(F.tolist())
print('    det F =', int(fF.det()), ' Smith F =', snf_diag(F))
print('    H = 6 F^-1 == (R4):', bool((H_R4 @ F == 6 * np.eye(6, dtype=np.int64)).all()))


def r3(x):
    u, v, a, b, c, d = x
    return (u - v) % 2 == 0 and (a + b - u) % 3 == 0 and (c + d - u) % 3 == 0 and (a - c - v) % 3 == 0


assert all(r3(F[:, k]) for k in range(6))
# the (R3) lattice has index 2*3^3 = 54 = |det F| (conditions independent: count solutions in (Z/6)^6)
cnt = sum(1 for x in itertools.product(range(6), repeat=6) if r3(x))
print('    columns of F satisfy (R3); |(R3) lattice mod 6| =', cnt, ' => index', 6 ** 6 // cnt,
      '(= |det F| = 54, so image == (R3) lattice)')
Q = np.array([[1, 1], [1, 0]], dtype=np.int64)
Qi = np.array([[0, 1], [1, -1]], dtype=np.int64)
assert (Q @ Qi == np.eye(2, dtype=np.int64)).all()
ok = True
for x in itertools.product(range(3), repeat=6):
    u, v, a, b, c, d = x
    Wp = Qi @ np.array([[a, b], [c, d]]) @ Q
    flag = (Wp[1, 0] % 3 == 0) and ((Wp[0, 0] - u) % 3 == 0) and ((Wp[1, 1] - v) % 3 == 0)
    r3mod3 = (a + b - u) % 3 == 0 and (c + d - u) % 3 == 0 and (a - c - v) % 3 == 0
    ok &= (flag == r3mod3)
print('    on all of (Z/3)^6: (R3) mod 3  <=>  Q^-1 W Q = [[u,*],[0,v]] mod 3 :', ok)
print('    symbolic: Q^-1 W Q = [[c+d, c],[a+b-c-d, a-c]]  ->  c+d=u, a+b=c+d, a-c=v (mod 3) == (R3) mod 3')


# ---------------- per-class machinery ----------------
def big(Mblock, r):
    out = np.zeros((6 * r, 6 * r), dtype=np.int64)
    for o in range(r):
        out[6 * o:6 * o + 6, 6 * o:6 * o + 6] = Mblock
    return out


gidx = {g: k for k, g in enumerate(S3)}
ID = (0, 1, 2)


def analyse(A, want_snf=False, label=''):
    n = A.shape[0]
    M, bases, tup = m3_matrix(A)
    r = len(bases)
    # run symmetry check
    for s in S3:
        perm = [6 * (k // 6) + gidx[compose(S3[k % 6], s)] for k in range(6 * r)]
        assert (M[np.ix_(perm, perm)] == M).all(), 'run symmetry fails'
    # k_ij(g), sectors from the definition (not from N)
    K = {}
    P = np.zeros((r, r), dtype=np.int64)
    S = np.zeros((r, r), dtype=np.int64)
    L = np.zeros((r, r, 2, 2), dtype=np.int64)
    for i in range(r):
        for j in range(r):
            for g in S3:
                kv = M[6 * i + gidx[ID], 6 * j + gidx[g]]
                P[i, j] += kv
                S[i, j] += sign(g) * kv
                L[i, j] += kv * rho(inverse(g))
    # N = F M F^-1, exact
    Fb = big(F, r)
    Hb = big(H_R4, r)
    X = Fb @ M @ Hb
    assert (X % 6 == 0).all(), 'F M F^-1 not integral'
    N = X // 6
    # expected block form: coordinates per orbit 0=u, 1=v, 2=a(W00), 3=b(W01), 4=c(W10), 5=d(W11)
    Nexp = np.zeros_like(N)
    for i in range(r):
        for j in range(r):
            Nexp[6 * i, 6 * j] = P[i, j]
            Nexp[6 * i + 1, 6 * j + 1] = S[i, j]
            for rr in range(2):
                for rp in range(2):
                    for cc in range(2):
                        Nexp[6 * i + 2 + 2 * rr + cc, 6 * j + 2 + 2 * rp + cc] = L[i, j, rr, rp]
    left_ok = bool((N == Nexp).all())
    # right-multiplication alternative (W_i -> sum_j W_j L_ij^T style) must FAIL unless degenerate
    Nalt = np.zeros_like(N)
    for i in range(r):
        for j in range(r):
            Nalt[6 * i, 6 * j] = P[i, j]
            Nalt[6 * i + 1, 6 * j + 1] = S[i, j]
            for rr in range(2):
                for cc in range(2):
                    for cp in range(2):
                        Nalt[6 * i + 2 + 2 * rr + cc, 6 * j + 2 + 2 * rr + cp] = L[i, j, cp, cc]
    right_ok = bool((N == Nalt).all())
    # adapted basis
    Lp = np.einsum('xy,ijyz,zw->ijxw', Qi, L, Q)
    L11, L12, L21, L22 = Lp[:, :, 0, 0], Lp[:, :, 0, 1], Lp[:, :, 1, 0], Lp[:, :, 1, 1]
    cong = bool(((L11 - P) % 3 == 0).all() and (L21 % 3 == 0).all() and ((L22 - S) % 3 == 0).all())
    D1, D2, D3 = (L11 - P) // 3, L21 // 3, (L22 - S) // 3
    Z = np.zeros((r, r), dtype=np.int64)
    I_ = np.block([[P, Z, Z], [D2, S, Z], [D1, L12, P]])
    II = np.block([[S, Z, Z], [L12, P, Z], [D3, D2, S]])
    # exact matrix of N on the Z_(3)-basis of the image, ordered (u*, y21*, y11*, v*, w12*, y22*)
    # basis vectors in (u, v, w11, w12, w21, w22) per orbit, then W = Q W' Q^-1
    def vec(o, kind):
        x = np.zeros(6 * r, dtype=np.int64)
        u = v = 0
        Wp = np.zeros((2, 2), dtype=np.int64)
        if kind == 'u':
            u, Wp[0, 0] = 1, 1
        elif kind == 'v':
            v, Wp[1, 1] = 1, 1
        elif kind == 'y11':
            Wp[0, 0] = 3
        elif kind == 'w12':
            Wp[0, 1] = 1
        elif kind == 'y21':
            Wp[1, 0] = 3
        elif kind == 'y22':
            Wp[1, 1] = 3
        W = Q @ Wp @ Qi
        x[6 * o:6 * o + 6] = [u, v, W[0, 0], W[0, 1], W[1, 0], W[1, 1]]
        return x
    order = ['u', 'y21', 'y11', 'v', 'w12', 'y22']
    B = np.array([vec(o, k) for k in order for o in range(r)], dtype=np.int64).T   # columns
    assert all(r3(B[6 * o:6 * o + 6, c]) or True for o in range(r) for c in range(B.shape[1]))
    # every basis vector satisfies the mod-3 conditions orbitwise
    for c in range(B.shape[1]):
        for o in range(r):
            u, v, a, b, cc, d = B[6 * o:6 * o + 6, c]
            assert (a + b - u) % 3 == 0 and (cc + d - u) % 3 == 0 and (a - cc - v) % 3 == 0
    fB = flint.fmpq_mat([[int(x) for x in row] for row in B])
    fN = flint.fmpq_mat([[int(x) for x in row] for row in N])
    Mb = fB.inv() * fN * fB
    Mb_int = np.zeros((6 * r, 6 * r), dtype=object)
    three_integral = True
    for a_ in range(6 * r):
        for b_ in range(6 * r):
            q = Mb[a_, b_]
            num, den = int(q.p), int(q.q)
            if den % 3 == 0:
                three_integral = False
            Mb_int[a_, b_] = (num * pow(den, -1, 3)) % 3
    Mb3 = Mb_int.astype(np.int64)
    blockdiag = np.zeros((6 * r, 6 * r), dtype=np.int64)
    blockdiag[:3 * r, :3 * r] = I_ % 3
    blockdiag[3 * r:, 3 * r:] = II % 3
    reduction_exact = three_integral and bool((Mb3 == blockdiag).all())
    rM, rN = rank_mod_p(M, 3), rank_mod_p(N, 3)
    rI, rII = rank_mod_p(I_, 3), rank_mod_p(II, 3)
    rP, rS = rank_mod_p(P, 3), rank_mod_p(S, 3)
    rPL = rank_mod_p(np.block([[P, L12], [Z, S]]), 3)
    spur3 = rM - rN
    formula = rI + rII - rP - rS - 2 * rPL
    res = dict(n=n, r=r, left_ok=left_ok, right_ok=right_ok, cong=cong, reduction_exact=reduction_exact,
               rankM_eq=(rM == rI + rII), rankN_eq=(rN == rP + rS + 2 * rPL), spur3=spur3, formula=formula,
               rP=rP, rS=rS)
    # (d) S invertible mod 3
    if rS == r:
        Si = inv_mod_p(S, 3)
        E = (D1 - L12 @ Si @ D2) % 3
        Kp = nullspace_mod_p(P, 3)          # rows: basis of ker P
        if Kp.shape[0]:
            img = (E @ Kp.T) % 3
            conn = rank_mod_p(np.concatenate([P % 3, img], axis=1), 3) - rP
        else:
            conn = 0
        res['E_conn'] = conn
    if want_snf:
        dM = snf_diag(M)
        dN = snf_diag(N)
        div3 = lambda d: sum(1 for x in d if x % 3 == 0)
        res['snf_spur3'] = div3(dN) - div3(dM)
        res['snfM'] = fmt_snf(snf_summary(dM))
        res['snfN'] = fmt_snf(snf_summary(dN))
    res.update(P=P, S=S, L=L, L12=L12, D1=D1, D2=D2, M=M, N=N, bases=bases)
    return res


# ---------------- T5 ----------------
print('\n(d) T5 (arcs i -> i+1, i+2 mod 5)')
A5 = circulant(5, {1, 2})
t5 = analyse(A5, want_snf=True)
print('    Smith M3 =', t5['snfM'], ' (R7 says 1^52 3^2 6^2 0^4)')
print('    Smith N  =', t5['snfN'], ' (R7 unglued says 1^48 3^2 6^6 0^4)')
print('    sectors: P', fmt_snf(snf_summary(snf_diag(t5['P']))), '| S', fmt_snf(snf_summary(snf_diag(t5['S']))),
      '| L (one column of W)', fmt_snf(snf_summary(snf_diag(t5['L'].transpose(0, 2, 1, 3).reshape(20, 20)))))
print('    left-mult block form', t5['left_ok'], '| right-mult alternative', t5['right_ok'], '| congruences', t5['cong'],
      '| exact mod-3 reduction to diag(I,II)', t5['reduction_exact'])
print('    rank_F3: P=%d S=%d (S invertible mod 3: %s, det S = %d)' % (t5['rP'], t5['rS'], t5['rS'] == 10,
      int(flint.fmpz_mat(t5['S'].tolist()).det())))
print('    spur3 (ranks) = %d, formula = %d, Smith-count spur3 = %d, connecting rank of E on ker P = %d'
      % (t5['spur3'], t5['formula'], t5['snf_spur3'], t5['E_conn']))
# stars
bases = t5['bases']
bidx = {b: k for k, b in enumerate(bases)}
P5 = t5['P']
stars = []
for w in range(5):
    Qs = [x for x in range(5) if x != w]
    v = np.zeros(10, dtype=np.int64)
    for pos, y in enumerate(Qs):
        comp = tuple(sorted(set(range(5)) - {w, y}))
        v[bidx[comp]] = (-1) ** pos
    stars.append(v)
stars = np.array(stars)
print('    P v_w = 0 over Z for all five stars:', bool(((P5 @ stars.T) == 0).all()))
print('    rank of the stars: over Q', int(flint.fmpz_mat(stars.tolist()).rank()), ', mod 3', rank_mod_p(stars, 3),
      ', dim ker P mod 3 =', 10 - t5['rP'])
Si = inv_mod_p(t5['S'], 3)
E5 = (t5['D1'] - t5['L12'] @ Si @ t5['D2']) % 3
img = (E5 @ stars.T) % 3
print('    rank(E on span(stars) -> coker P) mod 3 =',
      rank_mod_p(np.concatenate([P5 % 3, img], axis=1), 3) - t5['rP'], '(4 = isomorphism onto 4-dim coker)')
# the relation among the stars
print('    star relation: sum_w (-1)^w v_w =', ((np.array([(-1) ** w for w in range(5)]) @ stars)).tolist())

# ---------------- census ----------------
nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
for n in range(4, nmax + 1):
    cls = gentourng(n)
    stats = dict(classes=len(cls), left_ok=0, right_ok=0, cong=0, reduction_exact=0, rankM_eq=0, rankN_eq=0,
                 formula_eq=0, spur3_pos=0, S_inv=0, E_eq=0, snf_eq=0)
    spur_list = []
    for s in cls:
        A = adj_from_ascii(n, s)
        res = analyse(A, want_snf=(n <= 6))
        for key in ('left_ok', 'right_ok', 'cong', 'reduction_exact', 'rankM_eq', 'rankN_eq'):
            stats[key] += int(bool(res[key]))
        stats['formula_eq'] += int(res['formula'] == res['spur3'])
        stats['spur3_pos'] += int(res['spur3'] > 0)
        if 'E_conn' in res:
            stats['S_inv'] += 1
            stats['E_eq'] += int(res['E_conn'] == res['spur3'])
        if 'snf_spur3' in res:
            stats['snf_eq'] += int(res['snf_spur3'] == res['spur3'])
        spur_list.append((s, res['spur3']))
    print('\nn=%d' % n, stats)
    print('   spur3 distribution:', dict(sorted({v: sum(1 for _, x in spur_list if x == v)
                                                  for v in set(x for _, x in spur_list)}.items())))
    if n == 5:
        print('   n=5 spur3 by class:', spur_list)
