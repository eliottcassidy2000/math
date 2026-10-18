#!/usr/bin/env python3
"""Independent mod-P check of the note's 4d-table row for N=4 SU(2) SYM:
simple small N=4 at c=-9: (i) dim L vs the gauge-theory Schur index (contour integral),
(ii) classical freeness through a weight bound, (iii) nilpotency index of T in R_V.
The maximal ideal is computed from scratch:  w in M  iff  X w in M for the positive generators
X in {L_1, L_2, J^{+,3,-}_1, G1/G2/Gb1/Gb2_{1/2}} (these generate the positive part).
Usage: python n4_audit.py maxw2 [P]      (maxw2 = 2*max weight, e.g. 14 for weight 7)
"""
import sys
import time
import random
from fractions import Fraction
import flint
from alin import rref_vectors, rank, nullspace, leading_forms
from n4mod import (N4Vac, bracket_raw, TYPES, ODD, DELTA2, CH2, modekey, smul, sd, wt2, ch2, li,
                   pbw_basis, creation_modes)

maxw2 = int(sys.argv[1]) if len(sys.argv) > 1 else 12
P = int(sys.argv[2]) if len(sys.argv) > 2 else 2147483647
t0 = time.time()


def log(*a):
    print(*a, flush=True)


# ---------------- (0) super-Jacobi check of the bracket table (exact, generic c) -------------
def jac_check(ntests=4000, seed=1):
    rnd = random.Random(seed)
    c = Fraction(7, 3)
    one = Fraction(1)

    def br(X, Y):
        o, ce = bracket_raw(X, Y, c, one)
        return o, ce

    def br_lin(X, vec):  # [X, sum c_Z Z] (central parts of vec ignored: central commutes)
        out = {}
        cen = Fraction(0)
        for Z, cz in vec.items():
            o, ce = br(X, Z)
            for k, v in o.items():
                out[k] = out.get(k, 0) + cz * v
            cen += cz * ce
        return out, cen

    def lin_br(vec, Z):
        out = {}
        cen = Fraction(0)
        for X, cx in vec.items():
            o, ce = br(X, Z)
            for k, v in o.items():
                out[k] = out.get(k, 0) + cx * v
            cen += cx * ce
        return out, cen
    bad = 0
    for _ in range(ntests):
        X, Y, Z = [(rnd.choice(TYPES), 0) for _ in range(3)]
        X = (X[0], rnd.randint(-4, 4) * 2 + (1 if X[0] in ODD else 0))
        Y = (Y[0], rnd.randint(-4, 4) * 2 + (1 if Y[0] in ODD else 0))
        Z = (Z[0], rnd.randint(-4, 4) * 2 + (1 if Z[0] in ODD else 0))
        px, py = (X[0] in ODD), (Y[0] in ODD)
        yz, _ = br(Y, Z)
        lhs, lc = br_lin(X, yz)
        xy, _ = br(X, Y)
        r1, c1 = lin_br(xy, Z)
        xz, _ = br(X, Z)
        r2, c2 = br_lin(Y, xz)
        s = -1 if (px and py) else 1
        rhs = dict(r1)
        for k, v in r2.items():
            rhs[k] = rhs.get(k, 0) + s * v
        rc = c1 + s * c2
        diff = {k: lhs.get(k, 0) - rhs.get(k, 0) for k in set(lhs) | set(rhs)}
        if any(v != 0 for v in diff.values()) or lc != rc:
            bad += 1
    return bad


log("super-Jacobi violations in 4000 random triples:", jac_check())

# ---------------- (i) Schur index of N=4 SU(2) SYM by contour integral ----------------------
M = maxw2 + 2


def schur_index(M):
    # series in x = q^{1/2}; Laurent in z stored as dict exponent -> list of coeffs
    # single letter f = 2x/(1+x) * (z^2 + 1 + z^-2)
    sl = [0] * (M + 1)
    for n in range(1, M + 1):
        sl[n] = 2 * (-1) ** (n - 1)
    from collections import defaultdict
    # log PE = sum_k f(x^k, z^k)/k ; use exact Fractions
    logpe = defaultdict(lambda: [Fraction(0)] * (M + 1))
    for k in range(1, M + 1):
        for n in range(1, M // k + 1):
            for e in (2, 0, -2):
                logpe[e * k][n * k] += Fraction(sl[n], k)

    def mul(A, B):
        C = defaultdict(lambda: [Fraction(0)] * (M + 1))
        for ea, sa in A.items():
            for eb, sb in B.items():
                tgt = C[ea + eb]
                for i, x in enumerate(sa):
                    if x == 0:
                        continue
                    for j in range(0, M + 1 - i):
                        if sb[j] != 0:
                            tgt[i + j] += x * sb[j]
        return C
    # exp of logpe (no constant term in x): exp(L) = sum L^m/m!
    res = {0: [Fraction(1)] + [Fraction(0)] * M}
    term = {0: [Fraction(1)] + [Fraction(0)] * M}
    for m in range(1, M + 1):
        term = mul(term, logpe)
        term = {e: [x / m for x in s] for e, s in term.items()}
        for e, s in term.items():
            r = res.setdefault(e, [Fraction(0)] * (M + 1))
            for i in range(M + 1):
                r[i] += s[i]
    meas = {0: Fraction(1), 2: Fraction(-1, 2), -2: Fraction(-1, 2)}
    out = [Fraction(0)] * (M + 1)
    for e, cm in meas.items():
        s = res.get(-e)
        if s:
            for i in range(M + 1):
                out[i] += cm * s[i]
    return out


schur = schur_index(M)
log("Schur index coefficients (x = q^1/2):", [str(x) for x in schur[:maxw2 + 1]])

# ---------------- module ----------------------------------------------------------------------
one = flint.nmod(1, P)
c = -9 * one
V = N4Vac(c, one)
basis = pbw_basis(maxw2)
sector = {}
for w, lst in basis.items():
    for m in lst:
        sector.setdefault((w, ch2(m)), []).append(m)
log("dim V_{w/2}:", [(w / 2, len(basis.get(w, []))) for w in range(maxw2 + 1)])

posgens = [('L', 2), ('L', 4), ('J+', 2), ('J3', 2), ('J-', 2), ('G1', 1), ('G2', 1), ('Gb1', 1), ('Gb2', 1)]
Kf = {(0, 0): flint.nmod_mat(1, 1, [1], P)}
dimL = {0: 1}
Mbasis = {}
for w in range(1, maxw2 + 1):
    tot = 0
    for (ww, ch), Bs in sorted(sector.items()):
        if ww != w:
            continue
        rows = []
        for X in posgens:
            key = (w - X[1], ch + CH2[X[0]])
            if key not in Kf or Kf[key].nrows() == 0:
                continue
            tgt = sector[key]
            idx = {m: i for i, m in enumerate(tgt)}
            A = flint.nmod_mat(len(tgt), len(Bs), P)
            for j, b in enumerate(Bs):
                for mono, cc in V.act(X, b).items():
                    A[idx[mono], j] = cc
            Mb = Kf[key] * A
            ent = Mb.entries()
            nc = Mb.ncols()
            for i in range(Mb.nrows()):
                r = {Bs[j]: ent[i * nc + j] for j in range(nc) if ent[i * nc + j] != 0}
                if r:
                    rows.append(r)
        rr = rref_vectors(rows, Bs, P) if rows else []
        Km = flint.nmod_mat(len(rr), len(Bs), P)
        col = {m: j for j, m in enumerate(Bs)}
        for i, (_, row) in enumerate(rr):
            for mono, cc in row.items():
                Km[i, col[mono]] = cc
        Kf[(w, ch)] = Km
        tot += len(rr)
        Mbasis[(w, ch)] = nullspace([row for _, row in rr], Bs, P) if rr else \
            [{b: one} for b in Bs]
    dimL[w] = tot
log("dim L_{w/2}:", [(w / 2, dimL[w]) for w in range(maxw2 + 1)])
signed = [(-1) ** w * dimL[w] for w in range(maxw2 + 1)]
log("matches Schur index (with (-1)^F):", [Fraction(s) == schur[w] for w, s in enumerate(signed)])

# ---------------- grI, I_0, J, ker(pi) ---------------------------------------------------------
gr = {}
for key, vecs in Mbasis.items():
    Bs = sector[key]
    gr[key] = leading_forms(vecs, Bs, degree=lambda m: -li(m), modulus=P) if vecs else []
J = {}
allfree = True
for w in range(1, maxw2 + 1):
    for (ww, ch), Bs in sorted(sector.items()):
        if ww != w:
            continue
        I0 = [g for g in gr[(w, ch)] if li(next(iter(g))) == 0]
        jg = list(I0)
        for g in J.get((w - 2, ch), []):
            jg.append(sd(g))
        for X in creation_modes(maxw2):
            xw = -X[1]
            key = (w - xw, ch - CH2[X[0]])
            for g in J.get(key, []):
                jg.append(smul({(X,): one}, g))
        J[(w, ch)] = [row for _, row in rref_vectors(jg, Bs, P)] if jg else []
        rg, rj = len(gr[(w, ch)]), len(J[(w, ch)])
        if rj:
            assert rank(gr[(w, ch)] + J[(w, ch)], Bs, P) == rg, ("J not in grI", w, ch)
        if rg != rj:
            allfree = False
            log(f"   ker(pi) at weight {w / 2}, charge {ch / 2}: {rg - rj}")
    log(f"weight {w / 2}: dim I = {sum(len(gr[k]) for k in gr if k[0] == w)}, "
        f"dim J = {sum(len(J[k]) for k in J if k[0] == w)}  [{time.time() - t0:.0f}s]")
log("classically free through weight", maxw2 / 2, ":", allfree)

# ---------------- nilpotency index of T -------------------------------------------------------
T = {(('L', -4),): one}
T2 = smul(T, T)
T3 = smul(T2, T)
for name, poly, w in (('T', T, 4), ('T^2', T2, 8), ('T^3', T3, 12)):
    if w <= maxw2:
        I0w = [g for g in gr.get((w, 0), []) if li(next(iter(g))) == 0]
        r0 = rank(I0w, sector[(w, 0)], P) if I0w else 0
        r1 = rank(I0w + [poly], sector[(w, 0)], P)
        log(f"{name} in I_0 (i.e. zero in R_V): {r1 == r0}")
log(f"time {time.time() - t0:.0f}s")
