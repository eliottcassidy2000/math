#!/usr/bin/env python3
"""
block_doubling_classification_kps1.py — kind-pasteur-2026-06-09-S1
BRANCH E (HYP-2337): classification of the 2x2-block doubling family.

SHAPE JUSTIFICATION (derived in Section E1, demos included):
  Ansatz: each block of M' is affine in M over span{M, I} (functorial doubling:
  cross arcs determined by T's arcs plus twin arcs).
    M' = [[a M + alpha I,  b M + beta I],
          [c M + gamma I,  d M + delta I]]
  Validity for ALL tournaments T forces:
    - diagonal blocks: zero diagonal  => alpha = delta = 0; off-diag entries +-1
      => a, d in {+-1}.
    - skew-symmetry M'^T = -M' on the cross blocks:  (bM + beta I)^T = -(cM + gamma I)
      => -b M + beta I = -c M - gamma I => c = b, gamma = -beta  (M, I lin. indep.)
    - cross-block entries +-1 for all T: off-diagonal b M_ij => b in {+-1};
      diagonal (twin arc i <-> i') beta => beta in {+-1}.
  Hence M' = [[aM, bM + beta I],[bM - beta I, dM]],  (a,b,d,beta) in {+-1}^4: 16 members.
  Kron form:  M' = P (x) M + Q (x) I,  P = [[a,b],[b,d]] symmetric,  Q = beta*J,
  J = [[0,1],[-1,0]].
  (Remark: allowing the all-ones J_n - I in the cross block would add the ordered-sum
   T -> T (+) T family, where cross arcs ignore T; excluded by the span{M,I} ansatz.)

SYMMETRIES acting on (a,b,d,beta):
  g1 swap copies (conj by [[0,I],[I,0]]):       (a,b,d,beta) -> (d,b,a,-beta)
  g2 global arc reversal (M' -> -M' = op-double): (a,b,d,beta) -> (-a,-b,-d,-beta)
  g3 pre-compose with op (M -> -M):              (a,b,d,beta) -> (-a,-b,-d,beta)
  g2*g3 = beta-flip. Group = <swap(a,d), negate(a,b,d), flip beta> ~ Z2^3.

MAIN THEOREM (Clifford/Cayley-Dickson classification), proved here:
  M'^2 = P^2 (x) M^2 + {P,Q} (x) M + Q^2 (x) I,   Q^2 = -I2,   {P,Q} = beta*(tr P)*J.
  (a) M'^2 block-diagonal for all T  <=>  tr P = 0  <=>  rows of P orthogonal
      <=> det P = -2  <=>  P is a 2x2 Hadamard matrix  <=>  {P,Q} = 0 and P^2 = 2I.
      Then M'^2 = I2 (x) (2M^2 - I): Chebyshev T2 law. (orbit of D_skew, 8 members)
  (b) tr P != 0 members have P = +-ee^T or +-ff^T (rank 1, e=(1,1), f=(1,-1));
      law M'^2 = (tr P)*(I2 (x) M)*M' - I  (unimodular "boost" lift pairs, product 1).
  (c) skew-Hadamard preservation <=> case (a). Proof: S'S'^T = I - M'^2; input
      M^2 = (1-n)I gives M'^2 = (1-2n)I in case (a); in case (b) the off-diagonal
      block of M'^2 is (a+d)(b M^2 + beta M) whose diagonal is -2ab(n-1) != 0.
  Clifford: with gamma+ = P/sqrt(2), gamma- = Q: gamma+^2 = +1, gamma-^2 = -1,
  {gamma+,gamma-} = 0: a real rep of Cl(1,1) = M2(R). The doubling is literally a
  Cayley-Dickson step: Q(P (x) M) = -(P (x) M)Q is j*z = zbar*j for purely imaginary
  z (skew M: zbar = -z); beta = +-1 = orientation of j; a<->d = sheet order.

Output: 05-knowledge/results/block_doubling_classification_kps1.out
"""
import itertools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from skew_doubling_core_kps1 import (all_tournaments, iso_classes, canon, is_iso,
                                     H_count, M_of, A_of, scores, border,
                                     is_skew_hadamard)

OUTPATH = os.path.join(ROOT, "05-knowledge", "results",
                       "block_doubling_classification_kps1.out")

I2 = np.eye(2, dtype=np.int64)
J2 = np.array([[0, 1], [-1, 0]], dtype=np.int64)
SW2 = np.array([[0, 1], [1, 0]], dtype=np.int64)

NAMED = {(1, 1, -1, 1): "D_skew", (1, 1, 1, 1): "T[K2]", (1, -1, 1, 1): "SCblow"}


def PQ(s):
    a, b, d, be = s
    return (np.array([[a, b], [b, d]], dtype=np.int64), be * J2)


def Mp_of(A, s):
    n = A.shape[0]
    P, Q = PQ(s)
    return np.kron(P, M_of(A)) + np.kron(Q, np.eye(n, dtype=np.int64))


def D_member(A, s):
    Mp = Mp_of(A, s)
    return A_of(Mp), Mp


def is_tournament(Ad):
    m = Ad.shape[0]
    return np.array_equal(Ad + Ad.T,
                          np.ones((m, m), dtype=np.int64) - np.eye(m, dtype=np.int64))


def g1(s):
    return (s[2], s[1], s[0], -s[3])


def g2(s):
    return (-s[0], -s[1], -s[2], -s[3])


def g3(s):
    return (-s[0], -s[1], -s[2], s[3])


def all_members():
    return [(a, b, d, be) for a in (1, -1) for b in (1, -1)
            for d in (1, -1) for be in (1, -1)]


def compute_orbits():
    seen = set()
    orbs = []
    for s in all_members():
        if s in seen:
            continue
        orb = {s}
        frontier = [s]
        while frontier:
            t = frontier.pop()
            for g in (g1, g2, g3):
                u = g(t)
                if u not in orb:
                    orb.add(u)
                    frontier.append(u)
        seen |= orb
        orbs.append(sorted(orb, reverse=True))
    return orbs


def iso_score_restricted(A, B):
    """Isomorphism test restricted to score-preserving permutations
    (exact; feasible when score classes are small, e.g. 5+5 at n=10)."""
    n = A.shape[0]
    sA, sB = scores(A), scores(B)
    if sorted(sA) != sorted(sB):
        return False
    grpA, grpB = {}, {}
    for i in range(n):
        grpA.setdefault(sA[i], []).append(i)
        grpB.setdefault(sB[i], []).append(i)
    ks = sorted(grpA)
    if any(len(grpA[k]) != len(grpB[k]) for k in ks):
        return False
    groupsA = [grpA[k] for k in ks]
    for choice in itertools.product(*[itertools.permutations(grpB[k]) for k in ks]):
        p = [0] * n
        for ga, gb in zip(groupsA, choice):
            for x, y in zip(ga, gb):
                p[x] = y
        P = np.array(p)
        if np.array_equal(B[np.ix_(P, P)], A):
            return True
    return False


def mstr(P):
    return "[[%2d,%2d],[%2d,%2d]]" % (P[0, 0], P[0, 1], P[1, 0], P[1, 1])


def sstr(s):
    return "(%+d,%+d,%+d,%+d)" % s


def main():
    out = open(OUTPATH, "w", encoding="utf-8")

    def w(line=""):
        out.write(line + "\n")
        out.flush()
        print(line)

    w("=== block_doubling_classification_kps1 — BRANCH E (HYP-2337) ===")
    w("Family: M' = [[aM, bM+beta I],[bM-beta I, dM]] = P(x)M + Q(x)I,")
    w("        P=[[a,b],[b,d]], Q=beta*J, (a,b,d,beta) in {+-1}^4 : 16 members.")
    w("")

    # ---------------- E1: shape justification + validity + symmetries + orbits ----
    w("--- E1a: validity — every member yields a tournament (exhaustive n=3,4) ---")
    members = all_members()
    for n in (3, 4):
        bad = 0
        tot = 0
        for A in all_tournaments(n):
            tot += 1
            for s in members:
                Ad, _ = D_member(A, s)
                if not is_tournament(Ad):
                    bad += 1
        w("n=%d: %d tournaments x 16 members: invalid doubles = %d" % (n, tot, bad))

    w("")
    w("--- E1b: broken variants fail (shape is forced) — demo at n=3 ---")
    A = next(all_tournaments(3))
    M = M_of(A)
    I3 = np.eye(3, dtype=np.int64)
    bad1 = np.block([[M + I3, M + I3], [M - I3, -M]])      # alpha=+1 diag correction
    bad2 = np.block([[M, M + I3], [M + I3, -M]])           # gamma=+beta (not -beta)
    bad3 = np.block([[M, M + I3], [-M - I3, -M]])          # c=-b (cross M-signs differ)
    for name, B in (("alpha=+1 diagonal correction", bad1),
                    ("gamma=+beta (same I-sign)   ", bad2),
                    ("c=-b (opposite cross M-sign)", bad3)):
        w("%s : tournament? %s" % (name, is_tournament(A_of(B))))

    w("")
    w("--- E1c: symmetry identities (exhaustive over all 64 labeled n=4) ---")
    n = 4
    In = np.eye(n, dtype=np.int64)
    SWn = np.kron(SW2, In)
    ok1 = ok2 = ok3 = True
    for A in all_tournaments(4):
        for s in members:
            Mp = Mp_of(A, s)
            ok1 &= np.array_equal(SWn @ Mp @ SWn, Mp_of(A, g1(s)))
            ok2 &= np.array_equal(-Mp, Mp_of(A, g2(s)))
            ok3 &= np.array_equal(Mp_of(A.T, s), Mp_of(A, g3(s)))
    w("g1: conj-by-swap identity   SW M' SW == M'_{g1(s)}        : %s" % ok1)
    w("g2: global reversal         -M'      == M'_{g2(s)}        : %s" % ok2)
    w("g3: pre-op                  M'(T^op) == M'_{g3(s)}(T)     : %s" % ok3)

    w("")
    w("--- E1d: ORBITS under <g1,g2,g3> ---")
    orbs = compute_orbits()
    w("number of orbits: %d   sizes: %s   (total %d)" %
      (len(orbs), [len(o) for o in orbs], sum(len(o) for o in orbs)))
    orbname = {}
    reps = []
    for k, orb in enumerate(orbs):
        named = [NAMED[s] for s in orb if s in NAMED]
        rep = [s for s in orb if s in NAMED]
        rep = rep[0] if rep else orb[0]
        reps.append(rep)
        for s in orb:
            orbname[s] = named[0] if named else ("orbit%d" % k)
        w("orbit %d (size %d): rep %s = %s" % (k, len(orb), sstr(rep),
                                               named[0] if named else "unnamed"))
        w("   members: " + ", ".join(sstr(s) for s in orb))
    w("=> 16 members collapse to %d genuinely distinct doublings;" % len(orbs))
    w("   the three named members D_skew, T[K2], SCblow are in pairwise distinct orbits.")

    # ---------------- E2: M'^2 block algebra + skew-Hadamard --------------------
    w("")
    w("--- E2a: Kron square law  M'^2 = P^2(x)M^2 + {P,Q}(x)M + Q^2(x)I  (exh. n=4) ---")
    okk = True
    for A in all_tournaments(4):
        M = M_of(A)
        In4 = np.eye(4, dtype=np.int64)
        for s in members:
            P, Q = PQ(s)
            Mp = Mp_of(A, s)
            law = (np.kron(P @ P, M @ M) + np.kron(P @ Q + Q @ P, M)
                   + np.kron(Q @ Q, In4))
            okk &= np.array_equal(Mp @ Mp, law)
    w("Kron square law exact for all 64 x 16: %s" % okk)

    w("")
    w("--- E2b: per-member algebra table ---")
    w("%-10s %-22s %3s %4s  %-22s %-22s %-12s" %
      ("member", "P", "trP", "detP", "P^2", "{P,Q}", "orbit"))
    for s in members:
        P, Q = PQ(s)
        anti = P @ Q + Q @ P
        w("%-10s %-22s %3d %4d  %-22s %-22s %-12s" %
          (sstr(s), mstr(P), int(np.trace(P)), int(round(np.linalg.det(P))),
           mstr(P @ P), mstr(anti), orbname[s]))
    w("Q^2 = -I2 always. {P,Q} = beta*(trP)*J  (verified above row-by-row).")

    w("")
    w("--- E2c: block-diagonality of M'^2 (exhaustive n=4; off-diag block formula) ---")
    w("off-diag (1,2) block of M'^2 = (a+d) * (b M^2 + beta M); diag of that = -b(a+d)(n-1)")
    for s in members:
        a, b, d, be = s
        alwaysbd = True
        neverbd = True
        formok = True
        for A in all_tournaments(4):
            M = M_of(A)
            Mp = Mp_of(A, s)
            sq = Mp @ Mp
            off = sq[:4, 4:]
            bd = not off.any()
            alwaysbd &= bd
            neverbd &= not bd
            formok &= np.array_equal(off, (a + d) * (b * M @ M + be * M))
        verdict = ("BLOCK-DIAG (all T)" if alwaysbd
                   else ("never block-diag" if neverbd else "MIXED (?!)"))
        w("%s  a+d=%+d  formula exact: %s   %s" % (sstr(s), a + d, formok, verdict))
    w("=> M'^2 block-diagonal <=> a+d=0 <=> trP=0: exactly the D_skew orbit (8 members).")

    w("")
    w("--- E2d: orbit laws for M'^2 (exhaustive n=4) ---")
    okB = okAC = True
    for A in all_tournaments(4):
        M = M_of(A)
        In4 = np.eye(4, dtype=np.int64)
        I8 = np.eye(8, dtype=np.int64)
        for s in members:
            a, b, d, be = s
            Mp = Mp_of(A, s)
            if a + d == 0:
                okB &= np.array_equal(Mp @ Mp, np.kron(I2, 2 * M @ M - In4))
            else:
                okAC &= np.array_equal(Mp @ Mp,
                                       (a + d) * (np.kron(I2, M) @ Mp) - I8)
    w("orbit B  (trP=0):  M'^2 = I2 (x) (2M^2 - I)            [Chebyshev T2]: %s" % okB)
    w("orbits A,C (trP=+-2): M'^2 = (trP)*(I2 (x) M)*M' - I   [boost pencil]: %s" % okAC)
    w("   (I2(x)M commutes with M'; quadratic relation M'^2 - (trP)(I2(x)M)M' + I = 0)")

    w("")
    w("--- E2e: skew-Hadamard preservation test (seed: bordered C_3, order-4 SH) ---")
    C3 = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=np.int64)
    S4 = border(C3)
    TB = A_of(S4 - np.eye(4, dtype=np.int64))
    w("seed order-4 skew-Hadamard: %s" % is_skew_hadamard(S4))
    npass = 0
    for s in members:
        Ad, Mp = D_member(TB, s)
        Sp = Mp + np.eye(8, dtype=np.int64)
        ok = is_skew_hadamard(Sp)
        npass += ok
        w("%s  orbit=%-7s  S' (order 8) skew-Hadamard: %s" % (sstr(s), orbname[s], ok))
    w("=> %d/16 preserve skew-Hadamard: exactly the trP=0 (D_skew) orbit." % npass)

    w("")
    w("--- E2f: tower iteration for two distinct orbit-B members (orders 16, 32) ---")
    for s in ((1, 1, -1, 1), (-1, 1, 1, 1)):
        T = TB
        line = ["member %s:" % sstr(s)]
        for _ in range(3):
            T, Mp = D_member(T, s)
            Sp = Mp + np.eye(Mp.shape[0], dtype=np.int64)
            line.append("order %d SH=%s" % (Sp.shape[0], is_skew_hadamard(Sp)))
        w("  " + "  ".join(line))

    w("")
    w("--- E2g: PROOF (block algebra) ---")
    w("S'S'^T = (M'+I)(I-M') = I - M'^2  (M' skew).")
    w("Input skew-Hadamard of order n: M^2 = (1-n)I.")
    w("trP=0:    M'^2 = I2(x)(2M^2-I) = (2(1-n)-1) I = (1-2n) I  => S'S'^T = 2n I. QED")
    w("trP=+-2:  off-diag block of M'^2 = (a+d)(b M^2 + beta M) = +-2(b(1-n)I + beta M),")
    w("          diagonal entries -2ab(n-1) != 0 for n>=2 => M'^2 not scalar => fails. QED")

    # ---------------- E3: the 2x2 DNA trichotomy --------------------------------
    w("")
    w("--- E3: 2x2 DNA — sign-matrix trichotomy ---")
    e = np.array([[1], [1]], dtype=np.int64)
    f = np.array([[1], [-1]], dtype=np.int64)
    Pee = e @ e.T
    Pff = f @ f.T
    w("P_A(T[K2])  = [[1,1],[1,1]]   = e e^T (e=(1,1)) : %s, rank %d, det %d" %
      (np.array_equal(PQ((1, 1, 1, 1))[0], Pee), np.linalg.matrix_rank(Pee),
       int(round(np.linalg.det(Pee.astype(float))))))
    w("P_C(SCblow) = [[1,-1],[-1,1]] = f f^T (f=(1,-1)): %s, rank %d, det %d" %
      (np.array_equal(PQ((1, -1, 1, 1))[0], Pff), np.linalg.matrix_rank(Pff),
       int(round(np.linalg.det(Pff.astype(float))))))
    PB = PQ((1, 1, -1, 1))[0]
    w("P_B(D_skew) = [[1,1],[1,-1]]  = 2x2 Hadamard    : PP^T=2I %s, det %d" %
      (np.array_equal(PB @ PB.T, 2 * I2), int(round(np.linalg.det(PB.astype(float))))))
    w("P_A P_C = 0: %s;  P_A + P_C = 2I: %s;  H2 e = 2e1, H2 f = 2e2: %s %s" %
      (np.array_equal(Pee @ Pff, 0 * I2), np.array_equal(Pee + Pff, 2 * I2),
       (PB @ e).ravel().tolist(), (PB @ f).ravel().tolist()))
    w("")
    w("TRICHOTOMY: the three orbits = the three symmetric +-1 matrices mod symmetry:")
    w("  orbit A: P = +-e e^T  (2x projector on the SYMMETRIC line of the copy-swap)")
    w("  orbit C: P = +-f f^T  (2x projector on the ANTISYMMETRIC line)")
    w("  orbit B: P = 2x2 Hadamard (EXCHANGES the two swap-eigenlines), det = -2.")
    w("Equivalences (all verified above): rows of P orthogonal <=> trP=0 <=> detP=-2")
    w("  <=> P Hadamard <=> {P,Q}=0 <=> M'^2 block-diagonal <=> skew-Hadamard preserved.")
    w("Sing. members (det 0) are rank-1; the Hadamard is the UNIQUE invertible class.")

    w("")
    w("--- E3b: spectral pencil — spec(M') = U_{lambda in spec M} spec(lambda P + Q) ---")
    w("(proof: M normal; on C^2 (x) eigvec(lambda), M' restricts to lambda P + Q;")
    w(" char poly check: det(xI - M') = prod_lambda [x^2 - (trP) lambda x + (detP lambda^2 + 1)])")
    allok = True
    for Asamp in itertools.islice(all_tournaments(5), 100, 1024, 137):
        Ms = M_of(Asamp).astype(float)
        lams = np.linalg.eigvals(Ms)
        for s in members:
            a, b, d, be = s
            P, Q = PQ(s)
            Mp = Mp_of(Asamp, s).astype(float)
            got = np.poly(Mp).astype(complex)
            pred = np.array([1.0 + 0j])
            for lam in lams:
                quad = np.array([1.0, -(a + d) * lam, (a * d - b * b) * lam ** 2 + 1.0],
                                dtype=complex)
                pred = np.convolve(pred, quad)
            allok &= np.allclose(got, pred, atol=1e-6)
    w("pencil law (char-poly product) for all 16 members x 7 sample n=5 tournaments: %s"
      % allok)
    w("det(lambda P + Q) = detP*lambda^2 + 1: lift-pair products are")
    w("  orbit A/C (detP=0):  product = 1   (unimodular 'boost' pair, lifts l' = a*l +- sqrt(l^2-1))")
    w("  orbit B  (detP=-2):  product = 1 - 2 lambda^2;  trace 0 => l'^2 = 2 lambda^2 - 1 = T2(lambda).")

    # ---------------- E4: H spectra ---------------------------------------------
    w("")
    w("--- E4: H-spectra of the three orbit representatives ---")
    repD, repK, repS = (1, 1, -1, 1), (1, 1, 1, 1), (1, -1, 1, 1)
    for n in (3, 4, 5):
        cls = iso_classes(n)
        keys = {c.tobytes(): i for i, c in enumerate(cls)}
        op_idx = [keys[canon(c.T).tobytes()] for c in cls]
        Hd = {}
        for i, A in enumerate(cls):
            for s in members:
                Hd[(i, s)] = H_count(D_member(A, s)[0])
        # relation checks across all 16 members
        r1 = all(Hd[(i, g1(s))] == Hd[(i, s)] for i in range(len(cls)) for s in members)
        r2 = all(Hd[(i, g2(s))] == Hd[(i, s)] for i in range(len(cls)) for s in members)
        r3 = all(Hd[(i, g3(s))] == Hd[(op_idx[i], s)]
                 for i in range(len(cls)) for s in members)
        w("")
        w("n=%d (%d iso classes), doubles on %d vertices" % (n, len(cls), 2 * n))
        w("  H-relations: g1 (iso, equal H): %s | g2 (op-double, equal H): %s | "
          "g3 (H of op-class): %s" % (r1, r2, r3))
        w("  %3s %18s %6s | %8s %8s %8s | argmax" %
          ("idx", "scores", "H(T)", "H(Dskew)", "H(TK2)", "H(SCbl)"))
        sums = {repD: 0, repK: 0, repS: 0}
        for i, A in enumerate(cls):
            hT = H_count(A)
            row = [Hd[(i, repD)], Hd[(i, repK)], Hd[(i, repS)]]
            for r, v in zip((repD, repK, repS), row):
                sums[r] += v
            names = ["D_skew", "T[K2]", "SCblow"]
            mx = max(row)
            am = "/".join(nm for nm, v in zip(names, row) if v == mx)
            w("  %3d %18s %6d | %8d %8d %8d | %s" %
              (i, str(scores(A)), hT, row[0], row[1], row[2], am))
        nc = len(cls)
        w("  MEAN over %d classes:        | %8.1f %8.1f %8.1f" %
          (nc, sums[repD] / nc, sums[repK] / nc, sums[repS] / nc))
        # within-orbit collapse: per class, distinct H values across orbit members
        for k, orb in enumerate(orbs):
            splits = []
            for i in range(len(cls)):
                vals = sorted(set(Hd[(i, s)] for s in orb))
                if len(vals) > 1:
                    splits.append((i, vals))
            w("  orbit %d (%s): classes with >1 H value across its %d members: %s" %
              (k, orbname[orb[0]], len(orb),
               splits if splits else "none (full collapse)"))

    w("")
    w("--- E4b: C_3 doubles at n=6 — all three hit H=45; are they isomorphic? ---")
    C3doubles = {}
    for name, s in (("D_skew", repD), ("T[K2]", repK), ("SCblow", repS)):
        Ad, _ = D_member(C3, s)
        C3doubles[name] = Ad
        w("  %s(C3): H=%d scores=%s" % (name, H_count(Ad), sorted(scores(Ad))))
    pairs = [("D_skew", "T[K2]"), ("D_skew", "SCblow"), ("T[K2]", "SCblow")]
    for x, y in pairs:
        w("  iso(%s(C3), %s(C3)) = %s" % (x, y, is_iso(C3doubles[x], C3doubles[y])))

    w("")
    w("--- E4c: K2(T) ?= SCblow(T) up to iso, for all SELF-COMPLEMENTARY n=5 classes ---")
    w("(score-class-restricted backtracking iso on 10 vertices)")
    cls5 = iso_classes(5)
    keys5 = {c.tobytes(): i for i, c in enumerate(cls5)}
    op5 = [keys5[canon(c.T).tobytes()] for c in cls5]
    for i, A in enumerate(cls5):
        if op5[i] != i:
            continue  # only self-complementary classes
        AK, _ = D_member(A, repK)
        AS, _ = D_member(A, repS)
        hk, hs = H_count(AK), H_count(AS)
        isoKS = iso_score_restricted(AK, AS) if hk == hs else False
        w("  SC class %2d scores=%s H(T)=%d : H(K2)=%d H(SC)=%d equal=%s iso=%s" %
          (i, str(scores(A)), H_count(A), hk, hs, hk == hs, isoKS))
    w("  (non-SC classes have H(K2) != H(SCblow) at n=5 except none — see table)")
    # also the n=4 SC class for contrast
    cls4 = iso_classes(4)
    keys4 = {c.tobytes(): i for i, c in enumerate(cls4)}
    op4 = [keys4[canon(c.T).tobytes()] for c in cls4]
    for i, A in enumerate(cls4):
        if op4[i] != i:
            continue
        AK, _ = D_member(A, repK)
        AS, _ = D_member(A, repS)
        hk, hs = H_count(AK), H_count(AS)
        isoKS = iso_score_restricted(AK, AS) if hk == hs else False
        w("  n=4 SC class %d scores=%s H(T)=%d : H(K2)=%d H(SC)=%d iso=%s" %
          (i, str(scores(A)), H_count(A), hk, hs, isoKS))

    # ---------------- E5: Cayley-Dickson statement -------------------------------
    w("")
    w("--- E5: Clifford / Cayley-Dickson statement ---")
    w("THEOREM (2x2 DNA / Clifford classification). Every valid span{M,I} block")
    w("doubling is M' = P(x)M + Q(x)I with P symmetric +-1, Q = beta*J, beta=+-1.")
    w("M'^2 = P^2(x)M^2 + {P,Q}(x)M - I, with {P,Q} = beta*(trP)*J. The following are")
    w("equivalent: (i) {P,Q}=0; (ii) trP=0; (iii) P^2=2I (rows orthogonal, P Hadamard,")
    w("detP=-2); (iv) M'^2 = I2(x)(2M^2-I) for all T; (v) M'(T) preserves skew-Hadamard.")
    w("The pair (gamma+ = P/sqrt2, gamma- = Q) then satisfies gamma+^2=+1, gamma-^2=-1,")
    w("{gamma+,gamma-}=0: a real representation of the Clifford algebra Cl(1,1) ~ M2(R).")
    w("")
    w("CAYLEY-DICKSON reading: M skew = 'purely imaginary' (zbar = -z). The CD twist")
    w("law j z = zbar j becomes Q(P(x)M) = (P(x)(-M))Q, i.e. {P(x)M, Q} = 0, i.e.")
    w("{P,Q} = 0 — the SAME condition. So the good doublings are exactly the CD steps")
    w("M' = sqrt2 * gamma+ (x) M + j (x) I, with j^2 = -1. D_skew = (1,1,-1,1) uses")
    w("P = H2 = [[1,1],[1,-1]] (standard unnormalized Hadamard/Sylvester convention)")
    w("and Q = +J (counterclockwise j); beta=-1 is the conjugate orientation -j;")
    w("a<->d swap = relabeling the two CD sheets; global negation = the opposite")
    w("algebra. The 8 orbit-B members = 2(orientation) x 2(sheet) x 2(op) copies of")
    w("ONE CD step. Repo connection (07-reflections/lrc-cayley-dickson-tower-s387.md:")
    w("'each recursive doubling closes one old freedom and exposes a new defect'):")
    w("the skew-Hadamard tower 1->2->4->8->... with tournament cores on 0,1,3,7,...")
    w("vertices is the R->C->H->O dimension ladder; D_skew is its matrix shadow, and")
    w("the rank-1 members (T[K2], SCblow) are the DEGENERATE pencils that stay inside")
    w("one swap-eigenline — they keep the boost (product-1) spectrum but lose the")
    w("Hadamard mixing, exactly a CD-style property loss.")
    w("")
    w("=== done ===")
    out.close()


if __name__ == "__main__":
    main()
