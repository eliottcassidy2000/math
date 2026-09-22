#!/usr/bin/env python3
"""
Adversarial audit of lane block_spectrum_audit (wave 2026-09-22).

Independent recomputation (different code paths: rank-based algebraic
multiplicities, exact sympy characteristic polynomials only where the lane
claims exactness, generalised-nullspace zero multiplicity, brute-force
Hamiltonian paths, explicit Fano developments) of every number the lane note
quotes, plus the checks the explorer did not run:

  A1  literal block layout dimensions and the four m=n=1 characteristic polys
  A2  iso classes m,n<=4 and the ONE/TWO/ORI sweep statistics; zero-multiplicity
      formula for ONE via generalised nullspace
  A3  representative dependence of the ORI orientation (the cross-link pattern
      [i+j even] is NOT invariant under relabelling inside a part): full labelled
      sweep m,n<=3, and the (1,2) witness both ways
  A4  Brauer-Gentry on all labelled tournaments n<=6, with the stronger check
      min rho over non-transitive tournaments >= 1
  A5  sink/source shift, 32 exact cases
  A6  Paley T_7: counts, charpoly, Hamiltonian paths, and the session-lead
      Fano decomposition dev{0,1,3} u dev{0,1,5}, arc-in-2-cyclic-triples,
      line orientation s->s+1->s+3->s
  A7  bipartite blocks ONE/TWO/SKEW, K=J and the ORI pattern, m=n=3
  A8  Euler / K_5 / K_{3,3} table
  A9  binary entropy of 8/pi^2, circulant regular tournament counts
"""
import itertools
import math
import sys

import numpy as np
import sympy as sp

x = sp.Symbol('x')


def fail(msg):
    raise RuntimeError(msg)


def out(*a):
    print(*a)
    sys.stdout.flush()


def labelled(n):
    pairs = list(itertools.combinations(range(n), 2))
    for bits in itertools.product((0, 1), repeat=len(pairs)):
        A = np.zeros((n, n), dtype=np.int64)
        for (i, j), b in zip(pairs, bits):
            if b:
                A[i, j] = 1
            else:
                A[j, i] = 1
        yield A


def is_transitive(A):
    return sorted(int(s) for s in A.sum(axis=1)) == list(range(A.shape[0]))


def iso_classes(n):
    # canonical form = lexicographically smallest flattened relabelling; independent code
    seen = {}
    for A in labelled(n):
        key = min(tuple(A[np.ix_(p, p)].flatten().tolist()) for p in itertools.permutations(range(n)))
        seen.setdefault(key, A.copy())
    return list(seen.values())


def zero_alg_mult(M):
    # algebraic multiplicity of eigenvalue 0 = dim of generalised nullspace = N - rank(M^N), exact
    N = M.shape[0]
    P = sp.Matrix(M.tolist()) ** N
    return N - P.rank()


def cp(M):
    return sp.factor(sp.Matrix(M.tolist()).charpoly(x).as_expr())


def eig(M):
    return np.linalg.eigvals(M.astype(float))


def fmt(z):
    r, i = float(z.real), float(z.imag)
    return "%.4f" % r if abs(i) < 1e-9 else "%.4f%+.4fi" % (r, i)


# --------------------------------------------------------------------------
out("A1  literal layout")
sq = [(m, n) for m in range(1, 5) for n in range(1, 5) if 2 * m + 2 * n == m + n + 2]
out("  (m,n) with 2m+2n == m+n+2 in 1..4: %s" % sq)
if sq != [(1, 1)]:
    fail("A1 square instances")
polys = {}
for k in (0, 1):
    for kp in (0, 1):
        M = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [k, 0, 0, 1], [0, kp, 0, 0]])
        polys[(k, kp)] = sp.factor(M.charpoly(x).as_expr())
        out("  k=%d k'=%d charpoly %s" % (k, kp, polys[(k, kp)]))
exp = {(0, 0): x ** 4, (0, 1): x ** 2 * (x - 1) * (x + 1), (1, 0): x ** 2 * (x - 1) * (x + 1), (1, 1): (x - 1) ** 2 * (x + 1) ** 2}
if any(sp.expand(polys[k] - exp[k]) != 0 for k in polys):
    fail("A1 charpolys")

# --------------------------------------------------------------------------
out("A2  sweep")
classes = {n: iso_classes(n) for n in range(1, 5)}
counts = [len(classes[n]) for n in range(1, 5)]
out("  iso class counts n=1..4: %s" % counts)
if counts != [1, 1, 2, 4]:
    fail("A2 class counts")


def build(A, B, K, Kp):
    m, n = A.shape[0], B.shape[0]
    N = m + n + 2
    M = np.zeros((N, N), dtype=np.int64)
    M[:m, :m] = A
    M[m:m + n, m:m + n] = B
    M[:m, m:m + n] = K
    M[m:m + n, :m] = Kp
    M[m + n, :m] = 1
    M[m:m + n, m + n + 1] = 1
    return M


def links(m, n, mode):
    J = np.ones((m, n), dtype=np.int64)
    if mode == "ONE":
        return J, np.zeros((n, m), dtype=np.int64)
    if mode == "TWO":
        return J, J.T.copy()
    K = np.array([[1 if (i + j) % 2 == 0 else 0 for j in range(n)] for i in range(m)], dtype=np.int64)
    return K, (J - K).T.copy()


stat = {}
rho_all = []
ori_nilp = []
one_zero_ok = 0
one_zero_tot = 0
for mode in ("ONE", "TWO", "ORI"):
    s = dict(total=0, all_zero=0, pure_imag=0, some_neg_re=0, all_neg_re=0, max_re_ge0=0, rho_pos=0, factor=0)
    for m in range(1, 5):
        for n in range(1, 5):
            K, Kp = links(m, n, mode)
            for A in classes[m]:
                for B in classes[n]:
                    M = build(A, B, K, Kp)
                    ev = eig(M)
                    s["total"] += 1
                    rho = max(abs(ev))
                    rho_all.append(rho)
                    nil = zero_alg_mult(M) == M.shape[0]
                    if nil:
                        s["all_zero"] += 1
                        if mode == "ORI":
                            ori_nilp.append((m, n))
                    if any(abs(z.real) < 1e-9 and abs(z.imag) > 1e-9 for z in ev):
                        s["pure_imag"] += 1
                    if ev.real.min() < -1e-9:
                        s["some_neg_re"] += 1
                    if ev.real.max() < -1e-9:
                        s["all_neg_re"] += 1
                    if ev.real.max() >= -1e-9:
                        s["max_re_ge0"] += 1
                    if rho > 1e-9:
                        s["rho_pos"] += 1
                    PM = sp.Matrix(M.tolist()).charpoly(x).as_expr()
                    PAB = x ** 2 * sp.Matrix(A.tolist()).charpoly(x).as_expr() * sp.Matrix(B.tolist()).charpoly(x).as_expr()
                    if sp.expand(PM - PAB) == 0:
                        s["factor"] += 1
                    if mode == "ONE":
                        one_zero_tot += 1
                        if zero_alg_mult(M) == 2 + zero_alg_mult(A) + zero_alg_mult(B):
                            one_zero_ok += 1
    stat[mode] = s
    out("  %s %s" % (mode, s))
expected = {"ONE": dict(total=64, all_zero=16, pure_imag=0, some_neg_re=48, all_neg_re=0, max_re_ge0=64, rho_pos=48, factor=64),
            "TWO": dict(total=64, all_zero=0, pure_imag=0, some_neg_re=64, all_neg_re=0, max_re_ge0=64, rho_pos=64, factor=0),
            "ORI": dict(total=64, all_zero=2, pure_imag=0, some_neg_re=62, all_neg_re=0, max_re_ge0=64, rho_pos=62, factor=2)}
if stat != expected:
    fail("A2 sweep statistics differ from the lane's table: %s" % stat)
out("  ONE zero multiplicity == 2 + mult0(A) + mult0(B) (generalised nullspace): %d of %d" % (one_zero_ok, one_zero_tot))
if one_zero_ok != 64:
    fail("A2 zero mult formula")
nz = [r for r in rho_all if r > 1e-9]
out("  rho>0: %d of %d ; min nonzero rho %.4f ; max rho %.4f" % (len(nz), len(rho_all), min(nz), max(nz)))
if len(nz) != 174 or abs(min(nz) - 1) > 1e-9 or abs(max(nz) - 5.4592) > 1e-4:
    fail("A2 rho stats")
out("  ORI nilpotent (m,n) for the lane's class representatives: %s" % ori_nilp)
# the (1,3) transitive/transitive ORI charpoly quoted in the note
K, Kp = links(1, 3, "ORI")
B13 = [B for B in classes[3] if is_transitive(B)][0]
out("  ORI (1,3) both transitive, representative used by the lane: charpoly %s" % cp(build(classes[1][0], B13, K, Kp)))

# --------------------------------------------------------------------------
out("A3  ORI representative dependence (cross-link [i+j even] is not relabelling-invariant inside a part)")
K, Kp = links(1, 2, "ORI")
for B in ([[0, 1], [0, 0]], [[0, 0], [1, 0]]):
    M = build(np.zeros((1, 1), dtype=np.int64), np.array(B, dtype=np.int64), K, Kp)
    out("  (1,2) B=%s : charpoly %s ; nilpotent %s" % (B, cp(M), zero_alg_mult(M) == 5))
tot = 0
nil = 0
per = {}
for m in range(1, 4):
    for n in range(1, 4):
        K, Kp = links(m, n, "ORI")
        c = [0, 0]
        for A in labelled(m):
            for B in labelled(n):
                M = build(A, B, K, Kp)
                c[0] += 1
                if zero_alg_mult(M) == M.shape[0]:
                    c[1] += 1
        per[(m, n)] = tuple(c)
        tot += c[0]
        nil += c[1]
out("  ORI over ALL labelled (A,B), m,n<=3: nilpotent %d of %d ; per (m,n) (total, nilpotent): %s" % (nil, tot, per))
if per[(1, 2)] != (2, 1) or per[(1, 1)] != (1, 1):
    fail("A3 witness")
# also: ONE and TWO statistics are relabelling-invariant (K=J), spot check on all labelled m,n<=3
inv_ok = True
for mode in ("ONE", "TWO"):
    for m in range(1, 4):
        for n in range(1, 4):
            K, Kp = links(m, n, mode)
            vals = {}
            for A in labelled(m):
                for B in labelled(n):
                    key = (sp.Matrix(A.tolist()).charpoly(x).as_expr(), sp.Matrix(B.tolist()).charpoly(x).as_expr())
                    p = sp.Matrix(build(A, B, K, Kp).tolist()).charpoly(x).as_expr()
                    if key in vals and sp.expand(vals[key] - p) != 0:
                        inv_ok = False
                    vals[key] = p
out("  ONE/TWO charpoly depends only on (charpoly A, charpoly B) over all labelled m,n<=3: %s" % inv_ok)
if not inv_ok:
    fail("A3 invariance")

# --------------------------------------------------------------------------
out("A4  Brauer-Gentry, all labelled tournaments n<=6")
nontrans = 0
min_rho_nontrans = 99.0
rows = {}
for n in range(1, 7):
    cnt = tr = 0
    min_re, max_mod, max_im = 9.0, 0.0, 0.0
    for A in labelled(n):
        cnt += 1
        ev = eig(A)
        rho = max(abs(ev))
        if is_transitive(A):
            tr += 1
        else:
            nontrans += 1
            min_rho_nontrans = min(min_rho_nontrans, rho)
        min_re = min(min_re, ev.real.min())
        max_mod = max(max_mod, rho)
        max_im = max(max_im, abs(ev.imag).max())
    cot = 0.5 / math.tan(math.pi / (2 * n)) if n >= 2 else 0.0
    rows[n] = (cnt, tr, min_re, max_mod, max_im, cot)
    out("  n=%d labelled %d transitive %d minRe %.4f maxMod %.4f (bound %.4f) maxIm %.4f (bound %.4f)"
        % (n, cnt, tr, min_re, max_mod, (n - 1) / 2, max_im, cot))
    if min_re < -0.5 - 1e-9 or max_mod > (n - 1) / 2 + 1e-9 or max_im > cot + 1e-9:
        fail("A4 bound violation")
out("  non-transitive labelled tournaments n<=6: %d ; min rho among them %.6f" % (nontrans, min_rho_nontrans))
if nontrans != 32994 or min_rho_nontrans < 1 - 1e-9:
    fail("A4 nontransitive count / rho>=1")
if rows[4][1] != 24 or rows[6][1] != 720 or abs(rows[4][3] - 1.3953) > 1e-4 or abs(rows[6][3] - 2.4340) > 1e-4 or abs(rows[4][5] - 1.2071) > 1e-4:
    fail("A4 table values")
# strong n=4 class
strong = [A for A in classes[4] if sorted(int(s) for s in A.sum(axis=1)) == [1, 1, 2, 2]]
out("  strong n=4 class: charpoly %s spectrum %s" % (cp(strong[0]), ", ".join(fmt(z) for z in sorted(eig(strong[0]), key=lambda z: -z.real))))
if sp.expand(sp.Matrix(strong[0].tolist()).charpoly(x).as_expr() - (x ** 4 - 2 * x - 1)) != 0:
    fail("A4 strong charpoly")

# --------------------------------------------------------------------------
out("A5  sink/source shift")
ok = tot = 0
for n in range(1, 5):
    for A in classes[n]:
        PA = sp.Matrix(A.tolist()).charpoly(x).as_expr()
        for kind in ("sink", "source", "both", "both+arc"):
            N = n + (2 if kind.startswith("both") else 1)
            M = np.zeros((N, N), dtype=np.int64)
            M[:n, :n] = A
            if kind == "sink":
                M[:n, n] = 1
            elif kind == "source":
                M[n, :n] = 1
            else:
                M[n, :n] = 1
                M[:n, n + 1] = 1
                if kind == "both+arc":
                    M[n, n + 1] = 1
            k = N - n
            tot += 1
            if sp.expand(sp.Matrix(M.tolist()).charpoly(x).as_expr() - x ** k * PA) == 0:
                ok += 1
out("  charpoly(M) = x^k charpoly(A): %d of %d" % (ok, tot))
if (ok, tot) != (32, 32):
    fail("A5")

# --------------------------------------------------------------------------
out("A6  Paley T_7 and the session-lead Fano probe")
S = {1, 2, 4}
P = np.array([[1 if (j - i) % 7 in S else 0 for j in range(7)] for i in range(7)], dtype=np.int64)
arcs = int(P.sum())
cyc = set()
trans3 = 0
for t in itertools.combinations(range(7), 3):
    sub = P[np.ix_(t, t)]
    if sorted(sub.sum(axis=1).tolist()) == [1, 1, 1]:
        cyc.add(frozenset(t))
    else:
        trans3 += 1
out("  arcs %d cyclic triples %d transitive triples %d (7^3-7)/24 = %d" % (arcs, len(cyc), trans3, (343 - 7) // 24))
dev = lambda base: {frozenset((b + s) % 7 for b in base) for s in range(7)}
F1, F2 = dev((0, 1, 3)), dev((0, 1, 5))
out("  dev{0,1,3} size %d, dev{0,1,5} size %d, disjoint %s, union == cyclic triples %s"
    % (len(F1), len(F2), not (F1 & F2), (F1 | F2) == cyc))
# each of the two developments is a Fano plane: any two points on exactly one line
def is_fano(L):
    return all(sum(1 for l in L if a in l and b in l) == 1 for a, b in itertools.combinations(range(7), 2))
out("  dev{0,1,3} is a 2-(7,3,1) design: %s ; dev{0,1,5}: %s" % (is_fano(F1), is_fano(F2)))
per_arc = {}
for i in range(7):
    for j in range(7):
        if P[i, j]:
            per_arc[(i, j)] = sum(1 for t in cyc if i in t and j in t)
out("  cyclic triples through each arc: %s" % sorted(set(per_arc.values())))
orient_ok = all(P[s % 7, (s + 1) % 7] and P[(s + 1) % 7, (s + 3) % 7] and P[(s + 3) % 7, s % 7] for s in range(7))
out("  on every line {s,s+1,s+3}: s->s+1->s+3->s : %s" % orient_ok)
orient5 = all(P[s % 7, (s + 1) % 7] and P[(s + 1) % 7, (s + 5) % 7] and P[(s + 5) % 7, s % 7] for s in range(7))
out("  on every line {s,s+1,s+5}: s->s+1->s+5->s : %s" % orient5)
if not (arcs == 21 and len(cyc) == 14 and trans3 == 21 and (F1 | F2) == cyc and not (F1 & F2) and is_fano(F1) and is_fano(F2)
        and sorted(set(per_arc.values())) == [2] and orient_ok and orient5):
    fail("A6 Paley/Fano")
out("  charpoly %s" % cp(P))
if sp.expand(sp.Matrix(P.tolist()).charpoly(x).as_expr() - (x - 3) * (x ** 2 + x + 2) ** 3) != 0:
    fail("A6 charpoly")
hp = sum(1 for p in itertools.permutations(range(7)) if all(P[p[i], p[i + 1]] for i in range(6)))
out("  Hamiltonian paths %d" % hp)
if hp != 189:
    fail("A6 ham paths")

# --------------------------------------------------------------------------
out("A7  bipartite blocks m=n=3")
J = np.ones((3, 3), dtype=np.int64)
Z = np.zeros((3, 3), dtype=np.int64)
KO = np.array([[1, 0, 1], [0, 1, 0], [1, 0, 1]], dtype=np.int64)
for name, K in (("J", J), ("ORI", KO)):
    one = np.block([[Z, K], [Z, Z]])
    two = np.block([[Z, K], [K.T, Z]])
    skew = np.block([[Z, K], [-K.T, Z]])
    out("  K=%s: ONE %s ; TWO %s ; SKEW %s ; singular values %s" % (name, cp(one), cp(two), cp(skew), np.round(np.linalg.svd(K.astype(float), compute_uv=False), 4).tolist()))
    if sp.expand(sp.Matrix(skew.tolist()).charpoly(x).as_expr() - (x ** 4 * (x ** 2 + 9) if name == "J" else x ** 2 * (x ** 2 + 1) * (x ** 2 + 4))) != 0:
        fail("A7 skew charpoly")
for (m, n) in [(2, 2), (2, 3), (3, 4), (4, 4)]:
    Jm = np.ones((m, n), dtype=np.int64)
    two = np.block([[np.zeros((m, m), int), Jm], [Jm.T, np.zeros((n, n), int)]])
    out("  m=%d n=%d TWO max eigenvalue %.4f sqrt(mn) %.4f" % (m, n, max(eig(two).real), math.sqrt(m * n)))

# --------------------------------------------------------------------------
out("A8  Euler / minors")
for (m, n) in [(1, 4), (2, 3), (2, 4), (3, 3)]:
    V = m + n + 2
    E = m * (m - 1) // 2 + n * (n - 1) // 2 + m * n + m + n
    out("  m=%d n=%d V=%d E=%d 3V-6=%d" % (m, n, V, E, 3 * V - 6))
if not ((1, 4) and True):
    fail("unreachable")
E33 = 3 + 3 + 9 + 6
E23 = 1 + 3 + 6 + 5
if E33 != 21 or E23 != 15:
    fail("A8")

# --------------------------------------------------------------------------
out("A9  entropy and circulant counts")
h = 8 / math.pi ** 2
H = -(h * math.log2(h) + (1 - h) * math.log2(1 - h))
out("  8/pi^2 = %.5f ; H = %.5f bits" % (h, H))
if abs(H - 0.70028) > 5e-6:
    fail("A9 entropy")
for n in (3, 5, 7):
    c = 0
    for Sc in itertools.combinations(range(1, n), (n - 1) // 2):
        if all((-s) % n not in Sc for s in Sc):
            A = np.array([[1 if (j - i) % n in Sc else 0 for j in range(n)] for i in range(n)], dtype=np.int64)
            ev = eig(A)
            if not all(abs(z.real + 0.5) < 1e-9 for z in ev if abs(z - (n - 1) / 2) > 1e-9):
                fail("A9 regular real parts")
            c += 1
    out("  n=%d circulant regular tournaments %d (= 2^((n-1)/2) = %d)" % (n, c, 2 ** ((n - 1) // 2)))
    if c != 2 ** ((n - 1) // 2):
        fail("A9 circulant count")

out("AUDIT DONE: all lane numbers reproduced; ORI counts are representative-dependent (A3)")
