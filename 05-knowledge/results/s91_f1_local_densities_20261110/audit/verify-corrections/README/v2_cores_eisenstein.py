#!/usr/bin/env python3
"""Verifier of the corrections to the S91 README, section 3 ("Quadratic forms in the join calculus").

Question: the corrected sentence says that the lattice Z[omega] enters "only with a cyclic triangle in the
interface".  The session's output (bridges/b7_selfjoin_eisenstein.out, B10) shows this for |Z| = 3 only.  Here the
cores with |Z| <= 4 whose diagonal form is equivalent to a multiple of x^2 + xy + y^2 are listed with the
structure of Z (cyclic triangle or not).

The core polynomial is recomputed from its definition in bridges/b7_selfjoin_eisenstein.py (docstring B): the
left piece is a chain C1 with apex a, the right piece a chain C2 with apex b, C1 -> everything outside the left
piece, everything outside the right piece -> C2, an arbitrary tournament K on the core {a, b} u Z.  A Hamiltonian
path visits the core in some order; consecutive core vertices u, v need the arc u -> v unless the chain C1 can be
inserted after u = a or the chain C2 before v = b.  Only the quadratic part (c0, c1 + c2, c3) is needed: it counts
the core orders by (a first or not) x (b last or not), with weight x1 or y1 and x2 or y2.  The quadratic part is
checked here against direct counts of Hamiltonian paths of explicit gluings (second differences in the states).
"""
import itertools
import math

import numpy as np
from numba import njit


def pair_index(r):
    idx = np.zeros((r, r), np.int64)
    k = 0
    for i in range(r):
        for j in range(i + 1, r):
            idx[i, j] = idx[j, i] = k
            k += 1
    return idx, k


@njit(cache=False)
def quad_part(code, r, perms, pidx):
    """(c0, c1, c2, c3): numbers of admissible core orders with (a first, b last), (a first, b not last),
    (a not first, b last), (neither).  Admissible: every consecutive pair u, v has u -> v in K, except that the
    step out of a (C1 may be inserted) and the step into b (C2 may be inserted) are free."""
    K = np.zeros((r, r), np.int64)
    for i in range(r):
        for j in range(i + 1, r):
            if (code >> pidx[i, j]) & 1:
                K[i, j] = 1
            else:
                K[j, i] = 1
    c = np.zeros(4, np.int64)
    for t in range(perms.shape[0]):
        ok = True
        for i in range(r - 1):
            u = perms[t, i]
            v = perms[t, i + 1]
            if u != 0 and v != 1 and K[u, v] == 0:
                ok = False
                break
        if not ok:
            continue
        afirst = perms[t, 0] == 0
        blast = perms[t, r - 1] == 1
        if afirst and blast:
            c[0] += 1
        elif afirst:
            c[1] += 1
        elif blast:
            c[2] += 1
        else:
            c[3] += 1
    return c


@njit(cache=False)
def cyclic_triangle_in_Z(code, r, pidx):
    K = np.zeros((r, r), np.int64)
    for i in range(r):
        for j in range(i + 1, r):
            if (code >> pidx[i, j]) & 1:
                K[i, j] = 1
            else:
                K[j, i] = 1
    for i in range(2, r):
        for j in range(i + 1, r):
            for k in range(j + 1, r):
                if (K[i, j] and K[j, k] and K[k, i]) or (K[j, i] and K[k, j] and K[i, k]):
                    return True
    return False


def reduce_form(a, b, c):
    """reduced positive definite binary quadratic form equivalent to (a, b, c)."""
    while True:
        if abs(b) > a:
            k = (b + a) // (2 * a) if b > 0 else -((-b + a) // (2 * a))
            b, c = b - 2 * k * a, a * k * k - b * k + c
        if a > c:
            a, b, c = c, -b, a
            continue
        if abs(b) <= a <= c:
            if b < 0 and (abs(b) == a or a == c):
                b = -b
            return a, b, c


RES = []


def check(name, cond, detail=""):
    RES.append(bool(cond))
    print(("PASS  " if cond else "FAIL  ") + name + (("   [" + str(detail) + "]") if detail != "" else ""), flush=True)


out = {}
for r in (3, 4, 5, 6):
    pidx, npairs = pair_index(r)
    perms = np.array(list(itertools.permutations(range(r))), dtype=np.int64)
    forms = {}
    for code in range(1 << npairs):
        c = quad_part(code, r, perms, pidx)
        key = (int(c[0]), int(c[1] + c[2]), int(c[3]))
        forms.setdefault(key, []).append(code)
    eis = {}
    for (al, be, ga), codes in forms.items():
        D = be * be - 4 * al * ga
        if D >= 0 or al <= 0:
            continue
        g = math.gcd(math.gcd(al, be), ga)
        red = reduce_form(al // g, be // g, ga // g)
        if red == (1, 1, 1):
            cyc = sum(1 for code in codes if cyclic_triangle_in_Z(code, r, pidx))
            eis[(al, be, ga)] = (len(codes), cyc, g)
    out[r] = (len(forms), eis)
    print("|Z| = %d: %d cores, %d distinct quadratic parts; forms equivalent to a multiple of x^2 + xy + y^2: %s"
          % (r - 2, 1 << npairs, len(forms),
             {k: "%d cores, %d with a cyclic triangle in Z, content %d" % v for k, v in eis.items()}), flush=True)

check("A  numbers of distinct diagonal forms at |Z| = 1, 2, 3, 4 are 6, 20, 93, 628 (as in b7_selfjoin_eisenstein.out)",
      [out[r][0] for r in (3, 4, 5, 6)] == [6, 20, 93, 628], [out[r][0] for r in (3, 4, 5, 6)])
check("B  no form equivalent to a multiple of the Eisenstein form at |Z| <= 2", not out[3][1] and not out[4][1])
check("C  |Z| = 3: only (3, 9, 7), 24 cores, all with Z a cyclic triangle", out[5][1] == {(3, 9, 7): (24, 24, 1)}, out[5][1])
e4 = out[6][1]
check("D  |Z| = 4: only the form (5, 25, 35) = 5 (x^2 + 5xy + 7y^2)", list(e4) == [(5, 25, 35)], e4)
if e4:
    (n4, c4, g4), = e4.values()
    print("     |Z| = 4, form (5, 25, 35): %d cores, of which %d have a cyclic triangle inside Z" % (n4, c4))
    check("E  |Z| = 4: every core with that form has a cyclic triangle inside Z", n4 == c4, (n4, c4))

print("ALL CHECKS PASSED (%d checks)" % len(RES) if all(RES) else "SOME CHECKS FAILED")
raise SystemExit(0 if all(RES) else 1)
