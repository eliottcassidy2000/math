#!/usr/bin/env python3
"""Independent referee engine (audit2, 2026-10-18).  Written from scratch; imports nothing from the author's code.

Vacuum module of a mode (super)algebra with PBW basis, computed mod a prime P:
  * states = dict {sorted tuple of creation modes: coeff mod P}; a mode is (t, n2) with n2 = 2 * (mode index);
  * apply(X, mono) by recursive (super)commutation, memoised;
  * contravariant (Shapovalov) form via <X u, v> = <u, omega(X) v>, block-diagonal by (weight, charge);
  * maximal submodule I = radical of the form (no generating-set assumption: every leftmost creation mode is used);
  * Li filtration = PBW filtration; gr I by Li-ordered echelon form; arc ideal J = (gr I in Li 0)_d inside the
    supercommutative jet algebra; minimal generators of gr I as a differential ideal by graded Nakayama.
"""
import itertools
from fractions import Fraction
import flint

P_DEFAULT = 2147483647


class ModeVA:
    """Subclasses define: gens = list of (name, H2 = 2*weight, parity, charge); omega(mode); bracket_apply(X, Y, rest)."""

    def __init__(self, P=P_DEFAULT):
        self.P = P
        self.memo = {}
        self.inv2 = pow(2, P - 2, P)

    # ---------- basic data ----------
    def F(self, x):
        x = Fraction(x)
        return x.numerator % self.P * pow(x.denominator % self.P, self.P - 2, self.P) % self.P

    def H2(self, t):
        return self.gens[t][1]

    def par(self, X):
        return self.gens[X[0]][2]

    def chg(self, X):
        return self.gens[X[0]][3]

    def creation(self, X):
        return X[1] <= -self.H2(X[0])

    def li(self, X):  # Li degree of a creation mode: x_{t,j} <-> u^t_{-h_t-j}
        return (-X[1] - self.H2(X[0])) // 2

    def wt2(self, mono):
        return sum(-X[1] for X in mono)

    def mono_chg(self, mono):
        return sum(self.chg(X) for X in mono)

    def mono_li(self, mono):
        return sum(self.li(X) for X in mono)

    # ---------- module action ----------
    def apply(self, X, mono):
        key = (X, mono)
        r = self.memo.get(key)
        if r is not None:
            return r
        P = self.P
        if not mono:
            res = {(X,): 1} if self.creation(X) else {}
        elif self.creation(X) and X < mono[0]:
            res = {(X,) + mono: 1}
        elif self.creation(X) and X == mono[0] and self.par(X) == 0:
            res = {(X,) + mono: 1}
        else:
            Y, rest = mono[0], mono[1:]
            if X == Y:  # odd creation mode squared: X X rest = 1/2 {X,X} rest
                br = self.bracket_apply(X, X, rest)
                res = {m: c * self.inv2 % P for m, c in br.items() if c}
            else:
                sign = P - 1 if (self.par(X) and self.par(Y)) else 1
                inner = self.apply(X, rest)
                res = {}
                for m, c in self.apply_vec(Y, inner).items():
                    res[m] = (res.get(m, 0) + c * sign) % P
                for m, c in self.bracket_apply(X, Y, rest).items():
                    res[m] = (res.get(m, 0) + c) % P
                res = {m: c for m, c in res.items() if c}
        self.memo[key] = res
        return res

    def apply_vec(self, X, vec):
        P = self.P
        out = {}
        for m, c in vec.items():
            for m2, c2 in self.apply(X, m).items():
                out[m2] = (out.get(m2, 0) + c * c2) % P
        return {m: c for m, c in out.items() if c}

    def apply_word(self, word, vec):
        for X in reversed(word):
            vec = self.apply_vec(X, vec)
        return vec

    def lie_bracket_apply(self, terms, central, rest):
        """helper: sum_k co_k * apply(mode_k, rest) + central * rest."""
        P = self.P
        out = {}
        for (Z, co) in terms:
            if co % P == 0:
                continue
            for m, c in self.apply(Z, rest).items():
                out[m] = (out.get(m, 0) + c * co) % P
        if central % P:
            out[rest] = (out.get(rest, 0) + central) % P
        return {m: c for m, c in out.items() if c}

    # ---------- PBW basis ----------
    def creation_modes_upto(self, W2max):
        modes = []
        for t, g in enumerate(self.gens):
            n2 = -g[1]
            while -n2 <= W2max:
                modes.append((t, n2))
                n2 -= 2
        return sorted(modes)

    def basis(self, W2):
        """sorted monomials of weight W2 (doubled)"""
        if not hasattr(self, '_basis_cache'):
            self._basis_cache = {}
        if W2 in self._basis_cache:
            return self._basis_cache[W2]
        modes = self.creation_modes_upto(W2)
        out = []

        def rec(i, rem, cur):
            if rem == 0:
                out.append(tuple(cur))
                return
            for k in range(i, len(modes)):
                X = modes[k]
                w = -X[1]
                if w > rem:
                    continue
                if cur and X == cur[-1] and self.par(X):
                    continue
                cur.append(X)
                rec(k, rem - w, cur)
                cur.pop()
        rec(0, W2, [])
        self._basis_cache[W2] = out
        return out

    def blocks(self, W2):
        d = {}
        for m in self.basis(W2):
            d.setdefault(self.mono_chg(m), []).append(m)
        return d

    # ---------- contravariant form ----------
    def gram_all(self, W2max, verbose=False):
        """G[(W2,chg)] = (basis list, flint matrix) for W2 = 0..W2max (only those weights with states)"""
        P = self.P
        G = {(0, 0): ([()], flint.nmod_mat([[1]], P))}
        idx = {(0, 0): {(): 0}}
        for W2 in range(1, W2max + 1):
            for ch, B in sorted(self.blocks(W2).items()):
                n = len(B)
                Gm = [[0] * n for _ in range(n)]
                # group rows by leftmost mode
                groups = {}
                for i, b in enumerate(B):
                    groups.setdefault(b[0], []).append(i)
                for X, rows in groups.items():
                    oX = self.omega(X)
                    lowW2 = W2 + X[1]   # weight of u (X[1] negative)
                    lowch = ch - self.chg(X)
                    Bl, Gl = G[(lowW2, lowch)]
                    il = idx[(lowW2, lowch)]
                    # A[j] = omega(X) b_j expressed in lower basis
                    cols = []
                    for b in B:
                        v = self.apply(oX, b)
                        cols.append(v)
                    nl = len(Bl)
                    A = flint.nmod_mat(nl, n, [0] * (nl * n), P)
                    for j, v in enumerate(cols):
                        for m, c in v.items():
                            A[il[m], j] = c
                    # rows u_i
                    U = flint.nmod_mat(len(rows), nl, [0] * (len(rows) * nl), P)
                    for r, i in enumerate(rows):
                        u = B[i][1:]
                        k = il[u]
                        for col in range(nl):
                            U[r, col] = Gl[k, col]
                    R = U * A
                    for r, i in enumerate(rows):
                        for j in range(n):
                            Gm[i][j] = int(R[r, j])
                G[(W2, ch)] = (B, flint.nmod_mat(Gm, P) if n else None)
                idx[(W2, ch)] = {m: i for i, m in enumerate(B)}
            if verbose:
                print('gram weight', Fraction(W2, 2), 'done', flush=True)
        return G

    def radical(self, G):
        """I[(W2,ch)] = list of vectors (dict) spanning the radical (right kernel)."""
        P = self.P
        I = {}
        for key, (B, Gm) in G.items():
            if Gm is None:
                I[key] = []
                continue
            X, nul = Gm.nullspace()
            vecs = []
            for k in range(nul):
                v = {B[i]: int(X[i, k]) for i in range(len(B)) if int(X[i, k])}
                vecs.append(v)
            I[key] = vecs
        return I


# ---------------- linear algebra helpers ----------------
def rowspace_basis(vectors, cols, P):
    """vectors: list of dicts over keys in cols (list). returns list of dicts forming an echelon basis (rref rows)."""
    if not vectors or not cols:
        return []
    ci = {m: i for i, m in enumerate(cols)}
    M = flint.nmod_mat(len(vectors), len(cols), [0] * (len(vectors) * len(cols)), P)
    for r, v in enumerate(vectors):
        for m, c in v.items():
            M[r, ci[m]] = c % P
    R, rk = M.rref()
    out = []
    for r in range(rk):
        out.append({cols[j]: int(R[r, j]) for j in range(len(cols)) if int(R[r, j])})
    return out


def rank_of(vectors, cols, P):
    if not vectors or not cols:
        return 0
    ci = {m: i for i, m in enumerate(cols)}
    M = flint.nmod_mat(len(vectors), len(cols), [0] * (len(vectors) * len(cols)), P)
    for r, v in enumerate(vectors):
        for m, c in v.items():
            M[r, ci[m]] = c % P
    return M.rank()


# ---------------- associated graded / arc space machinery ----------------
class GradedAnalysis:
    """Given a ModeVA and the radical I (by (W2,ch)), compute gr I, the arc ideal J and the minimal generators."""

    def __init__(self, va, I, W2max):
        self.va = va
        self.P = va.P
        self.W2max = W2max
        self.I = I
        self.vars = va.creation_modes_upto(W2max)

    # supercommutative multiplication: variable x times monomial
    def mul_var(self, x, mono):
        va = self.va
        # position
        pos = 0
        while pos < len(mono) and mono[pos] < x:
            pos += 1
        if va.par(x) and pos < len(mono) and mono[pos] == x:
            return None, 0
        sign = 1
        if va.par(x):
            nodd = sum(1 for y in mono[:pos] if va.par(y))
            if nodd % 2:
                sign = -1
        return mono[:pos] + (x,) + mono[pos:], sign

    def mul_var_vec(self, x, vec):
        P = self.P
        out = {}
        for m, c in vec.items():
            m2, s = self.mul_var(x, m)
            if m2 is None:
                continue
            out[m2] = (out.get(m2, 0) + c * s) % P
        return {m: c for m, c in out.items() if c}

    def sort_signed(self, lst):
        """sort a list of variables with Koszul sign; return (tuple, sign) or (None, 0) if odd repeated"""
        va = self.va
        a = list(lst)
        sign = 1
        # bubble sort
        n = len(a)
        for i in range(n):
            for j in range(n - 1 - i):
                if a[j] > a[j + 1]:
                    if va.par(a[j]) and va.par(a[j + 1]):
                        sign = -sign
                    a[j], a[j + 1] = a[j + 1], a[j]
        for j in range(n - 1):
            if a[j] == a[j + 1] and va.par(a[j]):
                return None, 0
        return tuple(a), sign

    def d_mono(self, mono):
        va = self.va
        P = self.P
        out = {}
        for i, X in enumerate(mono):
            co = va.li(X) + 1
            X2 = (X[0], X[1] - 2)
            lst = list(mono)
            lst[i] = X2
            m2, s = self.sort_signed(lst)
            if m2 is None:
                continue
            out[m2] = (out.get(m2, 0) + co * s) % P
        return {m: c for m, c in out.items() if c}

    def d_vec(self, vec):
        P = self.P
        out = {}
        for m, c in vec.items():
            for m2, c2 in self.d_mono(m).items():
                out[m2] = (out.get(m2, 0) + c * c2) % P
        return {m: c for m, c in out.items() if c}

    def mono_deg(self, m):
        va = self.va
        return (va.wt2(m), va.mono_chg(m), va.mono_li(m))

    def run(self, verbose=True):
        va, P = self.va, self.P
        # monomial lists per multidegree
        self.monos = {}
        for W2 in range(0, self.W2max + 1):
            for m in va.basis(W2):
                self.monos.setdefault(self.mono_deg(m), []).append(m)
        # gr I per multidegree
        grI = {}
        for (W2, ch), vecs in self.I.items():
            if not vecs:
                continue
            B = va.blocks(W2)[ch]
            order = sorted(B, key=lambda m: (va.mono_li(m), m))
            basisI = rowspace_basis(vecs, order, P)
            for v in basisI:
                lead = min(v.keys(), key=lambda m: (va.mono_li(m), order.index(m)))
                l = va.mono_li(lead)
                sym = {m: c for m, c in v.items() if va.mono_li(m) == l}
                grI.setdefault((W2, ch, l), []).append(sym)
        self.grI = {D: rowspace_basis(v, self.monos[D], P) for D, v in grI.items()}
        degs = sorted(self.monos.keys())
        J = {}
        Dec = {}
        self.table = []
        for D in degs:
            W2, ch, l = D
            cols = self.monos[D]
            candJ = []
            candDec = []
            if l == 0:
                candJ += self.grI.get(D, [])
            for x in self.vars:
                xd = self.mono_deg((x,))
                Dl = (W2 - xd[0], ch - xd[1], l - xd[2])
                if Dl[0] < 0 or Dl[2] < 0:
                    continue
                for b in J.get(Dl, []):
                    v = self.mul_var_vec(x, b)
                    if v:
                        candJ.append(v)
                for b in self.grI.get(Dl, []):
                    v = self.mul_var_vec(x, b)
                    if v:
                        candDec.append(v)
            Dd = (W2 - 2, ch, l - 1)
            for b in J.get(Dd, []):
                v = self.d_vec(b)
                if v:
                    candJ.append(v)
            for b in self.grI.get(Dd, []):
                v = self.d_vec(b)
                if v:
                    candDec.append(v)
            J[D] = rowspace_basis(candJ, cols, P)
            Dec[D] = rowspace_basis(candDec, cols, P)
            gI = self.grI.get(D, [])
            # sanity: J subset grI, Dec subset grI
            if J[D]:
                assert rank_of(gI + J[D], cols, P) == len(gI), ('J not in grI', D)
            if Dec[D]:
                assert rank_of(gI + Dec[D], cols, P) == len(gI), ('Dec not in grI', D)
            self.table.append((D, len(cols), len(gI), len(J[D]), len(gI) - len(Dec[D])))
        self.J = J
        self.Dec = Dec
        return self.table

    def summary_by_weight(self):
        out = {}
        for (D, nV, nI, nJ, ngen) in self.table:
            W2 = D[0]
            o = out.setdefault(W2, {'V': 0, 'I': 0, 'J': 0, 'ker': 0, 'gens': []})
            o['V'] += nV
            o['I'] += nI
            o['J'] += nJ
            o['ker'] += nI - nJ
            if ngen:
                o['gens'].append((D[2], D[1], ngen))
        return out
