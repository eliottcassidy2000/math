#!/usr/bin/env python3
"""The quantum transfer from the derived arc space to the classical limit of a Virasoro minimal model (S44).

Builds on kind-pasteur-2026-10-18-S1 (singular_support_arc_spaces_classical_freeness_20261018.md, section 1.3 and open
problem 6) and reuses its engine vir_ss.py (imported, unchanged).

Setting.  V^c the universal Virasoro VOA at c = c_(p,q), v its singular vector of weight 2N, N = (p-1)(q-1)/2, I the
maximal ideal, gr_F V^c = C[x_0, x_1, ...] (x_j <-> L_(-j-2)), J = differential ideal of x_0^N (the arc ideal),
ker pi = grI / J.  K = Koszul complex of f_n = d^n x_0^N on C[x] (the derived arc space of Spec C[x]/(x^N)).

Transfer.  For a Koszul 1-cycle z = sum_n g_n e_n of bidegree (weight D, length r), i.e. sum g_n f_n = 0, put
    Phi(z) = sum_n PBW(g_n) L_(-1)^n v   in I_D,
where PBW(g) is the ordered PBW lift.  Its length-r component is sum g_n f_n = 0, so Phi(z) has length <= r-1.  The
first-order transfer T_1(z) is its length-(r-1) component; it lies in grI_(D,r-1) and is well defined modulo J (changing
lifts or ordering changes it by elements of J).  If T_1(z) lies in J, write T_1(z) = d_1(z') and continue:
T_2(z) = length-(r-2) component of Phi(z) - Phi(z'), etc.  Every element of ker pi arises as some T_k (descending
induction on the naive length), so ker pi is filtered by transfer order; the question is how much order 1 already gives.

Output per bidegree (D, s): dim ker pi_(D,s), dim H_1(D, s+1), rank of T_1 (from cycles of length s+1) modulo J, and,
when T_1 does not fill ker pi, the extra rank from T_2 (cycles of length s+2).  Also checks that Koszul boundaries
transfer into J (so T_1 factors through H_1) and that every T_1 lies in grI.

Usage: python vir_quantum_transfer_20261018.py p q Dmax [--prime P]
Session opus-2026-10-18-S44.
"""
import os, sys, time, argparse, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'voa_singular_support'))
import flint
from vir_ss import SingularSupport, block, block_index, full_order, apply_neg, deriv, central_charge


def add_into(acc, vec, s, P):
    for k, v in vec.items():
        acc[k] = (acc.get(k, 0) + s * v) % P


def apply_L(m, vec, P):
    """L_(-m) applied to a vector {PBW mono: coeff}."""
    out = {}
    for mono, co in vec.items():
        for mo2, co2 in apply_neg(m, mono):
            out[mo2] = (out.get(mo2, 0) + co * co2) % P
    return {k: v for k, v in out.items() if v}


def by_length(vec):
    out = {}
    for mono, co in vec.items():
        out.setdefault(len(mono), {})[mono] = co
    return out


def to_rows(vecs, D, r, P):
    idx = block_index(D, r)
    n = len(idx)
    rows = []
    for v in vecs:
        row = [0] * n
        for lam, co in v.items():
            row[idx[lam]] = (row[idx[lam]] + co) % P
        rows.append(row)
    return rows, n


def rank(rows, n, P):
    if not rows or not n:
        return 0
    return flint.nmod_mat(len(rows), n, [x for r in rows for x in r], P).rank()


def left_nullspace(rows, ncols, P):
    """Vectors c (len = #rows) with sum_i c_i rows[i] = 0, as a list of coefficient lists."""
    m = len(rows)
    if m == 0:
        return []
    if ncols == 0:
        return [[1 if i == j else 0 for i in range(m)] for j in range(m)]
    # A^T x = 0 with A^T of size ncols x m
    AT = flint.nmod_mat(ncols, m, [rows[i][j] for j in range(ncols) for i in range(m)], P)
    X, nul = AT.nullspace()
    return [[int(X[i, k]) for i in range(m)] for k in range(nul)]


class Transfer:
    def __init__(self, p, q, Dmax, P):
        self.p, self.q, self.Dmax, self.P = p, q, Dmax, P
        self.ss = SingularSupport(p, q, Dmax, P)
        self.ss.run(verbose=False)
        self.N = self.ss.N
        v = self.ss.vec_mod(self.ss.v)
        self.u = [v]                         # u_n = L_(-1)^n v
        self.f = [{(2,) * self.N: 1}]        # f_n = d^n x_0^N
        while 2 * self.N + len(self.u) - 1 < Dmax:
            self.u.append(apply_L(1, self.u[-1], P))
            self.f.append({k: x % P for k, x in deriv(self.f[-1]).items()})
        # check sigma(u_n) = f_n
        for n, (un, fn) in enumerate(zip(self.u, self.f)):
            top = by_length(un).get(self.N, {})
            assert {k: x for k, x in top.items() if x} == {k: x for k, x in fn.items() if x}, f"symbol of L_-1^{n} v"
        self.cache = {}

    def phi_basis(self, n, m):
        """PBW(m) applied to u_n, m a partition (k1 >= k2 >= ...): L_(-k1)(L_(-k2)(... u_n))."""
        key = (n, m)
        if key in self.cache:
            return self.cache[key]
        if not m:
            res = self.u[n]
        else:
            res = apply_L(m[0], self.phi_basis(n, m[1:]), self.P)
        self.cache[key] = res
        return res

    def k1_basis(self, D, r):
        out = []
        for n in range(len(self.f)):
            w = D - 2 * self.N - n
            if w < 0 or r - self.N < 0:
                continue
            for m in block(w, r - self.N):
                out.append((n, m))
        return out

    def d1_rows(self, basis, D, r):
        idx = block_index(D, r)
        rows = []
        for (n, m) in basis:
            row = [0] * len(idx)
            for mono, co in self.f[n].items():
                lam = tuple(sorted(mono + m, reverse=True))
                row[idx[lam]] = (row[idx[lam]] + co) % self.P
            rows.append(row)
        return rows, len(idx)

    def d2_rows(self, B1, D, r):
        """Koszul boundaries in K_1(D, r) coordinates."""
        i1 = {b: i for i, b in enumerate(B1)}
        rows = []
        nf = len(self.f)
        for n1 in range(nf):
            for n2 in range(n1 + 1, nf):
                w = D - 4 * self.N - n1 - n2
                if w < 0 or r - 2 * self.N < 0:
                    continue
                for m in block(w, r - 2 * self.N):
                    row = [0] * len(B1)
                    for mono, co in self.f[n1].items():
                        lam = tuple(sorted(mono + m, reverse=True))
                        j = i1[(n2, lam)]
                        row[j] = (row[j] + co) % self.P
                    for mono, co in self.f[n2].items():
                        lam = tuple(sorted(mono + m, reverse=True))
                        j = i1[(n1, lam)]
                        row[j] = (row[j] - co) % self.P
                    rows.append(row)
        return rows

    def phi_of(self, coeffs, basis):
        acc = {}
        for c, (n, m) in zip(coeffs, basis):
            if c:
                add_into(acc, self.phi_basis(n, m), c, self.P)
        return {k: v for k, v in acc.items() if v}

    def run(self):
        P = self.P
        ss = self.ss
        report = []
        for D in range(2 * self.N, self.Dmax + 1):
            t0 = time.time()
            # per length s: collect transfer images
            img1 = {}          # s -> list of T_1 vectors (block (D, s))
            img2 = {}          # s -> list of T_2 vectors
            h1 = {}
            bnd_ok = True
            in_grI = True
            for r in range(self.N, D // 2 + 1):
                B1 = self.k1_basis(D, r)
                if not B1:
                    continue
                rows1, nc = self.d1_rows(B1, D, r)
                cyc = left_nullspace(rows1, nc, P)
                rows2 = self.d2_rows(B1, D, r)
                rkb = rank(rows2, len(B1), P)
                h1[r] = len(cyc) - rkb
                if not cyc:
                    continue
                s = r - 1
                Js = ss.J.get((D, s), [])
                T1 = []
                lower = []
                for z in cyc:
                    ph = self.phi_of(z, B1)
                    bl = by_length(ph)
                    assert not bl.get(r), f"top length did not cancel at {(D, r)}"
                    t1 = bl.get(s, {})
                    T1.append(t1)
                    lower.append((z, ph, t1))
                img1.setdefault(s, []).extend(T1)
                # boundaries -> J ?
                if rows2:
                    for brow in rows2[:min(len(rows2), 25)]:
                        ph = self.phi_of(brow, B1)
                        t1 = by_length(ph).get(s, {})
                        if t1:
                            rJ, n = to_rows(Js, D, s, P)
                            rt, _ = to_rows([t1], D, s, P)
                            if rank(rJ + rt, n, P) != rank(rJ, n, P):
                                bnd_ok = False
                # T_1 inside grI ?
                g = ss.grI.get((D, s), [])
                rg, n = to_rows(g, D, s, P)
                rt, _ = to_rows([t for t in T1 if t], D, s, P)
                if rt and rank(rg + rt, n, P) != rank(rg, n, P):
                    in_grI = False
                # second order: cycles whose T_1 lies in J
                if s - 1 >= 0 and Js is not None:
                    B1s = self.k1_basis(D, s)
                    if B1s:
                        rows1s, ncs = self.d1_rows(B1s, D, s)
                        # solve: find combinations of T_1 lying in span(d1 rows of K_1(D,s)) = J_(D,s)
                        rT, n = to_rows(T1, D, s, P)
                        # left nullspace of [rows1s ; rT] restricted: c_T * T1 = -c_B * rows1s
                        stack = rows1s + rT
                        ns = left_nullspace(stack, n, P)
                        T2 = []
                        for vec in ns:
                            cB = vec[:len(rows1s)]
                            cT = vec[len(rows1s):]
                            if not any(cT):
                                continue
                            # z_comb = sum cT_i z_i ; z' = -sum cB_j (basis of K_1(D,s))
                            ph = {}
                            for c, (z, phz, t1) in zip(cT, lower):
                                if c:
                                    add_into(ph, phz, c, P)
                            phz2 = self.phi_of(cB, B1s)   # Phi of sum cB_j e_j with d1 = sum cB_j rows1s_j = -sum cT T1
                            add_into(ph, phz2, 1, P)
                            bl = by_length({k: x for k, x in ph.items() if x})
                            assert not bl.get(s), f"second-order cancellation failed at {(D, r)}"
                            T2.append(bl.get(s - 1, {}))
                        img2.setdefault(s - 1, []).extend(T2)
            # assemble per s
            for s in range(0, D // 2 + 1):
                g = ss.grI.get((D, s), [])
                Js = ss.J.get((D, s), [])
                ker = len(g) - len(Js)
                if ker == 0 and not img1.get(s) and not img2.get(s):
                    continue
                rJ, n = to_rows(Js, D, s, P)
                rJ_rank = rank(rJ, n, P)
                r1, _ = to_rows([t for t in img1.get(s, []) if t], D, s, P)
                rk1 = rank(rJ + r1, n, P) - rJ_rank if n else 0
                r2, _ = to_rows([t for t in img2.get(s, []) if t], D, s, P)
                rk12 = rank(rJ + r1 + r2, n, P) - rJ_rank if n else 0
                report.append(dict(D=D, s=s, li=D - 2 * s, ker=ker, H1_from=h1.get(s + 1, 0), T1=rk1, T1T2=rk12))
            line = '; '.join('s=%d(Li %d): ker %d, H1 %d, T1 %d, T1+T2 %d' % (x['s'], x['li'], x['ker'], x['H1_from'],
                                                                         x['T1'], x['T1T2'])
                             for x in report if x['D'] == D and (x['ker'] or x['T1'] or x['T1T2']))
            print('D=%2d %s  [boundaries->J: %s, T1 in grI: %s] (%.1fs)' % (D, line if line else '-', bnd_ok, in_grI,
                                                                           time.time() - t0), flush=True)
        return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('p', type=int)
    ap.add_argument('q', type=int)
    ap.add_argument('Dmax', type=int)
    ap.add_argument('--prime', type=int, default=2147483647)
    a = ap.parse_args()
    print('=== Vir_(%d,%d), c = %s, N = %d, prime %d' % (a.p, a.q, central_charge(a.p, a.q), (a.p - 1) * (a.q - 1) // 2, a.prime))
    T = Transfer(a.p, a.q, a.Dmax, a.prime)
    rep = T.run()
    tot = {}
    for x in rep:
        tot['ker'] = tot.get('ker', 0) + x['ker']
        tot['T1'] = tot.get('T1', 0) + x['T1']
        tot['T1T2'] = tot.get('T1T2', 0) + x['T1T2']
    print('TOTAL through D=%d: dim ker pi = %d, rank T1 = %d, rank T1+T2 = %d' % (a.Dmax, tot.get('ker', 0), tot.get('T1', 0), tot.get('T1T2', 0)))
    out = os.path.join(HERE, 'vir_quantum_transfer_%d_%d_D%d.json' % (a.p, a.q, a.Dmax))
    json.dump(dict(p=a.p, q=a.q, Dmax=a.Dmax, prime=a.prime, report=rep), open(out, 'w'), indent=1)


if __name__ == '__main__':
    main()
