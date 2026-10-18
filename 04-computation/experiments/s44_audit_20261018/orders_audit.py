#!/usr/bin/env python3
"""S44 audit: independent computation of the transfer-order filtration of ker pi = grI/J for Vir_(p,q).

Independent of vir_ss.py and of the S44 scripts:
  * own Virasoro engine (vcore.py), ASCENDING PBW order, ascending lift PBW(m) = L_-m1 (L_-m2 (... L_-mr u)) with
    m1 <= m2 <= ... (smallest mode outermost), whereas S44 uses descending order with the largest mode outermost;
  * a DIFFERENT algorithm: no incremental reduction of cycles.  For each weight D and each r,
        A_r = span{ PBW(m) L_-1^n v : N + len(m) <= r }  (= Phi of Koszul 1-chains of length <= r),
    and F_j(D,s) := tops at length s of A_(s+j) cap V^(<= s) is read off from one RREF of A_(s+j) with columns
    ordered by length (descending).  dim F_0(D,s) must equal dim J_(D,s) (computed separately from m * f_n);
    dim F_max(D,s) = dim grI_(D,s).  Order-j count = dim F_j - dim F_(j-1).
Outputs per bidegree: dim J, dim grI, ker, order counts.  Optionally compares with an S44 JSON.
Usage: python orders_audit.py p q Dmax [--prime P] [--compare path.json]
"""
import sys, time, json, argparse
import flint
from fractions import Fraction
from vcore import (parts, basis_by_length, index_of, block, block_index, apply_neg_vec, dpoly, polymul_mono,
                   central_charge)
from singvec import singular_vector_modp, singular_vector_exact


class Model:
    def __init__(self, p, q, Dmax, P):
        self.p, self.q, self.Dmax, self.P = p, q, Dmax, P
        self.N = N = (p - 1) * (q - 1) // 2
        self.c = central_charge(p, q)
        if 2 * N <= 20:
            ve = singular_vector_exact(p, q)
            self.v = {k: x.numerator % P * pow(x.denominator % P, -1, P) % P for k, x in ve.items()}
            self.v = {k: x for k, x in self.v.items() if x}
        else:
            self.v = singular_vector_modp(p, q, P)
        self.u = [self.v]
        while len(self.u) <= Dmax - 2 * N:
            self.u.append(apply_neg_vec(1, self.u[-1], P))
        self.f = [{(2,) * N: 1}]
        while len(self.f) <= Dmax - 2 * N:
            self.f.append(dpoly(self.f[-1], P))
        # symbol check: top-length part of u_n equals f_n
        for n, (un, fn) in enumerate(zip(self.u, self.f)):
            top = {k: x for k, x in un.items() if len(k) == N}
            assert top == fn, ('symbol of L_-1^n v', n)
            assert all(len(k) <= N for k in un), 'u_n has length > N'
        self._w = {}

    def w(self, n, m):
        """PBW_asc(m) u_n."""
        key = (n, m)
        if key not in self._w:
            if not m:
                self._w[key] = self.u[n]
            else:
                self._w[key] = apply_neg_vec(m[0], self.w(n, m[1:]), self.P)
        return self._w[key]

    def pairs(self, D):
        N = self.N
        out = []
        for n in range(0, D - 2 * N + 1):
            for m in parts(D - 2 * N - n):
                out.append((n, m))
        return out

    def J_dim(self, D, s):
        N, P = self.N, self.P
        rows = []
        idx = block_index(D, s)
        if not idx:
            return 0
        for n in range(0, D - 2 * N + 1):
            for m in parts(D - 2 * N - n):
                if len(m) != s - N:
                    continue
                row = [0] * len(idx)
                for mono, co in polymul_mono(self.f[n], m).items():
                    row[idx[mono]] = (row[idx[mono]] + co) % P
                rows.append(row)
        if not rows:
            return 0
        return flint.nmod_mat(len(rows), len(idx), [x for r in rows for x in r], P).rank()

    def filtration(self, D):
        """Returns dict r -> {s: #pivots at length s in RREF(A_r)}, plus rmax."""
        N, P = self.N, self.P
        basis = basis_by_length(D)
        idx = index_of(D)
        lens = [len(t) for t in basis]
        prs = self.pairs(D)
        if not prs:
            return {}, N - 1
        rmax = max(N + len(m) for n, m in prs)
        res = {}
        for r in range(N, rmax + 1):
            sel = [(n, m) for (n, m) in prs if N + len(m) <= r]
            if not sel:
                res[r] = {}
                continue
            data = [0] * (len(sel) * len(basis))
            for i, (n, m) in enumerate(sel):
                off = i * len(basis)
                for mono, co in self.w(n, m).items():
                    data[off + idx[mono]] = co
            M = flint.nmod_mat(len(sel), len(basis), data, P)
            R, rk = M.rref()
            cnt = {}
            for row in R.tolist()[:rk]:
                # first nonzero entry (pivot)
                for j, x in enumerate(row):
                    if int(x):
                        cnt[lens[j]] = cnt.get(lens[j], 0) + 1
                        break
            res[r] = cnt
        return res, rmax

    def rc_char(self, D):
        """dim of the simple quotient L_D from the Rocha-Caridi vacuum character of M(p,q)."""
        p, q = self.p, self.q
        num = {}
        K = 0
        while True:
            hit = False
            for k in ([0] if K == 0 else [K, -K]):
                e1 = p * q * k * k + k * (q - p)
                e2 = (p * k + 1) * (q * k + 1)
                if 0 <= e1 <= D:
                    num[e1] = num.get(e1, 0) + 1; hit = True
                if 0 <= e2 <= D:
                    num[e2] = num.get(e2, 0) - 1; hit = True
            if not hit and K > 3:
                break
            K += 1
        # divide by prod_{n>=1}(1-q^n): partition numbers
        part = [1] + [0] * D
        for n in range(1, D + 1):
            for m in range(n, D + 1):
                part[m] += part[m - n]
        return sum(num.get(i, 0) * part[D - i] for i in range(D + 1))

    def orders(self, D):
        N = self.N
        F, rmax = self.filtration(D)
        if F:
            dimI = sum(F[rmax].values())
            assert len(parts(D)) - dimI == self.rc_char(D), ('character check failed', D)
        else:
            assert len(parts(D)) == self.rc_char(D), ('character check failed', D)
        out = []
        for s in range(0, D // 2 + 1):
            if not F:
                break
            grI = F.get(rmax, {}).get(s, 0)
            Jd = self.J_dim(D, s) if s >= N else 0
            if s >= N and s <= rmax:
                F0 = F[s].get(s, 0)
                assert F0 == Jd, ('F_0 != J', D, s, F0, Jd)
            od = {}
            prev = Jd
            for j in range(1, rmax - s + 1):
                cur = F.get(s + j, {}).get(s, 0)
                if cur - prev:
                    od[j] = cur - prev
                prev = cur
            ker = grI - Jd
            assert sum(od.values()) == ker, ('order sum', D, s)
            if ker or grI:
                out.append(dict(D=D, s=s, li=D - 2 * s, J=Jd, grI=grI, ker=ker, orders=od))
        return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('p', type=int)
    ap.add_argument('q', type=int)
    ap.add_argument('Dmax', type=int)
    ap.add_argument('--prime', type=int, default=2147483647)
    ap.add_argument('--compare', default=None)
    ap.add_argument('--char', action='store_true', help='also check dim V - dim I against the Rocha-Caridi character')
    a = ap.parse_args()
    t0 = time.time()
    M = Model(a.p, a.q, a.Dmax, a.prime)
    print('=== audit Vir_(%d,%d) c=%s N=%d prime=%d (setup %.1fs)' % (a.p, a.q, M.c, M.N, a.prime, time.time() - t0), flush=True)
    allrows = []
    for D in range(2 * M.N, a.Dmax + 1):
        t1 = time.time()
        rows = M.orders(D)
        allrows.extend(rows)
        txt = '; '.join('s=%d(Li %d) ker %d %s' % (x['s'], x['li'], x['ker'], dict(sorted(x['orders'].items())))
                        for x in rows if x['ker'])
        print('D=%2d %s (%.1fs)' % (D, txt if txt else '-', time.time() - t1), flush=True)
    tot = {}
    for x in allrows:
        for k, v in x['orders'].items():
            tot[k] = tot.get(k, 0) + v
    print('TOTAL by order through D=%d: %s; dim ker pi = %d' % (a.Dmax, dict(sorted(tot.items())), sum(x['ker'] for x in allrows)))
    out = dict(p=a.p, q=a.q, Dmax=a.Dmax, prime=a.prime, rows=[dict(x, orders={str(k): v for k, v in x['orders'].items()}) for x in allrows])
    json.dump(out, open('audit_orders_%d_%d_D%d_P%d.json' % (a.p, a.q, a.Dmax, a.prime), 'w'), indent=1)
    if a.compare:
        S = json.load(open(a.compare))
        theirs = {(x['D'], x['s']): (x['ker'], {int(k): v for k, v in x['orders'].items()}) for x in S['report']}
        mine = {(x['D'], x['s']): (x['ker'], x['orders']) for x in allrows if x['ker']}
        keys = sorted(set(theirs) | set(mine))
        diff = [(k, theirs.get(k), mine.get(k)) for k in keys if k[0] <= min(S['Dmax'], a.Dmax) and theirs.get(k) != mine.get(k)]
        print('COMPARE with %s: %d bidegrees, %d differences' % (a.compare, len(keys), len(diff)))
        for d in diff[:20]:
            print('   DIFF', d)


if __name__ == '__main__':
    main()
