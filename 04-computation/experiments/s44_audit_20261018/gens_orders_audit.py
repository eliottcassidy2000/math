#!/usr/bin/env python3
"""S44 audit, HYP-9315(b)(c) and kind-pasteur's HYP-9294(a): minimal differential-ideal generators of grI and THEIR
transfer orders, defined invariantly:
    Dec(D,s) = sum_j x_j grI(D-j-2, s-1) + d grI(D-1, s)          (decomposables)
    #generators of order j at (D,s) = dim(F_j + Dec) - dim(F_(j-1) + Dec),   F_j = transfer-order filtration of grI
(F_0 = J).  Order-0 generators are the generators lying in J + Dec (only x_0^N at (2N, N)).
Own engine (orders_audit.Model).  Usage: python gens_orders_audit.py p q Dmax [P]
"""
import sys, time, json
import flint
from vcore import parts, basis_by_length, index_of, block, block_index, dpoly, asc
from orders_audit import Model


def rref_basis(rows, n, P):
    if not rows or not n:
        return []
    R, rk = flint.nmod_mat(len(rows), n, [x for r in rows for x in r], P).rref()
    return [list(map(int, r)) for r in R.tolist()[:rk]]


def rank(rows, n, P):
    if not rows or not n:
        return 0
    return flint.nmod_mat(len(rows), n, [x for r in rows for x in r], P).rank()


def main():
    p, q, Dmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    P = int(sys.argv[4]) if len(sys.argv) > 4 else 2147483647
    M = Model(p, q, Dmax, P)
    N = M.N
    grI = {}      # (D, s) -> list of dict polys (basis)
    gens = []
    t0 = time.time()
    for D in range(2 * N, Dmax + 1):
        basis = basis_by_length(D)
        idx = index_of(D)
        lens = [len(t) for t in basis]
        prs = M.pairs(D)
        rmax = max(N + len(m) for n, m in prs)
        Fb = {}   # (r, s) -> list of top polys
        for r in range(N, rmax + 1):
            sel = [(n, m) for (n, m) in prs if N + len(m) <= r]
            data = [0] * (len(sel) * len(basis))
            for i, (n, m) in enumerate(sel):
                off = i * len(basis)
                for mono, co in M.w(n, m).items():
                    data[off + idx[mono]] = co
            R, rk = flint.nmod_mat(len(sel), len(basis), data, P).rref()
            for row in R.tolist()[:rk]:
                piv = next(j for j, x in enumerate(row) if int(x))
                s = lens[piv]
                top = {basis[j]: int(x) for j, x in enumerate(row) if int(x) and lens[j] == s}
                Fb.setdefault((r, s), []).append(top)
        for s in range(N, rmax + 1):
            bidx = block_index(D, s)
            nb = len(bidx)
            to_row = lambda poly: [0] * nb if not poly else _row(poly, bidx, P)
            G = [_row(t, bidx, P) for t in Fb.get((rmax, s), [])]
            grI[(D, s)] = [dict(t) for t in Fb.get((rmax, s), [])]
            if not G:
                continue
            # decomposables
            dec = []
            for j in range(0, D):
                for poly in grI.get((D - j - 2, s - 1), []):
                    dec.append(_row({asc(k + (j + 2,)): v for k, v in poly.items()}, bidx, P))
            for poly in grI.get((D - 1, s), []):
                dp = dpoly(poly, P)
                if dp:
                    dec.append(_row(dp, bidx, P))
            rdec = rank(dec, nb, P)
            ngen = len(G) - rdec
            if ngen <= 0:
                continue
            # orders of generators
            prev = rank(dec + [_row(t, bidx, P) for t in Fb.get((s, s), [])], nb, P)  # F_0 + Dec
            od = {}
            if prev - rdec:
                od[0] = prev - rdec
            for jj in range(1, rmax - s + 1):
                cur = rank(dec + [_row(t, bidx, P) for t in Fb.get((s + jj, s), [])], nb, P)
                if cur - prev:
                    od[jj] = cur - prev
                prev = cur
            assert sum(od.values()) == ngen
            gens.append(dict(D=D, s=s, li=D - 2 * s, n=ngen, orders=od))
            print('generator(s) at (D=%d, s=%d, Li %d): %d, transfer orders %s' % (D, s, D - 2 * s, ngen, od), flush=True)
    print('Vir_(%d,%d) N=%d through D=%d: generators %s  (%.1fs)' % (p, q, N, Dmax, [(g['D'], g['s']) + ((g['n'],) if g['n'] > 1 else ()) for g in gens], time.time() - t0))
    json.dump(dict(p=p, q=q, Dmax=Dmax, P=P, gens=[dict(g, orders={str(k): v for k, v in g['orders'].items()}) for g in gens]),
              open('gens_orders_%d_%d_D%d.json' % (p, q, Dmax), 'w'), indent=1)


def _row(poly, bidx, P):
    row = [0] * len(bidx)
    for k, v in poly.items():
        row[bidx[k]] = (row[bidx[k]] + v) % P
    return row


if __name__ == '__main__':
    main()
