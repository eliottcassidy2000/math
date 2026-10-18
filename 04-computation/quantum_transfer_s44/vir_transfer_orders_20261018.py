#!/usr/bin/env python3
"""Transfer-order filtration of ker pi for Virasoro minimal models, all orders (S44).

For each weight D, Koszul 1-cycles are processed in increasing length r.  For a cycle z, the element Phi(z) in I is
reduced top-down: at length s = r - j its component t is reduced modulo the span W_s of tops of elements already
recorded at length s (initially the J-generators Phi(y), y in K_1(D,s), whose tops span J_s); if t lies in W_s the
corresponding combination is subtracted and the next length is examined; otherwise t is a new vector of grI at length
s, recorded with transfer order j (and the reduced element is added to the list at length s).  Every vector of ker pi
arises this way (descending induction on naive length), so the recorded orders give dims of the successive pieces of
ker pi_(D,s) by order; their sum must equal dim grI_(D,s) - dim J_(D,s) (checked).
Also reports the C[x][d]-span of the first-order class b_N mod J against the order-1 dimension.
Usage: python vir_transfer_orders_20261018.py p q Dmax [--prime P]
Session opus-2026-10-18-S44.
"""
import os, sys, time, argparse, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'voa_singular_support'))
sys.path.insert(0, HERE)
import flint
from vir_ss import block, block_index, central_charge, mult_x, deriv
from vir_quantum_transfer_20261018 import Transfer, by_length, add_into, left_nullspace, to_rows, rank


class Reducer:
    """Incremental row-echelon basis of tops at one (D, s), remembering for each basis row the element of V_D it is the
    top of (as an explicit vector), so a top can be cancelled by subtracting elements."""
    def __init__(self, D, s, P):
        self.D, self.s, self.P = D, s, P
        self.idx = block_index(D, s)
        self.cols = block(D, s)
        self.rows = []      # (pivot, row list, element dict)

    def reduce(self, t, elem):
        """Reduce top t (dict over block(D,s)) together with element elem; returns (t_reduced_row, elem_reduced)."""
        P = self.P
        row = [0] * len(self.cols)
        for lam, co in t.items():
            row[self.idx[lam]] = (row[self.idx[lam]] + co) % P
        e = dict(elem)
        for piv, brow, belem in self.rows:
            if row[piv]:
                f = row[piv]
                row = [(a - f * b) % P for a, b in zip(row, brow)]
                add_into(e, belem, -f, P)
        return row, {k: v for k, v in e.items() if v}

    def add(self, row, elem):
        P = self.P
        piv = next(j for j, x in enumerate(row) if x)
        inv = pow(row[piv], -1, P)
        row = [x * inv % P for x in row]
        elem = {k: v * inv % P for k, v in elem.items()}
        # keep echelon: eliminate this pivot from existing rows (not required for correctness)
        self.rows.append((piv, row, elem))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('p', type=int)
    ap.add_argument('q', type=int)
    ap.add_argument('Dmax', type=int)
    ap.add_argument('--prime', type=int, default=2147483647)
    a = ap.parse_args()
    P = a.prime
    T = Transfer(a.p, a.q, a.Dmax, P)
    N = T.N
    ss = T.ss
    print('=== Vir_(%d,%d), c = %s, N = %d: transfer-order filtration of ker pi' % (a.p, a.q, central_charge(a.p, a.q), N))
    report = []
    bN_span = {}   # (D,s) -> list of dicts: C[x][d]-span of the first-order class
    for D in range(2 * N, a.Dmax + 1):
        t0 = time.time()
        red = {}
        def R(s):
            if s not in red:
                red[s] = Reducer(D, s, P)
            return red[s]
        # order-0 elements: Phi(y) for y in a basis of K_1(D, s) (tops = d_1(y) span J_s)
        for s in range(N, D // 2 + 1):
            B = T.k1_basis(D, s)
            for (n, m) in B:
                el = T.phi_basis(n, m)
                top = by_length(el).get(s, {})
                if not top:
                    continue
                row, e = R(s).reduce(top, el)
                if any(row):
                    R(s).add(row, e)
        orders = {}    # s -> {order: count}
        for r in range(N + 1, D // 2 + 1):
            B1 = T.k1_basis(D, r)
            if not B1:
                continue
            rows1, nc = T.d1_rows(B1, D, r)
            cyc = left_nullspace(rows1, nc, P)
            for z in cyc:
                el = T.phi_of(z, B1)
                assert not by_length(el).get(r), 'top did not cancel'
                j = 1
                while r - j >= 0:
                    s = r - j
                    t = by_length(el).get(s, {})
                    if t:
                        row, el2 = R(s).reduce(t, el)
                        if any(row):
                            R(s).add(row, el2)
                            orders.setdefault(s, {})
                            orders[s][j] = orders[s].get(j, 0) + 1
                            break
                        el = el2
                        assert not by_length(el).get(s), 'reduction failed'
                    j += 1
        # compare with ker pi
        for s in range(0, D // 2 + 1):
            g = ss.grI.get((D, s), [])
            J = ss.J.get((D, s), [])
            ker = len(g) - len(J)
            od = orders.get(s, {})
            tot = sum(od.values())
            if ker or tot:
                report.append(dict(D=D, s=s, li=D - 2 * s, ker=ker, orders=od))
                if tot != ker:
                    print('  MISMATCH at D=%d s=%d: ker %d, recorded %d' % (D, s, ker, tot))
        line = '; '.join('s=%d(Li %d) ker %d %s' % (x['s'], x['li'], x['ker'], dict(sorted(x['orders'].items())))
                         for x in report if x['D'] == D)
        print('D=%2d %s (%.1fs)' % (D, line if line else '-', time.time() - t0), flush=True)
    # totals by order
    tot = {}
    for x in report:
        for k, v in x['orders'].items():
            tot[k] = tot.get(k, 0) + v
    print('TOTAL by order through D=%d: %s; dim ker pi = %d' % (a.Dmax, dict(sorted(tot.items())), sum(x['ker'] for x in report)))
    json.dump(dict(p=a.p, q=a.q, Dmax=a.Dmax, report=[dict(x, orders={str(k): v for k, v in x['orders'].items()}) for x in report]),
              open(os.path.join(HERE, 'vir_transfer_orders_%d_%d_D%d.json' % (a.p, a.q, a.Dmax)), 'w'), indent=1)


if __name__ == '__main__':
    main()
