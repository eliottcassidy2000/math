#!/usr/bin/env python3
"""Exact THM-4247 audit of the THM-4241 iota-projection degrees.

We enumerate the normalized full Hom lattice of THM-4241 through degree 42
and record, for full vectors of degree 34 or 42,
the degrees of v=m+iota*m and ell=m-iota*m.  The identity

    q(v)+q(ell)=4q(m)

is checked vector by vector.  Eisenstein integers are pairs a+b*omega with
omega^2+omega+1=0.
"""

from collections import Counter, defaultdict


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def mul(x, y):
    a, b = x
    c, d = y
    return a * c - b * d, a * d + b * c - b * d


def conj(x):
    a, b = x
    return a - b, -b


def norm(x):
    return mul(x, conj(x))[0]


OMEGA2 = (-1, -1)
H = (
    ((6, 0), (-4, -2)),
    ((-2, 2), (6, 0)),
)


def hidden_degree(a, b):
    x = (a, b)
    total = (0, 0)
    for i in range(2):
        for j in range(2):
            total = add(total, mul(mul(x[i], H[i][j]), conj(x[j])))
    assert total[1] == 0
    return total[0]


def residue2(x):
    return x[0] % 2, x[1] % 2


def enumerate_rows():
    # The maintained certificate proves these boxes do not touch their
    # coordinate boundary.  We recheck that independently here.
    eis = []
    for a in range(-8, 9):
        for b in range(-8, 9):
            x = (a, b)
            if norm(x) <= 42:
                eis.append((x, norm(x)))
    assert all(max(abs(x[0]), abs(x[1])) < 8 for x, _ in eis)

    hidden = []
    for am in range(-12, 13):
        for an in range(-12, 13):
            for bm in range(-12, 13):
                for bn in range(-12, 13):
                    a, b = (am, an), (bm, bn)
                    q = hidden_degree(a, b)
                    if q <= 168:
                        hidden.append((a, b, q))
    assert all(
        max(abs(a[0]), abs(a[1]), abs(b[0]), abs(b[1])) < 12
        for a, b, _ in hidden
    )
    by_parity = defaultdict(list)
    for a, b, q in hidden:
        by_parity[residue2(a) + residue2(b)].append((a, b, q))
    return eis, by_parity


def main():
    eis, by_parity = enumerate_rows()
    full_theta = Counter()
    projection_hist = {34: Counter(), 42: Counter()}
    projection_pairs = {34: Counter(), 42: Counter()}

    for d, nd in eis:  # coefficient of the mixed glue basis vector h
        parity = residue2(mul(OMEGA2, d)) + residue2(d)
        for ell_a, ell_b, qell in by_parity[parity]:
            base = nd + qell // 4
            assert qell % 4 == 0
            if base > 42:
                continue
            for a, na in eis:  # coefficient of the visible basis vector u
                qfull = base + 4 * na
                if qfull > 42:
                    continue
                qvis = 4 * nd + 16 * na
                assert qvis + qell == 4 * qfull
                full_theta[qfull] += 1
                if qfull in projection_hist:
                    projection_hist[qfull][qell] += 1
                    projection_pairs[qfull][(qvis, qell)] += 1

    expected = {34: 36288, 42: 16992}
    for degree, count in expected.items():
        assert full_theta[degree] == count
        assert sum(projection_hist[degree].values()) == count
        assert all(q % 12 == 0 and q > 0 for q in projection_hist[degree])

    eliminated = {
        34: Counter({12: projection_hist[34][12], 132: projection_hist[34][132]}),
        42: Counter({12: projection_hist[42][12], 168: projection_hist[42][168]}),
    }
    remaining = {
        degree: projection_hist[degree] - eliminated[degree]
        for degree in (34, 42)
    }
    assert sum(eliminated[34].values()) == 2112
    assert sum(remaining[34].values()) == 34176
    assert sum(eliminated[42].values()) == 864
    assert sum(remaining[42].values()) == 16128

    print("THM-4247 W=0 IOTA PROJECTION HISTOGRAM")
    for degree in (34, 42):
        print(f"full_degree={degree} full_count={full_theta[degree]}")
        print("hidden_projection_degree:count")
        for qell, count in sorted(projection_hist[degree].items()):
            print(f"  {qell}:{count}")
        print("(visible_projection_degree,hidden_projection_degree):count")
        for pair, count in sorted(projection_pairs[degree].items()):
            print(f"  {pair}:{count}")
        print(f"eliminated_hidden_projection_hist={dict(sorted(eliminated[degree].items()))}")
        print(f"remaining_hidden_projection_hist={dict(sorted(remaining[degree].items()))}")
        print(f"remaining_raw_vectors={sum(remaining[degree].values())}")


if __name__ == "__main__":
    main()
