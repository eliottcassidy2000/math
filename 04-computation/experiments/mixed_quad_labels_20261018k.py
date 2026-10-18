"""THM-4700 mixed-distance quadrilaterals, replayed as exact copy labels.

No colouring claim: this supplies four points and checks all four side norms.
The incoming proof supplies the conic point; this implementation is deliberately
separate from unit_path_ring_20261018/mixed_quadrilaterals.py.
"""
from fractions import Fraction as Q


def mixed_quad(geo, k):
    f, bar = geo.f, geo.bar
    s3 = f.add(f.scale(f.power(f.gen(0), 5), 2), f.one)
    assert f.mul(s3, s3) == f.scalar(-3)
    norm = lambda x: f.mul(x, bar(x))
    N = norm(k)
    A = f.scale(f.add(k, bar(k)), Q(1, 2))
    B = f.mul(f.scale(f.sub(k, bar(k)), Q(1, 2)), f.inv(s3))
    assert bar(A) == A and bar(B) == B and B != f.zero
    t = f.mul(f.add(f.one, A), f.inv(f.scale(B, 3)))
    sigma = f.mul(f.add(A, N), f.inv(f.scale(B, 3)))
    assert f.sub(f.power(sigma, 2), f.mul(N, f.power(t, 2))) == f.scale(f.sub(N, f.one), Q(1, 3))
    den = f.add(f.one, f.scale(f.power(t, 2), 3))
    z = f.scale(f.inv(f.mul(bar(k), f.add(f.one, f.mul(t, s3)))), 2)
    def split(x, s):
        assert f.power(s, 2) == f.scale(f.mul(norm(x), f.sub(f.scalar(4), norm(x))), Q(1, 3))
        r = f.mul(f.mul(s3, s), f.inv(norm(x)))
        return [f.scale(f.mul(x, f.add(f.one, r)), Q(1, 2)),
                f.scale(f.mul(x, f.sub(f.one, r)), Q(1, 2))]
    u1, u2 = split(z, f.mul(f.scale(sigma, 4), f.inv(f.mul(N, den))))
    e1, e2 = split(f.scale(f.mul(k, z), -1), f.mul(f.scale(t, 4), f.inv(den)))
    assert all(norm(x) == f.one for x in (u1, u2, e1, e2))
    points = [f.zero, u1, z, f.add(z, f.mul(e1, f.inv(k)))]
    assert len(set(points)) == 4
    expected = [f.one, f.one, f.inv(N), f.inv(N)]
    assert [norm(f.sub(points[(i+1) % 4], points[i])) for i in range(4)] == expected
    return points
