"""New Fano/binomial/endpoint composite. Run with Python 3."""
from pathlib import Path
import importlib.util
from math import comb
import json

ROOT = Path(__file__).parent
spec = importlib.util.spec_from_file_location("s265_fano", ROOT / "doubling/construct.py")
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


def tail(T, m):
    n = len(T)
    return tuple([row | sum(1 << j for j in range(n, n+m)) for row in T] +
                 [sum(1 << j for j in range(i+1, n+m)) for i in range(n, n+m)])


def wrap(T):
    n = len(T)
    return tuple([row | (1 << (n+1)) for row in T] + [(1 << n)-1, 1 << n])


def start_histogram(T, root, gate):
    """Forward endpoint DP times independent suffix-start DP at a fixed gate."""
    n = len(T); full = (1 << n)-1
    dp = [[0]*n for _ in range(1 << n)]
    dp[1 << root][root] = 1
    for S in range(1 << n):
        for v, value in enumerate(dp[S]):
            if not value: continue
            opts = T[v] & ~S
            while opts:
                b = opts & -opts; opts -= b
                dp[S | b][b.bit_length()-1] += value
    after = [c.hp(T, v) for v in range(n) if T[gate] >> v & 1]
    histogram = [0]*(n+1)
    for S in range(1 << n):
        if not dp[S][gate]: continue
        rest = full ^ S
        suffix = sum(h[rest] for h in after) if rest else 1
        histogram[S.bit_count()] += dp[S][gate]*suffix
    assert sum(histogram) == sum(dp[full])
    return histogram


def run():
    A = c.rebase(c.fano(0), 2); B = c.converse(A)
    rows = []
    for m in range(4):
        X, Y = tail(A, m), tail(B, m)
        expected = [0]*(len(X)+1)
        for j in range(m+1):
            expected[2+j] += 2*comb(m, j)
            expected[7+j] += comb(m, j)
        difference = [x-y for x,y in zip(c.rooted(X), c.rooted(Y))]
        assert difference == expected
        TX, TY = wrap(X), wrap(Y)
        hx, hy = c.hp(TX)[-1], c.hp(TY)[-1]
        assert hx == hy
        px = start_histogram(TX, 0, len(X)+1)
        py = start_histogram(TY, 0, len(Y)+1)
        assert [x-y for x,y in zip(px,py)] == [0] + expected + [0]
        assert sum(px) == 260*(1 << m) and sum(py) == 257*(1 << m)
        rows.append(dict(m=m, order=len(TX), H=hx, start_A=sum(px),
                         start_B=sum(py), start_gap=sum(px)-sum(py),
                         gate_histogram_difference=[x-y for x,y in zip(px,py)]))
    return dict(status="PASS", composites=8, rows=rows)


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
