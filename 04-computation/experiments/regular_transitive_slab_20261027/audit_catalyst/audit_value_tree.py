"""Verify a value-level merge certificate (derivation DAG) for the catalytic minimum g (THM-4820 audit).

Input lines 'D xa xb = ya yb + za zb' (x = y + z in coordinates a + b w); each non-unit value has one derivation.
Expanding the DAG gives a binary tree whose leaves are the six units.  The verifier checks y + z = x, that every
non-unit value reached has a derivation, that the expansion is well founded (no value depends on itself), that all
leaves are units, and recomputes the total cost sum beta+ with the norm computed in the sqrt(-3) model
(value = (X + Y sqrt(-3))/2, X = 2a - b, Y = b, norm = (X^2 + 3Y^2)/4).
Usage: python3 audit_value_tree.py FILE xa xb
"""
import re
import sys

sys.setrecursionlimit(100000)
fn, xa, xb = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
UNITS = {(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)}
pat = re.compile(r"^D\s+(-?\d+)\s+(-?\d+)\s*=\s*(-?\d+)\s+(-?\d+)\s*\+\s*(-?\d+)\s+(-?\d+)$")
der = {}
for line in open(fn):
    m = pat.match(line.strip())
    assert m, line
    v = tuple(map(int, m.groups()))
    x, y, z = v[0:2], v[2:4], v[4:6]
    assert (y[0] + z[0], y[1] + z[1]) == x, line
    assert x not in der, f"two derivations for {x}"
    der[x] = (y, z)


def norm(v):
    X, Y = 2 * v[0] - v[1], v[1]
    n4 = X * X + 3 * Y * Y
    assert n4 % 4 == 0
    return n4 // 4


memo = {}
state = {}


def cost(v):
    """returns (cost, number of leaves, total slack) of the expanded tree below v"""
    if v in memo:
        return memo[v]
    if v in UNITS:
        memo[v] = (0, 1, 0)
        return memo[v]
    assert v in der, f"value {v} has no derivation and is not a unit"
    assert state.get(v) != "open", f"cyclic derivation at {v}"
    state[v] = "open"
    y, z = der[v]
    cy, ly, sy = cost(y)
    cz, lz, sz = cost(z)
    beta = norm(v) - norm(y) - norm(z)
    bp = max(beta, 0)
    sl = bp - (norm(v) // 3 - norm(y) // 3 - norm(z) // 3)
    assert sl >= 0
    state[v] = "done"
    memo[v] = (cy + cz + bp, ly + lz, sy + sz + sl)
    return memo[v]


c, L, s = cost((xa, xb))
q = norm((xa, xb))
print(f"{fn}: root {(xa, xb)} norm {q}; {len(der)} distinct derived values; expanded tree has {L} unit leaves; "
      f"cost = {c}; Phi(root) = {q // 3}; total slack = {s}; cost - Phi = {c - q // 3}")
assert c == q // 3 + s
