#!/usr/bin/env python3
"""
Independent audit and projective-matroid reformulation of the arithmetic
Kakeya forcing-certificate game.

This file intentionally does not import the in-repo AK engine.  It checks the
literal and equal-suffix readings directly from the six-line certificate
data, after the invertible coordinate change

    (a,b) |-> (s,d) = (a+b,a-b).

Because every legal nonzero label has s != 0, row scaling reduces a label to
the single rational projective parameter rho=d/s.  If B is the signed
edge/seed incidence matrix and D_rho is the diagonal matrix of generator
parameters, the generator matrix is

    A = [ B | D_rho B ],

with paired columns (s_v,d_v) for every vertex v.

After the columns of the current wildcard set T are deleted, v fires exactly
when d_v is a coloop of the represented column matroid, equivalently when the
row space contains a vector supported only on d_v.  Thus forcing is a
distinguished-coloop peeling process.

The script:
  * verifies the score-14/9 literal-rule certificate independently;
  * prints exact generator combinations witnessing every firing round;
  * rejects the same certificate under the equal-suffix reading;
  * checks the complete score->1 family for several sizes;
  * performs deterministic hostile random comparisons against a direct
    untransformed support test implemented below.

Status: FINITE-EXACT.  No arithmetic-Kakeya theorem is claimed; the literal
rule is incompatible with the intuitive glued-copy definition.
"""

from fractions import Fraction
from itertools import product
import random


ZERO = (0, 0)
X = [ZERO, (1, 0), (0, 1), (1, 1)]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rref_with_provenance(rows, ncols):
    """Fraction RREF, retaining coefficients in the original input rows."""
    m = len(rows)
    aug = [
        [Fraction(x) for x in row] +
        [Fraction(int(i == j)) for j in range(m)]
        for i, row in enumerate(rows)
    ]
    rr = 0
    pivots = []
    for col in range(ncols):
        pivot = next((i for i in range(rr, m) if aug[i][col] != 0), None)
        if pivot is None:
            continue
        aug[rr], aug[pivot] = aug[pivot], aug[rr]
        scale = aug[rr][col]
        aug[rr] = [x / scale for x in aug[rr]]
        for i in range(m):
            if i != rr and aug[i][col] != 0:
                scale = aug[i][col]
                aug[i] = [x - scale * y for x, y in zip(aug[i], aug[rr])]
        pivots.append((col, rr))
        rr += 1
        if rr == m:
            break
    return [(col, aug[i][:ncols], aug[i][ncols:]) for col, i in pivots]


def rref(rows, ncols):
    """Fraction RREF without the much larger provenance augmentation."""
    work = [[Fraction(x) for x in row] for row in rows]
    rr = 0
    pivot_positions = []
    for col in range(ncols):
        pivot = next(
            (i for i in range(rr, len(work)) if work[i][col] != 0), None
        )
        if pivot is None:
            continue
        work[rr], work[pivot] = work[pivot], work[rr]
        scale = work[rr][col]
        work[rr] = [x / scale for x in work[rr]]
        for i in range(len(work)):
            if i != rr and work[i][col] != 0:
                scale = work[i][col]
                work[i] = [
                    x - scale * y for x, y in zip(work[i], work[rr])
                ]
        pivot_positions.append((col, rr))
        rr += 1
        if rr == len(work):
            break
    return [(col, work[i]) for col, i in pivot_positions]


class Instance:
    def __init__(self, dims, fs, seeds, t0=()):
        self.dims = tuple(dims)
        self.fs = tuple(dict(f) for f in fs)
        self.seeds = tuple((tuple(v), tuple(x)) for v, x in seeds)
        self.t0 = tuple(tuple(v) for v in t0)
        self.vertices = tuple(product(*[range(1, d + 1) for d in dims]))
        self.index = {v: i for i, v in enumerate(self.vertices)}

    @property
    def n(self):
        return len(self.vertices)

    @property
    def m(self):
        total = 0
        for i, f in enumerate(self.fs):
            suffix_count = 1
            for d in self.dims[i + 1:]:
                suffix_count *= d
            total += suffix_count * sum(x != ZERO for x in f.values())
        return total

    @property
    def r(self):
        return len(self.seeds)

    @property
    def t(self):
        return len(set(self.t0))

    @property
    def score(self):
        return Fraction(self.m + self.r, self.n - self.t)

    def labelled_generators(self, mode):
        """Return full-space generators as (name, sparse {vertex:(a,b)})."""
        ans = []
        for j, (v, x) in enumerate(self.seeds):
            ans.append((f"seed{j}:{v}:{x}", {self.index[v]: x}))
        for i, f in enumerate(self.fs):
            suffixes = tuple(product(
                *[range(1, d + 1) for d in self.dims[i + 1:]]
            ))
            for pre, x in sorted(f.items()):
                if x == ZERO:
                    continue
                right_pre = pre[:i] + (pre[i] + 1,)
                if mode == "strict":
                    pairs = ((s, s) for s in suffixes)
                elif mode == "literal":
                    pairs = product(suffixes, suffixes)
                else:
                    raise ValueError(mode)
                for left_suffix, right_suffix in pairs:
                    u = self.index[pre + left_suffix]
                    v = self.index[right_pre + right_suffix]
                    ans.append((
                        f"axis{i+1}:{pre}:{left_suffix}->{right_suffix}:{x}",
                        {u: x, v: (-x[0], -x[1])},
                    ))
        return ans


def transformed_rows(instance, mode, live, provenance=False):
    """Rows in interleaved live columns (s_v,d_v), plus generator names."""
    live_order = sorted(live)
    col = {
        (v, kind): 2 * j + kind
        for j, v in enumerate(live_order)
        for kind in (0, 1)
    }
    names = []
    rows = []
    for name, sparse in instance.labelled_generators(mode):
        row = [Fraction(0)] * (2 * len(live_order))
        for v, (a, b) in sparse.items():
            if v not in live:
                continue
            row[col[(v, 0)]] = a + b
            row[col[(v, 1)]] = a - b
        if any(row):
            names.append(name)
            rows.append(row)
    return live_order, col, names, rows


def projective_force(instance, mode, want_witnesses=False):
    """Distinguished-coloop peeling in (s,d) coordinates."""
    T = {instance.index[v] for v in instance.t0}
    trace = []
    witnesses = []
    while len(T) < instance.n:
        live = set(range(instance.n)) - T
        order, col, names, rows = transformed_rows(instance, mode, live)
        if not rows:
            break
        reduced = (
            rref_with_provenance(rows, len(rows[0]))
            if want_witnesses
            else [(p, row, None) for p, row in rref(rows, len(rows[0]))]
        )
        fired = []
        round_witnesses = {}
        for pivot, row, coeffs in reduced:
            if sum(x != 0 for x in row) != 1:
                continue
            j, kind = divmod(pivot, 2)
            if kind != 1:
                continue
            vertex = order[j]
            fired.append(vertex)
            if want_witnesses:
                round_witnesses[vertex] = [
                    (names[i], c) for i, c in enumerate(coeffs) if c
                ]
        if not fired:
            break
        trace.append(tuple(instance.vertices[v] for v in sorted(fired)))
        witnesses.append(round_witnesses)
        T.update(fired)
    return len(T) == instance.n, T, tuple(trace), tuple(witnesses)


def direct_force(instance, mode):
    """Independent untransformed rank test for hostile cross-checking.

    For a candidate v, let C contain every live coordinate column away from
    v together with the sum column x_v+y_v.  A legal target exists exactly
    when the difference column x_v-y_v is not in span(C).
    """
    T = {instance.index[v] for v in instance.t0}
    generators = instance.labelled_generators(mode)
    while len(T) < instance.n:
        live = sorted(set(range(instance.n)) - T)
        full_rows = []
        for _, sparse in generators:
            row = []
            for v in live:
                x = sparse.get(v, ZERO)
                row.extend((Fraction(x[0]), Fraction(x[1])))
            full_rows.append(row)
        if not full_rows:
            break
        fired = []
        for j, v in enumerate(live):
            rows_without_d = []
            rows_with_d = []
            for row in full_rows:
                base = []
                for k in range(len(live)):
                    if k == j:
                        continue
                    base.extend((row[2 * k], row[2 * k + 1]))
                sum_col = row[2 * j] + row[2 * j + 1]
                diff_col = row[2 * j] - row[2 * j + 1]
                rows_without_d.append(base + [sum_col])
                rows_with_d.append(base + [sum_col, diff_col])
            rank_without = len(rref(
                rows_without_d, len(rows_without_d[0])
            ))
            rank_with = len(rref(rows_with_d, len(rows_with_d[0])))
            if rank_with > rank_without:
                fired.append(v)
        if not fired:
            break
        T.update(fired)
    return len(T) == instance.n, T


def exploit_family(D, q):
    label = lambda j: (1, 0) if j % 2 else (0, 1)
    f1 = {(j,): label(j) for j in range(1, D)}
    f2 = {
        (side, e): (1, 1)
        for side in (1, D)
        for e in range(1, q)
    }
    seeds = [((1, 1), (1, 0)), ((1, 1), (0, 1))]
    for j in range(2, D + 1):
        previous = label(j - 1)
        complement = (0, 1) if previous == (1, 0) else (1, 0)
        seeds.append(((j, 1), complement))
    return Instance((D, q), (f1, f2), seeds)


def strict_13_over_7():
    """The independently reported equal-suffix [2,2,2] positive control."""
    fs = (
        {(1,): (1, 2)},
        {(1, 1): (1, 0), (2, 1): (1, 0)},
        {(2, 1, 1): (2, 1), (2, 2, 1): (0, 1)},
    )
    seeds = (
        ((1, 2, 2), (0, 1)),
        ((2, 2, 1), (1, 1)),
        ((1, 1, 1), (1, 0)),
    )
    return Instance((2, 2, 2), fs, seeds, (((1, 1, 2)),))


def random_instance(rng, dims):
    pool = [(1, 0), (0, 1), (1, 1), (1, 2), (2, 1)]
    fs = []
    for i, d in enumerate(dims):
        domain = product(
            *([range(1, x + 1) for x in dims[:i]] + [range(1, d)])
        )
        f = {}
        for key in domain:
            if rng.random() < 0.65:
                f[key] = rng.choice(pool)
        fs.append(f)
    vertices = tuple(product(*[range(1, d + 1) for d in dims]))
    seeds = []
    for _ in range(rng.randrange(0, min(5, len(vertices)) + 1)):
        seeds.append((rng.choice(vertices), rng.choice(pool)))
    t0 = rng.sample(vertices, rng.randrange(0, min(3, len(vertices)) + 1))
    return Instance(dims, fs, seeds, t0)


def main():
    print("PROJECTIVE NORMAL FORM")
    print("label (a,b) -> rho=(a-b)/(a+b); generator matrix A=[B|D_rho B]")
    print("firing = delete T columns, then peel v iff d_v is a coloop")

    inst = exploit_family(3, 3)
    require(inst.score == Fraction(14, 9), "wrong exploit score")
    for mode in ("literal", "strict"):
        p = projective_force(inst, mode, want_witnesses=True)
        d = direct_force(inst, mode)
        require((p[0], p[1]) == d, f"direct/projective mismatch in {mode}")
        print(f"3x3 mode={mode}: forced={p[0]} final={len(p[1])}/9 "
              f"score={inst.score} trace={p[2]}")
        if mode == "literal":
            for round_no, (vertices, ws) in enumerate(zip(p[2], p[3]), 1):
                print(f"  round {round_no}: {vertices}")
                for v in sorted(ws):
                    terms = " + ".join(f"({c})*{name}" for name, c in ws[v])
                    print(f"    witness d_{inst.vertices[v]}: {terms}")
    require(projective_force(inst, "literal")[0],
            "literal certificate should force")
    require(not projective_force(inst, "strict")[0],
            "strict certificate should fail")

    positive = strict_13_over_7()
    require(positive.score == Fraction(13, 7),
            "wrong strict positive-control score")
    positive_p = projective_force(positive, "strict")
    positive_d = direct_force(positive, "strict")
    require(positive_p[0] and (positive_p[0], positive_p[1]) == positive_d,
            "strict 13/7 positive control failed")
    print(f"STRICT POSITIVE CONTROL: score={positive.score} "
          f"trace={positive_p[2]}")

    print("EXPLOIT FAMILY")
    for D, q in ((2, 4), (3, 3), (4, 4), (6, 6)):
        g = exploit_family(D, q)
        expected = 1 + Fraction(1, D) + Fraction(1, q) - Fraction(1, D * q)
        require(g.score == expected, f"family score formula failed at {D,q}")
        literal = projective_force(g, "literal")[0]
        strict = projective_force(g, "strict")[0]
        require(literal and not strict, f"family forcing boundary failed at {D,q}")
        print(f"  {D}x{q}: score={g.score} literal={literal} strict={strict}")

    rng = random.Random(20260728)
    checked = 0
    for dims in ((2, 2), (2, 3), (2, 2, 2), (3, 2, 2)):
        for _ in range(20):
            g = random_instance(rng, dims)
            for mode in ("literal", "strict"):
                p = projective_force(g, mode)
                d = direct_force(g, mode)
                require((p[0], p[1]) == d,
                        f"random direct/projective mismatch {dims} {mode}")
                checked += 1
    print(f"HOSTILE CROSSCHECKS={checked}: all exact agreements")
    print("AUDIT_OK")


if __name__ == "__main__":
    main()
