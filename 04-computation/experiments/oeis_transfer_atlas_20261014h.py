"""Exact fixed-weight run compiler and scoped repository OEIS reference census.

Run: python -B -X utf8 04-computation/experiments/oeis_transfer_atlas_20261014h.py
No network calls, submissions, or changes to other research artifacts.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations
from math import factorial
from pathlib import Path
import json
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]
STEM = "oeis_transfer_atlas_20261014h"
SNAPSHOT = "de52220938d0"


def natural(x):
    if type(x) is not int or x < 0:
        raise ValueError("expected an exact nonnegative integer")
    return x


def rational(x):
    if type(x) not in (int, F):
        raise ValueError("weights must be exact int or Fraction, not bool/float")
    return F(x)


def add(a, b):
    return tuple((a[i] if i < len(a) else 0) +
                 (b[i] if i < len(b) else 0) for i in range(max(len(a), len(b))))


def mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return tuple(out)


def deriv(a):
    return tuple(i * a[i] for i in range(1, len(a)))


def ode(u, v):
    """L Y' + M Y + R = 0; exact fixed weight pair remains in the receipt."""
    u, v = rational(u), rational(v)
    a, b = u - 1, v - 1
    N, D = (F(0), F(1), F(0), -a*b), (F(1), -a-b, a*b)
    C = add(mul(deriv(N), D), tuple(-x for x in mul(N, deriv(D))))
    return mul(D, mul(N, N)), mul(C, add(N, tuple(-x for x in D))), mul(C, D)


def run_values(n, u, v):
    """P_j(u,v), 0<=j<=n, by a width-at-most-seven exact recurrence.

    O(n) rational arithmetic for FIXED u,v; bit cost and all coefficients of
    the bivariate polynomial are different complexity claims.
    """
    natural(n)
    L, M, R = ode(u, v)
    if L[0] != 0 or L[1] != 0 or M[0] != -1:
        raise ArithmeticError("invalid normalized kernel")
    y = []
    for k in range(n + 1):
        total = R[k] if k < len(R) else F(0)
        for i in range(1, min(k, len(M)-1)+1):
            total += M[i] * y[k-i]
        for i in range(2, min(k+1, len(L)-1)+1):
            total += L[i] * (k-i+1) * y[k-i+1]
        y.append(total)
    return tuple(y)


def w_values(n):
    natural(n)
    w = []
    forcing = {0: 1, 2: 3, 4: -5, 6: 1}
    for k in range(n+1):
        co = (k, -3, k+2, 5, 8-k, -1, 6-k)
        w.append(forcing.get(k, 0) + sum(co[i-1]*w[k-i]
                 for i in range(1, min(7, k)+1)))
    return tuple(w)


def brute_histogram(n):
    natural(n)
    ans = Counter()
    for p in permutations(range(n)):
        up = sum(y == x+1 for x, y in zip(p, p[1:]))
        down = sum(y == x-1 for x, y in zip(p, p[1:]))
        ans[up, down] += 1
    return ans


def composition_values(n, u, v):
    """Independent triangular polynomial-power calculation, intentionally slow."""
    a, b = rational(u)-1, rational(v)-1
    block = [F(0)] + [F(1)] + [a**(j-1)+b**(j-1) for j in range(2, n+1)]
    powers = [F(1)] + [F(0)]*n
    ans = powers[:]
    for k in range(1, n+1):
        powers = [sum(powers[j]*block[i-j] for j in range(i))
                  if i else F(0) for i in range(n+1)]
        ans = [x+factorial(k)*y for x, y in zip(ans, powers)]
    return tuple(ans)


def tournament_moments(n):
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    source_paths = tuple(permutations(range(n)))
    counts = []
    for mask in range(1 << len(edges)):
        out = [[False]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            out[i][j] = bool(mask & (1 << k))
            out[j][i] = not out[i][j]
        counts.append(sum(all(out[x][y] for x, y in zip(p, p[1:])) for p in source_paths))
    return len(counts), sum(counts), sum(x*x for x in counts), sorted(set(counts))


def reference_census():
    """Lexical references, never an attribution or accepted-contribution census.

    Explicit universe: tracked .md/.txt/.py/.out/.tex in SNAPSHOT under the
    five named research directories; no agents/navigation or JSON artifacts.
    The fixed Git object prevents concurrent writers changing this census.
    """
    roots = ("01-canon", "03-artifacts", "04-computation", "05-knowledge", "07-reflections")
    suffixes = {".md", ".txt", ".py", ".out", ".tex"}
    refs = {}
    candidates = []
    scanned = 0
    files = subprocess.check_output(["git", "ls-tree", "-r", "--name-only", SNAPSHOT,
                                     "--", *roots], cwd=ROOT, encoding="utf-8").splitlines()
    files = [p for p in files if Path(p).suffix in suffixes and
             Path(p).name not in ("INDEX.md", "MISTAKES.md") and STEM not in Path(p).name]
    scanned = len(files)
    specs = [f":(glob){root}/**/*{suffix}" for root in roots for suffix in sorted(suffixes)]
    specs += [":(exclude)**/INDEX.md", ":(exclude)**/MISTAKES.md"]
    matches = subprocess.run(["git", "grep", "-I", "-o", "-E",
                              r"\bA[0-9]{6}\b", SNAPSHOT, "--", *specs], cwd=ROOT,
                              capture_output=True, encoding="utf-8", check=False)
    if matches.returncode not in (0, 1):
        raise RuntimeError(matches.stderr)
    by_file = {}
    for line in matches.stdout.splitlines():
        _, rel, name = line.split(":", 2)
        rel = rel.replace("\\", "/")
        by_file.setdefault(rel, set()).add(name)
    for rel in sorted(by_file):
        for name in sorted(by_file[rel]):
            refs.setdefault(name, []).append(rel)
    for raw in sorted(files):
        rel = raw.replace("\\", "/")
        name = Path(rel).name.lower()
        if "oeis" in rel.lower() or "new_sequences" in name:
            candidates.append({"path": rel, "references": sorted(by_file.get(rel, ())),
                "status": "DISCOVERY_ROUTE_ONLY; title/content is not publication evidence"})
    return {"type": "LEXICAL_REFERENCE_CENSUS_NOT_ATTRIBUTION", "snapshot": SNAPSHOT,
            "roots": list(roots),
            "suffixes": sorted(suffixes), "excluded": ["INDEX.md", "MISTAKES.md", STEM+".*"],
            "scanned_files": scanned, "distinct_A_ids": len(refs),
            "reference_files": len(set(p for paths in refs.values() for p in paths)),
            "ids": refs, "candidate_filename_routes": candidates}


def verify():
    checks = 0
    def check(ok):
        nonlocal checks
        checks += 1
        if not ok:
            raise AssertionError(f"control {checks} failed")
    weights = (-1, 0, F(1, 2), 1, 2, 3)
    hist = [brute_histogram(n) for n in range(9)]
    for u in weights:
        for v in weights:
            vals = run_values(8, u, v)
            for n in range(9):
                expected = sum(count*F(u)**up*F(v)**down
                               for (up, down), count in hist[n].items())
                check(vals[n] == expected)
            check(vals == run_values(8, v, u))
    for u, v in ((2,0), (1,0), (0,0), (1,1), (F(2,3), F(-1,4))):
        check(run_values(30,u,v) == composition_values(30,u,v))
    known = (1,1,2,8,32,158,928,6350,49752,439670,4327904,46963358,
             556953448,7166360054,99428495088,1479600188798,
             23506712352248,397095175477430,7107209383674112,
             134345623603516190,2674426516381764744)
    w = w_values(1000)
    check(w[:len(known)] == known)
    check(tuple(map(F,w[:201])) == run_values(200,2,0))
    for n in range(101):
        check(run_values(n,1,1)[-1] == factorial(n))
    nud = run_values(100,1,0)
    for n in range(3,101):
        check(nud[n] == (n-1)*nud[n-1]+(n-2)*nud[n-2])
    check(w[3] != 2*w[2]+w[1])
    check(w[8] == 49752 and w[8] != 49760)
    moment_rows = []
    for n in range(1,6):
        size, total, squares, support = tournament_moments(n)
        check(F(total,size) == F(factorial(n),2**(n-1)))
        check(F(squares,size) == F(factorial(n)*w[n],4**(n-1)))
        moment_rows.append({"n":n,"labelled_tournaments":size,"sum_H":total,
                            "sum_H_squared":squares,"H_support":support})
    check(moment_rows[2]["H_support"] == [1,3])
    # Native parameter identities do not permit bool/float aliases.
    rejects = 0
    for n,u,v in ((True,2,0),(-1,2,0),(3.0,2,0),(3,True,0),(3,2,0.0)):
        try:
            run_values(n,u,v)
        except ValueError:
            rejects += 1
        else:
            raise AssertionError("malformed source accepted")
    check(rejects == 5)
    L,M,R = ode(2,0)
    check(L == (0,0,1,0,1,0,-1,0,-1))
    check(M == (-1,1,-3,5,5,3,-1,-1))
    check(R == (1,0,3,0,-5,0,1))
    return {"status":"FINITE-EXACT; mathematical proof in companion note",
            "checks":checks,"permutation_controls_through_n":8,
            "rational_weight_pairs":36,"W_n_0_through_20":list(w[:21]),
            "W1000_bits":w[1000].bit_length(),
            "W1000_big_endian_sha256":sha256(w[1000].to_bytes((w[1000].bit_length()+7)//8,"big")).hexdigest(),
            "moment_controls":moment_rows,
            "ODE_W":{"L":list(map(int,L)),"M":list(map(int,M)),"R":list(map(int,R))}}


if __name__ == "__main__":
    print(json.dumps({"verification":verify(),"repository_reference_census":reference_census()},
                     indent=2,sort_keys=True))
