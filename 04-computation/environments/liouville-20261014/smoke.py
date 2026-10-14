#!/usr/bin/env python3
"""Bounded, deterministic toolchain checks; not a mathematical proof audit."""
import importlib.metadata as metadata
import json
import platform
from fractions import Fraction
from time import perf_counter
import mpmath as mp
import numpy as np
import scipy.linalg
import sympy as sp

started = perf_counter()
x, y = sp.symbols("x y")
checks = {}
checks["rational_arithmetic"] = Fraction(1, 3) + Fraction(1, 6) == Fraction(1, 2)
checks["symbolic_factorization"] = sp.expand((x - 1)*(x + 1)*(x*x + 1)) == x**4 - 1
basis = sp.groebner([x*y - 1, y*y - 1], x, y, domain=sp.QQ)
checks["exact_groebner_reduction"] = basis.reduce(x*x - 1)[1] == 0
checks["hostile_nonzero_remainder"] = basis.reduce(x + 1)[1] != 0
matrix = sp.Matrix([[2, 1], [1, 3]])
checks["exact_linear_algebra"] = matrix.det() == 5 and matrix.inv()*matrix == sp.eye(2)
a = np.array([[2., 1.], [1., 3.]])
b = np.array([1., 2.])
solution = scipy.linalg.solve(a, b)
checks["numerical_linear_algebra"] = bool(np.linalg.norm(a @ solution - b) < 1e-12)
with mp.workdps(80):
    checks["multiprecision_arithmetic"] = bool(abs(mp.sqrt(2)**2 - 2) < mp.mpf("1e-75"))
result = {
  "scope": "Toolchain smoke checks only; no repository theorem validated.",
  "python": platform.python_version(),
  "packages": {p: metadata.version(p) for p in ["mpmath", "numpy", "scipy", "sympy"]},
  "checks": checks,
  "passed": all(checks.values()),
  "seconds": round(perf_counter() - started, 6),
}
print(json.dumps(result, indent=2, sort_keys=True))
raise SystemExit(0 if result["passed"] else 1)
