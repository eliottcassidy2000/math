"""Independent, diagnostic controls for the recurrence-bound proof.

This is not the main experiment. Exact checks use explicit exceptions rather
than assert, so optimized Python runs retain every check. Numerical checks
are diagnostics; the audit note gives the proof-level justification.
"""
from fractions import Fraction as F
from math import floor
import json
import mpmath as mp
import sympy as sp

checks = 0


def require(condition, label):
    global checks
    if not condition:
        raise RuntimeError(label)
    checks += 1


# Symbolic differential and both ordinary inverse-chart exclusions.
u, z, k = sp.symbols("u z k", nonzero=True)
quartic = k**2 * (u**2 + 1)**2 - 4*u*(u**2 - 1)
x = (k*z + k**2 + 2*u + k**2*u**2)/(2*u**2)
y = (x-k**2)*(x*u-1)/k
dxdu = sp.diff(x, u) + sp.diff(x, z)*sp.diff(quartic, u)/(2*z)
num = sp.fraction(sp.together(z*dxdu + 2*y))[0]
require(sp.rem(num, z**2-quartic, z) == 0, "dX/(2Y)=-du/z")
for target, expected in [
    (0, 4*u**2*(2*k**2*u+1)/k**2),
    (k**2, -4*u**2*(k**4-1)/k**2),
]:
    replacement = sp.solve(x-target, z)[0]
    require(sp.factor(replacement**2-quartic-expected) == 0,
            "inverse chart exceptional locus")
t = sp.symbols("t", positive=True)
spread = 4*t*(t**2-1)/(t**2+1)**2
require(sp.simplify(sp.diff(spread, t) +
                   4*(t**4-6*t**2+1)/(t**2+1)**3) == 0,
        "spread derivative")
require(sp.factor(1-spread**2 -
                   (t**4-6*t**2+1)**2/(t**2+1)**4) == 0,
        "absolute angular derivative identity")

# Exact circle controls, including floor equality and circular-wrap gaps.
# A rational rotation of prime denominator suffices to test this finite
# combinatorial argument whenever the prefix has distinct points.
circle_cases = 0
wrap_cases = 0
for prime in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 43]:
    for numerator in range(1, prime):
        alpha = F(numerator, prime)
        for nmax in range(2, min(prime-1, 20)+1):
            points = sorted((i*alpha % 1, i) for i in range(nmax+1))
            require(len({p for p, i in points}) == nmax+1,
                    "distinct circle prefix")
            gaps = []
            for j, (position, index) in enumerate(points):
                next_position, next_index = points[(j+1) % len(points)]
                gap = (next_position-position) % 1
                gaps.append((gap, next_index-index, j == len(points)-1))
            gap, difference, wrap = min(gaps)
            wrap_cases += int(wrap)
            for ell in [F(1, nmax), F(2, 2*nmax+1)]:
                require(floor(1/ell) == nmax, "floor convention")
                require(0 < gap < ell, "strict short gap")
                require(0 < abs(difference) <= nmax, "index size")
                require(difference*alpha % 1 == gap, "wrap-safe difference")
                for direction in [-1, 1]:
                    signed_index = direction*difference
                    require(direction*signed_index*alpha % 1 == gap,
                            "chosen orientation")
                circle_cases += 1
require(wrap_cases > 0, "wrap controls exercised")

# Numerical controls use direct smooth cubic-period integrals rather than
# relying solely on a library elliptic-integral call.
mp.mp.dps = 70
numerical_cases = []
for epsilon in [-1, 1]:
    for text_f in ["0.49", "0.25", "0.01", "1e-8"]:
        f = mp.mpf(text_f)
        d = 2 + (f if epsilon == 1 else 0)
        period = mp.quad(lambda v: 2/mp.sqrt((v*v+f)*(v*v+d)),
                         [0, mp.sqrt(f), 1, mp.inf])
        library_period = 2*mp.ellipk(1-f/d)/mp.sqrt(d)
        require(abs(period-library_period) < mp.mpf("1e-60"),
                "full cubic period")
        # Direct oval integration x=-1+(e2+1) sin(theta)^2.
        oval_period = mp.quad(
            lambda theta: 2/mp.sqrt(d-(d-f)*mp.sin(theta)**2),
            [0, mp.pi/2])
        require(abs(period-oval_period) < mp.mpf("1e-60"),
                "equal component periods")
        density = lambda s: 1/(4*mp.sqrt(1-s*s)*
                                mp.sqrt(1+epsilon*(f-s)))
        lower = mp.quad(density, [0, f])
        upper = mp.quad(density, [f, mp.mpf("0.5")])
        require(max(lower, upper) >= mp.pi/(24*mp.sqrt(6)),
                "global angular arc lower bound")
        require(lower+upper < period/2, "arcs below half period")
        elementary = (4*mp.asinh(1)+mp.log(d/f))/mp.sqrt(d)
        require(period < elementary, "three-piece upper bound")
        ff = F(text_f)
        dd = F(2)+(ff if epsilon == 1 else 0)
        j = 0
        while 2**j*ff < dd:
            j += 1
        w = F(5, 7)*(4+j)
        require(period < mp.mpf(w.numerator)/w.denominator,
                "rational period certificate")
        require(lower > f/5 and upper > (mp.mpf("0.5")-f)/5,
                "rational arc certificate")
        numerical_cases.append({"epsilon": epsilon, "f": text_f,
                                "J": j,
                                "global_even_budget": 2*floor(F(100,7)*(4+j))})

constants = []
for epsilon in [-1, 1]:
    c = mp.quad(lambda s: 1/(4*mp.sqrt(1-s*s)*mp.sqrt(1-epsilon*s)),
                [0, mp.mpf("0.5")])
    if epsilon == 1:
        formula = (mp.atanh(mp.sqrt(3)/2)-mp.atanh(1/mp.sqrt(2)))/(2*mp.sqrt(2))
    else:
        formula = (mp.atanh(1/mp.sqrt(2))-mp.atanh(mp.mpf("0.5")))/(2*mp.sqrt(2))
    require(abs(c-formula) < mp.mpf("1e-65"), "cusp arc constant")
    f = mp.mpf("1e-20")
    d = 2 + (f if epsilon == 1 else 0)
    period = 2*mp.ellipk(1-f/d)/mp.sqrt(d)
    residual = period-mp.log(32/f)/mp.sqrt(2)
    require(abs(residual)/(f*mp.log(1/f)) < 1, "cusp period diagnostic")
    constants.append({"epsilon": epsilon, "C": mp.nstr(c, 35),
                      "leading_coefficient": mp.nstr(mp.sqrt(2)/c, 35)})

print(json.dumps({"explicit_checks": checks, "circle_cases": circle_cases,
                  "wrap_cases": wrap_cases, "numerical_cases": numerical_cases,
                  "cusp_constants": constants}, indent=2, sort_keys=True))
