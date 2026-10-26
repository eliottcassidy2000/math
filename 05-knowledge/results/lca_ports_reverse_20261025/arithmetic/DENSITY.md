# Exact inverse-depth densities for the reset parameters

2026-10-25. **PROVED:** the natural-density statements below, with explicit
finite-cylinder tail control. **FINITE-EXACT:** the modular checks in the
standalone block. **OPEN:** whether any reset belongs to the marked 27
ray, and whether the guarded orbit from `(6,27)` halts. Density zero does
not mean absence, and does not answer either question.

This is a short continuation of [ANCESTRY.md](ANCESTRY.md). Use its odd
integer map `U`, terminating predecessor algorithm, and reset parameter

    X(t)=(128*64^t+7)/9,    t>=0,    R=6t+5.

Write `L(t)` for the number of legal predecessor steps from `X(t)` to
its integer root. For every integer `j>=0`, the exact conclusions are

    dens{t>=0 : L(t)>=j} = 2^(-j),
    dens{t>=0 : L(t)=j}  = 2^(-j-1).                (D1)

Here `dens` is ordinary natural density in the parameter `t`. Equivalently
these are relative densities among `R=5 (mod 6)`. For every fixed odd
integer `b>1`,

    dens{t>=0 : root(X(t))=b} = 0.                  (D2)

In particular, the possible reset parameters belonging to the unguarded
27 ray have density zero. This concerns the reset parameters, not the
frequency of resets along that ray or the natural-number density of
their exponentially spaced scalar values.

## One inverse branch retains a full ternary residue coordinate

For a 3-adic target `y`, the predecessor with `r=v_3(y)` is legal on
the cylinder

    C_r = (-1)^r*3^r + 3^(r+1) Z_3,    r>=0.

These cylinders are disjoint and have Haar masses `3^(-r-1)`. Writing
`y=3^r*((-1)^r+3z)`, the predecessor becomes

    P_r(y)=(2^(r+1)*(-1)^r+1)/3 + 2^(r+1)*z.       (D3)

The constant is an integer and the coefficient of `z` is a 3-adic unit.
Thus this affine map carries `C_r` bijectively onto `Z_3`: after a legal
inverse branch, the remaining ternary residue coordinate is unrestricted.

For any prescribed finite word `(r_0,...,r_(j-1))`, its exact domain is
one residue class modulo `3^A`, where `A=sum_i(r_i+1)`. This follows by
induction from `(D3)`: each next cylinder condition chooses one extension
of the current residue class. Different words of the same length have
disjoint domains, since the first differing valuation cannot hold twice.

The map `t -> X(t)` permutes residues modulo every power of 3, as proved
in ANCESTRY by `v_3(X(t)-X(s))=v_3(t-s)`. Therefore each of these finite
word domains pulls back to exactly one parameter class modulo `3^A`,
of natural density `3^-A`. This transports finite congruences; no general
identification of Haar measure with natural density is being assumed.

## Finite truncation controls the countably many branches

Fix `j>=1` and a cutoff `M>=1` for the inverse valuations. Put

    s_M = sum_(0<=r<M) 3^(-r-1) = (1-3^-M)/2.

Let `D_(j,M)` mean that the first `j` inverse steps are legal and all
their valuations satisfy `r<M`. It is a finite union of disjoint
parameter classes, with exact density

    dens D_(j,M) = s_M^j.                          (D4)

Any parameter with at least `j` legal inverse steps but outside this
truncated event has a first large valuation at some position `i<j`.
Cover it by the event `E_(i,M)` that the preceding `i` steps are legal
with valuations below `M`, and the current target is divisible by `3^M`.
Whether the large-valuation inverse itself is legal is not required for
this cover. Each `E_(i,M)` is again a finite union of residue classes;
by `(D3)`,

    dens E_(i,M) = s_M^i * 3^-M.

The covers for different `i` are disjoint, because they specify the
first large valuation. Hence the lower and upper natural densities of
`D_j={t:L(t)>=j}` satisfy

    s_M^j <= lower_dens D_j <= upper_dens D_j
           <= s_M^j + 3^-M*sum_(0<=i<j)s_M^i.       (D5)

The error is less than `2*3^-M`. Sending `M` to infinity proves the
first part of `(D1)`; subtracting the nested events `D_(j+1)` from `D_j`
proves the second part, including `j=0`.

For the actual sources `X(t)`, every legal inverse is an odd positive
integer greater than 1 and strictly smaller than its target. Thus the
finite-word cylinder events coincide with the ordinary integer inverse
algorithm, and every `L(t)` is finite. The 3-adic point 1 has infinitely
many legal inverses, all equal to 1, but no `X(t)>1` can reach it:
`U(1)=1` and an integer predecessor of a target greater than 1 is
greater than 1. This exceptional cycle is not silently counted as a
finite-depth integer root.

## A fixed root receives a density-zero set of parameters

If `b` is not an integer root, the set in `(D2)` is empty. Otherwise,
the points with root `b` are exactly its unique forward ray. Among
points at inverse depth less than `j`, at most the `j` scalars
`b,U(b),...,U^(j-1)(b)` have that root. Since `X(t)` is strictly
increasing, these account for at most `j` individual parameters. Every
remaining parameter with root `b` belongs to `D_j`. Its upper natural
density is therefore at most `2^-j`, for every `j`, proving `(D2)`.

This is a distributional separation of integer components. It loses the
actual root label and chronological credit. An infinite exceptional set
of density zero could still contain a reset on the marked ray. Neither
the root membership predicate nor the guarded debt inequality has been
removed from the fixed-source problem.

## Independent finite modular controls

The block below enumerates every parameter residue modulo `3^(j*M)`
for six `(j,M)` pairs, classifying it as truncated success, first large
valuation, or failed inverse. It independently replays ordinary integer
inverse chains for the first 80 representatives in each universe. The
modular state loses exactly `r+1` ternary digits per successful step;
it never supplies an unobserved digit. The literal 1-cycle and the
positive-integer boundary are explicit controls. No producer code is
imported, and all gates remain active under optimized Python.
Normal and optimized executions both returned `PASS`, with identical
complete JSON SHA-256
`3940c5f7e4c2536d57ad30988f3e896d1b8ce571cadaf99bf14aea75022d9436`.

```python
import json
from fractions import Fraction


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def v3(n):
    require(n > 0, "positive valuation domain")
    r = 0
    while n % 3 == 0:
        n //= 3
        r += 1
    return r


def x_residue(t, modulus):
    return ((128*pow(64, t, 9*modulus)+7)//9) % modulus


def modular_classification(t, j, cutoff):
    precision = j*cutoff
    y = x_residue(t, 3**precision)
    for i in range(j):
        require(precision >= cutoff*(j-i), "remaining precision budget")
        if y % (3**cutoff) == 0:
            return "LARGE", i
        r = v3(y)
        require(r < cutoff, "observed small valuation")
        m = y//(3**r)
        if m % 3 != pow(-1, r, 3):
            return "FAIL", i
        numerator = (1 << (r+1))*m+1
        require(numerator % 3 == 0, "integer inverse residue")
        precision -= r+1
        y = (numerator//3) % (3**precision)
    return "PASS", j


def literal_inverse(y):
    require(y > 1 and y % 2 == 1, "integer target boundary")
    r = v3(y)
    m = y//(3**r)
    numerator = (1 << (r+1))*m+1
    if numerator % 3:
        return None
    x = numerator//3
    require(1 < x < y and x % 2 == 1, "positive decreasing inverse")
    # Independent replay by actual odd/even scalar macros.
    z, k = (3*x-1)//2, 1
    while z % 2 == 0:
        z, k = 3*z//2, k+1
    require((z, k) == (y, r+1), "literal forward return")
    return x, r


def integer_classification(t, j, cutoff):
    y = (128*pow(64, t)+7)//9
    for i in range(j):
        if v3(y) >= cutoff:
            return "LARGE", i
        inverse = literal_inverse(y)
        if inverse is None:
            return "FAIL", i
        y, _ = inverse
    return "PASS", j


cases = [(1, 1), (1, 5), (2, 3), (3, 2), (4, 2), (5, 1)]
receipts = []
integer_replays = 0
for j, cutoff in cases:
    modulus = 3**(j*cutoff)
    successes, failures = 0, 0
    first_large = [0]*j
    residues = set()
    for t in range(modulus):
        residues.add(x_residue(t, modulus))
        result, index = modular_classification(t, j, cutoff)
        if result == "PASS":
            successes += 1
        elif result == "FAIL":
            failures += 1
        else:
            first_large[index] += 1
        if t < 80:
            require((result, index) == integer_classification(t, j, cutoff),
                    "modular versus literal integer classification")
            integer_replays += 1
    require(len(residues) == modulus, "parameter residue permutation")
    s = (1-Fraction(1, 3**cutoff))/2
    require(Fraction(successes, modulus) == s**j,
            "exact truncated inverse-depth density")
    for i, count in enumerate(first_large):
        require(Fraction(count, modulus) == s**i/Fraction(3**cutoff),
                "exact first-large tail cover density")
    require(successes+failures+sum(first_large) == modulus,
            "complete disjoint classification")
    lower = Fraction(successes, modulus)
    upper = Fraction(successes+sum(first_large), modulus)
    require(lower <= Fraction(1, 2**j) <= upper,
            "proved limiting density lies in finite bracket")
    require(upper-lower < Fraction(2, 3**cutoff), "uniform tail bound")
    receipts.append({
        "j": j, "cutoff": cutoff, "modulus": modulus,
        "successes": successes, "failures": failures,
        "first_large": first_large,
        "density_lower": str(lower), "density_upper": str(upper),
        "proved_limit": str(Fraction(1, 2**j)),
    })

require((2*1+1)//3 == 1, "exceptional inverse fixed point")
require(all((128*pow(64, t)+7)//9 > 1 for t in range(80)),
        "reset family excludes fixed point1")
require(literal_inverse(27) is None, "marked integer root")

print(json.dumps({
    "status": "PASS", "scope": "FINITE_MODULAR_CONTROLS_FOR_PROVED_DENSITIES",
    "parameter_residues_checked": sum(c["modulus"] for c in receipts),
    "integer_replays": integer_replays, "cases": receipts,
    "fixed_root_consequence": "Natural density zero; absence remains open",
}, sort_keys=True, separators=(",", ":")))
```
