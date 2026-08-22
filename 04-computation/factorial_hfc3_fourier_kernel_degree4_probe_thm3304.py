"""Exact Fourier--Dirichlet kernel and a finite-modular degree-four probe."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
import factorial_hfc3_symmetry_cells_support_thm3304 as probe


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


u, v, omega = sp.symbols("u v omega")
denominator = sp.expand((1 - u - v)
                        * (1 - omega*u - omega**2*v)
                        * (1 - omega**2*u - omega*v))
reduced = sp.rem(sp.Poly(denominator, omega),
                 sp.Poly(omega**2 + omega + 1, omega)).as_expr().expand()
require(reduced == 1 - u**3 - v**3 - 3*u*v,
        "cyclotomic denominator")


def kernel_coefficient(r, s):
    value = 0
    for ell in range(min(r, s) + 1):
        if (r - ell) % 3 or (s - ell) % 3:
            continue
        i = (r - ell) // 3
        j = (s - ell) // 3
        value += probe.multinomial((i, j, ell)) * 3**ell
    return value


# Selection and strict positivity on a complete finite hostile box.  The
# general statement follows directly from the displayed trinomial sum.
for r in range(25):
    for s in range(25):
        coefficient = kernel_coefficient(r, s)
        require((coefficient > 0) == ((r - s) % 3 == 0),
                f"selection/positivity at {(r, s)}")
require(kernel_coefficient(1, 0) == 0, "hostile nonzero character")
require(kernel_coefficient(1, 1) == 3, "positive balanced control")

# Degree-four omega-eigenbasis on e1=1:
#   b, a^2, a*b^2, b^4, a^3*b.
# Complete modular forms are evidence only; no projective elimination is
# claimed.  p>4*15+2 keeps every Dirichlet denominator invertible.
prime = 103
omega_p, variables, forms, _ = probe.cyclic_degree_four_forms_fast(
    prime=prime, exponents=(3, 6, 9, 12, 15))
require(omega_p == 46, "chosen cube root mod 103")
term_counts = {m: len(poly.terms()) for m, poly in forms.items()}
require(term_counts == {3: 35, 6: 210, 9: 715, 12: 1820, 15: 3874},
        "complete form term counts")

# Independent comparison at M3: rehomogenize the five Fourier monomials to
# degree four and integrate by the factorial numerator.  The forms may differ
# only by the common nonzero simplex denominator.
_, homogeneous_basis = probe.cyclic_degree_four_basis(prime)
slow3 = probe.form_to_sympy(probe.moment_form(homogeneous_basis, 3, prime),
                            variables, prime)
fast3 = forms[3]
slow_terms = dict(slow3.terms())
fast_terms = dict(fast3.terms())
key = next(key for key in slow_terms if fast_terms.get(key))
scale = int(slow_terms[key]) * pow(int(fast_terms[key]), prime - 2, prime)
scale %= prime
require(scale == 36, "common degree-12 simplex scale")
require(all((int(slow_terms.get(key, 0))
             - scale * int(fast_terms.get(key, 0))) % prime == 0
            for key in set(slow_terms) | set(fast_terms)),
        "fast/rehomogenized M3 comparison")

serialized = []
for exponent in sorted(forms):
    for powers, coefficient in forms[exponent].terms():
        serialized.append(str(exponent) + ":" + ",".join(map(str, powers))
                          + "=" + str(int(coefficient) % prime))
forms_hash = hashlib.sha256("\n".join(serialized).encode("ascii")).hexdigest()

payload = "\n".join([
    "kernel_EGF=1/(1-u^3-v^3-3*u*v)",
    ("simplex_moment=<a^r*b^s>=2*r!*s!/(r+s+2)!*"
     "[u^r*v^s]kernel_EGF"),
    "selection=zero iff r-s is nonzero mod 3; otherwise strictly positive",
    ("degree4_basis_on_e1=1=(b,a^2,a*b^2,b^4,a^3*b)"),
    ("phase_obstruction=the cubic moment has strictly positive coefficients; "
     "a null vector needs coefficient phase diameter at least pi/3"),
    "modular_prime=103; omega=46",
    "degree4_term_counts=(M3:35,M6:210,M9:715,M12:1820,M15:3874)",
    "degree4_forms_sha256=" + forms_hash,
    "degree4_projective_elimination=OPEN; no elimination certificate claimed",
    "controls=25x25 selection box; hostile (1,0); balanced (1,1); independent M3 rehomogenization scale 36",
]) + "\n"
print(payload, end="")
print("payload_sha256=" + hashlib.sha256(payload.encode("ascii")).hexdigest())
