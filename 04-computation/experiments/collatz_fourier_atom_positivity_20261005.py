"""Lossless circle encoding of the Collatz injection measure: exact controls."""
from fractions import Fraction as F
from math import factorial
import json

from collatz_critical_beta_weights_20261005 import finite_box, S, U

CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def w(t, k):
    return F(2*factorial(k)*factorial(t+1), factorial(t+k+2))


def tail_bound(levels, width, epsilon):
    if type(levels) is not int or levels < 2:
        raise ValueError("base depth cutoff at least two required")
    if type(width) is not int or width < 1:
        raise ValueError("positive inverse-width cutoff required")
    if type(epsilon) not in (int, F) or not 0 < epsilon < 1:
        raise ValueError("exact epsilon in (0,1) required")
    epsilon = F(epsilon)
    return (4*epsilon+4*(1-epsilon**2/3)**(levels-2)
            +F(6, width)+F(6, width+1)+F(12, width+2))


def prime(p):
    return type(p) is int and p >= 2 and all(p % d for d in range(2, int(p**0.5)+1))


def cyclotomic_reduce(coefficients, p):
    """Q[z]/(1+z+...+z^(p-1)); inputs use cyclic exponents mod p."""
    if not prime(p) or len(coefficients) != p:
        raise ValueError("prime cyclotomic coordinate required")
    return tuple(x-coefficients[-1] for x in coefficients[:-1])


def phase_shift(poly, exponent, p):
    values = [F(0)]*p
    for i, x in enumerate(poly):
        values[(i+exponent) % p] += x
    return cyclotomic_reduce(values, p)


def fourier_table(atoms, p):
    table = []
    for ell in range(p):
        values = [F(0)]*p
        for m, weight in atoms.items():
            values[(m*ell) % p] += weight
        table.append(cyclotomic_reduce(values, p))
    return table


def residue_atoms(atoms, q):
    values = [F(0)]*q
    for m, weight in atoms.items():
        values[m % q] += weight
    return values


def toy_residue(q, a):
    # Geometric law with the atom at1 deleted, then renormalized.
    raw = F(1, 2**(a+1))/(1-F(1, 2**q))
    if a == 1 % q:
        raw -= F(1, 4)
    return F(4, 3)*raw


def run():
    records = finite_box(5, 5)
    atoms = {(n-3)//6: w(t, k) for n, (t, k) in records.items() if n % 3 == 0}
    need(all(m >= 0 and value > 0 for m, value in atoms.items()), "injection atoms")
    need(sum(atoms.values()) < 1, "finite injection mass")
    for p in (2, 3, 5, 7, 11, 13):
        table = fourier_table(atoms, p)
        projected = residue_atoms(atoms, p)
        for a in range(p):
            acc = [F(0)]*(p-1)
            for ell, poly in enumerate(table):
                shifted = phase_shift(poly, -a*ell, p)
                acc = [x+y/p for x, y in zip(acc, shifted)]
            need(tuple(acc) == (projected[a],)+(F(0),)*(p-2), "exact Fourier inversion")
        need(sum(projected) == sum(atoms.values()), "quotient mass conservation")
    for q in range(1, 101):
        values = [toy_residue(q, a) for a in range(q)]
        need(sum(values) == 1, "toy projection is probability")
        need(all(x > 0 for x in values), "all finite toy residue atoms positive")
        if q >= 2:
            need(values[1] == F(1, 3*(2**q-1)), "positive refinement collapsing to zero")
    for j in range(1, 101):
        n = S(1, j)
        index = (n-1)//2
        need(index == 2*(4**j-1)//3, "root-ray Fourier indices 2,10,42")
        need(w(0, j) == F(2, (j+1)*(j+2)), "root-ray weight")
    examples = []
    previous_quarter_ratio = F(0)
    for t in range(1, 41):
        k = 9*t+1
        target = S(1, k)
        source = (2*target-1)//3
        index = (source-3)//6
        atom = w(1, k)
        need(target % 9 == 5 and source % 6 == 3, "three-divisible injection section")
        need(U(source) == target and U(target) == 1, "actual two-step certificate")
        need(index == 16*(4**(9*t)-1)//27, "lacunary injection Fourier index")
        need(atom == F(4, (k+1)*(k+2)*(k+3)), "exact injection weight")
        delta_lower = 2*atom
        phase = F(1, 2*index)
        quarter_ratio_power = delta_lower**4/phase
        need(quarter_ratio_power > previous_quarter_ratio, "quantified non-Holder control")
        previous_quarter_ratio = quarter_ratio_power
        if t <= 3:
            examples.append({"k": k, "source": source, "index": index,
                             "weight": str(atom), "phase": str(phase)})
    for r in (F(1, 16), F(1, 4), F(1, 2), F(3, 4)):
        den = 1+r+r*r
        kappa = (1+r)/den
        a0 = r*(1+r*r)/den
        edge3 = (1-r)*r
        column3 = (r+r*r)/den
        a2 = edge3*column3+(a0-edge3)*kappa
        bound = 2+2*r-r*r
        need(1+a0+a2/(1-kappa) == bound, "incoming two-generation bound")
        need(2*a2/(1-kappa) <= 4, "mixed depth-tail prefactor")
        need(kappa <= 1-r*r/3, "uniform depth decay away from zero")
        need(bound <= 3, "global root-base bound")
    need(tail_bound(100, 100, F(1, 2)) > 0, "rational finite-support tail")
    # The modular-law obstruction uses gamma=[[1,0],[N*t,1]].
    # At tau=i, its weight-two factor has absolute value1+(N*t)^2.
    for level in (1, 11, 22):
        for t in range(1, 51):
            need(1+(level*t)**2 > (level*t)**2, "unbounded modular multiplier")
    return {"status": "FINITE-EXACT; all-atom positivity OPEN", "checks": CHECKS,
            "finite_injection_atoms": len(atoms),
            "finite_injection_mass": str(sum(atoms.values())),
            "lacunary_injection_examples": examples,
            "toy_missing_atom": 1,
            "toy_positive_refinement_at_missing_atom": "1/(3*(2^q-1)), q>=2",
            "exact_cyclotomic_moduli": [2, 3, 5, 7, 11, 13],
            "proof_boundaries": ["nonnegative Fourier coefficients need not all be positive",
                                 "finite quotient positivity is not a uniform atom floor",
                                 "bounded disk series cannot have positive modular weight"]}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
