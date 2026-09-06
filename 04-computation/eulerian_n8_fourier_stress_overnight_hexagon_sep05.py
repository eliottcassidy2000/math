#!/usr/bin/env python3
"""Exact compiled n8 stress, with independent Python full-universe controls.

Temporary compilation is local and disposable. All semantic checks survive -O.
A nonzero modular determinant certifies a rational determinant is nonzero;
no inference is made from a zero modular determinant alone.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

import eulerian_boolean_fourier_repair_overnight_hexagon_sep05 as previous


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def modular_det(a, prime):
    require(all(len(row) == len(a) for row in a), "square modular matrix")
    b = [[x % prime for x in row] for row in a]
    det = 1
    for j in range(len(b)):
        pivot = next((i for i in range(j, len(b)) if b[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            b[pivot], b[j] = b[j], b[pivot]
            det = -det
        value = b[j][j]
        det = det * value % prime
        inverse = pow(value, -1, prime)
        for i in range(j + 1, len(b)):
            multiplier = b[i][j] * inverse % prime
            if multiplier:
                for k in range(j + 1, len(b)):
                    b[i][k] = (b[i][k] - multiplier * b[j][k]) % prime
            b[i][j] = 0
    return det % prime


def matrix_product(a, b):
    cols = list(zip(*b))
    return [[sum(x*y for x, y in zip(row, col)) for col in cols] for row in a]


def independent_control(data):
    n = data["n"]
    reps, orbit, sizes, index = previous.native.generator_orbits(n)
    orbits, complement, gauge = previous.switching_orbits(n)
    compressed = [sum(((r >> index[e]) & 1) << j for j, e in enumerate(gauge)) for r in reps]
    psi = [[sum(1-2*((h & x).bit_count() % 2) for h in members) for members in orbits] for x in compressed]
    triangles = previous.native.old.simple_cycles(n, 3, index)
    weighted = previous.native.old.weighted_operator(reps, orbit, triangles)
    require(data["reps"] == reps and data["sizes"] == sizes, "independent primal representatives and sizes")
    require(data["dual_reps"] == [min(o) for o in orbits] and data["dual_sizes"] == list(map(len, orbits)), "independent dual representatives and sizes")
    require(data["M"] == weighted and data["psi"] == psi and data["complement"] == complement, "independent full matrices and complement")


def main():
    source = Path(__file__).with_suffix(".cpp")
    reports = []
    primes = (101, 1009, 1000003)
    # Literal determinant positive/hostile controls, including nonzero integers
    # whose reduction at one prime vanishes.
    for p in primes:
        require(modular_det([[2, 3], [5, 7]], p) == p-1, "negative determinant control")
        require(modular_det([[1, 2], [2, 4]], p) == 0, "singular control")
    require(modular_det([[101]], 101) == 0 and modular_det([[101]], 1009) == 101, "modular-zero is not rational-zero")
    with tempfile.TemporaryDirectory(prefix="math-n8-fourier-") as directory:
        binary = str(Path(directory) / "enumerate")
        subprocess.run(["g++", "-O3", "-std=c++17", "-Wall", "-Wextra", str(source), "-o", binary], check=True)
        for n in range(3, 9):
            data = json.loads(subprocess.check_output([binary, str(n)], text=True))
            if n <= 7:
                independent_control(data)
            reps, psi, complement, weighted = (data[k] for k in ("reps", "psi", "complement", "M"))
            q = len(reps)
            require(q == previous.native.fixed.OrbitSampler(n).orbit_count, "independent Burnside count")
            require(sum(data["sizes"]) == sum(data["dual_sizes"]) == data["states"] == 2**((n-1)*(n-2)//2), "complete labelled universes")
            require(all(complement[complement[j]] == j for j in range(q)), "complement involution")
            even = [i for i, r in enumerate(reps) if r.bit_count() % 2 == 0]
            odd = [i for i, r in enumerate(reps) if r.bit_count() % 2]
            fixed = [j for j, k in enumerate(complement) if j == k]
            pairs = [(j, k) for j, k in enumerate(complement) if j < k]
            require(len(fixed) == len(even)-len(odd) and len(pairs) == len(odd), "parity/Fourier dimensions")
            require(all(psi[i][j] == 0 for i in odd for j in fixed), "selfdual modes even-supported")
            require(all(psi[i][a] == (-1)**reps[i].bit_count()*psi[i][b] for i in range(q) for a,b in pairs), "free-pair parity relation")
            # Independent character eigenvalue check uses the literal full matrix.
            product = matrix_product(weighted, psi)
            require(all(product[i][j] == data["eigenvalues"][j]*psi[i][j] for i in range(q) for j in range(q)), "complete Fourier diagonalization of multiplicity only")
            require(all(weighted[i][j] == 0 for side in (even, odd) for i in side for j in side), "native bipartition")
            c = [[int(weighted[i][j] > 0) for j in even] for i in odd]
            s = [[psi[i][j] for j in fixed] for i in even]
            pmat = [[psi[i][a]+psi[i][b] for a,b in pairs] for i in even]
            odd_basis = [[psi[i][a]-psi[i][b] for a,b in pairs] for i in odd]
            a = matrix_product(c, pmat)
            dets = {str(prime): modular_det(a, prime) for prime in primes}
            basis_dets = {str(prime): modular_det([ss+pp for ss,pp in zip(s,pmat)], prime) for prime in primes}
            odd_basis_dets = {str(prime): modular_det(odd_basis, prime) for prime in primes}
            require(any(basis_dets.values()), "complete even Fourier basis certified modulo a prime")
            require(any(odd_basis_dets.values()), "complete odd Fourier basis certified modulo a prime")
            # The test deliberately reports failure without declaring it a proof
            # of rational singularity; a positive value is an exact certificate.
            success = any(dets.values())
            report = {"n":n,"labelled_states":data["states"],"q":q,"even":len(even),"odd":len(odd),
                      "selfdual":len(fixed),"weighted_zero_modes":data["eigenvalues"].count(0),
                      "free_block_modular_determinants":dets,"even_basis_modular_determinants":basis_dets,
                      "odd_basis_modular_determinants":odd_basis_dets,
                      "free_block_nonsingular_certified":success,"full_data_sha256":digest(data),
                      "free_block_sha256":digest(a),"independent_full_matrix_python_control":n<=7,
                      "literal_permutation_representative_size_control":data["literal_permutations"],
                      "unrepaired_boolean_defect_rank":previous.exact_rank(previous.sp.Matrix(matrix_product(c,s))) if fixed else 0}
            if success:
                report.update(boolean_rank=2*len(odd), boolean_nullity=len(fixed))
            print("COMPLETE FOURIER STRESS", json.dumps(report,sort_keys=True),flush=True)
            reports.append(report)
    print("SEMANTIC SHA256",digest(reports))
    print("FINITE-EXACT n3..8; all-order transversality and native spectral gap remain OPEN")


if __name__ == "__main__":
    main()
