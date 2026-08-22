#!/usr/bin/env python3
"""Exact simple-spectrum audit for THM-3660's signed exceptional detector."""

from __future__ import annotations

import ast
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_SCRIPT = ROOT / "04-computation/lrc_exceptional_leakage_boundary_thm3660.py"
PARENT_OUTPUT = ROOT / "05-knowledge/results/lrc_exceptional_leakage_boundary_thm3660.out"
EXPECTED_PARENT_HASHES = (
    "37a7c775096d80c4d8acd9d39fc0812e32545a963ee572ea773157d0ef6804c6",
    "d24cda381808b5ec6009a47fc967f1293486e3ffb95d9153f4e4cd0bd3ed0730",
)
EXPECTED_SEMANTIC_SHA256 = "c9ffab0b789d087cc094e137274c6a17a2ef4c06aaa165a6969d77e404bc8ec9"

P = 13
N = P * P
MOD = 755373809845391722745761
ZETA = 123453826432109539554819
OMEGA = 298763986285447441216949
X_PLUS = frozenset(((12, 0), (0, 11), (6, 5), (9, 3)))
X_MINUS = frozenset(((0, 12), (12, 1), (6, 7), (3, 9)))
EXPECTED_PRODUCT = 669422013050837354847410
EXPECTED_CYCLIC_PRODUCT = 540112521324494664145058


def require(condition: bool, payload: object) -> None:
    if condition is not True:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def add(left, right):
    return ((left[0] + right[0]) % P, (left[1] + right[1]) % P)


def subtract(left, right):
    return ((left[0] - right[0]) % P, (left[1] - right[1]) % P)


def main() -> None:
    parent_hashes = (lf_sha256(PARENT_SCRIPT), lf_sha256(PARENT_OUTPUT))
    require(parent_hashes == EXPECTED_PARENT_HASHES,
            ("THM-3660 parent hashes", parent_hashes))
    require(pow(ZETA, P, MOD) == 1 and ZETA != 1, "zeta order")
    require(pow(OMEGA, N, MOD) == 1 and pow(OMEGA, P, MOD) != 1,
            "omega order")
    labels = tuple((r0, r1) for r0 in range(P) for r1 in range(P))
    detector = {
        label: (int(label in X_PLUS) - int(label in X_MINUS)) % MOD
        for label in labels
    }
    require(sum(detector.values()) % MOD == 0, "detector mean")

    fourier = {}
    for frequency in labels:
        total = 0
        for label in labels:
            exponent = -(frequency[0] * label[0] + frequency[1] * label[1]) % P
            total = (total + detector[label] * pow(ZETA, exponent, MOD)) % MOD
        fourier[frequency] = total
    zero_frequencies = tuple(frequency for frequency in labels
                             if fourier[frequency] == 0)
    require(zero_frequencies == ((0, 0),),
            ("Fourier zero set", zero_frequencies))
    nontrivial_values = tuple(fourier[frequency] for frequency in labels
                              if frequency != (0, 0))
    require(len(set(nontrivial_values)) == N - 1,
            "nontrivial Fourier eigenvalues not simple")
    product = 1
    for value in nontrivial_values:
        product = product * value % MOD
    require(product == EXPECTED_PRODUCT, ("augmentation determinant", product))

    cyclic_fourier = {}
    for frequency in range(N):
        total = 0
        for label in labels:
            assembled = label[0] + P * label[1]
            total = (total + detector[label]
                     * pow(OMEGA, -frequency * assembled % N, MOD)) % MOD
        cyclic_fourier[frequency] = total
    cyclic_zero_frequencies = tuple(
        frequency for frequency in range(N)
        if cyclic_fourier[frequency] == 0
    )
    require(cyclic_zero_frequencies == (0,),
            ("cyclic Fourier zero set", cyclic_zero_frequencies))
    cyclic_nontrivial_values = tuple(cyclic_fourier[k] for k in range(1, N))
    require(len(set(cyclic_nontrivial_values)) == N - 1,
            "cyclic Fourier eigenvalues not simple")
    cyclic_product = 1
    for value in cyclic_nontrivial_values:
        cyclic_product = cyclic_product * value % MOD
    require(cyclic_product == EXPECTED_CYCLIC_PRODUCT,
            ("cyclic augmentation determinant", cyclic_product))

    # Characteristic polynomial of convolution by the detector, stored with
    # coefficients in ascending degree order.
    characteristic = [1]
    for frequency in labels:
        eigenvalue = fourier[frequency]
        updated = [0] * (len(characteristic) + 1)
        for degree, coefficient in enumerate(characteristic):
            updated[degree] = (updated[degree] - eigenvalue * coefficient) % MOD
            updated[degree + 1] = (updated[degree + 1] + coefficient) % MOD
        characteristic = updated
    require(len(characteristic) == N + 1
            and characteristic[-1] == 1
            and characteristic[0] == 0,
            "characteristic polynomial shape")
    derivative = tuple(
        (degree * characteristic[degree]) % MOD
        for degree in range(1, len(characteristic))
    )
    # A direct eigenvalue check certifies square-freeness.
    require(all(sum(derivative[degree] * pow(value, degree, MOD)
                    for degree in range(len(derivative))) % MOD
                for value in fourier.values()),
            "characteristic polynomial repeated root")

    cyclic_characteristic = [1]
    for frequency in range(N):
        eigenvalue = cyclic_fourier[frequency]
        updated = [0] * (len(cyclic_characteristic) + 1)
        for degree, coefficient in enumerate(cyclic_characteristic):
            updated[degree] = (updated[degree] - eigenvalue * coefficient) % MOD
            updated[degree + 1] = (updated[degree + 1] + coefficient) % MOD
        cyclic_characteristic = updated
    cyclic_derivative = tuple(
        (degree * cyclic_characteristic[degree]) % MOD
        for degree in range(1, len(cyclic_characteristic))
    )
    require(len(cyclic_characteristic) == N + 1
            and cyclic_characteristic[-1] == 1
            and cyclic_characteristic[0] == 0
            and all(sum(cyclic_derivative[degree] * pow(value, degree, MOD)
                        for degree in range(len(cyclic_derivative))) % MOD
                    for value in cyclic_fourier.values()),
            "cyclic characteristic polynomial shape/squarefreeness")

    # The inverse multiplier on the augmentation ideal.
    inverse_n = pow(N, -1, MOD)
    inverse_kernel = {}
    for label in labels:
        total = 0
        for frequency in labels:
            if frequency == (0, 0):
                continue
            phase = pow(ZETA,
                        (frequency[0] * label[0] + frequency[1] * label[1]) % P,
                        MOD)
            total = (total + pow(fourier[frequency], -1, MOD) * phase) % MOD
        inverse_kernel[label] = total * inverse_n % MOD
    require(tuple(label for label in labels if inverse_kernel[label] == 0)
            == ((7, 7),), "inverse-kernel zero")
    require(len(set(inverse_kernel.values())) == N,
            "inverse-kernel address separation")

    convolution = {}
    for label in labels:
        convolution[label] = sum(
            inverse_kernel[left] * detector[subtract(label, left)]
            for left in labels
        ) % MOD
    target = {
        label: (int(label == (0, 0)) - inverse_n) % MOD
        for label in labels
    }
    require(convolution == target, "augmentation inverse identity")

    cyclic_detector = {
        assembled: detector[(assembled % P, assembled // P)]
        for assembled in range(N)
    }
    cyclic_inverse_kernel = {}
    for assembled in range(N):
        total = sum(
            pow(cyclic_fourier[frequency], -1, MOD)
            * pow(OMEGA, frequency * assembled % N, MOD)
            for frequency in range(1, N)
        ) % MOD
        cyclic_inverse_kernel[assembled] = total * inverse_n % MOD
    require(tuple(a for a in range(N) if cyclic_inverse_kernel[a] == 0)
            == (85,), "cyclic inverse-kernel zero")
    require(len(set(cyclic_inverse_kernel.values())) == N,
            "cyclic inverse-kernel address separation")
    cyclic_convolution = {
        assembled: sum(
            cyclic_inverse_kernel[left]
            * cyclic_detector[(assembled - left) % N]
            for left in range(N)
        ) % MOD
        for assembled in range(N)
    }
    cyclic_target = {
        assembled: (int(assembled == 0) - inverse_n) % MOD
        for assembled in range(N)
    }
    require(cyclic_convolution == cyclic_target,
            "cyclic augmentation inverse identity")

    vertical = {
        label: (detector[add(label, (0, 1))] - detector[label]) % MOD
        for label in labels
    }
    vertical_fourier = {}
    for frequency in labels:
        direct = sum(
            vertical[label]
            * pow(ZETA,
                  -(frequency[0] * label[0] + frequency[1] * label[1]) % P,
                  MOD)
            for label in labels
        ) % MOD
        multiplier = (pow(ZETA, frequency[1], MOD) - 1) % MOD
        require(direct == multiplier * fourier[frequency] % MOD,
                ("vertical multiplier", frequency))
        vertical_fourier[frequency] = direct
    vertical_support = tuple(frequency for frequency in labels
                             if vertical_fourier[frequency])
    require(vertical_support == tuple((u, v) for u in range(P)
                                      for v in range(1, P)),
            "vertical Fourier support")

    # Translation span: a vector is orthogonal to every translate of g iff
    # its Fourier support is confined to the trivial character.
    translate_rows = tuple(
        tuple(detector[subtract(label, shift)] for label in labels)
        for shift in labels
    )
    # Rank through character support is 168; check a direct modular RREF too.
    matrix = [list(row) for row in translate_rows]
    rank = 0
    for column in range(N):
        pivot = next((row for row in range(rank, N)
                      if matrix[row][column] % MOD), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, MOD)
        matrix[rank] = [value * inverse % MOD for value in matrix[rank]]
        for row in range(N):
            if row != rank and matrix[row][column]:
                factor = matrix[row][column]
                matrix[row] = [
                    (left - factor * right) % MOD
                    for left, right in zip(matrix[row], matrix[rank])
                ]
        rank += 1
        if rank == N:
            break
    require(rank == N - 1, ("translate rank", rank))

    semantic = digest_json((
        MOD, P, ZETA, OMEGA, parent_hashes,
        tuple(sorted(X_PLUS)), tuple(sorted(X_MINUS)),
        tuple((frequency, fourier[frequency]) for frequency in labels),
        zero_frequencies, product, tuple(characteristic),
        tuple((frequency, cyclic_fourier[frequency]) for frequency in range(N)),
        cyclic_zero_frequencies, cyclic_product, tuple(cyclic_characteristic),
        tuple((label, inverse_kernel[label]) for label in labels),
        tuple((label, convolution[label]) for label in labels),
        tuple((assembled, cyclic_inverse_kernel[assembled])
              for assembled in range(N)),
        tuple((frequency, vertical_fourier[frequency]) for frequency in labels),
        vertical_support,
        rank,
    ))
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source.read_text(encoding="utf-8")))),
            "Python assert node present")
    print("== THM-3661 LRC exceptional detector simple-spectrum rigidity ==")
    print(f"field=p:{MOD};groups:F13^2,C169;characters:169;zeta:{ZETA};omega:{OMEGA}")
    print(f"parent_sha256_lf={parent_hashes}")
    print("detector=signed indicator X_plus-X_minus;support:8;mean:0")
    print(f"split_fourier_zero_set={zero_frequencies};nontrivial_support:168;distinct:168;determinant:{product}")
    print(f"split_fourier_sha256={digest_json(tuple((frequency, fourier[frequency]) for frequency in labels))}")
    print(f"split_characteristic_degree=169;squarefree=True;sha256={digest_json(tuple(characteristic))}")
    print(f"cyclic_fourier_zero_set={cyclic_zero_frequencies};nontrivial_support:168;distinct:168;determinant:{cyclic_product}")
    print(f"cyclic_fourier_sha256={digest_json(tuple((frequency, cyclic_fourier[frequency]) for frequency in range(N)))}")
    print(f"cyclic_characteristic_degree=169;squarefree=True;sha256={digest_json(tuple(cyclic_characteristic))}")
    print("both_convolution_spectra=simple;centralizers=F[C_g],degree<169")
    print(f"split_translate_span=augmentation_ideal;rank:{rank}")
    print(f"vertical_derivative_fourier_support=156;sha256={digest_json(tuple((frequency, vertical_fourier[frequency]) for frequency in labels))}")
    print(f"split_inverse_kernel=support:168;distinct_values:169;zero:((7,7),);sha256={digest_json(tuple((label, inverse_kernel[label]) for label in labels))}")
    print(f"cyclic_inverse_kernel=support:168;distinct_values:169;zero:(85,);sha256={digest_json(tuple((assembled, cyclic_inverse_kernel[assembled]) for assembled in range(N)))}")
    print("inverse_identities=h*g=delta_0-(1/169)1 in both group laws")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("scope=finite group algebra and cyclotomic detector;not physical translation/current/LRC14")
    print("PASS")


if __name__ == "__main__":
    main()
