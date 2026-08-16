#!/usr/bin/env python3
"""Independent audit of the repaired THM-2594 r=5 theta contraction.

This companion deliberately separates two questions.

1. Reproducibility: run the primary exact interval program in normal and
   optimized modes and compare both transcripts with the stored output after
   masking only elapsed-time fields.
2. Consequence-object audit: use the primary program only to materialize its
   exact common-base table, then independently reimplement ANOVA centring,
   the toothpick transform, Phi_91 reduction, the THM-2512 factorization, all
   5,184 primitive coefficients, and the three attribution controls.

The linked-node identities are also checked directly with Fraction arithmetic.
They put the owner, deep/word, and source factors over one outer base y and one
finite ancestry sum.  They do not put those factors at one circle point and do
not construct a physical current.
"""

from __future__ import annotations

import hashlib
import importlib.util
import io
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from contextlib import redirect_stdout
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / "04-computation" / "lrc14_stage2_theta_contraction_opus_20260728.py"
STORED = ROOT / "05-knowledge" / "results" / "lrc14_stage2_theta_contraction_opus_20260728.out"

PRIMARY_LF_SHA256 = "09c43af0a0a56c7a0833bbfd13ed6a96bc5a7a3718aa1bc6b77a144bde101a06"
STORED_LF_SHA256 = "bef4ee9a18ff3e2f455bad66a95252dd9989b2f60953e26e8ea0c2dc6ae7f5df"
REPLAY_SEMANTIC_SHA256 = "184c159dd4aedb02b38b2029846b9e88f301fcc63dd163682bbabddc438b825e"
HISTORICAL_COMMIT = "c04aca388"
HISTORICAL_SCRIPT_REF = (
    HISTORICAL_COMMIT
    + ":04-computation/lrc14_stage2_theta_contraction_opus_20260728.py"
)
HISTORICAL_OUTPUT_REF = (
    HISTORICAL_COMMIT
    + ":05-knowledge/results/lrc14_stage2_theta_contraction_opus_20260728.out"
)
HISTORICAL_SCRIPT_LF_SHA256 = "aa03adf52f72379ede708e2040925b088548a7cc05fa999e54078e1a8e06e9e1"
HISTORICAL_OUTPUT_LF_SHA256 = "6cf260afdfc41e0e560adf294b2cf31b50148a67bbb4b4085195bce9c43dd149"
HISTORICAL_REPLAY_SEMANTIC_SHA256 = "d13f0b23768b74d41d4696c6e79e1597ff6e1e8accded7b504897f47245bbb46"

P = 13
Q = 7
D = 13**5
R_PACKET = 13**2
C3 = 2 * 13**5
T_DEN = 297836897838480
DEN_S_EXPECTED = 631451675259549179940893068080


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_bytes(path: Path) -> bytes:
    return path.read_bytes().replace(b"\r\n", b"\n")


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


ELAPSED_RE = re.compile(br"\[\s*\d+\.\d+s\]")


def normalize_replay(data: bytes) -> bytes:
    return ELAPSED_RE.sub(b"[ELAPSED]", data.replace(b"\r\n", b"\n"))


def load_primary_module():
    spec = importlib.util.spec_from_file_location("thm2594_primary", PRIMARY)
    require(spec is not None and spec.loader is not None, "cannot load primary module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_optimized() -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        [sys.executable, "-B", "-O", str(PRIMARY)],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def git_blob(ref: str) -> bytes:
    return subprocess.check_output(["git", "show", ref], cwd=ROOT)


def run_historical(optimized: bool) -> subprocess.CompletedProcess[bytes]:
    """Execute the immutable pre-repair script directly from its Git blob."""
    runner = (
        "import subprocess,sys;"
        "code=subprocess.check_output(['git','show',sys.argv[1]]);"
        "exec(compile(code,sys.argv[1],'exec',optimize=sys.flags.optimize),"
        "{'__name__':'__main__'})"
    )
    command = [sys.executable, "-B"]
    if optimized:
        command.append("-O")
    command.extend(["-c", runner, HISTORICAL_SCRIPT_REF])
    return subprocess.run(
        command,
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def poly_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def divide_by_monic(dividend: list[int], divisor: list[int]) -> list[int]:
    """Fresh high-to-low exact division, independent of the primary helpers."""
    require(divisor and divisor[-1] == 1, "divisor is not monic")
    rem = dividend[:]
    degree = len(divisor) - 1
    quotient = [0] * (len(dividend) - degree)
    for top in range(len(rem) - 1, degree - 1, -1):
        coefficient = rem[top]
        quotient[top - degree] = coefficient
        if coefficient:
            for j, value in enumerate(divisor):
                rem[top - degree + j] -= coefficient * value
    require(not any(rem), "nonzero polynomial remainder")
    return quotient


def phi91() -> list[int]:
    numerator = poly_mul([-1] + [0] * 90 + [1], [-1, 1])
    denominator = poly_mul([-1] + [0] * 12 + [1], [-1] + [0] * 6 + [1])
    value = divide_by_monic(numerator, denominator)
    require(len(value) == 73 and value[-1] == 1, "bad Phi_91")
    return value


PHI91 = phi91()


def reduce_phi91(poly: list[int]) -> tuple[int, ...]:
    work = poly[:] + [0] * max(0, 73 - len(poly))
    for top in range(len(work) - 1, 71, -1):
        coefficient = work[top]
        if coefficient:
            for j, value in enumerate(PHI91):
                work[top - 72 + j] -= coefficient * value
    return tuple(work[:72])


def exponent(i13: int, i7: int) -> int:
    return (7 * (i13 % 13) + 13 * (i7 % 7)) % 91


def cyclic_convolution_91(a: list[int], b: list[int]) -> list[int]:
    out = [0] * 91
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    out[(i + j) % 91] += ai * bj
    return out


def interaction_numerators(array: list[list[int]]) -> list[list[int]]:
    """Return numerators of 91 times the doubly centred 7 by 13 table."""
    require(len(array) == 7 and all(len(row) == 13 for row in array), "bad array shape")
    row_sums = [sum(row) for row in array]
    col_sums = [sum(array[ell][theta] for ell in range(7)) for theta in range(13)]
    total = sum(row_sums)
    result = [
        [
            91 * array[ell][theta]
            - 7 * row_sums[ell]
            - 13 * col_sums[theta]
            + total
            for theta in range(13)
        ]
        for ell in range(7)
    ]
    require(all(sum(row) == 0 for row in result), "centred row sum is nonzero")
    require(
        all(sum(result[ell][theta] for ell in range(7)) == 0 for theta in range(13)),
        "centred column sum is nonzero",
    )
    return result


def psi_raw(
    defect: list[list[int]], tau: int, a0: int, alpha: int, beta: int
) -> list[int]:
    """Directly expand THM-2512 (12)--(14) in the x=zeta_91 basis."""
    raw = [0] * 91
    for cut in range(7):
        for output in range(13):
            value = sum(
                defect[ell][
                    (output - tau * ((a0 * ell + cut) % 7)) % 13
                ]
                for ell in range(7)
            )
            raw[exponent(-alpha * output, -beta * cut)] += value
    return raw


def dtilde_raw(defect: list[list[int]], alpha: int, gamma: int) -> list[int]:
    raw = [0] * 91
    for ell in range(7):
        for output in range(13):
            raw[exponent(-alpha * output, -gamma * ell)] += defect[ell][output]
    return raw


def k_raw(u: int, beta: int) -> list[int]:
    raw = [0] * 91
    for j in range(7):
        raw[exponent(-u * j, -beta * j)] += 1
    return raw


def reduced_psi(
    defect: list[list[int]], tau: int, a0: int, alpha: int, beta: int
) -> tuple[int, ...]:
    return reduce_phi91(psi_raw(defect, tau, a0, alpha, beta))


def is_nonzero(vector: tuple[int, ...]) -> bool:
    return any(vector)


def audit_all_primitive(defect: list[list[int]]) -> tuple[int, str, int]:
    count = 0
    nonzero_coordinate_floor = 72
    digest = hashlib.sha256()
    for tau in range(1, 13):
        for a0 in range(1, 7):
            for alpha in range(1, 13):
                for beta in range(1, 7):
                    direct = reduced_psi(defect, tau, a0, alpha, beta)
                    factored = reduce_phi91(
                        cyclic_convolution_91(
                            k_raw(alpha * tau, beta),
                            dtilde_raw(defect, alpha, -beta * a0),
                        )
                    )
                    require(direct == factored, "factorization mismatch")
                    if is_nonzero(direct):
                        count += 1
                    nonzero_coordinate_floor = min(
                        nonzero_coordinate_floor, sum(value != 0 for value in direct)
                    )
                    digest.update(
                        (
                            f"{tau},{a0},{alpha},{beta}:"
                            + ",".join(map(str, direct))
                            + "\n"
                        ).encode()
                    )
    return count, digest.hexdigest(), nonzero_coordinate_floor


def distance_to_integer(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def danger(value: Fraction) -> bool:
    return distance_to_integer(value) < Fraction(1, 14)


def expected_deep(theta: int, y: Fraction) -> bool:
    if theta == 0:
        return Fraction(0) <= y < Fraction(13, 28)
    if theta == 1:
        return Fraction(1, 28) < y < Fraction(27, 28)
    if theta == 2:
        return Fraction(15, 28) < y < 1
    return False


def same_stalk_audit() -> tuple[int, int]:
    """Check the exact linked-node descents and deep-window partition."""
    y_samples = (
        Fraction(1, 100),
        Fraction(17, 101),
        Fraction(1, 2),
        Fraction(99, 100),
    )
    sheet_samples = (0, 1, 2, 12345, D - 1)
    checks = 0
    for y in y_samples:
        for u in range(13):
            w_u = (y + u) / 13
            require(13 * w_u - y == u, "owner node does not descend to y")
            for a in sheet_samples:
                x_ua = (w_u + a) / D
                require(C3 * x_ua - 2 * (y + u) / 13 == 2 * a, "deep sheet remains")
                require(13**6 * x_ua - y == u + 13 * a, "word sheet remains")
                for t in range(13):
                    theta = (t - 2 * u) % 13
                    require(
                        danger(C3 * x_ua - Fraction(t, 13))
                        == danger(Fraction(2 * y - theta, 13)),
                        "theta=t-2u descent fails",
                    )
                    checks += 1
            for q in range(13):
                w_q = (y + q) / 13
                for eprime in (0, 1, D - 1):
                    y_qe = (w_q + eprime) / D
                    require(D * y_qe - w_q == eprime, "source sheet identity fails")

    breakpoints = sorted(
        {Fraction(0), Fraction(1, 28), Fraction(13, 28), Fraction(15, 28), Fraction(27, 28), Fraction(1)}
    )
    cells = [(a + b) / 2 for a, b in zip(breakpoints, breakpoints[1:])]
    deep_checks = 0
    for y in cells:
        values = [danger(Fraction(2 * y - theta, 13)) for theta in range(13)]
        require(
            values == [expected_deep(theta, y) for theta in range(13)],
            "deep-window description fails",
        )
        require(sum(values) == 2 - int(danger(2 * y)), "deep root-count law fails")
        deep_checks += 13
    return checks, deep_checks


def split_prime_audit(reduced: tuple[int, ...]) -> tuple[int, ...]:
    rows = ((547, 64, 295), (911, 113, 213), (1093, 817, 380), (2003, 522, 1405), (2549, 266, 883))
    values = []
    for prime, root, expected in rows:
        require(pow(root, 91, prime) == 1, "root does not have order dividing 91")
        require(pow(root, 7, prime) != 1 and pow(root, 13, prime) != 1, "root order is proper")
        value = 0
        for coefficient in reversed(reduced):
            value = (value * root + coefficient) % prime
        require(value == expected, "split-prime value mismatch")
        values.append(value)
    return tuple(values)


def main() -> None:
    require(sha256(lf_bytes(PRIMARY)) == PRIMARY_LF_SHA256, "primary hash drift")
    require(sha256(lf_bytes(STORED)) == STORED_LF_SHA256, "stored hash drift")
    stored_normalized = normalize_replay(STORED.read_bytes())
    require(sha256(stored_normalized) == REPLAY_SEMANTIC_SHA256, "stored semantic hash drift")
    historical_script = git_blob(HISTORICAL_SCRIPT_REF).replace(b"\r\n", b"\n")
    historical_output = git_blob(HISTORICAL_OUTPUT_REF).replace(b"\r\n", b"\n")
    require(sha256(historical_script) == HISTORICAL_SCRIPT_LF_SHA256, "historical script drift")
    require(sha256(historical_output) == HISTORICAL_OUTPUT_LF_SHA256, "historical output drift")
    historical_normalized = normalize_replay(historical_output)
    require(
        sha256(historical_normalized) == HISTORICAL_REPLAY_SEMANTIC_SHA256,
        "historical semantic hash drift",
    )

    with ThreadPoolExecutor(max_workers=3) as pool:
        optimized_future = pool.submit(run_optimized)
        historical_normal_future = pool.submit(run_historical, False)
        historical_optimized_future = pool.submit(run_historical, True)
        primary = load_primary_module()
        normal_buffer = io.StringIO()
        with redirect_stdout(normal_buffer):
            state = primary.main()
            primary.stage2(state)
        optimized = optimized_future.result()
        historical_normal = historical_normal_future.result()
        historical_optimized = historical_optimized_future.result()

    normal_output = normal_buffer.getvalue().encode()
    require(normalize_replay(normal_output) == stored_normalized, "normal replay mismatch")
    require(optimized.returncode == 0, "optimized replay failed")
    require(not optimized.stderr, "optimized replay wrote stderr")
    require(normalize_replay(optimized.stdout) == stored_normalized, "optimized replay mismatch")
    for label, replay in (
        ("historical normal", historical_normal),
        ("historical optimized", historical_optimized),
    ):
        require(replay.returncode == 0, f"{label} replay failed")
        require(not replay.stderr, f"{label} replay wrote stderr")
        require(
            normalize_replay(replay.stdout) == historical_normalized,
            f"{label} replay mismatch",
        )

    (
        _e,
        _len_e,
        _qb,
        _tb,
        _cells,
        _deep,
        _d2,
        _qw,
        _gw,
        _g,
        _gcd2,
        pair_totals,
        joint,
        _joint_word,
        _joint_word_d2,
        colours,
        _fourier_colours,
        common_denominator,
        total_service,
    ) = state

    require(common_denominator == R_PACKET * D * D * T_DEN, "common denominator typing")
    require(common_denominator * 91 == DEN_S_EXPECTED, "decisive denominator drift")
    require(all(value >= 0 for row in pair_totals for value in row), "negative pair mass")
    require(all(pair_totals[u][u] == 0 for u in range(13)), "diagonal mass survives")
    require(colours[0] == 0 and sum(colours) == total_service, "colour/service mismatch")

    response = [[0] * 13 for _ in range(7)]
    for u in range(13):
        for q in range(13):
            for ell in range(7):
                for theta in range(3):
                    value = joint[u][q][ell][theta]
                    require(value >= 0, "negative Boolean-fibre mass")
                    response[ell][theta] += value
    require(any(response[ell][theta] for ell in range(7) for theta in range(3)), "zero response")

    defect = interaction_numerators(response)
    decisive = reduced_psi(defect, 1, 1, 1, 1)
    require(is_nonzero(decisive), "decisive coefficient vanished")
    split_values = split_prime_audit(decisive)
    primitive_count, primitive_digest, primitive_coordinate_floor = audit_all_primitive(defect)
    require(primitive_count == 5184, "not all primitive coefficients survive")

    beta_zero = reduced_psi(defect, 1, 1, 1, 0)
    require(not is_nonzero(beta_zero), "beta-zero hostile survives")

    flat_response = [[sum(response[ell])] * 13 for ell in range(7)]
    flat_defect = interaction_numerators(flat_response)
    require(not any(value for row in flat_defect for value in row), "flat interaction survives")
    require(not is_nonzero(reduced_psi(flat_defect, 1, 1, 1, 1)), "flat hostile survives")

    absolute_response = [[0] * 13 for _ in range(7)]
    for u in range(13):
        for q in range(13):
            for ell in range(7):
                for absolute_root in range(13):
                    theta = (absolute_root - 2 * u) % 13
                    if theta < 3:
                        absolute_response[ell][absolute_root] += joint[u][q][ell][theta]
    absolute_defect = interaction_numerators(absolute_response)
    absolute_decisive = reduced_psi(absolute_defect, 1, 1, 1, 1)
    require(is_nonzero(absolute_decisive), "fixed-root hostile vanished")
    absolute_primitive_count, absolute_digest, absolute_coordinate_floor = audit_all_primitive(
        absolute_defect
    )
    require(absolute_primitive_count == 5184, "fixed-root primitive bundle has a zero")

    fibre_count = 0
    fibre_by_shift = [0] * 13
    fibre_sum = [[0] * 13 for _ in range(7)]
    for u in range(13):
        for shift in range(13):
            q = (u - shift) % 13
            fibre_response = [
                [joint[u][q][ell][theta] if theta < 3 else 0 for theta in range(13)]
                for ell in range(7)
            ]
            fibre_defect = interaction_numerators(fibre_response)
            if any(value for row in fibre_defect for value in row):
                fibre_count += 1
                fibre_by_shift[shift] += 1
            for ell in range(7):
                for theta in range(13):
                    fibre_sum[ell][theta] += fibre_defect[ell][theta]
    require(fibre_sum == defect, "fibrewise centring does not sum to aggregate")
    require(fibre_count == 34 and fibre_by_shift[0] == 0, "fibre census drift")

    stalk_checks, deep_checks = same_stalk_audit()

    print("THM-2594 R=5 THETA CONTRACTION -- INDEPENDENT AUDIT")
    print("status=ACCEPT_NARROW_REALIZED_CANDIDATE_CONTRACTION")
    print(f"primary_lf_sha256={PRIMARY_LF_SHA256}")
    print(f"stored_lf_sha256={STORED_LF_SHA256}")
    print(f"replay_semantic_sha256={REPLAY_SEMANTIC_SHA256}")
    print(f"historical_commit={HISTORICAL_COMMIT}")
    print(f"historical_primary_lf_sha256={HISTORICAL_SCRIPT_LF_SHA256}")
    print(f"historical_stored_lf_sha256={HISTORICAL_OUTPUT_LF_SHA256}")
    print(f"historical_replay_semantic_sha256={HISTORICAL_REPLAY_SEMANTIC_SHA256}")
    print("historical_normal_replay=PASS")
    print("historical_optimized_replay=PASS")
    print("historical_stored_replay=PASS")
    print("historical_causal_wording=SUPERSEDED_BY_MISTAKE_295")
    print("normal_replay=PASS")
    print("optimized_replay=PASS")
    print("stored_replay=PASS")
    print(f"common_denominator={common_denominator}")
    print(f"decisive_denominator={common_denominator * 91}")
    print(f"decisive_nonzero={is_nonzero(decisive)}")
    print(f"decisive_reduced_sha256={sha256(','.join(map(str, decisive)).encode())}")
    print(f"decisive_nonzero_coordinates={sum(value != 0 for value in decisive)}/72")
    print(f"split_prime_values={split_values}")
    print(f"primitive_nonzero={primitive_count}/5184")
    print(f"primitive_bundle_sha256={primitive_digest}")
    print(f"primitive_nonzero_coordinate_floor={primitive_coordinate_floor}/72")
    print(f"fibrewise_interaction_nonzero={fibre_count}/169")
    print(f"fibrewise_by_shift={tuple(fibre_by_shift)}")
    print(f"same_stalk_fraction_checks={stalk_checks}")
    print(f"deep_window_checks={deep_checks}")
    print("owner_node=13*w_u-y=u_in_Z")
    print("deep_node=c3*X_u_a-2*(y+u)/13=2*a_in_Z")
    print("word_node=13^6*X_u_a-y=u+13*a_in_Z")
    print("linked_nodes_not_one_circle_point=PASS")
    print("h1_beta_zero=PASS_ZERO")
    print("h2_constant_column=PASS_ZERO")
    print(f"h3_fixed_absolute_root_raw_nonzero={sum(value != 0 for row in absolute_response for value in row)}/91")
    print(f"h3_fixed_absolute_root_centred_nonzero={sum(value != 0 for row in absolute_defect for value in row)}/91")
    print(f"h3_fixed_absolute_root_decisive_nonzero={is_nonzero(absolute_decisive)}")
    print(f"h3_fixed_absolute_root_primitive_nonzero={absolute_primitive_count}/5184")
    print(f"h3_fixed_absolute_root_bundle_sha256={absolute_digest}")
    print(f"h3_fixed_absolute_root_coordinate_floor={absolute_coordinate_floor}/72")
    print("causal_slaving_uniqueness=REFUTED_BY_H3")
    print("thm2512_generic_bridge=NOT_PROVED")
    print("thm2449_one_point_response=NOT_PROVED")
    print("physical_current=NOT_PROVED")
    print("row_exclusion=0")
    print("lrc14=OPEN")


if __name__ == "__main__":
    main()
