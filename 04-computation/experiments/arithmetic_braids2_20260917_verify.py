"""Replay all arithmetic-braids-II finite universes in normal and optimized Python.

Pins LF-normalized proof/source/certificate/output bytes in a single manifest.
No external Python packages are required.
"""
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
PREFIX = "arithmetic_braids2_20260917"
LANES = ("inverse_completion", "inverse_audit", "floor_reciprocity",
         "signed_cycles", "squarefree_symmetry", "fano_code")


def normalized(payload):
    return payload.replace(b"\r\n", b"\n")


def main():
    # Exact cross-lane lift exception; no floating-point or extrapolation.
    modulus = 23 ** 2
    ratio = 3 * pow(2, -1, modulus) % modulus
    if 3 ** 11 - 2 ** 11 != modulus * 331:
        raise RuntimeError("23-square factorization")
    if pow(ratio, 11, modulus) != 1 or any(pow(ratio, j, modulus) == 1 for j in range(1, 11)):
        raise RuntimeError("ratio order modulo529")
    for source in range(modulus):
        state = source
        for _ in range(11):
            state = (3 * state + 1) * pow(2, -1, modulus) % modulus
        if state != source:
            raise RuntimeError("affine branch period modulo529")
    records = []
    paths = {Path(__file__).relative_to(ROOT).as_posix(),
             f"05-knowledge/results/{PREFIX}_synthesis.md"}
    for lane in LANES:
        script = f"04-computation/experiments/{PREFIX}_{lane}.py"
        certificate = script[:-3] + ".json"
        transcript = f"05-knowledge/results/{PREFIX}_{lane}.out"
        regular = subprocess.run([sys.executable, script], cwd=ROOT,
                                 check=True, capture_output=True)
        cert_regular = normalized((ROOT / certificate).read_bytes())
        optimized = subprocess.run([sys.executable, "-O", script], cwd=ROOT,
                                   check=True, capture_output=True)
        cert_optimized = normalized((ROOT / certificate).read_bytes())
        if cert_regular != cert_optimized:
            raise RuntimeError(f"normal/optimized certificate differs: {lane}")
        if normalized(regular.stdout) != normalized(optimized.stdout):
            raise RuntimeError(f"normal/optimized transcript differs: {lane}")
        (ROOT / certificate).write_bytes(cert_regular)
        text = normalized(regular.stdout).decode("utf-8")
        text = text.replace(str(ROOT), ".").replace(ROOT.as_posix(), ".")
        (ROOT / transcript).write_text(text, encoding="utf-8", newline="\n")
        data = json.loads(cert_regular)
        if "source_sha256" in data:
            actual = sha256((ROOT / script).read_bytes()).hexdigest()
            if data["source_sha256"] != actual:
                raise RuntimeError(f"certificate source hash mismatch: {lane}")
        records.append({"lane": lane, "command": f"python {script}",
                        "ordinary_pass": True, "optimized_pass": True,
                        "certificates_identical_after_lf_normalization": True,
                        "certificate": certificate, "transcript": transcript})
        paths.update((script, certificate, transcript))
        if lane != "inverse_audit":
            paths.add(f"05-knowledge/results/{PREFIX}_{lane}.md")
        print(f"{lane}: normal + optimized + certificate comparison PASS", flush=True)
    artifacts = []
    for relative in sorted(paths):
        payload = normalized((ROOT / relative).read_bytes())
        artifacts.append({"path": relative, "sha256_lf": sha256(payload).hexdigest(),
                          "bytes_lf": len(payload)})
    manifest = {
        "status": "FINITE-EXACT reproductions; universal proofs in separately audited notes",
        "hash_basis": "UTF-8 artifact bytes normalized to LF",
        "runs": records, "artifacts": artifacts,
        "cross_lane_exact_control": {
            "identity": "3^11-2^11=23^2*331", "ratio_order_mod529": 11,
            "affine_branch_period_checks": 529,
            "scope": "The exponent-one affine branch modulo529, not integer orbit periodicity"},
        "independence": {
            "inverse_completion": "separate brute t search and closed-form carries in inverse_audit",
            "floor_reciprocity": "direct integer root and floor sums versus boundary formulas",
            "signed_cycles": "word enumeration versus direct forward trajectories and rational affine composition",
            "squarefree_symmetry": "direct square tests, sieve controls, HK/permutation/OCF Hamiltonian counts",
            "fano_code": "explicit finite code and lattice checks; theoretical audit recorded in note"},
        "open_boundaries": ["Global Collatz convergence", "Completeness of signed-cycle catalogs",
                            "Jointly squarefree terminating completions from the two separate constructions",
                            "A dynamical bridge from Fano/E8 or Bott periodicity to Collatz"]}
    destination = ROOT / f"05-knowledge/results/{PREFIX}_manifest.json"
    destination.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(f"manifest: {destination.relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()
