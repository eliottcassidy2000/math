# Bounded mathematical computation environment

Tested 2026-10-14 using Python 3.12.14, SymPy 1.14.0, mpmath 1.3.0,
NumPy 2.3.5 and SciPy 1.17.0. The exact package versions are pinned in
`requirements.txt`. These tools support symbolic, rational, multiprecision and
numerical experiments; successful smoke tests do not validate any repository
theorem or replace exact mathematical certificates.

## Reproduce

From the repository root, using a Python 3.12 installation on a POSIX shell:

```bash
python3 -m venv "$HOME/.venvs/math-computation-20261014"
. "$HOME/.venvs/math-computation-20261014/bin/activate"
python -m pip install -r 04-computation/environments/liouville-20261014/requirements.txt
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0 \
  python 04-computation/environments/liouville-20261014/smoke.py
```

The smoke script prints JSON and exits nonzero on a failed check. It exercises
rational arithmetic, symbolic factorization, exact Groebner reduction, a
hostile nonzero-remainder control, exact matrix inversion, numerical linear
algebra and 80-digit multiprecision arithmetic. All seven checks passed in the
observed environment. Numerical tolerances only test the numerical tools; they
are not proof certificates.

The observed run reused already installed packages through a virtual environment
created with `--system-site-packages`. A clean installation of the pinned
requirements above has not yet been replayed. The Python patch release and
package versions were recorded; this is not a full operating-system image or
wheel-hash lock.

## Optional tools and formal verification

GCC/G++ 14.2.0, GNU Make 4.4.1, GMP 6.3.0 and a LaTeX installation were available
in the observed environment. No SageMath, Julia, R, Lean, `lake` or `elan`
executable was found on PATH. They are not prerequisites for `smoke.py`.

The separate `04-computation/lean/TournamentH7` project pins Lean v4.30.0 and
Mathlib commit `c5ea00351c28e24afc9f0f84379aa41082b1188f`. Its documented
Mathlib cache estimate is approximately 6 GB. Lean installation, cache download,
`lake build TournamentH7.GMC2Main` and the endpoint `#print axioms` checks were
not performed for this environment checkpoint. Source inspection is not a
kernel build.

## Research hygiene

Review a computation's source and resource bounds before running it. Record its
exact universe, filters, positive and hostile controls, command, outputs and
source hashes separately from these toolchain smoke tests. Keep environment
files, caches and generated binaries outside versioned research artifacts.
