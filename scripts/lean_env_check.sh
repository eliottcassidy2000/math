#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

export ELAN_HOME="${ELAN_HOME:-/root/.elan}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/root/.cache}"
export LEAN_SCRATCH="${LEAN_SCRATCH:-/workspace/.scratch/lean}"
export PATH="${ELAN_HOME}/bin:${PATH}"

SMOKE_TEST=0
if [[ "${1:-}" == "--smoke-test" ]]; then
  SMOKE_TEST=1
fi

echo "Repo root: ${REPO_ROOT}"
echo "ELAN_HOME: ${ELAN_HOME}"
echo "XDG_CACHE_HOME: ${XDG_CACHE_HOME}"
echo "Scratch: ${LEAN_SCRATCH}"
echo "PATH: ${PATH}"
echo

missing=0
for tool in elan lean lake; do
  if ! command -v "${tool}" >/dev/null 2>&1; then
    echo "${tool}: MISSING"
    missing=1
  else
    echo "${tool}: $(command -v "${tool}")"
    "${tool}" --version
  fi
  echo
done

if [[ "${missing}" -ne 0 ]]; then
  echo "Missing Lean tools. Run scripts/bootstrap_lean.sh first." >&2
  exit 1
fi

if [[ "${SMOKE_TEST}" -eq 1 ]]; then
  mkdir -p "${LEAN_SCRATCH}"
  smoke_file="${LEAN_SCRATCH}/smoke.lean"
  printf 'theorem lean_bootstrap_smoke : True := by\n  trivial\n' > "${smoke_file}"
  echo "Running smoke test: lean ${smoke_file}"
  lean "${smoke_file}"
fi
