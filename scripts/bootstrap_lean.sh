#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

export ELAN_HOME="${ELAN_HOME:-/root/.elan}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/root/.cache}"
export LEAN_SCRATCH="${LEAN_SCRATCH:-/workspace/.scratch/lean}"
export PATH="${ELAN_HOME}/bin:${PATH}"

echo "==> Preparing persistent Lean directories"
mkdir -p "${ELAN_HOME}" "${XDG_CACHE_HOME}" "${LEAN_SCRATCH}"

if ! command -v elan >/dev/null 2>&1; then
  echo "==> Installing elan"
  curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --no-modify-path
else
  echo "==> elan already installed"
fi

if ! command -v lean >/dev/null 2>&1 || ! command -v lake >/dev/null 2>&1; then
  echo "==> Installing default stable Lean toolchain"
  elan default stable
fi

ENV_FILE="${REPO_ROOT}/.scratch/lean/env.sh"
mkdir -p "$(dirname "${ENV_FILE}")"
cat > "${ENV_FILE}" <<EOF
export ELAN_HOME="${ELAN_HOME}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME}"
export LEAN_SCRATCH="${LEAN_SCRATCH}"
export PATH="${ELAN_HOME}/bin:\${PATH}"
EOF

echo "==> Verifying Lean tools"
elan --version
lean --version
lake --version

echo
echo "Lean environment ready."
echo "Scratch directory: ${LEAN_SCRATCH}"
echo "Reusable shell env: ${ENV_FILE}"
