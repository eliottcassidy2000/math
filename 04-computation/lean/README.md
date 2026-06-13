# Lean Environment Bootstrap

This directory documents the Lean setup for future theorem-formalization work.
The repository does not yet define a full Lean/mathlib project; use the
bootstrap scripts to make `elan`, `lean`, and `lake` available first.

## Windows Workspace

From the repository root:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/bootstrap_lean.ps1
powershell -ExecutionPolicy Bypass -File scripts/lean_env_check.ps1 -SmokeTest
```

The local scratch directory is `.scratch/lean/`. It is intentionally ignored by
git and is safe for throwaway Lean experiments.

The Windows bootstrap installs Lean 4.30.0 by default. It uses `aria2c` for a
segmented download when available, which is much faster for the large Windows
Lean release archive.

## Linux / Cloud Codex Images

From the repository root:

```bash
bash scripts/bootstrap_lean.sh
bash scripts/lean_env_check.sh --smoke-test
```

The cloud bootstrap uses these persistent locations by default:

- `ELAN_HOME=/root/.elan`
- `XDG_CACHE_HOME=/root/.cache`
- `LEAN_SCRATCH=/workspace/.scratch/lean`

It also writes `.scratch/lean/env.sh`, which future shell sessions can source.

## Approved Install Policy

Lean tooling installs through `elan` are explicitly authorized for this project.
Future cloud images may persist `/root/.elan`, `/root/.cache`, and
`/workspace/.scratch` so Lean toolchains, Lake caches, and scratch files survive
between theorem-proving sessions.

## Project Policy

Do not create a full mathlib project until a first formal theorem target is
chosen. For now, start with standalone smoke files in `.scratch/lean/`.

A minimal Lean smoke file:

```lean
theorem lean_bootstrap_smoke : True := by
  trivial
```

Compile it with:

```bash
lean .scratch/lean/smoke.lean
```

Add `lean-toolchain`, `lakefile.lean`, and tracked `.lean` files only once the
formalization target and Lean version are chosen.
