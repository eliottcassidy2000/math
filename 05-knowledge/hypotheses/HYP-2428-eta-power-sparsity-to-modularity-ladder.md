# HYP-2428 - Eta powers form a sparsity-to-modularity cancellation ladder

**Status:** CONFIRMED computationally for the pilot identities; OPEN as a general
proof lens.
**Source:** codex-2026-06-11-P2.
**Artifacts:** `04-computation/cancellation_gate_atlas_codex.py`,
`05-knowledge/results/cancellation_gate_atlas_codex.out`.

## Statement

The powers of Euler's eta product form a ladder of cancellation mechanisms:

```text
eta^1  : Euler pentagonal sparsity
eta^3  : Jacobi triangular sparsity
eta^24 : Ramanujan Delta, dense but modularly controlled
```

The lesson is that cancellation gates need not remain sparse. They can move from
literal lacunarity to modular control.

## Evidence

Through `q^120`:

```text
eta^1  support  18/121, max |coeff| = 1
eta^3  support  16/121, max |coeff| = 31
eta^8  support  68/121, max |coeff| = 11914
eta^24 support 121/121, max |coeff| = 225755128648
```

The script verifies Jacobi's identity

```text
prod(1-q^m)^3 = sum_{k>=0} (-1)^k (2k+1) q^{k(k+1)/2}
```

through `q^120`. For `eta^24 = Delta/q`, every coefficient through the window is
nonzero; sparsity is gone, but modularity remains.

## Connection

This extends HYP-2424 and THM-485: Euler's pentagonal signs are not the only
cancellation gate. They are the first rung. Higher powers show how a gate can
stop being sparse while preserving arithmetic control.

The Type II code analogy is direct: the length-72 Gleason gate is not sparse in
the full coefficient list, but it is controlled by an invariant ring that forces
the dangerous low weights to vanish.
