# Threshold42 extension: a longer prescribed square-sum prefix, sooner

**Status: PROVED, explicit construction; independent audit accompanies this packet.**
This is a separate extension of the immutable
[`session124_recursive_square_prefix_20261031`](../session124_recursive_square_prefix_20261031/README.md)
core, not a replacement of its certificates. The ordinary square-sum Hamiltonian
path/cycle existence results remain prior work; this result concerns one specified
nested ray, not arbitrary prescribed paths.

## Theorem

Put \(P_0=(1,3)\), and, writing \(p_j\) in one-based order, set
\[
P_{l+1}=(1,3,6,19)\mathbin\Vert
 (49p_j+(-1)^j4)_{j=1}^{4l+2},\qquad
S_l=\frac{85\cdot49^l-1}{2}.
\]
For every integer \(l\ge1\) and **every** integer \(N\ge S_l\), the square-sum
graph on \([1,N]\) has a Hamiltonian cycle whose oriented traversal starts with
**the original unscaled \(P_l\)**. The bound improves the first checkpoint's
\((107\cdot49^l-1)/2\). The first thresholds are
\(S_0=42,S_1=2082,S_2=102042,S_3=5000082\).

Consequently every finite prefix of this one explicit infinite square-sum ray
occurs in such a cycle for all sufficiently large interval sizes.
The formulas for the ray and its square roots are unchanged from the linked core.
No claim of arbitrary prescribed-path extension or optimal cycle threshold is made.

## Proof and the changed induction boundary

A weak B pair consists of square-sum Hamiltonian paths on \([1,n]\) and
\([1,n+1]\), each beginning \((1,3)\), ending at 8 for an even-length member
and 10 for an odd-length member, with the same positional parity for each common
vertex. `SMALL_BASES.json` explicitly supplies all eleven such pairs for
\(42\le n\le52\). The existing frozen base data supply \(53\le n\le2081\).
The checker verifies each interval, every square edge, endpoint, prefix, and
common-vertex positional parity exactly.

Use these bases **only through 2081**, then recurse for **every \(n\ge2082\)**:
\[
m=\lfloor(n-24)/49\rfloor,\quad d=n-49m\in[24,72].
\]
Now \(42\le m<n\). Apply the unchanged paired certificate indexed by
\((m\bmod2,d)\). Its affine interval-ownership and square-edge identities apply
also at \(m=42\): blocks and tails are nonempty, endpoints distinct, and no
extra base-specific edge is used. The first checkpoint's complete symbolic
proof therefore still maps a weak B pair to a weak B pair. Strong induction
constructs B(n) for every \(n\ge42\), using no old transition tables.

Every paired word starts with the singleton sequence \((1,3,6,19)\) followed
by intact, forward block 4. If \(n\ge S_l\), then
\[
\lfloor(n-24)/49\rfloor\ge S_{l-1},
\quad S_l=49S_{l-1}+24.
\]
Induction on \(l\) therefore gives prefix \(P_l\) in both members.
The same quotient inequality and the unchanged terminal-cycle certificates
produce a cycle beginning with \(P_l\) for every \(N\ge S_l\).
Its last vertex is 8 if N is even and 15 if N is odd, closing to 1 via 9 or 16.
Terminal outputs need no paired-parity requirement because they are not reused.

**Important boundary:** retaining the old cutoff 2620 would invalidate the
new prefix induction for sizes 2082–2620. `core42.py` uses cutoff **2081**.
This is tested directly at 2081,2082,2083,2620,2621 and deeper boundaries.

The first checkpoint's original wider corollaries (2303 paths and depth-width
amplification) are not part of this minimal extension's proof or dependency set.
No claim is made that 42 is the smallest possible seed for a different invariant.

## Portable reproduction

Place this directory beside the frozen core directory and the existing
`square_sum_B_base_interval_20261019` base directory. Python standard library only:

```sh
python verify_extension.py
python -O verify_extension.py
python construct_cycle.py 2082
```

For another layout, pass both `--core-dir /path/to/frozen/core` and
`--base-dir /path/to/old/base/data`. `DEPENDENCIES.json` pins every immutable
first-checkpoint file and the eleven small bases. Existing large base data and
196 transition certificates are referenced, not duplicated. The new constructor
uses only old base pairs through 2081, the eleven new small pairs, and the frozen
new paired/direct tables. It uses no search solver, network, Git database,
old A/B/final transition table, or third-party C source.

Controls include full streamed cycles at 2082,2083,2620,2621,102042,102043,
and 5000082. Every vertex, edge, closing edge and prescribed prefix is checked;
no giant cycle file is needed. SHA256 hashes cover the decimal-lines stream.

## Scope of the other overnight target

Tournament Hamiltonian-path count infinite completeness remains **OPEN**:
this square-sum construction does not prove that all positive odd raw directed
Hamiltonian-path counts other than 7 and 21 are attained.
