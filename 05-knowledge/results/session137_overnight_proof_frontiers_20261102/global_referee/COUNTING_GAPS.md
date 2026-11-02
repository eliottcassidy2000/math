# Four tournament determinant gaps are #P functions

**PROVED; independently executable exact checks pass.** This is a counting-complexity consequence of the audited vertex-cycle identities, not an assertion that arbitrary nonnegative GapP functions belong to #P. The main determinant proof, equality cases, and literature audit are in `REPORT.md`.

## 1. Domain and theorem

An input is a binary adjacency matrix of a finite labelled tournament. Let `H(T)` count directed Hamiltonian vertex words and put

```
u(T)=det(I+A),          v(T)=(-1)^|T| det(A),
P(T)=u(T)+|v(T)|,
C_all(T)=number of simple directed cycles,
C_proper(T)=number of such cycles not spanning T.
```

Cycles are identified under cyclic rotation; they have length at least three. Define

```
delta(T)=H(T)-u(T),
e(T)=H(T)-u(T)-C_all(T),
g(T)=H(T)-P(T),
Gamma(T)=H(T)-P(T)-C_proper(T).
```

**Theorem.** All four functions are in #P on nonempty tournament inputs. They can be made total #P functions on matrix encodings by assigning zero to malformed/non-tournament inputs and to the empty input. For `delta` and `e` the empty value is naturally zero under `H(empty)=u(empty)=v(empty)=1`. For `g` and `Gamma`, zero at empty is an explicit extension, **not** the value of the displayed subtraction formula (which would be -1).

The certificates below have length `O(n^2 log n)` and a deterministic polynomial-time verifier. The construction is an exact count of new determinant-labelled objects. It does not supply an efficient Hamiltonian-path counter or a direct injection into Hamiltonian words.

## 2. Positive terminal coefficients and exact recurrences

For every nonempty induced tournament `R`, algebraic skew-matrix positivity proves

```
u(R)+v(R)>0,     u(R)-v(R)>0.
```

These are integers. Consequently `u(R)>=1` and `u(R)±v(R)>=1`. At empty, `u=v=1`, so `u>=1` still holds and `u-epsilon v=1-epsilon` is 0 or 2 for `epsilon in {+1,-1}`. The candidate's independent determinant/odd-cycle identities yield the following exact recurrences; their positivity is not a guess based on computed values.

Fix `x` to be the least original vertex label currently present. Let `C` range over cycles through `x`, let `R_C=T-V(C)`, and call a cycle odd or even according to its number of vertices. Then

```
delta(T)=delta(T-x)
 +sum_(C odd through x) 2 delta(R_C)
 +sum_(all C through x) u(R_C).                         (D)
```

Every cycle either avoids `x` or contains it. Subtracting this partition of `C_all` gives

```
e(T)=e(T-x)
 +sum_(C odd through x) 2 delta(R_C)
 +sum_(all C through x) (u(R_C)-1).                     (E)
```

Both recurrences have base zero at empty. All terminal coefficients are nonnegative integers, so they themselves prove `delta,e>=0` inductively.

For a fixed nonempty input `T`, compute

```
epsilon=-1 if v(T)>0, and epsilon=+1 otherwise.
```

Thus `epsilon v(T)=-|v(T)|`, including the zero case. The signed-margin identity gives

```
g(T)=delta(T-x)
 +sum_(C odd through x) 2 delta(R_C)
 +sum_(all C through x) (u(R_C)-epsilon v(R_C)).         (G)
```

All cycles avoiding `x` are proper cycles of `T`; a cycle through `x` is proper exactly when `R_C` is nonempty. Hence

```
Gamma(T)=e(T-x)
 +sum_(C odd through x) 2 delta(R_C)
 +sum_(all C through x)
     (u(R_C)-epsilon v(R_C)-1_[R_C nonempty]).           (Q)
```

The terminal coefficients in `(G)` and `(Q)` are nonnegative integers by the positivity facts above. At empty complement they are exactly `1-epsilon`, so they are explicitly 0 or 2, with no hidden use of a false empty strong bound. The sign `epsilon` is computed once from the original input and used at the one top-level `(G)` or `(Q)` step. All continuations use `(D)` or `(E)`, not another signed-margin step.

## 3. Abstract accepting traces

There are four mode codes: `delta=0`, `e=1`, `g=2`, `Gamma=3`. The initial mode is part of the verifier's specification/input parameter, not a nondeterministic choice. Given current vertex set `U`, choose `x=min U`.

A cycle through `x` is encoded by its old path `P=(p1,...,pell)` after removing `x`. Require `ell>=2`, distinct vertices of `U\{x}`, and all arcs in `x,p1,...,pell,x`. Because `x` is fixed, this encodes every cycle through it once, with no rotation factor. The cycle is odd exactly when `ell` is even.

The verifier accepts traces made from three actions:

1. **Skip:** remove `x`. A delta or e mode continues in the same mode. A g mode continues in delta mode. A Gamma mode continues in e mode.
2. **Colored recurse:** choose an odd cycle through `x`, choose a color `b in {0,1}`, delete its vertices, and continue in delta mode. The two colors are the factor 2 in the recurrences.
3. **Leaf:** choose any cycle through `x`, delete its vertices to obtain `R`, and supply an integer label `a` in the half-open interval `[0,B_mode(R))`, where

   ```
   B_delta(R)=u(R),
   B_e(R)=u(R)-1,
   B_g(R)=u(R)-epsilon v(R),
   B_Gamma(R)=u(R)-epsilon v(R)-1_[R nonempty].
   ```

   This action terminates and accepts.

At empty current vertex set, the verifier rejects; there is no empty accepting trace in delta or e mode. A nonterminating skip/recursion sequence is not accepted. A zero leaf bound has no admissible label. Every recursive action deletes at least one vertex, so a trace uses at most `n` records.

For delta and e modes, induction on `|U|` identifies the accepting-trace counts with `(D)` and `(E)`. For g and Gamma, the first action supplies exactly `(G)` and `(Q)` and then invokes the already-established recursive modes. This proves the counting identities at the abstract trace level.

## 4. Canonical fixed-length binary certificate

Counting complexity requires exactly one binary witness for every accepting trace, not an ambiguous variable-length description. Let the original labels be `[n]`, and set

```
L=ceil(log2(n+1)),
B=ceil(log2(2 n!+1)),
record_bits=5+(n+1)L+B.
```

Use a binary witness with **exactly n records**, each of this fixed width. Each record has the following fields in order:

| Field | Width | Meaning |
|---|---:|---|
| mode | 2 | Expected current mode, codes 0 through 3 |
| action | 2 | skip=0, recurse=1, leaf=2; code 3 rejected |
| old-path length | L | The number `ell` of old path vertices |
| old-path vertex slots | nL | Exactly n slots of width L; labels in `[n]`, zero padding |
| color | 1 | 0 or 1 only on colored recursion; zero otherwise |
| leaf label | B | Fixed-width unsigned integer; zero unless action is leaf |

The verifier checks the exact bit length and rejects any nonbinary symbol. For every active record it checks that the mode is the expected one. Further canonical rules are:

* Skip has zero path length, zero vertex slots, zero color, and zero label.
* Recurse and leaf use precisely the first `ell` vertex slots; all later slots are zero. Distinctness and membership tests ensure that active slots cannot be zero and cannot repeat or contain the root.
* Recurse requires an odd cycle and label zero. Its color is meaningful and can take either value.
* Leaf requires color zero and a label strictly smaller than its bound.
* After a leaf, **all remaining records are all-zero bit strings**. A nonzero suffix is rejected.
* A sequence reaching empty before a leaf, or exhausting all `n` records without a leaf, is rejected.

An all-zero record may represent a delta-mode skip before the terminal leaf, but it is unambiguously padding after the leaf: the terminal action determines where the active trace stops. There is no option to append extra active actions after it. Fixed widths also remove leading-zero ambiguity in integers. These rules give a bijection between accepted bit strings and the abstract accepting traces, including their color choices and leaf labels.

## 5. Polynomial bounds and verification time

For a `k`-vertex induced tournament, both `A` and `I+A` have entries in `{0,1}`. The Leibniz formula therefore bounds each absolute determinant by `k!<=n!`. Each nonnegative leaf bound is at most `2n!`, so the specified `B` bits suffice. The entire witness length is

```
n[5+(n+1)ceil(log2(n+1))+ceil(log2(2n!+1))]
 =O(n^2 log n).
```

Tournament validation, path membership/distinctness tests, cycle-arc checks, and state updates are polynomial-time. Integer determinants and their signs are computable in polynomial time with polynomial-bit intermediates, for example using Bareiss elimination. In the executable verifier only the original sign and the reached terminal coefficient require determinants; no enumeration of other branches is performed. Thus the verifier is polynomial-time in the adjacency-matrix input size and certificate length.

A nondeterministic machine guesses exactly this fixed witness length and applies the verifier. Its number of accepting branches equals each stated gap. This establishes membership in #P without invoking a general closure under subtraction.

## 6. Executable independent verification

`counting_gaps.py` is standalone Python using only the standard library. Its `GapSystem.verify(mode,witness)` is the polynomial-time verifier. Its separate recursive `traces`/`certificates` methods enumerate accepting traces for small diagnostic inputs; **the enumerator is not claimed to be polynomial-time**. Its independent expected counts use literal Hamiltonian permutations, separately minimum-rooted simple-cycle enumeration, and integer Bareiss determinants.

Completed checks:

* All **1,100 labelled tournaments of orders 0 through 5**.
* Exact matches for all four functions, with **9,882 accepting certificates** across that census.
* Every enumerated certificate accepted; no duplicate binary certificates.
* The six-vertex ordered sum `C3=>C3` independently gives `(delta,e,g,Gamma)=(5,3,4,2)` and exactly 14 accepted certificates, including nonzero Gamma and the two-color recursive branches.
* The seven-vertex `source=>C3=>C3` gives `(delta,e,g,Gamma)=(5,3,5,3)` and 16 certificates; this exercises the Gamma skip-to-e transition with a positive e remainder.
* **48 rejection controls**: wrong active mode, wrong initial mode, repeated path vertex, out-of-range leaf label, nonzero inactive color, nonzero path padding, nonzero suffix after a leaf, reserved action tag, short/long/nonbinary witness, malformed tournament matrices, and empty input.
* Normal and optimized Python executions produce byte-identical output. Test checks are explicit `require` calls and are **not removed** under `python -O`.

Reproduce:

```bash
python3 -B counting_gaps.py > counting_gaps.out
python3 -B -O counting_gaps.py > counting_gaps.optimized.out
cmp counting_gaps.out counting_gaps.optimized.out
```

The transcript prints the script's SHA-256. The universal proof is Sections 2–5; small certificates are independent diagnostics of its encoding and implementation.
