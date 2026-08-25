# Hostile audit: exact Sun 2-4-6-8 residue-class block census

## Verdict

**PASS for the displayed census parameters.**

For

\[
M=1{,}062{,}347=11\cdot13\cdot17\cdot19\cdot23,
\qquad R=459{,}490,
\qquad B=1{,}001{,}999{,}999{,}999{,}999,
\]

the current `04-computation/sun_2468_exact_residue_class_block_scan.py`
implements an exact finite census of the canonical representation number

\[
a(n)=\#\{(w,x,y,z)\in\mathbb Z_{\ge0}^4:
n={w+2\choose2}+{x+3\choose4}+{y+5\choose6}+{z+7\choose8}\}.
\]

The frozen nine-block replay was run from the current embedded source and is
bound to it by the printed SHA-256

```text
cedc808e7665d6ae98b687c8b0c152a9ebb3680afa4785ef6066e58f27fbba92.
```

It proves precisely

\[
\{n\in[1,B]:n\equiv R\pmod M,\ a(n)=0\}
=\{896{,}315{,}812{,}331{,}399\}.
\]

Thus the displayed integer is the **first and unique hole in this one residue
class through `B`**.  This is not a global least-counterexample proof, not a
census of the other `M-1` residue classes, and not a classification of all
counterexamples beyond `B`.

No in-scope correctness defect was found.  The general-purpose CLI does have
three out-of-scope hardening issues described below.

## Exact universe and canonical-domain audit

The engine enumerates the OEIS/canonical nonnegative-parameter convention:

- `c4`: \({x+3\choose4}\), `x >= 0`;
- `c6`: \({y+5\choose6}\), `y >= 0`;
- `c8`: \({z+7\choose8}\), `z >= 0`;
- triangular lane: with `i=w+1 >= 1`,
  \({w+2\choose2}=i(i+1)/2\).

There is exactly one zero atom in each higher-binomial role, and the vectors
are role-labelled, so equal numerical values in different roles neither merge
nor create spurious multiplicity.  The strict `value < high` cutoff is exact:
every triangular contribution is at least one, hence no higher triple with sum
at least `high` can represent a target at most `high`.

`choose_u64` uses the exact recurrence

\[
{n\choose j}={n\choose j-1}\frac{n-j+1}{j}.
\]

At the actual bound, all values and recurrence intermediates are far inside
`uint64_t`; summed higher atoms and discriminants are evaluated in `u128`.

## Triangular congruence lanes

For `T(i)=i(i+1)/2`, the scanner uses period `2M`.  This is always a valid
period because

\[
T(i+2M)-T(i)=M(2i+2M+1).
\]

For odd `M`, the minimal period may be `M`; using `2M` is harmless.  The root
classes in `[0,2M)` are disjoint, and stepping each by `2M` gives a bijective
partition of all positive indices satisfying the required congruence.  It does
not double-count indices.

For a fixed higher triple `s`, the target condition becomes

\[
T(i)\equiv R-s\pmod M,
\qquad
\max(1,\mathrm{start}-s)\le T(i)\le \mathrm{high}-s.
\]

The floating square root is only an initial estimate.  Both inverse-triangular
bounds are repaired by exact `u128` inequalities, so rounding mode and an
off-by-one square-root estimate cannot change the admitted index interval.
The subsequent residue-lane ceil/floor arithmetic is exact for the census
parameters.

## Address and interval boundaries

For each inclusive block `[start, high]`, the engine computes the first target
`base >= start` with `base == R (mod M)` and the exact target count

\[
1+\left\lfloor\frac{\mathrm{high}-\mathrm{base}}M\right\rfloor.
\]

Every mark is independently checked to be inside the interval, in the selected
class, and at address `(target-base)/M`.  The nine supplied blocks are adjacent
without overlap or gap and partition `[1,B]`:

| block | inclusive interval | first class target | class targets | zeros |
|---:|---|---:|---:|---:|
| 1 | `[1, 1,999,999,999,999]` | 459,490 | 1,882,624 | 0 |
| 2 | `[2,000,000,000,000, 11,999,999,999,999]` | 2,000,000,418,018 | 9,413,120 | 0 |
| 3 | `[12,000,000,000,000, 101,999,999,999,999]` | 12,000,000,210,658 | 84,718,082 | 0 |
| 4 | `[102,000,000,000,000, 201,999,999,999,999]` | 102,000,000,469,112 | 94,131,202 | 0 |
| 5 | `[202,000,000,000,000, 301,999,999,999,999]` | 202,000,000,520,206 | 94,131,202 | 0 |
| 6 | `[302,000,000,000,000, 401,999,999,999,999]` | 302,000,000,571,300 | 94,131,202 | 0 |
| 7 | `[402,000,000,000,000, 601,999,999,999,999]` | 402,000,000,622,394 | 188,262,404 | 0 |
| 8 | `[602,000,000,000,000, 801,999,999,999,999]` | 602,000,000,724,582 | 188,262,404 | 0 |
| 9 | `[802,000,000,000,000, 1,001,999,999,999,999]` | 802,000,000,826,770 | 188,262,404 | 1 |

The total is `943,194,644` selected targets.  The final selected target is
`1,001,999,999,866,611`.  The claimed hole has zero-based global class address
`843,712,847` and lies in block 9.

For blocks 7--9, an independent Python path reproduced the reported higher-box
sizes and feasible-triple counts:

| block | `#c4/#c6/#c8` | feasible triples |
|---:|---:|---:|
| 7 | `10963/868/262` | 2,220,090,999 |
| 8 | `11778/911/272` | 2,594,651,419 |
| 9 | `12452/945/279` | 2,927,893,995 |

These checks exercise the largest and most expensive enumeration boxes in the
manifest independently of the C++ triple-loop counters.

## Saturation and exact zero detection

Counts are stored in saturating `uint16_t` cells.  This does **not** compromise
zero detection:

- cells start at zero;
- each representation applies a monotone increment;
- saturation maps every count at least 65,535 to the positive sentinel 65,535;
- no operation ever maps a positive cell back to zero.

Therefore `count == 0` is exact, and every positive cell certifies at least one
canonical representation.  Counts below 65,535 are exact.

The generic option `--low-cutoff=65535` is misleading: a saturated cell whose
true count is larger still prints `65535`, so the predicate “true count <=
cutoff” is not exact at that endpoint. The census uses the zero predicate with
diagnostic cutoff zero, so this does not affect the singleton
hole result.  Likewise a reported minimum is exact whenever it is below
65,535; all supplied block minima are 0 or 4.

## Independent controls

`check_small_exact.py` compiles the embedded C++ source on a separate path and
compares full cell counts and metadata against direct Python enumeration.  It
covers odd and even moduli, target-free intervals, singleton intervals, both
interval endpoints, and the stated `[1,10000], M=33, R=20` control.  It also
checks the manifest partition, target addresses, last-three higher-box sizes
and feasible-triple counts, and actual-scope integer bounds.

Normal and optimized Python runs gave byte-identical output:

```text
small_exact_checks=92
manifest_arithmetic_checks=43
manifest_blocks=9
manifest_union=[1,1001999999999999]
manifest_class_targets=943194644
claimed_zero_vector=0,0,0,0,0,0,0,0,1
claimed_unique_zero=896315812331399
all_checks=PASS
```

Reproduce the independent audit from the repository root:

```powershell
python -B 04-computation/sun_2468_target_crt_class_independent_audit_thm4040.py
python -B -O 04-computation/sun_2468_target_crt_class_independent_audit_thm4040.py
```

## What the manifest proves, and what it does not

The implication from completed outputs to the class theorem is purely finite:
the blocks partition the full integer interval; within each block, target
addresses partition exactly the selected residue class; the enumeration is a
bijection with canonical representations; and exact zero cells across the
nine outputs are `0,0,0,0,0,0,0,0,1`, with the sole address decoding to the
claimed integer.

It does **not** bridge to other congruence classes.  In particular:

- “unique hole in `R mod M` through `B`” is justified;
- “smallest global hole” requires an independent no-hole census covering every
  integer below the candidate (or the other `M-1` classes);
- “all counterexamples” requires an unbounded theorem or a separately stated
  larger finite universe;
- an earlier global prefix check through `2e12`, if canonically certified, is a
  separate fact and does not fill the enormous off-class gap above that bound.

## General-CLI hostiles and recommended hardening

None of these fires for `M=1,062,347` and `B≈10^15`.

1. **Unchecked 64-bit offset addition.**  The accepted CLI range allows
   moduli close to `UINT64_MAX/2`, while expressions such as
   `root + triangular_period - first_index % triangular_period` and
   `first_index + offset` can overflow.  Restrict the documented modulus/input
   universe, or calculate and range-check these expressions in `u128`.
2. **Cutoff at the saturation sentinel.**  Reject `low_cutoff == 65535`, or
   label saturated output as `65535+` and promise exact low-count classification
   only for cutoffs below 65,535.
3. **Diagnostic-counter overflow.**  `feasible_higher_triples` and
   `representation_marks` are unchecked `uint64_t` diagnostics in arbitrary
   CLI scopes.  Add checked increments or document bounds.  The manifest totals
   are below seven billion and safe.

## Evidence packaging

The frozen output
`05-knowledge/results/sun_2468_target_crt_class_singleton_census_thm4040.out`
contains the embedded-source hash, all nine ordered block summaries, the exact
manifest cardinality, the singleton hole address, and the terminal `PASS`.
The wrapper compiles once with its displayed C++20 flags and returns success
only after all expected summaries and low-count rows match. The separate audit
compiles the same embedded source on a new temporary path with a different
optimization profile and checks direct Cartesian controls and manifest
arithmetic. A future second-compiler replay would be useful redundancy, but is
not a missing logical premise of the exact enumeration.
