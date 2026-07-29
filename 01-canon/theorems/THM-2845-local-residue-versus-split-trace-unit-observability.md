---
id: THM-2845
title: "Local residue versus split trace unit observability"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  A
  finite-dimensional associative algebra admits a linear scalar
  that is nonzero exactly on its units if and only if it is local, its
  residue division algebra is the base field, and the scalar is a nonzero
  multiple of the residue map.  On a product of local factors, every
  unit-safe scalar factors through the split residue algebra.  Outside
  F_2 it can use only one component; over F_2 it may use any odd number.
  Exact detection still requires one factor.  This unifies the modular
  p-group augmentation of THM-2839 with the split trace-zero-unit boundary
  of THM-2815/2842.
source: root/local-versus-split-unit-observability-2026-07-28
depends_on: []
related:
  - THM-2815-optimal-finite-laguerre-carrier-and-radial-selector-access-boundary
  - THM-2839-prime-power-unit-mass-full-spectrum-and-q11-response-provenance
  - THM-2842-ordered-positive-cone-vandermonde-multiplier-observability
script: 04-computation/algebra_local_split_unit_observability_thm2845.py
output: 05-knowledge/results/algebra_local_split_unit_observability_thm2845.out
script_sha256: 296abc5461741f96d491defaa744fe2fe78fbcb14f4018cbeba675889dc78d44
output_sha256: 807a59ef8329419fb61d012b0bdf6400f7276e4d3e2808622c19a30b8de4799e
hash_basis: LF-normalized bytes
---

# THM-2845 -- local residue versus split trace observability

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let `A` be a finite-dimensional unital associative algebra over a field
`k`, and let

```text
ell:A -> k                                               (1)
```

be `k`-linear.  Distinguish two properties:

```text
exact unit detector:
 ell(a)!=0 iff a is a unit;

unit-safe scalar:
 ell(u)!=0 for every unit u.                            (2)
```

No commutativity is assumed.

## 1. Exact scalar unit detectors are residue maps

The following are equivalent.

1. `ell` is an exact unit detector.
2. `A` is local, its residue division algebra is the base field,

   ```text
   A/J(A) isomorphic to k,                              (3)
   ```

   and

   ```text
   ell=c pi,                         c in k^*,          (4)
   ```

   where `pi:A->A/J(A)~=k` is the residue map.

Suppose `(1)` is exact.  Then

```text
ker ell={nonunits of A}.                                (5)
```

If `x` is a nonunit, then `ax` and `xa` are nonunits for every `a`.
Otherwise `x` would have a one-sided inverse; in a finite-dimensional
algebra, the corresponding injective multiplication map is surjective,
so `x` would be invertible.  Thus `ker ell` is a two-sided ideal.

Every element outside this ideal is a unit.  Hence the ideal is the
unique maximal ideal and equals `J(A)`.  Since `ker ell` has codimension
one,

```text
A/J(A)~=k,                                             (6)
```

and two nonzero linear maps to `k` with this kernel differ by a scalar,
giving `(4)`.  Conversely, in a local algebra an element is a unit
exactly when its residue is nonzero, proving the reverse implication.

## 2. Unit-safety is already exact on a local algebra

More generally, let `A` be local with residue division algebra

```text
K=A/J(A).                                              (7)
```

If `ell` is unit-safe and `j in J(A)`, then `ell(j)=0`.  Otherwise

```text
1-ell(1)ell(j)^(-1)j                                  (8)
```

would be a unit killed by `ell`.  Therefore `ell` descends to a
`k`-linear functional on `K`.  A functional nonzero on every element of
`K^*` has zero kernel and is injective, so

```text
[K:k]=1.                                               (9)
```

Thus on a local algebra, unit-safety, exact detection, and nonzero scalar
multiples of a `k`-valued residue map are the same property.

## 3. Split products and the sharp `F_2` exception

Now suppose

```text
A=A_1 x ... x A_D,
A_i/J(A_i)~=k.                                        (10)
```

Put `J=product_i J(A_i)`.  The argument `(8)` shows that every unit-safe
scalar kills `J`, so it has the unique form

```text
ell(a_1,...,a_D)=sum_(i=1)^D w_i pi_i(a_i).            (11)
```

Let

```text
r=#{i:w_i!=0}.                                         (12)
```

Then

```text
ell is unit-safe
 iff
   r=1,                         when k!=F_2,
   r is odd,                    when k=F_2.             (13)
```

Indeed, on a unit every residue coordinate is nonzero.  After rescaling
by the nonzero weights, cancellation is exactly the existence of

```text
y_1+...+y_r=0,                 y_i in k^*.             (14)
```

For `k!=F_2`, every `r>=2` has such a tuple: use opposite pairs, and for
odd `r` use one nonzero zero-sum triple plus pairs.  Over `F_2`, every
`y_i=1`, so cancellation occurs exactly when `r` is even.

Exact detection in `(10)` occurs if and only if `D=1` and `r=1`, by
Section 1.  In particular, the full trace on `F_2^3` is unit-safe but not
exact: it is nonzero on the sole unit `(1,1,1)`, yet also on the nonunit
`(1,0,0)`.

## 4. Exact finite-field census

For full nonzero weights over `F_q`, scaling coordinates reduces to the
number of nonzero solutions of one zero-sum equation:

```text
N_D(q)
 =#{(y_1,...,y_D) in (F_q^*)^D:sum_i y_i=0}
 =[(q-1)^D+(q-1)(-1)^D]/q.                            (15)
```

Additive-character orthogonality proves the second equality: the trivial
character contributes `(q-1)^D`, and each of the other `q-1` characters
contributes `(-1)^D`.

Thus a full weighted trace has a trace-zero unit for every `D>=2`,
except exactly when

```text
q=2 and D is odd.                                      (16)
```

If only `r` weights are nonzero, the number of killed units is

```text
(q-1)^(D-r) N_r(q).                                    (17)
```

Over any infinite field and `D>=2`, a full weighted trace always has a
trace-zero unit.

## 5. The two canonical poles

### 5.1. Modular p-group algebras are local

For every finite `p`-group `G`,

```text
F_p[G]/J(F_p[G])~=F_p,                                 (18)
```

and the residue map is augmentation.  Section 1 therefore recovers
THM-2839's exact modular statement:

```text
a in F_p[G] is a unit
 iff augmentation(a)!=0.                              (19)
```

The `p`-unit integral hypothesis in THM-2839 is the sidecar that makes
reduction land at this local pole.  It is load-bearing.  In
`Q[C_p]`, the group sum

```text
1+g+...+g^(p-1)                                       (20)
```

has nonzero augmentation `p` but vanishes in every nontrivial Fourier
component, so it is not a unit.

### 5.2. The Laguerre carrier is split

After complex scalar extension, THM-2815's carrier is

```text
A_D tensor_Q C ~= C^D,                                 (21)
```

and its positive readout is a full weighted trace.  For `D>=2`,
Sections 3--4 force trace-zero units; THM-2842 gives the rational explicit
example `ell_(D-1)`.  A Christoffel cardinal selector is a coordinate
projection, hence unit-safe, but it is an exact unit detector only when
`D=1`.

This local-versus-split distinction explains why a single augmentation
works in the modular group algebra while a positive Gauss--Laguerre
splitting cannot make scalar nullity nodewise.  Positivity of weights is
orthogonal to locality.

## 6. Minimal hostile controls

The boundaries already occur in the smallest algebras.

1. In `k[epsilon]/(epsilon^2)`, the residue coefficient is exact, while
   every radical-sensitive scalar kills some unit `1+c epsilon`.
2. On `C/R` or `F_4/F_2`, every base-field linear scalar kills a nonzero
   field element, hence a unit.
3. On `k^2`, the full weighted trace kills the unit `(w_2,-w_1)`.
4. On `F_2^2`, full trace kills the sole unit; on `F_2^3`, it is
   unit-safe but not exact.
5. The rational group sum `(20)` separates characteristic-zero split
   augmentation from THM-2839's modular local detector.

The theorem classifies scalar unit observability.  It does not claim that
unit detection alone supplies a physical Wick multiplier, an ancestry
phase, or a Gaussian moment.

## 7. Exact companion

The exact companion exhausts:

1. five truncated local algebras and all `73` scalar functionals, with
   the `10` unit-safe functionals exactly the residue maps;
2. `58/58` radical-sensitive killed-unit hostiles;
3. every `F_2`-linear scalar on `F_4`, with no unit-safe example;
4. all scalar functionals on `M_2(F_2)`, `M_2(F_3)`, and `M_2(F_5)`,
   with no unit-safe or exact detector;
5. `1,269` weighted split-product scalars, including the sharp `F_2`
   parity and exactness boundaries;
6. the formula `(15)` in `24` direct finite-field universes through
   `q=7,D=6`; and
7. the characteristic-zero group-sum hostile and modular delta controls
   for `p=2,3,5,7`.

All truth-bearing gates are explicit exceptions and all arithmetic is
exact.  Reproduce with

```text
python 04-computation/algebra_local_split_unit_observability_thm2845.py
python -O 04-computation/algebra_local_split_unit_observability_thm2845.py
```

Both modes byte-match the stored transcript.

## 8. Independent hostile audit

An independent audit checked the noncommutative one-sided-inverse step,
the exact-detector/local-residue iff, local unit-safety over a residue
division algebra, the split-product radical descent, the `F_2` parity
exception, and the finite-field character count.  It also replayed every
dual-number, `F_4/F_2`, matrix, split-product, and rational group-ring
hostile in normal and optimized modes against the stored transcript and
declared LF hashes.

**QED.**
