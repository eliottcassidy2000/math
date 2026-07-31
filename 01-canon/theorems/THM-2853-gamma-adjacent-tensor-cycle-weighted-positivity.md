---
id: THM-2853
title: "Gamma adjacent-tensor cycle-weighted positivity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  real alpha>0, every mixed Gamma moment of two or more normalized
  adjacent differences is strictly positive.  After clearing natural
  rising-factorial denominators, the tensor is exactly the
  cycle-weighted enumerator of permutations in which each distinguished
  mark is neither a singleton nor preceded by a label in its private
  set.  A one-marked-cycle family gives a coefficientwise sharp lower
  bound.  Signed or complex recombinations inside interlaced cone spans
  remain outside the result.
source: root/gamma-adjacent-tensor-cycle-positivity-2026-07-28
depends_on: []
related:
  - THM-2841-all-order-adjacent-difference-factorial-tensor-positivity
  - THM-2846-arbitrary-positive-cone-moment-three-transverse-boundary
  - THM-2865-gamma-transverse-null-holotopy-and-uniform-fourth-exit
script: 04-computation/gmc_gamma_adjacent_tensor_cycle_weight_thm2853.py
output: 05-knowledge/results/gmc_gamma_adjacent_tensor_cycle_weight_thm2853.out
script_sha256: cc36e93b4cf68bc7b9ee23fa119ecf8cb8e85ced4441d471fee3315c93beba1a
output_sha256: 3193522f00588d3c840ec42ae862ae8b963464f3442a9117ab10d9de3cbe4478
hash_basis: LF-normalized bytes
---

# THM-2853 -- Gamma adjacent-tensor cycle-weighted positivity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

For a real parameter `alpha>0`, define

```text
L_alpha(s^m)=(alpha)_m,
f_n=s^n/(alpha)_n,
d_n=f_(n+1)-f_n,                                      (1)
```

where

```text
(alpha)_m=alpha(alpha+1)...(alpha+m-1).
```

For every `k>=2` and all `n_1,...,n_k>=0`,

```text
L_alpha(product_(i=1)^k d_(n_i))>0.                  (2)
```

The positivity is polynomial and coefficientwise before `alpha` is
specialized.  Put

```text
E_n=(alpha)_(n+1)d_n
   =s^(n+1)-(alpha+n)s^n,
M=sum_i(n_i+1),                                      (3)
```

and

```text
W_alpha(n_1,...,n_k)
 :=L_alpha(product_i E_(n_i)).
```

Then `W_alpha` is a polynomial in `alpha` with nonnegative integer
coefficients and

```text
W_alpha(n_1,...,n_k)
 >= (k-1)! alpha (alpha)_(sum_i n_i)                 (4)
```

coefficientwise.  Equality in `(4)` holds exactly when

```text
k in {2,3} and n_1=...=n_k=0.                        (5)
```

Consequently,

```text
L_alpha(product_i d_(n_i))
 >= (k-1)! alpha (alpha)_(sum_i n_i)
       /product_i (alpha)_(n_i+1)
 >0.                                                 (6)
```

## 1. Algebraic expansion

Expanding `(3)` and applying the Gamma moments gives

```text
W_alpha(n_1,...,n_k)
 =sum_(S subset [k]) (-1)^|S|
    product_(i in S)(alpha+n_i)
    (alpha)_(M-|S|).                                 (7)
```

The alternating expression `(7)` does not display its sign.  The
following cycle-weighted model removes the cancellation.

## 2. Marked-cycle model

Take `M` labelled elements.  Distinguish marks

```text
x_1,...,x_k
```

and assign to `x_i` a private set `O_i` of `n_i` ordinary labels.  The
private sets are pairwise disjoint and, together with the marks, exhaust
the `M` labels.  Weight a permutation `pi` by

```text
alpha^(number of cycles of pi).                       (8)
```

The weighted enumerator of all permutations on `m` labels is
`(alpha)_m`: after adjoining a label, it either starts a new cycle, with
weight multiplier `alpha`, or is inserted after one of the `m` existing
labels.

Call `x_i` bad if either

```text
x_i is a singleton cycle, or
the predecessor of x_i lies in O_i.                   (9)
```

Fix a set `S` of marks required to be bad and delete those marks.  Each
deleted `x_i` is reconstructed independently in exactly two types:

```text
as a singleton cycle:                    weight alpha;
after one of the n_i labels in O_i:      n_i choices. (10)
```

The insertion sites in `(10)` are disjoint because the `O_i` are
disjoint.  Therefore the weighted intersection count is

```text
product_(i in S)(alpha+n_i)(alpha)_(M-|S|).           (11)
```

Inclusion-exclusion in the bad events proves the literal identity

```text
W_alpha(n_1,...,n_k)
 =sum_(pi: no x_i is bad) alpha^(cycles(pi)).          (12)
```

This proves coefficientwise nonnegativity.

## 3. Strict lower family and equality

Put all `k` marks in one directed `k`-cycle and let the ordinary labels
carry an arbitrary disjoint permutation.  Each marked predecessor is
another mark, so every such permutation is admitted by `(12)`.  There
are `(k-1)!` marked cycles, while the cycle-weighted enumerator on the
ordinary labels is `(alpha)_(sum_i n_i)`.  This proves `(4)`.

It remains to identify when this family exhausts `(12)`.

- If `k=2` or `3` and there are no ordinary labels, an admitted
  permutation cannot contain a singleton mark; it is therefore one full
  marked cycle.  Equality holds.
- If an ordinary label exists, insert it into the marked cycle
  immediately before a mark whose private set does not contain it.  The
  resulting permutation is admitted but is not in the disjoint
  marked-cycle family counted in `(4)`.  The inequality is strict.
- If `k>=4`, even with no ordinary labels, a product of two marked cycles
  of lengths at least two is admitted and is not counted in `(4)`.
  Again the inequality is strict.

This proves `(5)`, hence `(6)` and `(2)`.

## 4. Positive-cone consequence

For finite nonzero positive Gamma-adjacent cone elements

```text
U_j=sum_n lambda_(j,n)d_n,       lambda_(j,n)>=0,
```

multilinearity and `(2)` give, for every `k>=2`,

```text
L_alpha(U_1...U_k)>0.                                (13)
```

In particular every nonzero `U` in this cone satisfies

```text
L_alpha(U)=0,
L_alpha(U^k)>0 for every k>=2.                        (14)
```

At `alpha=1`, this is the factorial functional and the sign theorem
overlaps THM-2841.  THM-2841's forbidden-board argument has the stronger
index-sensitive derangement floor there; `(4)` is the uniform floor that
survives for every real Gamma shape.

## 5. Scope and the holotopy boundary

The marked-cycle formula is an all-order atomic tensor theorem.  It does
not imply positivity after signed or complex linear recombination of the
`d_n`, nor does it orient the complex projective span of two interlaced
positive cones.  THM-2846 already shows at `alpha=1` that two positive
interlaced cones can contain a nonzero complex line whose first three
moments vanish.  The proof-candidate THM-2865 asks whether that transverse
cancellation persists as `alpha` moves; THM-2853 neither assumes nor uses
that claim.

Thus tensor positivity and plane cancellation live on different levels:
every coordinate tensor in the positive basis is positive, while complex
coefficients can cancel their multilinear sum.  No multiplier
observation, arbitrary-radial coefficient theorem, SFC conclusion, NC2
closure, or GMC2 proof follows.

The first-order boundary is exact:

```text
L_alpha(d_n)=0.                                      (15)
```

No assertion is made for `alpha<=0`.

## 6. Exact evidence and independent audit

The exact companion:

1. checks nonnegative coefficients and the coefficientwise lower bound
   on `360` deterministic multi-index profiles;
2. verifies the exact equality classification `(5)`;
3. compares `(7)` with literal cycle-weighted permutation enumeration on
   `136` profiles with at most eight labels;
4. checks `6,534` positive rational Gamma specializations; and
5. verifies `90` first-order zero controls.

Every truth-bearing gate uses an explicit exception, so normal and
optimized executions are identical.  An independent audit rederived the
deletion/reconstruction bijection, the strict lower family, and the
equality boundary, and replayed the exact companion against the stored
transcript.

**QED.**
