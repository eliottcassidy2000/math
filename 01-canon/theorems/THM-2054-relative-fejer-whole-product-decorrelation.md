---
id: THM-2054
title: Relative Fejer whole-product decorrelation along a character line
status: >
  PROVED. If every scalar resonance of height at most H along an integer line
  lifts to the corresponding vector-character relation, then the line average
  of a product and its lifted-character torus average differ by at most twice
  the sum of the one-factor Fejer L1 errors. The statement extends to finite
  signed atom expansions. For two-scale characters (b_i,c_i) evaluated on
  (1,M), the lift condition is automatic when |M|>H sum_i|b_i|. For the exact
  seven-sector inclusion-exclusion algebra, after dropping the pinned zero
  factor, the complete H=2^20 error is <90816/1048577; rowwise H=2^19 is below
  every recorded pinned-base cap-Q margin. This last comparison is only a
  numerical budget compatibility result. The theorem proves the abstract
  nonresonant filter and its coefficient budget, but not a model-specific LRC
  plateau, bounded-resonance classification, or LRC(14).
source: codex-2026-07-21-LRC-unrelated-transfer
depends_on: []
related:
  - THM-546
  - THM-548
  - THM-699-total-variation-far-element-contraction
  - THM-2051
  - THM-2052
  - THM-2053
  - HYP-2684
script: 04-computation/lrc14_relative_fejer_atom_budget_referee_codex_20260721.py
output: 05-knowledge/results/lrc14_relative_fejer_atom_budget_referee_codex_20260721.out
script_sha256: 510d455298aa5f89c54c8c61c295b3a3f032f3830bacbad7e5a36f765f96c711
output_sha256: f9621dd60916ff64893ec2a9005b4f48687102faeb2644b4c5d0c0244edbd79d
---

# THM-2054 -- relative Fejer whole-product decorrelation

Let `T=R/Z` carry Haar probability measure. Fix characters

```text
chi_i in Z^d       (1<=i<=m)
```

and an integer direction `lambda in Z^d`. Put

```text
a_i=chi_i dot lambda,
```

and assume every `a_i` is nonzero. For bounded measurable functions
`f_i:T->C`, let

```text
p_i=F_H*f_i,
B_i=||f_i||_infinity,
epsilon_i=||f_i-p_i||_1,                               (1)
```

where, to fix the order convention used throughout this theorem,

```text
F_H(x)=sum_(|n|<=H) (1-|n|/(H+1)) exp(2 pi i n x).     (1a)
```

Thus `F_H` is the normalized positive Fejer kernel with Fourier support
`[-H,H]`; below, `N=H+1`.

## 1. The bounded-resonance lift theorem

Assume the height-`H` resonance-lift condition

```text
for every (n_1,...,n_m) in Z^m with |n_i|<=H,

sum_i n_i(chi_i dot lambda)=0
        implies
sum_i n_i chi_i=0 in Z^d.                              (2)
```

Zero frequency coordinates are allowed in (2); only the all-zero tuple is
vacuous. This point is essential for noncentered coverage atoms.

Define the actual line average and the lifted-character torus average by

```text
I_lambda(f)=integral_T product_i f_i(a_i t) dt,
J_chi(f)=integral_(T^d) product_i f_i(chi_i dot x) dx.  (3)
```

The latter need not be a product of independent one-factor averages: the
vector relations among the `chi_i` are deliberately retained.

**Theorem.** Under (2),

```text
|I_lambda(f)-J_chi(f)|
 <=2 sum_i epsilon_i product_(j!=i) B_j.               (4)
```

### Proof

Fejer convolution is a positive contraction, so

```text
||p_i||_infinity<=B_i,
```

and `p_i` has Fourier support in `[-H,H]`. Since multiplication by the
nonzero integer `a_i` preserves Haar measure on `T`, product telescoping gives

```text
|I_lambda(f)-I_lambda(p)|
 <=sum_i epsilon_i product_(j!=i) B_j.                 (5)
```

The character `x->chi_i dot x` also pushes Haar measure on `T^d` to Haar
measure on `T` (it is nonzero because `a_i!=0`). The same telescope gives

```text
|J_chi(f)-J_chi(p)|
 <=sum_i epsilon_i product_(j!=i) B_j.                 (6)
```

It remains to compare the two finite Fourier polynomials. Expanding them,

```text
I_lambda(p)
 =sum_(|n_i|<=H, sum_i n_i a_i=0) product_i p_i_hat(n_i),

J_chi(p)
 =sum_(|n_i|<=H, sum_i n_i chi_i=0) product_i p_i_hat(n_i).  (7)
```

Every vector relation in the second sum gives a scalar relation after dotting
with `lambda`; condition (2) supplies the converse. Thus the two index sets in
(7) agree exactly and

```text
I_lambda(p)=J_chi(p).
```

Combining this equality with (5)--(6) proves (4). QED.

## 2. Finite signed atom expansions

The result is stable under the Boolean/inclusion--exclusion expansions used by
coverage functions. Suppose

```text
Phi(x)=sum_alpha c_alpha product_i f_(alpha,i)(chi_i dot x),
```

with finitely many atoms. For every retained nonconstant factor, require its
character to be nonzero along `lambda` and require (2) for that atom; one global
condition (2), with zero frequencies inserted for omitted factors, suffices.
Omitted factors are simply dropped. A genuinely constant factor can be pulled
out exactly regardless of its character. A nonconstant factor evaluated on a
zero character is instead the point value `f(0)` and must be factored out
before applying the theorem, because an `L1` Fejer error does not control point
evaluation. Put

```text
B_(alpha,i)=||f_(alpha,i)||_infinity,
epsilon_(alpha,i)=||f_(alpha,i)-F_H*f_(alpha,i)||_1.
```

Under these inherited hypotheses, applying (4) atom by atom gives

```text
|integral_T Phi(lambda t)dt-integral_(T^d)Phi(x)dx|
 <=2 sum_alpha |c_alpha| sum_i epsilon_(alpha,i)
                         product_(j!=i)B_(alpha,j).     (8)
```

This is deliberately a whole-product estimate. No absolute sum over the
infinite relation lattice appears.

## 3. Explicit BV cost

If `f` has circular total variation `Var(f)`, then with `N=H+1`,

```text
||f-F_H*f||_1
 <=Var(f)(1+2 log N)/(4N).                             (9)
```

Indeed,

```text
||f-f(.-y)||_1<=Var(f)||y||
```

and, on `[-1/2,1/2]`,

```text
F_H(y)<=min(N,1/(4N|y|^2)).
```

Splitting at `|y|=1/(2N)` gives

```text
integral_T ||y||F_H(y)dy<=(1+2 log N)/(4N),
```

which proves (9). For interval indicators `Var(f)=2`; this recovers the
constant used in THM-2051.

## 4. Two-scale separation

Let

```text
e_i=b_i+M c_i,
chi_i=(b_i,c_i),       lambda=(1,M),                  (10)
```

with integer `b_i,c_i,M`, and assume `e_i!=0`. Then (2) holds whenever

```text
|M|>H sum_i |b_i|.                                    (11)
```

To see this, write

```text
B=sum_i n_i b_i,       C=sum_i n_i c_i.
```

A scalar resonance says `B+MC=0`. If `C=0`, then also `B=0`, so the vector
relation holds. If `C!=0`, integrality gives

```text
|M|<=|MC|=|B|<=H sum_i|b_i|,
```

contradicting (11). Thus every bounded scalar resonance lifts to the
two-dimensional cluster relation, and (4) or (8) applies.

The strict cutoff in this uniform sufficient criterion cannot be weakened to
`>=`. At `H=M=1`, take

```text
chi_1=(1,0),       chi_2=(0,1),       lambda=(1,1).
```

Then `(n_1,n_2)=(1,-1)` is a scalar resonance along `lambda` but not a vector
relation. Allowing zero coordinates in the resonance tuple is equally
necessary: adjoining any noncentered third factor and setting its frequency
to zero preserves this alias.

The same statement in `d` scales is already contained in Section 1: choose
the character vectors `chi_i` to retain every cluster coefficient and require
that the chosen scale direction `lambda` introduce no new height-`H` scalar
aliases.

## 5. Exact seven-sector atom budget

The signed coefficient budget can be made explicit for the seven-sector
inclusion--exclusion identity of THM-534. Let

```text
E={0,e_1,...,e_r},       7<=r<=11,
```

and for `A subset {1,...,6}` let `f_A` be the indicator of the complement of
the union of the sectors indexed by `A`. Then

```text
p_0(E)=sum_(A subset {1,...,6})(-1)^|A|
             integral_T product_(i=1)^r f_A(e_i t)dt.  (12)
```

The pinned zero offset has `f_A(0)=1`, since `A` contains only inner sectors,
and has been removed from the product. This removal is load-bearing: point
evaluation at a zero character is not controlled by an `L1` Fejer error.

The circular variation satisfies

```text
Var(f_A)<=2|A|,
sum_(A subset {1,...,6}) |A|=192.                      (13)
```

Suppose each nonzero offset has a lifted character `chi_i in Z^d` and
`e_i=chi_i dot lambda`, with the height-`H` resonance lift (2). Apply (8),
(9), and (13). If `p_0^lift` denotes the same signed atom expansion averaged
over the full lifted torus, then, for `N=H+1`,

```text
|p_0(E)-p_0^lift|
 <=384 r epsilon_H,
epsilon_H=(1+2 log N)/(2N).                            (14)
```

At the common relation cutoff `H=2^20` from THM-2051, `r<=11` and
`log(2^20+1)<21` give the fully explicit estimate

```text
|p_0(E)-p_0^lift|<90816/1048577.                       (15)
```

This is smaller than the smallest number recorded in HYP-2675's proposed
`cap-Q` comparison,

```text
129643/980980-90816/1048577
 =46851988331/1028633065460>0.                         (16)
```

Equation (16) is an error-budget comparison, not a plateau theorem. MISTAKE-
080/082 prohibit importing `Q(k-1)` solely from cardinality or calling a
decorrelated limit a finite-scale majorant. An LRC application must prove that
its particular full-torus model `p_0^lift` is the pinned-base model to which
the claimed plateau bound applies, or establish the correct model-specific
bound directly.

The scale cutoff can be halved if the five recorded margins are used row by
row. At `H=2^19`, `log(2^19+1)<20`, and (14) gives

```text
k             8          9          10         11         12
error <   18368      20992       23616      26240      28864
          -----      -----       -----      -----      -----
          174763     174763      174763     174763     174763
```

Each entry is below its corresponding HYP-2675 `cap-Q` number. The tightest
remaining numerical margin is the `k=9` value

```text
2064067449/171439007740>0.
```

Again, this conclusion is conditional on correct model identification; it is
not a repair by numerical slack alone.

For two-scale offsets `e_i=b_i+Mc_i`, the lift used in (14) is automatic under

```text
|M|>H sum_i |b_i|                                      (17)
```

by Section 4. Thus the factorwise analytic and signed-coefficient costs are
closed; model identification, the decorrelated plateau, and the bounded-
resonance branch remain.

## 6. Exact scope for the LRC wide branch

HYP-2684 proposed controlling

```text
integral_T H(x,Mx)dx-integral_(T^2)H(x,phi)dx dphi
```

through a mixed-variation bound on the full Boolean function `H`. For the
actual finite LRC coverage algebra, Section 2 supplies a different rigorous
route: expand `H` into finitely many products of one-dimensional interval
predicates, smooth each factor, use (11) to identify the finite constant terms
exactly, and pay only the two physical-space telescopes (8). The coefficient
sum and the interval variations give an explicit, possibly large, budget.

This proves the abstract **nonresonant** half of HYP-2684 without assuming a
mixed-variation estimate, and Section 5 pays the precise inclusion--exclusion
coefficient budget. It does not by itself prove the required `p0<=cap_k`
inequality: one must still identify and bound the correct full-torus plateau,
without the pinned-zero error of MISTAKE-080, and route bounded resonances
through HYP-2682/HYP-2676 or the newer THM-2052/2053 finite atlas. The theorem
therefore advances the wide-branch glue but does not prove LRC(14).

## 7. Assumption challenge and Tournament Analysis

The faithful combinatorial object is the height-filtered signed resonance
hypergraph on the character factors. A binary tournament orientation loses
the integer coefficient tuple, including which coordinates are zero, and
there is no canonical antisymmetric switch preserving (2). Tournament scores,
cycles, SCCs, and Hamiltonian paths are therefore not applicable to the proof.

The challenged assumption is that decorrelation requires pointwise Fourier
decay of the already-composed multivariable Boolean function. It is enough to
regularize its one-dimensional factors, preserve the whole finite product,
and separate scalar from vector resonances before taking absolute values.

## 8. Exact referee

The stored referee is an arithmetic regression, not an independent proof of
the Haar-pushforward, BV, or lifted-model arguments. It checks
`sum_A|A|=192`, every rational comparison in Section 5, the rowwise `H=2^19`
minimum, and the strict-cutoff alias. It passes both

```text
python 04-computation/lrc14_relative_fejer_atom_budget_referee_codex_20260721.py
python -O 04-computation/lrc14_relative_fejer_atom_budget_referee_codex_20260721.py
```

and byte-matches the stored output. The frozen SHA-256 hashes are recorded in
the frontmatter.
