---
id: THM-2348
title: "Prime-type rectangularity, target-token conditioning, and two-prime owner factorization"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For a
  prime-token allocation energy, robust factorization by oriented prime
  type is equivalent to vanishing every ANOVA component meeting two
  types, to additive type separation, to vanishing all conditional
  rectangle differences on the equal-prime quotient, and to Cartesian
  tied optima after every typewise perturbation. Fixing labelled target
  tokens contracts ANOVA components and Boolean Mobius coefficients by
  explicit subset sums; uniform deletion preserves the surviving ANOVA
  tensors. For two distinct prime targets, an exact constant-inside /
  strict-outside criterion characterizes product owner sets and gives the
  composite lift formula with one cohabitation correction theta. The
  correction has a universal Gordian triangle bound and equals the
  negative THM-2176 interaction defect at the unknot. A fixed Cartesian
  optimum alone does not imply decoupling, individual two-block Mobius
  coefficients are gauge-dependent, and no unconditional connected-sum
  additivity or owner-block deletion theorem is asserted.
source: codex-2026-07-25-prime-type-rectangularity
depends_on:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2339-prime-owner-deletion-and-target-allocation-hypergraph
  - THM-2346-global-allocation-anova-normal-form
related:
  - THM-2248-higher-interaction-defect-complex-and-tropical-trace-spectrum
  - THM-2336-prime-target-gordian-owner-diagram-and-bypass-split
script: 04-computation/prime_type_rectangularity_conditioning_thm2348.py
output: 05-knowledge/results/prime_type_rectangularity_conditioning_thm2348.out
script_sha256: d9a2b9a244497bf8c90173e67b4eae0dcd2b22eeb304fd1860752fcb9dac9b80
output_sha256: 3a6c6614519eb20450688f835e1931f3b4c467c6104512c17e7062ea03f7b281
hash_basis: working-tree bytes (LF)
---

# THM-2348 -- robust prime-type factorization and conditioning

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The independent audit checked the ANOVA grouping, the conditional
rectangle telescope for arbitrarily many prime types, and the exact
opposite-corner perturbation proving the difficult converse. It also
checked the labelled conditioning law, the two-colour score-table gauge,
the constant-inside / strict-outside owner criterion, and every metric
inequality below. The exact companion replays identically under ordinary
and optimized Python and matches the stored transcript after LF
normalization.

THM-2346 gives a unique interaction tensor for every subset of labelled
target tokens. This theorem identifies exactly when those tensors permit
the target allocation problem to split by **oriented prime type**, rather
than merely appearing to split at one accidental optimum. It also gives
the functorial law for fixing target tokens and the sharp two-prime
Gordian specialization.

## 1. Prime-type quotient and mixed ANOVA tensors

Let the current packet partition have block set

```text
C=pi,                         m=|C|>=2,
```

and let `S` be the labelled prime-occurrence set of the target. Partition
it by oriented prime type:

```text
S=disjoint_union_(P in Pcal) S_P.                  (1)
```

The allocation energy

```text
E:C^S->R
```

is invariant under

```text
Gamma=product_P Sym(S_P).                           (2)
```

It therefore descends to

```text
E_bar:Y->R,

Y=C^S/Gamma
  isomorphic to product_P Y_P,                     (3)
```

where `Y_P=C^(S_P)/Sym(S_P)` is the finite count-table space for copies
of type `P`.

Write THM-2346's unique uniform ANOVA expansion as

```text
E=sum_(T subset S) H_T.                            (4)
```

Call `T` **mixed-type** when it meets at least two distinct sets `S_P`.

> **Robust prime-type rectangularity theorem.** The following are
> equivalent.
>
> **(A) Mixed ANOVA vanishing.**
>
> ```text
> H_T=0 for every mixed-type T.                     (5)
> ```
>
> **(B) Additive type separation.** There are unique mean-zero,
> type-symmetric functions `E_P:C^(S_P)->R` such that
>
> ```text
> E=H_emptyset+sum_P E_P,
>
> E_P=sum_(emptyset!=T subset S_P)H_T.              (6)
> ```
>
> **(C) Conditional rectangle vanishing.** For all distinct `P,Q`, all
> `p_0,p_1 in Y_P`, `q_0,q_1 in Y_Q`, and every fixed configuration `z`
> of the remaining type coordinates,
>
> ```text
> E_bar(p_0,q_0,z)+E_bar(p_1,q_1,z)
> -E_bar(p_0,q_1,z)-E_bar(p_1,q_0,z)=0.             (7)
> ```
>
> **(D) Perturbation-robust tied factorization.** For every family of
> arbitrary functions `phi_P:Y_P->R`, putting
>
> ```text
> Phi(y)=sum_P phi_P(y_P),
> ```
>
> gives
>
> ```text
> argmin_Y(E_bar+Phi)
>  =product_P pr_P(argmin_Y(E_bar+Phi)).            (8)
> ```

Here “arbitrary” in (D) means arbitrary functions on the whole finite
type-configuration spaces, not merely tokenwise unary fields. Rational
perturbations suffice when `E` is integer-valued. The product in (8) is
the product of its nonempty projections.

### Proof

The orthogonal ANOVA summands in (4) split into the constant, the summands
supported inside one `S_P`, and the mixed-type summands. Uniqueness gives
(A) iff (B), including the formula and uniqueness in (6).
Type-symmetry follows because `Gamma` permutes the ANOVA summands inside
each displayed sum.

Condition (B) makes every rectangle difference (7) vanish. Conversely,
assume (C), fix a base point `y^0 in Y`, and replace its type coordinates
one at a time. Equation (7) says that the increment caused by replacing
the `P` coordinate is independent of all other coordinates. Hence

```text
E_bar(y)
 =E_bar(y^0)
  +sum_P [E_bar(y_P,y^0_(-P))-E_bar(y^0)],          (9)
```

which is (B) on the quotient and therefore on the labelled space.

Condition (B) also makes (8) immediate. For the last converse, suppose a
rectangle in (7) has nonzero value

```text
Delta=e_00+e_11-e_01-e_10.                         (10)
```

Use large typewise penalties to freeze every other type at `z` and to
exclude all `P,Q` states except the selected two. Row and column
potentials do not change `Delta`. If `Delta<0`, put
`epsilon=-Delta/2` and use

```text
phi_P(p_0)=0,
phi_Q(q_0)=-e_00,
phi_Q(q_1)=epsilon-e_01,
phi_P(p_1)=epsilon-e_10+e_00.                      (11)
```

The surviving table becomes

```text
[ 0        epsilon ]
[ epsilon  0       ].                              (12)
```

If `Delta>0`, put `epsilon=Delta/2` and use

```text
phi_P(p_0)=0,
phi_Q(q_1)=-e_01,
phi_Q(q_0)=epsilon-e_00,
phi_P(p_1)=-e_10-epsilon+e_00.                     (13)
```

The table becomes

```text
[ epsilon  0       ]
[ 0        epsilon ].                              (14)
```

In either case the tied optimum consists of two opposite corners, while
the product of its projections contains all four corners. This contradicts
(D), proving (D) implies (C). QED.

### Score-table form and the two-colour guardrail

For a mixed `T`, THM-2346 writes

```text
H_T=sum_(c in C)lambda_(c,T) z_c^(tensor |T|).      (15)
```

Its exact Gram-rank calculation gives

```text
m>=3:
  H_T=0 iff lambda_(c,T)=0 for every c;

m=2, C={a,b}:
  H_T=0 iff
  lambda_(a,T)+(-1)^|T| lambda_(b,T)=0.            (16)
```

Thus (5) has an exact test in the supplied score-table coordinates.
Individual blockwise Mobius coefficients are not invariant consequences
of rectangularity. In particular, for two blocks,

```text
mu_a({p,q})=-mu_b({p,q})!=0
```

can be globally unary. The intrinsic objects are `H_T`, (7), and the
combination in (16), not the two coefficients separately.

## 2. Why one rectangular optimum is insufficient

For one unperturbed energy let

```text
O=argmin_Y E_bar.                                  (17)
```

The following are equivalent:

```text
O=product_P pr_P(O);

O is closed under replacing one prime-type coordinate
by that coordinate from any other point of O;

RectGap(E)
 :=max_(y in product_P pr_P(O))
      [E_bar(y)-min E_bar]
 =0.                                               (18)
```

These are exact tests for this **one** tied optimum. They do not imply the
robust conditions of Section 1.

Indeed, on two two-point type spaces the allocation table

```text
E_M=
[ 0  M ]
[ M  M ],                     M>0,                 (19)
```

has the singleton Cartesian optimum `{(a,a)}`, while

```text
H_{p,q}=-(M/4) h tensor h,

h(a)=1, h(b)=-1.                                   (20)
```

The hidden mixed interaction is unbounded as `M` grows. One realization
by block score tables is

```text
w_a({p})=w_a({q})=w_a({p,q})=0,

w_b({p})=w_b({q})=w_b({p,q})=M.                    (21)
```

Dominance can therefore hide arbitrarily large coupling behind a fixed
product optimum. The robustness quantifier in (D) is essential.

## 3. Exact target-token conditioning

Let `D subset S` be a set of **labelled target tokens**, fix an allocation

```text
a_D in C^D,
```

and put

```text
E^a(x_(S\D))=E(x_(S\D),a_D).                       (22)
```

If `H^a_U` denotes the ANOVA component of the restricted energy, then

> **Conditioning law.** For every `U subset S\D`,
>
> ```text
> H^a_U(x_U)
>  =sum_(V subset D)H_(U union V)(x_U,a_V).         (23)
> ```

Each term on the right retains zero marginal in every active coordinate
of `U`; grouping (4) by `U=T\D` and ANOVA uniqueness proves (23).
Sequential conditioning commutes. Uniformly averaging over all values of
the deleted labelled tokens kills every summand with `V!=emptyset`, so

```text
uniform marginalization over D sends H_U to H_U.   (24)
```

Equation (23) is labelled conditioning. Conditioning only a quotient count
uses the induced multinomial weights and is not silently identified with
(23).

### Block-score and Mobius contraction

For each block `c`, put

```text
D_c={d in D:a_d=c},

kappa=sum_c w_c(D_c),

w^a_c(A)=w_c(A union D_c)-w_c(D_c).                (25)
```

Then

```text
E^a=kappa+sum_c w^a_c(A_c),                        (26)
```

and the Boolean Mobius coefficients contract exactly as

```text
mu^a_c(T)
 =sum_(V subset D_c)mu_c(T union V)
                    for emptyset!=T subset S\D.    (27)
```

For `emptyset!=U subset S\D`, one transported representative in
THM-2346's centred-indicator coordinates is

```text
lambda^a_(c,U)
 =sum_(V subset D)
    lambda_(c,U union V)
    product_(v in V) z_c(a_v).                     (28)
```

The resulting tensor is gauge-invariant even when the representatives in
(28) are not unique.

This operation fixes target-token coordinates. THM-2339's deletion theorem
instead deletes a **source packet block** that is itself a nontrivial prime
target. They agree on the corresponding singleton self-owner face, but
they are not the same map. For a composite target, deleting an owner block
also requires that no other target token be sent to it; mere conditioning
does not supply that exclusivity. Equations (23), (27), and (28) are the
exact surviving obstruction.

## 4. Two distinct prime targets

Let `P,Q` be nonisomorphic oriented prime knots and let the source packet
blocks be `K_B`, `B in pi`. Define

```text
a_B=d_G(K_B,P)-u(K_B),

b_B=d_G(K_B,Q)-u(K_B),

mu_B
 =d_G(K_B,P#Q)-d_G(K_B,P)-d_G(K_B,Q)+u(K_B).
                                                               (29)
```

An ordered owner pair `(B,C)` has normalized allocation energy

```text
E(B,C)=a_B+b_C+1_[B=C]mu_B.                         (30)
```

Put

```text
alpha=min_B a_B,          A=argmin_B a_B=Own_pi(P),

beta=min_C b_C,           D=argmin_C b_C=Own_pi(Q),

g_B=a_B-alpha,            h_C=b_C-beta.             (31)
```

> **Exact product-owner criterion.**
>
> ```text
> Opt_pi(P#Q)=A times D                              (32)
> ```
>
> if and only if there is a scalar `theta` such that
>
> ```text
> 1_[B=C]mu_B=theta             for every (B,C) in A times D,
>
> g_B+h_C+1_[B=C]mu_B>theta
>                                for every (B,C) outside A times D.
>                                                               (33)
> ```

This is direct subtraction of `alpha+beta` from (30). Equality in the
second line creates an extra tied allocation; reversal of the strict
inequality beats the proposed product.

Let

```text
U_0=sum_(B in pi)u(K_B).
```

Under (33),

```text
Lambda_x(pi;P#Q)
 =U_0+alpha+beta+theta

 =Lambda_x(pi;P)+Lambda_x(pi;Q)-U_0+theta.          (34)
```

If `A times D` contains any off-diagonal pair, then `theta=0`. The only
way it contains no off-diagonal pair is

```text
A=D={B_0}.                                         (35)
```

Thus every product-owner case except a common unique owner forces exact
normalized additivity. In the exceptional case, `theta=mu_(B_0)` may be
nonzero precisely while every competitor in (33) stays strictly above it.

This is a conditional theorem about the displayed allocation energy and
its exact optimum. It is not unconditional additivity of unknotting number,
and THM-2339's one-owner deletion law does not imply (33).

## 5. The cohabitation term is a bounded continuation defect

For an arbitrary source knot `K`, define

```text
mu_K(P,Q)
 =d_G(K,P#Q)-d_G(K,P)-d_G(K,Q)+u(K).               (36)
```

Then

```text
-2 min(u(P),u(Q))
 <=mu_K(P,Q)
 <=2 min(u(K),u(P),u(Q)).                          (37)
```

For the lower bound, triangle inequality and connected-sum
nonexpansion give

```text
d_G(K,P#Q)>=d_G(K,P)-u(Q),

d_G(K,Q)<=u(K)+u(Q).
```

Substitution gives `mu_K(P,Q)>=-2u(Q)`; interchange `P,Q`.
For the upper bound,

```text
d_G(K,P#Q)<=d_G(K,P)+u(Q),

d_G(K,Q)>=|u(K)-u(Q)|,
```

so

```text
mu_K(P,Q)
 <=u(K)+u(Q)-d_G(K,Q)
 <=2 min(u(K),u(Q)).
```

Again interchange `P,Q`.

At the unknot `U`,

```text
mu_U(P,Q)
 =u(P#Q)-u(P)-u(Q)
 =-sigma(P,Q),                                     (38)
```

where `sigma` is THM-2176's nonnegative connected-sum interaction
cocycle. Thus `mu_K` is the target-row continuation of that old root
defect, not an unrelated coefficient. The bounds (37) are universal; no
sharp knot equality classification is claimed.

## 6. Scope and exact companion

The theorem supplies:

```text
global robust factorization  <->  mixed ANOVA vanishing;

target-token restriction     ->   exact tensor contraction;

two-prime product owners      <->  constant-inside / strict-outside gap;

root interaction sigma       ->   target-row continuation mu_K.     (39)
```

It deliberately does not turn a symmetric cohabitation tensor into a
tournament, infer decoupling from one optimum, recover gauge-dependent
block coefficients, delete a composite owner without exclusivity, or
prove connected-sum additivity.

The exact companion exhausts all `81` integer `2x2` tables in
`{-1,0,1}`, checks both opposite-corner normal forms, verifies the full
conditioning identity on a three-colour three-token bank, audits the
two-colour gauge hostile and an unbounded family of fixed-optimum
interactions, exhausts the two-prime owner criterion on a three-block
integer bank, and checks the triangle-bound algebra over a finite hostile
metric universe. Every load-bearing check raises explicitly under ordinary
and optimized Python.

Reproduce with

```bash
python3 04-computation/prime_type_rectangularity_conditioning_thm2348.py
python3 -O 04-computation/prime_type_rectangularity_conditioning_thm2348.py
```

Both transcripts must match

```text
05-knowledge/results/prime_type_rectangularity_conditioning_thm2348.out
```

byte-for-byte after LF normalization.
