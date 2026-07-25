---
id: THM-2346
title: "Global allocation ANOVA normal form and tournament boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Every
  finite allocation energy has a unique uniform Hoeffding/ANOVA expansion
  into zero-marginal token interactions. For the target-allocation energies
  of THM-2339, the global tensors are explicit centred-indicator transforms
  of the blockwise Boolean Möbius coefficients. This gives the exact global
  score-table gauge kernel, unary and pairwise iff tests, equivariance under
  equal-prime token permutations, and the exact two-colour Ising interaction
  in THM-2339's word-metric hostile. With blocks as relation vertices and an
  ordered token pair sampling a skew block observable, every allocation
  pair tensor is symmetric, so exact antisymmetric-relation reduction occurs
  only in the unary case. Other tournament encodings are not ruled out. No
  knot distance, unknotting number, catalyst, or classification is asserted.
source: codex-2026-07-25-global-allocation-anova
depends_on:
  - THM-2339-prime-owner-deletion-and-target-allocation-hypergraph
related:
  - THM-2336-prime-target-gordian-owner-diagram-and-bypass-split
  - THM-2340-owner-word-anova-target-landing
script: 04-computation/global_allocation_anova_thm2346.py
output: 05-knowledge/results/global_allocation_anova_thm2346.out
script_sha256: 1627b380344b0e61fead9d406cb26a23fcadc880c141a75b6f8f0b226b6b6b38
output_sha256: 01068f405e51f20c00fb7f561608f4cb989d4932df9c2585fd809dffef9aa939
hash_basis: working-tree bytes (LF)
---

# THM-2346 -- global allocation ANOVA normal form

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The independent audit rederived the commuting-projection expansion and
the centred-indicator Möbius bridge, including the triangular
lower-order terms.  It checked the Gram eigenvalues and their exceptional
`k=1` and two-colour kernels, the complete score-table gauge iff, the
equal-prime quotient, and the symmetric-versus-alternating tournament
boundary.  The two-token hostiles and all rank controls replay identically
under ordinary and optimized Python; the stored transcript and both LF
hashes match.

THM-2339 gives a geometrically canonical Boolean Möbius table for every
individual packet block, but it deliberately does not claim uniqueness after
the coloured hyperedges are summed over all allocations. This theorem takes
that final quotient.

The answer is the uniform Hoeffding/ANOVA expansion of the **global**
allocation energy. Its components are unique zero-marginal tensors. For the
special allocation carrier, they occupy a small diagonal Veronese, or
colour-anisotropic Potts, subspace. In particular, its binary interaction is
symmetric cohabitation energy, not an antisymmetric orientation.

## 1. The unique global zero-marginal normal form

Let `S` be a finite token set and let `C` be a finite colour set of size

```text
m=|C|>=2.
```

In the application, `S` is THM-2339's labelled list of oriented prime-target
occurrences and `C=pi` is the set of current partition blocks. Start more
generally with an arbitrary energy

```text
E:C^S -> R.                                             (1)
```

For `s in S`, let `A_s` be uniform averaging in coordinate `s`:

```text
(A_s E)(x)
 =1/m sum_(c in C) E(x with x_s replaced by c).         (2)
```

The operators `A_s` are commuting orthogonal projections for the uniform
inner product on `C^S`. For `T subset S`, define

```text
Pi_T
 =prod_(s in T)(I-A_s) prod_(s notin T)A_s,

H_T=Pi_T E.                                             (3)
```

Although written on `C^S`, `H_T` depends only on `x_T`. It has zero marginal
in every active coordinate:

```text
A_s H_T=0                    for every s in T.           (4)
```

For an explicit formula, put

```text
M_A E(x_A)
 =m^(-(|S|-|A|)) sum_(y in C^S, y_A=x_A) E(y).          (5)
```

Then

```text
H_T(x_T)
 =sum_(A subset T)(-1)^(|T|-|A|) M_A E(x_A).            (6)
```

> **Global ANOVA theorem.**
>
> ```text
> E(x)=sum_(T subset S)H_T(x_T).                        (7)
> ```
>
> The summands in (7) are mutually orthogonal, and (4) makes them unique.
> The `T`-summand space has dimension `(m-1)^|T|`.

### Proof

Because the projections commute,

```text
I
 =prod_(s in S)[A_s+(I-A_s)]
 =sum_(T subset S)Pi_T,                                 (8)
```

which proves (7). If `T!=U`, choose `s` in their symmetric difference.
One component lies in the image of `A_s` and the other in the image of
`I-A_s`, so they are orthogonal. Expanding the factors `I-A_s` gives (6).

The one-coordinate function space is

```text
R^C=<1> orthogonal-sum V,

V={f:C->R:sum_(c in C)f(c)=0},                          (9)
```

where `dim(V)=m-1`. The `T`-summand is the tensor product of one copy of `V`
for each `s in T` and constants elsewhere. Its dimension is
`(m-1)^|T|`. The direct-sum dimensions total

```text
sum_(T subset S)(m-1)^|T|=m^|S|,
```

so (4) and (7) leave no further freedom. QED.

The constant `H_emptyset` is retained in (7), so this is lossless for the
whole energy, not only for its minimizers. If only the owner set is retained,
adding a constant is harmless; many additional nonlinear changes can also
preserve an `argmin`. This theorem does not identify that much coarser
owner-only equivalence.

## 2. Exact interaction-order tests

Define the global interaction degree

```text
deg_A(E)=max({|T|:H_T!=0}),                             (10)
```

with the evident convention for a constant energy.

> **Order criterion.**
>
> ```text
> E is constant plus unary terms
>     iff H_T=0 for every |T|>=2;
>
> E is constant plus unary and pair terms
>     iff H_T=0 for every |T|>=3.                       (11)
> ```

Here the pair terms are arbitrary functions of the two indicated colours;
no symmetry or orientation is built into the phrase "pairwise".

There is also a base-free finite-difference test. For two colours `a,a'`,
let

```text
Delta_s^(a,a')E
 =E(x with x_s=a)-E(x with x_s=a').                    (12)
```

Then

```text
unary
 iff Delta_s Delta_t E=0
     for every distinct s,t and every colour choice;

pairwise
 iff Delta_s Delta_t Delta_u E=0
     for every distinct s,t,u and every colour choice.  (13)
```

The implications from the decompositions are immediate. Conversely, fix one
base colour in every coordinate and expand `E` by inclusion-exclusion over
the coordinates changed from that base point. The terms supported on two,
respectively three, changed coordinates are precisely the mixed differences
in (13). Vanishing of all differences of the stated order also kills every
higher difference containing one of them, leaving exactly the indicated
orders. Thus (11) and (13) are iff statements, not sufficient tests only.

## 3. From blockwise Möbius coefficients to global ANOVA

Return to THM-2339's allocation form. For every colour `c in C`, let

```text
w_c:2^S->R,                 w_c(emptyset)=0,            (14)
```

and put

```text
A_c(x)={s in S:x_s=c},

E(x)=sum_(c in C)w_c(A_c(x)).                           (15)
```

Write the Boolean Möbius coefficients as

```text
mu_c(B)
 =sum_(A subset B)(-1)^(|B|-|A|)w_c(A),

w_c(A)=sum_(nonempty B subset A)mu_c(B).                (16)
```

For a colour `c`, use the centred indicator

```text
z_c(d)=1_[d=c]-1/m.                                     (17)
```

It belongs to `V`, and

```text
1_[x_s=c]=1/m+z_c(x_s).                                 (18)
```

For nonempty `T subset S`, define

```text
lambda_(c,T)
 =sum_(B superset T)m^(-(|B|-|T|))mu_c(B).              (19)
```

> **Möbius-to-ANOVA bridge.**
>
> ```text
> H_emptyset
>  =sum_(nonempty B subset S)m^(-|B|)
>      sum_(c in C)mu_c(B),                             (20)
>
> H_T(x_T)
>  =sum_(c in C)lambda_(c,T)
>      prod_(s in T)z_c(x_s)              (T!=emptyset). (21)
> ```

### Proof

Equations (15) and (16) give

```text
E(x)
 =sum_(c in C)sum_(nonempty B subset S)
    mu_c(B) prod_(s in B)1_[x_s=c].                    (22)
```

Substitute (18). Choosing the constant `1/m` on `B\T` and the centred
indicator on `T` produces

```text
m^(-(|B|-|T|))mu_c(B) prod_(s in T)z_c(x_s).            (23)
```

This term depends on `T` and has zero marginal in each coordinate of `T`,
so uniqueness in Section 1 places it in `H_T`. Summing (23) over `B`
containing `T` proves (19)--(21). The choice `T=emptyset` proves (20). QED.

Thus each order-`k` allocation tensor lies in

```text
D_k=span{z_c^(tensor k):c in C}
     subset Sym^k(V).                                  (24)
```

The inclusion in `Sym^k(V)` identifies the colour spaces in the active
token coordinates. The coefficient can depend on the token subset `T`, but
within that tensor the colour arguments are symmetric.

Equation (19) records an important triangular effect. A blockwise hyperedge
`mu_c(B)` contributes not only to global order `|B|`; uniform marginalization
also sends it into every `H_T` with `T subset B`. Consequently a nonzero
blockwise pair coefficient is not by itself an iff test for a nonzero global
pair interaction.

Combining (11) and (21) gives exact allocation-specific tests:

```text
global unary
 iff sum_c lambda_(c,T) z_c^(tensor |T|)=0
     for every |T|>=2;

global pairwise
 iff sum_c lambda_(c,T) z_c^(tensor |T|)=0
     for every |T|>=3.                                 (25)
```

## 4. Rank and the exact score-table gauge kernel

Use the ordinary Euclidean inner product on `R^C` for this rank calculation.
The simplex contrasts satisfy

```text
<z_c,z_d>
 = (m-1)/m          if c=d,
 = -1/m             if c!=d.                          (26)
```

The Gram matrix of `{z_c^(tensor k):c in C}` therefore has diagonal

```text
a^k,             a=(m-1)/m,
```

and off-diagonal

```text
b^k,             b=-1/m.
```

Its two eigenvalues are

```text
eta_perp=a^k-b^k                 with multiplicity m-1,

eta_one =a^k+(m-1)b^k            with multiplicity 1.  (27)
```

It follows that

```text
dim(D_1)=m-1;

dim(D_k)=m           for m>=3 and k>=2;

dim(D_k)=1           for m=2 and k>=1.                (28)
```

For `m>=3,k>=2`, both eigenvalues in (27) are positive. For `m=2`,
`z_b=-z_a`, which gives the last line directly.

Now compare two supplied block-score families and prefix their differences
by `delta`. Equations (20)--(28) give the complete global gauge kernel.

> **Score-table gauge theorem.** The two families give exactly the same
> global function `E:C^S->R` iff `delta H_emptyset=0` and:
>
> - for `m>=3`, `delta lambda_(c,{s})` is independent of `c` for every
>   token `s`, while
>
>   ```text
>   delta lambda_(c,T)=0
>       for every c and every |T|>=2;                  (29)
>   ```
>
> - for `m=2`, writing `C={a,b}`,
>
>   ```text
>   delta lambda_(a,T)
>    +(-1)^|T| delta lambda_(b,T)=0
>       for every nonempty T.                          (30)
>   ```

For `|T|=1`, the only relation among the `z_c` is `sum_c z_c=0`, yielding
the first condition in (29). For larger `T`, (27) proves independence when
`m>=3`. Equation (30) is exactly
`z_b^(tensor |T|)=(-1)^|T|z_a^(tensor |T|)`. Together
with (20), these conditions kill every unique ANOVA component, and hence are
both necessary and sufficient.

The smallest hostile control to blockwise uniqueness is already
two-colour and two-token. Let `S={p,q}`, `C={a,b}`, and take

```text
mu_a({p,q})=1,          mu_b({p,q})=-1,
```

with every other coefficient zero. Both blockwise pair coefficients are
nonzero, but

```text
E=
  [ 1   0 ]
  [ 0  -1 ]

 =1/2 h(x_p)+1/2 h(x_q),                              (31)
```

where rows and columns are ordered `a,b` and

```text
h(a)=1,                 h(b)=-1.
```

Thus `H_{p,q}=0`: opposite blockwise pair hyperedges have become a purely
unary global field. This abstract finite-energy example is not asserted to
be realized by Gordian distances.

The binary allocation subspace is especially concrete:

```text
D_2 subset Sym^2(V),

dim(D_2)=1                       for m=2,
dim(D_2)=3=dim(Sym^2(V))         for m=3,
dim(D_2)=m<dim(Sym^2(V))
                                  for m>=4.            (32)
```

Thus two blocks give an Ising coupling, three blocks permit every symmetric
zero-marginal pair tensor, and four or more blocks give a strict
colour-anisotropic Potts subspace.

## 5. Equal-prime symmetry and the endpoint quotient

Let the token type of `s` be the oriented prime-knot type of `P_s`, and put

```text
G_J=product_P Sym(S_P),

S_P={s:P_s isomorphic to P}.                          (33)
```

For `g in G_J`, act on allocations by

```text
(g*x)_s=x_(g^(-1)s).                                  (34)
```

Because `J_(gA)` is isomorphic to `J_A`, THM-2339's geometrically supplied
score tables obey

```text
w_c(gA)=w_c(A),

mu_c(gA)=mu_c(A).                                     (35)
```

Uniform coordinate averaging commutes with this action. Therefore

```text
H_(gT)(g*x_T)=H_T(x_T),

lambda_(c,gT)=lambda_(c,T).                           (36)
```

The global normal form is equivariant, not merely invariant after summing
its orders.

An allocation orbit is determined exactly by the contingency table

```text
n_(P,c)=#{s in S_P:x_s=c},

sum_(c in C)n_(P,c)=|S_P|.                            (37)
```

Hence

```text
|C^S/G_J|
 =product_P binom(|S_P|+m-1,m-1).                     (38)
```

This is the same endpoint quotient as THM-2339, now visible at every ANOVA
order. Both `E` and its tied optimum set descend losslessly to these tables.
The quotient does not permute packet blocks: the colours remain the
geometrically distinct blocks of `pi`.

If every token has the same prime type, there is one tensor family for each
order `k`, transported among all `k`-subsets, and the full energy depends
only on the colour-count vector.

## 6. The exact tournament boundary

This section uses one precise tournament typing. The vertices are the packet
blocks, equivalently the colours `C`. Choose an order on each token pair
`{s,t}` only to say which assigned block is the first argument. A weighted
binary orientation observable is a zero-marginal matrix

```text
R_(s,t):C x C->R,

R_(s,t)(c,d)=-R_(s,t)(d,c).                          (39)
```

Nonzero signs orient the block pair; zeros remain ties. Reversing the token
pair reverses the sign gauge. This is a genuine antisymmetric relation, not
the transitive preorder obtained by sorting scalar energy values.

For a general energy, the canonical pair tensor has the orthogonal split

```text
H_{s,t}=Sym(H_{s,t})+Alt(H_{s,t}),

Sym(Q)(c,d)=[Q(c,d)+Q(d,c)]/2,

Alt(Q)(c,d)=[Q(c,d)-Q(d,c)]/2.                       (40)
```

> **Antisymmetric-relation criterion.** A general finite energy is exactly
> a constant plus unary terms plus observables of the form (39) iff
>
> ```text
> H_T=0                         for every |T|>=3,
>
> Sym(H_{s,t})=0                for every token pair.  (41)
> ```

Necessity follows by projecting the proposed representation. Conversely,
under (41), use the unique `H_emptyset`, `H_{s}`, and the already skew
`H_{s,t}` in (7). If every off-diagonal entry of a nonzero pair observable
is nonzero, its signs give a weighted tournament; otherwise it is an
antisymmetric relation with honest ties.

If one **single** block relation is required for all token pairs, there is
one further exact condition: all nonzero tensors `H_{s,t}` must be scalar
multiples of one fixed skew tensor. Without that specified common-relation
sidecar, pair-dependent orientations do not constitute one tournament.

For every allocation energy (15), however, (21) gives

```text
H_{s,t}
 =sum_(c in C)lambda_(c,{s,t}) z_c tensor z_c
 in Sym^2(V).                                         (42)
```

The symmetric and alternating tensor squares intersect only at zero. Hence:

> **Allocation tournament boundary.**
>
> ```text
> an allocation energy is exactly reducible to unary fields
> plus antisymmetric block relations of type (39)
>
> iff it is globally unary.                            (43)
> ```

In particular, every genuine nonzero allocation pair coupling is
cohabitation rather than orientation. When two equal-prime occurrences can
be transposed by `G_J`, (36) independently forces their pair tensor to be
symmetric, so indistinguishability itself also forbids orienting that pair.

Equation (43) does **not** rule out a tournament built from different
vertices, a separately supplied asymmetric observable, or a lossy
orientation of scalar comparisons. Such a proposal owes its own vertices,
pairwise observable, sign gauge, tie semantics, preserved target, and loss
sidecar.

## 7. The THM-2339 word-metric hostile in global coordinates

Use THM-2339's abstract conical word-metric hostile with tokens `p,q` and
blocks `a,b`. Subtract the common root baseline

```text
u(a)+u(b)=2
```

from its four lift costs. With the owner of `p` indexing rows and the owner
of `q` indexing columns, the two normalized allocation energies are

```text
E_0=
  [ 1  2 ]
  [ 2  1 ],

E_1=
  [ 0  2 ]
  [ 2  1 ].                                           (44)
```

For `h(a)=1,h(b)=-1`, exact uniform ANOVA gives

```text
E_0
 =3/2-1/2 h(x_p)h(x_q),                               (45)

E_1
 =5/4-1/4 h(x_p)-1/4 h(x_q)
      -3/4 h(x_p)h(x_q).                              (46)
```

Thus the canonical pair matrices are

```text
H^0_{p,q}
 =
 [ -1/2   1/2 ]
 [  1/2  -1/2 ],

H^1_{p,q}
 =
 [ -3/4   3/4 ]
 [  3/4  -3/4 ],                                     (47)
```

and

```text
H^1_{p,q}-H^0_{p,q}
 =-1/4 h tensor h.                                    (48)
```

Their squared norms under the uniform measure on `C^2` are

```text
||H^0_{p,q}||_2^2=1/4,

||H^1_{p,q}||_2^2=9/16.                               (49)
```

Both are pure symmetric Ising interactions. Their alternating/tournament
parts vanish exactly.

In THM-2339's blockwise coordinates,

```text
metric       mu_a({p,q})       mu_b({p,q})
------------------------------------------------
d_0                -1                 -1
d_1                -2                 -1.             (50)
```

For two colours, `z_b=-z_a=-h/2`, so (21) sends the row of (50) to

```text
1/4[mu_a({p,q})+mu_b({p,q})] h tensor h,              (51)
```

which is exactly (47).

The original hostile keeps every one-token block score fixed, but changing
the two-token block score also changes the uniform global mean and effective
unary fields:

```text
E_1-E_0
 =-1/4-1/4h(x_p)-1/4h(x_q)-1/4h(x_p)h(x_q).           (52)
```

This is not a contradiction. "One-token score" and "global unary ANOVA
component after averaging over the other token" are different coordinates.
The irreducible failure of global unary reduction is precisely the last
term in (52).

The owner information remains the same as in THM-2339:

```text
argmin(E_0)={aa,bb},

argmin(E_1)={aa}.                                     (53)
```

The ANOVA form explains the change as a strengthened same-block interaction
together with an `a`-favouring effective field. It asserts no realization of
this abstract hostile by knots.

## 8. Scope and the refined knot owner object

For each supplied packet, partition, and composite target, the exact chain is

```text
Gordian subset-score tables
 -> blockwise Boolean Mobius presentation
 -> global zero-marginal ANOVA tensors
 -> full energy on equal-prime contingency tables
 -> tied minimizing endpoint orbits.                  (54)
```

The map from the first to the second arrow is invertible block by block. The
second to third arrow takes the exact global gauge quotient described in
Section 4. The third stage is still lossless for the whole energy. Passing
only to interaction norms, an optimum, or a tournament shadow loses data.

The order profile

```text
I_k(E)=sum_(|T|=k)||H_T||_2^2                         (55)
```

is a gauge-invariant diagnostic, but the tensors themselves, not merely the
numbers `I_k`, are the lossless normal form. Equations (11), (21), and (43)
replace the informal questions "is the target interaction unary, pairwise,
or tournament-like?" by exact tests.

No value in (54) is computed without the Gordian-distance table required by
THM-2339. This theorem proves no new knot distance, unknotting number,
connected-sum law, catalyst, or knot classification.

## 9. Exact companion

The optimization-safe exact companion verifies:

- projection, zero marginals, orthogonality, and reconstruction;
- the Möbius-to-ANOVA bridge for deterministic exact score banks;
- the `m=2` blockwise-pair/global-unary gauge hostile;
- all rank cases in (28) by exact rational Gaussian elimination;
- equal-prime orbit counts and energy descent to contingency tables;
- the symmetric allocation pair sector and an independent skew
  three-colour tournament positive control;
- the matrices, decompositions, norms, and minimizers in (44)--(53).

Reproduce with

```bash
python3 04-computation/global_allocation_anova_thm2346.py
python3 -O 04-computation/global_allocation_anova_thm2346.py
```

and compare both transcripts with

```text
05-knowledge/results/global_allocation_anova_thm2346.out
```

This finite companion checks the exact identities and hostile controls. The
general proof is the commuting-projection argument, the centred-indicator
expansion, and the Gram-rank calculation above.
