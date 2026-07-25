---
id: THM-2340
title: "Owner-word ANOVA and exact target-support landing criteria"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For any
  complex 13 by 13 target-character twist matrix, normalized Fourier
  inversion splits its target coefficient energy orthogonally into the
  zero fibre, pure-a axis, pure-b axis, and mixed locus. These four
  energies are exactly the squared grand mean, ordered row main effect,
  ordered column main effect, and two-way interaction. Thus a pure-a
  target survives iff the row-mean profile is nonconstant; pure-b iff the
  column-mean profile is nonconstant; and a fork target survives iff the
  twist matrix is not additively separable. Applied to the full semantic
  current of THM-2334/2337, this is a lossless 169-current criterion for
  matching a terminal word to its target-support locus. It proves no
  component positive, excludes no scalar row, and does not prove LRC(14).
source: codex-2026-07-25-owner-word-anova
depends_on:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
related:
  - THM-2195-transitive-quotients-exactly-control-universal-substitution-products
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
script: 04-computation/lrc14_owner_word_anova_thm2340.py
output: 05-knowledge/results/lrc14_owner_word_anova_thm2340.out
script_sha256: e96788e30e648794c861c2617f8b66fbf46a67492f7662970ecbe0bc01a96946
output_sha256: 4055c3c280ce152928dc87f4f82d81576410ed1e4fd96ea32a5a934b475920f2
hash_basis: working-tree bytes (LF)
---

# THM-2340 -- pure owners are main effects; a fork is interaction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The independent audit rederived the normalized Fourier conventions, the
ordered main-effect identities, the double-centring/fork equivalence, the
complete orthogonal energy invoice, and all three support-mask Gram kernels.
It also checked target-swap covariance and replayed the exact companion under
ordinary and optimized Python with byte-identical output.

THM-2334 turns the full marked LRC current into a `13 x 13` matrix of
target-character twists.  Its total nonzero-target energy is the variance
of that matrix.  THM-2337 refines the target plane into the three
word-support loci

```text
pure a axis:       q_a!=0, q_b=0;
pure b axis:       q_a=0, q_b!=0;
fork locus:        q_a!=0, q_b!=0.                 (1)
```

The variance has an exact lossless refinement.  It is ordinary two-way
ANOVA on the twist matrix:

```text
pure a  <-> ordered row main effect,
pure b  <-> ordered column main effect,
fork    <-> irreducible row/column interaction.     (2)
```

This is the intrinsic replacement for forcing the fork
`j -> {a,b}` into a tournament edge.  Choosing one head marginalizes the
matrix and destroys its interaction.  ANOVA retains all four orthogonal
pieces and gives an iff target-landing criterion for each word type.

## 1. Normalized target Fourier transform

Let

```text
F=F_13,
zeta=exp(2*pi*i/13),
H:F x F -> C.                                      (3)
```

Think of `H(s,t)` as the full semantic current with target character
`(s,t)`.  Define its normalized inverse transform

```text
A(x,y)
 =1/13^2 sum_(s,t in F)zeta^(-sx-ty)H(s,t).        (4)
```

In the LRC application, `(x,y)` are THM-2309's ordered target coordinates
`(q_a,q_b)`, and (4) is exactly THM-2334's target-fibre inversion.  The
theorem itself is the abstract finite Fourier identity (3)--(4), so it
does not depend logically on the current construction.

Define the row means, column means, and grand mean

```text
R(s)=1/13 sum_t H(s,t),
C(t)=1/13 sum_s H(s,t),
mu=1/13^2 sum_(s,t)H(s,t).                         (5)
```

The doubly centred interaction is

```text
I(s,t)=H(s,t)-R(s)-C(t)+mu.                        (6)
```

It has zero mean on every row and every column:

```text
sum_t I(s,t)=0,
sum_s I(s,t)=0.                                    (7)
```

The exact reconstruction is

```text
H(s,t)
 =mu+[R(s)-mu]+[C(t)-mu]+I(s,t).                   (8)
```

The four summands in (8) are pairwise orthogonal for the normalized
Hermitian inner product on `C^(F x F)`.

## 2. The two pure axes are the ordered main effects

Set `y=0` in (4).  Summing first in `t` gives

```text
A(x,0)=1/13 sum_s zeta^(-sx)R(s).                  (9)
```

One-dimensional Parseval and `A(0,0)=mu` therefore give

```text
sum_(x!=0)|A(x,0)|^2
 =1/13 sum_s |R(s)|^2-|mu|^2
 =1/13 sum_s |R(s)-mu|^2.                         (10)
```

Consequently

```text
some pure-a target coefficient survives
 iff R is nonconstant.                             (11)
```

The ordered transpose argument gives

```text
sum_(y!=0)|A(0,y)|^2
 =1/13 sum_t |C(t)-mu|^2,                          (12)

some pure-b target coefficient survives
 iff C is nonconstant.                             (13)
```

The order matters.  Transposing `H` exchanges (10) and (12), exactly as
swapping blocker labels `a,b` exchanges the two pure words.

## 3. The fork is precisely the interaction

Let `A_I` be the inverse transform of `I`.  The three non-interaction
terms in (8) have Fourier support respectively at

```text
(0,0),
{(x,0):x!=0},
{(0,y):y!=0}.                                      (14)
```

Equations (7) say conversely that `A_I` vanishes on both coordinate axes.
Therefore

```text
A_I(x,y)=A(x,y)       when x!=0 and y!=0,
A_I(x,y)=0            when x=0 or y=0.             (15)
```

Two-dimensional Parseval yields

```text
sum_(x!=0,y!=0)|A(x,y)|^2
 =1/13^2 sum_(s,t)|I(s,t)|^2.                     (16)
```

Hence

```text
some mixed/fork target coefficient survives
 iff I is not identically zero
 iff H(s,t) is not of the form f(s)+g(t).           (17)
```

The last equivalence is exact: if `I=0`, take
`f=R` and `g=C-mu`; conversely every additive matrix has zero double
centering.

This is why the fork is not a tie to be oriented.  It is the part of the
current which no choice of one-target marginal can retain.

## 4. The complete orthogonal energy invoice

Combining (8), (10), (12), and (16) gives

```text
1/13^2 sum_(s,t)|H(s,t)|^2

 =|mu|^2
  +1/13 sum_s |R(s)-mu|^2
  +1/13 sum_t |C(t)-mu|^2
  +1/13^2 sum_(s,t)|I(s,t)|^2                      (18)

 =|A(0,0)|^2
  +sum_(x!=0)|A(x,0)|^2
  +sum_(y!=0)|A(0,y)|^2
  +sum_(x!=0,y!=0)|A(x,y)|^2.
```

The dimension ledger agrees term by term:

```text
grand mean:          1,
row main effect:    12,
column main effect: 12,
interaction:       144,

1+12+12+144=169.                                   (19)
```

In particular THM-2334's exact nonzero-target variance decomposes as

```text
sum_((x,y)!=(0,0))|A(x,y)|^2
 =row main-effect energy
  +column main-effect energy
  +interaction energy.                             (20)
```

Thus "`H` is nonconstant" proves that at least one of the three word
types lands, while (11), (13), and (17) say exactly which types.

## 5. Equivalence with the support-mask Gram forms

THM-2337 defines the three masks

```text
P_a(x,y)=1_(x!=0,y=0),
P_b(x,y)=1_(x=0,y!=0),
P_ab(x,y)=1_(x!=0,y!=0).                           (21)
```

Let

```text
h(d)=12 if d=0,
h(d)=-1 if d!=0.                                  (22)
```

The unnormalized finite transforms of the masks are

```text
P_a_hat(ds,dt)=h(ds),
P_b_hat(ds,dt)=h(dt),
P_ab_hat(ds,dt)=h(ds)h(dt).                        (23)
```

Therefore, for `sigma in {a,b,ab}`,

```text
sum_q P_sigma(q)|A(q)|^2
 =1/13^4 sum_(ell,ell')
   H(ell)conjugate(H(ell'))
   P_sigma_hat(ell'-ell).                          (24)
```

Substituting (23) into (24) and summing the unused character coordinate
gives exactly (10), (12), or (16).  Thus the `169 x 169` mask Gram form
can be evaluated without forming it:

```text
average the rows and test their variance;
average the columns and test their variance;
double-centre and test the residual norm.           (25)
```

This is the cheapest exact decisive test once the `169` twisted currents
are available.

## 6. Sharp equality and hostile controls

Every equality boundary is inhabited.

### Zero only

If

```text
H(s,t)=c!=0,
```

then `A(0,0)=c` and every nonzero target coefficient vanishes.  A nonzero
full current alone therefore forces none of the three landing energies.

### One pure owner only

If `H(s,t)=f(s)` with `f` nonconstant, then only the pure-a axis can
survive.  If `H(s,t)=g(t)`, only pure-b can survive.

### Both pure owners but no fork

If

```text
H(s,t)=f(s)+g(t)
```

with both centred functions nonzero, both main effects are positive but
`I=0`, so every mixed target coefficient vanishes.  Two active pure
transfers do not imply a fork.

### Fork only

If `f,g` are nonzero centred functions and

```text
H(s,t)=f(s)g(t),
```

then `R=C=mu=0` and `I=H`.  Only the mixed locus survives.

These are exact abstract matrices, not assertions that every one is
realized by a canonical LRC rectangle.  They prove that none of the three
iff tests can be weakened to total variance, two marginal variances, or
nonzero full current.

## 7. Tournament and information-loss audit

The intrinsic observable is not a pairwise comparison between `a` and
`b`.  It is the ordered two-coordinate response matrix.

```text
object:
  H(s,t), the full target-character current;

operations:
  row averaging, column averaging, and double centring;

preserved:
  target order, zero/pure/mixed support type, exact energy, and target
  swap covariance;

destroyed by choosing one head for a fork:
  the interaction I and all 144 mixed Fourier coordinates;

destroyed by forgetting target order:
  which 12-dimensional main effect belongs to pure a versus pure b;

cheapest decisive tests:
  the three norms in (10), (12), and (16).          (26)
```

This is analogous to the substitution boundary in THM-2195: an
unmarked quotient can look symmetric while the marked response carries
the obstruction.  Here the fork is an interaction component, not a
cosmetic tournament tie.

## 8. LRC specialization and remaining boundary

For the full marked word/deepest-comb/bare current of THM-2334 and
THM-2337, take

```text
H(s,t)
 =the boundary value of the target-character twist (s,t),
A(x,y)
 =the full-word target-fibre Abel aggregate.        (27)
```

Then:

```text
terminal word sigma={a} lands in its support type
 iff the row means of H are nonconstant;

terminal word sigma={b} lands in its support type
 iff the column means of H are nonconstant;

terminal word sigma={a,b} lands in its support type
 iff H has nonzero two-way interaction.             (28)
```

These are target-support statements.  They do not yet retain an
all-`91`-unit aggregate, a bounded visible address, or the
terminal-component phase of THM-2303.  Most importantly, the theorem does
not prove the relevant right side of (28) for a canonical strict-row
current.  It replaces a vague polarization problem by one exact finite
functional obstruction.

No scalar profile is excluded.  The exact ledger remains `165`, and
LRC(14) remains open.

## 9. Exact companion

The companion checks the ANOVA reconstruction, zero row/column means,
pairwise norm invoice, all three exact mask Gram identities, target-swap
law, and dimension ledger using rational arithmetic.  Its controls
separately inhabit zero-only, pure-a-only, pure-b-only, both-pure/no-fork,
and fork-only boundaries.  Every load-bearing check raises explicitly
under ordinary and optimized Python.

Reproduce with

```bash
python3 04-computation/lrc14_owner_word_anova_thm2340.py
python3 -O 04-computation/lrc14_owner_word_anova_thm2340.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_owner_word_anova_thm2340.out
```

byte-for-byte after LF normalization.
