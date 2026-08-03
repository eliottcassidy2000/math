# Constrained norm-phase transfer atlas and seam Cayley deletion

**Status: PROMOTED AS
[THM-3274](../01-canon/theorems/THM-3274-norm-fibre-constrained-phase-transfer-and-refinement-invoice.md),
VERIFIED-EXACT, AND INDEPENDENTLY HOSTILE-AUDITED.** A separate integer-coded
finite-field engine reproduced the refinements, spectra, Rabin certificate,
Fourier resultants, derivative identity, recurrences and hostiles without
importing the target companion. This note records an abstract finite-field
walk calculation. It does not construct a physical LRC ancestry walk,
endpoint current, safety certificate, row exclusion, or decrement for
`LRC(14)`.

## 1. Inheritance pass and concept board

- Closest proved mechanism: THM-3268 proves that freely choosing any nonzero
  increment gives the equitable quotient `14J_12-I_12` and a two-mode closed
  form.
- Canonical hostile: repeating the single translation `(1,0)` agrees with the
  free walk at one aggregate step but returns in characteristic 13 and does
  not iterate the quotient.
- Corrected near miss: THM-3267 proves that full `q` retains the chosen norm
  phase, `[q]` retains only parity, and `det(q,R)` retains none. Those are
  factorization statements, not a physical transition law.
- Least-used sidecar: THM-3246's negative second-corrector seam consists of
  the twelve exponents `0,...,5,162,...,167`, one in every norm phase and in
  twelve distinct projective directions.

The live concept board was:

| object / representation | question | lost coordinate | cheapest test |
|---|---|---|---|
| one multiplicative norm fibre as an additive increment set | is norm phase equitable? | the exact point in a phase fibre | compare all 14 sources in every phase |
| one projective scalar line | can parity or determinant repair phase? | affine-line coset | refine by the squared transverse determinant |
| the THM-3246 seam as a 12-letter increment alphabet | does one-per-phase imply a small transfer? | exact Singer placement | one outgoing colour-refinement step |
| full additive Cayley graph versus punctured graph | can the 168-state seam spectrum still be described exactly? | deleted zero | compare a direct characteristic polynomial with Fourier orbit factors and a principal-minor identity |
| one fixed increment | which varying-set conclusions survive? | increment word | repeat through lengths `2` and `13` |

The fixed-versus-varying distinction is load-bearing throughout:

```text
varying-set walk: choose a fresh t in S at every step, x+t != 0;
fixed-increment walk: choose one t once, repeat it, and stop on hitting 0.
```

## 2. Candidate A: a fixed norm fibre is equitable but not two-mode

Work in the chosen model

```text
F_169=F_13[u]/(u^2-2),       alpha=1+2u,
phi(alpha^j)=j mod 12,
X=F_169^*.
```

Let

```text
S_0={t in X:phi(t)=0},       |S_0|=14.
```

At every step choose a fresh `t in S_0`, subject to `x+t!=0`. The partition
of `X` by `phi(x)` is outgoing-equitable. In phase order `0,...,11`, its exact
quotient is

```text
Q_0 =
(0,0,2,2,2,0,2,2,0,2,1,0)
(0,0,2,0,2,0,0,2,2,2,2,2)
(2,2,1,2,0,2,1,2,2,0,0,0)
(2,0,2,2,2,0,0,2,0,0,2,2)
(2,2,0,2,0,2,2,0,1,2,1,0)
(0,0,2,0,2,2,2,2,2,0,2,0)
(2,0,1,0,2,2,2,0,1,0,2,2)
(2,2,2,2,0,2,0,0,0,2,0,2)
(0,2,2,0,1,2,1,0,2,2,0,2)
(2,2,0,0,2,0,0,2,2,2,0,2)
(1,2,0,2,1,2,2,0,0,0,2,2)
(0,2,0,2,0,0,2,2,2,2,2,0).
```

It is symmetric. Its row sums are `13` in source phase zero and `14` in the
other eleven phases. Its characteristic polynomial is

```text
p(lambda)=lambda^12-13lambda^11-77lambda^10+884lambda^9
          +1508lambda^8-19811lambda^7-2325lambda^6
          +160510lambda^5-96880lambda^4-341392lambda^3
          +243616lambda^2+27968lambda-23104.             (1)
```

The exact Rabin--Frobenius test proves that `(1)` is irreducible modulo `331`:

```text
x^(331^12)=x mod p,
gcd(x^(331^6)-x,p)=gcd(x^(331^4)-x,p)=1 mod 331.
```

Therefore `(1)` is irreducible over the rationals. Since `Q_0` is symmetric,
its spectrum consists of twelve distinct real algebraic roots. This is the
first sharp contrast with THM-3268: constraining to one norm fibre preserves
phase lumpability but expands the exact response from two modes to twelve.

If `C_n(d)` counts varying-`S_0` paths with terminal/source phase difference
`d`, every `C_n(d)` obeys the order-twelve recurrence defined by `(1)`. A
nonzero `12 by 12` Hankel determinant for each of the twelve `d` proves that
the scalar order is exactly twelve, not merely at most twelve. The determinant
tuple has digest

```text
92dc06c80eebcc5059a1c83719b987b3b4223b5b0aa3fbc842216ec3b0f0c934.
```

The first three profiles are

```text
n=0: (168,0,0,0,0,0,0,0,0,0,0,0),
n=1: (182,196,196,196,196,196,196,196,196,196,196,196),
n=2: (4522,2548,2548,2548,2548,2548,2548,2548,2548,2548,2548,2548).
```

The mechanism is multiplicative rescaling. For source `x` of phase `a`, write
`t=xz`; then `z` runs through the norm fibre of phase `-a`, and the target
phase is `a+phi(1+z)`. It depends on `a`, not on the point `x` within that
phase. Multiplying the whole walk by `alpha^{-r}` proves that every absolute
norm fibre `S_r` has a phase-shift conjugate of `Q_0`.

The fixed-increment hostile fails much earlier. With `t=1`, the same-phase
sources `alpha^0` and `alpha^12` land in phases `10` and `9`, respectively;
even the pointwise phase partition is not equitable. Repeating `t` gives

```text
n=2:  (12,14,14,14,14,14,14,14,14,14,14,14),
n=13: (156,0,0,0,0,0,0,0,0,0,0,0).
```

Thus `(1)` belongs to a fresh-choice norm-fibre operation, not to one fixed
translation.

## 3. Candidate B: a projective line needs exactly a squared transverse sidecar

Take the additive projective direction

```text
L^*=F_13^*(1,0)={(lambda,0):lambda!=0}.
```

Under varying-set semantics, norm phase alone is not equitable. Every source
phase has four distinct outgoing phase profiles. For example, the two
phase-zero points `alpha^0=(1,0)` and `alpha^12=(8,8)` see

```text
(1,0,2,0,2,0,2,0,2,0,2,0),
(1,2,0,0,0,1,2,0,2,2,0,2).
```

Write a point as `x=(a,b)`. An increment in `L^*` preserves `b`, while the
norm `a^2-2b^2` sees the transverse coordinate only through `b^2`. The exact
coarsest outgoing-equitable refinement of the phase partition is therefore

```text
(phi(x),b^2).                                             (2)
```

One colour-refinement step from phase produces precisely the cells in `(2)`,
and `(2)` is already equitable. There are 48 nonempty cells: 36 of size four
and 12 of size two.

The mechanism makes the spectrum transparent. The block `b^2=0` merges one
punctured affine line and has six phase cells; its quotient is a cell split of
the loopless `K_12`, with eigenvalues `11,-1`. Each of the six nonzero square
values merges the two lines `b=+-sqrt(s)`. Its seven-cell quotient is a cell
split of the loopless `K_13`, with eigenvalues `12,-1`. Hence the 48-state
quotient has spectrum

```text
12 with multiplicity 6,
11 with multiplicity 1,
-1 with multiplicity 41,                                (3)
```

and minimal polynomial

```text
(lambda-12)(lambda-11)(lambda+1).
```

Every fixed linear observable of the transfer obeys

```text
C_(n+3)=22C_(n+2)-109C_(n+1)-132C_n.                    (4)
```

For relative phase difference, the three exact closed forms are

```text
C_n(0)=(300*12^n+26*11^n+1858*(-1)^n)/13,
C_n(d)=(168*12^n-168*(-1)^n)/13,             d odd,
C_n(d)=(144*12^n+26*11^n-170*(-1)^n)/13,    d even, d!=0.
```

Their total is

```text
156*12^n+12*11^n,
```

as forced by twelve nonzero affine lines of size 13 and the punctured base
line of size 12. These formulas compute an exponential path family with two
integer exponentiations and one parity bit.

The result is direction-independent. For `L=dF_13`, multiplication by
`d^{-1}` carries the walk to the displayed model and preserves relative norm
phase. In original coordinates the normalized sidecar is

```text
(phi(x)-phi(d), (det(d,x)/Norm(d))^2).                   (5)
```

This is a concrete repair of THM-3267's determinant blindness: determinant
alone carries no global norm phase, but its normalized square is exactly the
missing state when the consumer is a walk constrained to one scalar line.
The consumer matters.

Again, fixing one scalar increment instead of varying over `L^*` discards the
48-state transfer and returns the characteristic-13 translation hostile.

## 4. Candidate C: the seam forces all 168 points, but its spectrum is still exact

Now reinterpret the THM-3246 negative-owner seam as the abstract increment
alphabet

```text
N={alpha^j:j=0,...,5,162,...,167}.                       (6)
```

This reinterpretation is not present in THM-3246 and carries no physical seam
semantics. It is only the constrained-set experiment requested after
THM-3268. Multiplying by `alpha^6` normalizes `(6)` to the consecutive Singer
arc `{alpha^0,...,alpha^11}`.

The set has one increment in each norm phase, but that balance supplies no
phase quotient. Within each of the twelve source phases, all fourteen points
have different outgoing phase profiles. Thus the first outgoing refinement
of the phase partition has 168 singleton cells. Consequently the coarsest
outgoing-equitable refinement refining phase is the discrete point/Singer-
exponent partition. An exhaustive `GL_2(F_13)` check also finds that the set
stabilizer of `N` is only the identity.

This is a precise stopping obstruction, not merely a failed guess: no proper
phase-refining equitable transfer exists for the varying-seam walk.

There is nevertheless a closed spectral representation. Restore the deleted
zero and let `B` be the 169-vertex additive Cayley adjacency with increment
set `N`; let `A` be its principal submatrix on `F_169^*`. Additive Fourier
characters diagonalize `B`. The 168 nontrivial characters split into fourteen
projective Galois orbits of size twelve, so

```text
chi_B(lambda)=(lambda-12)*product_[c] F_[c](lambda),     (7)
```

where `[c]` ranges over the fourteen projective character directions and each
`F_[c]` is an explicitly frozen monic degree-twelve integer polynomial. There
are thirteen distinct factors: the directions `(1,7)` and `(0,1)` have the
same factor.

Translation symmetry makes all 169 principal characteristic minors of
`lambda I-B` equal. Since the derivative of a determinant is the sum of its
principal minors,

```text
chi_A(lambda)=chi_B'(lambda)/169.                        (8)
```

The script derives `(7)` in `Z[zeta_13]`, applies `(8)`, and independently
computes the characteristic polynomial of the explicit `168 by 168`
adjacency. The coefficient tuples agree exactly. Their digests are

```text
chi_B: b8025e76edb12075836f20791c13418030c09e73e87b2011c25dba1f0b8a27e3,
chi_A: 653c5c1d8251e610ab744b6a72f5bba5e782cf2bb74738c0f1d33dbba1be1dd4.
```

The repeated Fourier-orbit factor survives differentiation. The exact
rational factor degrees of `chi_A` are `12` and `156`; the degree-twelve
factor is

```text
lambda^12-lambda^11+27lambda^10-66lambda^9+66lambda^8
+1039lambda^7-4029lambda^6+5784lambda^5+9062lambda^4
-55849lambda^3+110475lambda^2-102233lambda+41003.         (9)
```

The degree-156 cofactor digest is

```text
ccee063f7dd21de43d7502f4d1db02153cd78bf0dac92a65307824238cbfbf98.
```

An exact gcd with the derivative proves `chi_A` squarefree. Hence the
minimal polynomial of the full seam transfer has degree 168. Matrix entries
and all scalar observables satisfy the recurrence specified by `(8)`, but no
proper phase-compatible transfer state exists, and no low-order recurrence is
claimed for an arbitrary scalar observable.

Two hostiles sharpen the boundary.

1. Another one-per-phase transversal with exponents
   `(0,13,2,15,4,17,6,19,8,21,10,23)` also refines immediately to 168
   singletons, but its punctured characteristic-polynomial digest is
   `58aaacff1ebc624f20e95541ba3d214356389739137c0683feff88722c512609`.
   One-per-phase balance does not determine the spectrum.
2. Repeating the single seam increment `1` gives the length-13 profile
   `(156,0,...,0)`. It is not the varying-alphabet matrix `A`.

## 5. Connection contract and information-loss ladder

The three constrained sets form a clean information ladder:

| source walk and quotient | predicate preserved | information destroyed | exact repair / obstruction |
|---|---|---|---|
| varying fixed norm fibre, `x -> phi(x)` | every phase-to-phase path count | point within a norm fibre | no extra state; exact 12-state quotient, but irreducible order 12 |
| varying projective line, `x -> phi(x)` | only coarse totals | transverse affine-line coset | normalized squared determinant `(5)`; exact 48-state, three-mode transfer |
| varying seam alphabet, `x -> phi(x)` | phase labels but not transition counts | complete Singer placement | first refinement is discrete; retain all 168 points |
| full seam Cayley graph `B -> A=B without 0` | all punctured paths | translation diagonalization | derivative identity `(8)` restores an exact spectral description |
| any varying-set walk -> one fixed `t` | the alphabet contains `t` and may match one aggregate row | the increment word and branching | no quotient transfer; characteristic-13 hostile |

This is not a monotone relation between set sizes: the 14-element norm fibre
compresses to 12 states, the 12-element scalar line needs 48, and the
12-element seam needs 168. The operative invariant is the increment set's
interaction with multiplicative norm phase and additive stabilizers, not its
cardinality.

The projective result is the most useful reframe. THM-3267 correctly says
that `det(q,R)` is globally blind to the absolute norm phase. The constrained
line consumer changes the question: its increments preserve a transverse
determinant, and phase only needs the square of that determinant after
normalization. A globally blind coordinate can be the exact *dynamic* sidecar
for a specific operation.

## 6. Validation and reproduction

The companion freezes the exact universe `F_169^*`, the increment sets, the
zero exclusion, and both increment semantics. It includes:

- the THM-3268 all-increment quotient as a positive control;
- direct point-level dynamic programs independent of the norm and projective
  quotients;
- all twelve norm fibres, not only the displayed representative;
- an exact mod-331 irreducibility certificate and twelve nonzero Hankel
  determinants;
- one-step colour refinement and an independent equitability check for the
  48 projective cells;
- direct `168 by 168` characteristic computation independent of the seam's
  cyclotomic Fourier reconstruction;
- exhaustive `GL_2(F_13)` seam stabilizers;
- fixed-increment and phase-matched-transversal hostiles; and
- no assertion nodes, floating literals, randomness, or fitted recurrence.

Run

```text
python3 04-computation/lrc_constrained_norm_phase_transfer_atlas_20260803.py
python3 -O 04-computation/lrc_constrained_norm_phase_transfer_atlas_20260803.py
```

Both executions byte-match

```text
05-knowledge/results/lrc_constrained_norm_phase_transfer_atlas_20260803.out.
```

LF-normalized hashes:

```text
script: 00ac273344020cdc569b75ba314438497d5462f8871f3343e83419ceba95387a
output: 680359f4fd51bdfb89bfb66b7bfae709ebb1e96a453f3ac6ec73707cf93a7597
```

## 7. Honest frontier and next decisive tests

Nothing here supplies a physical increment set on a literal LRC ancestry
sheet. The chosen norm gauge is external, the seam-to-increment
reinterpretation is abstract, and path multiplicities do not control a
loneliness maximum, positive current, or row safety.

The next tests are now sharper:

1. Construct the missing same-ancestry `(q,R)` record from THM-3267 and list
   the *actual* lawful increments. Test whether they are a norm fibre, a
   scalar-line union, a seam-like rigid set, or none of these.
2. If the physical increments preserve one projective direction, carry the
   normalized squared determinant `(5)` before any phase quotient. It is the
   cheapest exact sidecar exposed here.
3. For unions of norm fibres, compute whether the twelve-state quotients add
   and whether their characteristic polynomials collapse. This is the natural
   cyclotomic-number interpolation between THM-3268's two modes and Candidate
   A's twelve.
4. For time-dependent prescribed increment sets, keep the bounded transfer
   product rather than fitting one stationary recurrence. Fixed and varying
   semantics already show why this distinction is necessary.
5. Independently reconstruct the seam factors and audit the exact claim that
   the rational factorization has degrees `12` and `156` before any canonical
   promotion.

Tournament language is not useful here: the norm-fibre and scalar-line
graphs are symmetric, while the seam graph has both missing and asymmetric
arcs. Forcing a total orientation would discard precisely the transition
counts being studied. The native objects are equitable partitions, additive
Cayley graphs, and transfer matrices.
