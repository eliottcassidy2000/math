# Tournament cut switching: centered coronal and total-walk compiler

**Status:** structural synthesis for
[THM-3315](../01-canon/theorems/THM-3315-tournament-cut-switching-centered-coronal-walk-compiler.md).
This note gives a loss-aware sequence compiler for a native tournament operation:
cut switching.  It is a theorem-scale algebraic identity, not a literature
novelty claim.  The finite computation is an independent hostile audit of
the formulas, not their proof.

Companions:

- [`tournament_switching_centered_coronal_walk_compiler_20260803.py`](../04-computation/tournament_switching_centered_coronal_walk_compiler_20260803.py);
- its [frozen exact output](../05-knowledge/results/tournament_switching_centered_coronal_walk_compiler_20260803.out).

The outcome is that the centered characteristic polynomial of a tournament
is common to its entire switching class, but it is not enough to compile the
ordinary total directed-walk sequence of a selected switched tournament.  A
single load-bearing sidecar repairs the loss: the signed-cut observer
numerator.

## Inheritance pass

**Closest proved mechanism.**  [THM-3248 -- Q4 paired-owner Stirling
compiler](../01-canon/theorems/THM-3248-q4-paired-owner-stirling-compiler.md)
retains both the determinant and the full walk numerator.  Its
determinant-only hostile is the correct warning here: an operator invariant
does not determine a sequence observed through a chosen vector.

**Canonical hostile.**  [THM-1415 -- switching is the canonical star
quotient](../01-canon/theorems/THM-1415-switching-is-the-canonical-star-quotient.md)
already shows that the transitive triangle and the directed triangle lie in
one switching class.  They have radically different long-walk behavior, so
any switching-class-only sequence claim must fail at order three.

**Corrected near miss.**  Treating switching as preservation of the ordinary
adjacency spectrum is false.  Switching conjugates the *centered* adjacency
operator; the all-ones rank-one term used to recover ordinary adjacency is
not conjugation-invariant.  Equally, the skew tournament matrix is not the
symmetric switching form of [THM-2501 -- switching fourth moment, signed
four-cycles, and Gram energy](../01-canon/theorems/THM-2501-switching-fourth-moment-signed-c4-and-gram-energy.md).

**Least-used relevant sidecar.**  The signed cut itself, viewed as an
observer vector rather than merely an operation, appears not to have been
used in the inherited join/substitution compilers.  Keeping it turns the
violent quotient of THM-1415 into a small exact interface.

## Exact identity

Let `T` be a labelled tournament on `n` vertices.  Use row-to-column
adjacency

```text
A_ij = 1 if i -> j, and A_ij = 0 otherwise,
```

and let `J` be the all-ones matrix.  Define the centered adjacency

\[
 C=2A-J.                                                       \tag{1}
\]

Thus `C_ii=-1`, while each off-diagonal entry records the intrinsic arc
orientation by a sign.  Choose a cut sign vector
`d in {+1,-1}^n`, put `D=diag(d)`, and identify `d` with `-d`.  Reversing
every arc between the two sign parts produces a tournament `T^d` with

\[
 C(T^d)=DCD,
 \qquad
 2A(T^d)=DCD+J.                                               \tag{2}
\]

Define the switching-invariant centered denominator and the signed observer
coronal

\[
 \begin{aligned}
 P(z)&=\det(I-zC),\\
 U_d(z)&=d^T(I-zC)^{-1}d=\frac{N_d(z)}{P(z)},
 \end{aligned}                                                \tag{3}
\]

where

\[
 N_d(z)=d^T\operatorname{adj}(I-zC)d                           \tag{4}
\]

has degree at most `n-1`.  Finally let

\[
 G_{T^d}(x)={\bf1}^T(I-xA(T^d))^{-1}{\bf1}
            =\sum_{k\geq0}w_k(T^d)x^k,                        \tag{5}
\]

where `w_k` is the number of directed walks of length `k`, with arbitrary
initial and terminal vertices.

Then

\[
 \boxed{
 \det(I-2zA(T^d))=P(z)-zN_d(z),
 \qquad
 G_{T^d}(2z)=\frac{N_d(z)}{P(z)-zN_d(z)}.}                    \tag{6}
\]

Consequently the pair `(P,N_d)` is a complete constant-coefficient compiler
for the total-walk sequence of the selected switch.  The common polynomial
`P` is the switching-class coordinate; `N_d` is the observer/cut sidecar.
Both are unchanged by the gauge `d -> -d`, and both transform covariantly
under simultaneous relabelling of `C` and `d`.

### Proof

Conjugating the first matrix in (6) by `D` and using `D1=d` gives

\[
 \begin{aligned}
 D(I-2zA(T^d))D
   &= I-zC-zdd^T.                                             \tag{7}
 \end{aligned}
\]

The matrix determinant lemma therefore yields

\[
 \det(I-2zA(T^d))
   =P(z)\left(1-zd^T(I-zC)^{-1}d\right)
   =P(z)-zN_d(z).                                             \tag{8}
\]

For the observer, Sherman--Morrison gives

\[
 \begin{aligned}
 G_{T^d}(2z)
 &=d^T(I-zC-zdd^T)^{-1}d\\
 &=U_d(z)+\frac{zU_d(z)^2}{1-zU_d(z)}
  =\frac{U_d(z)}{1-zU_d(z)}
  =\frac{N_d(z)}{P(z)-zN_d(z)}.                              \tag{9}
 \end{aligned}
\]

The identities are first valid in the rational-function field.  Equations
(4), (8), and (9) clear denominators and make them polynomial/formal-series
identities.  No invertibility assumption at a numerical value of `z` is
needed.  This proves (6).

## Finite interface and coefficient compiler

Write

\[
 U_d(z)=\sum_{k\geq0}u_kz^k,
 \qquad u_k=d^TC^kd,                                         \tag{10}
\]

and scale the target sequence by

\[
 s_k=2^kw_k(T^d),
 \qquad G_{T^d}(2z)=\sum_{k\geq0}s_kz^k.                     \tag{11}
\]

Equation (9) is equivalent coefficientwise to

\[
 \boxed{
 s_0=u_0=n,
 \qquad
 s_k=u_k+\sum_{i=0}^{k-1}u_i s_{k-1-i}\quad(k\geq1).}        \tag{12}
\]

This is an online convolution form of the compiler.  It also has a genuinely
finite interface.  If `P(z)=sum_i p_i z^i`, then the first `n` moments give

\[
 [z^m]N_d(z)=\sum_{i=0}^{m}p_i u_{m-i}
 \qquad(0\leq m<n),                                         \tag{13}
\]

and Cayley--Hamilton makes every later coefficient of `P U_d` vanish.
Thus `P` together with `u_0,...,u_(n-1)` determines `N_d`, hence the entire
sequence.  If

\[
 Q_d(z)=P(z)-zN_d(z)=\sum_{j=0}^{n}q_jz^j,                   \tag{14}
\]

then

\[
 \sum_{j=0}^{n}q_js_{k-j}=0\qquad(k\geq n),                 \tag{15}
\]

with zero padding if `Q_d` has smaller degree.  Dividing by `2^k` converts
this immediately to a constant-coefficient recurrence for the unscaled
walk counts.

This is complementary to the inherited sequence compilers:

- [THM-3181 -- tournament half-grid reciprocity and repeated-join
  recurrence](../01-canon/theorems/THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md)
  obtains C-finiteness at fixed path-cover depth under repeated order join;
- [THM-3202 -- C3 repeated-join moving-jet formula and C-finite
  obstruction](../01-canon/theorems/THM-3202-c3-repeated-join-moving-jet-formula-and-cfinite-obstruction.md)
  shows why a fixed jet fails when the observer depth moves;
- THM-3248 keeps the full quotient walk numerator through substitution.

Here the tournament order is fixed and the operation varies over its cut
cube.  The numerator sidecar is again indispensable, but diagonal
conjugation plus one rank-one update makes its transport exact.

## Uniform switching-cube averages

Average uniformly over the `2^(n-1)` cut signs modulo global negation.  This
is a labelled cut-cube average, not an average over distinct isomorphism
classes.  Since

\[
 \mathbb E_d[dd^T]=I,                                       \tag{16}
\]

we obtain

\[
 \boxed{
 \mathbb E_d U_d(z)=\operatorname{tr}(I-zC)^{-1}.}           \tag{17}
\]

Jacobi's determinant identity gives

\[
 \operatorname{tr}(I-zC)^{-1}=n-z\frac{P'(z)}{P(z)}.        \tag{18}
\]

Multiplying by `P` and then using (6) yields the two polynomial averages

\[
 \boxed{
 \mathbb E_d N_d(z)=nP(z)-zP'(z),}                           \tag{19}
\]

and

\[
 \boxed{
 \mathbb E_d\det(I-2zA(T^d))
   =(1-nz)P(z)+z^2P'(z).}                                   \tag{20}
\]

Thus the average observer numerator and the average ordinary adjacency
characteristic denominator collapse back to the switching-invariant
polynomial `P` and its derivative.  Individual target walk sequences do not.

## Sharp order-three hostile

Let `T` be the transitive tournament `1 -> 2 -> 3`, with `1 -> 3`.  Its
centered denominator is

\[
 P(z)=1+3z+6z^2+4z^3.                                      \tag{21}
\]

Switching at the source and switching at the middle vertex are cuts of the
same cardinality.  The first target remains transitive; the second target is
the directed triangle.  Their complete data are

| singleton switched | `N_d(z)` | `det(I-2zA(T^d))` | `(w_k)_(k>=0)` |
|---|---|---|---|
| source | `3+6z+4z^2` | `1` | `3,3,1,0,0,...` |
| middle | `3+6z+12z^2` | `1-8z^3` | `3,3,3,3,3,...` |

Therefore the first failed implication is already exact at three vertices:

```text
centered spectrum (indeed the entire switching class)
  + cut cardinality
does not determine the ordinary total-walk sequence.          (22)
```

The missing coordinate is `N_d`.  This example is minimal, since there is
only one tournament up to isomorphism below order three and no directed
cycle can occur.

## Source, target, map, and loss ledger

**Source.**  A tournament's centered switching class together with a signed
cut observer `(C,d)`, where `d~-d`.

**Target.**  The complete total directed-walk sequence, its rational ordinary
generating function, and a degree-at-most-`n` recurrence for the selected
switched tournament `T^d`.

**Map.**  Diagonal conjugation `C -> DCD`, followed by the rank-one recovery
`2A(T^d)=DCD+J`; equivalently, after conjugating back, the update is
`I-zC -> I-zC-zdd^T`.

**Preserved.**  `P(z)=det(I-zC)`, the centered spectrum with multiplicity,
and every function of it such as the switching-cube averages (19)--(20).

**Restoring sidecar.**  `N_d=d^T adj(I-zC)d`, equivalently the first `n`
signed centered moments.  It remembers how the selected cut observer sits
inside the centered operator.

**Destroyed when the sidecar is forgotten.**  The fixed all-ones observer,
the selected adjacency characteristic polynomial, and the ordinary walk
sequence.  The hostile (21)--(22) shows that cut size does not repair this.

**Not restored even by `(P,N_d)`.**  Hamiltonian path counts, path-cover
profiles, SCC order, order-join factorization, and substitution run content.
In particular this note supplies no converse to
[THM-1936 -- signed Redei join multiplicativity](../01-canon/theorems/THM-1936-signed-redei-join-multiplicative.md)
and no refinement of
[THM-1862 -- order-join reduction principle](../01-canon/theorems/THM-1862-order-join-reduction-principle.md).

**Cheapest decisive test.**  The two singleton cuts of the transitive
triangle above.  Any proposed switching compiler that omits observer data
must separate `3,3,1,0,...` from `3,3,3,3,...`; if it cannot, it has already
failed.

## Exact audit

The companion script uses only exact integer polynomial and matrix
arithmetic.  It exhausts all labelled tournaments through order five and
one representative of every cut sign modulo global negation:

```text
orders                 1 through 5
labelled tournaments   1,099
switch cases           16,933
observer moment cells  252,811
```

For every case it independently checks:

- centered conjugation and the `d~-d` gauge;
- direct matrix-power moments against `U_d=N_d/P`;
- polynomial truncation and Cayley--Hamilton closure;
- the direct determinant of `I-2zA(T^d)` against `P-zN_d`;
- direct target walks against the convolution (12);
- the recurrence (15); and
- the cut-cube averages (17)--(20).

Normal and optimized runs agree with the frozen transcript.  The three
frozen digests are

```text
case     f5edbe1308c1f43503d0f790acbdabb9b09dbd741d6eed7c3df2ce1eb0057ad4
average  8289c5c82f0194084af60048cf055219ada6739fb6d76887603f5477c7fcde6f
semantic 282a0c09458e9b539e349135f3800cd8bfc850971b4e8782fc14de75ea36c1be
```

The script AST contains no `assert` node and no floating literal.

## Named next probes

1. **Moment-side switching distributions.**  Determine which coefficients of
   `N_d` have closed Rademacher moment formulas beyond the mean (19).  The
   next hostile check should distinguish genuine centered spectral
   information from signed-cut placement information.
2. **Low-rank substitution compatibility.**  Test whether `(P,N_d)` composes
   through any nontrivial tournament substitution family.  The expected
   obstruction is that substitution needs multi-owner/run data, exactly as
   THM-3248 warns.
3. **Collision search.**  Find the least order at which distinct cut
   observers have identical `(P,N_d)` but different Hamiltonian or
   path-cover profiles.  Such a witness would sharply delimit the compiler
   without weakening its exact walk claim.
