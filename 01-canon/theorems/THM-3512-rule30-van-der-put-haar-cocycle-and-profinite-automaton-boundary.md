---
id: THM-3512
title: "Rule 30 van der Put--Haar cocycle and profinite-automaton boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  After
  removing the permanent low state bit, the Rule 30 odometer embedding has
  flat exact van der Put valuation on every dyadic shell.  Its normalized
  one-step displacement is an everywhere-odd coboundary whose deterministic
  dyadic block averages and Haar lifting details recover the innovation
  height.  The embedding has uniform strict derivative zero but is nowhere
  locally analytic.  Every fixed modular projection of the reduced van der
  Put word is eventually zero, while no finite least-significant-bit-first
  Mealy transducer computes the full embedding; every exact infinite
  synchronous realization, measured after N input digits, needs at least
  2^(N-E_(N+1))>=2^floor(N/3) reachable states.  A projective
  sibling ratio obeys
  an exact closed three-local renormalization.  No Rule 30 prize consequence
  is claimed.
source: root/rule30-normalized-displacement/van-der-put-haar/2026-08-16
audit: >
  PASS (2026-08-16), independent proof, scope, hostile, and replay audit.
  The auditor rederived the shell law, strict derivative and analyticity
  boundary, Haar lifting, exact prefix fibers, synchronous-state invoice,
  codimension identities, quasisymmetry equivalence, projective recurrence,
  unrestricted quotient hostile, and Mahler boundary.  An independent
  exhaustive prefix check covered all pairs through depth 12.  Ordinary,
  optimized, and stored transcripts are byte-identical.
depends_on:
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
  - THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary
  - THM-3503-rule30-odometer-ultrametric-regrading-and-orbit-closure-dimensions
  - THM-3507-rule30-normalized-dyadic-displacement-sibling-trace-and-assouad-spectrum
related:
  - THM-3468-rule30-radial-green-fold-innovation-discrepancy-and-fixed-seed-carrier-boundaries
  - THM-3500-rule30-dyadic-section-cut-defect-and-cross-depth-valuation-carrier
  - THM-3501-rule30-universal-cover-green-potential-and-slack-holonomy-seam
external:
  - Vladimir Anashin, "Automata finiteness criterion in terms of van der Put series of automata functions," arXiv:1112.5089 (2011; CITED for context, not needed for the internal necessary-direction proof)
  - Eric Rowland and Reem Yassawi, "Profinite automata," Advances in Applied Mathematics 85 (2017), 60--83, arXiv:1403.7659 (CITED for the profinite-automaton analogy only)
script: 04-computation/rule30_van_der_put_haar_automaton_thm3512.py
output: 05-knowledge/results/rule30_van_der_put_haar_automaton_thm3512.out
script_sha256: 1773f563166b4f0e36db62e7c82eaf791e8640be2f140bc6ca4e8cbce6f2a909
output_sha256: 9351016e53b914e27e0514b3127e68e86e25cd0cdf2a1afa6b6bff98d3b3a71d
hash_basis: raw bytes
---

# THM-3512 -- Rule 30 van der Put--Haar cocycle and profinite-automaton boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This theorem gives a literal non-Archimedean and automata interpretation of
THM-3507's normalized displacement.  The useful object is not a metaphorical
wave: it is first a van der Put coefficient of the actual odometer embedding,
then a lossless deterministic Haar lifting of consecutive orbit blocks.

## 1. Inheritance and clocks

Retain THM-3507's notation

```text
Phi(x)=x xor ((2x) or (4x)),      R_t=Phi^t(1),
iota:Z_2 -> X,                    iota(t+1)=Phi(iota(t)),
q_m=2^m,                          v_m=nu_2(iota(t+q_m)-iota(t)),
d_m=v_(m+1)-v_m,
U_m(t)=(iota(t+q_m)-iota(t))/2^(v_m).                 (1)
```

The valuation in (1) is independent of `t`, and every `U_m(t)` is an odd
2-adic unit.  The first innovation is `v_0=1`.  Every packed state is odd, so
define the desuspended embedding and its innovation height by

```text
J(t)=(iota(t)-1)/2 in Z_2,
h_m=v_m-m-1.                                          (2)
```

THM-3507 gives

```text
nu_2(J(t)-J(s))=v_m-1=m+h_m
when m=nu_2(t-s),                                     (3)

h_0=0,                 h_(m+1)-h_m=d_m-1>=0.          (4)
```

THM-3463's no-`111` law says that two consecutive gaps cannot both equal one.
Equivalently, the nonnegative increments in (4) have no consecutive zeros.
Therefore

```text
boxed: h_m>=floor(m/2).                               (5)
```

The inheritance pass is: THM-3507's signed sibling trace is the closest
mechanism; fixed-precision automaticity is the corrected near miss; the
canonical hostile is a synchronous transducer forced to buffer every
unemitted phase bit; and the least-used sidecars are the van der Put owner,
the signed Haar detail, and the first unit after projective cancellation.

## 2. The normalized displacement is a van der Put shell

For `n>=1`, let `n^-` be the integer obtained by deleting the highest nonzero
binary digit of `n`.  The van der Put coefficient of a continuous
`F:Z_2->Z_2` at `n` is

```text
C_n(F)=F(n)-F(n^-).                                   (6)
```

If `n=2^m+r` with `0<=r<2^m`, then `n^-=r`.  Applying (1)--(2) gives the
literal identity

```text
boxed:
C_(2^m+r)(J)=J(r+2^m)-J(r)
             =2^(v_m-1) U_m(r)
             =2^(m+h_m) U_m(r).                      (7)
```

Thus every coefficient on the complete dyadic shell
`[2^m,2^(m+1))` has the same exact valuation `m+h_m`; after the standard
1-Lipschitz normalization,

```text
b_(2^m+r):=C_(2^m+r)(J)/2^m=2^(h_m)U_m(r),
nu_2(b_(2^m+r))=h_m.                                  (8)
```

The entire van der Put sequence reconstructs `J`.  Passing only to the shell
valuation `h_m` destroys the oriented unit and phase owner; `U_m(r)` is the
necessary sidecar.  The maps `J` and `-J` are the cheapest abstract hostile:
they have the same shell valuations and opposite owners.

Equation (3) also proves directly that `J` is 1-Lipschitz.  The reduced
coefficients in (8) are therefore the standard reduced van der Put
coefficients of a compatible 2-adic map.

## 3. Uniform zero strict derivative, but no local analyticity

For distinct `s,t`, with `m=nu_2(t-s)`, equation (3) gives

```text
nu_2((J(t)-J(s))/(t-s))=h_m.                           (9)
```

By (5), the right side tends to infinity uniformly as `s` and `t` approach
one another.  Hence `J` has uniform strict derivative

```text
boxed: J'(t)=0 for every t in Z_2.                    (10)
```

This does not make `J` locally constant: `iota` is a homeomorphism and the
affine desuspension in (2) is injective.  It follows that `J` is nowhere
locally analytic.  Indeed, if `J` were represented by a convergent power
series on a nonempty ball, its derivative series would vanish at every point
of that ball by (10).  Uniqueness of a characteristic-zero power series would
then kill every nonconstant coefficient, making `J` constant there, contrary
to injectivity.

This is a genuine ultrametric boundary.  A nonconstant function may be more
than first-order flat at every point because `Z_2` has no connected intervals.
No real mean-value principle may be imported.

## 4. Exact coboundary and deterministic Haar lifting

Since `v_0=1`, the base normalized velocity is the exact odometer coboundary

```text
U_0(t)=J(t+1)-J(t).                                   (11)
```

Telescoping any dyadic orbit block, first at scale zero and then at arbitrary
scale `m`, gives

```text
boxed:
2^(-r) sum_(j=0)^(2^r-1) U_m(t+j2^m)
 =2^(h_(m+r)-h_m) U_(m+r)(t).                         (12)
```

The quotient on the left is an integer.  Its exact 2-adic valuation is
`h_(m+r)-h_m`.  In particular,

```text
2^(-m) sum_(j=0)^(2^m-1) U_0(t+j)=2^(h_m)U_m(t).      (13)
```

At `t=0`, (13) is the dyadic Volkenborn average of `U_0`.  It converges to
zero by (5), and the convergence is uniform in the translated starting phase.
The everywhere-odd reduction `U_0 mod 2=1` does not see this cancellation;
reduction modulo two and division by the growing block length do not commute.

For the lossless two-channel lifting, put

```text
L_m(t)=(U_m(t)+U_m(t+q_m))/2
      =2^(d_m-1)U_(m+1)(t),
V_m(t)=(U_m(t+q_m)-U_m(t))/2.                         (14)
```

Then

```text
boxed:
U_m(t)=L_m(t)-V_m(t),
U_m(t+q_m)=L_m(t)+V_m(t),                             (15)

V_m(t) is even iff d_m=1,
V_m(t) is odd  iff d_m>1.                             (16)
```

For the multilevel detail, subtract the left block sum from the adjacent
right block sum and normalize in the usual deterministic Haar fashion:

```text
D_(m,ell)(t)
 =2^(-ell-1) [
    sum_(j=2^ell)^(2^(ell+1)-1) U_m(t+j2^m)
   -sum_(j=0)^(2^ell-1) U_m(t+j2^m)] .                (17)
```

Applying (12) to the two halves gives

```text
boxed:
D_(m,ell)(t)=2^(h_(m+ell)-h_m)V_(m+ell)(t).           (18)
```

Equations (12), (15), and (18) are a lossless block lifting scheme once the
gap and detail are retained.  They are **not** claimed to be a Kozyrev
`L^2(Q_2)` wavelet expansion: the samples are consecutive points of an
odometer tower, the values lie in `Q_2`, and no positive energy identity is
used.  Calling (17) a deterministic Haar detail records only its exact
sum/difference algebra.

The fractal dimensions acquire a cancellation-rate form.  From
`v_m=m+1+h_m` and THM-3503,

```text
dim_H X=1/(1+limsup_(m->infinity) h_m/m),
dim_P X=1/(1+liminf_(m->infinity) h_m/m),              (19)
```

with `1/(1+infinity)=0`.  Thus the ambient compression of the orbit closure
is exactly the asymptotic 2-adic cancellation of the coboundary averages.

## 5. Fixed-modulus automaticity is a diagonal false friend

Fix `K>=1`.  Equations (5) and (8) imply

```text
b_n=0 mod 2^K for every n>=2^(2K).                    (20)
```

Every fixed modular projection of the reduced coefficient word is therefore
eventually zero, hence trivially ultimately periodic and 2-automatic.  But
the full `Q_2`-valued set `{b_n}` is infinite, because its exact valuations
`h_m` are unbounded.

This is precisely the useful part of the profinite-automaton analogy: each
fixed quotient is finite, while no uniform finite state set survives the
inverse limit.  Rowland--Yassawi are cited only for that established general
framework.  No algebraicity or automaticity of the adaptively normalized
Rule 30 owner word `U_m(r)` is asserted.

The destroyed coordinate in (20) is the first nonzero digit after the moving
height `h_m`.  The faithful adaptive object is the pair

```text
(h_m, 2^(-h_m)b_(2^m+r))=(h_m,U_m(r)).                (21)
```

This is the coefficient version of THM-3458's moving-observer warning.  A
fixed precision becomes silent exactly because the informative precision
moves outward with the scale.

## 6. Finite synchronous Mealy automata are impossible

Consider a deterministic synchronous Mealy transducer reading and writing
binary digits least significant first.  If it computes `J`, then after the
`m`-bit input prefix `r`, the two inputs `r` and `r+2^m` enter the same state
before receiving the respective tails `0...` and `10...`.  Their normalized
output difference is exactly

```text
[J(r+2^m)-J(r)]/2^m=b_(2^m+r).                        (22)
```

It depends only on the reached state.  A finite state set would make the
range of the reduced coefficients finite, contradicting (8) and (5).  This
is also the necessary direction of the cited van der Put automata criterion,
but (22) is a complete internal proof for the present map.

There is a sharper quantitative invoice.  After `N` input digits, associate
to each of the `2^N` phase prefixes the pair

```text
(first N output digits, reached transducer state).    (23)
```

The pair map is injective.  If two prefixes gave the same pair, appending the
same infinite tail would give the same full output, contradicting injectivity
of `J`.  The number of possible output prefixes is

```text
|J(Z_2) mod 2^N|=|X mod 2^(N+1)|=P_(N+1)=2^(E_(N+1)). (24)
```

In fact the fibers are exact.  For `0<=r,s<2^N`,

```text
J(r)=J(s) mod 2^N
 iff iota(r)=iota(s) mod 2^(N+1)
 iff r=s mod P_(N+1).                                 (24a)
```

The last equivalence is the exact width-`N+1` seed period.  Thus every output
prefix has exactly `2^(N-E_(N+1))` input prefixes above it.  Equivalently, the
output is a bijective code of the low `E_(N+1)` phase bits, while the remaining
`N-E_(N+1)` bits are an exact information backlog.

Consequently every exact infinite synchronous realization has, after reading
`N` input digits,

```text
boxed:
# reachable states at depth N >=2^(N-E_(N+1)).        (25)
```

Here `E_(N+1)` counts the ones in the length-`N` word
`epsilon_1...epsilon_N`.  The no-`111` law gives
`E_(N+1)<=ceil(2N/3)`, so

```text
# reachable states at depth N >=2^(floor(N/3)).       (26)
```

The meaning is literal: the output prefix has emitted only `E_(N+1)` phase
bits, and a synchronous device must buffer the remaining dimension-defect
bits in its state.  An asynchronous weighted transducer can escape this
no-go: consume one phase bit and emit a block of `d_m` state levels.  This is
why THM-3507's finite owner/borrow graph must have edge lengths; (25) rules
out the unweighted synchronous substitute, not the proposed weighted gap
certificate.

This proves no lower bound for computing the one fixed Rule 30 center bit.
It concerns the complete injective phase-to-packed-state embedding.

The backlog lower bound has an exact fractal-codimension exponent.  Using
THM-3503 and ignoring the harmless one-index shift,

```text
limsup_N (N-E_(N+1))/N = 1-dim_H X,
liminf_N (N-E_(N+1))/N = 1-dim_P X.                  (26a)
```

This identifies the exponent of the proved fiber/state lower bound, not the
minimal state complexity of a uniform transducer: no matching uniform upper
construction is asserted.

## 7. Bounded gaps are exactly quasisymmetric coding

The radial distance law makes the quasisymmetric boundary exact:

```text
boxed:
iota:Z_2->X is quasisymmetric
 iff J:Z_2->(X-1)/2 is quasisymmetric
 iff sup_m d_m<infinity.                              (27)
```

If `d_m<=G`, a ratio of input distances `2^k` becomes at most `2^(Gk)`;
for ratios at most one, `d_m>=1` gives the identity bound.  Thus one may use

```text
eta(t)=t       for 0<=t<=1,
eta(t)=t^G     for t>=1.                              (28)
```

Conversely, apply a quasisymmetry control to the triple
`(0,2^m,2^(m+1))`.  Its input distance ratio is `2`, while its output distance
ratio is `2^(d_m)`, so `2^(d_m)<=eta(2)` uniformly.

Together with THM-3507, (27) adds the equivalent condition

```text
quasisymmetric embedding
 iff dim_L X>0
 iff X is uniformly perfect.                         (29)
```

Thus a finite weighted graph with bounded edge cost would be not only a
dimension certificate but a bounded-delay, quasisymmetric coding certificate.

## 8. Projective sibling ratio and closed three-local renormalization

Remove the common odd amplitude of a sibling pair by defining

```text
G_m(t)=-U_m(t+q_m)/U_m(t) in Z_2^*.                  (30)
```

The sibling trace becomes

```text
boxed:
1-G_m(t)=2^(d_m) U_(m+1)(t)/U_m(t),
d_m=nu_2(1-G_m(t)).                                   (31)
```

The first unit after cancellation is

```text
Z_m(t)=(1-G_m(t))/2^(d_m)=U_(m+1)(t)/U_m(t).          (32)
```

It is always odd.  The filtration index `d_m`, not its mod-two graded value,
is the information: reducing (32) modulo two always gives one.

Let `G_j=G_m(t+jq_m)` and `A_j=U_m(t+jq_m)`.  Then
`A_(j+1)=-G_j A_j`, while the two parent units are
`A_0(1-G_0)/2^d` and `A_2(1-G_2)/2^d`.  Since
`A_2=G_0G_1A_0`, their projective ratio gives the closed recurrence

```text
boxed:
G_(m+1)(t)
 =-G_m(t)G_m(t+q_m)
   [1-G_m(t+2q_m)]/[1-G_m(t)].                       (33)
```

Equation (33) preserves every future gap without retaining the common
amplitude.  To reconstruct the signed displacements, one must restore one
base amplitude on each `q_m`-translation residue and then use (30); odd
rescaling is the cheapest hostile to the projective quotient.

### 8.1 Adaptive precision is unavoidable

Knowing `G mod 2^K` cannot distinguish `d=K` from any larger cancellation:
all such values equal one modulo `2^K`.  More strongly, recovering
`Z=(1-G)/2^d mod 2^K` requires `G mod 2^(d+K)`.  On unrestricted odd-unit
inputs, (33) does not descend to residues modulo `2^K` for any `K>=2`;
modulo two it is trivially constant.  This abstract nonclosure does not prove
that every physical Rule 30 residue quotient fails.  A quotient that resolves
the exact cancellation must retain

```text
d (or an explicit overflow state) + Z mod 2^K,        (34)
```

or equivalent information.  Whether a smaller enriched quotient closes on
the physical Rule 30 subshift remains open.  This is the exact adaptive-
precision tariff behind THM-3507's requested owner/borrow state.

There is a useful torsion split:

```text
d_m=1  iff G_m=-1 mod 4,
d_m>=2 iff G_m= 1 mod 4.                              (35)
```

On the second branch, the 2-adic logarithm is faithful and
`nu_2(log G_m)=d_m`.  It is not faithful on the first branch:
`log(-1)=0`.  A logarithmic quotient therefore needs the mod-four sign as a
sidecar.  Equivalently, for odd siblings `A,B`,

```text
AB mod 8 is 1 or 5 for d=1,
            3      for d=2,
            7      for d>=3.                          (36)
```

### 8.2 Scalar ratio collapse is impossible

If `G_m(t)` were independent of phase at any scale, (33) would reduce to

```text
g -> -g^2.                                            (37)
```

For every odd `g`, `nu_2(1+g^2)=1`.  Hence (37) would force
`d_(m+1)=d_(m+2)=1`, contradicting no-`111`.  Thus every scale has genuine
projective phase variation.  Any scalar-ratio graph is a demonstrated
underquotient.

## 9. Mahler is a lossy scale mixer for this target

The van der Put result does not transfer coefficientwise to the Mahler basis.
Already

```text
J(0)=0,          J(1)=3,          J(2)=12,
C_2(J)=J(2)-J(0)=12,              nu_2(C_2)=2,
a_2=J(2)-2J(1)+J(0)=6,            nu_2(a_2)=1.         (38)
```

Thus the second Mahler coefficient mixes the lower shell into the level-one
increment and destroys its exact innovation valuation.  Mahler expansions
remain valid global coordinates for continuous 2-adic functions, but they do
not carry the flat shell invariant without a triangular lower-level sidecar.

## 10. Preservation, loss, and scope

| source -> target | preserved | destroyed | required sidecar / hostile |
|---|---|---|---|
| `J` -> full van der Put word | all values of `J` | nothing | full coefficient sequence |
| van der Put word -> shell height | innovation metric and dimension slopes | phase and signed owner | `U_m(r)`; hostile `J` versus `-J` |
| sibling pair -> low/high Haar lift | lossless signed pair | nothing | gap plus signed detail |
| reduced word -> fixed `mod 2^K` word | fixed output precision | moving first live digit | `h_m` and adaptive unit |
| `J` -> synchronous Mealy state | exact online digits | none if state is unbounded | depth-`N` backlog `N-E_(N+1)` |
| sibling pair -> `G_m` | every gap and projective evolution | common odd amplitude | one base amplitude per residue |
| `G_m` -> fixed residue | cancellations below cutoff | exact large gap and normalized unit | overflow plus `Z_m` |
| van der Put -> Mahler | the continuous function globally | flat shell typing coefficientwise | triangular lower-shell data |

No tournament is intrinsic.  The faithful vertices are phase-tree edges,
consecutive orbit blocks, and projective sibling ratios.

Run

```bash
python3 04-computation/rule30_van_der_put_haar_automaton_thm3512.py
python3 -O 04-computation/rule30_van_der_put_haar_automaton_thm3512.py
```

and compare both byte-for-byte with
`05-knowledge/results/rule30_van_der_put_haar_automaton_thm3512.out`.
The companion contains no `assert` gates.  Its finite universes are
implementation controls for the exact proofs above.

This theorem proves no eventual nonperiodicity or balance of the Rule 30
center column, no fixed-seed query lower bound, no unbounded or bounded actual
innovation gap, no algebraicity or nonalgebraicity of a Rule 30 generating
function, and no Rule 30 prize.  It makes no literature novelty claim.
