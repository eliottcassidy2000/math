> **STATUS:** session synthesis and idea provenance.  Current mathematical
> truth is proved THM-3471, together with THM-3458, THM-3463,
> and THM-3468.  None of the three Rule 30 prizes is claimed here.

# Rule 30 after scalar randomness: transverse slack, marked carries, and macro renormalization

The session began with three tempting but incomplete pictures:

1. the nonperiodic valuation clock inside the center should resist
   cancellation;
2. independent innovation coordinates should push the marked center toward
   balance; and
3. large Hankel/digit/block ranks should make prediction slow.

Each picture failed for a precise reason.  The repair in all three cases was
to retain the operation applied **after** the attractive invariant: Green
slack, a marked phase arc, or charged scale growth.

## 1. The source-depth filtration is not the right scalar filtration

THM-3468 decomposed every nonlinear correction by source depth `u` and Green
slack `v`.  It proved each fixed `(u,v)` channel periodic and left open whether
the valuation clock could survive the infinite sum.  The natural next guess
was to sum every `v` at fixed `u` and hope a bounded set of source depths stayed
too simple to cancel the clock.

That guess fails as strongly as possible.

The complete fixed-`u` strip has one Motzkin diagonal compiler

```text
R_u(z,q)=z^(u+1) A_u(W_q)/(1+qz),
q^2W_q^2+(1+qz)W_q+z^2=0.
```

At `q=1`, all finite source-depth models live in one quadratic field.  The
first three forced strips satisfy

```text
R_0+R_1+R_2=0.
```

More surprisingly, `R_2` alone is exactly the Frobenius tail of the valuation
clock:

```text
R_2(z)=z^3G(z^2),
B+R_2=z^3/(1+z^2).
```

So the depth-two strip is the perfect correction which would turn the
density-`1/3` backbone into the balanced odd-time word.  It has density `1/6`
and is pointwise disjoint from the backbone.  The outer two strips reproduce
that same correction; the two copies cancel in characteristic two.

This is why another scalar nonperiodic component was the wrong target.  The
first three boundary branches are informative only before transport slack is
forgotten.

## 2. The first transverse jet is already nonzero

With slack marker `q`, the boundary circuit becomes

```text
R_0+R_1+R_2
 =(1+q) times a nonzero unit-series factor.
```

Its order at `q=1` is exactly one, and

```text
d/dq(R_0+R_1+R_2)|_(q=1)=R_1.
```

The first jet weights a source event by `v mod 2`; at a fixed target this is
just a checkerboard source-depth color.  It is the cheapest possible repair,
although genuinely finer transport information will need higher Hasse jets or
the full marker.

This is a literal instance of the repository rule “test filtration--observer
commutation before scalarizing.”  It also sharpens two older analogies:

- Berggren's three siblings need labels to reconstruct their parent; here the
  three boundary siblings become zero only after their transport labels are
  removed.
- Factorial Stokes/boundary-current work likewise warns that an unlabelled
  boundary can vanish while a face-labelled first response survives.

Neither analogy imports the other problem's theorem.  The common operation is
the labelled circuit.

## 3. Innovation independence hides a maximally complicated carry

The innovation atlas is a measure-space homeomorphism

```text
Gamma:Z_2 -> F_2^N.
```

This makes Haar-random phase look like independent fair coordinates.  But
phase `+1` becomes a triangular Boolean odometer.  At the `r`th innovation,
its carry `q_r` has odd truth-table weight.  One parity fact then forces three
global consequences:

```text
top ANF coefficient = 1,
degree(q_r)=r-1,
every Walsh coefficient is nonzero for r>=3.
```

Thus the coordinates are independent while the dynamics in those coordinates
has maximal algebraic carry degree and full Fourier support.  Independence of
the chart is not simplicity of the shift.

The cochain formulation makes the dichotomy exact:

- even current holonomy gives a coboundary and two primitives, differing by
  the owner constant;
- odd holonomy has no primitive on the base cube and creates the next cyclic
  double cover.

This is a cleaner explanation of “innovation” than period doubling alone.

## 4. The origin is not free, but it is a moving terminal arc

An unpointed odometer cannot distinguish `Gamma(0)` from any other codeword.
The physical light-cone boundary supplies the missing calibration:

```text
T_k(-k)=0,
c_k = xor_(h=-k)^(-1) Q_k(h).
```

At every finite cutoff, the boundary conditions force phase zero modulo the
exact seed period; the inverse limit selects phase zero uniquely.  So the
owner is not arbitrary.

But the prize-density debt remains.  The center is the parity of a **moving
terminal arc** of a current whose innovation-depth members are spectrally
full.  Total current weight, Haar balance, Walsh support, and innovation
holonomy do not determine that arc.
The physical left-front version says precisely which adjacent-parent OR events
are integrated.  A serious density probe must control their ordered terminal
discrepancy, not another phase average.

## 5. A real prediction advance, in a declared model

The complexity lane changed direction from lower bounds to charged batching.
For a block of `s` output cells and height `h`, precompute the map from its
`s+2h` raw cone.  Charging table construction and every lookup gives

```text
T=O(h2^(s+2h)+n^2/(sh)+n),
S=O(2^(s+2h)+n/log n) words.
```

On a `Theta(log n)` word-RAM, `s=2h=Theta(log n)` at half-word address size
gives

```text
T=O(n^2/log^2 n),
S=O(n/log n) words.
```

The table can emit the whole center prefix at the same asymptotic cost.  This
is a Four-Russians-style upper bound with preprocessing charged, not a Prize-3
solution and not a Turing-machine lower-bound statement.

The instructive hostile is internal: for every fixed tail, left permutivity
makes the universal block map a full-rank permutation on its owner prefix.
Perfect universal block rank coexists with fast batched evaluation of one
marked seed.  Rank is not time.

## 6. Three blocks simultaneously: an exact scale object

Group the binary row into `h`-bit symbols and sample every `h` time steps:

```text
C_h=pi_h F^h pi_h^(-1).
```

Then `C_h` is a radius-one CA over alphabet `F_2^h`; its local rule consumes
three adjacent `h`-bit parents simultaneously.  It is left-permutive in the
first parent and preserves uniform Bernoulli measure.  Across scales,

```text
C_(mh)=G_(m,h) C_h^m G_(m,h)^(-1).
```

This is an exact renormalization identity, but compression grows the alphabet
from `2^h` to `2^(mh)`.  It is neither a three-child arithmetic tree nor a
fixed finite renormalization.  The cumulative one-seed windows nevertheless
visit every local three-block address for `h<=6` in the audited finite
universe—another warning that universal local richness and marked temporal
complexity are different predicates.

## 7. What genuinely moved

The session did not solve a prize.  It did replace three vague strategies by
three exact targets:

1. **Nonperiodicity:** work in the `u>=3` completion while retaining slack
   jets; source depth without transport is provably too coarse.
2. **Density:** estimate the ordered terminal-arc discrepancy of `Q_k`, with
   the physical cone basepoint fixed.
3. **Prediction:** compress the width-changing block cocycle, or prove growth
   in a declared model which charges preprocessing and state alphabet.

The cross-prize invariant is now visible: every finite-scale quotient looks
simple—quadratic, Bernoulli, or triangular—but the physical observable moves
through an unbounded marked extension.  The likely next advance will not be a
better scalar statistic.  It will be a law for how that marked extension grows.
