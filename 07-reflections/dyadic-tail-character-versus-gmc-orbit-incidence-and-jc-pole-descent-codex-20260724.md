# Dyadic tail character versus GMC orbit incidence and JC pole descent

**Status:** proved mechanism comparison, stratified-completion transfer, and
scoped negative audits.
The fair-extraction theorem is [THM-2160](../01-canon/theorems/THM-2160-dyadic-checksum-extracts-a-fair-bit-under-the-critical-run-deadline.md).
The GMC and Jacobian conclusions below identify what that mechanism preserves
and loses; they are not new implications between the three problems.

## 1. What the coin proof actually uses

For critical run length `n`, let `h` be the largest power of two with `h<=n`.
The first `h` bits are constant, and the first `2h` bits are not. On the tail
`y=(y_1,...,y_h)`, THM-2160 uses

```text
S_h(y)=sum_(i=1)^h i y_i mod h
```

and separates the lower and upper halves of `Z/hZ`.

This is not one universal prefix involution. On the Hamming-weight-`j` shell,
write `j=2^a u` with `u` odd and rotate the tail by

```text
delta=h/2^(a+1).
```

The rotation preserves the Bernoulli monomial weight and changes `S_h` by
`j delta=h/2 mod h`. It therefore bijects the two outputs. The chosen rotation
depends on the sidecar `j`; the checksum alone does not determine the proof.
At total weight `h`, the two remaining words `0^h1^h` and `1^h0^h` have
opposite checksums.

The checksum's one-flip improvement has a different cause. Position `h` has
coefficient zero in `Z/hZ`, so it globally ignores the last tail bit.
THM-2160 proves that a composition-exact shell coloring cannot globally ignore
two fixed tail coordinates: its signed layer enumerator would have a double
root at `-1`, while the required endpoint polynomial has a simple root.

The faster construction escapes this no-go by **stratified obliviousness**.
Writing `h=2t` and relative tail `z`, it decides the exact-`n=h` stratum
`z_1=1` from only

```text
z_2+...+z_t mod 2,
```

then completes each Hamming layer on `z_1=0` by its exact signed binomial
defect. It stops at `3h/2` on the dyadic boundary, and a weight-`h-1`
obstruction proves that this half-tail use is shell-optimal. Thus the sample
`00001` needs only flip six, not flips six and seven.

The reusable object is therefore an affine cyclic cocycle on fixed-composition
shells:

```text
Hamming shell --rotation--> Hamming shell
       |                         |
       S_h                       S_h+h/2.
```

Exchangeability makes all atoms in a shell equiprobable. Without that
constant-weight property, the same bijection need not preserve the target.
The stronger lesson is that a quotient may ignore many coordinates on one
stratum only if a complementary stratum retains enough information to pay
the exact defect.

## 2. GMC: the faithful generalization is incidence, not pairing

The additive DvdK carrier has a finite root set `Omega`, a transitive Galois
action, an equivariant weight

```text
w(alpha)=alpha^k/p'(alpha),
```

and a selected root packet. Unlike Bernoulli words of fixed composition, the
root weights are conjugate field elements, not equal scalar masses. There is
in general no output-reversing involution pairing opposite weights.

The smallest clean control is

```text
p(X)=X^3-2,                  k=0.
```

If the roots are `alpha,zeta alpha,zeta^2 alpha`, then

```text
1/p'(alpha)=alpha/6
```

and the three weights sum to zero around a three-cycle. No two are negatives
of one another, and an odd three-element root set has no fixed-point-free
involution. Any proposed universal pair cancellation therefore fails before
the analytic or t-adic packet-selection step.

The exact group-theoretic boundary is simple. For a transitive `G`-set
`G/H` and a sign character `chi:G->{+/-1}`, a `chi`-equivariant bicoloring
exists exactly when `H subset ker(chi)`: the color of `gH` must be
`chi(g)` times the base color, and this is well-defined precisely under that
stabilizer condition. When it exists, it is balanced. Dyadic coin rotations
satisfy the criterion on their relevant orbits. For `X^3-2`,
`G=S_3`, `H=S_2`, and the sign character is nontrivial on `H`, so even the
abstract coloring is obstructed.

The correct replacement is already formalized in
[`GMC2RootPacketAlgebra.lean`](../04-computation/lean/TournamentH7/TournamentH7/GMC2RootPacketAlgebra.lean):
`packetSum_eq_zero_of_fixed` and `galoisPacket_baseValue_eq_zero` sum all group
translates with uniform incidence. This handles cyclic cancellation of any
order rather than demanding pairs.

More importantly, GMC(2) no longer has a `DvdK1` formalization leaf. The
front-door theorem is
[`GMC2Main.gmc2`](../04-computation/lean/TournamentH7/TournamentH7/GMC2Main.lean),
which invokes the proved
[`GMC2DvdKOmegaWiring.singlePolyCrux_holds`](../04-computation/lean/TournamentH7/TournamentH7/GMC2DvdKOmegaWiring.lean).
Both are imported by
[`TournamentH7.lean`](../04-computation/lean/TournamentH7/TournamentH7.lean).
THM-2101's monodromy, transcendental, and Tate-packet wrappers remain useful
alternate formalization projects, but they do not block the root theorem.

The typed transfer ledger is:

```text
coin j-subsets -> GMC root packets:
map: no canonical map;
shared mechanism: a finite group acts on selected subsets;
coin predicate: equal Bernoulli mass of two checksum halves;
GMC predicate: a base-fixed packet sum must vanish;
lost coordinate: equivariant root weights and their full-orbit sum;
restoring sidecar: Galois action + Lagrange full-root identity;
decisive hostile test: X^3-2.
```

## 3. JC: deck invariance stops at the fraction field

For the nonmonic quartic of
[THM-2158](../01-canon/theorems/THM-2158-quartic-quadratic-deck-parity-and-exact-finite-pole-criterion.md),

```text
P=V^2 z^4+beta z^3+gamma z^2+delta z+epsilon,
```

the quadratic deck makes the canonical approximate root invariant, but only
as an element of `C(x)[z]`:

```text
H_0=V z^2+(beta/(2V))z+(4 gamma V^2-beta^2)/(8V^3).
```

The deck involution fixes the entire base fraction field. Hence invariance
cannot imply polynomial regularity: for any nonunit `V`, the fixed element
`1/V` is not in `C[x]`. The apparent two-coordinate valuation sidecar
actually compresses to the single effective divisor

```text
sum_pi max(0,
  3nu_pi(V)-nu_pi(4 gamma V^2-beta^2))[pi].
```

The cubic congruence automatically forces `V|beta` in the UFD `C[x]`, and
this divisor is invariant under polynomial fibre translations. THM-2180
goes further in the reduced Keller branch: the normalized hostile face
`(Z-1)^3(Z+3)` proves `V|beta` directly, leaving only

```text
V divides 4gamma-(beta/V)^2.
```

The generic hostile quartic is already sharp for this last pole:

```text
P=x^2 z^4+x z^3,
H_0=xz^2+z/2-1/(8x).
```

It has polynomial linear coefficient and is deck invariant, while the constant
coefficient retains a finite pole. Thus the coin's output-toggling character
has no analogue on the live JC obstruction: the pole lies in the deck's
trivial representation. THM-2181 compresses the exact square prefix and closes
the already monic depressed case, but a polynomial root with leading term
`Vz^2` still needs a terminal monicization or quadratic-member lemma. The
useful future operation must first lower or exclude the remaining pole divisor,
not apply the involution again.

The typed ledger is:

```text
coin checksum -> quartic deck:
map: character parity suggests even/odd decomposition only;
preserved predicate: fraction-field descent of deck-invariant expressions;
destroyed coordinate: finite-prime valuation and polynomial integrality;
restoring sidecar: the single translation-invariant pole divisor above;
decisive hostile test: x^2 z^4+x z^3.
```

## 4. LRC: stratified completion becomes phase-danger geometry

THM-2167's good-prime digit equation also cuts each unrestricted layer into
equal affine fibres of size `q^11`. That superficial resemblance to a Hamming
shell is real but insufficient: lonely-runner success is not exchangeable
inside one carry fibre. THM-2174 exhibits millions of primitive distinct rows
on one carry-owner path realizing every divisor-mask target for the original
row and its single deletions.

The faithful finite-denominator quotient is the **phase-danger clutter**

```text
Min_subset{{i:14|a v_i|_Q<Q}:a in (Z/QZ)^*}.
```

At `Q<=14` it collapses to the divisor mask. At `Q=25`, aligned two-digit
words matter and every separate zero-digit mask can agree while exact phase
success differs.

Here the coin's stratified completion suggests the right positive move, not a
false equivalence. A repeated radix state may ignore an arbitrarily long
middle block on the stratum where its two boundary `h`-digit windows agree;
that pump preserves every denominator dividing `q^h`. The complementary
unbounded stratum must pay an endpoint-current defect. THM-2162/2174 make
that payment exact:

```text
sum [H(Wr)-H(Wl)]=6W mu/7,
W<=K/(7mu)
```

for a positive-mass deletion core extended to a zero-safe row.

Thus the LRC analogue of the coin theorem is not a global involution. It is a
two-lane carrier:

```text
bounded phase depth -> aligned digit-window completion;
unbounded phase scale -> signed endpoint current.
```

The divisor masks and danger clutters are hypergraphs under deletion
containment. Forcing them into a tournament would lose the relevant arity;
THM-2168 instead passes to the generated character lattice and uses the
labelled blocker incidence graph. That transverse quotient eliminates the
all-independent and scalar-`4+3` three-blocker lanes; the truly rank-one
scalar `5+3` tail is exactly where the quotient loses dimension.

## 5. Research consequence

The puzzle contributes a useful diagnostic rather than a new cross-problem
theorem:

1. identify the orbit action;
2. retain the shell/charge that chooses the action element;
3. ask whether weights are constant, equivariant, or valuation-sensitive;
4. locate any ignored coordinate as the kernel of an explicit observable.

For GMC the answer is full orbit incidence, already formalized. For LRC it is
the danger clutter plus endpoint current. For the quartic Jacobian branch it
is the one translation-invariant finite-pole divisor.
