# Dyadic tail character versus GMC orbit incidence and JC pole descent

**Status:** proved mechanism comparison and scoped negative transfer audit.
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

The one-flip improvement has a different cause. Position `h` has coefficient
zero in `Z/hZ`, so the checksum ignores the last tail bit. That bit is read
only when it is itself the first change. Thus the sharper deadline is a
kernel-coordinate fact, not a stronger cancellation estimate.

The reusable object is therefore an affine cyclic cocycle on fixed-composition
shells:

```text
Hamming shell --rotation--> Hamming shell
       |                         |
       S_h                       S_h+h/2.
```

Exchangeability makes all atoms in a shell equiprobable. Without that
constant-weight property, the same bijection need not preserve the target.

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
`1/V` is not in `C[x]`. The exact missing sidecar is the valuation vector

```text
nu_pi(beta)-nu_pi(V),
nu_pi(4 gamma V^2-beta^2)-3nu_pi(V)
```

at every irreducible `pi|V`.

The hostile quartic is already sharp:

```text
P=x^2 z^4+x z^3,
H_0=xz^2+z/2-1/(8x).
```

It passes the first divisibility and is deck invariant, while the constant
coefficient retains a finite pole. Thus the coin's output-toggling character
has no analogue on the live JC obstruction: the pole lies in the deck's
trivial representation. A useful future operation must lower or exclude the
two valuation defects, not apply the involution again.

The typed ledger is:

```text
coin checksum -> quartic deck:
map: character parity suggests even/odd decomposition only;
preserved predicate: fraction-field descent of deck-invariant expressions;
destroyed coordinate: finite-prime valuation and polynomial integrality;
restoring sidecar: the two local divisibility defects above;
decisive hostile test: x^2 z^4+x z^3.
```

## 4. Research consequence

The puzzle contributes a useful diagnostic rather than a new cross-problem
theorem:

1. identify the orbit action;
2. retain the shell/charge that chooses the action element;
3. ask whether weights are constant, equivariant, or valuation-sensitive;
4. locate any ignored coordinate as the kernel of an explicit observable.

For GMC the answer is full orbit incidence, already formalized. For the
quartic Jacobian branch the ignored coordinate is not harmless: it is exactly
the finite-pole valuation sidecar still missing from descent.
