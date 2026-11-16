# What a target wheel needs that a residue wheel does not supply

## A sufficient physical descent certificate

Take a finite set of states, each specifying extra readiness of a
tournament (for example, a rooted or cut profile), and a well-founded
rank. A wheel row supplies a congruence and finite exceptions, an
inverse taking the **actual target** N to a smaller ranked target q
with a declared state, and a forward construction on every witness
in that input state that gives the required output state and count N.
If these rows cover every allowed target-state pair outside checked boundary
cases, induction proves coverage. A Riesel divisor row in the
[7-free wheel](descent/RIESEL_DESCENT.md) is an exact instance with
q=N/p and ordinal sum for the forward construction. On unrestricted
odd targets, prime and 7p inputs require other rows.

For a forward rule N=mq+a and readiness condition q≡r mod K, the
target condition is N≡a+mr mod mK (unless one has already shown that
the ready predicate needs less information). A wheel modulo 120
for N alone does not generally authorize division followed by
another independent condition modulo 120.

## A scalar nonzero offset cannot come from an opaque substitution

Fix a tournament skeleton, uniformly replace k≥1 of its vertices
by copies of an arbitrary nonempty T (converses are allowed), and
leave all other vertices single. Uniform means every outside
vertex beats all or loses to all vertices of a given copy, as
dictated by the skeleton. Suppose the resulting Hamiltonian count
is aH(T)+b for **every** T, without choosing its port state.
Then for all T the output is H(skeleton) H(T)^k.
Consequently k=1, a=H(skeleton), and b=0.

Proof: compress each Hamiltonian path to its skeleton walk with
maximal same-block runs coalesced. A skeleton walk that visits
replacement block i in r_i≥1 runs, and every unreplaced vertex once,
contributes a product of ordered path-cover counts with exactly
r_i nonempty slots. When T=TT_t, each contribution is a product
of surjection counts from t labeled vertices to r_i labeled slots:
within a given slot all its vertices have a unique increasing
order. The quantities are nonnegative. For any walk with r_i≥2,
its contribution grows without bound as t grows (all other factors
are positive for t large). But the assumed scalar affine function
takes the constant value a+b on TT_t. Thus no possible walk can
revisit a replacement block. The remaining walks are precisely
Hamiltonian paths of the original skeleton, and each lifts in
H(T)^k ways. Finally H(T) takes unboundedly many values, e.g.
ordinal triangle powers, so equality with an affine function forces
the stated conclusion. There is no nonzero constant and no uniform
even scalar radix in this architecture. This is not a no-go for
selective contacts or for rules that prepare T using its target.

This is the same reason a count-only compression loses continuation:
passive transitive padding keeps H=1 but exposes any hidden multiple
visit. A useful affine row has to constrain additional information, or
alter the prepared representative rather than treat its ordinary
count as a complete interface.

## A real +2 seam whose entire support mod120 is known

For the [one-arc reverse bridge](../unit_shields_wheel_followup_20261116/bridge/ODD_PROJECTOR.md),
each lower anchor is H0=85Q−14 and upper H1=H0+2, with
2Q−3h=1. The strong right word before reversal has the unique form
x(10)^t, where x starts and ends in 1 and t≥1 is fixed by x.
Thus Q≡2 mod3; H0≡1 mod5. For t≥2, the last four bits are 1010,
so Q≡3 mod8 and H0≡1 mod8. CRT forces

    H0 ≡81 mod120,   H1 ≡83 mod120   (t≥2).

Only x=101 yields t=1: the defect identity and h≤(A+1)/2 for
an ending-one head bound its odd address A by 5, leaving 3 and 5
to check. Its pair 1941,1943 is 21,23 modulo120. Thus the primitive
projector really gives infinitely many **absolute** counts in
the unit class 83, with a certified +2 from the class 81.
But the anchor theorem for this seam gives only (2/3) ln X+
O(sqrt(ln X)) lower counts through X. It cannot fill even the single
class 83; a finite union of polynomially many-per-log tilted rays
cannot give coverage. The earlier fixed dyadic reverse bridge
illustrates the same distinction: for t≥3 its lower rays have
residues {69,21,45,93} and {117} modulo120, respectively;
the paired upper residues are obtained by adding 2.

The right next condition is not “can this residue be hit by an
unbounded family?” but “does selecting this residue place the exact
count at one endpoint for every target in a decidable input class?”
The new [length and amplitude wheel](PORTS.md) proves an analogous
separation inside the marked 24:13 family.
