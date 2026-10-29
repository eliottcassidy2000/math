# Beyond 76: odd-period ring charts and the minimal reflected double

2026-10-29. Source pin `88c969273f2b6bc4ea2a8c5d29bb36c558eb59d2`.
**PROVED** below means an elementary proof is supplied, not an external novelty
claim. **FINITE-EXACT** claims have the exhaustive standalone checker
`check_crt.py`; normal and `-O` receipts agree. No full-word census at period 55
is claimed.

## Inheritance and scope

- Closest inherited mechanism: the all-odd-period legal-word reader in
  `session74_eleven_mechanisms_20261027/session74-clocks/RESULT.md` and
  `THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md`.
- Canonical hostile: `lucas_prime_reflection_towers_20261028/COMPOSITE_76.md`:
  the F19 coordinate is not reflection closed, although F4 is.
- Corrected near miss: prime-degree permutation arguments do not apply to
  a composite carrier merely because the reflection fixed-point count fits.
- Least-used sidecar: the prime-power lift in the new
  `reciprocal_lucas_carriers_20261028/THEORY.md`, P2. Prime splitting and
  local Lucas valuation lifting are inherited there, not claimed anew here.
- Session93 `incoming/CONCURRENT_COMPARISON.md` supplies the 76-state
  three-observation repair. We test its exact scope beyond 76.

Anchor: compatible composite charts. Niche: nonreduced prime-power clocks.
Wildcard: replace a failing quotient with its reflected double rather than
restoring every state. Concept board: CRT ideals; word dihedral action;
future-observation equivalence; nilpotent clocks; source-preserving weights.
All source paths above except the theorem are below `05-knowledge/results/`.

## 1. PROVED: every odd-period ring has an explicit one- or two-chart normal form

Let O=Z[phi], phi²=phi+1, n positive odd, and

    a=F_(n-1)-1, b=F_n, alpha=a+b phi=phi^n-1,
    L=L_n, g=gcd(a,b), M=L/g².

Then

    g = 1 if 3 does not divide n, and 2 if 3 divides n;
    R_n=O/(alpha) ≅ Z/L                         (3 does not divide n),
    R_n             ≅ F4 × Z/(L/4)             (3 divides n).

For M>1 the explicit odd-coordinate root is

    r = -(a/g) (b/g)^(-1) mod M,
    u+v phi ↦ u+r v mod M.

In the second case append the coordinate u+v omega in F4, omega²=omega+1.
Thus complete CRT factor charts are obtained by factoring M into *prime
powers*, not merely distinct primes.

**Proof.** Multiplication by a+b phi has integer matrix
`[[a,b],[b,a+b]]` and determinant a²+ab-b²=-L. If an odd prime divides
both a and b, phi^n=1 in O/(p). Taking the determinant norm gives -1=1,
a contradiction. Thus g is a power of two. Fibonacci numbers modulo 2 have
period 3. Modulo 4 the pair has period 6, and when n is odd and divisible
by 3, b=2 mod4 and a=0 mod4. Therefore g=2 exactly in that case.
Write beta=alpha/g=A+B phi; gcd(A,B)=1 and norm(beta)=-M. Hence
`gcd(B,M)=gcd(B,A²)=1`. Substituting r=-A/B gives r²-r-1=0 modulo M.
The displayed evaluation is onto Z/M, kills beta, and both quotient and
target have M elements, so it is an isomorphism. For g=2, M is odd;
beta is a unit modulo 2, and the ideals (2),(beta) are comaximal.
CRT gives O/(2 beta)=O/(2)×O/(beta), with O/(2)=F4. The n=1 ring is
the one-element zero ring and is interpreted separately. QED.

An immediate clock warning: if 3|n, the additive group is C2×C_(L/2),
and +1 has order L/2, not L. Named additions +1 and +phi together generate
it. If 3 does not divide n, +1 alone generates the additive group.

Useful exact factorizations and chosen roots:

| n | full ring | roots of phi in the odd factors |
|---|---|---|
|9|F4 × F19|5 mod19|
|15|F4 × F11 × F31|4 mod11; 19 mod31|
|21|F4 × F29 × F211|24 mod29; 179 mod211|
|25|F11 × F101 × F151|4 mod11; 79 mod101; 124 mod151|

The whole odd-coordinate roots for n=15,21,25 are respectively 81 mod341,
1445 mod6119, and 135722 mod167761. The checker tests the normal form
identities for all odd n through 201 as controls; the proof is all-order.

## 2. PROVED: a reflected double is the exact minimal dihedral carrier

Let W_n be cyclic binary words with no adjacent ones and let
`q(w)=sum_i w_i phi^i` in R_n. Fix **any ring quotient** pi:R_n→A. Define
C by moving position i to i+1 and R by moving i to -i. Write f=pi q and

    D(w) = (f(w), f(Rw)).

Then the image D(W_n), with its actual realization constraints, is the
coarsest future-observation quotient of W_n for observation f and the named
operation menu {C,R}. Its exact transitions are

    C:(x,y) ↦ (phi x, phi^(-1) y),
    R:(x,y) ↦ (y,x).

**Proof.** q(Cw)=phi q(w), and RC=C^(-1)R. Thus the formulas hold. Every
word in C,R equals C^j or C^jR, so D determines every future f-observation.
Conversely the empty operation and R already observe the two entries of D.
Two words with different D therefore cannot be identified by any exact
future-observation carrier. QED.

This works even for even n without asserting that q is injective there.
For odd n, q is the inherited bijection, so one can equivalently speak of
states in R_n and a transported word reflection rho. The theorem does not
pretend rho is a ring homomorphism.

The state bound is |A|², often much smaller than L_n, and the actual image
can be smaller still. Multiplicities may be retained on image states when
counting source words. Since C,R permute source fibres bijectively, fibre cardinalities
are constant along dihedral quotient orbits. Arbitrary source weights need
not be constant: their named source values must be preserved or transported
separately. A later named operation not
on the menu requires a new compatibility test; fibre counts alone do not
repair a failure.

### Why inversion is not an odd-characteristic algebra repair

For n≥3 the singleton words at positions 0,1,2 are legal. If transported
reflection on a nonzero odd-characteristic quotient were additive, these
words would force f(1)=1, f(phi)=phi^(-1), f(phi²)=phi^(-2). Additivity and
phi²=phi+1 would give phi^(-2)=phi^(-1)+1, whereas the inverse quadratic
relation gives phi^(-2)=1-phi^(-1). Hence 2 phi^(-1)=0, impossible when
2 is a unit. Reflection cannot be implemented by that additive shortcut.
In characteristic 2 the obstruction vanishes: for every odd n divisible
by 3, the F4 projection is genuinely closed under R by x↦x², because
omega^(-1)=omega² and Frobenius is additive.

## 3. FINITE-EXACT: arithmetic addition can destroy the double's compression

For odd n let T be **ring addition by 1**, transported back through the
legal-word bijection. It is not a bit toggle. Partition refinement starts
from observation f and repeatedly assigns each state the triple
`(current_class, class_after_R, class_after_T)`. The fixed partition is the
exact coarsest future-observation carrier for {R,T}; every finite future
word is covered at stabilization.

| n | observed modulus | minimal {C,R} states | minimal {T,R} states | refinement rounds |
|---|---:|---:|---:|---:|
|9|19|73|76|2|
|15|11|111|1364|6|
|15|31|751|1364|4|
|15|341|1361|1364|2|
|21|29|785|24476|12|
|21|211|16906|24476|7|
|21|6119|24473|24476|2|
|25|11|121|167761|34|
|25|101|10201|167761|18|
|25|151|22651|167761|18|

At n=25, modulus 11, the two odd-coordinate source states 135721 and
88564 both have D=(3,4). Adding 1 produces D=(4,3) and D=(4,5), respectively.
Thus even the provably optimal dihedral quotient cannot process T.
The source pin and coordinate formula fix the meaning of these states;
the receipt gives more hostile pairs for every proper chart.

The shortcut `(f(w), f(Rw), f(RTw))` reconstructs all states for the full
odd projection at n=9,15,21. But for n=25 modulo 11 it has only 1331 profiles,
with maximum fibre 654. It is not a universal three-test reconstruction
rule. The observed full {T,R} closure in the table is **not** asserted for
all odd periods or all prime powers. F4 is an explicit counterweight to
any slogan that reflection plus addition always forces faithful storage:
that projection remains a closed four-state carrier.

## 4. PROVED: translation closure reduces general repairs to subgroups

Here is an all-size organizing theorem, separate from the finite census.
Let G be any finite abelian group, pi:G→A a homomorphism, and rho:G→G any
permutation. Suppose the operation menu includes translations by an
additive generating set and rho. The coarsest exact pi-observation carrier
is G/H, where H is the largest subgroup contained in ker(pi) satisfying

    rho(x+h)-rho(x) in H for every x in G, h in H.

One exact descending algorithm starts H_0=ker(pi) and sets

    H_(j+1)={h in H_j : rho(x+h)-rho(x) in H_j for every x in G}.

Each H_j is a subgroup: addition follows by telescoping differences and
inverses by changing x. The sequence stabilizes, and its limit satisfies
the displayed condition. Every compatible subgroup survives every stage.
Finally, any observation equivalence preserved by all translations is a
coset equivalence: x~y iff 0~y-x; its zero class is a subgroup. This proves
both minimality and the algorithm. Positive translation words suffice,
since every translation has finite order and therefore its inverse is a
positive power.

For R_n with 3∤n, T alone supplies all translations. For 3|n, use both
+1 and +phi before invoking this subgroup theorem; the table's {T,R}
minimality was computed directly and does not smuggle in +phi. This is a
precise boundary between a ring-ideal quotient (initial chart), an additive
subgroup quotient (compatible carrier), and a general automaton partition.
They must not be identified by definition.

### Sparse certificate: primitivity at period 25 without a state census

The subgroup theorem makes the n=25 full-carrier claim much cheaper to
verify than to discover. Since the additive carrier is cyclic of order
N=167761=11·101·151, every translation-invariant nontrivial block system
is reduction modulo one of the six proper divisors d>1. The following
six pairs are witnesses that reflection preserves none of those systems:

| d | x | y | rho(x) | rho(y) |
|---:|---:|---:|---:|---:|
|11|1|12|1|158418|
|1111|0|1111|0|167186|
|1661|0|1661|0|4740|
|101|0|101|0|8039|
|151|2|153|103683|164330|
|15251|1|15252|1|118933|

In each row x=y modulo d but rho(x) is not rho(y) modulo d. The file
`sparse_certificate.json` supplies legal source words as 25-bit integers
for both states of each row. `verify_sparse_certificate.py` checks only
these twelve word occurrences, their reflections, the selected root,
the ideal norm and complete divisor list; it does **not** enumerate the
167761 states. The all-odd reader bijection is the stated inherited input.

Consequently **the group generated by arithmetic +1 and transported word
reflection on period-25 words is primitive**. This finite-specific theorem
has an elementary sparse proof, independently verifiable without the
original exhaustive search. Every nonconstant observation modulo a divisor
of N therefore has the full N-state future carrier. This strengthens the
individual chart census, but by itself does not identify the group as alternating
or symmetric and does not prove primitivity for arbitrary composite periods.
The companion `NINE_CARRIES_AND_ALTERNATING_GROUPS.md` supplies the further
classification step for periods 25 and 55.

This is a concrete example of the distinction the project seeks: an exact
future representation requires 167761 states, while a proof of that
requirement needs only six hostile pairs. The short certificate proves a
universal obstruction over all candidate block systems; it is not a
compressed replacement for the carrier itself.

## 5. PROVED/EXACT: the first 11-adic repeat changes the clock

At n=55,

    L_55=312119004989=11²·199·331·39161.

The selected 11-primary component is Z/121, with phi=37. It is not F121
and not two copies of F11. Direct arithmetic gives

    ord_121(37)=55, ord_11(4)=5,
    conjugate 1-37=84 mod121 has order 110.

The residue 11 is nonzero but its square is zero: this component is
nonreduced. Dropping its nilpotent layer shortens the local clock by 11.
The other conjugate cannot replace the selected root in a period-55 chart,
since its 55th power is -1 rather than 1. These are exact local calculations,
not enumeration of the 312-billion-state word carrier. In particular,
there is no reflection-closure conclusion at n=55 from this local test.

## Loss ledger and remaining questions

| Map | Preserved | Lost | Repair/test |
|---|---|---|---|
|R_n→odd coordinate|ring arithmetic, selected local clocks|F4 tag when n is divisible by 3|normal-form CRT coordinate; test R fibres|
|word→D|all future C,R observations|states indistinguishable by this menu|exact minimal; named T hostile above|
|Z/121→F11|residue arithmetic|11-step local lift in phase|retain full121 residue for exact return clock|
|G→G/H|specified pi and translation/rho futures|only behaviorally invisible differences|subgroup fixed-point theorem|

OPEN: classify all proper odd-primary projections whose {T,R} carrier is
faithful; prove or refute the three-test full-odd repair for all odd n
multiple of 3; derive double-image cardinalities without enumerating words.
These are narrower, testable questions than identifying every reflection
with algebraic conjugation or claiming a composite analogue of a prime-
degree group theorem.
