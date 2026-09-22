# Glued number lines: exact conjugacies, XOR, and signed clock budgets

**2026-09-21. Status: PROVED scoped results; CITED classical sparsity input;
FINITE-EXACT controls; independently audited. Collatz remains OPEN.**
No novelty claim. No new general theorem in this packet is claimed to be
Lean-formalized. The pasted blueprint's assertion of such a formalization
is corrected in the [blueprint audit](glued_xor_20260921_blueprint.md).

The user's ideas have three exact survivors: an even-scale version of
signed Collatz conjugacy that retains its valuation shell, a two-bit
description of a tournament with a directed square fixed, and a lattice
model with two three-cycles exchanged by sign. The decisive additional
result is a sharper distinction between the possible clock behaviors on
the positive and negative-parameter Collatz sheets.

## Inheritance and the working board

The closest proved mechanisms are the signed parameter/content laws in
[arithmetic braids II](arithmetic_braids2_20260917_signed_cycles.md), the
[previous Catalan gate](catalan_elliptic_20260921_catalan.md), and the
[bounded-discrepancy argument](collatz_guards_20260921_discrepancy.md).
The canonical hostiles are a fixed point that an alleged conjugacy must
preserve, a transitive tournament changed by vertex switching, and a
density-zero set whose reciprocal sum can still diverge. The corrected
near miss is applying an identity after forgetting its division rule,
sampling measure, or directed structure. The least-used sidecars are the
two-adic shell, equal-time noncoalescence, and the distinction between
three components and a single period-three orbit.

Anchor: signed Collatz maps and actual clock behavior. Niche: the complete
four-vertex XOR model. Wildcard: exact CRT transport to primitive lattice
pairs and a sixfold geometry. Relevant method cards are **type every
analogy**, **find the hidden second coordinate**, and **controlled
forgetting needs a sidecar** in the maintained research protocol.

| Live concept | Operation | Preserved predicate | Decisive missing coordinate |
|---|---|---|---|
| Signed Collatz | Negation / dilation | Forward arrows when guards match | Parameter, shell and parity branch |
| Clock discrepancy | Product/carry normalization | Exact orbit identity | Distinctness or repetition |
| Known cycle core | Commuting relabelling | Least period | Component count is not cycle length |
| Order-four tournament | Two diagonal flips | Chosen directed square | Vertex labels and the third path chord |
| Squarefree sieve | CRT digit splitting | Prime-square exclusion | Ordinary height and ring operations |
| Hexagonal rotation | Sixfold lattice action | gcd | Time average need not equal spatial density |

## 1. The proposed sign-and-offset pairs depend on the division rule

Write A_b(n)=3n+b for the raw map, O_b(n)=oddpart(3n+b) on all valid
integer inputs, and U_b for its restriction to odd inputs and odd b.
The latter always removes at least one factor of two.

**PROVED raw classification.** For H(n)=an+c with integer a!=0,

```text
H A_b=A_d H  iff  d=ab-2c.
```

Thus the user's requested opposite-sign pairs do have exact raw embeddings:

```text
3n-1 on positive integers --H(n)=3-4n--> 3n-2 on negative n=-1 mod4;
3n+1 on negative integers --H(n)=-3-4n--> 3n+2 on positive n=1 mod4.
```

These are conjugacies onto invariant sublattices, not onto the entire
opposite half-line. They fail after ordinary odd-part acceleration:
the first sends the source fixed point 1 to -1, but O_(-2)(-1)=-5.
Indeed even b makes O_b raw on every odd input, so its only possible
cycle is the fixed point -b/2 when that point is odd. The requested
opposite-sign target sectors have no fixed points or cycles at all.

**PROVED accelerated rigidity.** An integer-affine embedding intertwining
odd-accelerated maps on an entire odd ray must have

```text
c=0, a odd, d=ab.
```

Two exact valuation cells force away every translation. In particular
U_b(-n)=-U_(-b)(n) is a forward conjugacy, not arrow reversal.

Even scaling can be restored by retaining the shell. On v2(m)=1 define

```text
V_(2b,1)(m)=(3m+2b)/2^(v2(3m+2b)-1).
```

Then V_(2b,1)(2n)=2U_b(n). Combining this with negation gives genuine
opposite-sign correspondences, with corrected parameter signs:

```text
positive U_(-1) --n -> -2n--> negative V_(+2,1);
negative U_(+1) --n -> -2n--> positive V_(-2,1).
```

Alternatively a translation gives 3m-2 or 3m+2 in a **parity-swapped
shortcut map**, provided its other branch is translated too. The full
classification, first-image exceptions and controls are in the
[conjugacy lane](glued_xor_20260921_conjugacies.md).

## 2. The dual three-cycle families are conjugate copies, with distinct periods

There are three known positive U_(-1) cycles, of odd least periods 1,2,7.
Negation pairs them with the three known negative U_(+1) cycles. Completeness
at all periods remains open. Within either family, these are three
components of different sizes, not a period-three orbit.

On their combined 20 tagged states, the iteration permutation has cycle
type 1^2 2^2 7^2. A commuting bijection preserves least period. Its entire
centralizer has order

```text
2*(2*2^2)*(2*7^2)=1568=2^5*7^2.
```

Consequently this particular known core has **no commuting symmetry of
order three**. Sign reflection is a commuting involution matching three
pairs. A threefold action on another object may still be useful, but
cannot be imported as a period-preserving symmetry of this core.

## 3. A stronger clock dichotomy, and why the 3n-1 sheet is worth studying

For a positive odd orbit with parameter b=+/-1, put
K_L=sum of its first L actual halving exponents,
D_L=K_L-L log_2(3), and q_L=2^D_L. Then

```text
n_L q_L = n_0 product_(i<L)(1+b/(3n_i))
        = n_0+(b/3)sum_(i<L)q_i.
```

An elementary packing proof gives a useful independent result: on any
nonperiodic positive U_1 orbit, the distinct images are prime to six,
so rearranging their reciprocals gives
`n_L q_L <= C n_0 L^(1/9)`. This forces liminf D_L=-infinity and even
`D_L<=-(8/9)log_2 L+O(1)` infinitely often. In particular a divergent
orbit cannot have D_L tend to +infinity, contrary to the pasted blueprint.

A classical input strengthens this to a full limit. Garcia--Tal's
[1999 paper](https://matwbn.icm.edu.pl/ksiazki/aa/aa90/aa9033.pdf), equation (6)
with Proposition 1, supplies a power-saving bound on the number of
distinct orbit values in each interval. Dyadic shell summation then
proves summability of their reciprocals. The paper and hypotheses were
read directly and independently reviewed. The detailed note distinguishes
this quantitative input from mere Banach-density zero.

The resulting dichotomy is:

| Positive odd dynamics | Eventually periodic | Hypothetical infinite orbit |
|---|---|---|
| U_(+1) | D_L -> +infinity; sum q_L diverges | D_L -> -infinity; sum q_L converges |
| U_(-1) | D_L -> -infinity; sum q_L=3n_0 | D_L -> -infinity; sum q_L<3n_0 |

For U_(-1), summability and D_L -> -infinity already follow from the
additive identity and positivity, without a literature input. The strict
budget inequality in the infinite case uses reciprocal summability.
Its positive boundary term is
`lim n_L q_L=n_0 product_(i>=0)(1-1/(3n_i))`.

This answers the user's preference for the negative-parameter sheet with
an exact invariant: **budget exhaustion detects eventual periodicity**.
The known cycles exhaust the budget; proving that every positive start
does so remains open. This is not a proof that there are only three cycles.
Proofs, source qualifications and exact known-cycle budgets appear in the
[blueprint/discrepancy lane](glued_xor_20260921_blueprint.md).

## 4. The four-vertex XOR model is exact, with its surrounding cycle fixed

Fix 0->1->2->3->0. The two diagonal orientations are independent bits
(a,b). Their flips form C2^2. Rotating the vertex labels acts as
`(a,b)->(1-b,a)`, cycling through all four completions. Therefore the
four completions are all isomorphic and each has five Hamiltonian paths.

For the skew sign matrix S, the three perfect matchings give

```text
pf(S)=s01*s23-s02*s13+s03*s12.
```

With that directed square fixed, this reduces to `-s02*s13`: the sign
is exactly a diagonal XOR, while det(S)=pf(S)^2=1 discards it. Fixing only
the Hamiltonian path leaves **three** bits and includes all four
order-four tournament types. Vertex switching is another operation again;
it can change a transitive tournament with one Hamiltonian path into a
strong tournament with five.

The [tournament lane](glued_xor_20260921_tournaments.md) proves these
distinctions and exhausts all 64 labelled tournaments, all relabellings,
and all switching classes. It also exhibits the failure of the one-type
cycle-fixed phenomenon at five vertices.

## 5. There is a genuine two-triangle geometry behind the sieve, with a limit

The lattice rotation R(a,b)=(-b,a+b) preserves gcd and satisfies R^3=-I.
Every nonzero lattice orbit has six points, split by R^2 into two
three-cycles exchanged by sign. Primitive lattice points have limiting
density 6/pi^2. This gives the proposed shape a precise mathematical home.

It does not explain the density merely by counting six points. For the
hexagon max(|a|,|b|,|a+b|)<=N, remove the origin. Its total and primitive
counts are `6 sum r` and `6 sum phi(r)`. For the centered square they are
`8 sum r` and `8 sum phi(r)`. Hence the two models have the **same exact
primitive fraction for every N**, despite their different rotational
symmetries. The multiplicity cancels in the normalized fraction.

The actual sieve bridge is explicit: split each residue modulo p^2 into
two base-p digits and assemble them by CRT. This is a bijection
`Z/M^2 Z -> (Z/M Z)^2`, M a product of distinct primes, under which
`p^2|n` iff both coordinates are divisible by p. Both probabilities are
therefore `product_p(1-p^(-2))=1/zeta(2)=6/pi^2` in the justified limit.

Adjoining a sign also doubles both numerator and denominator. Moreover,
gcd is constant on each six-orbit, so its primitive-indicator time average
is 0 or 1, not 6/pi^2. A Collatz orbit interpretation still needs a map
preserving legal steps and the relevant sampling law.
Proofs and all finite/global boundaries:
[density lane](glued_xor_20260921_density.md).

## Reproduction, audit and remaining frontier

```text
python 04-computation/experiments/glued_xor_20260921_verify.py
```

The [master replay](../../04-computation/experiments/glued_xor_20260921_verify.py)
runs four normal/optimized pairs, verifies source hashes, pins the source
archive and proof notes, and records the exact universes in the
[manifest](glued_xor_20260921_manifest.json). The blueprint lane checks
120,097 exact controls over 18,264 actual bounded prefixes and audits the
inherited Lean inventory; it does not formalize the new limit theorems.
Conjugacy and tournament notes passed separate peer reviews; the
discrepancy argument and primary-source implication passed two reviews.
The density lane's maps and count proofs were independently checked by
the parent, with optimized-safe exhaustive controls.

The surviving research targets are concrete: universal budget exhaustion
for positive 3n-1; exclusion of infinite positive 3n+1 orbits despite their
necessarily negative discrepancy limit; and a genuine trajectory-preserving
map beyond raw affine or finite group analogies. None is discharged by
duplicating known cycles, taking a spatial density, or forgetting a guard.
