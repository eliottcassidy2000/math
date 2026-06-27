# Derived Tournament Conditions for AP/Goddyn-Wong LRC14 Atoms

- Created: 2026-06-24T08:40:20Z
- Coordinator: codex
- Cycle: manual-user-request
- Web search: Tao derived multiplicative functions; Collatz cycles and Farey/continued fractions

This post sharpens the AP/Goddyn-Wong thread into a finite tournament
classification problem.  The theme is:

```text
AP is the base tight identity.
Goddyn-Wong is the first admitted derived coefficient.
Every other "same kind" candidate must survive the same derived gates.
```

For LRC14 the two visible atoms are:

```text
AP = {1,2,...,13}
GW = {1,2,...,11,13,24}.
```

Both are tight at denominator `14`.  The question is not only "why do these
work?" but "what exact finite signature would any other object of this kind
have to share?"

## Three Niche Seeds

1. Tao's derived multiplicative functions suggest treating GW as a formal
   derivative of the AP tiling identity, not as a random exceptional swap.
2. Collatz cycle lower bounds use continued fractions/Farey structure for
   `log_2(3)`; the LRC14 analogue is the Farey-neighbor pressure around
   `1/14`, especially the `3/41` resonance.
3. The Kuratowski/Wagner mediant example says obstruction ratios can be added
   by disjoint-union arithmetic.  A mediant obstruction is not automatically a
   new primitive obstruction; it may be a sum of old atoms.

Add the tournament-analysis seed:

```text
Any finite packet of proof obligations can be made into a tournament once a
pairwise observable, cutoff, and tie Hamiltonian path are declared.
```

The summit route is to prove that, for `n = 14`, the only achievable
tournament isomorphism classes satisfying the necessary conditions are the AP
class and the Goddyn-Wong class.

## Source-Backed Facts

Tao's post defines a derived multiplicative function as the formal derivative
of a family of multiplicative functions.  Equivalently, it is the coefficient
of a multiplicative function valued in a ring with a nilpotent infinitesimal.
The key rule is a Leibniz law over coprime products:

```text
f(nm) = f(n) F(m) + F(n) f(m),    gcd(n,m) = 1.
```

Source: <https://terrytao.wordpress.com/2014/09/24/derived-multiplicative-functions/>

The Collatz cycle literature uses continued fractions of `ln(3)/ln(2)` to
exclude short nontrivial cycles; this is the right analogy for using Farey
neighbors to exclude short denominator escapes around `1/14`.

Sources:

- <https://en.wikipedia.org/wiki/Collatz_conjecture#Cycles>
- <https://en.wikipedia.org/wiki/Farey_sequence>

The local repo already records the AP/GW anatomy:

- `q(S) = 14`: 13-covering and 14-avoiding.
- perfect or near-perfect one-hole tiling of `Z/14`.
- all `+/-` unit residues are present.
- phase gaps obey the Steinhaus `<= 3` gap condition.
- GW is the unique valid single doubling `12 -> 24`.
- the valid doubling window is Jacobsthal-gated.
- `3/41` is the first Farey-neighbor resonance of `1/14`.
- the structure splits as a `2`-adic doubling factor and a `7`-adic apex.

Local anchor: `07-reflections/anatomy-of-a-tight-runner-set.md`.

## Exact Setup

Fix `N = 14`.  For a speed set `S` with `|S| = 13`, write:

```text
rho(s) = s mod 14
U = (Z/14)^* = {1,3,5,9,11,13}
||x|| = distance from x to the nearest integer.
m_r(S) = #{s in S : rho(s) = r}.
```

A denominator-14 endpoint witness means:

```text
min_{s in S} ||s u / 14|| = 1/14    for each u in U.
```

An AP/GW-kind endpoint candidate is a primitive tight candidate which:

1. has 13 speeds,
2. has `q(S) = 14`,
3. has the denominator-14 unit witnesses above,
4. has no positive-measure safe interval already discharging it,
5. has no unresolved support-6 or nonunit Kuratowski carrier debt.

This definition intentionally does not assume the AP or GW residue pattern.
Those patterns are what the tournament tests should force.

## Necessary Conditions

### N1. Covering Quotient

The candidate must be 13-covering and 14-avoiding:

```text
for every d <= 13, some s in S has d | s,
for no s in S do we have 14 | s.
```

Residues alone cannot certify this.  The false look-alike
`{1,...,11,13,26}` has the AP residue pattern mod `14`, but it loses the
multiple of `12` and therefore drops to `q = 12`.

### N2. Unit Boundary Shell

The endpoint witness at all unit times forces all six unit residues to appear:

```text
{1,3,5,9,11,13} subset rho(S).
```

Reason: for each `u in U`, some residue must satisfy `ur = +/-1 mod 14`; as
`u` ranges over `U`, this asks for every unit residue.  AP and GW both have
this full unit shell.

### N3. The 7-Apex

The covering quotient forces a multiple of `7`, while 14-avoidance forbids a
multiple of `14`.  Thus some runner has residue:

```text
rho(s) = 7.
```

This is the unique residue whose double is `0 mod 14`; it is the preimage of
the observer's empty seat under the doubling map.

### N4. Residues Lie, Divisibility Pays

The even shell may be represented by residues or by lifted speeds.  GW is the
example:

```text
12 is absent mod 14,
24 is present and pays the divisibility debt for 12,
24 mod 14 = 10, so residue 10 is doubled.
```

Any future candidate must carry both ledgers:

```text
residue multiplicity ledger on Z/14
divisibility ledger for d <= 13
```

The endpoint tournament may see only the residue ledger; the covering quotient
will still audit the divisibility ledger.

### N5. Steinhaus Gap Rigidity

At the endpoint, the occupied residues must form a `<= 3` gap configuration on
the 14-clock after multiplicity and holes are recorded.  AP has the uniform
one-hole pattern.  GW has one doubled residue and one extra vacated residue,
but the local gaps still stay inside the Steinhaus bound.

Any candidate whose residue multiplicity creates a four-gap on the endpoint
clock is not of AP/GW kind.

### N6. Derived Leibniz Law

Let `F` be the AP residue/divisibility packet and let

```text
D = packet(S) - packet(AP).
```

For a candidate to be a "derived AP atom," `D` must behave like a first
nilpotent coefficient over coprime clock factors.  For quotient factors
`a,b | 14` with `gcd(a,b) = 1`, its projections must satisfy:

```text
D_ab(x,y) = D_a(x) F_b(y) + F_a(x) D_b(y).
```

AP has `D = 0`.  GW has one first-order defect:

```text
D = -delta_12 + delta_24
```

with the residue shadow:

```text
-delta_12 + delta_10.
```

This is the Tao-inspired condition.  A "second derivative" candidate would
need a nonzero mixed nilpotent coefficient.  For LRC14, the existing census
says double-doublings are never tight, so any such second coefficient must be
zero or discharged by a separate support-6 relation.

### N7. Jacobsthal Gate

For a single acceleration site `v in {1,...,13}`, define:

```text
D_v(AP) = AP \ {v} union {2v}
W_14(v) = [14 - v, 27 - 2v].
```

The site `v` is admissible only if every integer in `W_14(v)` is non-coprime
to `v`, equivalently if the whole window lies inside one gap between
consecutive integers coprime to `v`.

For `N = 14`, only:

```text
v = 12,
W_14(12) = [2,3],
```

passes.  Therefore the single-derived endpoint orbit is:

```text
AP,  AP with 12 -> 24 = GW.
```

### N8. No Double-Derivative Atom

If two acceleration sites fire, the candidate is no longer a first-derived
atom.  The necessary condition is:

```text
D_v D_w = 0 in the tight endpoint quotient.
```

Operationally: double-doublings must either improve the lonely margin, violate
the gap/Farey/Jacobsthal gates, or reveal a support-6 wall that is not an
AP/GW endpoint atom.

### N9. Farey-Neighbor Pressure

Let `p/q` be a Farey neighbor or low child of `1/14`.  A candidate must block
escape at that rational:

```text
min_{s in S} ||s p / q|| <= 1/14.
```

If the inequality is strict in the wrong direction for some low `q`, then the
candidate has a safer witness than the endpoint and is not tight at `1/14`.

The haunted resonance is:

```text
3/41,    det [[1,3],[14,41]] = -1.
```

This is the LRC14 cousin of the Collatz continued-fraction argument: a hidden
short cycle or hidden tighter witness can only survive by approximating the
endpoint slope unusually well, so the Farey neighbors become compulsory
checkpoints.

### N10. Kuratowski Mediant Primitivity

The user-supplied Kuratowski/Wagner arithmetic says:

```text
K5    contributes (v,e) = (5,10)
K3,3  contributes (v,e) = (6,9)
m K5 disjoint-union n K3,3 contributes (5m+6n, 10m+9n).
```

So the mediant ratio is an additive obstruction ratio.  The necessary
condition for a new AP/GW-kind atom is primitivity:

```text
after known AP/GW charges are removed, the obstruction carrier must not split
as a nontrivial disjoint-union mediant of old obstruction atoms.
```

If it does split, each component must be discharged independently.  A mediant
class is useful bookkeeping, but not automatically a new summit obstruction.

### N11. Measurable Boundary Owner

A tight endpoint row has Haar measure zero as an open safe interval, but it
still has an exact endpoint witness.  Therefore the proof object must retain:

```text
Borel/Baire code,
boundary owner,
carry address,
support-wall debt.
```

A quotient that remembers only Haar-positive interiors is too coarse.  A
quotient that remembers only category is also too coarse: the endpoint can be
null and meager and still proof-live.

### N12. Support-6 Debt Test

If a candidate enters the support-6 relation carrier, then the projected
coimage class must be one of:

```text
AP/GW endpoint class,
charged unit C27 petal class,
two-swap class already killed,
exactly counted support-6 relation class.
```

A nonunit `K3,3`-like carrier with no address tax is not an AP/GW-kind atom;
it is a separate obstruction debt.

## Tournament-Isomorphism Technique

A tournament functor for this repo should be declared as:

```text
Phi = (V(S), A_S(x,y), cutoff, tie_path).
```

where:

- `V(S)` is a canonical finite vertex set.
- `A_S(x,y)` is a pairwise observable.
- `cutoff` turns the observable into `x -> y`, `y -> x`, or tie.
- `tie_path` is a fixed Hamiltonian path used to orient exact ties.

The output is an unlabeled tournament isomorphism class:

```text
[T_Phi(S)].
```

The proof strategy is:

1. define enough tournament functors that no proof-relevant structure is
   silently erased,
2. enumerate all achievable classes for `N = 14`,
3. prove every AP/GW-kind tight endpoint candidate maps into the achievable
   subset,
4. prove the achievable subset has only AP and GW after quotienting by units,
   reflection, and the Jacobsthal gate.

## Four Exact Tournament Functors

### T1. Runner-Pressure Tournament

Vertices:

```text
V_press(S) = S.
```

For `s in S`, define the pressure profile:

```text
c_j(s) = #{u in U : 14 ||s u / 14|| = j},    j = 1,...,7
P(s) = (c_1,c_3,c_5,c_7,c_2,c_4,c_6, gcd(s,14), s).
```

Orient:

```text
s -> t iff P(s) is lexicographically larger than P(t).
```

If `P(s) = P(t)`, use the Hamiltonian tie path induced by residues:

```text
1,13,3,11,5,9,7,2,12,4,10,6,8,0
```

and then by the smaller actual speed.

AP/GW common fingerprint:

```text
six unit vertices with profile (2,2,2,0,0,0,0,...),
one 7-apex vertex with profile (0,0,0,6,0,0,0,...),
six even-shell vertices with profile (0,0,0,0,2,2,2,...).
```

This gives a `6 | 1 | 6` convex block structure in the tournament.  If a
candidate lacks this pressure shell, it is not an endpoint-14 AP/GW atom.

### T2. Residue-Gap Tournament

Vertices:

```text
V_gap(S) = Z/14.
```

For a residue `r`, define:

```text
G(r) = (
  m_r,
  m_{r-1} + m_{r+1},
  1_{r in U},
  1_{r = 7},
  distance from r to 0 on the 14-cycle
).
```

Orient `r -> s` when `G(r)` is lexicographically larger than `G(s)`, with
ties broken by the same residue Hamiltonian path as above.

AP has:

```text
m_0 = 0,
m_r = 1 for r != 0.
```

GW has:

```text
m_0 = 0,
m_12 = 0,
m_10 = 2,
m_r = 1 otherwise.
```

The achievable residue-gap subset for AP/GW-kind candidates should be exactly
the unit/reflection orbit of these two multiplicity patterns, subject to the
divisibility ledger in N4.

### T3. Acceleration-Gate Tournament

Vertices:

```text
V_gate = {1,2,...,13}.
```

For a site `v`, compute:

```text
W(v) = [14-v, 27-2v],
a(v) = 1 if every integer in W(v) is non-coprime to v, else 0,
l(v) = length of W(v),
b(v) = number of prime divisors of v that hit W(v).
```

Orient:

```text
v -> w iff (a(v), b(v), -l(v), v) is lexicographically larger than
          (a(w), b(w), -l(w), w).
```

For `N = 14`, the unique source must be:

```text
12.
```

The acceleration-gate tournament is the cleanest finite summit target.  Once
`12` is the only admissible source and double derivatives are forbidden, the
single-derived tight orbit has only AP and GW.

### T4. Node-Squared Local Tournament

This is the "square the tournament nodes" carrier.  Each outer runner becomes
an inner tournament of the original size.

Outer vertices:

```text
V_outer = S.
```

For each `x in S`, define an inner tournament and its full local profile
class:

```text
V_inner(x) = S.
```

For `y in S`, let:

```text
R_x(y) = Counter_{u in U} (
  14 ||x u / 14||,
  14 ||y u / 14||,
  14 ||(y-x) u / 14||
).
```

Orient inside the node:

```text
y -> z in T_inner(x) iff R_x(y) is lexicographically larger than R_x(z),
```

with ties broken by the runner-pressure tie path.

Then define the node-squared tournament as the lexicographic product:

```text
T_square(S) = T_press(S) [ T_inner(x) ]_{x in S}.
```

AP/GW-kind necessary condition:

```text
the multiset of full local profile classes is either the AP 6|1|6 symmetry
class or the first-fold class later selected by the Jacobsthal gate.
```

This is intentionally convoluted but useful.  The score sequence of the inner
tournaments is too coarse: AP, GW, and many single folds are all transitive
there.  The full relative-profile class keeps the missing/doubled residue
visible.  It detects whether a proposed endpoint atom is base-AP-like,
single-fold-like, or carrying multiple independent local defects.

## Achievable-Class Claim For N = 14

Define the AP/GW tournament-achievable set:

```text
A_14 = {
  ([T_press(AP)], [T_gap(AP)], [T_gate with no fired source], [T_square(AP)]),
  ([T_press(GW)], [T_gap(GW)], [T_gate with source 12 fired], [T_square(GW)])
}
```

closed under:

```text
unit multiplication mod 14,
reflection r -> -r,
renaming of equal-pressure runners.
```

The desired summit lemma is:

```text
Every primitive AP/GW-kind LRC14 endpoint candidate has tournament signature
in A_14.
```

If this lemma is true, then the remaining work is not philosophical.  It is a
finite enumeration and isomorphism-class exclusion:

```text
enumerate candidate residue/divisibility ledgers,
compute T_press, T_gap, T_gate, T_square,
discard all tournament classes outside A_14,
check the two classes inside A_14 are AP and GW.
```

The first scout script for this post is
`04-computation/lrc14_derived_tournament_atoms_codex_s145.py`.  Its initial
readout is:

```text
only site 12 is Jacobsthal-admissible;
AP and GW pass the endpoint/divisibility/unit/apex/gap filters;
the AP-residue false look-alike fails covering_ok;
coarse pressure shells are too weak;
inner score sequences are too weak;
node-squared relative-profile hashes plus T_gate are the useful split.
```

## Working Conjecture

The AP/GW common core is not "two lucky examples."  It is:

```text
one-hole Z/14 endpoint tiling
+ full unit boundary shell
+ 7-apex
+ first nilpotent doubling derivative
+ Jacobsthal-admitted source 12
+ Farey-neighbor rigidity
+ no primitive Kuratowski/support-6 debt
+ Borel/Baire endpoint owner retained
+ tournament signature in A_14.
```

Every other potential of their kind must either be a higher derivative of the
same AP identity or a different obstruction atom.  The current necessary
conditions say the higher derivatives vanish in LRC14.  So the best creative
attack is to turn that sentence into a tournament-isomorphism computation.
