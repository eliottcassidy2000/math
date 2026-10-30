# Independent audit: ramified Pisano attachment

**Verdict:** the formulas and the 16-lift result are correct, provided
“orientation descends” is restricted to differences with **nonzero reduction**.
A nonzero kernel difference reduces to zero, so a literal demand that its
orientation descend is impossible. Nor does an arbitrary invariant tournament
have a well-defined quotient tournament.

## Scope and independent evidence

Read the fresh source `session109_source/05-knowledge/results/unit_lift_clock_20261015a.md`.
Its Theorem U1 explicitly assumes the prime does not divide the quadratic
discriminant and explicitly excludes ramified primes. Its odd-order orbit
pairing and classical binomial order-lifting method are prior scope; the
ramified filtration and the attachment/coherence calculation below extend
that scope. No claim of discovering classical valuation lifting is justified.

The independent stdlib checker uses integral-coordinate permutations, builds
each ramified ideal as an actual image set, and enumerates all points and
orbits. It imports no producer code. Every test raises an exception on failure;
none uses `assert`, so Python optimization does not remove checks.

| k | points | nonzero g-orbits | sign bits |
|---|---:|---:|---:|
| 1 | 25 | 8 | 4 |
| 2 | 625 | 48 | 24 |
| 3 | 15625 | 248 | 124 |

Normal and `python3 -O` runs pass and produce byte-identical JSON.

## All-depth proof

Put pi=2phi−1, so pi²=5, and work in the ramified quadratic extension of
Z_5. Its integral ring is Z_5[phi]=Z_5[pi], since 2 is invertible, with
uniformizer pi and residue field F_5. The quotient R_k is modulo pi^(2k).
For x nonzero modulo pi^(2k), write s=v_pi(x) and r=2k−s.
The shell has 4·5^(r−1) points.

For g=phi⁴=2+3phi, g−1=1+3phi has norm −5, hence valuation 1.
If z has positive valuation t, binomial expansion shows

    v_pi((1+z)^5−1)=t+2.

Indeed 5z has valuation t+2; the other mixed terms have valuations at
least 2+2t, and z⁵ has valuation 5t>t+2. The least term is unique.
For e coprime to 5, the leading term ez similarly gives valuation t.
Consequently for each positive integer d,

    v_pi(g^d−1)=1+2v_5(d).

The condition g^d x=x is therefore 1+2v_5(d)≥r. Its least positive
solution is 5^floor(r/2). In particular every orbit has odd length.
A nonzero orbit cannot contain its negative: an odd-order permutation
of the orbit cannot take x to −x (or square that equality and invert 2
modulo the orbit length).

The shells r=2j−1 and r=2j each contribute 4·5^(j−1) orbits.
Summing gives 2(5^k−1) nonzero orbits, hence B_k=5^k−1 antipodal
pairs. Independently choosing one orbit from each pair gives precisely
2^B_k labelled Cayley tournaments invariant under the specified multiplier.
These are not isomorphism-class counts and not counts of all tournaments.

## Reduction and the exact meaning of attachment

For k≥2, a nonzero downstairs element of valuation s retains that valuation
in every upstairs lift. Its remaining depth grows from r to r+2, so its
orbit length multiplies by 5. The inverse image of a downstairs orbit of
length L has 25L points; each upstairs orbit has length 5L and maps onto
the downstairs orbit. Hence exactly five upstairs orbits lie above it.

The reduction kernel has 25 points and consists of the upstairs shells
r=1,2 plus zero. Those two shells supply 8 orbits, paired into 4 antipodal
pairs. Thus

    B_k=5B_(k−1)+4,  B_1=4.

For each downstairs antipodal pair, label its five upstairs pairs by their
projection. Their five signs are independent in an arbitrary invariant
tournament: 32 possibilities. Exactly two patterns are constant and thus
descend; exactly one matches a specified downstairs sign. All four kernel
signs remain free. Therefore:

* a fixed downstairs tournament has exactly 2⁴=16 coherent lifts;
* all coherent upstairs tournaments number 2^(B_(k−1)+4);
* arbitrary upstairs tournaments number 2^(5B_(k−1)+4);
* coherence imposes 4B_(k−1) independent binary constraints.

Equivalently an arbitrary five-bit pattern can be parametrized by one
reference sign and four relative signs. There is no canonical downstairs
orientation without specifying coherence or an additional selection rule.
The coherent fraction among all upstairs choices is 2^(−4B_(k−1)).
Reduction is not a loopless tournament homomorphism on within-fiber edges:
those edges collapse. The statement concerns only edges between fibers.

The checker explicitly constructs noncoherent invariant tournaments:
phi and 6phi upstairs modulo 25 reduce to the same nonzero phi modulo 5,
but receive opposite signs. Likewise phi and 26phi modulo 125 reduce to
phi modulo 25 with opposite signs. All antipodal sign constraints still
hold, so these are genuine counterexamples to unrestricted descent.

## Q² reversal

Here Q is multiplication by phi in the (1,phi) basis, so Q² is multiplication
by h=phi² and h²=g. For an element with remaining depth r let
L=5^floor(r/2), its odd g-orbit length. Then h^(2L)=1 modulo pi^r.
Also h≡−1 modulo pi, so h^L−1 is a unit. Factoring
(h^L−1)(h^L+1)=0 modulo pi^r gives h^L=−1 modulo pi^r.
Since h^L=h·g^((L−1)/2), multiplication by h takes the g-orbit of x
onto its negative orbit. Thus **every** g-invariant Cayley tournament is
reversed by Q²; no additional sign compatibility is required. Conversely,
Q² reversal implies g=Q⁴ invariance by applying reversal twice.
Q² cannot preserve a tournament as well as reverse it. These statements
hold for both coherent and incoherent choices and were checked at every
nonzero point for k≤3.

## Reproduction

```sh
python3 session109_pisano_symmetry/audit/check_ramified_attachment.py \
  session109_pisano_symmetry/audit/ramified_attachment_audit.json
python3 -O session109_pisano_symmetry/audit/check_ramified_attachment.py \
  /tmp/ramified_attachment_optimized.json
cmp session109_pisano_symmetry/audit/ramified_attachment_audit.json \
  /tmp/ramified_attachment_optimized.json
```

No publication or external mutation was performed by this audit.

## Structural addendum: ordinary lexicographic composition

A coherent lift is (abstractly) the ordinary lexicographic product

    T_k ≅ T_(k−1)[H_k],

where H_k is the Cayley tournament on the reduction kernel K, with
connection set S_k∩K. To prove this, choose any set-theoretic section
sigma:R_(k−1)→R_k. The bijection (u,z)↦sigma(u)+z identifies the vertex
set with R_(k−1)×K. If u≠v, coherence makes the orientation between
(sigma(u)+z) and (sigma(v)+w) depend only on v−u, hence on the edge
u→v in T_(k−1), independently of z,w. If u=v, the difference is w−z,
so the induced tournament is H_k in every fibre. This is exactly the
lexicographic product definition. An additive section need not exist and
is not required. Each fibre is a module (all outside vertices have uniform
orientation to the entire fibre); translation identifies its internal
orientation with the same H_k.

The additive bijection R_1→K, x↦5^(k−1)x, intertwines multiplication
by g. Thus H_k has exactly 16 labelled invariant Cayley possibilities.
Conversely any lower tournament and any such H_k define a unique coherent
connection set: pull back the nonzero lower connection set and adjoin
its specified kernel connection set. The product assertion therefore
expresses both necessary structure and sufficiency.

Starting from 16 possibilities at level 1, fully coherent towers through
level k number 16^k=2^(4k). Their top tournaments uniquely determine every
lower tournament by descent and every kernel restriction, so this also
counts the labelled top tournaments satisfying coherence at every stage.
The unrestricted top count is 2^(5^k−1). These are counts on fixed labelled
rings with the specified symmetry and reductions, not isomorphism counts;
possible isomorphisms between products do not invalidate these counts.
Iteration gives an abstract iterated product of k invariant 25-vertex
kernel tournaments.

The notation 2^(k(p−1)) versus 2^(p^k−1) is valid conditionally in an
analogous setting where the established kernel bit count is p−1 and
the unrestricted bit count is p^k−1. Here p=5; this audit alone does not
establish such hypotheses at arbitrary primes.

The checker now also verifies all 625·624=390000 ordered distinct vertex
pairs at k=2 against the lexicographic product adjacency rule, using the
actual independently enumerated upper orbit-sign tournament and a
coordinate section. Normal and optimized reruns pass identically. The
all-depth argument above applies to every sign choice, not only the
particular exhaustive adjacency example.
