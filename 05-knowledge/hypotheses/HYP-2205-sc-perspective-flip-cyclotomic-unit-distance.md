# HYP-2205: SC perspective flips are the right cyclotomic carrier for the unit-distance H-gap bridge

**Status:** OPEN synthesis, with THM-409 proved and S629 finite evidence.

**Claim.** The shared structure behind the unit-distance `H=7`/`H=21` echoes
and self-converse tournament symmetry is not raw scalar equality.  It is a
carrier-with-conjugation pattern:

1. the tournament side has forbidden evaluations `H=7` and `H=21`;
2. the unit-distance side has additive edge-carrier echoes at the Eisenstein
   angle, including exact `u(5)=7` and the lattice echo `n=11 -> 21`;
3. the Moser `n=21` row is an explicit retained-side-channel object
   `P_2^-` with `57` unit edges and a unit spine, not a tournament with
   `H=21`;
4. the self-converse side is controlled by the anti-automorphism coset and its
   canonical involution on rooted perspectives, not by a canonical involution
   on individual vertices.

Thus the right common object is a conjugation carrier whose useful invariants
live after one observer/rooted lift.

## Evidence

S629 adds `04-computation/sc_perspective_flip_cyclotomic_s629.py` and stores
`05-knowledge/results/sc_perspective_flip_cyclotomic_s629.out`.

The SC atlas verifies through `n=6`:

```text
n | U(n) | rooted P(n) | SC classes | SC rooted
1 | 1    | 1           | 1          | 1
2 | 1    | 2           | 1          | 2
3 | 2    | 4           | 2          | 4
4 | 4    | 12          | 2          | 8
5 | 12   | 48          | 8          | 32
6 | 56   | 296         | 12         | 60
```

It also proves the exact theorem behind the computation: THM-409 says
`Anti(T)` is a coset over `Aut(T)`, so edge reversal induces a canonical
involution on `Aut(T)`-vertex orbits.  At `n=6`, some anti-automorphisms have
vertex cycle type `(6,)`; therefore the older sentence "the anti-map is an
involution on vertices" is false in general.  The corrected sentence is:

```text
edge flip is an involution on rooted perspectives.
```

That is precisely the observer-coupled correction from HYP-2121 in a
self-converse setting.

## H=7, H=21, and Unit Distance

S627/HYP-2204 already separated scalar equality from carrier equality:

```text
u(5)=7 = 4 spine edges + 3 tile/bulk edges,
literal unit-tournament H at n=5 is 15, not 7.

Harborth n=11 lattice echo = 21 = 10 spine + 11 tile/bulk,
exact planar u(11)=23.
```

S629 keeps that separation but adds the cyclotomic label:

```text
Phi_3(2) = 7,
Phi_3(4) = 21,
3*Phi_3(2) = 21.
```

The Eisenstein/60-degree lattice makes these scalars visible as edge-count
echoes, while the tournament Hamiltonian-path invariant keeps them forbidden.
The equality is a warning that a quotient is collapsing side channels, not a
proof that the two carriers are the same.

## The n=21 Case

THM-408 gives a concrete non-scalar reading of `21`:

```text
P_2^- has 21 vertices, 57 unit edges, and a 20-edge unit spine.
P_2^+ has 22 vertices, 60 unit edges, and a 21-edge unit spine.
```

Each added Moser slab contributes:

```text
+8 vertices, +27 unit edges = +8 spine steps +19 bulk edges.
```

So beyond the first small Eisenstein echoes, the relevant recursion is a
side-channel-rich Moser slab/ear recursion.  The number `21` has moved from a
forbidden Hamiltonian-path count to a carrier size with `57` additive unit
edges.

This is why the initial "just place Eisenstein integers" picture should be
treated as a low-shell model, not as the whole recursive law.  The next exact
proof target is not raw `21`; it is classification of endpoint-compatible
extensions of retained `21`-cores.

## Proposed Tests

1. **SC perspective fingerprint test.** For self-converse classes beyond
   `n=6`, compute the rooted-perspective flip profile: fixed perspectives,
   transposed perspective pairs, anti-cycle types, and `sigma^2` automorphism
   order.  Test whether H-maximizers and near-H-gap classes have special
   perspective-flip profiles.

2. **Moser 21-core extension test.** For the exact `n=21`, `57`-edge cores,
   classify all endpoint-compatible ears by gain packet, direction support,
   unit-spine compatibility, and totally-unfaithful obstruction labels.

3. **Cyclotomic side-channel test.** Track when `Phi_3` values appear as
   edge-count echoes, spine/bulk splits, LRC `C=2n-1` composite burdens, and
   tournament forbidden-H guardrails.  The prediction is that coincidences
   become useful only when the retained side channels are named.

4. **Rooted-conjugation transfer test.** Replace raw complement merging by
   THM-409's perspective involution in scripts that currently use SC classes
   as scalar fixed points.  The predicted gain is fewer false symmetries and
   better detection of observer-coupled payload.

## Tournament Analysis

Vertices in S629 are proof routes and carrier obligations, not points:

- SC anti-coset theorem;
- rooted perspective flip atlas;
- unit-spine carrier section;
- raw `H=7`/`H=21` scalar match;
- Eisenstein edge-count echo;
- `n=21` Moser slab recursion.

The proof-route tournament is transitive under the chosen gauge:

```text
SC anti-coset theorem
> unit-spine carrier section
> rooted perspective flip atlas
> n=21 Moser slab recursion
> Eisenstein edge-count echo
> raw H=7/H=21 scalar match.
```

This is a useful negative result: raw scalar matching is the least reliable
route once side channels are retained.

## Guardrails

The `Cl_2(pi/3)` and `1.014` similarities remain inspirational rather than
proved equalities.  HYP-2183 already records the honest status: Sawin's
unit-distance exponent is a proven arithmetic-carrier amplification, while
the tournament-side entropy exponent has not been measured in a theorem.  In
this hypothesis, the rigorous shared object is the primitive-cube-root
carrier/conjugation pattern, not an asserted equality of exponents.

## Assumption Challenge

Alternate vertices considered: tournament vertices, rooted perspectives,
anti-automorphism coset elements, Moser slab layers, unit spines, unit-edge
packets, Eisenstein shells, literal H-evaluations, and proof obligations.

Chosen vertices: proof routes/carrier obligations.  This preserves the
predicate "does the quotient retain conjugation/perspective payload and
H-gap side channels?" and destroys literal point coordinates, continuous
embeddings, and raw Hamiltonian-path equality.

Challenged assumption: that a self-converse tournament supplies a canonical
vertex-level involution, or that a unit-distance scalar `7`/`21` means a
tournament with `H=7`/`H=21`.
