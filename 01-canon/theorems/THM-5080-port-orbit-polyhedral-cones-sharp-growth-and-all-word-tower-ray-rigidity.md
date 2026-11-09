---
id: THM-5080
title: "Ternary port orbits: exact minimal cones, sharp growth and all-word tower-ray rigidity"
status: "PROVED by exact rational certificates and all-word induction; independently audited. Orbit cofiniteness and binary/ternary Chan-Pak completeness remain OPEN."
session: Session236, 2026-11-09
source: 05-knowledge/results/session236_chanpak_offsets_orbits_20261109/README.md
proof:
  - 05-knowledge/results/session236_chanpak_offsets_orbits_20261109/orbit/PROOF.md
  - 05-knowledge/results/session236_chanpak_offsets_orbits_20261109/last_letter/PROOF.md
related: [THM-5070, HYP-9640]
---

# THM-5080: exact cones and all-word rigidity, not orbit coverage

Let M1,...,M15 be the frozen base-3 generators of THM-5070, let s0=(1,0,1),
T=M1=[[1,3,0],[0,0,1],[1,1,0]], and v=(3,1,2). Let R be the entire orbit
and R+ its nonempty-word part. Its first coordinate is I(Omega(Tournament),3),
not the binary Hamiltonian-path count. Numbering follows the inherited file.
Words below are chronological: (4,6) means M6 M4 s0.

## A. Exact minimal cones

The conical hull of R has four extreme rays:

    (1,0,1), (97,32,33), (25,17,17), (1,8,8),

attained by (), (4,6), (4,2), (4). Its inward facets are

    (-1,2,1), (-17,-824,849), (0,-1,1), (8,7,-8).

The conical hull of R+ has five extreme rays:

    (193,33,129), (292,65,162), (97,32,33), (25,17,17), (1,8,8),

attained by (4,6,1), (4,6,5), (4,6), (4,2), (4). Its inward facets are

    (-3039,6402,2909), (-1,2,1), (-17,-824,849), (0,-1,1), (768,1415,-1511).

Exact facet checks prove generator invariance (240 and 375 inequalities).
Initial states are contained; every extreme ray is attained. These two directions
prove equality of the cones for the entire infinite orbit by induction. They do
not characterize reachable integer points.

## B. Sharp expansion and short certificates

On R+, E/h>=33/193, E/h<=8, O/h<=8 and h/O<=97/33, all sharply. Every generator
satisfies h'>=h+3E, hence h'/h>=292/193 after the first step. Equality occurs at
s=(193,33,129) followed by T. Consequently every word of length r>=1 ending
at h satisfies r<=1+floor(log_(292/193)(h)). An attained h has a tournament
witness on at most 1+4[1+floor(log_(292/193)(h))] vertices.

This gives polynomially checkable logarithmic-length certificates, not a
polynomial-time inverse algorithm or a proof that every large h occurs.

## C. Arbitrary-word tower rigidity

A product P in the fifteen generators satisfies Pv=lambda v, lambda>0, iff
every letter is T. Then P=T^r and lambda=2^r. The singleton-seeded orbit
never hits that ray.

For the complete 57-matrix ternary alphabet of marked blocks of orders 1 through4,
the same conclusion holds with allowed letters T,T^2,T^3,T^4. This upgrades
THM-5070(5)'s bounded-word test to a theorem at every word length.

Proof: inspect the final letter. In the 15-alphabet every non-tower has no
nonnegative preimage of v with positive h. In the 57-alphabet there are four
tower powers, 52 separating-form certificates and one exceptional matrix

    [[7,18,12],[2,6,7],[5,7,6]],

whose positive-h nonnegative preimage ray is (45,1,5). All 57 matrices map the
positive-h nonnegative cone into h<=3O, excluding that exception. Backwards
induction now proves the claim. The exception really maps (45,1,5) to 131v;
the invariant cone is essential. Breaking port-covering odd cycles into odd
paths also proves h<=xO for every physical port state and real x>=1.

## D. Offset boundary

With L(s)=h+2E+O, T^3=I+v(1,2,1) gives

    T^(3k)s = s + ((8^k-1)/7)L(s)v.

If a fixed word maps infinitely many points of a line parallel to v into a
finite union of such lines, subtracting two images gives Pv parallel to v.
Thus it is a tower power in either alphabet. This excludes non-tower digit
branches only in that finite-parallel-line scheme. Different directions,
nonlinear sets, larger alphabets and fixed-port operations remain available.

## Scope and checks

The state (10,4,5) satisfies all exact nonempty cone inequalities and h=1 mod3,
but is unreachable because 10 is a permanent ternary gap. Exact conical geometry
therefore does not settle integer reachability. The 56 binary matrices and the
all-odd-except-7,21 conjecture remain open.

The independent audit physically reconstructs all 57 matrices from 1,105 marked
blocks on three independent port-state inputs, then verifies all 615 cone
inequalities, sharp growth and all last-letter certificates. Unsaturated closure
to h=3000 reproduces 53,732 states and ten inherited gaps. No coverage extension
is claimed. Complete proofs, scripts and audit reports are in the linked packet.
