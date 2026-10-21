# Exact mixed old/new zonotope refinement

**PROVED**, 2026-10-21. Separate corollary to `README.md`, preserving that
packet's audited theorem, producer and receipts unchanged. This improves the
quadratic/mixed part of its general unfinished-old-source gate, not its four
native examples, whose old sources were already complete.

Use the same original unsigned unfinished charge multisets F_o,F_n. Let Z_o
and Z_n denote their signed-sum zonotopes (convex hulls of all full signed
sums). Set

    G(F_o,F_n)=max_(o in Z_o,n in Z_n) [Q(n)+beta(o,n)].

Each zonotope has a vertex list contained in its twelve vectors
2C_i-Vsum defined in the base packet. Hence G is exactly the maximum of at
most144 integer expressions. The lists can have duplicates or zeros.

This vertex reduction follows directly from the geometry, not from assuming
that a separate Q-maximizer also maximizes the mixed objective. A vertex of
a planar zonotope is exposed by a generic linear functional; every source
copy then chooses the sign of that functional on its vector. The resulting
signs form one of the twelve angular cuts. The objective is affine in o
for fixed n and convex in n for fixed o. Successive replacement by vertices
therefore shows its global maximum occurs on a pair of listed vertices.

## Dominance with arbitrary destinations, signs and nonarrival

Let z_o,u and z_n,u be any old/new contributions assigned to terminals, using
arbitrary subsets of the unfinished sources and arbitrary signs. Then

    sum_u [Q(z_n,u)+beta(z_o,u,z_n,u)] <= G(F_o,F_n).     (1)

Proof: choose independent symmetric signs sigma_u for the terminals and put
O=sum sigma_u z_o,u and N=sum sigma_u z_n,u, using the **same** sigma_u for
the old and new group at a given terminal. All cross-terminal terms have
expectation zero, so

    E[Q(N)+beta(O,N)]
       =sum_u [Q(z_n,u)+beta(z_o,u,z_n,u)].

Independently add each unused old source with a symmetric sign. This does
not change the expectation. Independently add each unused new source with
a symmetric sign. Adding v increases expected Q(N) by Q(v) and changes
expected beta(O,N) by zero. Thus the expected objective of a full signed
assignment is at least the left side of(1). At least one full signed
assignment attains at least its expectation and is bounded above by G.
All random signs are proof devices; no physical clock is changed or guessed.
This also handles empty groups and sources that never arrive.

Together with the base packet's linear support estimates, equation(1) gives

    Delta(d) <= Delta(h)+G(F_o,F_n)+L_o+L_n,
    for every d>=h.                                  (2)

All original guards, source multiplicities, actual past clocks and residual
accounting remain as in the base theorem. In particular the physical target
remains equivalent to R'<=M-J'; this corollary supplies no new pair credit.

## Comparison and a strict hostile to scalar replacement

Let D_o=max_o Q(o), D_n=max_n Q(n). Cauchy--Schwarz gives

    G <= D_n+2 sqrt(D_o D_n).

The refinement is exact integer arithmetic with at most144 evaluations,
rather than a rounded irrational bound. If F_o is empty, G=D_n. If F_n is
empty, G=0 even when D_o>0.

Strictness already occurs with one old charge(1,0) and one new charge(0,1):
D_o=D_n=1 but |beta|=1, so G=2 whereas D_n+2sqrt(D_oD_n)=3.
The lost coordinate in the scalar bound is relative angular support. One
must not replace G by D_n+max beta between separately chosen Q-maximizers;
the 144-pair joint maximum is the justified expression.

## Independent finite controls

`mixed_gate.py` uses the unchanged base charge algebra but separately
enumerates all independent sign choices for every pair of multisets of at
most two of the twelve charge vectors. It compares these full sign searches
against the144-vertex evaluator, checks the real Cauchy inequality in integer
squared form, and includes empty-old, empty-new and strict angular controls.
It also enumerates all absence/three-terminal/sign assignments for a mixed
four-source packet, including an initially empty terminal, and verifies both
(1) and the complete energy bound(2).

    python3 mixed_gate.py --output mixed_receipt.json
    python3 -O mixed_gate.py --verify mixed_receipt.json

This is a structural refinement, not a new claim of uniform native coverage.
The previously certified finite native examples have D_o=0 and are unchanged.
