# Common-path cancellation, joint runs, and arithmetic clock scope

**Status: PROVED** L1--L4 in the stated languages. Finite coefficient and
literal-program certificates accompany them; no universal return closure is
claimed. This extends the inherited ordered-run substitution mechanism,
including **THM-3121**, `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`.

## L1. Two-complement common-frame substitution

For a legal Hamilton path P in T, let tau_P(T) preserve its consecutive arcs
and reverse every other arc. Let Q be a legal Hamilton path in tau_P(T).
Suppose a vertex v has the **same ordered predecessor and successor** in P
and Q, counting a missing predecessor/successor as an endpoint marker.
Replace v by a transitive block TT_m, m>=1. Expand its occurrence in both
paths into the increasing path through that block. Then both expanded paths
are legal, and

    tau_(Q[m]) tau_(P[m]) (T[v <- TT_m])
      = (tau_Q tau_P T)[v <- TT_m].                         (1)

The equality preserves every exterior name and the ordered clone addresses.
A uniform newborn attachment commutes with it too: use the old v-mask bit
on every clone, and the old bits on the exterior vertices.

**Proof.** After the first complement the block is the consecutive-edge
interval tournament, so its increasing path is still legal. Every cross-block
edge reverses except the incoming/outgoing path hinges. The second path uses
exactly those same hinges. Hence it is legal at the block boundary; away from
the block legality is inherited from Q. The second complement reverses the
same internal non-path edges again and cancels every cross-block exception.
Exterior arcs follow the original two-complement program. Uniform birth
relations commute with substitution entry by entry. This proves (1). No
assumption that the intermediate block remains a module is made. ∎

The binary symmetric-difference formula for two complements is inherited;
the useful extra condition is a **physical ordered frame** that makes the
substitution legal through the non-modular intermediate relation.

## L2. A two-sample all-size affine-return certificate

Assume v has exactly one in-neighbour or exactly one out-neighbour in T.
Use an L1 program, and choose the newborn's direction so that the same
one-neighbour condition holds in its final base target W. Define T_m,W_m by
transitive substitution at v. For any fixed integers a,b,

    H(W_m)=a H(T_m)+b for all m>=1

holds **if and only if it holds for m=1 and m=2**.

**Proof.** Collapse every maximal path segment in the clone block to a visit
to v. A Hamilton path can visit the block at most twice: each visit requires
an exterior entry (respectively exit), except possibly the first (last), and
there is only one exterior vertex capable of making that entry (exit).
If c_r counts the quotient walks with r visits, the number of assignments of
the m ordered clone labels to r nonempty increasing runs is r! S(m,r).
Consequently

    H(T_m)=c_1+(2^m-2)c_2,
    c_1=H(T),  c_2=(H(T_2)-H(T))/2.                        (2)

The same holds for W. Two evaluations determine both coefficients, proving
the equivalence. L1 proves actual realization of the program for every m.
If T and W are strong, every nonempty substitution remains strong: an old
directed cycle through v supplies paths between any two clones, and the
strong quotient supplies paths to and from exterior vertices. ∎

The six certificates in `lift_certificates.json` satisfy these hypotheses.
Their clone ports are respectively 0,0,7,6 for sources 1457,1459,1463,1465;
both programs for a given source use the same port. The 1463 port has one
in-neighbour; the other three have one out-neighbour. All use two complements
and one birth, without helpers, vertex deletions, or class-valued resets.

## L3. Joint response compiler and the complete retained golden observer

For k ordered path slots and endpoint constraints at named **exterior**
vertices, the same one-neighbour argument bounds the total number of module
runs by k+1. Empty slots may be allowed; that does not increase the bound.
For any fixed such query,

    F(m)=sum_(r=1..k+1) c_r r! S(m,r),                    (3)

where the nonnegative integer c_r counts actual collapsed path components,
including their order and endpoint incidences. Thus k+1 evaluations at
m=1,...,k+1 determine this entire query. Equivalently its recurrence divides
`product_(j=1..k+1)(X-j)` for m>=1. This is the inherited ordered-run method
with a new one-neighbour bound, not a formula for arbitrary forgotten joint
port tensors. Queries naming particular clone endpoints need their boundary
sidecars and are not covered by (3) as stated.
The determination is over the integers/rationals with realization guards:
the evaluation matrix has diagonal factorials. For k=2 the residual
differences must be divisible by 2 and 6. It is not an invertibility claim
over arbitrary small-characteristic coefficient rings.

Keep the original exterior root r=2. Write E for the number of Hamilton paths
ending there and S for the ordered two-slot response whose first path starts
there. The golden observer is `(H,X,Y)=(H,(S-E)/2,(S+E)/2)`.
Put `R_1=1, R_2=2^m-2, R_3=3^m-3*2^m+3`.
The following triples are the **exact coefficients in this R basis**:

| Source / offset | H coefficients | E coefficients | S coefficients |
|---|---|---|---|
| 1457 input | 1457,594,0 | 269,0,0 | 1457,594,0 |
| 1457 / +2 | 4373,1782,0 | 809,0,0 | 4373,1782,0 |
| 1457 / +4 | 4375,1782,0 | 811,0,0 | 4375,1782,0 |
| 1459 input | 1459,923,0 | 509,0,0 | 1805,2192,828 |
| 1459 / +2 | 4379,2769,0 | 1527,0,0 | 5419,6578,2484 |
| 1459 / +4 | 4381,2769,0 | 1501,0,0 | 5527,6684,2536 |
| 1463 input | 1463,620,0 | 97,42,0 | 4861,2066,0 |
| 1463 / +4 | 4393,1860,0 | 158,66,0 | 8524,3602,0 |
| 1465 input | 1465,538,0 | 4,0,0 | 14580,11552,2616 |
| 1465 / +4 | 4399,1614,0 | 115,0,0 | 35579,27292,6114 |

These finite coefficients are obtained by exact subset endpoint counts for
m=1,2,3, justified at all m by the run bound. Size m=4 is an independent
literal control. In particular, the two 1457 programs satisfy the all-size law

    (H,X,Y) -> (3H+c, 3X, 3Y+c), c in {2,4}.             (4)

For the 1459 / +2 program the law is

    (H,X,Y) -> (3H+2, 3X+2^m, 3Y+2^m),                 (5)
    2^m=(H+387)/923.

The last equality is a realization guard, not permission to divide arbitrary
observer data by 923. Over modulo 7^a, 923 is a unit, so this dyadic coordinate
is recoverable from the current H without adding a phase coordinate.

All five inherited Q4 row matrices act on these full observers. In the old
F49 chart `(a,b,g)=(4(H+2X-Y),Y,H-Y) mod7`, (4) becomes

    (a,b,g) -> (3a,3b+c,3g).                             (6)

This is a genuine affine observer transport on each stated physical family.
It does **not** prove that its output belongs to the input grammar for another
application. The paths and full target relation are still necessary data.

There is a literal inverse **on each certified image family**: delete its
designated newborn, undo Q[m], then undo P[m]. Both reverse path operations
are legal because a complement preserves its selected path. For (4) the
observer inverse is `((H'-c)/3,X'/3,(Y'-c)/3)`. Its three modulo-three guards
are necessary; recognition of the actual image relation and path addresses
is additionally required. Those arithmetic guards alone do not classify a
tournament or license a source reset. This gives a realization-preserving
inverse for six infinite image families, not closure of their union.

## L4. Clone clocks are split; retain the ancestral cubic phase separately

For any fixed finite path-slot query under a single transitive substitution,
the run expression is an integer linear combination of j^m. Modulo a prime p,
terms r! S(m,r) with r>=p vanish, and the remaining bases lie in F_p. Therefore
every such response has period dividing p-1 for m>=1. This conclusion does
not require a one-neighbour port. It does require transitive substitution and
a fixed context; an actual complement can destroy that context.

More generally modulo p^a the same finite exponential formula has period
dividing the unit-group exponent after m>=a: nonunit bases vanish and unit
bases have that exponent. This is an upper bound, not a claim of a minimal
period in all contexts.

For the table above, the exact integer-observer periods modulo 7^a are

    3*7^(a-1) for both sides of the 1457 and 1463 lifts;
    6*7^(a-1) for both sides of the 1459 and 1465 lifts.    (7)

**Proof of exactness.** H has a unit coefficient of 2^m in every row, forcing
the order `3*7^(a-1)`. In the second pair of families, S additionally has a
unit coefficient of 3^m; once the H mode returns this forces
`6*7^(a-1)`. The 3^m coefficient is zero in the first pair. The elementary
orders follow from `2^3-1=7` and `3^6-1=7*104`, followed by binomial lifting.
Changing between `(H,E,S)` and the golden coordinates only divides by 2,
which is a unit at these moduli. ∎

This distinguishes the earlier irreducible cubic clock of order171 in F343.
Its ancestral n=3 state is retained explicitly in every certificate, using
`F7[t]/(t^3+t^2+1)` and the inherited root `eta=4+4t+t^2`. The old n=4 states
are `(4,6,1)` and `(6,4,1)`; multiply by eta^(-1) to retain n=3. The new clone
operation does **not** tick that interval-growth clock. Normal-basis F49
coordinates and their Frobenius swaps are recorded separately from both
physical growth operations. Equal field sizes never substitute for a map.

## Remaining return problem

The two observed affine maps (4) would generate the desired ternary words
algebraically if their outputs returned to one admissible source grammar.
That implication remains unproved. In particular, a full observer return
does not reconstruct a path selector or an actual target tournament. The
new precise task is to classify which output families admit another L1/L2
program, or to find a finite set of larger joint interfaces closed under both
offsets. The bounded searches recorded in AUDIT are neither a classification
of all two-complement programs nor an all-size negative theorem.
