# An induced two-digit family from the old strong 25-core

**PROVED IN NOTES:** the all-module response identities, all-m induced
ternary children, and marked two-vertex extension below. **FINITE-EXACT:**
literal certificates. **OPEN:** repeatable rearming and Chan–Pak.

## Inheritance and the new choice of source

The old H25 bridge joined a singleton plus a strong six-core to a strong
seven-core; its endpoint was the M6² arithmetic reset. See
[H25-BRIDGE](../../strategy_fibres_absorption_20261113/fibres/H25-BRIDGE.md).
The new incoming
[CORE_REPLACEMENTS](../../prime_ternary_lift_20261113/CORE_REPLACEMENTS.md)
supplies both ternary digits from 15·2^m−1 by matching the full module response,
with a bounded change of the old core. That response method is inherited.

Here changing the source chart to the **old six-core itself** gives both
digits by induced extensions: no old arc is changed. The old bridge's failed
one- and two-port neutral inflations are not being reused as a return.

For a nonempty tournament M, let K_j(M) count ordered covers by j nonempty
directed paths. In a rooted frame A[p←M], compression of the maximal M runs
gives the exact identity H=Σ a_j K_j(M). The coefficient a_j counts legal
exterior words with j separated copies of the port. This is the inherited
run expansion; its literal quotient retains all uniform module incidences.

## One actual strong source, both children

Rows are outgoing masks; indices start at zero.

| Frame | Rows | Module port | Run vector (a₁,a₂,a₃,…) |
|---|---|---:|---|
| S | (44,17,2,6,13,30) | 2 | (25,12,0,…) |
| T₂ | (110,92,24,112,97,70,4) | 6 | (77,36,0,…) |
| T₄ | (2,8,11,97,15,87,151,63) | 2 | (79,36,0,…) |

The finite coefficient vectors have two independent controls: a complete
exterior-gap compiler, and literal Hamiltonian subset DP after substitution.
All three base frames are strong; their base H values also pass direct
permutation enumeration. Inflation of a vertex in a strong frame by a
nonempty tournament remains strong: exterior directed paths connect every
pair of vertices within the inflated set as well as the exterior vertices.

**Theorem.** For every nonempty M,

    H(T₂[M])=3H(S[M])+2H(M),
    H(T₄[M])=3H(S[M])+4H(M).

For every integer m≥1, taking M=TT_m therefore gives a strong source

    h_m=12·2^m+1

with **induced strong children** 3h_m+2 and 3h_m+4. The source has m+5
vertices, and its two children have m+6 and m+7 vertices.

**Proof.** The coefficient vectors differ by 3 times the source vector plus
2 or 4 times its first coordinate. K₁(M)=H(M). For transitive M,
K₁=1 and K₂=2^m−2. The exact embeddings below preserve every old vertex,
every old arc, and all module incidences. ∎

| m | Source | +2 child | +4 child |
|---:|---:|---:|---:|
| 1 | 25 | 77 | 79 |
| 2 | 49 | 149 | 151 |
| 3 | 97 | 293 | 295 |
| 4 | 193 | 581 | 583 |
| 5 | 385 | 1157 | 1159 |

These are actual infinite construction families, not an assertion of
infinitely many prime values or new spectrum records.

## The births in the original source labels

For +2, the embedding into T₂ is

    0→4, 1→2, module2→module6, 3→5, 4→3, 5→0.

The sole new vertex is target vertex 1. In source labels it beats
0,1,4 and all of M, and loses to 3 and 5. This one specified newborn cut
therefore gives 3h_m+2 for every m.

For +4, the embedding into T₄ is

    0→6, 1→3, module2→module2, 3→4, 4→5, 5→7.

The new vertices are target 0 and 1. Equivalently, append a sink a to the
whole old source, then append b with b→old vertex1; every other old vertex
and a point to b. This gives a completely specified two-birth constructor.

## The marked reason for +4, and its exact failure of automatic return

For a rooted tournament (D,u), define

    G(D,u)=Σ H(D−V(P)),

where the sum ranges over all nonempty directed paths P starting at u.
Equivalently, G counts ordered two-path covers whose first path starts at u
and whose second path may be empty. In the orbit coordinates G=E+O.

**Generic two-birth lemma.** Add a sink a, then b with b→u and every other
old vertex plus a pointing to b. The resulting tournament has H(D)+2G(D,u)
Hamiltonian paths.

**Proof.** If the path does not use b→u, it ends with a,b and contributes
H(D). If it does use that arc, either a immediately precedes b or a is the
last vertex. In both cases the remaining old vertices split into a path
starting at u and a complementary path, possibly empty. Each case contributes
G. These cases are disjoint and exhaustive. ∎

For S[M], keep old vertex1 as root. Its exact G response is

    G=27K₁(M)+12K₂(M)=H(S[M])+2H(M).

There is a structural two-run guard. The module's sole exterior out-neighbour
is root1. An unmarked Hamiltonian path can leave the module only through that
vertex, so it visits the module at most twice. In a G cover the root already
starts the first path; neither of its two paths can leave the module after
entering it. Thus the marked cover has at most two module runs as well. The
coefficient identities do not depend on cancellation of higher-run terms.

Thus the generic lemma independently gives the +4 identity. The complete
marked response is checked by an exterior-permutation/two-path-gap compiler,
independently of the unmarked run recurrence.

For TT_m the preparation is G−H=2. It is **not** automatically retained by
the chosen children: already at m=1, no vertex in either displayed child
has G−H=2. Their full lists are retained in the receipt. Scalar attainment
therefore does not turn the construction into an induction on its own.

Changing to a different witness is a real finite option: the strong cores

    H77: (22,124,88,113,96,69,1), root6,
    H79: (64,17,67,7,13,31,58), root0

have G−H=2. Their module response vectors differ from those of the displayed
children, so this observation does not assert a uniform repair of the entire
dyadic family. Further rearming work must preserve the full required response
vector, or explicitly change the parameter and source chart.

## The module-count boundary

For M=C3, the three counts are 147,447,453. The increments are 2H(M)=6
and 4H(M)=12, rather than 2 and 4. This is the same digit-scaling boundary as
the incoming 15·2^m−1 construction. Retaining a module does not allow arbitrary
recursive substitution while keeping the digit fixed.

Reproduce with

    python -B 05-knowledge/results/prime_terminal_20261113_1610/repairs/check.py

The all-m claim comes from the finite coefficient identity, not the finite
DP cutoff. Literal positive controls use m=1,…,8, through order15. The two
unmarked base response compilers and the marked compiler are exact finite
enumerations; no polynomial fitting is used.

The finite H79 rearm is now supplied explicitly in
[SOURCE_TRANSFER.md](SOURCE_TRANSFER.md). It enters a different source family
and is neutral only at this base module. That distinction is retained in the
composed construction.
