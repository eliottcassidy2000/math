# The omitted internal gap restores a factor of three

2026-11-10. **PROVED explicit observer and inherited divisibility applied
to the proposed variable-gap lane. FINITE-EXACT controls separate below.**
The five-core sufficient covering statement and Chan–Pak remain OPEN.

## Why this lane was opened

The five-core bank is not the complete H<=9 parent language. In addition
to its two strong H9 cores, that language includes

    C3 -> TT_r -> C3,

with arbitrary internal gap r>=0, arbitrary internal apex cut, and arbitrary
exterior transitive flanks. Session251 tested the gapless sixth diagnostic
but explicitly did not collapse nonempty gaps into exterior padding. The
incoming Session252 associative seven-state law provides an exact way to
retain the missing middle word. We tried it as an additional constructive
carrier rather than as another fixed finite core label.

The outcome is a precise stopping rule for the coprime target lanes: every
apex child with r>=1 has count divisible by3, regardless of the cuts or the
exterior flanks. The empty gap is genuinely exceptional. This divisibility
is already implicit in the inherited ordered-block capacity formula,
[Session193 block-divisibility, equation(8)](../../session193_triangulated_proof_20261105/block_divisibility/PROOF.md).
The broader 3^(t-2) triangle-bank divisor and its adjacent-triangle boundary
are also inherited from
[Session228 general divisibility](../../session228_recursive_constructions_20261108/mobile_triangle/GENERAL_DIVISIBILITY.md).
We claim a new explicit specialization/diagnosis here, not a new general
divisibility theorem.

## Exact seven-coordinate formula for a nonempty gap

Let the apex beat k vertices of the first triangle and l vertices of the
last, where k,l in{0,1,2,3}. Triangle rotation makes the cardinalities
sufficient for these two local cuts. Let the r>=1 internal singleton bits
in dominance order form an integer M with exactly r bits, allowing leading
zeroes. Put

    Z=2^r,
    R=Z-1-reverse_r(M),
    alpha=M mod2, beta=R mod2,
    delta_j=1 if j in{1,2}, and0 otherwise.

R is the complemented reverse middle word. It cannot be discarded: two
words with the same g-value need not have the same response.

For the entire coloured core C3 -> TT_r -> C3, the positive seven-tuple
(h,d,L,i,o,s,e), used in the outer-flank formula, is

    d = 9,
    L = 144 Z,
    i = 36(3-l)Z + 36R + 9(3-k),
    o = 36kZ + 36M + 9l,
    s = 3(3-k),
    e = 3l,
    h = 9 + 6(delta_k+delta_l) + 9(g(M)-1)
          +9kR+3k beta+9(3-l)M+9k(3-l)Z+3alpha(3-l).      (1)

The actual child with exterior addresses A,B is therefore

    H = h+9(g(A)-1+g(B)-1)+L AB+iA+oB+s(A mod2)+e(B mod2).

Every coefficient and h in (1) is divisible by3. This is a complete exact
compiler for the variable-gap source grammar; it is not merely a residue
test on an intermediate statistic.

**Derivation.** A coloured triangle has seven-state coordinates

    (d,L,O+,O-,I+,I-,q)
       =(3,12,3k,-k,3(3-k),-(3-k),2delta_k).

A nonempty transitive middle word has

    (1,Z,M,-alpha,R,-beta,g(M)-1).

Multiply these two states and the last triangle state using the exact
ordinal-join law in Session252 `algebra/PROOF.md`, equation(A4), then set
h=d+q, i=I+,o=O+,s=-I-,e=-O-. This gives (1) term by term.

## Empty gap and the broader separator mechanism

For r=0 the middle is the identity, not a nonempty state with a zero negative
middle channel. The exact gapless values are

    d=9, L=144,
    i=36(3-l)+9(3-k), o=36k+9l,
    s=3(3-k), e=3l,
    h=9+6(delta_k+delta_l)+10k(3-l).                    (2)

For k=l=1, (2) gives h=41, not divisible by3. Applying the r>=1 formula
at r=0 would incorrectly erase the signed endpoint cross term. This is a
hostile boundary control, not a negligible convention.

In fact the separator need not be transitive. For EVERY nonempty coloured
tournament D, every apex child of

    TT_a -> C3 -> D -> C3 -> TT_b

has count divisible by3. In the positive matrix channel, a triangle's
matrix modulo3 has only a possible (1,3) entry; sandwiching any middle
upper-triangular matrix between two such matrices gives zero. In the
negative channel the middle diagonal entry of every nonempty D is zero
by the inherited alternating-path identity. The only potentially surviving
cross term between the two triangle endpoint channels is therefore killed.
All seven response coordinates of C3 -> D -> C3 vanish modulo3; arbitrary
exterior flanks preserve this. The direct ordered-block capacity formula
gives the same proof: only the ADJACENT triangle-pair term can lose both
triangle factors, and a nonempty separator removes that adjacency term.

## What survives after division by three

Equation(1) is an exact positive integer formula after division by3. It
does not supply a tournament whose Hamilton count equals that quotient.
When k=0 and A=0, or when l=3 and B=0, an end triangle is an actual SCC
factor of the child; deleting that ordinal factor gives a genuine quotient
construction. In particular, with no exterior flanks, k=0 or l=3 suffices.
The flank condition matters because an apex-to-flank return route can
destroy this factorization. For mixed boundary contact k>0,l<3, such a
quotient has not been proved.
This is the same distinction between a divisible count and a realizable
quotient that invalidates naive arithmetic transport elsewhere in the repo.

For Chan–Pak, an infinite prime/7p target lane coprime to3 cannot be reached
by the nonempty-gap construction at all. Multiples of3 already admit the
ordinary multiplication route from smaller realized counts (with the
known exceptional seed63 handled separately). Thus the variable gap does
not fill the critical missing obligation. A productive escape would have
to change the old ordinal incidences, add another apex, or use boundary
cores whose path-cover divisor behaves differently. These operations are
outside this exact source grammar and need their own realization sidecar.

`check_gap.py` tests all k,l and all middle words for r=1..4 against literal
weighted-path signatures, recounts selected actual flank children, and
tests the nonempty-separator residue theorem for all labelled D of orders
1..3 and every cut of D. It separately checks the gapless H41 hostile.
