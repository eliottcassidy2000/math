# Physical-source barriers for the complete C2/C3 bilinear reductions

**Status: PROVED, architecture-specific.** These results challenge proposed ways of proving universal odd coverage. They do not refute coverage by the full C1/C2/C3 union and do not assert any new gap in the unrestricted tournament Hamilton-path spectrum.

## 1. Setup and inherited input

A physical source is a nonempty rooted tournament (V,c) with exactly one Hamilton path starting at c. Its profile is

    v=(D,E,U,d)=(H(V), e_c(V), R_c(V), H(V-c)),

where R counts ordered spanning two-path covers whose second path starts at c and whose first path may be empty; H(empty)=1. The singleton is o=(1,1,1,1).

The frozen Session197 proofs `bilinear-c2-inverse/PROOF.md` and `c3-two-source/PROOF.md` characterize the entire strong C2 and C3 scalar classes by physical source pairs. They also supply the exact source grammar

    O(D,E,U,d)=(D,0,2U,d),
    M_L(D,E,U,d)=(g(L)D-(L mod 2)E+LU,
                  D-[L=1]E, Q(L)D+B(L)U, D),

with L>=1 and the integer gadget constants of ROOTED.md. A source can be weak. Its deletion coordinate must not be discarded.

Our nearest proved mechanism is that grammar and its bounded exact inverse. The canonical hostile is replacing physical membership by a response cone. The less-used sidecar is the exact one-path slice D=1. The concept board is physical membership, bilinear image, dyadic padding, arithmetic content, and presentation-local descent.

## 2. The one-path slice is exactly dyadic

**Theorem 1.** The physical-source profiles with D=1 are precisely

    (1,1,1,1),  and  (1,0,2^k,1), k>=1.

**Proof.** Every nontransitive tournament contains a directed triangle. Each of the triangle's three Hamilton words injects into the Hamilton words of the whole tournament: insert each remaining named vertex deterministically before the first vertex it beats, or at the end if there is no such vertex. Deleting the named inserted vertices recovers the old word. Thus H=1 forces transitivity. A transitive tournament has one Hamilton word, and the marked source root must be its first vertex. Every other vertex may independently belong to either ordered cover path, each path being in transitive order. This gives U=2^(|V|-1). For |V|>=2, E=0 and d=1; the singleton is the exceptional displayed profile. Conversely the source-marked transitive tournament realizes every displayed state. QED.

**Corollary 1.1.** The physical source-profile set, and its subset of primitive profiles with gcd(D,E,U,d)=1, are not semilinear, and hence are not definable in Presburger arithmetic in these ordinary integer coordinates. The same holds after the invertible linear change to the positive coordinates in the C2 proof.

Indeed intersection with D=1,E=0,d=1 and projection to U would otherwise make the powers of two semilinear. Every profile in that slice is primitive. A one-dimensional semilinear set is eventually periodic. If p>=1 were an eventual period, choose 2^k>p beyond the threshold; then 2^k+p would have to be a power of two, although it lies strictly between 2^k and 2^(k+1).

This does not exclude a binary-digit automaton: powers of two do have a regular binary language. It does not imply that the bilinear scalar images are nonsemilinear or that they omit any infinite target class.

**Corollary 1.2.** No fixed finite collection of congruence tests can replace exact source membership, even when combined with the elementary physical cone and parity constraints. If m is the least common multiple of the tested moduli, choose 2^k>2m and put U=2^k+2m. The fake profile (1,0,U,1) agrees modulo every tested modulus with the actual profile (1,0,2^k,1), but U lies strictly between consecutive powers of two. Both satisfy D,d odd, U=E modulo 2, D-E>=d, U>=D+d-E, and U>=2d. The congruence obstruction applies to any predicates depending only on the finitely many residue profiles, not merely linear congruences.

## 3. A parity-respecting false C2 target and its actual C3 repair

Use the C2 form

    N(v,w)=(D-E+U)g-(D-E+d)s
             +(D+3U-d)L-(U+E)f,

where w=(g,s,L,f). The triangle source v=(3,1,3,1) is actual. The vector w=(1,0,6,1) satisfies every elementary cone and parity constraint in Corollary 1.2, but it is not physical by Theorem 1. Substitution gives

    N((3,1,3,1),(1,0,6,1)) = 5+11*6-4 = 67.

The complete physical C2 inverse excludes 67: its proved cutoff is U<=floor((67-1)/6)=11, so a finite exact source enumeration and all admissible N/M pairs are exhaustive. The supplied verifier independently implements that enumeration and checks the absence. Therefore source-cone substitution can produce a genuinely false C2 membership assertion, not merely an invalid graph recipe for an otherwise C2-realizable number.

Nevertheless 67 is not a counterexample to the C2/C3 union. The C3 B-template has the actual representation

    B((1,0,4,1),(1,0,2,1))=67.

Both inputs are source-marked transitive tournaments. Its left response is (11,4,30,4), so grafting the converse of the second source gives H=11+2*30-4=67.

The verifier constructs the seven-vertex tournament and checks H=67 and C=3 by full Hamilton-word enumeration. The response row used throughout verification is (11,4,30,4). This is an explicit scope separator: false in C2, true in C3, not an unrestricted gap.

## 4. Arithmetic content obstructs universal source-wise descent

Define c(v)=gcd(D,E,U,d), using all four coordinates.

**Theorem 2.** Every coordinate of any O/M_L descendant of v is divisible by c(v). For each of the two C2 forms and eight C3 forms, the scalar evaluated on descendants of v,w is divisible by c(v)c(w).

**Proof.** Each source move is an integer linear map on the four coordinates, so divisibility persists inductively. Every scalar form is bilinear with integer coefficients. Every monomial therefore has the asserted product as a divisor. QED.

All physical D are odd, so every source content is odd. This also follows internally from the grammar: g(L) and Q(L) are odd by their recurrences, while B(L) is even for L>=2. Induction then preserves D,d odd and U=E modulo 2. For M_L, its first coordinate is odd because L(U-E) is even modulo 2; its second and third coordinates have the same parity (check L=1 separately). O preserves these properties directly. Consequently, if c(v)c(w)>1 and a permitted source-pair form has value H, no combination of source-wise O/M_L moves and changes among the permitted forms can produce either 3H+2 or 3H+4. Both targets are nonzero modulo the odd divisor c(v)c(w). This remains true for arbitrary move depth and arbitrary gadget indices. It forbids a universal local child theorem for every pair presentation, not a theorem that first replaces the presentation with unrelated primitive sources.

For a concrete physical example, two M_2 moves from o produce

    v=(9,3,9,3).

With w=(1,0,2,1), the N-form has H=69. Every subsequent source-wise pair continuation has count divisible by three, so neither 209 nor 211 is obtainable in that continuation architecture. This is not a claim that those targets are absent from C1/C2/C3.

**Prime-target consequence.** A prime target represented by the C2/C3 source forms requires both sources to have content one. Divisibility alone permits a content equal to that prime, but the proved positive lower bounds exclude it: C2 has H>=3U_v U_w+1, C3 Type I has H>=2U_v U_w, nonboundary Type II has H>=3U_v U_w, and its singleton boundary has H>=6U. Since c(v)<=U_v, all these force H>c(v)c(w). The regular C2 exception 15 is not prime. Thus primitive-source filtering is a sound necessary test for prime targets; it is not sufficient for realization or coverage.

## 5. Why the deletion sidecar remains necessary

The inherited weak-source theorem gives actual profiles

    a=(9,0,18,1),  b=(9,0,18,3).

They share (D,E,U) but have different deletion counts. Against w=(1,0,2,1),

    N(a,w)=133,  N(b,w)=129.

Thus using gcd(D,E,U)=9 would falsely conclude that these outputs are divisible by nine. The correct four-coordinate contents are respectively one and three. Physical sources with these exact profiles are recovered and directly counted by the verifier; the quoted weak-source collision is supporting provenance, not a required arithmetic assumption.

## 6. Dyadic padding of a finite bank has zero-density output

**Theorem 3.** Fix m actual source profiles. Permit arbitrary independent sink padding O^k on either input, but no additional M_L moves. The number of distinct C2/C3 scalar values at most N obtained from these inputs is at most

    1 + 10 m^2 (1+floor(log_2 N))^2,   N>=1.

The extra one permits the exceptional C2 value 15. Thus this restricted construction has natural density zero and cannot cover all sufficiently large odd targets.

**Proof.** All admissible C2/C3 representations satisfy the common source cutoff U<=H/3 from the inherited positive bounds. Padding multiplies each seed's U by 2^k. Since every seed has U>=1, at most 1+floor(log_2 N) padding choices per seed can occur in an output H<=N. There are two C2 forms and eight C3 forms. Counting all ordered pairs and all ten forms bounds their distinct values by the displayed expression; rejecting inadmissible boundaries or identifying duplicate profiles can only reduce it. The bound is o(N). QED.

This barrier is only for a finite bank with O-padding. It does not apply to arbitrary M_L gadgets, to an increasing bank, or to the addition of the unrestricted C1 scalar image. Primitive restriction does not repair source recognition: by Corollary 1.1, even the primitive source set is nonsemilinear.

## 7. Reproduction and stopping boundary

Run `python3 verify.py` and `python3 -O verify.py` in this packet. The verifier checks the exact small source bank by an independent recurrence, the complete C2 negative at 67, actual C3 realization, parity and cone controls, content examples, and literal deletion-collision graphs. All checks use explicit exceptions.

The established barriers reject cone-only source selection, fixed finite congruence recognition, omission of the deletion sidecar, and universal continuation from arbitrary imprimitive source pairs. They leave open whether unrelated realizable source pairs cover every remaining odd target. Any successful all-target proof must discharge physical source membership, and must be allowed to change a trapped presentation.
