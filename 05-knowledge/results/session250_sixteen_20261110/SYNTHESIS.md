# From sixteen objects to generators that retain what matters

**Session250, 10 November 2026. COMPLETE mathematical deliverable.** New
proved-in-note mechanisms and exact implementation certificates; no
independent canon audit or mathematical-priority claim. **Chan–Pak remains
OPEN.** The original source of the user's five-block decomposition remains
unresolved. A precise candidate interpretation is constructed below, with
that limitation explicit.

## 1. The requested 1–5–4–5–1, reconstructed rather than misidentified

The inherited three-vertex digraph census has sixteen isomorphism classes
and arc-count polynomial

    D_3(z)=1+z+4z^2+4z^3+4z^4+z^5+z^6.

The five bins empty / sparse (one or two arcs) / half-full (three arcs) /
dense (four or five arcs) / complete have sizes exactly

    1, 5, 4, 5, 1.

This is a natural complement-symmetric reconstruction of the requested
pattern, **not confirmation that it was the historical referent**. The
fine arc grading is inherited from `digraph_metagraph_s20dm.out`; no source
for the exact historical grouping was recovered from the inspected
repository material or earlier conversation search. The bounded provenance
ledger records what was inspected, without an exhaustive absence claim.

All four middle types are self-complementary. Transpose fixes ten of the
sixteen types, complement fixes four, and their product fixes two. Their
Klein-four quotient has eight orbits. These are isomorphism-class fixed
points, not labelled graphs equal to their bitwise complements.

The genuinely useful information is the loss ledger. Transpose preserves
Hamiltonian-path count H; complement does not. The five bins also do not
determine H. True marked-arc classes have polynomial

    M_3(z)=z(1+z)^5,

with total32, whereas zD_3'(z) has total48. A marked class counts
automorphism-orbits of arcs, not raw arcs. The all-order twisted Burnside
formulas make the same distinction for every n.

[Candidate reconstruction, explicit sixteen rows and symmetry formulas](geometry/DIGRAPH_CANDIDATE_SIXTEEN.md).

## 2. A different sixteen-point geometry yields a new exact recurrence

The old Clebsch identification concerns the sixteen **cut switches of a
five-vertex tournament**, not sixteen digraph types. Its value is a genuine
map from the old even-path-forest expansion to vertex boundaries.

For any n-vertex tournament, the function s -> H(T^s) on the switch cube
has only even Fourier supports of sizes at most2 floor(n/3). Every surviving
path-forest component has even edge length, costs at least three vertices,
and contributes only its two endpoints to the switching character.

If A sums the n vertex-star switches, the uniform minimal annihilator is

    product_(j=0..floor(n/3)) (A-(n-4j)I).

Consequently switch-word sums have explicit rational generating functions.
At n=5 the entire eigenvalue−3 sector disappears:

    A H=H+30,
    sum_(k>=0) (A^k H)(T)z^k
      =(H(T)-15/2)/(1-z)+(15/2)/(1-5z).

There is also an exact inverse interface. The whole switching packet is a
quadratic signed cut energy

    H(T^s)=15/2+sum_(i<j)c_ij s_i s_j,

and the c_ij are explicit degree-two plus degree-four polynomials in the
ten tournament edge signs. Solving those sign equations is necessary and
sufficient for realizing a proposed packet. A Clebsch eigenfunction with
the right mean alone is not sufficient.

[Boundary proof and recurrence](geometry/SWITCHING_BOUNDARY_RECURRENCE.md);
[signed cut energy and exact inverse equations](geometry/ISING_INVERSE_INTERFACE.md).

## 3. Recover A000568 without losing switching stabilizers

A switching packet is a torsor K=F2^n/<1>. Permutations that preserve its
switching class form a group P and act affinely on K through a cocycle.
The ordinary tournament classes in the packet are exactly its affine
P-orbits. Weighted Burnside therefore recovers the full H-enumerator.

An even-order element of P has no affine fixed point, since that would
make it an automorphism of a tournament. An odd-order permutation with
c cycles has exactly2^(c-1) fixed switch states. Thus the number of
ordinary classes in the packet is

    (1/|P|) sum_(pi in P, odd order) 2^(c(pi)-1).

The unweighted count loses the cocycle; its H-weighted refinement needs
it. On a fixed affine coset, Fourier averaging keeps only endpoint sets
meeting **every permutation cycle evenly**. This intersects the new
boundary cutoff with the old odd-cycle Burnside mechanism.

The switching recurrence descends to ordinary classes, with loops and
switch multiplicities retained. At order five the global first-H walk GF
is90/(1-5z)+2/(1-z): its leading mode records twelve classes, while its
initial value records their H sum92. Assigning each class the labelled
packet mean would incorrectly give90.

This is an exact A000568-related alternative carrier, not a faster
all-order census or a recurrence in vertex number. Its missing practical
step is a birth rule updating both endpoint tensors and the affine
stabilizer cocycle.

[Affine Burnside and fixed-cycle filtering](orbits/AFFINE_BURNSIDE_BOUNDARY.md).

## 4. New explicit sequences from the inherited two-apex grammar

The positive36-state two-apex transfer and the value13 were already in
the repository. They are recovered here, not claimed as discoveries.
The new contribution is a closed first-moment formula for an arbitrary
fixed exceptional tournament Q, and an explicit certified second moment
for the two-apex case.

For k fixed apices, let asc_Q(pi) count forward adjacent steps of an apex
order. Define

    P_Q(E)=(1+E)^2/4^k
      *sum_pi (E+3)^asc_Q(pi)(E-1)^(k-1-asc_Q(pi)).

Summing H over all2^(kn) incidence words on a transitive spine of length n
has exponential generating function P_Q(exp(2^k x)). Thus its only
possible exponential bases are j2^k, 0<=j<=k+1.

For two apices, writing m_n=sum_w H(w) and s_n=sum_w H(w)^2,

    m_n=(0^n+3*4^n+3*8^n+12^n)/8,
    L^3(L-36)(L-24)^2(L-16)^2(L-12)^2(L-8)^4(L-4)^3 s=0.

The latter is not an unchecked prefix fit. After finding a candidate,
1144 consecutive residuals on the1144-state tensor carrier are verified
exactly; Cayley–Hamilton propagates them to every later index. The
certificate includes initial terms and the rational GF numerator.

For a uniform four-letter word, the relative variance converges to441/400.
Hence even exact first moments do not identify typical values. These words
are in bijection with tournament classes marked by an unordered deletion
pair whose complement is transitive. Forgetting that pair weights a class
by its automorphism-orbits of eligible pairs, not by a constant factor.

[Proofs, recurrences and finite witness receipts](sequences/PROOF.md).

### A new synthesis consequence: moment certificates themselves recur

Session249's finite alternating support bound now has a concrete recurrence
engine. For fixed j and odd L, set Y=(H-1)/2 and

    C_(n,j,L)=sum_(r=j)^(j+L)(-1)^(r-j) binom(r,j)
                sum_(|w|=n) binom(Y(w),r).

This is C-finite in n: each binomial is a polynomial in H and every fixed
raw moment is a matrix coefficient of a fixed tensor power. Moreover
C_(n,j,L)>0 certifies an actual word with H=2j+1, by the inherited
pointwise alternating bound. This is a rational certificate hierarchy,
not a target-uniform positivity theorem.

The first nontrivial example is especially instructive:

    C_(n,1,1)=(-s_n+6m_n-5*4^n)/4
      ~ -(841/102400)36^n.

It shares the displayed17th-order annihilator and eventually becomes
negative, even though H=3 persists at every spine length n>=1: use a
directed triangle on the two apices and one spine vertex, and place all
remaining vertices before it in a transitive prefix. Thus this fixed
moment truncation eventually fails to certify an existing value.
Higher moments or direct positive-language constructions are necessary;
the failure is not a missing H value. `run_all.py` checks the first cases.

## 5. The strongest constructive reduction: only strong deletion cores matter

For a rooted unit-offset seed, let B be the root deletion and O its
out-neighbours. The inherited tail mechanism has

    H(L o 0^t)=c+2^t M,     c=H(B)-T.

Here T is the selected-endpoint two-path-cover response. The new
seven-coordinate ordinal-gluing law retains H, the two-slot path count,
both endpoint counts and three rooted-cover responses. Cutting each
Hamiltonian path just before or after a selected vertex proves the two
boundary sums are each at least H(B).

This forces a sharp normal form. Every nonconstant unit-offset seed has

    B=TT_a ▷ C ▷ TT_b,    O nonempty and contained in C,

where C is one strongly connected component, itself a unit-offset
deletion core; a singleton is allowed. Conversely every such padding is
valid. Its slope is M_B=2^b M_C. **Deleting the padding preserves the odd
part of the ray slope.**

Thus joining separate marked components cannot generate new odd lines.
Search for genuinely new unit-offset lines can be confined to strongly
connected deletion cores. This is stronger than a finite search bound
and directly changes the constructive search space. It does not prove
that every required odd line has such a core.

[Exact gluing, core-normal-form theorem and examples](constructive/PROOF.md).

## 6. What to do next, now that the interfaces are explicit

1. Enumerate or construct **strong deletion cores** with unit offset and
   primitive odd slope, using endpoint responses rather than adding
   transitive padding. Retain finite below-ray coverage separately.
2. Seek a positive word substitution that updates the inherited two-apex
   boundary state and changes H by a useful arithmetic rule. Exact moments
   organize the family but do not replace this realization step.
3. At vertex birth, transport endpoint Fourier tensors **and** affine
   stabilizer cocycles. Their combined carrier is what the A000568/H
   bridge requires; neither scalar packet spectra nor coarse bins suffice.
4. Resolve the historical five-block source before claiming that the
   reconstructed digraph grouping is the user's original object.

This session uses the repository's connection-contract and quotient-loss
principles: exact operation first, lost coordinate second, useful transfer
third. No new general META-PATTERNS card is promoted from one session.

## 7. Reproduction and publication scope

Run `python3 run_all.py` here. It executes the new formula examples and
the exact moment certificate using only Python's standard library. The
checks include all64 three-vertex digraphs; all labelled tournaments
through five vertices for the switching and affine formulas; small
rooted gluing controls; and the positive word/moment examples. This is
verification of the new deliverable, not a renewed audit of inherited
research. No new canon ID is assigned. Hashes and replay output accompany
the final packet; publication is recorded in the session receipt.
