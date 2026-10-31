# Exact million-prefix tournament Hamiltonian-path construction

Status: PROVED elementary transfer and finite target pruning; FINITE-EXACT
coverage/obstruction, independently checked seed coordinates and actual children.
No global spectrum theorem and no novelty claim. The finite conclusion is
that every positive odd H≤999999 other than 7 and 21 is attained. This packet
does not re-prove the global exclusion of 7 and 21.

## Objects and coordinates

Take an arbitrary nonempty tournament A, a marked external vertex x with a fixed
cut, and an arbitrarily long transitive right spine v1,...,vk. Every vertex of A
beats every spine vertex; earlier spine vertices beat later ones. The bits of a
word w prescribe x→vi iff wi=1. The resulting tournament is (A▷spine)+x.

For the current parent A and cut retain:

* h = H(A);
* F = H(A+x);
* t = number of Hamiltonian paths of A whose last vertex beats x;
* q = number of ordered pairs (P,Q) of disjoint directed paths partitioning A,
  where P may be empty, Q is nonempty, and x beats the first vertex of Q.

Paths are vertex words, not unordered endpoint sets. A path of length one is
allowed. Thus 0≤t≤h and q≥0. The scalar F alone does not compose.

## Exact singleton transfer (PROVED)

Appending one parent sink v gives

    bit1: (h,F,q,t) -> (h,F,       2q+h,0),
    bit0: (h,F,q,t) -> (h,F+q+h-t,2q,  h).

Proof of q: if v belongs to Q and Q has other vertices, delete its terminal v;
this contributes q. If v belongs to P, delete its terminal v, allowing the
resulting P to be empty; this contributes another q. The leftover case Q=(v)
contributes h precisely for bit1. The new t is respectively0 orh because every
parent Hamiltonian path ends at its sink v.

For bit1, v is a sink of the entire augmented tournament, so F is unchanged.
For bit0, paths ending at v correspond to paths of A+x not ending at x; there
are F−t. Every other path has v immediately before x (v loses to all of A),
and is P,v,x,Q. When Q is nonempty these contribute q, and when Q is empty
these contribute h. Thus F'=F−t+q+h. These are disjoint exhaustive cases.

The singleton proof is self-contained and does not require a block-matrix theorem.

## All-word finite pruning (PROVED)

Fix target ceiling N. All realizable states have F≥1 and0≤t≤h. Appending bit1
preserves F; appending bit0 increases F by q+h−t≥q. Also q never decreases.
Therefore once q>N, no subsequent bit0 can produce F≤N. Subsequent all-one
words produce no new value. Likewise once F>N no extension can return belowN.
This proves that the finite state traversal in seeded_words.py exhausts the
image belowN over **all finite words**, not merely words up to a chosen depth.
Repeated states can be identified because the displayed coordinates compose
exactly. Zero-q neutral loops are harmless under state deduplication.

## Exact finite results

All labelled nonempty seed tournaments of orders1 through5, with all2^n cuts,
produce273 distinct (h,F,q,t) states. seed_roots.json retains one actual
(n,edge-code,cut-mask) witness per state. Edge-code convention: enumerate i<j
lexicographically; bit1 means i→j, bit0 means j→i. Cut bit1 means x→i.

The complete all-word image below1,000,000 has499,980 distinct values. Its
allowed-odd misses (besides7 and21) are exactly:

    34869,85885,103061,116587,153451,245923,535405,594793,
    595543,609873,618723,767135,774073,782103,824119,
    826373,868367,949777.

The first two are supplied by exact ordinal products:

    34869 = 3*59*197;
    85885 = 5*89*193.

Each factor occurs directly in the grammar. Consequently every positive odd
integer≤103059 except7,21 has an actual tournament witness. This advances the
specific THM-4104 finite prefix80405 but does not change the conjecture's
status. Witnesses can be reconstructed from the finite state traversal by
retaining predecessor addresses (or a second target-specific traversal).

The first obstruction to this grammar plus ordinal products is
103061=7*14723, with14723 prime. Because7 is excluded globally, any tournament
with H=103061 must have a strong component with that entire value. Thus its
absence from this grammar cannot be repaired by products of smaller grammar
values. The prime misses594793,774073,949777 are further such obstructions.
No absence from the unrestricted tournament spectrum is asserted.

All seeds through order4 yield only46 states; their first allowed missing
value is105. A single transitive seed class alone misses13 at every order.
These hostiles prevent treating binary word length or scalar growth as a
coverage proof.

## Concrete old-frontier witnesses

The file target_witnesses.jsonl gives actual seeds and words for:

*80407: order17, seed(5,9,27), word00111110100;
*80689: order16, seed(5,40,23), word1111011000;
*103059: order17, seed(5,8,14), word00101010000.

check_target_witnesses.py constructs the full adjacency matrices and uses
independent subset endpoint DP, without using the transfer law. A witness for
prime80407 is necessarily strong, and80689=7*11527 is necessarily strong by
the known exclusion7 and SCC multiplicativity.

## Verification and boundaries

The independent audit enumerates literal permutations and every cut between
paths for all labelled seeds≤5. It imports no primary profile computation;
its273 states agree exactly. It also literally constructs1092 deterministic
seed-plus-word tournaments and verifies H by independent subset DP. Normal
and Python-O outputs are compared. No check relies on assert.

Complete commands are in `README.md`; `run_checks.py --work-dir WORK`
repeats the primary and independent checks in normal and optimized modes.

The support traversal is reproducible and every state has an inductive actual
constructor; it is not an arbitrary numeric closure. The independent audit
checks the transfer on a finite universe, while the proof supplies its general
scope. The prefix computation is bounded; no global interval tiling, all-prime
construction, asymptotic density statement, or universal finite seed bank is
proved. The grammar's explicit misses are a reason to enlarge the carrier,
not a disproof of Chan–Pak Conjecture6.2.

## Four-coordinate positive block action (PROVED elementary reduction)

Set g=h−t, and order the state as(F,q,g,h). For an arbitrary nonempty marked
right block B, write its coordinates as(H,E,p,qB,s,tB,Fb), where E counts all
ordered two-path covers including an empty side, p counts those with nonempty
first path ending at a vertex that beatsx, qB counts those with nonempty second
path starting at a vertex beaten byx, and s counts full paths starting at a
vertex beaten byx. Then direct ordinal-cut cover decomposition
reduces to

    F' = H F + p q + (H−s)g + (Fb−H)h,
    q' = E q + qB h,
    g' = (H−tB)h,
    h' = Hh.

Every coefficient is nonnegative:0≤s,tB≤H and Fb≥H by the insertion injection.
The four-coordinate right action is sufficient for constructor output; it is
not claimed to recover the full seven-coordinate marked-cover monoid element.
Here is a self-contained derivation. Write I_A=F_A−s_A−t_A and similarly
for B, counting paths with x internal. In a path of (A▷B)+x, either no B
vertex precedes any A vertex, or x separates a nonempty prefix of B from a
nonempty suffix of A. The first case contributes

    (s_A+I_A)H_B + t_A s_B + h_A(t_B+I_B),

according as x is before/internal to A, between the blocks, or internal
to/after B. The second case is P_A P_B x Q_A Q_B, with P_B and Q_A
nonempty, and contributes q_A p_B. Expanding F_A and F_B gives the
displayed F' formula. This partition is disjoint and exhaustive because
all edges between blocks point from A to B.

For q', if Q contains a vertex of A, split P and Q at the A/B boundary.
The A cover is counted by q_A and the unrestricted B cover by E_B.
If Q contains no A vertex, P contains a full path of A and the B cover
is counted by q_B. Hence q'=q_A E_B+h_A q_B. All full paths of A▷B
end in B, giving t'=h_A t_B and h'=h_A H_B, so g'=h_A(H_B−t_B).
Every term uses uniquely determined path splits. This proves all four
coordinates without an external formula.

For a directed triangle B, specify its cut by j=number of vertices beaten byx.
Cyclic symmetry means j suffices. Its profile is

    H=3, E=12, p=9−3j, qB=3j, s=j, tB=3−j,
    Fb=3 ifj=0or3, and Fb=5 ifj=1or2.

Consequently the four Cj actions are

    C0: (3F+9q+3g,      12q,    0,  3h),
    C1: (3F+6q+2g+2h,   12q+3h, h,  3h),
    C2: (3F+3q+g+2h,    12q+6h, 2h, 3h),
    C3: (3F,            12q+9h, 3h, 3h).

Here C3 names a triangle block with all three vertices beaten byx, not a new
unmarked tournament name. Each block_target_witnesses.jsonl code Cj uses cycle
0→1→2→0 and letsx beat vertices0,...,j−1 within that block.

### Expanded-language pruning warning

The q>N pruning theorem above is an exact image theorem for singleton words.
For arbitrary blocks a p=0 extension may multiplyF byH without usingq. Indeed
p=0 iffx beats every vertex ofB, and then the augmented old tournament precedes
B in an ordinal join. Thus the only possibilities after q>N are further SCC
products; a generalized search must retain these products or explicitly factor
them out. The frozen block targets are witnesses, not an exhaustive
expanded-language decision procedure. Its successes are independently checked
on actual tournaments and do not depend on pruning completeness.

## Million-prefix closure (FINITE-EXACT, verified)

Four explicit single-triangle-block repairs close the four prime/7-prime misses
in the singleton-word image. The other14 misses are products of directly
attained smaller values. These exact factorizations are frozen in
check_mod3_and_products.py. In particular782103=3*260701 and
618723=49*9*23*61, so no nonexistent H=7 factor is used.

The four actual repairs, in block_target_witnesses.jsonl, pass a separate
full-adjacency subset-DP implementation that uses no profile transfer formula:

    order16: H=103061;
    order20: H=594793;
    order19: H=774073;
    order20: H=949777.

Normal and Python-O runs byte-match. Therefore, combining the all-word support
computation, the14 exact factorizations and these four independent DP receipts,

> Every positive odd m≤999999 other than7 and21 is the Hamiltonian-path count
> of an explicitly constructible finite tournament.

This is a bounded computational construction theorem, not the unbounded
conjecture. All million-prefix witnesses are obtained by the precise seed,
spine, triangle-block and ordinal-product operations in this packet. The
large traversal stores support rather than a separate adjacency certificate
for every value, but its state transitions preserve actual constructors and
target-specific witness reconstruction is available in `construct_witness.py`.

## Global mod3 restriction on the block language (PROVED)

This restriction applies to an arbitrary initial four-state vector, and hence
to every actual seed; it is not a bounded census inference.

Modulo3, C0 andC3 annihilate every coordinate. C1 andC2 send a state to

    (a*g+2h, 0, b*h, 0),

with(a,b)=(2,1)or(1,2). In such a state any singleton mapsq,g,h tozero, while
possibly changingF. Any following Cj then annihilatesF as well. Also three
consecutive nonconstant triangle blocks annihilate all coordinates: the first
has the displayed form, the second has onlyF possibly nonzero, and the third
killsF. The all-zero residue state is preserved by every later block.

It follows that any output not divisible by3 in the singleton/triangle-block
language contains at most two triangle blocks; if there are two, they must be
consecutive, and both must have nonconstant cuts. Thus every non-3 prime output
has one of the templates

    seed ; singleton-word ; (C1orC2)^r ; singleton-word ; apex,
    r∈{0,1,2}.

Forr=0 the two singleton words combine. Triangle blocks do not form a freely
iterable prime-producing amplifier. This does NOT imply the reduced templates
miss infinitely many primes, or prevent a universal construction using them;
the two singleton words remain unbounded. It gives a necessary grammar
normal form for any global prime-lane proof using this alphabet.

The symbolic proof is independently hostile-tested on all81 input residues,
all constant triangle actions, all8 triples of nonconstant triangle actions,
and all8 pairs separated by one singleton. The proof extends one intervening
singleton to any nonempty singleton word because the resultingq,g,h remainzero.

## A global p-adic restriction for arbitrary block sequences (PROVED)

Letp be prime, and let an arbitrary sequence of nonempty marked ordinal blocks
act on any actual initial state. Call an occurrence p-singular ifp divides both
H(B) andE(B). Let r be the number of such occurrences so far. No restriction
is placed on the intervening blocks. Then

    p^r divides h;
    p^max(0,r−1) divides q andg;
    p^max(0,r−2) divides F.

Proof: atr=0 the claims are integrality. A non-singular step preserves these
bounds: h' is a multiple ofh, q' combinesq andh, g' is a multiple ofh, andF'
is an integral linear combination ofF,q,g,h. At a singular step H andE each
supply a further factorp. Thus h gains one level; both summands ofq' have
valuation at leastr_old; g' has valuation at leastr_old; and the summands ofF'
have valuation at leastmax(0,r_old−1), using the extraH onF and the previously
proved valuations ofq,g,h. This is exactly the stated bound withr_new=r_old+1.
Forr_old=0 or1, read negative lower bounds as the trivial integrality bound.
QED.

In particular, a prime output other thanp permits at most two p-singular block
occurrences, even if arbitrary non-singular blocks intervene. More generally,
an output withv_p(F)=e permits at moste+2 such occurrences. This is a global
restriction on elimination/construction sequences with the required retained
state, not a claim about arbitrary vertex-deletion sequences of unrestricted
tournaments.

For triangle blocks, H=3,E=12, so every occurrence is3-singular. Strong order4
block code4 hasH=5,E=32 and gcd(H,E)=1; it evades this specific obstruction for
everyp. Thus retaining the two-path-cover countE separates two seemingly
similar odd-count multiplier blocks. Avoiding the obstruction is necessary for
some iterable prime-producing designs, not sufficient to prove universality.
