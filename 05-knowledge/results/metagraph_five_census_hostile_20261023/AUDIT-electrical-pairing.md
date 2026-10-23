# Independent fixed-root census, six-set recursion, and padding audit

October 23, 2026. **PASS** for the explicit order-23 pair and the stated
all-order extension of this particular vertex-refinement obstruction.
Auditor: the concurrent electrical-pairing lane.

## Independent finite calculation

`audit_independent.py` imports no producer code and starts directly from
the two connection sets. It enumerates all **7,315** five-subsets containing
vertex zero. Root-preserving permutations of the other four positions give
a separately constructed table of **48** rooted types. Both complete
48-entry count records agree, including every type. Translation invariance
then makes the equality hold at every vertex.

For the decisive six-set separator, the independent code does not test
six-vertex degree sequences. It counts transitive subsets by their unique
source. For a candidate mask M and requested size k, use

    F(M,k) = sum_{v in M} F(M intersect N^+(v),k-1),
    F(M,0)=1.

Each transitive subset has exactly one source, and that source dominates
all remaining members, so the recurrence counts each exactly once. It
reproduces **1,173** and **1,265** transitive six-sets. Controls include
TT_2=253 and TT_3=23*C(11,2)=1,265, using the actual regular degrees.
All transitive counts through order five agree. Different TT_6 counts
prove nonisomorphism without any assumption about circulant automorphisms.

Normal generation and optimized replay passed 15 and **16 checks**
respectively over these complete finite universes:

    python -B audit_independent.py
    python -B -O audit_independent.py

I also ran the producer's independent Python `check.py`; it passed **27
checks** and reproduced receipt SHA256
`e07d815325e112874961e6cd6943baf122052accc06231569092ac3e8527cbac`.
That program separately verifies all rooted and unrooted sizes <=5 and
the fixed-root six-incidences 306 versus 330. The peer program's TT_6
method is distinct from that score-vector enumeration.

## All-round proof

The initial rooted <=5 record is constant within each graph by vertex
transitivity and equal between graphs by the exact census. Both graphs
are 11-regular. An ordinary vertex refinement that retains old colors and
out- and in-neighbor color multisets therefore sees eleven copies of the
same color on each side at every round. This proves perpetual equality,
not merely agreement for a finite number of iterations.

The scope excludes pair-rooted flags, shared-neighbor contractions,
path-marked ports, and other higher-order refinement procedures. The
computation does not refute the new actual-port five-flag observer.

## Common-prefix extension

For TT_k -> T_i, partition every rooted motif of size <=5 into its prefix
and core parts. All cross arcs are forced. A root in the core uses the
equal rooted core census; a root in the prefix uses the equal unrooted
core census and the same prefix root position. This proves equal initial
records at corresponding prefix vertices and all core vertices.

Refinement preserves that correspondence: each core vertex sees the same
prefix colors as incoming neighbors and eleven core neighbors in each
direction; each prefix vertex sees the corresponding earlier and later
prefix colors and all 23 core vertices. Induction applies at every round.

There are two independent reasons the padded graphs remain nonisomorphic.
Their unique nontrivial strongly connected components are the original
nonisomorphic cores. Also the six-set separator persists quantitatively:

    TT_6(TT_k -> T_i) = sum_{r=0}^6 C(k,6-r) TT_r(T_i).

The terms r<=5 agree, so the difference is always **92**, at every k>=0.
Thus a six-vertex motif continues to separate the padded pair at each
order n>=23, although the padded tournaments are not strongly connected
for k>0. No infinite strong-family claim is needed.

## Bounded search and isomorphism scope

The prime-order multiplier-orbit proof in the report is valid: order-p
translation groups are Sylow p-subgroups of the target automorphism group,
because its permutation representation embeds it in S_p. Conjugating
them aligns an isomorphism with the normalizer of the regular p-cycle,
which is affine. This justifies the bounded cyclic search's class labels.
The peer audit did not rerun every discovery row at p=13,17,19,31; those
remain bounded discovery outputs and are not dependencies of the explicit
order-23 theorem.

## Structural extension to every fixed motif cutoff

The added all-fixed-k theorem follows by comparing at least
`2^((p-1)/2)/(p-1)` prime-cyclic isomorphism classes with at most
`(C(p,k)+1)^t_k` full unrooted k-census vectors. The exponential-polynomial
inequality guarantees collisions for all sufficiently large odd primes.
The root-orbit factor in the conversion to singly rooted records is
essential: a rooted type corresponding to an Aut(H)-orbit O occurs
`|O| N_H(T)/p` times at each ambient root. The report retains this factor
and the lower-order extension coefficient `C(p-l,k-l)` explicitly.

The finite audit verifies that 1063 is prime and that the displayed exact
inequality holds there for k=5 and t_5=12. The general asymptotic proof
does not depend on this deliberately loose numerical example. Ordinary
vertex refinement remains blind by regularity and transitivity. All
prime-cyclic cores are strong because every nonzero step generates a
directed Hamiltonian cycle. The same prefix argument gives failures at
every sufficiently large integer order, without falsely asserting that
padding preserves strong connectivity.

This is a structural boundary for fixed singly rooted motif records and
ordinary refinement, not a statement about path-marked flags or higher
joint incidence. It identifies which additional kind of information the
all-size program must retain.
