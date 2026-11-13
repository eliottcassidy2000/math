# An arbitrary-core seam compiler with nonnegative path-cover coefficients

**PROOF CANDIDATE, independent review requested.** This is a uniform full-H
composition rule for a fixed exterior and arbitrary middle tournament. Its
source is the direct splitting bijection in PROOF §3; path-cover substitution
is inherited from THM-1975 (`the-path-cover-polynomial-is-the-refined-compositional-invariant.md`, early verified formulation) and THM-3121
(`path-cover-walk-content-substitution-kernel.md`, explicit universal kernel),
as well as the newer marked interface calculus.
The new application makes every coefficient explicit from the exterior
Hamilton skeleton, including absent-edge repair obligations.

Let A,B be fixed nonempty disjoint tournaments. Fix a set E⊆A×B and reverse precisely
those cross arcs in the ordinal join A→M→B, leaving every arc between the
middle tournament M and the exterior unchanged. The middle M is arbitrary:
it need not be transitive, strong, or have a prescribed order. Call the result
G_E(M).

For a permutation π of V(A)∪V(B), call it admissible if every consecutive
step OTHER THAN an A→B step is an actual exterior arc. Thus a missing B→A
step cannot be repaired by the middle and excludes π. Every A→B gap is a
legal middle slot; there is also a leading legal slot iff π starts in B,
and a trailing legal slot iff π ends in A. Let s(π) be the number of legal
slots. Among A→B gaps, precisely those whose pair belongs to E are mandatory:
they must contain middle vertices to repair the missing arc. Let r(π) count
these gaps. The leading/trailing slots are optional.

Let C_k(M) count ordered covers of V(M) by k nonempty directed paths,
including the vertex set partition and internal path orders. Set C_0(empty)=1
and C_0(M)=0 when M is nonempty. Then

    H(G_E(M)) = Σ_{k≥0} d_k C_k(M),
    d_k = Σ_{admissible π} binom(s(π)−r(π), k−r(π)),       (1)

with out-of-range binomials zero. These are finite NONNEGATIVE integer
coefficients, depending only on A,B,E. In particular d_k=0 above |A|+1
and above |B|+1, so at most min(|A|,|B|)+1 path-cover coordinates are needed.

Proof. Delete all middle vertices from a Hamilton path of G_E(M), retaining
the order π of its exterior vertices. Each middle run is a nonempty path
of M. Its possible location is exactly a legal slot: entering M requires
an A vertex, leaving M requires a B vertex, with the stated endpoint
exceptions. Missing A→B arcs force occupancy. All other bad steps are
unrepairable. Conversely choose an admissible π, k of its legal slots
including its r mandatory slots, and an ordered k-path cover of M. Insert
the cover paths into the chosen slots in their linear order. This is an
inverse bijection. Counting the optional occupied slots proves (1).
The slot bound follows because every internal A→B slot consumes a distinct
A and B; leading/trailing exceptions cannot simultaneously exceed the
stated two bounds. No hypothesis on internal arcs of M was used. □

## Transitive specialization and exact growth obstruction

For M=TT_m, C_k(M)=k! S(m,k), because an ordered set partition has a unique
increasing path in each block. Therefore

    H(G_E(TT_m)) = Σ_k d_k k! S(m,k).                    (2)

Equivalently each skeleton contributes
Σ_(j=0)^r (−1)^j binom(r,j)(s−j)^m, by inclusion–exclusion of empty mandatory
slots. Formula (2) exposes positivity which that alternating power sum hides.
For each fixed k≥1, k!S(m,k) counts surjections [m]→[k]; if m≥k it is at least
k! k^(m−k), by assigning the first k elements bijectively then choosing freely.
Thus if d_k>0 for any k≥2 and m≥k, a fixed full target c enforces

    m ≤ k + log_k(c/(d_k k!)).

If all d_k for k≥2 vanish, every positive middle length gives the same
count d_1. This is an exact all-height dichotomy in this fixed exterior
language: either every additional length is neutral, or count growth is
unbounded with an explicit exponential witness. The zero-length boundary
is d_0, not d_1, and must be treated separately.

## Exact repair invoice for two exteriors

For two exteriors on the same named seam, let d and d' be their coefficient
vectors. Their full-count difference for every middle tournament is

    ΔH(M)=Σ_k (d'_k−d_k) C_k(M).                         (3)

Hence an exact +2 replacement is a concrete integer constraint on the
entire path-cover vector of M, not on H(M)=C_1(M) alone. A sufficient
coefficient identity yielding ΔH=2H(M) forces M to be transitive whenever
one demands unit +2 repair. More complicated signed coefficient vectors
can be tested by the same finite skeleton compiler without truncating
middle order. This turns the five-core search into a finite exterior
calculation and a still-open realization problem for positive path-cover
vectors; it does not solve that realization problem by formal subtraction.

The slots here are structural boundary slots, not arbitrary directed graph
edges. The theorem applies only while all A→M→B incidences remain uniform.
Allowing the seam to alternate around a port requires a richer state.

## Arbitrary middle blocks cannot hide the activated length cost

For any M of order n and every k, C_k(M)≥k!S(n,k): each ordered set partition
admits at least one directed Hamilton path in each block (Rédei). Since d_k≥0,

    H(G_E(M)) ≥ H(G_E(TT_n)).                            (4)

Thus the exponential size obstruction does not disappear by making the seam
nontransitive. The transitive middle minimizes H for EVERY fixed exterior
in this uniform-seam language, simultaneously and at every size. In fact it
is the unique minimizer up to isomorphism: d_1≥H(A)H(B)>0 (proved below),
so equality forces C_1(M)=H(M)=1, which is equivalent to M being transitive. This is a
comparison theorem, not a classification of all tournament substitution.

For the single-backarc exterior of PROOF §3, nonempty M gives the particularly
short exact invoice

    H(J(M)) = h_A h_B H(M) + L_a R_b (2H(M)+C_2(M)).     (5)

The avoiding sector traverses M once. The backarc sector splits M into two
possibly empty paths: the two one-nonempty choices give 2H(M), and the
ordered two-nonempty choices give C_2(M). Consequently

    H(J(M)) ≥ h_A h_B + 2^n L_a R_b ≥ 1+2^n.            (6)

The first bound is attained by TT_n. Since h_A h_B,L_aR_b are positive,
equality in that first bound forces H(M)=1, hence forces M transitive.
For fixed exterior and target c the sharp size budget is therefore

    n ≤ floor(log_2((c−h_A h_B)/(L_a R_b))).

Even arbitrary cyclic middle blocks cannot exceed the transitive budget.
To make an indefinitely renewable fixed-count feedback return, one must
change the exterior/feedback language or delete/replace old seam vertices;
changing only the middle's internal structure cannot work.


## Every genuine feedback set is visible, with a sharp universal budget

If E is nonempty, the growing case always occurs. Fix one reversed edge
b→a and Hamilton paths P_A,P_B. Split P_A just before a and P_B just after b.
The exterior skeleton

    A_left, B_left ending b, a starting A_right, B_right

is admissible and has exactly two legal slots, regardless of which of the
outer pieces are empty. Other missing A→B steps merely mark some of those
slots mandatory. This skeleton contributes one to d_2. Distinct selected
edges and path pairs give distinct skeletons (its sole B→A step recovers
the selected edge, and its A/B subsequences recover the paths). Therefore

    d_2 ≥ |E| H(A)H(B),                 d_1 ≥ H(A)H(B).

The latter inequality comes from concatenating each A/B path pair in the
original order. Thus no nonempty feedback set gives a constant-length
response. This avoids leaving the d_2 positivity condition unexplained.

There is a stronger OCF certificate with an especially simple size budget.
For n=|M|≥1 put

    W_odd(M)=Σ_{S⊆M, |S| odd} H(M[S])H(M\S),
    W_E=Σ_{(a,b)∈E} H(A−a)H(B−b).

The old packings wholly inside A,B,M contribute h_A h_B H(M). For every
(a,b)∈E and odd S, a path of M[S] gives a directed odd cycle
b→a→P→b. Complete it with arbitrary disjoint odd-cycle packings confined to
A−a, B−b, and M\S. These contributions are disjoint across all edges and S:
each has exactly one cross-module cycle, which recovers its selected edge,
its ordered middle path, and its middle vertex set. Therefore

    H(G_E(M)) ≥ h_A h_B H(M)+2 W_E W_odd(M)
               ≥ h_A h_B+2^n W_E ≥ 1+|E|2^n.           (7)

The second inequality uses Rédei positivity on every induced subtournament;
there are 2^(n−1) odd subsets. Thus for fixed exterior and activated count c,

    n ≤ floor(log_2((c−h_A h_B)/W_E)).

This holds for every internal middle structure and any nonempty feedback
set. It is sharp as a universal statement: singleton A,B with their unique
cross arc reversed and transitive M have H=1+2^n. There is no possibility
of preserving a fixed target through unbounded middle growth by adding
internal cyclic structure or several feedback arcs while this exterior
and its interface directions are fixed.
