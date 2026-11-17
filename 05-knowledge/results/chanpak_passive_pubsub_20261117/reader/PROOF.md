# Passive two-run readers: the cone, target steering, and exact cancellation

A passive reader inserts a nonempty tournament B into a homogeneous slot of a quotient tournament: all vertices of B have the same arcs to each exterior vertex. Write h=H(B) for Hamilton paths and p=P₂(B) for two labelled path covers, with empty paths allowed. We classify the quotients whose spanning paths can visit the slot in at most two runs, construct physical reader weights approaching every slope in the resulting cone, and exhibit exact all-height cancellation and a scalar collision with provably distinct continuations. These are not universal target placement.

## The sharp cone

Let I be the exterior vertices that beat the slot, and O the exterior vertices beaten by the slot. Assume the quotient is strong and the exterior is nonempty. Call o→i, for o in O and i in I, a return edge. There is always at least one return edge under these assumptions. A strong reader has at most two slot runs for every substitution **if and only if** the return edges form a nonempty star. If they do not, a pair of disjoint return edges gives a directed path

    b₁ → o₁ → i₁ → b₂ → o₂ → i₂ → b₃

after substituting at least three vertices. Extend this to a Hamilton path by successive vertex insertions; insertion in a tournament preserves the ordered subsequence of the existing path (prepend, append, or take the first orientation change among internal adjacencies). The three b vertices therefore remain in three different runs. Conversely every return to B requires crossing from O to I; two such returns would use disjoint edges. A nonempty bipartite edge set with no disjoint pair is a star.

Let L_i(I) count ordered covers of I by a path starting at i and an optional other path. Let R_o(O) analogously mark the terminal vertex o. The two-run skeleton count is

    α = sum_(o→i) L_i(I) R_o(O).

Indeed any such skeleton is an optional I prefix, slot run, nonempty O path ending o, nonempty I path starting i, slot run, optional O suffix. A one-run skeleton with no return edge gives H(I)H(O). If it uses a return edge, it has precisely one of the same two split patterns with a run omitted. In one pattern no new forward I-to-O edge can be missing; in the other it can be missing, only because one of its endpoints is the star hub. For instance, if the hub is in I, it is the first vertex of the second I path, so the jump from the first I path to the first O path is safe, while the forward jump from the second I path to the O suffix (if present) can fail. Let Z be the number of failed patterns. Then 0≤Z≤α; the coefficient of a single spanning B-path before subtracting double-counted empty slots is H(I)H(O)+2α−Z. Thus the exact read is

    α p + β h,   β = H(I)H(O)−Z.

For any tournament A and marked vertex i, L_i(A)≥H(A): cut a full path just before i and retain both pieces as a cover. The terminal version is the same. Hence α≥H(I)H(O)>0 and

    -α < β ≤ α.

In particular, no reader in this entire class can make a small output shift from an interface with |Δp|≫|Δh|: its absolute output shift is at least α max(0, |Δp|−|Δh|). This includes the current reflected +2 strands whose genuinely observed P₂ shift is large. This bound does not cover arbitrary three-run readers.

## Two full physical coefficient banks

Write u>0 in binary as b₀...b_(ℓ−1), leading bit 1, and tag the 1 positions in a TT_ℓ. Define the explicit even statistic

    U(u) = sum_(i>j) b_i b_j 2^(i-j-1) + sum_i b_i b_(i+1).

This counts optional ordered covers with both the end of the first path and start of the second path tagged (both paths then nonempty). If the first endpoint is later, the intervening positions are freely assigned to either path; if earlier, the endpoints must be adjacent. The weighted starting cover count is exactly u.

For the **negative-bank** reader take I=TT_ℓ and O a singleton q, with return edges q→tagged vertices. For the **positive-bank** reader take I to be that TT_ℓ with apex i beating the tags, O a singleton q, and the sole return edge q→i. Both quotients are strong: from the slot go to q, return through the tag at the first vertex in the negative bank or through apex i in the positive bank, then reach all of I. This gives the two literal reads

    Φ_-(u;B) = u p + (1−U(u)) h,
    Φ_+(u;B) = (u+1) p + (u+1−(u mod 2)−U(u)) h.

For the first formula the failed cross patterns are exactly the U(u) two-tag covers; for the second there is no failed pattern because the only missing forward I-O arc is the very return edge, whose endpoints are already used. The marked start at i counts a singleton path plus u covers through its tags. Its incoming block's ordinary H is u+1−(u mod 2)−U(u): an odd cycle through i uses a tagged earlier bit and untagged later bit, and contributes 2^(j−i−1)+[j=i+1]. Add these inversion weights to U(u) and sum over all later positions of each tagged bit. The sum is u−(u mod 2), since the adjacent corrections add one for each tagged bit other than the last.

Their difference is therefore **exactly** p+(u−(u mod 2))h. For fixed physical B a target D can choose this difference whenever (D−p)/h is a nonnegative even integer: choose u that even integer (or add one), with u=1 for the zero case. This is an exact difference between two actual counts, not a construction of a single tournament of count D.

Also, if v=rev_ℓ(u), then U(u)/u = v/2^ℓ + O(ℓ/2^ℓ): expand uv/2^ℓ, whose terms with later endpoint first are exactly the exponential sum in U; the remainder and adjacent terms are O(ℓ). This makes β/α asymptotically dense in [-1,0] and [0,1] in the two banks. Prescribing a fixed leading prefix of u does not remove this slope freedom, since the leading digits of v are the trailing digits of u. Keeping the error bound rather than asserting exact ratio equality is important.

## Target-selected compensated rails

For a fixed dyadic λ in [0,1), choose a bit width s and r with 0≤r<2^s and rev_s(r)=2^s λ. Give the reader the address

    u = 2^(s+t) a + r,    a>=1, t>=1.

The t blank places separate the high word and low suffix. In U, the cross-pair sum factors, with no adjacency correction, so

    U(u)=U(a)+U(r)+λ 2^(s+t) a.

For modules whose differences satisfy Δp=λ Δh, the negative-bank output shift is then exactly

    ΔΦ_- = (1−U(a)−U(r)+rλ) Δh,

independent of the absolute height t. For modules whose differences satisfy Δp=−(1−λ)Δh, the positive-bank output shift is exactly

    ΔΦ_+ = ((r+1)λ−(r mod 2)−U(a)−U(r)) Δh.

Thus a target gap may choose a finite motif a through its required U(a), while a target count on that motif's dyadic rail independently chooses t by a power-of-two test. This test is exact and includes the true physical coefficient; it is not a density claim. For a prescribed U>0, trim trailing zeros (this preserves U) and use the gap d between the first and last 1: 2^(d−1)≤U. Hence a finite motif search up to width 2+log₂ U is complete. Trailing zeros in a can be absorbed into t. For U=0 the canonical motif is the single bit 1.

The endpoint λ=1 has its own exact physical reader: choose u=2^ℓ−1, so all tags are marked, U=2^ℓ−2, and Φ_-=(2^ℓ−1)(p−h)+2h. This is the actual cyclic substitution C₃[TT_ℓ,B,1]. With strong eight-vertex sources in outgoing-mask notation

    (156,65,162,54,198,211,13,74)  : (h,p)=(527,7448),
    (30,28,208,36,200,23,43,107)   : (h,p)=(535,7456),

we have Δp=Δh=8. The all-height output pair is N,N+16 with N=6921·2^ℓ−5867. At ℓ=6 these are the actual primes 437077 and 437093 (trial division). For a genuine interior λ=1/2, take u=2^t+1 and

    (138,136,11,112,135,151,183,12) : (199,4120),
    (100,237,48,37,107,192,12,93)    : (247,4144).

These strong sources have Δp=24 and Δh=48. The negative reader yields 8041·2^(t−1)+4319 and that count plus 72, all t≥2. At t=6 the lower count 261631 is prime by trial division. None of these examples implies infinite primes.

## A hidden scalar collision that a second passive read always separates

Take the positive bank at u=2^t+1, t≥2. Here its exact (α,β)=(2^t+2,2^(t−1)+1) makes it read (2^(t−1)+1)(2p+h). Two other strong eight-vertex inputs have the following independently checked three-slot profiles:

    C=(216,29,81,52,32,135,186,30) : (h,p,P₃)=(387,6024,46899),
    D=(86,128,42,67,46,11,182,61)  : (419,6008,45603).

Both outputs have H=12435(2^(t−1)+1) at every height, despite distinct source profiles. Their new P₂, however, never agrees. For a positive reader let L₊ and L₋ count nonempty paths in I starting at its apex i, weighted by P₂ of their complement and by +1 or (-1)^path length. The odd-cycle packing identity gives

    P₂(output) = 2 P₂(I) P₂(B) + 2 L₊ P₃(B) + 2 L₋ H(B).

A cycle crossing blocks necessarily uses the unique o→i return and one ascending path through I, optionally one path in B, then the outgoing singleton. Its odd-length condition is picked by summing the signed and unsigned reads of the I and B subpaths; the outgoing singleton accounts for reversing the sign of the parity condition. For B, the optional selected path with a positive sign gives P₃(B); with a negative sign it gives P₁(B), since each odd cycle or singleton in the pack can be assigned either of two positive colors or one negative color. A direct edge i→o missing from the graph would contribute only a two-cycle and is excluded by this same odd parity.

Here L₊>0 and L₋<0. To see the latter without any pack argument, for a designated vertex i let F(I,i) be the displayed signed weighted starting sum. The selected singleton contributes -P₂(I−i); all longer paths by their next vertex j give -F(I−i,j). Inductively F(I,i) is minus the ordinary weighted *terminal* count at i with complement weighted by H, because P₂(I−i) is H(I−i) plus the sum of all those terminals there, and partitioning according to j→i gives exactly the ordinary weighted-terminal recurrence. It is strictly negative. Hence from C to D each contribution to the P₂ difference is negative: Δp=-16, ΔP₃=-1296, Δh=32. Every subsequent strong two-run reader from the cone therefore distinguishes the two outputs (their H is equal, but α>0). This is loss of a scalar observation with a physical way to expose the retained continuation data, not a scalar hole.

Run `python replay.py` for literal quotient-substitution checks, the six strong source profiles, the higher two-cover transfer directly on both collision hosts at t=2, pumped arithmetic and finite primality certificates. Normal and optimized Python agree.

## What physical nesting costs: additive visit rank

The same combinatorial argument removes the two-run restriction. For a general homogeneous slot let ν be the maximum matching size in its bipartite return graph O→I. The maximum number of slot runs in a Hamilton path after inserting any nonempty B is exactly

    min(|B|, 1+ν).

For the upper bound, choose a crossing from O to I between each two consecutive runs. They are disjoint by disjointness of the intervening path intervals. For the lower bound use any matching and choose distinct b vertices in B to spell the alternating subsequence as above, then extend to a Hamilton path. No internal arc in B is needed for this subsequence.

When one homogeneous-slot reader is inserted into the slot of a second, the new return graph for the innermost slot is the **disjoint union** of the inner and outer return graphs. An outer I vertex beats everything in the inner graph, and an outer O vertex loses to everything in it, so no mixed return arc is possible. Thus ν adds under actual nesting. In particular d nested nontrivial strong two-run readers require d+1 slot runs on large enough inputs.

This higher observation cannot be bypassed by just a different algebraic basis of the same two counts. If r=1+ν, the Hamilton count of the nested reader on B is a sum with nonnegative skeleton multiplicities over ordered covers of B into j nonempty paths for 1≤j≤r. The multiplicity at r is positive by the attained maximum. Inclusion-exclusion expresses the j-nonempty-cover count as a linear combination of P₁,...,P_j with coefficient 1 on P_j. The P_j are linearly independent functions of B, already on transitive TT_n with P_j(TT_n)=j^n. It follows that the coefficient of P_r really survives. Iterated physical passive composition therefore does not remain a two-coordinate scalar program even if each primitive reader does. The replay independently checks the predicted maximum run count on twelve literal nested cases, including small B that cap the rank.
