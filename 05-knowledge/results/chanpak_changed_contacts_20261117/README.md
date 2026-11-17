# Two new vertices contacting a whole growing slot: exact induced factor-seven return

**What this session established.** A growing slot can be contacted by two new vertices to return exactly from `H=p+8` to `7p` on an infinite ray. [Changing the contacts independently](HETEROGENEOUS.md) gives a second prime-capable ray and a physical strong tripler. An independent [single seam](MERSENNE.md) realizes every `7(2^k−1)` for `k≥3`, starting immediately after 7 and 21. The homogeneous braid has a sharp all-height classification and a short list of remaining changing-host interfaces. None of these results is universal prime placement.


Take a tournament D with at least two vertices, with marked vertex q. Add a transitive slot B=TT_k, k >= 1, such that B beats q and all vertices of D-q beat B. Call this source T. To make T', add a second nonempty transitive slot A=TT_l with D→A→B. The old arcs are kept. Define

* G=H(D), J=H(D-q), e=number of Hamiltonian paths of D ending q, L=G-e.
* W=number of ordered partitions of D into directed paths X,Y with X possibly empty, Y of length at least two vertices, and Y starting q.

For either tournament a Hamiltonian path visits B in at most two runs: every nonterminal B-run must exit to q, which occurs once. For k vertices, precisely one assignment to one B-run and 2^k−2 assignments to two ordered nonempty B-runs are possible. Each assignment supplies unique increasing orders within the runs.

In the source, one B-run can be last (L possibilities), or occur in X B Y with Y starting q (W+J possibilities, including exactly J with Y={q}). Two B-runs necessarily have X B Y B with Y starting q and length >=2 (W possibilities). Thus H(T)=L+J+W(2^k−1).

In the child, put s=2^l. For just one B-run, the terminal B-run requires Dpath A B, with G skeletons. Otherwise use X [A] B Y [A] with Y starting q, and A material can be apportioned arbitrarily between the two indicated slots (exactly s possibilities; empty slots are allowed). There are W+J choices of (X,Y). For two B-runs, apportion A between the slots immediately preceding them. If len(Y)>=2 the allocations are unrestricted; if Y={q}, the second A-slot must be nonempty, because q loses to B. Therefore

    H(T') = G+s(W+J) + (sW+(s−1)J)(2^k−2).

## Explicit strong family

Use D outgoing masks (26,4,25,18,2), q=4, indexed from zero. Direct small counting gives (G,J,e,W)=(9,5,4,5). At l=2 the source six rows before blowing up B are (58,36,57,50,2,16); the child eight rows are (250,228,249,242,194,16,32,96), both with the B slot at 5. Binary masks refer to arcs directed out of the row index. The quotient source and child are strong; replacing a vertex of a nontrivial strong tournament with any tournament preserves strength, since directed excursions through the rest of the quotient connect vertices of the replaced module. In fact for this D and all k,l >=1,

    H(T)=5·2^k+5,
    H(T')=(2^(l+1)−1)(H(T)−10)+14.

In particular l=2 gives the exact induced restoration H(T')=7(H(T)−8), for every k. Prime quotients include 17,37,157,317 at k=2,3,5,6; the entire ray is physical with or without primality. This does not assume a fixed normalized dyadic root menu and contacts a whole slot whose size grows.

For any tournament B in place of TT_k in the fixed l=2 fixture, the same skeleton count gives H(T')=7H(T)−56H(B), since one run contributes H(B) and two runs contribute the number of ordered covers by two nonempty paths. The constant affine intercept uses H(TT_k)=1.

## Uniqueness of the output ray within this architecture

Suppose varying k gives H(T')=7(H(T)−8) for fixed D,q,l. Match the coefficient of 2^k to get (7−s)W=(s−1)J. In particular W>0, so q has at least one outneighbor t. Split each Hamiltonian path of D-q before t, inserting q before the suffix; this injects the J paths into W, because one recovers them by deleting q and concatenating. Thus W>=J. The coefficient identity forces s=4 and W=J: s=2 requires J=5W, and s>=8 is impossible.

After matching slopes, the constant equation reduces to 7e=6(G+J)−56. Since e<=J just by deleting q, and G>=J by inserting q at a legal position in each Hamiltonian path of D-q, one gets J<=56/5. The positive even integer G+J is divisible by 7, and 6(G+J)/7−8=e<=11 forces G+J=14. Then e=4 and 5<=J<=7, with J=7 excluded by the known tournament hole. Consequently J=W=L=5, G=9, precisely the same ray above. In particular this is an infinite changed-contact lift but not a universal arbitrary-target selector.

Replay: `python replay.py` (and `python -O replay.py`). It counts small quotient skeletons directly by permutations and uses an independent subset DP on inflated tournaments through k=7. The arbitrary-host cutoffs below are proofs, not extrapolations from this check.

## The same rigidity for arbitrary homogeneous A

The formulas for the child remain valid replacing the term G by hG, and s by p elsewhere, with h=H(A) and p=P_A(2)=Σ_S H(A[S])H(A[A\S]), taking empty count 1. The unrestricted allocation uses p ordered path covers with an optional empty part; forcing the second part to be nonempty subtracts h. Matching 7W now requires (7−p)W=(p−h)J>0. Indeed p>h by the two endpoint splits, and J>0. Every subset of A contributes at least one by Rédei, so p>=2^{|A|}. It follows that A has at most two vertices, hence reduces to the earlier classification.

## What if the target chooses D and A?

Let h=H(A), P=P₂(A), t=2^k and e=G−L. The exact same computation gives the following particularly useful separation of the constant and the growing part:

    Δ := H(T')−7H(T)
       = (h−7)(G+J)+7e +(t−1)[(P−7)W+(P−h)J].

The desired constant is Δ=−56. Splitting a Hamiltonian path of D immediately before q shows L≤W. Also e≤J and J≤G by deletion/insertion. Hence G≤W+J. If h≥7 all terms in the original count difference
`(h−7)G+7e+(P−7)(W+J)+(t−2)[(P−7)W+(P−h)J]`
are nonnegative. If h<7 but P≥7 and P+h≥14, the first and third terms together are at least `(P+h−14)(W+J)`, and the others are nonnegative as well. As P≥2^|A|, any possible returning A is TT₁, TT₂ or TT₃ (the directed triangle has h=3, P=12).

If A=TT₁ or TT₂ and W=0, q is a sink, so G=e=J. The two deficits become (t−6)J and (3t−8)J, neither −56 for odd J and powers t≥2. Otherwise W≥J; put z=W−J≥0. For TT₁ or TT₂ the positive equations for −Δ are respectively

    56 = 6L−e +(4t+2)J+5(t−1)z,
    56 = 6L−e +6J+3(t−1)z.

Using e≤J and H(T)=L+tJ+(t−1)z≤(t+1)J+tz proves H(T)≤21 for TT₁ (and k≤3); since 21 is forbidden, H(T)≤19. For TT₂ with z>0 it proves k≤4 and H(T)≤71. These bounds follow by maximizing (t+1)J+tz over positive odd J and z≥0 or z≥1 with the resulting linear inequality, for t=2,4,8,16; the replay reproduces the tiny arithmetic maximizations. If TT₂ has z=0 it is exactly the already classified ray, even with a new choice of D at each target.

For A=TT₃,

    Δ = −6L+e +(7t−13)J+(t−1)W.

Since L≤W, t≥8 gives a strictly positive deficit. Therefore the **only potentially unbounded alternative in this contact geometry** is A=TT₃, k=1 or 2, with D changing. At a prescribed p the exact two tests on the host are

    L+J+(t−1)W = p+8,
    7[L−(t−2)J]−e = p+64.

These tests include the absolute count; they are not just congruence conditions.

### The transitive deletion cannot supply that last alternative

If D−q is transitive, J=1. Read the indicator of q's outneighbors along the transitive order as a binary integer w, earliest bit most significant. Then W=w: each second path is determined by its first neighbor after q and a subset of later vertices. Write b=w mod 2. We have e=1−b, G=1+2C(w), L=2C(w)+b. The weighted-cut expansion from [THM-4114](../../../01-canon/theorems/THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature.md), or directly the [odd-cycle expansion](https://arxiv.org/abs/2307.05569), gives

    C(0)=0,   C(2r+1)=C(r),   C(2r)=C(r)+(r+(r mod 2))/2  (r>0).

Indeed a new terminal zero creates even path segments from each earlier 1 into the new terminal vertex, with weight 1 if adjacent and 2^(distance−2) otherwise. A terminal one creates none. This keeps the extra paths that would be lost by assuming that deleting q always leaves a concatenated Hamiltonian path.

For t=2 the equation is w−12C(w)=−58+7b. If w=2r, the doubling rule forces r≤14. If w=2r+1, then r−6C(r)=−26 and r is even; putting r=2u gives u+6C(u)+3(u mod 2)=26, so w≤105. For t=4 the equation is 3w−12C(w)=−72+7b; reducing mod 3 gives b=0, and substituting w=2r gives r even with C(r)=6. Thus w=4u with u≤12. The replay finds no solutions within these **proved exhaustive bounds**. Consequently even this last alternative must change the nontransitive deletion, or relax the exterior contacts.

### A constant tail word cannot repeat the return

More generally fix J and grow D−q by ordinal singleton suffixes, recording a 1 if q beats the new vertex. The state (G,W,e) transforms as follows, using the [marked endpoint cancellation](../runner_aliquot_relay_20261117/NEAR_SINK_PEEL.md#endpoint-cancellation):

    1: (G, W, e) ↦ (G, 2W+J, 0)
    0: (G, W, e) ↦ (G+J+W−e, 2W, J).

For any nonempty word of length m, the linear coefficient of W in the net change in the TT₃ deficit is (t−1)(2^m−1)−6S, where S is the sum of 2^(i−1) over the zero positions i. It is odd, hence nonzero. On iterating a fixed word, e is fixed after the first repetition and W grows strictly unless it is identically zero (all zeros starting at W=0). Thus a fixed nonempty tail word cannot preserve the deficit at every repetition with growing W, or give the same exact deficit at unbounded repetition counts. This leaves target-dependent, nonperiodic contacts and changes of deletion genuinely distinct.

## Why a single opaque slot is limited in any fixed exterior

This is independent of the braid. Replace one unmarked vertex in any fixed tournament template by B with uniform exterior contacts. Write R_j(B) for ordered covers by j nonempty paths. Any Hamiltonian path count of such a template is a fixed nonnegative linear combination of R_j: collapse the module's consecutive runs to symbols. The same observation applies to a Hamiltonian path ending at an exterior mark and to a rooted ordered two-cover. Thus either a fixed affine comparison of two such counts, or the earlier near-sink inverse equations (using endpoint cancellation for E), has the form

    Σ_{j=1}^d a_j R_j(B) = constant.

For n=|B|, splitting any one path of a j-cover into two gives n−j choices, and a resulting cover has at most j preimages (one merge at an adjacent seam). Therefore `j R_(j+1) ≥ (n−j)R_j`. If the highest nonzero coefficient has j≥2, its term eventually dominates both every lower term and the fixed constant; all solutions have bounded n. If only R₁=H(B) occurs, it fixes H(B). For the aliquot deficits −28 or −56 that fixed odd H(B) can only be 1 since 7 is forbidden. Then B is transitive and any nonconstant output below X forces n=O(log X): if a term R_j with j≥2 is read, it grows at least as `j^(n−j)` by assigning one vertex to each of j paths and assigning the others freely. Hence a finite bank of fixed single-slot templates contributes only O(log X) exact returns. The new physical ray is at that scale, not a claim of coverage.
