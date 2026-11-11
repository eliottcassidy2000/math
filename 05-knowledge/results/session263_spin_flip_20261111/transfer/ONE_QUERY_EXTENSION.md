# A single physical future count encodes every bounded-depth source state

**Status: PROVED**, using THM-5130(B) at its inherited audited status and the new free-action theorem in `PROOF.md`. No scalar-target universality is claimed.

## 1. A mixed-radix continuation with an explicit decoder

Let s=(r,M,h,E) be a literal selected tournament and set q=h+r. Assume r is known, 0<=q<=B and 0<=E<=B, where B is an integer >=1. There is no need to bound M for this first theorem. Let m be any nonnegative integer such that

    3^m >= B 2^m.                                  (1)

Append the THM-5130 distance-1 transitive sign-flip gadget with first mark i=1 and m vertices following the adjacent marked pair. The transitive gadget has m+2 vertices and its closer adds one more. The resulting slope is exactly

    S=10·3^m M+3·2^m q−E,                         (2)

and the new offset is −r. Close this new pair with ONE transitive tail vertex and its closer. The single resulting tournament has

    N=H=−r+2S,                                    (3)

and its order is the original order plus m+5. All gadget sizes and final tail sizes are legal; m=0 means zero vertices FOLLOWING a nonempty two-vertex gadget, not an empty gadget.

**Theorem 6 (one-query state encoding).** Given N,r,B,m under these hypotheses, recover the complete original state uniquely as follows. Set S=(N+r)/2 and A=10·3^m, D=3·2^m. Then

    M = nearest integer to S/A,
    q = nearest integer to (S−AM)/D,
    E = AM+Dq−S,
    h = q−r.                                      (4)

There are no ties. Consequently the single physical continuation (3) is injective on all such states with the same r.

*Proof.* The remainder R=Dq−E satisfies

    |R| <= DB+B <=4·2^m B <=4·3^m < A/2.

Thus the first nearest-integer operation is exact. Also (1) implies 2^m>=B: for m>=0, 2^m >=(3/2)^m >=B. Therefore E<=B<=2^m<D/2, proving the second nearest-integer operation exact. The last two equations recover E and h exactly. This includes q=0 and negative remainder −E; ordinary floor division would fail there, so nearest-integer rounding is essential. ∎

The map is a genuine encoding of the entire interface, not just an assertion that a generic linear functional should separate finitely many points. It gives a uniform explicit continuation and an integer decoder. It retains the finite head and endpoint correction −E. In particular, erasing E would destroy the last radix digit.

## 2. Uniform finite-depth cost and exact scalar diversity

Use the initialized free grammar of `PROOF.md`: the marked TT2 source and all binary F3,F5 words of length k. Here r=(-1)^(k+1) is known, M<=3(57/2)^k, and q,E<=M/2. For simple integer bounds choose

    B_k=3·29^k,
    m_k=min{m>=0: 3^m>=B_k 2^m}.                 (5)

Use the SAME continuation m_k on every source of depth k. This gives EXACTLY 2^k distinct scalar Hamiltonian path counts, because Theorem 6 recovers the source and the free-action decoder recovers its word. The physical order is bounded by

    4k+2+m_k+5 <=14k+11.                           (6)

Indeed B_k<=2^(2+5k), and (3/2)^(2(2+5k))>=2^(2+5k)>=B_k, so m_k<=4+10k. More sharply m_k=ceil(log_(3/2)(B_k)). No logarithmic rounding is needed in implementation; condition (5) uses exact integers.

Every encoded count satisfies

    N <=1+26·3^m_k B_k.                            (7)

For B_k>1, minimality gives 3^m_k<3 B_k^p, where p=log(3)/log(3/2). Thus N<1+78 B_k^(1+p). This is an explicit exponential-in-depth bound, with linear-in-depth order overhead.

**Scope of the improvement.** The three-query argument had only 2^(k/3) guaranteed distinct scalar values from 2^k source states. This one-query extension removes that loss IN SOURCE DEPTH: all 2^k states now become distinct scalar counts. It does NOT automatically improve the numerical-cutoff exponent, because the radix gadget increases count magnitudes. The crude cutoff exponent here is log(2)/((1+p)log(29)), approximately 0.0555, whereas the three-query bound was about 0.0690. These are lower bounds for different explicitly described count families, not density records. The one-query theorem is primarily a lossless compiler and physical observability result.

## 3. Exact synthesis and recognition, not arbitrary target coverage

For fixed k, encode a binary word by the physical construction above. Conversely, given a proposed scalar N, apply (4), check the bounded-state hypotheses and exact equality (3), then use the free-action decoder to test whether the recovered state belongs to the initialized orbit with word length exactly k. This is an exact membership-and-construction algorithm for the encoded finite-depth family. Its integer arithmetic and recovered word length are polynomial in the bit length of N and k. No counting oracle is required for synthesis: the inherited matrices compute N from the chosen word, and the literal tournament is produced by the same operations.

If the sign r is unknown, this statement does not claim that one count distinguishes both offset fibres. Known depth parity is part of the input. If no bound on q,E is available, one fixed m cannot provide this radix separation argument. The required continuation is size-dependent, and this dependence is explicit rather than hidden in a limiting argument.

The odd slope of the original grammar remains odd after this continuation: q is even, 10·3^m M is even, and E is odd, so S is odd. Thus the encoder itself also supplies distinct odd-slope sources of the opposite sign at each fixed depth, and their entire rays −r+2^t S are realized for every t>=1. Distinctness of their slopes follows from (3), NOT from an unproved scalar-injectivity claim about the original short-letter grammar.

## 4. Controls and source pins

Run `python3 -B 05-knowledge/results/session263_spin_flip_20261111/transfer/one_query_controls.py`. The exact receipt tests all words through depth 10 with a common size-dependent continuation at each depth, single-count decoding, distinctness, graph-order and magnitude bounds, plus direct subset counts for four small new encoding constructions. A hostile q=0,E=1 example explicitly shows that floor division recovers M−1 rather than M.

The inherited formula is THM-5130(B), file `01-canon/theorems/THM-5130-sign-flip-towers-gadget-transfer-law-and-seeds-of-both-signs.md`, pinned at file commit `7c14d6741eab74550c9337f6143751399d36053e`; specialize delta=1, i=1, a=m+2. New freeness and source bounds are Theorems 1–3 of this lane's `PROOF.md`. The single-query map, exact radix inequalities, all-depth distinct-count construction, and bounded-depth cost are new deductions. Chan–Pak and all-odd-coefficient coverage remain OPEN.
