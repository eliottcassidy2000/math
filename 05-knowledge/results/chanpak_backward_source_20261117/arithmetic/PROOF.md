# Backward exact height: a unique tall candidate and finite fine-switch sources

**PROVED universal source-order bound; PROVED necessary inverse tests and density obstruction for the specified reader language.** Actual coverage remains open.

Use the branching reader of `pubsub_branch/PROOF.md`. Its legal state is t, its coefficients are a=12t-3, b=15t-4, and its read is N=aH+bP₂. Legal t are positive odd and their counting function is Θ(X^σ), where 82^(-σ)+5976^(-σ)+96264^(-σ)=1 and σ≈0.1450940571. We use the proved language count, not an extrapolation from word samples.

## 1. The target selects a single candidate once the reader is tall

As in `height_select/archive/PROOF.md`, let J=4H+5P₂. Since5a-4b=1,

    H=5N-bJ, P₂=aJ-4N.

Suppose the source must have a same-order companion with profile(H+10,P₂-8), preserving the physical output+2. The parity identity P₂=(-1)^n H mod3 forces n even. Order2 has no such pair, so n>=4. Necessarily P₂-8>=4(H+10). Consequently every genuine lower source lies in

    ceil((24N+48)/(a+4b)) <= J <= floor((5N-1)/b),
    J=12 mod24.                                             (1)

The congruence is the necessary even-order gate P₂=0 mod8 and P₂=H mod3. The real width of this interval is

    N/[b(a+4b)] - 1/b - 48/(a+4b).

If N<24b(a+4b), (1) contains at most one admissible J. An exact constant-number-of-divisions test suffices: take the largest J=12 mod24 below the upper bound and reject it if it misses the lower bound. There is no hidden need to search a huge source-parameter interval at tall height. Passing this test is not a proof of physical source realization.

The source H cannot equal1. A tournament with H=1 is transitive, has P₂=2^n, and every tournament at that order has P₂>=2^n by the all-singleton term in the odd-cycle expansion. Its proposed companion would violate this lower bound. Thus H>=3, so

    N >= (a+4b)H+48b >= 936t-249.                           (2)

Reader states above(N+249)/936 cannot support the native+2 at all. Known source holes7,21 and the forbidden companion H+10=21 give further rejection gates H not in{7,11,21}; these are inherited spectrum exclusions, not new proofs of those exclusions.

## 2. Fine source pairs have a uniform finite order bound

Let T,T' be tournaments of the same order n, with H(T)=h>1, H(T')=h'>1, and nonzero Δ=P₂(T')-P₂(T). Then

    n <= max(h,h') + v₂(Δ) - 2.                            (3)

Proof. The singleton-free SCC core C of T has order r<=h, by the all-order finite-core argument in `pubsub/publisher/PROOF.md`. Since h>1, r>=3. Its P₂ is divisible by4: in the odd-cycle expansion at q=2, the all-singleton term is2^r, and a nontrivial packing has at least one cycle and at least one component, hence a factor at least4. Reinsert the n-r singleton SCCs to get

    v₂(P₂(T)) >= n-r+2 >= n-h+2.

Apply the same argument to T'. A difference is divisible by the smaller of the two guaranteed powers, proving(3) whenever n-max(h,h')+2 is positive; otherwise(3) follows trivially because Δ is itself divisible by4.

For the native switch Δh=10, ΔP₂=-8 this gives the explicit **even-order** bound

    n <= H+11.                                             (4)

This is independent of how the sources are constructed. Hence a proposed low-H repair cannot evade inspection merely by moving to enormous source order. For any fixed candidate H, all physical source pairs are decidable by enumerating finitely many tournaments through H+11 and checking their profiles (or checking literal one-arc adjacency if the task additionally requires that edit).

Combining(1) and(4), a tall reader leaves at most one source profile and a proved finite source-order ceiling. This is a genuine exact backward reduction, albeit not a practical general enumeration algorithm.

## 3. Unlimited padding and bounded source H do not enlarge the reader exponent

Permit *every* source tournament with H<=K, not just known examples, and arbitrary singleton padding. The resulting output set below X has size

    O_K(X^σ).                                               (5)

Indeed every such source has a singleton-free core C of at most K vertices. There are finitely many such cores. Write(H,P) for its profile and let s be the number of singleton SCCs, so the full source has profile(H,2^sP). Since t>=1,

    N=(12t-3)H+(15t-4)2^sP >= 11t2^sP.

The number of allowed t at this padding is O((X/(2^sP))^σ). Summing over s is a convergent geometric series; summing over finitely many cores proves(5). Singleton positions can give different sources but not new(H,P₂) profiles. The proof does not assume the sources support the+2, so it remains valid after imposing that stricter requirement.

There is a quantitative version. The number of possible cores of order at most K is at most

    1+sum_(r=1)^K 2^(r(r-1)/2) = 2^(K²/2+O(K)).

Consequently if K=c sqrt(log₂ X), the output count is at most

    X^(σ+c²/2+o(1)).                                        (6)

For any c<sqrt(2(1-σ))≈1.3076 this is density zero, even relative to odd targets. Thus for density-one odd N, any representation in this reader grammar must use a source whose H exceeds c sqrt(log₂ N). A source library that grows only by padding, or whose H remains below this scale, cannot close Chan–Pak through these readers.

This is not a lower bound for arbitrary tournament constructions; it is a sharp statement about what changing sources must accomplish within this particular branching-reader language.

## 4. What source changes could actually help

A successful replacement must make the nontrivial source core grow with the target, not merely its singleton padding. Moreover it must meet the exact two-coordinate inverse, not just the required slope or displacement. The useful design target is a family of physically realizable source fibers

    J=4H+5P₂, with consecutive allowed H values differing by10,

whose J values intersect the target-selected intervals(1). The source-order theorem makes every proposed local H window testable without an arbitrary order cutoff. It does not guarantee that such fibers exist at the required J.

Allowing new reader invariants changes the coefficients and can evade these reader-specific density bounds; allowing more than two subscribed run coordinates can evade the two-coordinate source formulation. Neither freedom should be confused with a proof of physical coverage.
