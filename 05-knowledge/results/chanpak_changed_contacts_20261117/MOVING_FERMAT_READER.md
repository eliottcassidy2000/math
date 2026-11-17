# A changed right reader lands with a transitive left block

The fixed-reader obstruction in [the Fermat-target note](FERMAT_TARGET.md) does not require changing the left block to a nontransitive one. Changing the right tournament itself suffices at both 257 and 65537. The new 65537 witness is strong of order 16, one vertex smaller than the previous construction; no minimum-order claim is made.

In the left scaffold (8,1,3,6), replace vertex 0 by TT_m, mark vertex 1, and apply the same one-arc seam into the converse of the right reader. With P=2^m its count is

    N = (3P−1)d + (4P+1)e − (P−1)f.

Here d is the full right path count, e the ordered terminal cover count, and f the full starting count at its marked port. The following are literal choices. Outgoing masks refer to the right reader before taking its converse.

| target N | m | right reader rows | port | (d,e,f) | direct affine check |
|---|---:|---|---:|---|---|
| 7·257 | 1 | (34,124,9,17,37,76,29) | 6 | (101,146,20) | 867·2+65=1799 |
| 7·65537 | 4 | (34,460,9,177,167,326,157,293,93) | 2 | (1907,5711,139) | 28426·16+3943=458759 |

Both right readers have a Hamiltonian path ending at the port. After taking the converse the two frames close through the seam, so the outputs are strong. The replay verifies the profile reads and the whole assembled graphs by independent subset DP. These are constructive landings at counts already within the established numerical range, not new range coverage.

For a prescribed N and m the profile-selection problem is exact: `d≤N/(6P+1)` and `(P−1)f≡(3P−1)d−N (mod 4P+1)` with 0≤f≤d and e≥d. The gcd controlling this congruence is gcd(P−1,5); in particular when m≡0 (mod 4) this forces N≡2d (mod 5). These constraints prune the search but do not select a physically realizable right tournament for every target.

## Changing only this right reader still cannot cover all sevenfold targets

Allow an arbitrary nonempty B of h=H(B), P=P₂(B), along with any right reader. The exact count remains

    N=(3P−h)d+(4P+h)e−(P−h)f ≥ (6P+h)d,

since e≥d and f≤d. Also P≥(n+1)h and P≥2^n for n=|B|, and P≡0 (mod 4) if n≥2. If d=1 the reader is transitive and its profile must be (1,1,1) if the port is the source, or (1,2^j,0) with j≥1 predecessors. These elementary conditions leave **no** possible profile at N=35, and only the following at N=119:

    (h,P,d,e,f)=(1,2,5,11,5) or (1,2,9,9,7).

For transparency, no tournament census is needed for this pruning. At n≥3 the bound gives d=1; for d=1, the source case 6P+h has incompatible h modulo 24 and size bound when n≥2, the j=1 case 11P+h has incompatible h modulo 44, the j=2 case gives 79 at n=2 and exceeds 119 at n≥3, and j≥3 already exceeds 119 at n≥2. At n=2 and d=3 the equation at N=119 is 17e−3f=86 with f≤3, impossible. At n=1 the exact equation is 5d+9e−f=N with d≤N/13 and d≠7. This leaves the stated profiles after applying the transitive d=1 rule. The replay also verifies the finite arithmetic screening.

The first alleged profile at 119 has f=d: all Hamiltonian paths start at the marked port. By strong-component decomposition and Camion's theorem (rotating a Hamiltonian cycle gives every start in a nontrivial first component), this forces the port to be a source and hence e=d, contrary to e=11. The second has e=d, which forces every terminal cover to concatenate. Using the singleton port plus any path of the deletion shows f equals that deletion's full count, so f=7 is impossible. Thus this left scaffold misses both 7·5 and 7·17 even with a changing marked right tournament and arbitrary B. Both counts are realized elsewhere.

## Padding two fixed cores does not make a pure Fermat family

Fix a left module B₀ with profile (h,P₀) and a pointed right tournament with triple (d,e,f). Ordinal-pad B₀ by m singleton vertices, and prepend l new ordinal sources to the represented right tournament (retaining its port). For l≥1 the right starting count becomes zero and its terminal cover count is 2^l e: the transitive prefix can be split arbitrarily between the two paths. The total is

    N = 4eP₀·2^(m+l) + 3dP₀·2^m + eh·2^l − dh.

Suppose N=7(2^t+1). If a=min(m,l)≥1, then t≥a: the first term alone less the fixed dh is at least dh·(8·2^(2a)−1)≥7(2^a+1), using e≥d and P₀≥2h. Reducing the equality modulo 2^a now gives 2^a | dh+7. Thus one of m,l is bounded, even if both are chosen from t.

If infinitely many targets were hit, pass to a subsequence with one padding length fixed. When l grows and m is fixed, the expression is A·2^l+C where C=d(3P₀2^m−h). Since both l and t grow, comparison mod 2^min(l,t) forces C=7. But C is the product of two attained positive odd counts (the right count and the full scaffold count), and neither can be 7. If m grows and l is fixed (including l=0), the analogous constant comparison forces the affine intercept to be 7. The argument in [the fixed-reader obstruction](FERMAT_TARGET.md#why-no-fixed-right-reader-gives-the-full-pure-fermat-ray) then applies: the left module's fixed h must divide 7, so h=1, and the remaining normalized slope is either odd ≥5 or forces a deletion of count 7. Neither gives a power of two even if t is selected independently of m. Therefore this entire two-height padding scheme hits only finitely many pure Fermat-form targets.

## A general bound on a fixed carry

The earlier 65537 witness used P=27h−1 and P=2^s to turn 127P+27h into 128P+1. More generally a physical profile with P=A h−c for positive integer c has `n≤A−2`, because splitting gives P≥(n+1)h. Thus a carry with **fixed** coefficient A only uses finitely many block profiles. If P is a power of two, then `P<A(A−2)!` by h≤n!, so producing an unbounded power-of-two coordinate requires A to grow as well.

For A=27, c=1 we get n≤25 and P<2^89. The congruence 2^s≡−1 (mod 27) gives s≡9 (mod 18), leaving s+7=16,34,52,70,88. Only 16 is a power of two. Hence 65537 is the only *classical Fermat number* this exact carry can hit, although some other 2^t+1 forms are still arithmetically possible. A repeatable construction needs a changing coefficient, not just an inflated source.

Replay: `python replay_moving_reader.py` (also verified under `-O`).
