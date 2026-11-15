# Changing the anchor: a local quadratic atlas that still misses almost every prime

**All-order identities and scoped bounds.** Start with the [one-outlet replacement theorem](../chanpak_chronological_core_changes_20261114/ONE_OUTLET.md). A follow-up derivation noticed that a nonempty ordinal tail on the exterior makes an anchored parity sum vanish. This gives a new family with two independently replaceable cores. It does not give Chan–Pak. The prime-density theorem below concerns only a finite anchor menu with *equal normalized core ratios*.

## Exact physical gluing

For a tournament X let h(X) count Hamiltonian paths, with h(empty)=1. Let Q(X) count ordered two-slot path covers allowing an empty slot. Put P(X)=Q(X)−h(X). The bounds P(X)≥|X|h(X) and h(X)≤|X|! will be important. The first counts the distinct cuts of each full path.

Take R with marked q, and nonempty arbitrary cores W,U. Orient R→W→U, U→q, and R−q→U. Let A=h(R), K=h(R−q), and B be the number of two-slot covers of R with first slot starting at q. Then
```
N = (A+K) h(U)h(W) + P(U)[B Q(W) − K h(W)].
```
The result is strong exactly when q can reach all of R in R: the only port back into R is q, while R reaches both W and U. The case of q a source of R is an old singleton cyclic wrap: R−q and W combine into one ordinal block. Use a genuine nonsource to change that quotient.

For the proof, first replace U by one vertex v. The exterior in the one-outlet formula is R→W. The cycles through v correspond to a first path starting at q in an ordered two-cover of R→W, of length at least two in vertices. The first path supplies the cycle and the second supplies an ordering of its complement. Before the length restriction there are BQ(W) covers, by independently splitting the two slots at the ordinal boundary. The signed first-slot length sum vanishes, since the W side is nonempty and odd-cycle packings cancel the two slot assignments of any odd block. The length-one contribution is Kh(W); it has negative sign. Writing E and O for the even and odd cycle totals of the old host, the unrestricted even and odd first-slot totals are equal. Even cycles require odd first-slot size, so E=BQ(W)/2−Kh(W) and O=BQ(W)/2. The exterior has h(R→W)=Ah(W). Inserting those data into the one-outlet replacement proves the identity. See the [odd-cycle packing expansion](https://arxiv.org/html/2412.10572v3#S3.SS4.Thmcorollary20).

### Absolute inverse without synchronization

Fix U and put L=B P(U), M=(B−K)P(U)+(A+K)h(U). Both are positive since B≥K (take first slot just q). The unknown W's coordinates t=h(W), z=P(W) have to satisfy
```
N = L z + M t,     z≥|W| t,     t,z positive odd.
```
Thus choose t by the exact linear congruence Mt≡N (mod L), bounded by N/(M+L). Recover z=(N−Mt)/L and cap the order of W by floor(z/t). The gcd test and physical profile test are independent. For transitive U and W of positive orders m,s,
```
N−(A+K) = (2^m−1)(B·2^s−K).
```
This is a literal divisor decoder with only O(log N) candidate m for each anchor, independent of any bounded-order search for W. The transitive pair lane has at most O(log² N) distinct values through N for any fixed anchor. Do not infer coverage from exact inversion.

## A three-form atlas removes the target-prime splitting obstruction

Equal normalized ratios P(U)/h(U)=P(W)/h(W)=x/y in lowest terms give h(U)=ay, h(W)=by, and
```
N = ab F(x,y)
F(x,y) = B x² + (B−K)xy + (A+K)y².
D = (B−K)² − 4B(A+K).
```
Here x,y,a,b are positive odd and x≥max(|U|,|W|)y. A prime target has a=b=1. For target 7p, p>7 prime, the only possibilities are F=p with ab=7 or F=7p with ab=1. In fact F≥A+2B≥3, and F=7 would force x=y=1 since x=y is primitive and x>y gives x≥3, hence F≥11. At x=y=1 both cores have order one and a=b=1.

The following literal padded anchors are listed by outgoing bit rows. They are all strong after gluing, and none has q a source. The display includes outputs for U=W=TT_j, j=1,2,3.

| R | q | (A,K,B) | D | Literal outputs |
|---|---:|---|---:|---|
| (60,57,58,48,32,0) | 0 | (3,1,24) | 145 | 51,289,1341 |
| (56,49,51,54,32,0) | 2 | (5,3,36) | −63 | 77,431,2003 |
| (112,97,99,103,110,64,0) | 1 | (9,5,40) | −1015 | 89,479,2219 |

The square classes are 145, −7 and 145·(−7). For every odd target prime p outside the finite set dividing the leading coefficients or discriminants, at least one D is a quadratic residue mod p, since the product of the three Legendre symbols is +1. For that form, choose y=1 mod p and take a root of F(x,1)=0 mod p; the root is simple. Hensel lifts to F(x,1)=p in Z_p, or to F(x,1)=7p there. For the **same form** the derivative in x is odd when y=1, since all B are even and K odd. F(1,1) is odd. Therefore every odd target, including p or 7p, is represented with x,y odd in Z_2. The two congruence solutions can be combined by CRT to any prescribed precision at 2 and p simultaneously. This clears exactly these two local tests, not all places, and is not a claim of physically realizable local (x,y).

The pattern is a three-character circuit: a single quadratic norm forces splitting in one field, while a product relation on three characters ensures at least one of their split lanes is available. The distinction between a character/witness equation and a physically admissible positive inverse is the same useful distinction made by the old [fruit elliptic tree](../prime_shells_20260921_elliptic.md) and the [squares/doubles transport deck](../collatz_sqdbl_20260925_squares_doubles_foundry.md). These are analogies about what data a transport loses, not reductions of Collatz.

## Dyadic tails create infinitely many discriminant square classes

Padding an anchor R on the right by TT_s multiplies B by 2^s and leaves A,K unchanged (q stays in R, and its anchored first-slot condition does not see the tail). Hence with b=B2^s,
```
D_s = [b−(2A+3K)]² − 4(A+K)(A+2K).
```
For all sufficiently large s this lies strictly between consecutive nonnegative squares: if z=b−(2A+3K), require 2z−1>4(A+K)(A+2K). For example a transitive five-vertex R with q first has (A,K,B)=(1,1,16) and D=97, giving 191 with two TT_2 cores, but this example is a singleton wrap in disguise.

More strongly, every fixed rational square class occurs only finitely often. Fix its squarefree integer representative d and split s by parity. For either parity, b=c u² with c fixed positive and u a power of two. Square-class membership gives
```
d v² = c²u⁴ − (4A+6K)c u² + K².
```
The quartic has four distinct roots over C: the quadratic in u² has positive nonzero discriminant and nonzero constant term. The affine curve is smooth of genus one after completion, so [Siegel's theorem on integral points](https://www.math.ens.psl.eu/~obenoist/refs/Silverman.pdf#page=280) makes its integer points finite (this is a finiteness conclusion, not an effective search bound). D=0 can occur at most twice. In particular the multiquadratic field generated by the nonzero D_s has infinite degree, since finite-dimensional square-class subspaces over F_2 are finite.

These new arithmetic fields do *not* make new boundary wires. The tail can instead be absorbed into W; then P(W)/h(W)=x/y transforms to 2^s(x+y)/y−1 while U's ratio is unchanged. Varying the discriminant is one controlled way to break synchronization **relative to the original anchor**. It does not remove the need for physical profile choice.

## In spite of the local atlas, a fixed synchronized menu misses almost every prime

For every fixed finite set of anchors, the synchronized construction above produces at most O(X log log X/log² X) distinct prime targets p≤X, and at most that many prime cofactors p≤X of 7p targets. The implied constants can depend on the menu. Here is a sieve proof that covers real and imaginary discriminants alike.

For a fixed form, F has nonnegative coefficients with B>0 and A+K>0, so F≤7X bounds x,y≪sqrt(X). Discard y≤X^(1/4) at cost O(X^(3/4)) lattice pairs. Otherwise h(U)≥y and h(U)≤|U|! imply |U|≫log X/log log X. Since x≥|U|y, the remaining pairs lie in a rectangle of side lengths H≪sqrt(X) and J≪sqrt(X)log log X/log X.

If D is square, 4BF is the product of the two positive integer factors 2Bx+(B−K±sqrt(D))y. Both are at least 2Bx, since B−K≥sqrt(D) in the square case. For F=p or 7p at least one factor is bounded by 28B, hence x and y are bounded. Now suppose D nonsquare. Use the odd primes ℓ not dividing 7BD. The number ρ(ℓ) of zeros mod ℓ is 2ℓ−1 when D splits (a union of two distinct lines), and 1 otherwise. For squarefree d supported on these primes, the count of lattice pairs with d|F is HJρ(d)/d² plus O(ρ(d)((H+J)/d+1)), by counting each residue class in the rectangle.

Use the [Selberg upper-bound sieve](https://people.math.harvard.edu/~elkies/M229.20/sieve.pdf), with z=X^(1/40). Its main denominator G(z) is the sum over squarefree d<z of products of ρ(ℓ)/(ℓ²−ρ(ℓ)). In the split lane the factor is at least 2/ℓ. For completeness, the product of 1+2/ℓ over good split ℓ≤z^ε is ≫log z by Mertens in the fixed reduced residue classes (half the classes split). In this expansion weighted by 2^ω(d)/d, the mean of log d is at most 2∑_{ℓ≤z^ε} (log ℓ)/ℓ=O(ε log z); choose a fixed small ε so at least half the weight has d<z. Thus G(z)≫log z. In the error sum over d<z², 3^ω(d)ρ(d)≤6^ω(d)d≤d⁴, so the total error is O((H+J)z⁸+z¹⁰). Prime targets or prime cofactors p≤z add at most z distinct values. The bound now follows from HJ/G(z) and the small-y discard.

This is an explicit separation: there is **no simultaneous inert-prime exclusion for the atlas**, but even that atlas, or any fixed finite augmentation in equal-ratio mode, realizes a density-zero fraction of primes and 7p targets. The missed targets of the menu are not asserted to be missing in the tournament spectrum. A route that remains open is to let the target select an *unsynchronized* pair via the linear inverse, then solve the bounded physical profile problem uniformly. These bounds do not establish that such a choice always exists.

## Reproduction

`python replay.py` and `python -O replay.py` independently count the literal glued tournaments by subset/end-point DP. The checks include 600 seeded triples and all nine displayed values; [receipt](replay.txt). The infinite conclusions are proved above, not inferred from the trials.
