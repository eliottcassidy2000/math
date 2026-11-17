# Decimal scales: physical alternatives, a Cantor limit, and the nested prime-offset phase

The archived [decimal window work](../powers_two_decimal_windows_20261007.md) kept the exponent phase, rather than assuming far-apart powers have independent digits. Its suffix analysis studies powers of 2 modulo powers of 5. In the binary-apex seam below the roles reverse: powers of 5 modulo powers of 2 govern the exact repair. The archived [least-hole packet](../least_hole_near_twins_20261117/README.md) also distinguished the Farey-seam exclusions from actual targets and constructed 10^(2j)+11^j using a five-vertex tile with (H,P₂)=(11,100). All new limits are restricted to stated physical grammars.

## Hold the decimal scale while changing the source

The following literal outgoing rows give two-slot count 10000:

| H | rows |
|---:|---|
| 121 | two ordinal copies of (24,1,3,6,14) |
| 311 | (510,44,384,100,430,132,54,74,234) |
| 317 | (482,372,481,439,389,272,312,354,0) |
| 383 | (218,496,379,498,256,17,48,372,97) |

Take an ordinal chain of j choices and wrap it in a directed triangle with two singletons. This is a strong graph and its count is exactly P₂+H, by counting one or two runs through the chain. Therefore every

    10^(4j) + 121^a·311^b·317^c·383^d,        a+b+c+d=j,

is realized. At fixed j the choices are distinct by prime factorization (121=11²). This gives an exact target-side core selector: determine the only possible power 10000^j less than N, subtract it and factor the remainder over those four allowed atoms, requiring exactly j factors. The residual is smaller than the decimal base. This works for **100098587=10⁸+311·317**, independently certified prime. The older one-stock decimal read would instead use offset 121². Fixed prime offsets are not obtained for j>1 in this product grammar.

## An actual fractal of normalized Hamiltonian counts

Let C be the (11,100) tile, and take n ordinal copies of C. One apex either beats a full copy (bit 1) or loses to it (bit 0). Require first bit 1 and last bit 0, ensuring strength. This is a different grammar from the cyclic wrapper with two singletons.

For a bit word b, use the [marked interval law](../chanpak_archive_target_windows_20261117/FLANKS.md). Uniform contact gives marked weight 89 on the active side, signed marked weight −11 on that side (the optional signed total is zero and the empty term has weight 11), and no local extra count. Thus:

    H(b)=11^n
       +7921 Σ_(i<j,b_i=1,b_j=0) 11^(n−1−j+i)100^(j−i−1)
       +121·11^(n−2) · number of adjacent 10 transitions.

This identifies the exact entire grammar, not only approximate plots. Set rho=11/100 and define the compact set

    A={1+Σ_(i≥1) e_i rho^i : e_i∈{0,1}}.

The finite sets of strong reads H(b)/100^(n−2) converge in Hausdorff distance to **7921·{xy:x,y∈A}**, with error O(n rho^(floor(n/2))). For one direction, take the first bits of b and the reverse-complement of its last bits. The product of these two positional values counts the weighted inversions in the display; spurious pairs with overlapping or crossed positions have total O(n rho^(n−1)). The base and adjacent terms vanish at rate O(n rho^(n−2)). Conversely prescribe arbitrary two finite fringes of length floor(n/2) from two addresses in A, fill any leftover space, and the omitted tails cost O(rho^(floor(n/2))). This proves both directions. Covering independent fringes of length m gives at most 4^m intervals of width O(rho^m); the limiting product set has upper box dimension at most log(4)/log(100/11)<1. No claim of equality of the dimension is needed.

The absolute counts have a much stronger cutoff at the right edge. For fixed numbers of ones and zeros, exchanging adjacent 01 to 10 increases the read: it adds an inversion and only lengthens existing affected inversions, whose weights strictly increase with distance. Hence a maximum is a run of ones followed by zeros. With s ones and t zeros the exact read is

    100^(s+t)−100^s·11^t−100^t·11^s+3·11^(s+t).

This also follows directly from the cyclic two-block-plus-apex read (P_X−H_X)(P_Y−H_Y)+2H_XH_Y. The maximum has |s−t|≤1 by minimizing (100/11)^s+(100/11)^t, and is always strictly below 100^n. Conversely the extreme pair alone forces H≥7921·100^(n−2). Thus for k≥1 this grammar misses every count from 100^k through just below (7921/100)·100^k. In particular it misses 100^k+p for fixed positive p at all sufficiently large k. The fractal limit is a genuine structural recurrence, but its approach to the endpoint does not give an absolute placement.

## Every integer decimal exponent, and where its prime offsets can live

The [binary-apex seam](../target_core_selection_20261117/README.md) realizes the exact strong read

    be+F(b−1)F(e−1)−[b,e odd]       b,e≥2.

Here F is the binary inversion read from the flank decoder. Use b=2^r·5^a, e=2^(k−r)·5^(k−a), for a,r∈[0,k] and b,e≥2. Since F(2z+1)=F(z), and at least one arm is even, the count is

    10^k+F(5^a−1)F(5^(k−a)−1).

These are actual strong decimal-scale values for every k≥1, independent of how the powers of 2 are distributed. A prime offset requires a=0 or k: F(x)=1 only for a binary string of all ones (including x=0), and 5^j is not a power of 2 for j>0. The unique eligible offset is q_k=F(5^k−1), verified prime at k=1,2,10 (5,23,8573437). This asserts primality of the **offset**, not of the full decimal-centered count. Also F(5^j−1)≥(5^j−1)/2+1 for j≥1 by a first zero removal, so no fixed offset repeats at arbitrarily high k in this exact-product lane.

The normalized q_k/5^k has a nested 2-adic phase. If x=5^k−1 has low-to-high bits b_i, unrolling the binary recurrence gives

    q_k/5^k = Σ_(i≥0)(1−b_i)2^(−i−1) + O(k/5^k).

At a zero bit the exact summand before division is floor(x/2^(i+1)) plus the next bit; each rounding error is O(1), and the number of bits is O(k). The tail of the displayed series from the padded leading zeros is also O(1/5^k). The order of 5 modulo 2^m is 2^(m−2) for m≥2, because v₂(5^k−1)=2+v₂(k) by elementary squaring and the factor for odd exponents. Therefore as k varies, x runs uniformly over multiples of four modulo every 2^m. Reversing and complementing those low bits gives uniform dyadic subinterval frequencies, hence the normalized offsets are uniformly distributed on [3/4,1]. Conditionally on v₂(k)=s, the first s+2 complemented bits are ones and the next zero. The conditional uniform band is [1−2^(−s−2),1−2^(−s−3)], and the exponent class has density 2^(−s−1). This is the decimal phase recurrence, not an independence assumption.

Separate from those realized offsets, for any fixed prime p different from 2 and 5 a prime divisor ell of 10^r+p (ell coprime to 10) recurs on r modulo ord_ell(10). A finite list of factors gives periodic composite exclusions; for example 17 divides every 10^k+7 with k≡1 mod16, and 107 does so for k≡2 mod53. The divisor itself is the first exception in those two progressions. Combining such periodic tests does not prove infinitely many primes or physical placement for fixed p.

Replay: `python replay_decimal.py` checks the new literal equal-P₂ tiles, the prime example, the word formula by literal subset DP, and the finite and asymptotic ingredients' arithmetic controls.

There is an exact phase gap if both 5-power arms grow and the *total* decimal exponent has k≡2 mod4. Write a+b=k. If both arms are odd, their first phase bands force q_a q_b/5^k≤49/64+o(1) as min(a,b) grows. If they are both even, exactly one has 2-adic valuation 1 and the other has valuation at least 2, so q_a q_b/5^k≥105/128+o(1). Thus no accumulation in (49/64,105/128) is possible on this balanced-exponent lane. Fixing a short arm changes the problem and can bridge the gap. This is precisely the kind of correlation that independent phase models would miss.

## Why moving contacts alone cannot keep a fixed unit repair on the decimal stocks

The old two-exponential obstruction assumes fixed reads. A related obstruction here allows every tag within a repeated stock to change with height. Put j ordinal copies of C on the left of a fixed center D, and compare two fixed choices of apex contacts in D while holding the chosen left contacts common between them. Write L for the left weighted start and L^- for its signed part. From the adjacent interval law their exact finite difference is

    ΔH=11^j Δκ + L ΔR_D + L^- ΔR_D^-.

The signed left term L^- equals 11^(j−1) times the signed marked start on the final tile: a tag in any earlier tile leaves a nonempty optional suffix, whose parity sum cancels. If ΔR_D=0, ΔH is thus divisible by 11^(j−1). If ΔR_D≠0 and ΔH is a fixed constant, the formula forces L/11^j to stay bounded. A first tag in tile i gives L≥11^j(100/11)^(j−i) by the split injection and the ordinary suffix two-cover product. Therefore every tag must lie in the last O(1) tiles (or there are no tags). It follows that L and the whole difference are divisible by 11^(j−O(1)). In both cases an absolute nonzero constant gap fails at unbounded heights.

If both left and right decimal stacks grow, the conclusion is simpler: every difference term contains the ordinary H of at least one whole flank, since the center's own h,P are unchanged by its apex contacts, so the long cross-center term cancels. With stack lengths j,k this forces divisibility by 11^min(j,k). Thus the exact +2 repair cannot be exported to these decimal scales merely by programming the flanks. Changing the core or the way it reads the flanks is mathematically consequential.
