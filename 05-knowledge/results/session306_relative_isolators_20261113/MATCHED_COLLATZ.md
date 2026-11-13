# An ordinary minus Collatz edge in a fixed-order strong tournament chart

Take S305's crossing matched-core tournament X(a,b,c): a transitive order with four marked endpoints 0<1<2<3 and gaps of widths a,b,c, with only the edges 0–2 and 1–3 reversed. All widths below are positive. Its inherited count is

    H(a,b,c) = 1 + 2^(a+c+1)·3^b + 2(2^a+2^c)·2^b.

The two reverse edges close overlapping directed cycles; every gap lies in their strongly connected union.

## Exact edge

For u≥1, b≥1 the move

    X(u+1,b,u+1) → X(u,b+1,u+1)

is a physical ordinary-minus Collatz step on the count:

    3H(u+1,b,u+1)−1 = 2H(u,b+1,u+1).

The order, four named endpoints, and matched reversed arcs are unchanged. Move the last unmarked vertex of the first gap one place across endpoint 1, making it the first vertex in the second gap. The moved vertex and endpoint 1 are adjacent in the total order. Therefore the only arc changed is the arc between that vertex and endpoint 1. No other relative order or matched reversal changes. The ternary term triples then halves, because the central width increases by one while the outer total drops by one. The binary term does the same: the outer sum of powers changes from 2^(u+2) to 3·2^u, while increasing b doubles the term. The constants satisfy 3·1−1=2. Because the output count is odd, the valuation is **exactly one**, so this is an odd-to-odd edge of the unmodified minus map. Reversing the flip realizes the exact inverse arithmetic edge, not just two unrelated tournaments with the two target counts.

## Rigidity of the two-width palette

Suppose two charts X(a,b,c) and X(d,b+k,e) with fixed positive outer widths and fixed central offset k satisfy

    3H(a,b,c)+s = 2^v H(d,b+k,e),  s∈{−1,+1}, v≥1

for infinitely many positive integers b with b+k≥1. Since 1,2^b,3^b are independent (any nonzero linear combination has eventually the sign of its largest surviving base), coefficients must agree. The constants require (s,v)=(−1,1) or (+1,2). The ternary coefficients give k=1 and d+e=a+c−v by prime factorization (also for a negative offset, treating powers as rationals). Binary coefficients give

    3(2^a+2^c) = 2^(v+1)(2^d+2^e).

Sort a≤c and d≤e. If c−a≥2, the left side has four nonzero bits, impossible on the right. If c−a=1, its bits are at a and a+3, so d+e=2a+1−2v; the required outer total would be 2a+1−v, impossible. If a=c, the left bits are at a+1 and a+2, so {d,e}={a−v,a−v+1}; the outer-total condition forces v=1. This proves the displayed minus chart is the **only** infinite chart edge of the stated form, up to exchanging its outer widths, and explains why it cannot reenter: the output outer widths differ by one, while the input requires equality. It also rules out ordinary-plus chart edges with fixed outer widths, a fixed central shift, and fixed valuation within X. It is not a no-go theorem for nested matchings, other skeletons, nonfixed widths, or nonlinear composed returns.

## Design lesson

S305's 3H−2 width exchange and this 3H−1=2H' flip come from the same three positive modes, but place different budgets on the outer sum and product. The Collatz equation makes this visible coefficient by coefficient, and the paper's fixed-old-label view distinguishes an actual single flip from two scalar representatives. To chain inverse blocks, use a genuinely different marked core or a count-preserving change of chart. Repeating an exact numerical relation in the same two-width family cannot supply that change.

## Also allow the nested chart and zero outer widths

The companion nested matching N from S305 has coefficients, for b>0,

    H_N = A·3^b + (2A+1)·2^b + 1−A,   A=2^(a+c+1).

Use nonnegative outer widths for both charts. The leading coefficient still forces k=1 and A_src=2^v A_dst, hence A_src≥4. The binary coefficient condition is 3B_src=2^(v+1)B_dst. A nested source has odd B_src and is impossible. A crossing source with exactly one zero width has B_src=2(1+2^c) with c≥1 (up to exchanging widths); its valuation is one, also impossible. Both widths zero would give A_src=2, already excluded. Thus the source is crossing with positive widths and constant 1. A nested destination has constant 1−A_dst<0, and a crossing destination with at least one zero width has negative constant by the inherited S305 formula. Neither can match the positive source constant after applying 3H+s. This reduces to the positive crossing case, so the classification above holds over **both** matching charts with zero widths allowed. Central width zero is immaterial for infinitely many b with fixed shift.

This result concerns chartwise Collatz relations. It does not exclude arbitrarily large isolated matches with varying parameters, or multi-edge macros whose intermediate states leave the charts.

## No longer macros in the two matching charts

This last obstruction does not require intermediate tournaments to stay in the chart. Fix an arbitrary signed odd-step word of length ℓ≥2 with exact valuations summing to v≥ℓ. As a scalar map it has the form (3^ℓ H+C)/2^v, for a constant C independent of b. Suppose it relates two crossing or nested charts for infinitely many b, with fixed nonnegative outer widths and a fixed central shift. Comparison of the 3^b term forces central shift ℓ and outer sum drop v. The binary-coefficient comparison requires 3^ℓ B_src=2^(v+ℓ) B_dst. Hence a nested source (odd B_src) or a zero-width crossing source (valuation at most two) is impossible. For a crossing destination the comparison is

    3^ℓ (2^a+2^c) = 2^(v+ℓ)(2^d+2^e).

Sort a≤c and d≤e and put Δ=c−a, Δ′=e−d. If both are positive, comparison of 2-adic valuations and outer sums gives Δ′=Δ+v+2ℓ. Yet

    (1+2^Δ′)/(1+2^Δ) ≥ (1+2·8^ℓ)/3 > 3^ℓ,

contradicting the equality of odd parts. If Δ′=0 the right odd part is 1, impossible. For Δ=0<Δ′ the same comparisons give Δ′=v+2ℓ−2≥3ℓ−2; hence the required equality 3^ℓ=1+2^Δ′ is impossible for ℓ≥2. Thus no fixed signed macro reaches a crossing destination, even with zero outer widths. For a nested destination, let m be the smaller source width, let Δ be the difference, and note B_dst=1+2^(2m+Δ−v+2) is odd. If Δ>0 then m+1=v+ℓ, so the binary equation on odd parts becomes 3^ℓ(1+2^Δ)=1+2^(v+2ℓ+Δ), impossible by the same ratio inequality. If Δ=0 then m+2=v+ℓ, giving 3^ℓ=1+2^(v+2ℓ−2), also impossible. Thus neither destination type supports a fixed signed macro of length two or longer, even through arbitrary intermediate counts. The parameters (including the sign/valuation word) are fixed independently of b.
