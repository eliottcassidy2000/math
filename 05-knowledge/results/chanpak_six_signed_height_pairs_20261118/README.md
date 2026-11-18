# Signed height selectors for the six coefficients 4,2,5,3,7,13

The new result is on **absolute counts**, not congruence classes:

> For each a in {4,2,5,3,7,13} and integer k≥1, put p=a·2^k+1. There is a tournament with 7p Hamiltonian paths. There is also one with p paths except when p itself is 7 or 21. No primality assumption on p is needed.

As a companion, both signed dyadic rails 2^k±1 and 7(2^k±1) are classified exactly for k≥1, with only the known holes among those addresses. These are target-directed physical constructors; no large tournament enumeration is needed to generate an output. This does **not** cover arbitrary coprime targets.

## 1. The coprime-to-15 reduction, including the exceptional core

Let S be the attainable positive counts. Ordinal join multiplies Hamiltonian counts. Given odd N, write N=3^u5^v m with gcd(m,15)=1. If m≠7 and all such m have a witness, combine it with u copies of the directed triangle and v copies of the four-vertex strong count-five core. If m=7 and v≥1, use 35 and append the remaining factors; otherwise if u≥2 use 63 and append 3^(u−2). If neither condition holds the result is exactly 7 or 21. Literal 35 and 63 witnesses are in small.tsv. Thus proving the reduced coprime condition would suffice; it is not being assumed or proved by the construction below. In the opposite direction it is a necessary subcase, so the reduction is an equivalence. The earlier prefix certificate is [the background, not the basis of the infinite constructions](../chanpak_prime_frontier_selector_20261118/README.md).

## 2. A physical signed address, with the height read exactly

Use the [earlier prime-bearing wrapper](../chanpak_prime_frontier_selector_20261118/PRIME_RAIL.md). There is a small fixed tournament D, a chosen subset I⊆D, a singleton b, and a variable transitive module A of width w≥1. Orient b→A→D; orient I→b and b→D−I. Let d=H(D). Let G count ordered partitions D=P∪S into directed paths, with P nonempty ending in I and S allowed empty. Let Q count the same partitions with S empty or its first vertex outside I.

A spanning path has at most two A-runs: any A-run other than the first position of the whole path must start immediately after b. A single A-run gives d+G+Q templates (start at b, follow a nonempty prefix before b, or start in A and later cross through b). Two gives G templates: the first A-run is at the start, and the two A-runs are separated by a nonempty D-prefix followed by b. The transitive A fills the one slot uniquely, and its two ordered nonempty slots in 2^w−2 ways. Thus the literal count is

    H = G·2^w + β,       β = d+Q−G.

For an arbitrary nonempty module X the same count is H(X)·β+G·f_X(2), where f_X(2)=2H(X)+the ordered nonempty two-cover count. The 7-fold partner therefore can be exact on the *same variable X* by using **different** fixed cores and marked sets on the two sides: it is not the impossible operation of replacing a vertex by the entire lower tournament as an unchanged uniform module. Compare the [uniform 7-amplification no-go theorem](../seven_gap_atlas_20261117/README.md).

## 3. Twelve literal scaffolds and the finite starts

Rows and I bitmasks are constants in rails.py; the independent subset DP recomputes the three profile entries and whole-host instances at widths one and two. For p=a·2^k+1, use G=a·2^j on the left or G=7a·2^j on the right, and width w=k−j; β is 1 or 7, respectively.

| positive rail | G left | β left | j left | G right | β right | j right |
|---|---:|---:|---:|---:|---:|---:|
| a=1 (covers a=2,4 by shift) | direct | 1 | — | 896 | 7 | 7 |
| a=3 | 6 | 1 | 1 | 42 | 7 | 1 |
| a=5 | 20 | 1 | 2 | 35 | 7 | 0 |
| a=7 | 112 | 1 | 4 | 392 | 7 | 3 |
| a=13 | 208 | 1 | 4 | 91 | 7 | 0 |

The direct 2^k+1 construction reverses the long arc in the transitive tournament of order k+2. Its odd cycles all use that arc; an odd-size subset of k middle vertices gives each cycle, and two distinct cycles cannot coexist in one disjoint packing. Their positive contribution therefore gives exactly 1+2^k by [Grinberg–Stanley](https://arxiv.org/abs/2307.05569).

The right dyadic-plus row, for example, uses D=(2,132,9,35,15,23,63,125), I={1,2,5}, (d,G,Q)=(127,896,776). It gives 896·2^w+7=7(2^(w+7)+1). For the minus pair a=1 use G=1024,β=−1 on the left (j=10), and G=224,β=−7 on the right (j=5). Both sevenfold dyadic signs were already proved in earlier [positive](../chanpak_changed_contacts_20261117/PLUS_RAY.md) and [negative](../chanpak_changed_contacts_20261117/MERSENNE.md) scaffold work. We use independent witness profiles here and add the paired coefficient table and the explicit lower-minus model. In particular the positive 128 f(2)+H lift was already present even for arbitrary modules. The low values before a scaffold has positive width are all supplied explicitly in small.tsv from the earlier first-birth atlas. Its table of 26 distinct values is checked graph by graph. The only omitted requests are ones that reduce to the proved holes.

In particular for k≥4 there are exact tournaments for all four values 2^k−1,2^k+1,7(2^k−1),7(2^k+1). For the six coefficient word, the only holes on the lower positive rail are a=3,k=1 (7) and a=5,k=2 (21). On the upper rail even those two addresses are covered, as 49 and 147. The upper positive rail never asks for either forbidden count.

### The same incoming module gives both rails, with a nonconstant offset

Fix an arbitrary nonempty A. If one of the two rows for a given coefficient has smaller j, pad its incoming A by that many ordinal singletons. This keeps h=H(A) unchanged and multiplies its B=f_A(2) by the required power of two. The two resulting tournaments then have exactly

    p = α B + ε h,     7p = 7(α B + ε h),

where (α,ε) runs through (128,+1),(6,+1),(20,+1),(112,+1),(208,+1),(1024,−1). The paired function in rails.py takes the same arbitrary A and applies the indicated two different contexts. This does not use one output as a module of the other. In particular the additional coefficients are not just six fixed signed rays: their common incoming H(A) can vary as well.

For a prescribed numerator p the profile test is now explicit. At |A|≥2, B≡0 mod 4 by complementary bipartition and tournament oddness, so h≡εp mod 4α. Height requires B=(p−εh)/α exactly. Every assignment of vertices to two slots has a witness (at least one path on each side), so B≥2^|A|, and the two one-slot choices give B≥2h. Thus |A|≤floor(log₂[p/(α−1/2)]). This is a finite architectural test, not a proof that the required h,B pairs exist. Unit checks at 3 and 5 can be imposed before the profile search; for example with α=6, p inherits h modulo 3.

## 4. Outside the finite prime bank

Examples actually checked for primality, rather than just assuming the exponential expression is prime:

| prime p | 7p | witness orders |
|---:|---:|---:|
| 206,158,430,209 = 3·2^36+1 | 1,443,109,011,463 | 39,41 |
| 2,748,779,069,441 = 5·2^39+1 | 19,241,453,486,087 | 43,46 |
| 824,633,720,837 = 192·2^32+5 | 5,772,436,045,859 | 40,42 |
| 2,305,843,009,213,693,951 = 2^61−1 | 16,140,901,064,495,857,657 | 60,65 |

The first three primality claims use complete trial division through the square root; the third uses the shared count-five four-vertex module with 32 ordinal singletons, not a unit-offset ray. The last uses the [Lucas–Lehmer criterion](https://www.mersenne.org/various/math.php) with exponent 61 and ending residue zero. These are explicit construction orders, not claimed minimums or claims of new prime discoveries.

Input N on one of the specified rails determines its height rather than authorizing an overshoot. Divide by the requested scale (1 or 7), choose the indicated sign, and factor the shifted integer into a power of 2 times its odd part. The odd part selects a coefficient family and the exact valuation supplies k (the even coefficients 2,4 only shift it). This is a height-and-mantissa rule, not a residue-class fit. For arbitrary odd p one can still choose the sign with the larger v₂(p±1), guaranteeing at least two binary levels. The missing uniform theorem would force a suitable fixed or target-dependent signed scaffold G,β from its resulting odd part **within that height budget**. If the incoming module is only transitive, a finite list of fixed signed coefficients gives at most O(log X) distinct values up to X and cannot do that for all coprime p. The arbitrary-module pairing is not subject to that particular ray bound, but still lacks a target-to-profile theorem. No general constructive closure is claimed.

## 5. What transfers from the archived six numbers

The archived [Bernoulli/modular genealogy](../borwein_ladder_owner_numbers_vsc_twelve_20261028.md) identifies the unordered set as {d+1 : d|12}; its prime members give the von Staudt–Clausen denominator 2730 of B12. The Monster/Baby Monster reading of the *ordered* word is a separate mnemonic and does not explain the graph scaffold. In the [old metagraph eta-product observer](../joint_metagraph_dyadic_20261113/README.md), these are exactly the levels N for which the normalized square-root eta product has integral coefficients. More generally, for

    F_N(q)=product_(N does not divide d)(1−q^d)^(24/(N−1)),

its normalized kth root is integral iff k(N−1)|24, from the first q coefficient for necessity, and the integer-exponent product for sufficiency. For such a root put A=24/[k(N−1)]. At a prime ℓ with s=v_ℓ(A), Frobenius shifts the first nonconstant term to **exact** height ℓ^s, with nonzero coefficient −A/ℓ^s mod ℓ. For the six square-root levels in the user's order the first mod-2 heights are 4,4,1,2,2,1; mod-3 heights are 1,3,3,3,1,1. This is an older style of height-versus-residue mechanism, not an eta-to-tournament map.

As six independent clock periods, the [exact cover calculation](../harmonic_reader_residual_20261116/WHEELS.md) leaves at least 576 of 5460 phase addresses uncovered. The new graph theorem does something different: it uses the same six numbers as mantissas of proposed target values, and supplies literal paired reads at every height on those rays. No inference from a clock density or number coincidence was used in the graph proof.

## Replay

`python rails.py | diff -u replay.out -` checks the small literal graphs and prime examples. Use `python rails.py 3 36 --scale 7` or `python rails.py 1 61 --sign -1 --scale 7` to materialize the outgoing masks. The optional discovery program pairs.cpp reads an input table of small D tournaments (order, H, outgoing masks, no header) and searches for 7-scaled signed coefficients; the theorem only needs the rows in rails.py, not completeness of any discovery bank.
