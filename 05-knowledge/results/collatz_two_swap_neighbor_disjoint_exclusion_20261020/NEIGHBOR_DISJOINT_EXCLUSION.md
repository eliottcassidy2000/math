# All-length exclusion for same-direction disjoint edits of neighboring ones

20 October 2026. Status: **PROVED**, with separately audited content, cited analytic dependencies, and independently verified finite-exact certificates. See [the independent audit](audit/NEIGHBOR_DISJOINT_EXCLUSION_AUDIT.md). No claim about arbitrary disjoint two edits is made.

## 1. The exact family

Let p>=2, 1<=r<p, gcd(p,r)=1, q=p+r, and start from the zero-intercept lower mechanical word. Put n=p-r and use the source-one residue labels a_j=-jr modulo p. Perform either pair of disjoint, same-direction actual adjacent swaps

    plus:  move residues u,u+n right, 0<=u<r-n;
    minus: move residues u,u+n left,  n<=u<r.

These are precisely the same-direction disjoint edits of two consecutive source ones in this primitive k=1 universe. Both orders are legal and have the same output. The complete legality and actual-content identity are proved in the [neighbor-disjoint content theorem](../collatz_two_swap_neighbor_disjoint_content_20261020/NEIGHBOR_DISJOINT_CONTENT.md).

**Theorem.** No output in this family is an integer-cycle branch word, for either clock sign and at any length.

Unlike the overlapping theorem, this family has no p=2,q=3 boundary paths at all. The theorem excludes only this specified disjoint subclass. It does not exclude disjoint edits of nonneighboring ones or mixed directions, nor does it silently include repeated mechanical sources.

## 2. Reduction to p<8000

For an output w write

    D=2^q-3^p, B=B(w), g=gcd(|D|,B).

The inherited content theorem gives g<=121^(4 p^(2/3)). For p>=8000=20^3, use the exact inequalities

    log121<29/6, log3>13/12.

The first is certified by the first 21 nonnegative Taylor terms for exp(29/6), whose rational sum exceeds 121; the second follows from the first two positive terms of 2 artanh(1/2). Hence

    log(g/3^p)<(29/30-13/12)p=-7p/60<-p/9.          (2.1)

If the output were integral, then |D|=g, so the nonzero multiplicative form U=2^q3^(-p)-1 would satisfy |U|<exp(-p/9). Also |U|<1/2, so

    |q log2-p log3|<2 exp(-p/9).                    (2.2)

The [audited absolute-gap arithmetic certificate](../collatz_two_swap_endpoint_exclusion_20261020/ENDPOINT_TWO_EDIT_CLASSIFICATION.md) excludes (2.1) for every coprime p>=100 with p<q<2p. Its proof uses the stated explicit real Matveev theorem on U, obtains p<10^12, applies Legendre's criterion, and excludes every one of the 18 remaining convergents by exact rational interval bounds for both signs. Those steps depend on the gap inequality and clock range, not on the endpoint edit that first supplied it. Therefore every possible integral output in the present family has

    2<=p<8000.                                     (2.3)

This explicitly retains the cited Matveev and exact interval dependencies. No numerical estimate for an endpoint resultant is being transferred to the present family.

## 3. Independently checkable signed descent through 2^26

Put N=2^26=67,108,864. The new `verify_signed_descent.cpp` is the earlier checked unsigned-128-bit producer with the single change of its inclusive magnitude limit from 2^24 to 2^26. It proves:

* Every positive starting magnitude 1<=m<=N under the half-step map enters the cycle 1,2.
* Every negative starting magnitude 1<=m<=N, under m/2 for even m and (3m-1)/2 for odd m, enters one of the three magnitude cycles

      1;
      5,7,10;
      17,25,37,55,82,41,61,91,136,68,34.

It verifies the displayed cycle edges first. Every other odd start must reach a strictly smaller positive magnitude or a verified cycle minimum; even starts halve immediately. Strong induction proves the finite proposition. Trajectory arithmetic is checked before every multiplication/addition; overflow, a one-million-step cap, bad cycle edges, or incomplete coverage cause a fatal nonzero exit. No trajectory is assumed to stay below N.

The optimized and unoptimized runs agree exactly:

| Map | Odd non-base starts | Even starts covered | Accelerated odd steps | Largest descent length |
|---|---:|---:|---:|---:|
| positive | 33,554,431 | 33,554,432 | 117,191,811 | 237 |
| negative magnitude | 33,554,429 | 33,554,432 | 117,183,239 | 174 |

One positive and three negative known minima complete all N starts in each map. The largest prehalving numerators are 474,637,698,851,092 and 254,287,025,337,584, respectively.

This is a new executed finite certificate, not an appeal to an external global verification bound. It does not prove convergence for all signed integers.

## 4. The complete small-weight reduction

Any integer cycle other than the four verified cycles has an integer minimum absolute value M>N, hence M>=N+1. Multiplication over its p odd entries gives

    2^q=product(3+1/a_i) positive,
    2^q=product(3-1/a_i) negative,

where every a_i>=M. For alpha=log3/log2 this yields

    |q/p-alpha|<1/((3M-1)log2)
                <3/[2(3N+2)]<1/(2N),

using log2>2/3 and -log(1-t)<t/(1-t) for the negative sign. Since p<8000<8192 and N=8192^2, the error is strictly below 1/(2p^2). With gcd(p,q)=1, Legendre forces q/p to be an ordinary continued-fraction convergent of alpha.

The exact rational logarithm checker gives the complete list below weight 8000:

    (p,q)=(2,3),(5,8),(12,19),(41,65),
          (53,84),(306,485),(665,1054).

The next denominator is 15601, from 24727/15601, beyond the residual range. Eighty-term positive artanh-series bounds certify every continued-fraction digit, with no floating-point decision.

The known cycles must be handled separately, not assigned the false minimum premise M>N. Their shapes are (1,2),(1,1),(2,3),(7,11). Coprimality of p,q excludes multiple traversals. The first two are outside p>=2 and p<q<2p. The shape (2,3) is already on the list and has no paths of the present type. The -17 cycle has shape (7,11), not on this convergent list, and is added explicitly.

## 5. Every remaining actual path

`verify_small_neighbor_candidates.py` performs both orders of every legal disjoint neighbor pair in the eight shapes, including the original cut. It independently enumerates actual same-direction edge pairs moving consecutive source labels and checks equality with the residue-generated universe. Two carry algorithms agree, and full integrality is checked by gcd(B,|D|).

| p | q | Ordered paths | Maximum actual content |
|---:|---:|---:|---:|
| 2 | 3 | 0 | not applicable |
| 5 | 8 | 4 | 1 |
| 7 | 11 | 4 | 1 |
| 12 | 19 | 8 | 311 |
| 41 | 65 | 28 | 19 |
| 53 | 84 | 36 | 11 |
| 306 | 485 | 208 | 1 |
| 665 | 1054 | 452 | 37 |
| total | | 740 | |

Every one of the 740 ordered paths is nonintegral. The partial factor 311 at p=12 is retained, rather than confused with the full clock magnitude 7153. The four paths at the known -17 shape also fail integrality directly. Normal and optimized Python produce identical receipts with 2,996 explicit checks.

Sections 2 through 5 exhaust all possible integer outputs of the stated family and prove the theorem. The independent audit verified the complete proof application and both finite reductions.

## 6. Reproduction and scope

From this directory:

    g++ -std=c++17 -O2 -Wall -Wextra verify_signed_descent.cpp -o /tmp/neighbor-signed-descent
    /tmp/neighbor-signed-descent > descent.replay.json
    cmp descent.receipt.json descent.replay.json
    python3 verify_small_neighbor_candidates.py > candidates.replay.json
    python3 -O verify_small_neighbor_candidates.py > candidates.optimized.replay.json
    cmp candidates.receipt.json candidates.replay.json
    cmp candidates.receipt.json candidates.optimized.replay.json

The content theorem and its prime-7 argument remain a separate theorem. The present closure additionally requires the cited Matveev/continued-fraction arithmetic certificate, the new finite signed descent through 2^26, and the complete 740-path residual check. No arbitrary disjoint-edit theorem or general Collatz resolution follows.
