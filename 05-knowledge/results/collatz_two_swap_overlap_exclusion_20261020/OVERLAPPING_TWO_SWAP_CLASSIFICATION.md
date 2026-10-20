# All-length integer classification of overlapping mechanical two-swap paths

20 October 2026. Status: **PROVED**, with explicitly cited analytic dependencies and **FINITE-EXACT** computational dependencies. Both new finite computations and the complete proof application have passed a materially different independent audit; see [the audit report](audit/SIGNED_DESCENT_AND_SMALL_WEIGHT_AUDIT.md).

## 1. Precise result

Start from the nonconstant zero-intercept lower mechanical word

    c_i=floor((i+1)p/q)-floor(ip/q),   0<p<q.

An actual adjacent swap exchanges unequal neighbouring bits on the circle. A two-swap path is called **overlapping** when its two affected bit sets intersect. An **undo** returns the original word by using the same physical edge twice. At q=2, the two directed edge indices describe the same physical edge and count as undo, not as a distinct nonundo path.

**Primitive two-one theorem.** Suppose gcd(p,q)=1, q=p+r with 1<=r<p, and apply either legal overlapping two-distinct-one chain:

    plus: 110 -> 101 -> 011,
    minus: 011 -> 101 -> 110.

The output is an integer-cycle branch word if and only if (p,q)=(2,3). In the ordinary cut the two exceptional paths are

    011 -> 110 -> 101, fixed point -7;
    011 -> 101 -> 110, fixed point -5.

Both are rotations of the known negative cycle -5,-7,-10. The theorem includes both clock signs and every cyclic cut.

**Overlapping-path corollary.** Among all nonconstant mechanical sources, not necessarily primitive, the only integer outputs of a nonundo overlapping two-swap path are those two p=2,q=3 outputs and their rotations. In particular, none has exact cyclic adjacent-swap distance two from the balanced source necklace: the exceptional outputs are balanced rotations themselves.

This does not exclude arbitrary disjoint two-swap paths. Combined with the previously proved repetition classification, the remaining genuine two-edit problem is primitive disjoint edits.

## 2. Inherited large-weight reduction

For the primitive k=1 source q=p+r, label ones by residues a_j=-jr modulo p. The complete legal parameter ranges are

    plus:  0<=u<min(r,p-r), move u then u+r right;
    minus: max(0,p-2r)<=u<p-r, move u+r then u left.

These ranges are proved in the [overlap-content theorem](../collatz_two_swap_overlap_lattice_content_20261020/OVERLAPPING_TWO_ONE_CONTENT.md). For either output, put

    D=2^q-3^p,
    B=sum_{i:w_i=1} 2^i 3^(number of ones after i),
    g=gcd(|D|,B).

Integrality is equivalent to g=|D|, with fixed point B/D and sign equal to the sign of D.

The cited packet proves

    g_plus<=38^(4p^(2/3)),
    g_minus<=46^(4p^(2/3)).

It retains the actual common clock root, uses a prime-5 valuation to prove noncancellation of two independent relation polynomials, and applies Minkowski's second theorem to their exponent lattice. Its separate Matveev/continued-fraction corollary excludes every integral path with p>=4096. Therefore it remains to classify

    2<=p<4096, all coprime 1<=r<p, all legal u.

The previous packet counted this unreduced universe but did not check its nearly seven billion paths. The following new finite magnitude certificate replaces that pathwise enumeration by a complete reduced candidate argument.

## 3. A self-contained signed magnitude certificate

Let N=2^24=16,777,216. The new computation proves:

**Finite signed-descent proposition.**

* Every positive starting magnitude 1<=n<=N under T_+(n)=n/2 for even n and (3n+1)/2 for odd n enters the cycle 1,2.
* Every starting magnitude 1<=n<=N under T_-(n)=n/2 for even n and (3n-1)/2 for odd n enters one of the three magnitude cycles

      1;
      5,7,10;
      17,25,37,55,82,41,61,91,136,68,34.

The second map is exactly the original map on negative integers after taking absolute values. No externally reported global Collatz verification range is used.

### Why the finite algorithm proves the proposition

The producer `verify_signed_descent.cpp` first verifies every edge of the four listed cycles using the full half-step maps. It then tests every odd start n<=N outside the explicitly verified cycle minima.

For an odd value x it calculates 3x+1 or 3x-1, then divides by every factor of 2. These accelerated odd transitions are exact compositions of the half-step map. Iteration stops only on reaching a smaller positive value than the starting n or one of the verified cycle minima. Strong induction on the starting magnitude proves the claim. Even starts halve immediately to a smaller magnitude and need no trajectory loop.

Every mathematical trajectory operation uses checked unsigned 128-bit arithmetic. The program exits unsuccessfully on a possible overflow, a per-start cap of one million odd steps, an invalid cycle edge, or incomplete coverage. A failure is not discarded or counted as a completed case. The optional modulo-2^64 trace checksum is separate and has no mathematical role.

Both optimized and unoptimized builds completed successfully with identical receipts:

| map | non-base odd starts tested | even starts covered | total accelerated odd steps | largest descent length |
|---|---:|---:|---:|---:|
| positive | 8,388,607 | 8,388,608 | 29,285,622 | 181 |
| negative magnitude | 8,388,605 | 8,388,608 | 29,303,560 | 172 |

There is one positive base minimum and three negative base minima, completing all N starts in each map. The largest prehalving numerators encountered are respectively 60,342,610,919,632 and 30,541,433,029,400, far below the checked arithmetic limit.

A materially different independent implementation uses the full half-step maps, checked 64-bit arithmetic, and strong-induction descent. It does not import the producer implementation or use accelerated odd transitions. Its complete positive and negative runs take respectively 58,576,417 and 58,604,988 half-steps; the longest descents are 287 and 273. Both verify the same finite signed proposition.

This is a finite theorem through N, not a claim that all integers converge to the listed cycles.

## 4. Minimum magnitude forces a small convergent

Consider any integer cycle other than the four explicitly listed cycles. By Section 3, its minimum absolute value M is an integer satisfying

    M>N, hence M>=N+1.                               (4.1)

Let its p odd steps have absolute values a_1,...,a_p, all at least M. Multiplying step ratios around the cycle gives

    2^q=product_i(3+1/a_i) for a positive cycle,
    2^q=product_i(3-1/a_i) for a negative cycle.

Put alpha=log3/log2. For the positive sign,

    0<q/p-alpha<1/(3M log2).

For the negative sign, use -log(1-t)<t/(1-t) for 0<t<1:

    0<alpha-q/p<1/((3M-1)log2).

The elementary series log2=2 artanh(1/3) gives log2>2/3. By (4.1), the larger of the two bounds is strictly below

    3/[2(3N+2)] < 1/(2N).                            (4.2)

For p<4096, p^2<N=4096^2, so (4.2) is strictly less than 1/(2p^2). Since gcd(p,q)=1 in the primitive family, Legendre's criterion forces q/p to be an ordinary continued-fraction convergent of alpha.

It is important that the known negative cycles are not fed into this lower-minimum argument. In particular, the -17 cycle has (p,q)=(7,11), whose ratio is not one of these convergents.

## 5. Complete candidate list and direct actual swaps

`verify_small_convergent_overlaps.py` encloses log2 and log3 by exact rational positive-series bounds with 80 terms. At every continued-fraction step, both interval endpoints have the same integer part. It extracts all convergents with 2<=p<4096 and p<q<2p:

| p | q | sign of D |
|---:|---:|:---|
| 2 | 3 | negative |
| 5 | 8 | positive |
| 12 | 19 | negative |
| 41 | 65 | positive |
| 53 | 84 | negative |
| 306 | 485 | positive |
| 665 | 1054 | negative |

The next convergent is 24727/15601, whose denominator exceeds the entire residual range. The excluded denominator-one convergents cannot be a k=1 two-distinct-one source with p>=2.

The four known integer cycles have shapes

    positive 1,2: p=1,q=2;
    negative -1: p=1,q=1;
    negative -5 cycle: p=2,q=3;
    negative -17 cycle: p=7,q=11.

A word traversing a known cycle multiple times has a common factor in its p and q. Since the present source shape is coprime, only one traversal can occur. Thus only (2,3) and (7,11) require direct checks in the k=1 family. The former is already on the convergent list; the latter is added explicitly.

The eight resulting shapes contain exactly 906 legal ordered overlapping paths. The checker:

1. constructs the original mechanical word and both exact legal parameter ranges;
2. performs every swap on the current bit word, including the original cut;
3. independently enumerates actual neighbouring-edge paths and checks completeness;
4. computes B both by affine composition and by the positional power sum;
5. tests full gcd content, not a norm or a prime-only surrogate;
6. directly verifies every parity and closure in each integral output.

Its 3,663 explicit checks pass in normal and optimized Python, with identical receipts. All paths at p>=3 are nonintegral. The only outputs with g=|D| are

    p=2,q=3, plus: 101, B=7,D=-1, fixed point -7;
    p=2,q=3, minus:110, B=5,D=-1, fixed point -5.

At the hostile -17-cycle shape (7,11), all six overlapping paths have content 1 while |D|=139. This direct check is consistent with the earlier independently proved necklace distance three of the genuine -17 cycle, but does not need that distance result as an extra dependency.

Together with the complete signed minimum argument and the inherited p>=4096 exclusion, this proves the primitive two-one theorem.

## 6. Other overlapping paths and nonprimitive sources

The preceding result has a narrow k=1 two-distinct-one hypothesis. The following separate arguments justify the broader overlapping-path corollary.

A two-one chain requires adjacent source ones, hence a mechanical gap of one. Since mechanical gaps are floor(q/p) or ceil(q/p), this forces q<2p. Thus every primitive mechanical two-one overlap belongs to the k=1 family already classified above.

### A single labelled one moving twice

A nonundo overlapping path of this type is 100 -> 010 -> 001 or its reverse. The source has a gap between successive cyclic ones of at least three. For a mechanical source, gaps are floor(q/p) or ceil(q/p), so a gap at least three forces q>2p.

Every positive integer cycle satisfies

    2^q=product_i(3+1/a_i)<=4^p,

since each positive odd a_i>=1. Hence q<=2p. A negative integer fixed point requires D<0, hence q<p log_2(3)<2p. Thus the single-one two-step path is nonintegral for either sign. The q=2 coincident-edge path is undo and was excluded by definition.

### Repeated mechanical sources

If d=gcd(p,q)>=2, the mechanical source is the d-th power of a nonconstant primitive source v of shape p_0,q_0. The [audited repetition classification](../collatz_two_swap_repetition_classification_20261020/TWO_EDIT_REPETITION_CLASSIFICATION.md) gives the common factor

    H=sum_{j=0}^{d-1} (2^q_0)^(d-1-j)(3^p_0)^j>=7,

which divides both the source carry and D. Rotate the cut outside the three-bit affected block; rotation preserves gcd(H,B) because H divides D and is coprime to 6. In this cut a nonundo overlap changes the carry by a power-of-2-and-3 unit times 3 or times 5. Since H is coprime to 6 and exceeds 5, it cannot divide that perturbation. Therefore it cannot divide the edited carry, and an integer output is impossible.

The mirrored d=2 exception in the general two-edit repetition theorem uses edges separated by q_0>=2; those affected bit sets are disjoint. It is not an overlapping exception.

Finally, every two actual swaps with intersecting affected sets are either undo, a single-one two-step path, or a two-one chain, once the q=2 edge convention is handled. The three cases above prove the overlapping-path corollary, without asserting anything about disjoint edits.

## 7. Reproduction, dependencies, and remaining scope

The proof uses the frozen subexponential overlap-content theorem and its separately cited Matveev/continued-fraction large-weight reduction. The new finite dependencies are entirely supplied here:

    g++ -std=c++17 -O2 -Wall -Wextra verify_signed_descent.cpp -o /tmp/signed-descent
    /tmp/signed-descent > descent.replay.json
    g++ -std=c++17 -O0 -Wall -Wextra verify_signed_descent.cpp -o /tmp/signed-descent-o0
    /tmp/signed-descent-o0 > descent.unoptimized.replay.json
    cmp descent.receipt.json descent.replay.json
    cmp descent.receipt.json descent.unoptimized.replay.json
    python3 verify_small_convergent_overlaps.py > candidates.replay.json
    python3 -O verify_small_convergent_overlaps.py > candidates.optimized.replay.json
    cmp candidates.receipt.json candidates.replay.json
    cmp candidates.receipt.json candidates.optimized.replay.json

The independent code, overflow/termination controls, exact candidate checks, and audit receipts are in `audit/`. The source and receipt hashes are frozen in the packet manifest.

No seven-billion-path enumeration was performed or claimed. Its completeness obligation is discharged instead by the independently verified bounded-magnitude proposition, the signed product inequalities, exact convergent exhaustion, and direct checks of every remaining actual path.

The result is an all-length classification of overlapping two-swap outputs from mechanical sources. It is not a proof of the Collatz conjecture, not an exclusion of arbitrary two edits, and not a claim of historical novelty. Primitive disjoint two-edit patterns remain outside this theorem.
