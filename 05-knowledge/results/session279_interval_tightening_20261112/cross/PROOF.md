# Ground-state moment compression: remove the ubiquitous error first

**Status: PROVED identities and interval theorems; FINITE-EXACT rational
experiments. A certificate compression result, not a faster exact enumeration
algorithm or a proof of Chan–Pak.**

## Inheritance and the less-obvious connection

The direct input is Session278's [divisibility-defect gas](../../session278_signed_corrections_20261112/enumeration/PROOF.md).
For each odd cycle partition, its Burnside contribution is

    W product_e (1-q_e),       q_e=1-2^(-d_e), d_e>0.

Some cycle lengths may already be integrated exactly in a backbone B. The
remaining d_e are the differences `min(a,b)-gcd(a,b)` for incompatible pairs
outside B. Singleton cycles never contribute such edges.

A different, older source is the LRC reversed-peel moment hierarchy:
[THM-2210, adaptive moment LP](../../../../01-canon/theorems/THM-2210-nested-binomial-minorant-and-adaptive-moment-lp-hierarchy.md).
It already proves the exact LP/dual interpretation of factorial-moment bounds.
[THM-2237, truncated Boolean moments](../../../../01-canon/theorems/THM-2237-truncated-boolean-moment-interval-and-parity-top-atom-majorants.md)
already identifies parity as an ambiguity invisible to proper moments.
Neither theorem is new here. The new application is to **shift the defect
energy before optimizing those moments**, and retain the number of edges
whose energy really changed. Ordinary Bonferroni expands probabilities near
one; the shift removes a large uniform cost exactly and makes most remaining
coordinates zero.

The connection has a typed ledger. Source: LRC multiplicity laws. Target:
weighted artificial Bernoulli laws of arithmetic Burnside defects. Map: the
factorial-moment cone and its polynomial dual. Preserved: zero-atom bounds.
Lost: cycle lengths, divisibility graph and geometric LRC data. Required
sidecars here: baseline weight, minimum positive defect and exceptional-edge
count. No tournament path-count realization follows from this counting map.

## 1. Exact minimum-defect extraction

Fix one sector with M positive defects. For M>0 set

    delta=min_e d_e,       rho=2^(-delta),
    z_e=1-2^(-(d_e-delta)),
    R=#{e:d_e>delta}.

Then

    W product_e (1-q_e) = W rho^M product_e (1-z_e).       (1)

Exactly M-R of the z's are zero. The probabilities of the R other coordinates
are strictly between zero and one. The common contribution `delta M` is removed
from the energy before any approximation. For M=0 use delta=0, rho=1, R=0.
All equal-defect sectors become exact at degree zero, regardless of M.

This is not a claim that the minimum defect has a universal fixed value. For
example the decisive four-cycle core `(13,15,17,25)` has minimum defect10,
not12. The cell carries its actual minimum.

## 2. No-reference conversion from the old low moments

Suppose a producer gives, for each cell `(M,delta,R)`, the five numbers

    C_j=sum_(sectors in cell) W e_j(q),        j=0,...,4.   (2)

The consumer need not see the partitions, exact count, full defect list or
higher corrections. Put `a=1-1/rho`, `b=1/rho`. Since `z_e=a+b q_e`,

    G_k=rho^M sum_(j=0)^k binom(M-j,k-j) a^(k-j)b^j C_j.  (3)

For k>M the value is zero. This is a finite exact affine elementary-symmetric
transform: first choose which j of the k selected coordinates contribute the
`b q` term, then choose the remaining k-j positions. Thus

    G_k=sum_(sectors in cell) W rho^M e_k(z).

In particular `G_k=0` for k>R. The negative coefficients in (3) are harmless
when evaluated as exact rationals; floating cancellation would not be a
certificate.

Aggregate cells with equal R. Let X count independent active exceptional
coordinates in each sector, and mix those laws with their positive weights
`W rho^M`. Then `G_k=sum_x binom(x,k) p_x` for a nonnegative unnormalized law
on `{0,...,R}`, and its zero mass is exactly the desired contribution. This
weighted probability law is an analytical construction, not a probability
distribution on tournament isomorphism classes.

For comparison, forgetting R and keeping M yields a valid but weaker support
bound. A000568 at n64 with B={3} needs422 low residue bits under that support
bound, but362 when R is retained. The coefficient order is still four.

## 3. Exact degree-four polynomial interval in each scalar stratum

The following is a specialization of the inherited LP dual, made explicit so
no floating solver is used. For a nonnegative law on `{0,...,r}`, known through
moments `G_0,...,G_4`, write

    Q_S(x)=product_(s in S)(s-x)/s.

For r>4, the upper degree-four vertices have root sets

    S=(i,i+1,j,j+1),        1<=i, i+1<j, j+1<=r.           (4)

Every such polynomial equals one at zero and is nonnegative at every positive
integer: a negative interval lies only between two adjacent roots. Lower
vertices with Q(0)=1 have

    S=(1,i,i+1,r),          2<=i, i+1<r.                  (5)

They are nonpositive at all positive support nodes. The zero polynomial is
also a valid lower bound. Include lower-degree feasible polynomials if desired;
they cannot worsen the interval. The implementation includes all such
endpoint/adjacent-pair certificates through degree four.

To see completeness of (4)-(5), fix Q(0)=1. This leaves four coefficients.
At a vertex four distinct positive evaluation constraints bind, hence Q has
those four roots. For an upper polynomial the roots must pair adjacently;
otherwise an integer lies in a negative interval. For a lower polynomial the
outer positive intervals must be cut off by roots1 and r, with an adjacent
middle pair. A nonzero lower bound with positive objective can be scaled to
Q(0)=1; otherwise the zero lower bound is at least as good. An upper bound
with Q(0)>1 can be scaled down to one. The finite LP has an attained optimum;
its constraint matrix has no lineality on polynomials with fixed constant
when r>=4. This gives all optimal possibilities. At r<=4, ordinary binomial
inversion gives the exact zero atom directly.

If `Q(x)=sum_k alpha_k binom(x,k)`, its integral is `sum_k alpha_k G_k`.
The coefficients are computed by exact forward differences of Q at0,...,4.
Taking the maximum lower and minimum upper values is therefore the **optimal
scalar degree-four moment interval for that R stratum**. Sum stratum intervals.
This is not asserted optimal among every conceivable use of the original
arithmetic cell information: aggregation and the general moment cone forget
constraints of the dyadic Bernoulli sources.

Finally intersect with the integers by ceiling the lower endpoint and flooring
the upper. The consumer reports an exact value only if those integers coincide;
otherwise it returns the interval, not a guessed nearest integer. Its width w
requires at most `bit_length(w)` independently certified low residue bits for
unique recovery, as in Session277. No residue is manufactured by this consumer.

## 4. A structural exactness band, and its sharp first failure

**Theorem.** With backbone B={3,5,7,9,11}, the degree-four ground-moment
certificate is rationally exact at every order n<=69. The first residual
rank-five core occurs at70, and the moment interval then has positive width.

Proof. All remaining cycles have length at least13. Let t be their number.
For t<=3 there are at most3 defect edges. For t=4 with a repeated length,
there are at most5 edges, and extracting the minimum removes at least one,
so R<=4. For four distinct odd lengths of total at most69, the only possibilities
are `(13,15,17,19)`, `(13,15,17,21)`, `(13,15,17,23)` and `(13,15,19,21)`.
In each the three edges from13 have defect12, which is the minimum; thus R<=3.
For t=5, the only multisets of total at most69 are `13^5`, `13^4 15`,
`13^4 17`, and `13^3 15^2`. Their positive defects are all12, so R=0.
Six remaining cycles cost at least78. Adding backbone cycles or singleton cycles
cannot alter the residual defect list. This proves the all-n<=69 statement.

At70, `(13,15,17,25)` has defects

    12,12,12,14,10,16.

After extracting delta10, its five nonzero z values are

    3/4,3/4,3/4,15/16,63/64.

Their independent active-count law has positive mass at all six nodes0,...,5.
Four moments cannot pin its zero mass in the general nonnegative moment cone;
the explicit parity perturbation in Section5 proves strict width. All other
sectors can be held fixed, so the complete moment interval is also nontrivial.
The exact experiment at70 requires31 low residue bits. This is a sharp boundary
for this certificate, not for efficient enumeration of A000568 itself.

For B={3}, a separate finite exact core enumeration gives the earlier boundary39,
with core `(5,7,7,9,11)`. The first R>=5 cores are48 for B={3,5} and for
B={3,5,7}. These latter threshold claims are explicitly FINITE-EXACT controls;
the combinatorial proof above concerns B through11.

## 5. No fixed backbone and moment depth repairs absolute accuracy forever

The positive-sector obstruction of Session277 and fixed-Bonferroni obstruction
of Session278 have a stronger descendant: even **optimal scalar moment LPs**
at fixed depth after this minimum shift retain unbounded absolute ambiguity.
This is a limit on that certificate family, not on all methods using arithmetic.

Fix a finite backbone B and depth d. Choose odd primes a<b<c outside B and
put r=d+1. The core `(a,b^r,c)` has r edges a-b, one edge a-c and r edges b-c.
Their minimum positive defect is a-1; its extraction leaves exactly r equal
positive probabilities `z=1-2^(-(b-a))`. The resulting binomial law nu on
0,...,r has every atom positive. For

    v_k=(-1)^k binom(r,k),
    epsilon=min_k nu_k/binom(r,k)>0,

the laws `nu+epsilon v` and `nu-epsilon v` are nonnegative, have identical
factorial moments through degree d, and have zero masses differing by
`2 epsilon`. The identity `sum_k binom(k,j)v_k=0` for j<r is the r-th finite
difference of a degree-j polynomial, the inherited Boolean-parity mechanism.

Pad the physical cycle partition with f singleton cycles. Its positive
shifted weight is a fixed constant times

    2^[binom(f,2)+f(r+2)]/f!,

which tends to infinity. Adding all other sectors unchanged gives two feasible
laws for the aggregate retained moments whose zero masses differ by that
unbounded factor. Therefore the optimal scalar moment interval cannot have
width below one eventually at fixed B,d. It must adapt its backbone, depth or
other information. No implication about missing Hamiltonian counts is made.

For the order70 core above, epsilon=1/65536. Its shifted padded weight is

    2^[binom(f,2)+4f+59]/(82875 f!).

Thus the scalar degree-four ambiguity is at least
`2^[binom(f,2)+4f+44]/(82875 f!)`. The alternative moment laws need not be
realized by different Burnside arithmetic sectors. This explicit qualification
is essential: the theorem is about the moment-positivity relaxation.

## 6. Exact experiments and what they cost

All columns below use the same odd-partition producer for a given backbone,
and moments only through four. The legacy column reproduces Session278's
bounds from that packet. Entries are sufficient low residue bits, not runtimes.

| n | B | legacy | unshifted per-M moment LP | shifted, M support | shifted, R support |
|---|---|---:|---:|---:|---:|
| 48 | {3} | 281 | 264 | 114 | 83 |
| 56 | {3} | 479 | 462 | 239 | 193 |
| 64 | {3} | 738 | 721 | 422 | 362 |
| 64 | {3,5} | 493 | 468 | 210 | 166 |
| 64 | {3,5,7} | 310 | 277 | 206 | 160 |
| 64 | {3,5,7,9,11} | 121 | 71 | 0 | 0 |
| 69 | {3,5,7,9,11} | see results | see results | see results | 0 |
| 70 | {3,5,7,9,11} | see results | see results | see results | 31 |

The retained packet adds `(delta,R)` to the prior edge-count sidecar. Hence the
gain is **not** an improvement from the identical information with no added
coordinate. The full partition producer can already compute the exact answer
more directly; it is intentionally used here to establish what information
is sufficient in the compressed consumer. Efficient generation of these
stratified moments by the neighboring marked-response DP remains OPEN.
There is no asymptotic speedup claim.

`ground_moments.py --consume packet.json` is a no-reference mode: it reads only
the promised cell moments, transforms them and certifies the interval. It never
calls the producer or computes an exact count. Reference checks in experiment
mode use a separate direct gcd exponent formula, not the ground expansion.
The saved packets for n64 use B={3} and B={3,5,7,9,11}, respectively. The latter
collapses to the exact1721-bit integer without a supplied answer or residue.

Reproduce from repository root:

```bash
P=05-knowledge/results/session279_interval_tightening_20261112/cross
python3 "$P/ground_moments.py" --experiment > /tmp/session279-ground-results.json
python3 "$P/ground_moments.py" --consume "$P/packet_n64_B3.json"
python3 "$P/ground_moments.py" --consume "$P/packet_n64_B11.json"
```

The no-reference consumer checks consistency conditions, not authenticity of a
malicious or erroneous producer. As with any certificate, its mathematical
input obligations must be established. The accompanying controls check every
displayed dual pointwise through support40, affine transforms against explicit
edge lists, core thresholds, and oracle containment for every experiment row.
