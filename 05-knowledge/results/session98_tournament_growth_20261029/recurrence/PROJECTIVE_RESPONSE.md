# Fixed-degree projective symmetry states for A000568

## Result first

The inherited finite-horizon gcd-response recurrence admits a **strictly coarser exact quotient when the requested output is one coefficient**. Replace the future response `h_j` by

\[
s_j=h_j-jh_1,
\]

and absorb `2^(r h_1)` into its weight, where `r` is the degree still to fill. Histories whose response vectors differ by a multiple of the size vector `(1,3,5,...)` now merge. This gives a closed insertion recurrence and a normalized formal differential hierarchy. It is an improvement in state count, not a new Burnside formula, polynomial-time algorithm, or individual-tournament growth classification.

The attached independent checker compares this quotient to an unnormalized response DP through 100 and to direct odd partitions through 30. At target 100, the sum of stage-state counts falls from **16,330 to 13,544**, and peak state count from **2,953 to 2,721**. Both versions in this comparison compute the same *single target*, with the same descending elimination and horizon pruning. This is not a comparison against the prior all-coefficients implementation's running time.

## 1. Inheritance and provenance

Read-only source snapshot: `origin/main` at `5862c94555eafea4e9a0638a17f7c661f8e45d37`. Relevant files:

* `01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md`: unit tail, prime tail, formal generating-function decomposition.
* `01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md`: divisor-profile quadratic form and add-a-cycle increment.
* `05-knowledge/results/weighted_core_transport_20261020/README.md`: cycle-index boundary and full odd-cycle differential hierarchy; failure of cycle index to predict individual-core children.
* `05-knowledge/results/tournament_future_boundary_20261020/README.md`, Propositions 1–2: finite-horizon response equivalence, multiplicity insertion, and pruning; already computes through 200.
* `05-knowledge/results/session86_relations_attachment_carriers_20261028/A000568_COMPOSITION.md`: independent Burnside/composition verification and observer-relative limitations.

The Davis/Burnside count is classical, also recorded at [OEIS A000568](https://oeis.org/A000568), consulted 2026-10-29. The PDE and the unnormalized finite-response quotient below are explicitly inherited. The additional claim here is the fixed-degree **projective** quotient and its normalized implementation; no historical-priority claim beyond this repository comparison is made.

## 2. Classical input and inherited functional equation

For a finite odd-cycle multiset `mu`, write

\[
e(\mu)=\sum_{a\in\mu}(a-1)/2+\sum_{a<b\text{ occurrences}}\gcd(a,b),\qquad
z_\mu=\prod_{j\text{ odd}}j^{m_j}m_j!.
\]

Occurrence indexing in the pair sum includes equal cycle lengths. Then
`a_n = sum_(|mu|=n) 2^e(mu)/z_mu`. An even cycle reverses an antipodal pair under a permutation power and cannot fix a tournament. All-odd cycle types admit independent orientations on their pair orbits. Burnside gives the formula.

For a symbolic parameter `q`, define

\[
\mathcal F(\mathbf y;q)=\sum_{\mathbf m\ge0}
q^{\sum_jm_j(j-1)/2+\sum_j j\binom{m_j}{2}+\sum_{i<j}\gcd(i,j)m_im_j}
\prod_{j\text{ odd}}\frac{y_j^{m_j}}{j^{m_j}m_j!}.
\]

All operations take place in the completion by `deg(y_j)=j`; each degree has finitely many terms. The inherited cycle-insertion PDE is

\[
\partial_{y_p}\mathcal F(\mathbf y;q)
=\frac{q^{(p-1)/2}}p
\mathcal F((q^{\gcd(p,j)}y_j)_j;q).                 \tag{1}
\]

Proof: differentiate a monomial and remove the distinguished p-cycle. The denominator loses a factor p and the exponent loses `(p-1)/2+sum_j gcd(p,j)m_j`. This proves (1), including repeated p-cycles.

With `y_j=x^j` and `q=2`, `A(x)=F((x^j);2)` is the ordinary A000568 generating function. For `h=(h_j)`, let

\[
A_h(x)=\mathcal F((2^{h_j}x^j)_j;2),\qquad a_n(h)=[x^n]A_h(x).
\]

The inherited hierarchy is

\[
A'_h(x)=\sum_{p\text{ odd}}2^{(p-1)/2+h_p}x^{p-1}A_{h+g_p}(x),
\quad g_p(j)=\gcd(p,j).
\]

It has initial coefficient `a_0(h)=1` and determines all coefficients recursively in degree. Only `h_j` for `j<=n` can affect `a_n(h)`.

## 3. The fixed-degree quotient theorem

Fix remaining mass `r` and a set `L` of available odd lengths containing 1; discard lengths greater than r. Two response vectors h and h' multiply **every exact-r continuation** by the same nonzero scalar iff

\[
h'_j-h_j=j(h'_1-h_1)\qquad(j\in L,\ j\le r).       \tag{2}
\]

Here a continuation is a multiplicity vector k with `sum_j j k_j=r`; its cross-interaction multiplier is `2^(sum_j h_j k_j)`. Intrinsic future weights and factorial denominators are the same for the two histories and cancel in ratios.

**Proof.** If (2) holds, the ratio is `2^(r(h'_1-h_1))`, independent of k. Conversely compare the continuation consisting of r unit cycles with the continuation consisting of one j-cycle and r-j unit cycles. Equality of their ratios forces `(h'_j-h_j)-j(h'_1-h_1)=0`. This tests every retained coordinate. For r=0 the condition is vacuous. QED.

Thus `s_j=h_j-jh_1`, with `s_1=0`, is the coarsest projective response invariant for this exact-degree, individual-continuation observer. This does **not** prove minimality for arbitrary symbolic algorithms or sums that allow cancellations.

Strictness is witnessed before processing lengths 3 and 1 by histories `(9,9,9,9,9)` and `(15,15,15)`. Both have mass 45. Their `(h_1,h_3)` are `(5,15)` and `(3,9)`, different for the inherited quotient but both have `s_3=0`. At a fixed remaining mass r their future ratios differ only by `2^(2r)`, absorbed into the accumulated weight. These histories are legitimate descending-elimination states: all their lengths exceed 3. A smaller witness is `(5,9,9)` versus `(23)`, both mass 23, with normalized `s_3=-2`.

## 4. Closed multiplicity recurrence

Take a fixed desired size N, process odd p in descending order, and let r be remaining mass. A partial history has weight `w=2^e/z` and response `h_j=sum_(a in history) gcd(a,j)`. Store

\[
W=w2^{rh_1},\qquad (r,(s_j)_{j\in L,j\le r}),\quad s_j=h_j-jh_1.
\]

Initialize `r=N`, `W=1`, and all s zero. For each `0<=k<=floor(r/p)`, insert k cycles of length p and set

\[
\begin{aligned}
r'&=r-kp,\\
s'_j&=s_j+k(\gcd(p,j)-j),\\
W'&=W\,\frac{2^{k(p-1)/2+p\binom{k}{2}+k s_p+r'k}}{p^k k!}. \tag{3}
\end{aligned}
\]

Discard p and all coordinates j>r'. Add weights of equal successor states. After processing 1, the weight at r=0 is exactly a_N. States with r>0 do not contribute to the target.

**Proof.** Ordinary insertion multiplies w by
`2^(k(p-1)/2+p*C(k,2)+k*h_p)/(p^k*k!)`. Also `h'_1=h_1+k`. Substitute `h_p=s_p+p*h_1` and compare `W'=w'2^(r'h'_1)` to `W=w2^(rh_1)`. The h_1 coefficient cancels because `r=kp+r'`, leaving (3). Response update follows by the same subtraction. Equality of normalized states makes every transition identical, so merging is exact by induction. QED.

Some normalized exponents may be negative; arithmetic is rational, and the implementation deliberately uses exact Fractions. Nonnegative final integer counts do not justify rounding intermediate weights.

At the final length 1, s is identically zero: all histories of the same committed mass merge, regardless of their number of cycles. Earlier, with remaining lengths {1,3}, only `s_3=-2*(number of historical cycles not divisible by 3)` survives. This is a concrete partial-structure classification simpler than `(total cycles, cycles divisible by 3)`. With {1,3,5}, the two nondivisibility counts give `s_3` and `s_5`; composite j retain the full gcd-weighted deficiency `s_j=-sum_a(j-gcd(a,j))`.

## 5. Normalized size recursion and formal GF

Define `b_n(s)=a_n(s)` for `s_1=0`. The gauge identity is

\[
a_n(h)=2^{nh_1}b_n((h_j-jh_1)_j).
\]

Normalize the inherited hierarchy to obtain the exact degree recursion

\[
\boxed{\quad n b_n(s)=\sum_{\substack{p\le n\\p\text{ odd}}}
 2^{(p-1)/2+s_p+n-p}
 b_{n-p}((s_j+\gcd(p,j)-j)_{j\le n-p}).\quad}       \tag{4}
\]

The initial value is `b_0=1`, and `a_n=b_n(0)`. No unknown core bank occurs. A formal-series version, setting `B_s(x)=sum_n b_n(s)x^n`, is

\[
B'_s(x)=\sum_{p\text{ odd}}2^{(p-1)/2+s_p}x^{p-1}
 B_{s+g_p-\mathrm{id}}(2x).                         \tag{5}
\]

Only finitely many terms affect a given coefficient, so (4) and (5) are well-defined formal identities, despite infinite state dimension. This is a normalized hierarchy, **not a scalar finite-order differential equation**.

The analytic-convergence boundary is decisive: `a_n >= 2^(n(n-1)/2)/n!`, the identity-permutation contribution. Consequently `A(x)` has radius of convergence zero. The earlier phrase “fixed analytic background” for the unit-tail series must therefore be read as a formal background, not as a convergent power series near zero.

## 6. Verification and limitations

Run `python3 verify.py > receipt.json`. The checker uses explicit exceptions, so checks remain active under `python3 -O`.

* Direct independent partition summation equals normalized and raw response insertion for every n<=30.
* Equation (4), memoized independently of multiplicity insertion, agrees through n=20.
* The initial 11 terms agree with the A000568 prefix.
* Normalized and raw fixed-target recurrences agree at n=40,60,100.
* A strict reachable quotient collision is explicitly checked.

| N | Raw peak | Projective peak | Sum raw stage states | Sum projective stage states |
|---:|---:|---:|---:|---:|
|40|130|105|638|510|
|60|420|361|2,238|1,780|
|100|2,953|2,721|16,330|13,544|

These state counts include all stored remaining-mass states, including terminal non-target masses; both compared algorithms use the identical convention. They count dictionary keys after each processed length, not operations, memory bytes, or partition histories. The prior packet's n=200 computation remains the larger existing range.

The quotient sacrifices responses at other target degrees; to obtain several a_n, this implementation runs separately for each n. It does not recover individual child distributions, automorphism groups, or edge-flip metagraph neighbors. It offers a precise answer to “which partial symmetry data can be forgotten as size changes?”: under a fixed final mass, a uniform *per-vertex* response can be folded into weight; only deviations from that linear size response remain observable.
