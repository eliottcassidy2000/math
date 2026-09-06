# Rational orientation transport cannot produce an exactly fair deadline extractor

**Status: PROVED + FINITE-EXACT; [independent audit accepted](continuing8_20260906_coin_rational_handoff_audit.md).** The algebraic obstruction below is self-contained. Its application uses the proved spine normal form. It excludes an entire stationary/periodic polynomial-state architecture, at every finite linear deadline slope. It does **not** construct a slope below two, prove a positive uniform gap above one, or exclude arbitrary finite descriptions with dilation or nonlinear decoding.

## Inheritance and corrected frontier

The current mathematical baseline is `fe3ac5816`. The relevant proved statements are:

- **THM-2966, `01-canon/theorems/THM-2966-spine-normal-form-for-critical-run-fair-extractors.md`:** actual row polynomials have integer Bernstein coefficients in their exact binomial boxes; the two orientation sums characterize exact fairness.
- **THM-3340, `01-canon/theorems/THM-3340-single-donor-cyclic-rotation-proves-all-pointwise-AMM-floors.md`:** the pointwise deadline frontier is closed. Optimizing different extractors for different rows does not supply one uniform extractor.
- **THM-3343 and THM-3344, `.../THM-3343-shifted-donor-rotation-bisects-exactly-the-dyadic-annuli.md` and `.../THM-3344-orientation-splitting-saves-exactly-one-dyadic-donor-bit.md`:** rotation supplies the donor capacity; one orientation split saves one donor bit. The latter sharpness theorem keeps every other row at its pointwise floor. It does not prohibit a later interior row from paying across an annulus.
- **THM-3337, `.../THM-3337-cross-shell-compression-attains-the-T4-floor.md`:** an actual cross-shell replacement has deadlines `(2,3,5,5,7,7,8)`. Its two nonzero shell residuals cancel exactly. In particular row 3 reads bit 5, beyond its old shell endpoint 4. This is the inherited cheap hostile to compulsory shellwise closure.
- **THM-3342, `.../THM-3342-sublinear-deadline-excess-is-impossible-for-fair-critical-run-extractors.md`, and MISTAKE-368:** every *fixed* extractor fails to have sublinear deadline excess. This does not prove that the infimum over extractors is greater than one.
- **THM-3577, `.../THM-3577-amm-r512-offset-transition-and-causal-horizon.md`, THM-3373, `.../THM-3373-r8-slack-kernel-conformal-locality-width-five.md`, and THM-4086, `.../THM-4086-rule-a-transition-clock-and-phase-cocycle.md`:** a finite policy failure is not alternative-prefix infeasibility; causal/slack labels and the dyadic phase survive any legitimate state compression.

The earlier HYP-9075 two-adjacent-shell negative computation concerns colorings that stay inside their original endpoints. It must not be extended to the THM-3337 spill. MISTAKE-361 is the other load-bearing correction: unbounded transport can erase a deadline obstruction, so an abstract capacity flow is not an extractor.

Targeted searches for rational orientation sums, finite-state coin transport, rational gluing, and the content argument found no prior statement of the result below. No external priority claim is made. In particular the new proof does not import Fatou, Kronecker, or a literature theorem about finite automata.

The live concept board is: (1) exact Bernstein integrality, (2) signed cross-shell residuals, (3) rational matrix realization, (4) reduction modulo a denominator prime, (5) dyadic dilation and the retained phase/row capacity. The new result connects (1), (3), and (4), while the controls below show why (2) alone is insufficient and why (5) remains available.

## 1. An elementary rational-gluing theorem

Let

`F_i(X) in Q(X) intersect Z[[X]]`,

where the rational function is regular at zero and its Taylor coefficients there are integers. Let `phi_i(X)=a_i X+b_i`, with nonzero integers `a_i` and integers `b_i`.

**Theorem A.** If the identity of rational functions

`sum_i F_i(phi_i(X)) = c in Q`

holds, every prime dividing the reduced denominator of `c` divides `product_i a_i`. In particular, if all `a_i` are `+1` or `-1`, then `c` is an integer. Integer signs or multiplicities on the summands can be absorbed in `F_i`.

**Proof.** Write `F_i=A_i/B_i`, with integral polynomials, normalized so that the joint gcd of all coefficients of `A_i` and `B_i` is one. Every `B_i` is primitive. Indeed, if a prime divided all its coefficients, the formal identity `B_i F_i=A_i` in `Z[[X]]` would force that prime to divide every coefficient of `A_i` too. This contradicts the normalization. No assertion about `B_i(0)=1` is required.

Write `c=u/v` in lowest terms and suppose a prime `ell` divides `v` but no `a_i`. Reduction of a primitive `B_i` modulo `ell` is a nonzero polynomial. Substitution by `a_i X+b_i` is an automorphism of `F_ell[X]`, so each reduced `B_i(phi_i)` is also nonzero. Clearing denominators in the displayed identity and reducing modulo `ell` gives

`0 = u product_i B_i(phi_i(X))` in `F_ell[X]`.

The right side is nonzero, since `u` is a unit and this polynomial ring is a domain. Contradiction. QED.

The conclusion controls denominator **prime support**, not prime valuations. The slope hypothesis matters. Set `F_d(T)=T/(1-dT)`, whose Taylor coefficients are integers. Then

`F_3(4X+1)+F_6(-X) = -1/2`.

The prime 2 is precisely permitted by the nonunit affine slope 4. Both sides are interpreted as rational functions on their common domain; no convergence of a Taylor expansion at the substituted argument is being assumed.

For the coin problem the relevant immediate corollary is

`F(X)+G(1-X) != 1/2` whenever `F,G in Q(X) intersect Z[[X]]`.

## 2. Every finite-linear-deadline extractor has nonrational orientation sums

Use `p=P(0)`, `q=1-p`, and let the critical value `m` be the initial constant-run length. Assume an actual deterministic exactly fair extractor has a pathwise deadline

`T(m) <= C m+D`, for fixed finite `C,D`.

By THM-2966, after refining stopped subtrees, the row at `0^m1` has a polynomial `W_m(p)` and the row at `1^m0` has a polynomial `V_m(u)` in its native variable `u=P(1)`. Their Bernstein depths satisfy `d_m <= gamma m+B`, for some fixed nonnegative `gamma,B`, and their head-count coefficients lie between zero and the corresponding binomial coefficient. Define the *head*, rather than doubled signed, orientation sums

`F(p)=sum_(m>=1) p^m(1-p)W_m(p)`,

`G(u)=sum_(m>=1) u^m(1-u)V_m(u)`.

**Theorem B.** Neither `F` nor `G` is rational. Equivalently, their ordinary Taylor coefficient sequences at zero each have infinite Hankel rank over `Q`.

**Proof of the analytic and arithmetic premises.** Each summand is an integral polynomial divisible by the corresponding `m`th power. Thus the formal coefficient at order `n` receives contributions only from `m<=n` and is an integer. These formal sums are genuine analytic germs at zero. More generally, for complex `z`, the exact Bernstein box gives

`|W_m(z)| <= (|z|+|1-z|)^d_m`.

Fix a real `p_0` in `(0,1)`. At `z=p_0`, the base

`|z| (|z|+|1-z|)^gamma`

is `p_0<1`. On a sufficiently small complex neighborhood it stays strictly below one, so the defining series for `F` converges there normally. The same argument works near zero and for `G`. Consequently the physical orientation sums are analytic along the entire real interval, and their germs at zero belong to `Z[[X]]`.

Exact fairness is `F(p)+G(1-p)=1/2` on `(0,1)`. If either orientation sum were rational, analytic continuation along the interval identifies it with that rational function everywhere on the interval. A pole there is impossible because the physical sum is analytic. The identity then makes the other orientation sum rational too, including its analytic germ at zero. The corollary of Theorem A contradicts the identity. QED.

For the Hankel formulation, the columns of the infinite Hankel matrix are the successive shifts of the coefficient sequence. Finite rank makes their span finite-dimensional and shift-invariant; the characteristic polynomial of the shift supplies a constant-coefficient linear recurrence. Its generating function is rational. Conversely a rational generating function has an eventual such recurrence, hence finite Hankel rank after adjoining the finite prefix. This proves the stated equivalence without an external classification theorem.

The linear bound is used for analytic identification, not for a numerical slope estimate. No superlinear-deadline extension is asserted. In particular this theorem gives no quantitative lower bound on `C*`.

It also shortens the proof of THM-3342: after that theorem's Pólya--Carlson step establishes rationality of the two integer germs under sublinear excess, Theorem A gives the contradiction immediately. Its subsequent Fatou normalization, root-of-unity pole classification, sixth-root restriction, and backward integer-valued-polynomial argument are unnecessary for that conclusion. This is a proof simplification, not a correction of its valid theorem or a new proof of the cited Pólya--Carlson input.

## 3. Exact finite-state consequence and its boundary

Suppose, for even one orientation, the actual row polynomials eventually have a fixed-dimensional polynomial realization

`s_(m+1)=A(p)s_m`, `W_m(p)=b(p)^T s_m`,

with a polynomial initial state and fixed polynomial matrix/vector over `Q`. A finite preperiod or a periodic list of matrices and outputs is allowed. Grouping one period of length `h`, its tail sum is a finite sum of terms with matrix inverse

`(I-p^h A_period(p))^-1`.

Its determinant has constant term one, so the resulting orientation sum is rational. Theorem B rules out exact fairness whenever these are genuine rows obeying a finite linear deadline and the integer Bernstein boxes. Finite prefix surgery cannot repair this architecture: it adds only a polynomial to `F` or `G`.

For an integer matrix realization of the actual head rows, the denominator-content obstruction can be read directly from the determinant formula. It provides a symbolic impossibility certificate without a state census. The assertion concerns **orientation-resolved polynomial row realization**. It does not mean that all finite-state algorithms or all small signed boundary summaries are impossible.

Three controls specify this boundary.

1. **A legal but unfair stationary rule.** At critical value `m`, read `m-1` additional bits and use their parity, with opposite output conventions on the two initial orientations. Then `T(m)=2m`, and

   `W_1=1`, `W_(m+1)=(2p-1)W_m+(1-p)`, `V_m=1-W_m`

   are actual Bernstein-box polynomials. Their orientation sums are

   `F(p)=p(1+p)/(1+2p)`, `G(u)=u^2/(1+2u)`.

   At `p=1/3` the head probability is `16/35`, not `1/2`. It is fair at the single bias `p=1/2`, illustrating why that test is insufficient.

2. **A rational fair algorithm outside the premise.** For the ordinary von Neumann pair rule, the two head-orientation sums are exactly

   `F(p)=p-p^2/2`, `G(u)=u^2/2`.

   They sum to `1/2` after reflection but have half-integral Taylor coefficients. Indeed the words `0011(00)^k01` all have critical value two and first produce an unequal pair at arbitrarily late times. No finite `T(2)` exists; an infinite constant-pair continuation can also remain undecided. This hostile prevents conflating the theorem with a ban on ordinary finite-state fair-bit extraction.

3. **A signed aggregate can look legal while its row is impossible.** Set `K_1=0` and, for `N>=2`,

   `K_N(p)=(1-2p)^(2N)-p^N-q^N`.

   This is a three-exponential stationary boundary state. It vanishes at both endpoints, satisfies the necessary parity `K_N=1-p^N-q^N (mod 2)`, and obeys the necessary tail bound `|K_N(p)|<=p^N+q^N`. For the last bound, `0< p,q<1` gives `(p-q)^2<=p^2+q^2<=max(p,q)`.

   Nevertheless its first handoff cannot be an actual row. The row-one signed contribution must be `pq(Delta_1+Delta'_1)` with each conditional signed share in `[-1,1]`. But

   `K_2/(pq)=-6+16pq`,

   which equals `-3` at `p=1/4`, below the allowed `-2`. Thus even parity, endpoint values, geometric decay, and the full pointwise tail capacity do not replace the orientation-resolved row box. This is an explicit stopping boundary, not an assertion that every aggregate-state realization forces rational orientation sums.

## 4. Connection contract and the surviving next object

The source is the actual deadline-labelled Bernstein spine. The map takes its head-count rows to the two ordinary orientation generating functions, then takes a proposed rational matrix realization to its integral denominator content modulo two. It preserves exact fairness and integer row realization; the affine prime theorem explains exactly when a chart operation preserves the obstruction. It deliberately loses the numerical deadline slope, so it cannot answer `C*<2`.

The needed sidecar is orientation and the full causal row box. The three-exponential hostile demonstrates the first failure if one retains only a signed aggregate. The THM-3337 control demonstrates the opposite failure: forcing each shell residual separately to vanish deletes a genuine exact solution.

The strongest survivor is a finite description with a **nonstationary scale operation**, such as dyadic dilation or a degree/phase-dependent transition. A Mahler-type state that retains the change of variable is not rationalized by the fixed-matrix argument. Nor is a nonlinear state-to-row decoder covered. The inherited THM-3343 dyadic construction is a positive control: the companion independently reconstructs its legal donor rows and its fair prefixes, and both orientation Hankel matrices of order 30 have full rank modulo 101. This is finite corroboration of an existing safe architecture, not a new all-scale construction inferred from ranks.

One precise next task is to specify a dyadic boundary state *together with* its orientation-resolved row decoder and prove a capacity-preserving recursion. The scalar phase-free compression already fails in THM-4086, while a fixed polynomial matrix decoder is now excluded here. The open object is therefore a scale-aware, capacity-preserving action; a finite coefficient fit or a zero aggregate residual alone is not a sufficient certificate.

No LRC, Laurent sign, or Smith result is asserted from the vocabulary of this connection. The common operation is concrete: retain integral content before applying a quotient or chart. Here reduction modulo a denominator prime produces an impossibility, whereas the recent Laurent norm computations used a nonzero modular resultant to certify a unit. Neither predicate transports without its own exact source map.

## 5. Reproduction and finite scope

Run the adjacent `continuing8_20260906_coin_rational_handoff.py` with Python and Python `-O`. It imports no repository producer or numerical solver. Its universe includes all 782 head-count boxes through depth four, invertible affine substitutions over `F_2,F_3,F_5`, exact rational stationary/von-Neumann controls, every nonconstant eight-bit word of the inherited THM-3337 replacement, and independent THM-3343 donor and complete prefix identities through horizon 64. Modular Hankel ranks are auxiliary finite controls, not the proof of Theorem B.

Both modes pass **5,070 active exact gates** and have byte-identical LF output. No finite deadline-profile optimization or large-scale policy bank was rerun. The source, ordinary output, optimized output, and this report are frozen outside the repository for independent audit; the parent task owns filing and status promotion.
