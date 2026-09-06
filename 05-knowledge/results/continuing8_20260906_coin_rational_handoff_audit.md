# Independent referee: rational orientation handoff

**Status: PROVED / independently audited relative to the inherited spine normal form.** The producer's Theorems A and B, their stationary/periodic polynomial-state consequence, and their stated scope boundaries pass independent mathematical review. No mathematical repair was requested. The separately labelled polynomial-substitution corollary below is a parent-proposed extension, not a change to the frozen producer.

The reviewed primary is `continuing8_20260906_coin_rational_handoff.md`, frozen SHA256 `600e4dae1f4815c4a2ba2eb69c66140c8d6c8920e8e4a4a05db407fd86fe6eff`. Its source SHA256 is `359eba70c2e89226e44486441919d6a95a911188054de6fd19f9e09de7ead804`; its output SHA256 is `94e04c1f7b92b12f70952e5542370fec4ccb427f3e7890ddd550ea3d2e84c86e`. The independent companion pins those source/output bytes but does not import or execute that producer.

## 1. Algebraic content and chart scope

For each rational integral Taylor germ write `F=A/B` with integral polynomials and joint coefficient gcd one. If a prime divided the content of B, coefficientwise reduction of the formal identity `BF=A` would force it to divide every coefficient of A. Thus B is primitive. This is sufficient; no assertion that B has constant coefficient one, no integer-valued-polynomial theorem, and no Fatou normalization are required.

Let `c=u/v` be reduced and let a prime ell divide v. Clearing all denominators in an identity `sum F_i(a_i X+b_i)=c` and reducing modulo ell gives

`0 = u product B_i(a_i X+b_i)`.

If every a_i is nonzero modulo ell, each substituted denominator is a nonzero polynomial, since affine substitution is an automorphism of the polynomial ring over the field. The right side cannot vanish in a domain. This proves the claimed denominator-prime support restriction; it does not give denominator valuations. Signs and integer multiplicities can indeed be absorbed into the integral germs.

The exact hostile `F_3(4X+1)+F_6(-X)=-1/2`, where `F_d(T)=T/(1-dT)`, verifies that excluding nonunit chart slopes would be false. It is a rational-function identity, so no Taylor convergence at the substituted input is presumed.

## 2. Analytic passage from a genuine extractor

The application uses **THM-2966, `01-canon/theorems/THM-2966-spine-normal-form-for-critical-run-fair-extractors.md`**: actual head rows have integral Bernstein counts in their binomial boxes, separately for the two initial orientations. A pathwise deadline `T(m)<=Cm+D` gives Bernstein depth `d_m<=gamma m+B`, with fixed nonnegative gamma and B.

Each summand `X^m(1-X)W_m(X)` is an integral polynomial divisible by `X^m`. Hence every formal coefficient is a finite sum of integers. To identify that germ with the physical orientation probability, put `S(z)=|z|+|1-z|>=1`. The exact box yields

`|z^m(1-z)W_m(z)| <= |1-z| S(z)^B [|z| S(z)^gamma]^m`.

At every real point in `[0,1)`, the bracket is strictly less than one. Continuity gives a complex neighbourhood with a uniform geometric majorant, so the series is normally convergent there. This includes the origin and justifies both the integer analytic germ and analytic continuation along the physical interval. The argument for the other orientation uses its own native variable.

If one germ equals a rational function R, overlapping analytic neighbourhoods continue that equality along `(0,1)`. A purported first pole there is removable because the physical sum is analytic. Fairness then identifies the other physical orientation with `1/2-R(1-X)`. Near its origin, any endpoint pole of the reflected rational function is also removable, since the other orientation is analytic there. Thus one rational orientation would force two rational integral germs with reflected sum one half, contradicting Section 1. Both orientations are nonrational.

This proof uses finite linear depth growth for normal convergence. It proves no quantitative slope gap, no superlinear-deadline result, and no impossibility for ordinary von Neumann extraction.

For completeness, the infinite Hankel columns are the successive left shifts of the coefficient sequence. Their finite-dimensional span, if finite, is preserved by coordinatewise left shift; Cayley--Hamilton gives a constant-coefficient recurrence and a rational generating function. Conversely an eventual recurrence leaves only a finite prefix to adjoin to a finite column span. Thus the infinite-rank reformulation is exact. It does not say that every leading finite Hankel determinant is nonzero.

## 3. The architecture consequence retains the right source

An eventual fixed polynomial matrix recurrence for the actual row polynomials gives a rational orientation sum by a matrix resolvent. For period h, grouping phases yields denominators `det(I-X^h M(X))`, whose constant coefficient is one. A finite preperiod or finite row surgery adds only a polynomial. One such orientation already contradicts Section 2 under actual fairness and finite linear deadlines.

This is specifically an orientation-resolved polynomial row realization. An aggregate signed state can fail to reconstruct a row in its Bernstein box. A nonstationary substitution/dilation, a degree-dependent action, or nonlinear row decoder is not reduced to this resolvent argument. The producer keeps those exclusions explicit.

The inherited cross-shell control **THM-3337, `01-canon/theorems/THM-3337-cross-shell-compression-attains-the-T4-floor.md`**, is essential: its actual row 3 reads bit 5 beyond the old endpoint 4, and two nonzero shell residuals cancel. Thus forcing each shell separately to close would discard a genuine source. The positive dyadic control is **THM-3343, `01-canon/theorems/THM-3343-shifted-donor-rotation-bisects-exactly-the-dyadic-annuli.md`**; finite Hankel measurements here corroborate it but do not prove its all-scale fairness.

## 4. Independent exact companion

The verifier uses SymPy only for exact rational polynomial/matrix identities; its finite-field, literal word, donor, and determinant paths are separately implemented. All checks use a raising gate rather than removable assertions.

- It derives the affine coefficient-change determinant symbolically through degree six and checks the nonunit-slope exception.
- It constructs every one of the 782 Bernstein boxes through depth four by selecting actual binary suffix words in Hamming layers, then compares direct leaf probabilities to the reconstructed power polynomial at two rational biases.
- It sums the legal stationary example by a two-dimensional matrix resolvent and recovers `F=X(1+X)/(1+2X)`, `G=X^2/(1+2X)` and head probability `16/35` at bias one third.
- It checks the literal first-pair von Neumann decomposition and twenty words of critical value two with stopping time `6+2k`. Its fair rational germs have half-integral coefficients; it has no finite deadline at that critical value.
- It independently divides the aggregate handoff at the first row, obtaining `-6+16pq`, equal to `-3` at `p=1/4`, outside the allowed sum of two signed shares in `[-1,1]`. Endpoint, parity and tail controls alone therefore fail.
- It reconstructs each dyadic donor by subtracting literal prescribed heads from half of an entire annulus Hamming layer, rather than using the producer's signed-defect recursion. All six complete dyadic prefixes through horizon 64 satisfy the exact polynomial fairness identity and both donor boxes.
- It uses a fraction-free modular determinant recurrence, rather than the producer's rank elimination. At sizes `2,4,8,12,16,24,30`, the two shifted Hankel determinants modulo 101 are respectively `(100,0),(6,3),(21,87),(65,39),(80,77),(5,78),(75,20)`. The exceptional second orientation at size two has rank one; all other displayed determinants are nonzero.
- It enumerates all 254 nonconstant eight-bit words in the inherited spill. The head counts by weight are `(0,4,14,28,35,28,14,4,0)`. The two shell signed profiles are `(0,0,-2,-4,0,4,2,0,0)` and its negative.

These are finite controls accompanying the analytic proof, not a claim that large observed Hankel rank by itself establishes nonrationality.

## 5. Separately audited polynomial-substitution extension

**Corollary (parent-proposed, proved here).** Retain `F_i in Q(X) intersect Z[[X]]`, regular at zero, and replace the affine charts by arbitrary nonconstant integral polynomials phi_i. If `sum F_i(phi_i(X))=c in Q`, every prime ell dividing the reduced denominator of c makes at least one phi_i constant modulo ell. Equivalently, ell divides every nonconstant coefficient of at least one chart.

**Proof.** Primitive B_i reduce to nonzero polynomials. Composition by a nonconstant polynomial over a field is injective: for a nonzero polynomial of positive degree its composed degree is the product of the two positive degrees, and nonzero constants stay nonzero. Consequently, if every reduced phi_i is nonconstant, all reduced `B_i(phi_i)` are nonzero. The cleared-product contradiction of Section 1 applies unchanged. A reduction that lowers the degree but leaves it positive is harmless. QED.

In particular, charts whose nonconstant coefficient content is one all preserve the integer-constant obstruction, including `X^2` and `2X-X^2`. The companion checks all 15,546 pairs of a nonzero outer polynomial and a nonconstant inner polynomial of degree at most two over `F_2,F_3,F_5`, by direct coefficient convolution, and verifies the degree-product identity. It also checks the two displayed charts stay nonconstant in these characteristics.

This corollary still assumes rational F_i. It does not prohibit a Mahler functional equation involving genuinely nonrational series, and it does not convert a scale-dependent fair-coin architecture into a stationary one. The retained object in a useful Mahler/coin connection is the actual chart action together with the orientation row decoder and capacity, not the common word “dyadic.”

## 6. Frozen reproduction

Run `python3 -B continuing8_20260906_coin_rational_handoff_audit.py` and the same command with `-O`. The script locates the frozen producer source beside itself when filed and its output in `05-knowledge/results`; outside-workspace fallback paths are explicit. No producer import, state census, deadline optimization, or network dependency is used.

Both modes pass **18,256 exact always-active gates**, with byte-identical actual LF output. Audit source SHA256: `2fd45385c767a4c23d170ff9c36be1d531455ec7d435a7718cd2bfab44bb7c73`. Audit output SHA256: `474c15855caf0aea2e00ae0da9d0b136e89d91341a208e5f2d30f44bebf012fb`.

All primary source/output bytes remain unchanged. The parent task owns filing, current-status promotion and Git integration.
