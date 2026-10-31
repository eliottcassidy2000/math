# Session 127: GSMA inheritance, exact norm gate, and sparse-unit kernel

## Executive result

The user's 2026-10-31 clarification names **GSMA = Generative Source Morphism Architechture** (their spelling). Searches found no pre-existing formal definition under this name in the pinned repository or connected saved conversations. Accordingly the architecture below is an **explicit proposed instantiation**, not a recovered theorem under an established acronym.

There is an exact, uniform mathematical merge. For `k = 3, 4, 5`, set

\[
f_k(X)=X^k-\sum_{j=0}^{k-1}X^j,
\quad F_k(P,Q)=P^k-\sum_{j=0}^{k-1}P^jQ^{k-j},
\quad U_k(z)=\frac{z^{k-1}}{1-z-\cdots-z^k}.
\]

For coprime integers `p,q` with `p != 0`,

\[
\boxed{U_k(q/p)=\frac{p q^{k-1}}{F_k(p,q)},\qquad
\gcd(F_k(p,q),pq)=1.}
\]

Thus the reduced denominator of the rational value is **exactly** `|F_k(p,q)|`:

* `U_k(q/p)` is an integer iff `F_k(p,q) = +/-1`.
* For a **fixed finite set S of rational primes**, `U_k(q/p)` belongs to `Z[S^{-1}]` iff every prime dividing `F_k(p,q)` lies in S.

This is a genuine Thue–Mahler norm condition. It is not a statement about Mahler's three-halves map or Thue–Morse words. The rational-function identity holds for every `k >= 2`; irreducibility of `f_k` is verified here for the requested three orders.

The cyclic inverse-unit source has the exact finite kernel census:

| Order | All exponents `0 <= n <= 500` with `alpha_k^{-n}` in `Z + Z alpha_k` |
|---|---|
| 3 | `n=0`: `1`; `n=3`: `2-alpha_3`; `n=16`: `-103+56 alpha_3` |
| 4 | `n=0`: `1`; `n=4`: `2-alpha_4` |
| 5 | `n=0`: `1`; `n=5`: `2-alpha_5` |

**FINITE-EXACT**, not an all-exponent theorem at orders 4 and 5. The companion Tribonacci lane supplies a separately cited global cubic classification and certified continued-fraction analysis. The identity at exponent 16 is already in the literature; no novelty claim is made.

## 1. Exact inheritance and correction hierarchy

All repository evidence is pinned to commit
`522535a7b6f519fd0c1e4304a2c5b1c909fce4c9` of `eliottcassidy2000/math`.
`SOURCE_PINS.json` records full paths, Git blob IDs, SHA-256 hashes and byte sizes.

**Closest proved mechanism:** `THM-227`,
`01-canon/theorems/THM-227-knacci-mersenne.md`, proves the k-nacci companion-matrix trace law and explicitly records

\[
(X-1)f_k(X)=X^{k+1}-2X^k+1.
\]

Its exact polynomial identity is relevant here; the trace law alone does not imply a norm-collapse theorem.

**Historical observation recovered:**
`04-computation/rational_irrational_symbiosis.py` and
`05-knowledge/results/rational_irrational_symbiosis.out`, together with
`04-computation/free_exploration_s90g.py` and
`05-knowledge/results/free_exploration_s90g.out`, record the Tribonacci continued-fraction prefix `[1;1,5,4,2,305,...]`, its convergent `103/56`, and the next convergent `31461/17105`. These are floating-point historical computations, not exact certificates for their full printed tails. The new arithmetic must replace, not silently inherit, floating-tail claims.

**Canonical hostile / corrected near miss:** `01-canon/MISTAKES.md`, MISTAKE-028, retracts the inference from k-nacci/Mersenne trace coincidences to forbidden tournament Hamiltonian counts. In particular 31, 63 and 127 are attained, and the same Tribonacci trace sequence contains both forbidden and attained counts. Neither the value 305 nor a unit identity establishes a tournament-spectrum result. The old scripts contain speculative rhetoric outside their useful CF observation; that prose is not promoted.

**Least-used relevant sidecar:** the explicit connection contract in
`00-navigation/RESEARCH-PROTOCOL.md`, section 5: source, target, map, preserved predicate, information destroyed, restoration sidecar, decisive test. This is the closest authenticated repository architecture and is applied literally below.

**Negative recovery result:** targeted filename enumeration; targeted Git shows; connected GitHub content searches for GSMA, Generative Source, Morphism Architecture, Thue-Mahler, norm collapse, and source morphism; saved-conversation searches; and shared-note searches did not recover a pre-existing GSMA formalism or a named all-order Thue–Mahler norm-collapse theorem. This is a scoped search result, not proof that no such private source exists.

## 2. Fields, norm forms, and the denominator theorem

Let `alpha_k` be the unique positive real zero of `f_k`. Descartes' rule gives at most one positive root, while `f_k(1)=1-k<0` and `f_k(2)=1>0` put it in `(1,2)`. Work in the order `A_k=Z[alpha_k]`, with ordered basis `(1,alpha_k,...,alpha_k^{k-1})`; no maximal-order identification is needed.

The supplied checker certifies irreducibility modulo 3, 2, and 5 for `k=3,4,5`, respectively, using the Frobenius/Rabin criterion. Consequently these define number fields of degrees 3,4,5. Multiplication by `p-q alpha_k` has determinant

\[
N_{\mathbb Q(\alpha_k)/\mathbb Q}(p-q\alpha_k)
 =\prod_i(p-q\alpha_{k,i})=F_k(p,q).
\]

The forms are

\[
F_3=p^3-p^2q-pq^2-q^3,
\quad F_4=p^4-p^3q-p^2q^2-pq^3-q^4,
\]
\[
F_5=p^5-p^4q-p^3q^2-p^2q^3-pq^4-q^5.
\]

Clearing the rational-function denominator gives the displayed evaluation identity. If a prime divides `p`, then `F_k(p,q) = -q^k` modulo that prime; if it divides `q`, then `F_k(p,q)=p^k` modulo it. Coprimality of `p,q` therefore proves `gcd(F_k(p,q),pq)=1`, and hence the no-cancellation assertion. If `q=0`, primitivity forces `p=+/-1`, and both sides correctly give the value zero with denominator one.

No rational input is a pole: a rational zero of the monic `f_k` with constant term `-1` would be `+/-1`, and direct substitution excludes both for `k>=2`.

The exact S-integrality criterion is the special `a=+/-1` Thue–Mahler equation

\[
F_k(p,q)=\pm\prod_{\ell\in S}\ell^{e_\ell},\qquad
\gcd(p,q)=1,\quad e_\ell\ge0.
\]

This matches the definition in Gherga–Siksek, *Efficient resolution of Thue-Mahler equations*, arXiv:2207.14492v2 (2024), published Algebra & Number Theory 19 (2025), 667–714:
https://arxiv.org/abs/2207.14492v2 . The introduction, equation (2) and its following paragraph, also states the classical Thue–Mahler finiteness theorem. **CITED consequence:** for each of the three orders and each fixed finite S, only finitely many rational inputs have S-integral values. This gives finiteness, not an explicit list or an effective numerical bound in this packet; no solution list is imported from that paper.

Moreover `F_k(p,q)=+/-1` iff `p-q alpha_k` is a unit in `A_k`: determinant `+/-1` gives an integral inverse multiplication matrix. The corresponding localized statement is exact as well: the element is a unit of `A_k[S^{-1}]` iff its norm is an S-unit. The adjugate matrix proves sufficiency, multiplicativity of determinant proves necessity.

## 3. Sparse inverse-unit source and high-coordinate kernel

The constant term of `f_k` gives

\[
\alpha_k^{-1}=\alpha_k^{k-1}-\sum_{j=0}^{k-2}\alpha_k^j.
\]

Write `alpha_k^{-n}=sum_j c_j(n) alpha_k^j`. Multiplication by the inverse has the **integer, invertible** update

\[
(c_0,c_1,\ldots,c_{k-1})\longmapsto
(c_1-c_0,c_2-c_0,\ldots,c_{k-1}-c_0,c_0).
\]

To prove it, multiply a proposed next vector by `alpha_k` and reduce `alpha_k^k=sum_{j=0}^{k-1}alpha_k^j`; the constant coefficient is the old `c_0` and every subsequent coefficient is the corresponding old `c_j`.

Define the additive projection

\[
\pi_{\rm high}:A_k\longrightarrow\mathbb Z^{k-2},\qquad
\sum_j c_j\alpha_k^j\mapsto(c_2,\ldots,c_{k-1}).
\]

Then the rank-two collapse is **membership in the kernel** `Z+Z alpha_k`, not vanishing of the norm. The kernel is an additive lattice, generally **not a subring**. The unit powers form a cyclic multiplicative subgroup, not the entire unit group. The orbit is tested against the additive kernel; these structures must not be conflated.

The inherited telescoping identity gives the uniform baseline

\[
\boxed{\alpha_k^{-k}=2-\alpha_k.}
\]

It also gives the homogeneous residual identity

\[
\boxed{(p-q)F_k(p,q)=p^{k+1}-2p^kq+q^{k+1}.}
\]

For primitive `p,q`, a useful local sidecar is
`gcd(F_k(p,q),p-q) | (k-1)`, since reduction modulo `p-q` gives
`F_k(p,q) = (1-k)q^k` and `gcd(q,p-q)=1`.

The exceptional cubic point obeys

\[
\alpha_3^{-16}=-103+56\alpha_3,
\quad F_3(103,56)=-1,
\quad U_3(56/103)=-323008.
\]

The exact coefficient map first gives `(p,q)=(-103,-56)`, with norm `+1`. Normalizing to positive denominator gives `(103,56)` by multiplying the algebraic element by `-1`; the odd-degree norm changes to `-1`, while the rational input `q/p` is unchanged.

The scalar `-323008` is `-103*56^2`, not a Tribonacci term or a dimension-volume value.

## 4. Domain, loss, and hostile controls

### Analytic domain is load-bearing

`U_k` is the generating function of the fundamental k-step recurrence with coefficients beginning at degree `k-1`. Its Taylor series has radius `rho_k=1/alpha_k`. To see this directly, `D_k(t)=1-sum_{j=1}^k t^j` is decreasing for positive `t` and has its unique positive zero `rho_k`; no complex zero can lie at smaller modulus because `|sum z^j| <= sum |z|^j < 1`.

At the baseline `p/q=2`, the input `1/2` is inside the series disk and `U_k(1/2)=2`. At the cubic exceptional point, the positive identity `alpha_3^{-16}=56 alpha_3-103` gives `56/103 > 1/alpha_3`. Thus the negative value is a **rational continuation value outside the positive-series convergence disk**, not the sum of a convergent positive-coefficient series.

### Norm is not the sparse kernel

Every `alpha_k^{-n}` is a unit, with norm `(-1)^{(k+1)n}`. Yet `alpha_k^{-1}` has a nonzero coefficient of `alpha_k^{k-1}` for every `k>=3`, so it fails the rank-two kernel test. Keeping the norm but deleting the high coordinates destroys precisely the property we need to detect.

### A cyclic-unit census is not a full Thue census

At order four, `F_4(1,-1)=1`, so `1+alpha_4` is a unit and `U_4(-1)=-1`. But its selected real embedding is greater than one, so it is not any negative power `alpha_4^{-n}` with `n>=0`. Even a complete classification of negative powers on the rank-two plane would not classify all rational integer values. At orders four and five no all-height Thue or Thue–Mahler solution classification is asserted here.

### Scalar value does not preserve source identity

The norm observer is multiplicative, not additive; it is not a ring homomorphism. The high-coordinate projection is additive, not generally multiplicative. GSMA must use a typed diagram, not call all arrows homomorphisms of one unnamed structure.

## 5. Proposed GSMA instance and connection contract

| Stage | Type and map | Preserved | Lost / required sidecar |
|---|---|---|---|
| Generative source | `n -> alpha_k^{-n}` in a cyclic unit subgroup | multiplicative addition of exponents | k, selected generator and exponent |
| Coefficient realization | ordered integral power basis | exact algebraic element, all coefficient relations | basis and defining polynomial must remain fixed |
| Sparse gate | test `pi_high=0` | exact membership in `Z+Z alpha_k` | cannot replace this test by unit norm |
| Binary form | `a+b alpha -> (p,q)=(a,-b)` | exact element and signed norm | coefficient signs, primitive normalization |
| Rational GF observer | `(p,q) -> U_k(q/p)` for `p != 0` | integrality / S-integrality iff unit / S-unit norm | generally loses pair, exponent, embedding and convergence domain |
| CF observer | selected positive root plus a convergent and its predecessor | exact approximation / next quotient when certified | norm alone loses root choice and CF-adjacency data |

**Cheapest decisive tests:** baseline `(p,q)=(2,1)` at every order; cubic `(103,56)` with opposite convergence-domain behavior; inverse-power `n=1` as unit-but-not-sparse hostile; quartic `(1,-1)` as a full-unit versus negative-power-source hostile. All relevant integer identities are checked in both normal Python and `python -O`.

## 6. Reproduction and exact scope

Run from this directory:

```
python3 -B check_unit_kernel.py > unit_kernel.normal.json
python3 -B -O check_unit_kernel.py > unit_kernel.optimized.json
cmp unit_kernel.normal.json unit_kernel.optimized.json
```

The producer update is replayed independently using generic integer polynomial convolution, monic reduction and binary powering of the explicit inverse. A second cancellation check multiplies every produced vector by `alpha_k^n` and obtains 1. Exact multiplication-matrix determinants verify the norm forms at every primitive pair in `[-8,8]^2`, alongside the rational-function and local-gcd identities. Frobenius certificates prove the three irreducibilities. Saved checks remain active under optimization. Vector-sequence hashes and selected full checkpoint vectors make the finite universe auditable.

**Most promising bounded extension:** use the source-retaining sparse kernel to organize candidates, then solve the independent full Thue–Mahler problem for an explicitly fixed S. The architecture prevents a bounded orbit census from masquerading as a global Diophantine closure. The exact denominator theorem is already uniform in order; global arithmetic completeness, not a larger numerical coincidence, is the remaining issue.
