# Primitive prime-power Hecke traces: a marked extension of PR93

2026-10-18. Owned research packet; no change to the imported proofs.

**Result.** The primitive cyclic correspondence at level $p^r$ retains a transitive inert Cartan action. Its exact weight-two Newton slope gives a trace divisible by $p^{2r}$. At level $p^2$, adding the single imprimitive scalar branch destroys this conclusion: its trace contribution is $p^2$. This extends the *Hecke trace* mechanism; it does not establish the still-open prime-power truncated-series congruences.

**Status.** The finite-ring maps, formal-series norm identities, coefficient budget and support counterexample below are **PROVED**. The CM application is **PROVED given the stated classical full-level CM reciprocity input** used by the imported prime-level argument. Its modular source/frame hypotheses are retained explicitly. The six reconstructed $s_7$ trace polynomials and their CM evaluations are **FINITE-EXACT**. The higher Dwork/truncation comparison and the level-11 coefficient obstruction remain **OPEN**. No general literature priority is claimed for a prime-power Hecke/CM consequence.

## 1. Inheritance and the supplied paper

The preserved PR93 head was `5b7f276c3274046af1688a382b4cb183e3c23be5`. The parent session imported all nine files. The inputs here are:

- [PR power-structure and inert-CM proof](beukers_remark114_power_structure_and_inert_cm_20261018.md), especially its Hecke polynomial, norm-unit statement and inert-prime argument;
- [PR weight-one $s_7$ proof](cooper_s7_weight_one_root_frobenius_20261018.md), especially the retained Fricke sheet $s^2=(1+t)(1-27t)$ and its stated open higher-power/level-11 problems;
- [PR reflection proof and attribution repair](hyp9310_is_klmsy_corollary_reflection_proof_20261018.md): the harmonic CM law is attributed to KLMSY, not treated as newly established here.

The supplied 42-page attachment is Huber–Schultz–Ye, *Ramanujan–Sato series for $1/\pi$*, Acta Arithmetica 207 (2023), 121–160, [DOI](https://doi.org/10.4064/aa220621-19-12), local source [PDF](<C:/Users/Eliott/Downloads/Ramanujan_Sato series for 1__.pdf>). PDF pages 5–7, 13–22 were extracted; complete pages 5, 13 and 22 were rendered and visually read. The cover is PDF page 1, so printed page 12 is PDF page 13.

Three mechanisms matter here, as **CITED** supplied inputs:

1. Theorems 2.1–2.3 build a differential equation from a *marked modular form and Hauptmodul*, retaining the ramification factor $w$. A scalar modular value alone does not specify the differential frame.
2. Lemma 3.1 and Proposition 3.2 use **primitive** matrices: $ad=n$, $0\le b<d$, $\gcd(a,b,d)=1$. The degree is $\psi(n)=n\prod_{\ell\mid n}(1+1/\ell)$. This is the usable bridge to prime powers.
3. Theorems 3.5 and 3.8 retain a branch tangent and second derivative when recovering a Ramanujan–Sato coefficient. They exemplify why equality of a CM value does not authorize forgetting its marked local branch.

[Beukers, *Supercongruences using modular forms*, arXiv:2403.03301v3](https://arxiv.org/abs/2403.03301v3), and the PR's cited CM reciprocity argument supply the modular setting. Their proofs are accepted inputs, not re-audited in this packet.

Concept board: primitive cyclic support; Cartan projective frame; integral norm; valuation versus integer precision; missing Dwork comparison. The cheapest hostile is the scalar $pI$ branch, not another prime-level census.

## 2. A lossless projective-circle codec at every odd prime-power level

Let $p$ be odd, $R=\mathbb Z/p^r\mathbb Z$, and let $d$ be a nonsquare modulo $p$. Set $O=R[\omega]$, $\omega^2=d$, with conjugation $\bar\omega=-\omega$.

A vector $(a,b)$ is primitive iff at least one coordinate is a unit. Its norm $a^2-db^2$ is then a unit, since $d$ is nonsquare. Thus

$$\mathbf P^1(R)=O^*/R^*,\qquad |\mathbf P^1(R)|=p^{r-1}(p+1)=:n_r.$$

Multiplication by $O^*$ acts freely and transitively on these classes: a primitive $z$ is itself a unit, and the unique class taking $[z]$ to $[w]$ is $[w/z]$. Explicit projective representatives are

$$ (1,b),\ b\in R;\qquad (a,1),\ a\in pR.$$

The actual circle map is

$$[z]\longmapsto z/\bar z\in\{w\in O^*:w\bar w=1\}.\tag{1}$$

Its kernel is $R^*$: $z=\bar z$ forces the $\omega$ coefficient to vanish because 2 is invertible. An explicit inverse has two charts:

$$z=1+w\quad\text{if }1+w\text{ is a unit};\qquad
z=\omega(1-w)\quad\text{otherwise}.\tag{2}$$

In the second case $w\equiv-1\pmod p$, so $1-w$ is a unit. Both choices satisfy $z/\bar z=w$. This proves bijectivity without a cardinality guess. The second chart is essential at $w=-1$. Reduction from depth $r$ to $r-1$ has exactly $p$ preimages per projective point when $r>1$.

**Preservation statement.** The map retains the full primitive cyclic address and the Cartan multiplication action. Reduction modulo $p$ forgets $r-1$ digits. For $r>1$ this is a quadratic *ring*, not a finite field. No planar unit-distance realization or colouring conclusion follows from calling (1) a circle.

## 3. Primitive norms are integral units modulo $p$

For $F(q)\in1+q\mathbb Z[[q]]$, define the cyclic norm

$$\mathcal N_m F(q)=\prod_{b=0}^{m-1}F(\zeta_m^b q^{1/m}).$$

Translation invariance removes fractional powers; Galois invariance and algebraic integrality give $\mathcal N_mF\in1+q\mathbb Z[[q]]$. Reduction modulo a prime above $p$ gives

$$\mathcal N_{p^j}F(q)\equiv F(q)\pmod p.\tag{3}$$

Indeed all $p$-power roots of unity reduce to 1, and Frobenius turns the product into $F(q)$; the coefficients on both sides are rational integers, so the congruence is modulo $p$, not merely modulo a chosen uniformizer.

For determinant $p^r$, write the primitive upper triangular representatives as

$$\gamma_{j,b}=\begin{pmatrix}p^{r-j}&b\\0&p^j\end{pmatrix},\quad
0\le j\le r,\quad0\le b<p^j,\quad\gcd(p^{r-j},b,p^j)=1.$$

There are $c_0=1$, $c_r=p^r$, and $c_j=p^j-p^{j-1}$ for $0<j<r$. The primitive product is exactly

$$\prod_\gamma F(q^{p^{r-2j}}\zeta_{p^j}^b)
=F(q^{p^r})\mathcal N_{p^r}F(q)
\prod_{j=1}^{r-1}
\frac{\mathcal N_{p^j}F(q^{p^{r-j}})}
 {\mathcal N_{p^{j-1}}F(q^{p^{r-j-1}})}.\tag{4}$$

The denominators are integral unit series. Modulo $p$, (3) telescopes (4) to

$$F(q^{p^r})F(q^{p^{r-1}})\equiv F(q)^{n_r}.$$

Consequently

$$U_r(q):=\frac{\prod_\gamma F(q^{p^{r-2j}}\zeta_{p^j}^b)}{F(q)^{n_r}}
\in1+p q\mathbb Z[[q]].\tag{5}$$

No eta-product identity is needed for this unit statement. The eta identities in PR93 remain relevant to its separate low-degree power structure.

The exact norm algorithm avoids numerical roots of unity. If $qF'/F=\sum B_nq^n$ and $\mathcal N_mF=\sum C_nq^n$, then

$$C_0=1,\qquad nC_n=\sum_{i=1}^n B_{mi}C_{n-i}.\tag{6}$$

Computing through $q^K$ consumes the original series through $q^{mK}$. The implementation rejects an insufficient input precision. Independent normalized resultants of finite polynomials verify (6).

## 4. The framed inert prime-power trace theorem

Retain this entire source, not only a rational modular value:

- an Atkin–Lehner extension $\Gamma$ of $\Gamma_0(N)$ with normalized Hauptmodul $t=q+\cdots\in q\mathbb Z[[q]]$;
- a holomorphic weight-two form $\Theta(\tau)=F(t(\tau))\in1+q\mathbb Z[[q]]$, trivial character on $\Gamma_0(N)$, with the permitted Atkin–Lehner character and unique zero at the pole of $t$, of order $\rho$ as in the imported setting;
- an imaginary quadratic source $\alpha$, order discriminant $D$, $p\nmid ND$, $p$ odd and inert, and $t_0=t(\alpha)\in\mathbb Q\cap\mathbb Z_p$;
- the **CITED full-level CM reciprocity action** on the ratios below, with the marked away-$p$ level/differential frame retained. Local units in $O_{K,p}^*$ act on primitive cyclic $p^r$ subgroups through $\mathbf P^1(\mathbb Z/p^r)$ and through inertia at the chosen $p$-adic place. The action transports the normalized ratios exactly, or up to frame/character units, hence preserves their valuations.

The last item is the same local-reciprocity mechanism used in the PR's prime-level inert proof, now applied at level $p^r$. It is not inferred merely from $t_0\in\mathbb Q$. If this full-level arithmetic input has not been supplied for another proposed modular family, the application to that family remains **CONDITIONAL** on it.

Define

$$T_r(X,t)=\prod_{\gamma=(a,b;0,d)}
\left(d^2-\frac{\Theta(\gamma\tau)}{\Theta(\tau)}X\right),
\qquad f_\gamma=\frac{d^2\Theta(\alpha)}{\Theta(\gamma\alpha)}.\tag{7}$$

The modular polynomial argument gives $T_r\in\mathbb Z[t,X]$ and
$\deg_t[X^m]T_r\le\lfloor m\rho\rfloor$: coefficients are invariant modular functions, with poles only at the marked pole of $t$; their integral $q$-expansions become integral $t$-coefficients because $t=q+\cdots$. This is the imported coefficient argument for the primitive representatives from the supplied paper. Since $n_r$ is even, the top coefficient is $U_r(t)$. Equation (5) gives

$$U_r(t)\equiv1\pmod p,\qquad U_r(t_0)\in\mathbb Z_p^*.\tag{8}$$

The finite CM value ensures $\Theta(\alpha)\ne0$; (8) then ensures that every numerator factor $\Theta(\gamma\alpha)$ is nonzero as well. The local Cartan action is transitive by §2, and inertia preserves the chosen valuation. Hence all $f_\gamma$ have the same valuation. Their product, with multiplicities and frames retained, has valuation

$$L_r=2\sum_{j=0}^r j c_j
=2r n_r-\frac{2(p^r-1)}{p-1}.$$

Therefore the exact root valuation is

$$\lambda_r=2r-\frac{2(p^r-1)}{(p-1)p^{r-1}(p+1)}.\tag{9}$$

For odd $p$, the subtracted quantity is strictly between 0 and 1. All elementary symmetric functions lie in $\mathbb Q_p$, by (7)–(8) and rational $t_0$, so their valuations are integers. The monic polynomial of these marked roots has the coefficient bounds

$$v_p(e_k(f_\gamma))\ge\lceil k\lambda_r\rceil,
\quad\text{and in particular}\quad
\sum_\gamma f_\gamma\equiv0\pmod {p^{2r}}.\tag{10}$$

The rational/unramified coefficient-value condition is important for the ceiling. Mere algebraic $p$-integrality would only imply the real-valued valuation bound before rounding. The odd-prime condition also matters: the same counting formula at $p=2,r=2$ gives $\lambda=3$, not a precision-four conclusion.

For the five weight-two sources in the imported power-structure note, the modular/integrality hypotheses are the inherited ones. The arithmetic conclusion uses the stated full-level CM reciprocity input; it does not identify an unmarked rational $t_0$ with all its modular frames.

## 5. The sharp support boundary: scalar backtracking

All integral upper triangular matrices of determinant $p^r$ split by content $p^j$ into scalar multiples of primitive matrices of determinant $p^{r-2j}$. Thus

$$1+p+\cdots+p^r=\sum_{j=0}^{\lfloor r/2\rfloor}\psi(p^{r-2j}),
\qquad\psi(1)=1.$$

At $r=2$ the extra matrix is $pI$. It fixes the modular source but its marked ratio from (7) is $p^2$, not 1. Hence

$$\operatorname{Tr}_{\rm all,p^2}=\operatorname{Tr}_{\rm primitive,p^2}+p^2,
\qquad v_p(\operatorname{Tr}_{\rm all,p^2})=2.\tag{11}$$

The last equality follows from (10), not a generic lack-of-cancellation assumption. Replacing primitive support by all determinant matrices or all composed one-step paths is therefore an invalid precision-preserving operation. Composition also has multiplicities, so it must be specified separately from either support set.

The multiplicity is explicit. On periodic functions put
$A_1H(\tau)=H(p\tau)+p^2\sum_{b\bmod p}H((\tau+b)/p)$, and let $A_2$ sum over the primitive determinant-$p^2$ matrices with the same $d^2$ weight. Direct composition gives

$$A_1^2=A_2+(p+1)p^2 I.\tag{12}$$

There are $p+1$ backtracking paths, each with scalar weight $p^2$. On $q$-series, this follows independently from $A_1=V_p+p^3U_p$, $U_pV_p=I$, and
$A_2=V_p^2+p^3V_pU_p+p^6U_p^2-p^2I$.
Thus even the composed-path trace has exact valuation 2 at an inert source; it is not the primitive trace. Removing these branches preserves the specified cyclic support, not all Hecke paths.

For $s_7$, $t_0=-1$, $p=3$, the exact values are

$$\operatorname{Tr}_{\rm primitive,9}=\frac{603742572}{25150225},\quad v_3=4;
\qquad\operatorname{Tr}_{\rm all,9}=\frac{830094597}{25150225},\quad v_3=2.$$

The primitive value also shows that a uniform modulus $p^{2r+1}$ does not follow even at this actual CM source. At the generic cusp value $t=0$, the primitive level-nine trace is $748$, a 3-adic unit: the CM/inert hypothesis cannot be discarded.

## 6. Exact $s_7$ experiment and remaining obligations

Use the inherited normalization

$$\theta=\sum_{x,y\in\mathbb Z}q^{x^2+xy+2y^2},\quad
\Theta=\theta^2,\quad
t\theta^3=q\prod_{n\ge1}(1-q^n)^3(1-q^{7n})^3.$$

The program constructs $U_r$ using (4) and the trace using coefficient extraction from $1/\Theta$. It reconstructs $U_r(t)$ and $U_r(t)\sum f_\gamma$ within their proved degree bounds and checks four additional $q$-coefficients. Here $\rho=2/3$. The exact rational evaluations are recorded in JSON; zero traces, if encountered, are retained as zero rather than assigned a finite valuation.

| $p,r$ | primitive degree | exact root slope | $v_p(\operatorname{Tr})$ at $t=-1$ |
|---|---:|---:|---:|
| 3,1 | 4 | $3/2$ | 2 |
| 3,2 | 12 | $10/3$ | 4 |
| 3,3 | 36 | $95/18$ | 7 |
| 5,1 | 6 | $5/3$ | 2 |
| 5,2 | 30 | $18/5$ | 5 |

Additional valid inert evaluations use $t=1/2,D=-12$ and $t=1/8,D=-19$, plus $(p,r)=(11,1)$. Bad level primes, ramified primes and nonintegral parameters are excluded. Higher observed valuations are not promoted to a general improvement.

The next precise obligation is an independently proved comparison between these primitive Hecke coefficients and
$F_{p^r}(t)/F_{p^{r-1}}(t^\sigma)$, with its normalization and denominator controlled to $p^{2r}$. The PR supplies only the needed prime-level comparison. Equation (10) alone does **not** supply it, and does not settle $c(np^r)\equiv c(np^{r-1})\pmod {p^{2r}}$.

Likewise the level-11 weight-one problem remains a genus-one/Fricke-cover coefficient problem. Replacing it by the finite projective circle in §2 would forget the modular zero/pole data; none of the calculations here closes that gap.

**Computational universe.** Cartan/Hilbert-90 maps: $p=3,5,7,11$, depths 1–3, every projective point and three source orbits per level. Cyclic norms: three explicit integral polynomials, degrees $1,3,5,9$, compared with normalized exact resultants. Modular traces: the six listed prime/depth pairs, with source normalization and the three declared CM parameters. Strict type, bad-prime, split/ramified-algebra, nonprimitive-vector, missing-chart and insufficient-precision hostiles are included. This is not a search over all modular families or a validation of all PR theorem claims.

Reproduce with `python -B -X utf8 04-computation/experiments/pr93_modular_bridge_20261018a.py`; add `-O` for the independent assertion-disabled replay. Output is [the exact JSON receipt](pr93_modular_bridge_20261018a.json).

Normal and optimized runs agree at **11,698 exact checks**. LF JSON SHA256:
`c093dc5502914046f276a69292bf957d7bef18e525aced77e7a08840eae30035`.
Independent peer proof/API review and optimized-versus-saved replay: **PASS**.
The review covers the elementary extension and its stated CM premise; it does
not re-prove that cited arithmetic input or the missing truncation comparison.
