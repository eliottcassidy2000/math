# Session59: three parameters, two geometric obstructions, and a Frobenius repair

Status: **PROVED** by the elementary arguments below, with cited classical zeta rationality only in §2; **FINITE-EXACT** controls in §8. Derived 2026-10-25. These are new connections to the inherited Session58 packet, not claims of literature priority. No actual numerical gap-addition relation, hourglass, magic square, or Riemann-Hypothesis result is obtained.

## 1. Inheritance, scope, and concept board

Read-only source commit: `2b96df316ec3`, Session58. Closest proved mechanisms:

- `05-knowledge/results/session58_hourglass_complement_geometry_20261025/tropical-gap-triples/RESULT.md`: the counts `N_k`, `M_k`, and disjoint primitive ancestry.
- `.../projective-support-bridge/PROJECTIVE_SUPPORT_BRIDGE.md`: the squarefree Gaussian index bijection with `P^(k-1)(F_3)` and its support/divisibility sidecar.
- `.../absolute-geometry/ABSOLUTE_GEOMETRY_BRIDGES.md`: extension compatibility and Frobenius cycle data; polynomial absolute zeta does not manufacture descent.

Canonical arithmetic hostile: at center root 65, the index triple `(1,0),(0,1),(1,1)` passes the support filter but has gaps `4056,3000,3696`, not a Schur triple. Corrected near miss: a rational ordinary generating function was never a theorem that the sequence is a Hasse–Weil extension count. Least-used sidecar: the Frobenius permutation of the three geometric points in an unordered configuration.

Portfolio: **Anchor**—distinguish arithmetic support size from extension degree; **Niche**—test torification after distinct-point deletion; **Wildcard**—restore unordered geometric counts through Frobenius cycle types.

| Live concept | Object / operation | What this pass establishes | Lost datum / cheapest test |
|---|---|---|---|
| Arithmetic candidate sequence | `k` split primes; ordinary/exponential series | Fractional zeta pole orders obstruct a single finite-field realization | Numerical gaps; center 65 |
| Fixed-dimension support geometry | Hold `k` fixed, vary `q` | A finite-type stratified Z-model with polynomial counts | This is not variation of integer primes |
| Torification | Expand in `t=q-1` | Total obstruction at `k>=4`; primitive obstruction at `k>=3` | Negative coefficient is necessary obstruction, not a classification of F1 models |
| Primitive ancestry | Delete coordinates zero in all three columns | Polynomial binomial transform, matching Session58 at `q=3` | No creation of additive relations |
| Quotient and descent | Order forgetting versus Frobenius | Exact `111`, `12`, `3` correction | Split orbit counts are not all quotient rational points |
| Absolute zeta | Fixed-`k` integer counting polynomial | A Soulé polynomial transform exists despite failed torification | No descent, spectral positivity, or RH implication |

The meaningful parameter distinctions are:

1. `k`: arithmetic number of split primes / projective dimension plus one;
2. `q`: field size in a fixed-`k` auxiliary geometric family;
3. `m`: extension degree over a fixed finite base `F_Q`, with `q=Q^m`.

The formulas below must not exchange these parameters.

## 2. Rational ordinary series is not rational Hasse–Weil zeta

The inherited counts are

```text
N_k = (21^k - 6*7^k - 3*5^k + 14*3^k - 6)/48,  k>=0,
M_0 = 0,
M_k = (20^k - 6*6^k - 3*4^k + 14*2^k)/48,      k>=1.
```

Their ordinary series are rational, as Session58 proves. Define instead the formal exponential series, with constant term one,

```text
Z_N(z) = exp(sum_(k>=1) N_k z^k/k),
Z_M(z) = exp(sum_(k>=1) M_k z^k/k).
```

Using `exp(sum_(k>=1) a^k z^k/k)=(1-az)^(-1)` gives

```text
Z_N(z)^48 = (1-7z)^6 (1-5z)^3 (1-z)^6 / [(1-21z)(1-3z)^14],
Z_M(z)^48 = (1-6z)^6 (1-4z)^3 / [(1-20z)(1-2z)^14].
```

Neither `Z_N` nor `Z_M` is rational. If, for example, `Z_N=R` were rational, the order of `R^48` at `z=1/21` would be a multiple of 48, whereas the displayed rational function has a simple pole there. Equivalently, `Z_N'/Z_N` has residue `-1/48`, impossible for the logarithmic derivative of a rational function. For `M` use `z=1/20`.

**Corollary (with classical zeta rationality).** Neither `N_k` nor `M_k`, indexed by `k>=1`, is the extension-field point-count sequence of one finite-type scheme over any finite field. The same obstruction applies to a finite disjoint union of locally closed varieties. No fixed base field was assumed in this impossibility.

**Finite-initial-repair boundary.** Changing finitely many terms multiplies `Z` by `exp(P(z))` for a polynomial `P` with zero constant term. This adds a polynomial `P'` to its logarithmic derivative and cannot change the nonintegral residue. Thus even eventual agreement with either sequence is impossible.

Mechanism: the `1/48` is legitimate for free arithmetic signed/permutation orbits at each independent `k`, but it is not an integral multiplicity of a Frobenius eigenvalue. This does not prohibit geometric realization of the count for one `k`, a fixed-`k` family in `q`, virtual formulas, stacks with other counting conventions, or an unrelated infinite dynamical system.

## 3. A genuine fixed-k geometric carrier

Let `k>=1`. Consider ordered triples `(P_1,P_2,P_3)` of pairwise distinct points of `P^(k-1)` such that in each homogeneous coordinate, the number of nonzero entries among the three points is `0`, `2`, or `3`. Call this the **support-balanced** predicate. It is independent of scalar representatives.

To be precise about scheme structure, use the following *stratified model* `X_k` over Z. For each ordered triple of nonempty supports `(S_1,S_2,S_3)` satisfying the row predicate, take the product of the three standard projective support tori, then delete each pairwise diagonal that can occur. A diagonal can occur only when its two supports coincide; it is closed in that product. Let `X_k` be the finite disjoint union of these open subschemes of tori. It is a reduced finite-type Z-scheme, and its points over every field are exactly the displayed configurations with the support stratum retained.

This construction is deliberately stratified: the support-balanced constructible locus in the ambient product need not itself be locally closed, so no unproved identification with an induced ambient subscheme is used. Forgetting the disjoint-union topology retains the field-point set, not necessarily the topology. Before deleting diagonals, this model is literally a disjoint union of split tori.

Write `C_k(q)=#X_k(F_q)`, and set `C_0(q)=0` (empty object). For all finite prime powers `q`, put `t=q-1`. Then

```text
t^3 C_k(q)
  = (1+3t^2+t^3)^k - 3(1+t^2)^k - 3t(1+t+t^2)^k
    + t(2t+3)(1+t)^k + 2 - 2t^2.                 (C)
```

This numerator is divisible by `t^3` in `Z[t]`; thus `C_k(q)` is an integer polynomial. At `q=3`,

```text
C_k(3)=6N_k.
```

The factor six orders the three projective points. The `q=3` fiber, not all fibers, is the inherited squarefree arithmetic support filter.

### Proof of (C)

Count ordered triples of vector representatives in a `k` by `3` matrix, then divide by the free scalar action of size `t^3`.

- A permitted row has `1+3t^2+t^3` choices.
- A specified zero column leaves `1+t^2` row choices. Two zero columns force the full matrix zero. After excluding zero columns, the count is `(1+3t^2+t^3)^k-3(1+t^2)^k+2`.
- For a specified pair of proportional columns, fix their nonzero ratio (there are `t` choices). If their common first entry vanishes, the third entry must vanish; otherwise the third is arbitrary. The local count is `1+tq=1+t+t^2`. Requiring the third column to be nonzero removes `q^k` matrices; it also excludes a zero common column, since that would force everything zero. Thus the pair event has size `t[(1+t+t^2)^k-q^k]`.
- The intersection of any two pair events is the all-proportional locus, of size `t^2(q^k-1)`. Inclusion–exclusion adds twice this intersection.

The result is precisely (C). Pairwise distinct projective columns make both the scalar and column-permutation actions free on rational lifts; this statement alone does not count rational points on the geometric permutation quotient (§6).

Divisibility by `t^3` can alternatively be checked by the coefficients of degrees zero, one, and two in the displayed numerator, all of which vanish. Division by a monomial then leaves an integer polynomial.

## 4. The diagonal-removal torification obstruction

Expand `C_k(1+t)` in powers of `t`. Direct coefficient extraction from (C) gives

```text
[t^0] C_k(1+t) = 0,
[t^1] C_k(1+t) = k(k-1),
[t^2] C_k(1+t) = 3*binom(k,2)-7*binom(k,3)
                = k(k-1)(23-7k)/6.              (L)
```

For the last identity, the coefficient of `t^5` in the numerator of (C) is

```text
3k(k-1) - 3[binom(k,2)+3binom(k,3)+binom(k,4)]
           + 2binom(k,3)+3binom(k,4),
```

which is the claimed value. Therefore it is negative for every `k>=4`.

**Theorem.** For every `k>=4`, no finite split-torus stratification has uniform count polynomial `C_k(q)`. In particular the stratified model `X_k` is not torifiable over Z. More strongly, for any fixed finite base `F_Q`, no split-torus-stratified finite-type scheme has extension counts `C_k(Q^m)` for every `m>=1`.

Proof: a finite decomposition into split tori has a counting polynomial `sum_r a_r(q-1)^r` with nonnegative integers `a_r`. Equality for all finite fields, or merely for the infinitely many arguments `Q^m`, forces equality of polynomials. The negative coefficient in (L) contradicts this.

First hostile:

```text
C_4(1+t) = t^9+12t^8+54t^7+109t^6+102t^5+78t^4
           +56t^3-10t^2+12t.
```

Boundary/positive control: `k=0,1` are empty, and at `k=2` the support filter is automatic on three distinct projective-line points. The usual unstratified configuration space is `Conf_3(P^1)`, isomorphic to `PGL_2`; it has count `q^3-q=t^3+3t^2+2t` and the standard split-group torification. This is a positive realization of the same count and field-point functor, **not** an assertion that the finer disjoint support model is isomorphic to `PGL_2`. The coefficients of `C_3(1+t)` are nonnegative, but we do not infer a torification from that necessary numerical test.

The source support-only model is torified; deleting proportionality diagonals causes the obstruction. This is a concrete instance of a familiar non-closure phenomenon: torifiability is not automatic for arbitrary open complements. We exclude split-torus stratifications only—not nonsplit torus decompositions, blueprint models, or all notions of F1 geometry.

## 5. Primitive ancestry strengthens the obstruction

Let `X_k^prim` be the union of the same strata with every coordinate active in at least two of the three points, i.e. row weights `2` or `3` only. Write `D_k(q)` for its count; put `D_0(q)=0`. At `q=3`, this is precisely the Session58 primitive support predicate, so `D_k(3)=6M_k`.

Every configuration has a unique union `U=S_1 union S_2 union S_3` of active coordinates. Deleting all other coordinates gives an element of `X_|U|^prim`; reinsertion is inverse. Consequently the geometric support partition gives the **polynomial identity**

```text
C_k(q)=sum_(j=0)^k binom(k,j)D_j(q).             (A)
```

For `k>=1`, repeating the proof of §3 with the all-zero row removed yields

```text
t^3D_k(q) = (3t^2+t^3)^k - 3t^(2k) - 3t(t+t^2)^k
            + (3t+2t^2)t^k.
```

Equivalently, for `k>=2`,

```text
D_k(1+t)=t^(k-2){t^(k-1)[(t+3)^k-3]-3(t+1)^k+2t+3}.   (P)
```

`D_1=0` is handled directly; no negative-power expression is needed at that endpoint. `D_2=C_2`. For every `k>=3`, the bracket in (P) has constant term zero and linear term `-(3k-2)t`, while the first term in braces begins at degree at least two. Thus

```text
D_k(1+t)=-(3k-2)t^(k-1)+O(t^k),  k>=3.          (PL)
```

The first nonzero coefficient is negative. Hence no finite split-torus decomposition realizes this primitive count, over Z or over any fixed finite base, for **every `k>=3`**. The first hostile is

```text
D_3(1+t)=t^6+9t^5+24t^4+15t^3-7t^2.
```

This strengthens the total threshold from four to three; the lower-support ancestry contributions in (A) can hide a primitive obstruction. At `k=3`, adding three copies of `C_2` changes the coefficient `-7` to `2`.

## 6. Unordered rational configurations need Frobenius cycle data

`C_k(q)/6` counts unordered triples **of F_q-rational points**. It does not count all F_q-rational points of a geometric unordered configuration quotient. The latter also includes Frobenius-stable triples whose individual points are not rational.

Let `U_k(q)` count the Frobenius-stable unordered triples of distinct geometric points satisfying the same support predicate. Their Frobenius permutation has one of three types. Define

```text
T2_k(q)=sum_(s=1)^k binom(k,s)
        [(q^2-1)^(s-1)-(q-1)^(s-1)](1+q+...+q^(s-1)),

T3_k(q)=sum_(i=0)^(k-1)(q^(3i)-q^i).
```

Then

```text
U_k(q)=C_k(q)/6 + T2_k(q)/2 + T3_k(q)/3.         (F)
```

This is an exact Frobenius-stable-set count; no general quotient-existence theorem is needed for the enumeration. When the geometric quotient of the stratified model is formed, its rational points have precisely this interpretation. The ordinary unordered projective-line configuration space in the `k=2` example below is a standard explicit quotient.

### Proof by orbit types (independent of twisted Burnside)

**Type `111`:** every point is rational, giving `C_k(q)/6`.

**Type `12`:** write the triple `{P,Frob(P),Q}`, with `P` of exact degree two and `Q` rational. The two conjugates have the same support `S`, of size `s`. If `Q` had a nonzero coordinate outside `S`, that row would have weight one. Conversely `supp(Q) subset S` is sufficient, since the conjugate pair contributes two nonzero entries on each coordinate in `S`.

The number of exact-support `S` projective points defined over `F_(q^2)` but not `F_q` is `(q^2-1)^(s-1)-(q-1)^(s-1)`. Each degree-two orbit is counted twice. There are `1+q+...+q^(s-1)` allowed rational choices of `Q`. Sum over supports and divide by two, giving `T2/2`.

**Type `3`:** a degree-three orbit has a common support and hence every row has weight zero or three; the predicate is automatic. Because three is prime, exact-degree-three points are all `P^(k-1)(F_(q^3))` points except the rational ones. Dividing their number by three gives `T3/3`.

These types are disjoint and exhaustive, proving (F). The same formula is the `S_3` twisted-Frobenius average `(C+3T2+2T3)/6`, but the direct proof avoids assuming that quotient rational points equal rational-point orbits.

### Smallest nonempty hostile

At `k=2`, all three twisted terms equal `q^3-q`. Thus

```text
split rational triples = (q^3-q)/6,
degree-two-plus-one    = (q^3-q)/2,
degree-three orbit     = (q^3-q)/3,
U_2(q)                = q^3-q.
```

For `q=3`, these are `4+12+8=24`. The original arithmetic count is `N_2=4`, **not 24**. The 20 extra geometric points are descent twists, not extra integer-gap candidates.

For reference, `U_3(3)=824` versus `N_3=150`, and `U_4(3)=22248` versus `N_4=3736`. Finite polynomial expansions of `U_k(1+t)` for small `k` happen to have nonnegative coefficients; that is not a proof of torifiability and is not used here.

## 7. Absolute zeta: a valid transform, not a descent certificate

For each fixed `k`, `X_k` and `X_k^prim` are finite-type Z-schemes satisfying an integer polynomial point-count condition. If `C_k(q)=sum_j c_j q^j`, the classical polynomial Soulé transform is therefore

```text
zeta_Ck(s)=product_j (s-j)^(-c_j).
```

Likewise use the coefficients of `D_k`. This is the normalized polynomial counting transform; it does not assert an F1 descent. The obstruction in §4 shows why polynomial-counting and split-torus descent are different properties in this very family. Since `C_k(1)=D_k(1)=0`, the formal normalization exponent at `q=1` is zero. A zero value is not a claim that the varieties are empty, nor an identification with any raw Hom-point convention.

At `k=2`, `C_2(q)=q^3-q`, so the transform is `(s-1)/(s-3)`. In contrast, treating the *varying-k* arithmetic sequence `N_k` as an extension sequence is ruled out by §2. These are different constructions, not inconsistent zeta formulas.

The coefficient `k(k-1)` in §4 equals the number of ordered unequal coordinate labels. This alone supplies no root-system action or Weyl-group theorem; no numerological bridge is inferred.

## 8. Exact controls and reproduction

Run from this directory:

```bash
python3 verify.py > normal.log
python3 -O verify.py > optimized.log
cmp normal.log optimized.log
```

The verifier uses only the Python standard library and explicit exception gates. `normal.log`, `optimized.log`, and `verification.json` are saved matching outputs.

1. Directly enumerate all unordered triples of canonical projective representatives, multiply by six for ordered counts, and test both total/primitive predicates: `q=2,k=0..5`; `q=3,k=0..4`; `q=5,7,k=0..3`. The largest direct triple universe is `binom(57,3)=29260` at `q=7,k=3`.
2. Use independent exact polynomial arithmetic through `k=20` to check (C), (P), divisibility, (A), specialization to `6N_k/6M_k`, and the low/leading coefficient formulas.
3. Construct explicit finite fields `F_(q^2)` and `F_(q^3)` for prime `q=2,3,5`; check each nonzero element satisfies `a^(q^d-1)=1`. Enumerate projective Frobenius orbits for `k=1..3`, independently enumerate eligible `12` and `3` configurations, and compare with (F). Largest extension projective set: `P^2(F_125)`, 15751 points.
4. Retain positive control `Z(A^1/F_3,z)=1/(1-3z)` and its integral logarithmic residue; retain the `-1/48` fractional-pole hostile. The all-index nonrationality theorem rests on the proof in §2, not a finite numerical approximation.
5. Retain all three minimal nonempty failures: total `k=4` coefficient `-10`, primitive `k=3` coefficient `-7`, and quotient `q=3,k=2` split4 versus geometric24.

Parent independently checked total `C_k` by a separate direct projective-triple program for `q=2,3,k=1..4`, `q=5,k=1..3`, and independently expanded the low coefficients through `k=6`. That external audit is reported separately by the parent; the present directory does not claim to contain its transcript.

## 9. Connection and failure ledger

| Source → target / map | Preserved | Destroyed / sidecar | Decisive test |
|---|---|---|---|
| Arithmetic signed indices → `X_k(F_3)` / reduce and order | Distinctness, supports, necessary valuation predicate | Actual gaps and additive phase / Gaussian factors or integer gap list | Center65 hostile |
| Varying-k sequence → proposed fixed-base extension sequence | Numerical values only | Frobenius-integral multiplicities absent | Nonintegral zeta residue; proposed map rejected |
| Support-only tori → distinct configurations / remove diagonals | Coordinate supports and projective scalar classes | Torification may fail / diagonal arrangement | `-10t^2` at k4, `-7t^2` primitive k3 |
| Total → primitive / delete inactive coordinates | Reversible support ancestry | Ambient label set / chosen active subset | Binomial polynomial identity |
| Ordered → split unordered / divide by6 | Rational-point orbit count | Nontrivial descent / Frobenius permutation | 4 versus24 |
| Split unordered → geometric unordered / add cycles12,3 | Complete stable triples | Integer arithmetic interpretation does not extend | Explicit extension-field orbit enumeration |
| Fixed-k polynomial → absolute zeta / coefficient product | Counting-polynomial coefficients | Descent and RH positivity / actual geometric/spectral data | Polynomial-count model can still fail torification |

Hypothesis outcomes: `k` as extension degree **REFUTED**; fixed-`k` polynomial geometric carrier **PROVED**; torification survives distinctness **REFUTED** in the thresholds above; unordered quotient points equal split rational orbits **REFUTED**; Frobenius repair **PROVED**; any numerical gap-addition consequence **NOT CLAIMED**.

New productive gates: (i) classify which *coarser unstratified* carriers, if any, admit torifications at the non-obstructed small values; (ii) investigate the quotient's nonnegative coefficient phenomenon without confusing it with a construction; (iii) retain full Frobenius cycle data whenever the cyclotomic support lane applies power/permutation operators. None is needed to close the present proofs.

## 10. Primary-source pins and import boundaries

- Bernard Dwork, *On the rationality of the zeta function of an algebraic variety*, American Journal of Mathematics 82 (1960), 631–648, DOI `10.2307/2372974`. Primary publication/translation record: https://m.mathnet.ru/php/archive.phtml?jrnid=mat&option_lang=eng&paperid=201&wshow=paper . Imported only the classical rationality theorem; its proof was not re-audited. The theorem's consequence is applied to §2, not used to invent a geometric model.
- Javier López Peña and Oliver Lorscheid, *Torified varieties and their geometries over F1*, arXiv:0903.2173v3, §1.1 definitions and Lemma1.2, §1.2 counting/zeta, §1.3.6 split Chevalley examples: https://arxiv.org/html/0903.2173v3 . Imported finite split-torus decomposition and resulting field-point counting, and the standard split-group positive control. We do not use the converse “polynomial implies torification”; no such converse is claimed.
- Anton Deitmar, Shin-ya Koyama, Nobushige Kurokawa, *Counting and zeta functions over F1*, arXiv:1306.3681v2, §1: https://arxiv.org/html/1306.3681v2 . Imported the all-prime integer-polynomial condition and coefficient-product normalization. The candidate formulas, torification obstructions and Frobenius corrections here are not asserted as results of that paper.
- Weiyan Chen, *Twisted cohomology of configuration spaces and spaces of maximal tori via point-counting*, arXiv:1603.03931, https://arxiv.org/abs/1603.03931 . Bibliographic context only: twisted point-counting of configuration spaces is an established method. No theorem from its abstract is a dependency of the direct cycle-type proof in §6.

All primary links checked 2026-10-25. No novelty claim is made for torification, the zeta formalism, or Frobenius orbit counting; the local contribution is their exact application to the Session58 candidate carriers and the explicit failure thresholds.

Operational note: all new files were written only in `/tmp/session59-f1-arithmetic`. No repository publication was attempted by this lane. The inherited checkout had unrelated dirty inbox paths. The bounded inbox tool reported `.machine-id not found`; no registration file was fabricated and no inbox mutations were attempted.
