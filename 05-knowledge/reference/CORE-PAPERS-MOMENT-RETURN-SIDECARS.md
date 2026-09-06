# Moment-return and factorial sidecar source pins

These source records were moved intact from CORE-PAPERS.md on 2026-09-05
to keep the mandatory startup surface bounded. Original freshness labels
remain in force; the move does not revalidate a paper or upgrade its status.

<a id="strong-factorial"></a>

### Edo--van den Essen — *The Strong Factorial Conjecture*

- **Primary / freshness:** [arXiv:1304.3956v2](https://arxiv.org/abs/1304.3956);
  *J. Algebra* **397** (2014), 443--456, [DOI](https://doi.org/10.1016/j.jalgebra.2013.09.011).
  **PUBLISHED / stable; source and journal record checked 2026-07-24.**
- **Imported role / notation gate:** Definitions 2.1/2.7 and Conjectures 2.4/2.8 index `FC(m)`/`SFC(m)` by ambient variables; `N(f)` is support/window size. Thus `t` univariate terms mean `SlotSFC_1(t)`, not `SFC(t)` (MISTAKE-350).
  The Gaussian identity
  `E[(W+Zh)^(2j)]=binom(2j,j)L((sh)^j)` identifies its `n=1` two-charge case.
- **Repo consumers:** [THM-1790](../../01-canon/theorems/THM-1790-the-emp-floor-detection-depth-at-least-degree-plus-one.md)
  and [HYP-8765](../hypotheses/HYP-8765-gmc2-radial-channel-return-tower.md).
- **Does not prove:** SFC, the full HYP-8765 cutoff, or a uniform NC2 depth.
  THM-2022 instead uses a coefficient-dependent good prime; THM-1790 proves
  the complementary projective lower bound.

<a id="effective-returns"></a>

### Erman--Smith--Várilly-Alvarado — *Laurent polynomials and Eulerian numbers*

- **Primary / freshness:** [arXiv:0908.2609v2](https://arxiv.org/abs/0908.2609),
  *Journal of Combinatorial Theory A* **118** (2011), 396--402,
  [DOI 10.1016/j.jcta.2010.02.006](https://doi.org/10.1016/j.jcta.2010.02.006).
- **Imported role:** studies generic constant-term equations through toric
  geometry and Eulerian numbers and records the sharp Sturmfels effective
  question.  In the repo this is the target “first nonzero constant term by
  exponent-width `m+n`,” with Newton-polygon branch count as its geometry.
- **Repo consumers:**
  [THM-1630 effective redirect](../../01-canon/theorems/THM-1630-tnc-is-duistermaat-van-der-kallen-theorem-2.md),
  [THM-1650 Newton polygon](../../01-canon/theorems/THM-1650-newton-polygon-of-the-effective-dvdk-bound.md).
- **Does not prove:** the Sturmfels `m+n` effective bound, NC2, or GMC(2).
  Generic regular-sequence/degree information is not a uniform first-return
  theorem for every Laurent polynomial.

<a id="hyper-bessel"></a>

### Baricz--Singh — *Zeros of some special entire functions*

- **Primary / freshness:** [arXiv:1702.00626v2](https://arxiv.org/abs/1702.00626),
  *Proceedings of the AMS* **146** (2018), 2207--2216,
  [DOI 10.1090/proc/13927](https://doi.org/10.1090/proc/13927).
- **Imported role:** the positive-parameter hyper-Bessel / `0F_q` zero theorem
  places all zeros on the negative real axis.  Gauss multiplication makes the
  repo's `Phi_(p,q)(x)=sum x^k/((pk)!(qk)!)` an exact instance.
- **Repo consumer:**
  [THM-2023, Laguerre--Polya boundary theorem](../../01-canon/theorems/THM-2023-gmc2-hyperbessel-boundary-laguerre-polya.md).
- **Does not prove:** GMC or NC2.  It handles the `Phi` sharp boundary, not the
  opposite Mittag--Leffler/root-of-unity-filter `Psi` boundary, and it does not
  exclude cancellation at an allowed negative-real zero without further work.

<a id="legendre-radar"></a>

### Mangoubi--Kadets--Weller Weiser — common roots of Legendre polynomials

- **Primary / freshness:** [IAS seminar announcement, 2026-02-18](https://www.ias.edu/math/events/special-year-learning-seminar-41).
  **SEMINAR ONLY:** no public paper or proof was located on 2026-07-21.
- **Imported role:** the announcement states that only finitely many Legendre
  polynomials vanish at a fixed nonzero point, with quantitative bounds to be
  discussed, and attributes the work to Dan Mangoubi, Borys Kadets, and Adi
  Weller Weiser.
- **Repo consumer:**
  [THM-2021, Legendre refinement](../../01-canon/theorems/THM-2021-gmc2-legendre-finite-recurrence-closure.md).
- **Does not prove:** a citable published finite-recurrence
  theorem until a paper/proof is public.  It is no longer a dependency of
  NC2: [THM-2018's recurrence mechanism](../../01-canon/theorems/THM-2018-gmc2-resonance-algebraic-egf-and-proportional-shadow.md)
  and [THM-2022's full proof](../../01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md)
  supersede that need.  It also does not establish Stieltjes' stronger “at
  most one” claim.
