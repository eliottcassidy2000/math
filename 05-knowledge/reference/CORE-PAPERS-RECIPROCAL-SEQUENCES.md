# Core papers sidecar: reciprocal integer sequences

> Moved verbatim from `CORE-PAPERS.md` on 2026-08-03 (boxeph) to free shared
> startup bytes. Same maintenance rule as the main file: refresh primary
> link/version, name the first consumer, sharpen **does not prove**.

### Downey--Ong--Sellers — *Beyond the Basel Problem: Sums of Reciprocals of Figurate Numbers*

- **Primary / freshness:** [author-hosted preprint](https://www.d.umn.edu/~jsellers/downey_ong_sellers_cmj_preprint.pdf),
  [archived copy](https://web.archive.org/web/20130529032918/http://www.math.psu.edu/sellersj/downey_ong_sellers_cmj_preprint.pdf),
  *College Mathematics Journal* **39** (2008), 391--394,
  [JSTOR 27646686](https://www.jstor.org/stable/27646686).
- **Imported role:** supplies the classical telescoping program and exact
  reciprocal sums for figurate-number families, including the triangular
  identity `sum 1/T_n=2`.  These are the seed rows for the repo's polygonal,
  simplex, digamma, and trigamma mass surfaces.
- **Repo consumers:**
  [THM-2000 support-harmonic/figurate surface](../../01-canon/theorems/THM-2000-support-harmonic-abel-dini-figurate-surface.md),
  [THM-2005 support-Dirichlet atlas](../../01-canon/theorems/THM-2005-support-dirichlet-automatic-tournament-atlas.md).
- **Does not prove:** the repo's support-versus-multiplicity collision law,
  Abel--Stieltjes/Dini iff criterion, iterated Bertrand boundary, full
  support-Dirichlet profile, or tournament reciprocal atlas.

### Applegate--Pol--Sloane — *The Toothpick Sequence and Other Sequences from Cellular Automata*

- **Primary / freshness:** [arXiv:1004.3036v2](https://arxiv.org/abs/1004.3036),
  *Congressus Numerantium* **206** (2010), 157--191.
- **Imported role:** develops the toothpick cellular automaton and its dyadic
  recurrences and product-like generating functions.  The repo uses the exact
  A139250 formula as an automatic-sequence test case and as a model for
  operation-generated tournament-adjacent sequences.
- **Repo consumers:**
  [THM-2000](../../01-canon/theorems/THM-2000-support-harmonic-abel-dini-figurate-surface.md),
  [toothpick bridge tangent](../../00-navigation/TANGENTS.md),
  [Perron/toothpick experiment](../../04-computation/perron_toothpick_klein_S315.py).
- **Does not prove:** convergence or a closed form for the support reciprocal
  mass, a tournament interpretation, or novelty of any repo-generated
  sequence.  Cellular-automaton multiplicities must still be collapsed before
  applying the support-harmonic lens.
