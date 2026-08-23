# Core papers: Mahler's `3/2` problem and rational bases

**Literature status checked 2026-08-23.**  This is the routed source ledger for
the compact entry in [CORE-PAPERS](CORE-PAPERS.md).  The problem remains
**OPEN**.  Repo deductions are typed separately from cited results.

## Mahler — *An unsolved problem on the powers of `3/2`*

- **Primary / published:** [EMS reprint and bibliographic record](https://ems.press/books/dms/252/4994),
  original *J. Australian Math. Soc.* **8** (1968), 313--321.
- **CITED result:** Mahler defines a `Z`-number by
  `{xi(3/2)^n}<1/2` for every `n>=0`, proves that the candidate set is at most
  countable, and gives the quantitative counting bound `O(x^0.7)` up to
  height `x`.
- **Repo consumers:**
  [THM-2228](../../01-canon/theorems/THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization.md)
  and [THM-3848](../../01-canon/theorems/THM-3848-rational-base-prefix-atom-tree-and-lonely-runner-separation.md).
- **Does not prove:** existence or nonexistence of a `Z`-number.  A finite
  safe-prefix family, even at every horizon, is not one fixed positive orbit.

## Flatto--Lagarias--Pollington — fractional-part range

- **Primary / published:** [Acta Arithmetica record and paper](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/70/2/108521/on-the-range-of-fractional-parts-p-q),
  **70** (1995), 125--147, DOI `10.4064/aa-70-2-125-147`.
- **CITED result used here:** every `xi>0` has
  `limsup {xi(3/2)^n}-liminf {xi(3/2)^n}>=1/3`.
- **Does not prove:** the `1/2` range needed to exclude Mahler `Z`-numbers,
  density, or equidistribution for an individual `3/2` orbit.

## Akiyama--Frougny--Sakarovitch — rational-base numeration

- **Primary / published:** [author PDF](https://math.tsukuba.ac.jp/~akiyama/papers/3half_H65fullTH-060406.pdf),
  *Israel J. Math.* **168** (2008), 53--91,
  DOI `10.1007/s11856-008-1056-4`.
- **CITED result:** least-significant-first base-`p/q` expansions are unique
  for integers; the integer-expansion language is nonregular, while addition
  admits a finite right-transducer.
- **Repo separation:** THM-3848's normalized torus wall partition keeps an
  explicit conversion map.  Its `p^N` atom count is neither the count of
  rational-base integer words nor the count of THM-2228 safe carry words.
- **Does not prove:** Mahler's conjecture, normality of minimal words, or an
  identification of the distinct finite-prefix trees.

## Dubickas--Mossinghoff — *Lower bounds for Z-numbers*

- **Primary / published:** [AMS article record](https://www.ams.org/journals/mcom/2009-78-267/S0025-5718-09-02211-X/),
  *Math. Comp.* **78** (2009), 1837--1851,
  DOI `10.1090/S0025-5718-09-02211-X`.
- **CITED result:** exact algorithms prove that any `Z_(3/2)`-number would
  exceed `2^57` and connect safe prefixes to a ceiling iteration on integers.
- **Does not prove:** nonexistence.  THM-3848's all-horizon denominator-19
  family changes its initial integer and does not evade or strengthen the
  fixed-height lower bound.

## Andrieu--Eliahou--Vivion — rational-base normality

- **RADAR / PREPRINT v2:** [arXiv:2510.11723v2](https://arxiv.org/abs/2510.11723v2),
  checked 2026-08-23.
- **CITED conditional route:** the paper conjectures normality of all minimal
  and maximal rational-base words, proves equivalence with residue-class
  equidistribution for every ceiling orbit, and proves that the conjecture
  would imply nonexistence of `Z_(p/q)`-numbers for `1<q<p<q^2`.
- **Does not prove:** its normality/equidistribution conjecture or Mahler's
  problem.  For `3/2`, even the weaker finite-word occurrence needed by that
  route remains conjectural; numerical evidence is not a terminal sidecar.

## Typed repo synthesis

- **PROVED:** THM-2228 identifies ordinary stabilization as a separate gate.
- **PROVED + VERIFIED-EXACT:** THM-3848 solves the normalized finite-prefix
  atom tree and the associated mixed-power lonely-runner row.
- **PROVED:** THM-3848's closed formal safe-tail shift is nonsofic, with
  entropy `log(3/2)`; this is not the unknown stabilized Z-language.
- **CONDITIONAL:** THM-3848 derives radical growth from an explicitly assumed
  ABC schema.  It gives no safe-horizon or stabilization bound.
- **OPEN:** whether any positive `xi` satisfies the strict inequality for all
  powers of `3/2`.
