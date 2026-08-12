# Pythagorean-tree source sidecar

> **Freshness:** primary records checked 2026-08-12.  This is a routed
> extension of `CORE-PAPERS.md`, not an independent truth surface.

## Berggren descendant geometry and fixed-hypotenuse fibres

### Kőszegyová--Csókási--Hirjak / Janičková--Csókási

- **Primary:** Lucia Kőszegyová, Evelin Csókási, Juraj Hirjak,
  *Structure of Primitive Pythagorean Triples in Generating Trees*,
  *European Journal of Pure and Applied Mathematics* **17** (2024),
  2127--2141, [DOI](https://doi.org/10.29020/nybg.ejpam.v17i3.5323) and
  [open PDF](https://www.ejpam.com/ejpam/article/download/5323/1687).
  The matrix-power and area antecedent is Lucia Janičková and Evelin
  Csókási, *Metric properties in Berggren tree of primitive Pythagorean
  triples*, [arXiv:2304.05230](https://arxiv.org/abs/2304.05230).
- **Imported role:** the three descendants of a primitive parent are
  noncollinear, lie in `2x+2y-3z+c=0`, and form a triangle of area
  `2 sqrt(17) ab`; the sources give the U-power/side-vector formulas and
  prove the triangle is neither right nor isosceles.
- **Consumer:**
  [THM-3334](../../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md).
- **Does not prove:** THM-3334's fixed-cusp Farey fan, horocycle translation,
  `Q` encoding, complete `8/9`--`9/8` angle law, Gaussian torsor, unbounded
  CRT fibres, ancestry/XOR atlas, or any tournament/LRC/JC transfer.  The
  published paper explicitly leaves acute-versus-obtuse behaviour open.

### Tripathi -- *On Pythagorean Triples Containing a Fixed Integer*

- **Primary:** Amitabha Tripathi, *Fibonacci Quarterly* **46/47** (2008/09),
  331--340, [official PDF](https://www.fq.math.ca/Papers1/46_47-4/Tripathi.pdf)
  and [DOI](https://doi.org/10.1080/00150517.2008.12428142).
- **Imported role:** for odd `c`, the number of primitive Pythagorean triples
  with hypotenuse `c` is `2^(omega(c)-1)` exactly when `c>=3` has no prime
  divisor `3 mod 4`, and is zero otherwise.
- **Consumer:** THM-3334, which refines the count to an affine Gaussian
  factor-choice torsor and proves unbounded U-spine fibres by CRT.
- **Does not prove:** the count without its prime condition (`c=3` is the
  minimal hostile), Berggren ancestry words, XOR matchings, or a tournament.
