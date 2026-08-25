# Pythagorean-tree source sidecar

> **Freshness:** primary records checked 2026-08-12.  This is a routed
> extension of `CORE-PAPERS.md`, not an independent truth surface.

## Square-pyramidal / cannonball classification

### Bennett -- *Lucas' Square Pyramid Problem Revisited*

- **Primary:** Michael A. Bennett, *Acta Arithmetica* **105** (2002),
  341--347, [DOI 10.4064/aa105-4-3](https://doi.org/10.4064/aa105-4-3),
  [publisher record](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/105/4/83879/lucas-square-pyramid-problem-revisited),
  and [author PDF](https://personal.math.ubc.ca/~bennett/paper21.pdf).
- **Imported role:** Bennett's Theorem 2.1 bounds the positive solutions of
  `x(x+1)(x+2)=6y^2` by three; the paper's three displayed solutions saturate
  that bound.  Restricting them to `x=2s,y=2t` proves that the positive integer
  solutions of

  ```text
  1^2+2^2+...+s^2=t^2
  ```

  are exactly `(s,t)=(1,1)` and `(24,70)`.  The paper also records Watson's
  1918 proof and gives a modern algebraic route through quartic Pell-type
  equations.
- **Consumer:**
  [THM-3335](../../01-canon/theorems/THM-3335-square-triangular-pell-markov-pythagorean-selector.md),
  where the classification makes the selector/cannonball intersection global:
  the unique positive selector root that is square-pyramidal is `s_3=70`, at
  pyramid height `24`.
- **Does not prove:** THM-3335's Pell/Markov/Pythagorean compiler, its
  even/odd split, the separate question whether `q_k` is triangular, any
  skew-EW design, or an LRC/tournament/JC transfer.

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

## 2026-08-24 repo-derived rational-edge update

[THM-4057](../../01-canon/theorems/THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge.md)
is **PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED** and is derived
from elementary Stern/Calkin--Wilf recurrences plus prior repo carriers. It
identifies primitive-Pythagorean Calkin--Wilf ordinals as `3,5 mod 6`, gives
the exact three Berggren ordinal transducers and radix inverse, and classifies
marked branch-triangle closure gcds and the unique unordered-Pythagorean root
collision. The sources above do not supply those ordinal/tournament claims;
THM-4057 does not claim they prove a tournament, LRC, or Jacobian result.
