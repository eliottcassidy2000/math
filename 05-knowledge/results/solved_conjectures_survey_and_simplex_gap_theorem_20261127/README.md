# Transfer survey: from shifted Waring offsets to an exact simplex semigroup in every degree

I scanned the 200 entries returned from the [solved-conjectures results directory](https://github.com/1337isnot1337/solved-conjectures/tree/main/results) and read the relevant proof files, not treating the directory's title as proof of every claim. This note proves a new companion result for the repo's figurate program. It does not claim to solve a second named conjecture in the source repository or the fixed-five Pollock conjecture.

## 1. Exact theorem: the stable positive-simplex gaps are always symmetric

For \(d\ge2\), the positive simplex numbers are \(\binom{d+n}{d}\), \(n\ge0\). Let an offset \(N\ge0\) be **eventually representable** if \(j+N\) is a sum of exactly \(j\) positive simplex numbers for every sufficiently large \(j\) (depending on \(N\)). Subtracting the minimum term 1 and putting it back as padding proves that this is equivalent to membership in
\[
S_d=\left\langle B(n): n\ge1\right\rangle,\qquad
B(n)=\binom{d+n}{d}-1.
\]
This is the shifted-atom move from the [Benfield–Lippard parity analysis](https://github.com/1337isnot1337/solved-conjectures/blob/main/results/benfield-lippard-waring-gap-parity/PROOF.md) and its [separate symmetry obstruction](https://github.com/1337isnot1337/solved-conjectures/blob/main/results/benfield-lippard-waring-gap-symmetry/PROOF.md), applied to the [repo's figurate families](../figurates_local_global_carry_and_norm_20261126/README.md). Their shifted *powers* and our shifted *simplex numbers* have different arithmetic.

**Theorem.** List in strictly increasing order all prime powers \(t_i=p_i^{j_i}\) with \(j_i\ge1\) and \(t_i\mid d\), and put \(g_i=B(t_i)\). Then the minimal generators of \(S_d\) are exactly \(d,g_1,\ldots,g_\ell\). The least element in each residue modulo \(d\) is given by the \(d\) distinct values
\[
\sum_{i=1}^{\ell} b_i g_i,\qquad 0\le b_i<p_i .
\]
Consequently the largest gap and the total number of gaps are
\[
F=\sum_i(p_i-1)g_i-d,\qquad G=(F+1)/2 .
\]
More precisely, for every \(0\le N\le F\), exactly one of \(N,F-N\) is representable. This reflection symmetry holds for **every** simplex degree, including mixed-prime degrees; no empirical exponent cutoff is used.

## 2. Why the carry occurs at prime powers, and why it pays for itself

Set \(d_i=\gcd(d,g_1,\dots,g_i)\) with \(d_0=d\). The product identity
\[
1+B(n)=\prod_{u=1}^n(1+d/u)
\]
is read in \(\mathbb Q_q\) for each prime \(q\mid d\). For \(q^a\parallel d\) and \(n<q^a\), every factor is \(1\bmod q^{a-\lfloor\log_q n\rfloor}\). Thus the current divisibility can drop only when \(n\) passes a listed threshold. At \(t_i=p_i^{j_i}\) the term \(u=t_i\) is the only one of its valuation, giving \(v_{p_i}(B(t_i))=v_{p_i}(d)-j_i\), **also** when \(j_i=v_{p_i}(d)\): all other factors are \(1\bmod p_i\) and the exceptional factor is \(1+d/t_i\). It follows that \(d_i=d_{i-1}/p_i\). More generally if \(t_i\le n<t_{i+1}\) (no upper restriction at the last threshold) then \(B(n)\) is divisible by \(d_i\). If additionally \(n<2t_i\), then \(B(n)\equiv g_i\bmod d_{i-1}\), since \(u=t_i\) is still the only term of that valuation for \(p_i\), and both numbers have all the previous divisibilities at the other primes.

Each residue multiple of \(d_i\) has exactly one digit expression in \(g_j\) for \(j\le i\) with \(0\le b_j<p_j\): the new \(g_i\) generates the cyclic quotient of size \(p_i\) modulo \(d_{i-1}\), and the remaining residue is decoded by induction. Let \(W_i=\sum_{j\le i}(p_j-1)g_j\) and \(W_0=0\). We need these residue representatives to lie **below** the actual target, not just in the correct classes.

If \(s<t\) are adjacent listed thresholds, then \(B(t)\ge(p_s+1)B(s)\), where \(p_s\) is the prime of \(s\). If they have different bases, \(d\ge st\), and even the last factor \(1+d/t\) in the ratio of binomial coefficients is at least \(s+1\ge p_s+1\). If they have the same base, then \(t\ge p_s s\), and each factor in that ratio is at least 2, since \(t\le d\). Induction now gives \(W_{i-1}<g_i\), and for \(i>1\) even \(W_{i-1}<p_{i-1}g_i/(p_{i-1}+1)\).

To generate \(B(n)\), the initial range before \(t_1\) is divisible by \(d\). Otherwise use the current \(t=t_i\).
- At \(n=t\) we already have \(g_i\). For \(t<n<2t\), its top digit is 1. The jump \(B(t+1)-B(t)=d(1+B(t))/(t+1)\) exceeds \(W_{i-1}\), by the preceding bound and \(d\ge t\ge p_{i-1}\) (or \(W_0=0\)).
- For \(n\ge2t\), the ratio of binomial coefficients at \(2t,t\) is at least \((3/2)^t\ge t\ge p_i\). So \(B(n)\ge p_i g_i > (p_i-1)g_i+W_{i-1}\), enough for every digit choice in this range.

In both cases subtract the digit representative and pay the nonnegative residual with copies of \(d\). This proves generation; the gcd drop and increasing order prove minimality. Moreover \(p_i g_i\) is a nonnegative combination of earlier generators: it lies in their residue subgroup, and that subgroup's digit cost is at most \(W_{i-1}<g_i\). Thus every sum can be reduced to its unique digit form plus a nonnegative multiple of \(d\). This proves the claimed Apéry set (least representatives), not just an upper bound for it. Digitwise complement sends \(w\) to \(W_\ell-w\), so \(N=w+kd\) is representable exactly when \(k\ge0\), while \(F-N=(W_\ell-w)-(k+1)d\) is representable exactly when \(k\le-1\). It also gives the exact generating function
\[
\sum_{N\in S_d}z^N=\frac{1}{1-z^d}\prod_i\frac{1-z^{p_i g_i}}{1-z^{g_i}} .
\]
Unlike a prime-factor Euler product, these factors enumerate an additive, normalized carry.

| Degree | Generators | Last stable gap | Number of stable gaps |
|---:|---|---:|---:|
| 2 | 2, 5 | 3 | 2 |
| 3 | 3, 19 | 35 | 18 |
| 4 | 4, 14, 69 | 79 | 40 |
| 8 | 8, 44, 494, 12869 | 13399 | 6700 |
| 12 | 12, 90, 454, 1819 | 2805 | 1403 |
| 30 | 30, 495, 5455, 324631 | 1309899 | 654950 |

The tetrahedral row means \(j+N\) is eventually a sum of exactly \(j\) **positive** tetrahedra iff \(N\in\langle3,19\rangle\). In each residue mod 3 the first possible offsets are \(0,19,38\). This does not bound the required \(j\) by five. Indeed \(5+35=40\) is not a sum of five positive tetrahedra under this offset law but is a sum of two twenties, so the distinction from an at-most-five theorem is already real. The independent replay computes Apéry sets from all shifted atoms up to a certified stopping point, then compares each residue and each generator to the formula for every \(2\le d\le150\).

## 3. Other transferable ideas identified in the survey

| Source lane | What the proof preserves | Next repo test, not a claim of resolution |
|---|---|---|
| [Abelian voltage covers](https://github.com/1337isnot1337/solved-conjectures/blob/main/results/erdos-gyarfas-abelian-covers/PROOF.md) and [affine polarity edit barriers](https://github.com/1337isnot1337/solved-conjectures/blob/main/results/erdos-85-four-cycle-threshold/PROOF.md) | Tree gauges reduce a lifted route to a chord sum; potential vertex collisions are explicit failures, not just trace deficits. | Apply the chord-sum test to the remaining [affine 125-panel](../moore_125_affine_panels_and_nemenyi_inverse_20261127/README.md) stratum. A partial route either uses distinct endpoints or supplies an exact forbidden equality in the phase variables. |
| [One-intersection transversals](https://github.com/1337isnot1337/solved-conjectures/blob/main/results/blazsik-lemons-one-intersection/PROOF.md) and [multicolor sails](https://github.com/1337isnot1337/solved-conjectures/blob/main/results/multicolor-steiner-sail/PROOF.md) | A cover on a colored neighborhood is constrained by the incidence at its central point. | In the conditional [Moore matching design](../general_position_moore_design_and_pollock_separator_20261127/README.md), first ask for an independent transversal of endpoint pairs under actual block intersections. The no-one-point-intersection hypothesis is *not* automatic there, so do not import the coloring conclusion. |
| [Two-speed insertion](https://github.com/1337isnot1337/solved-conjectures/blob/main/results/lonely-runner-two-insertions/PROOF.md), [eight-speed spectrum](https://github.com/1337isnot1337/solved-conjectures/blob/main/results/eight-speed-spectrum/PROOF.md), and [fixed relation hyperplanes](https://github.com/1337isnot1337/solved-conjectures/blob/main/results/loneliness-on-a-hyperplane/PROOF.md) | The usable certificate is a native interval or torus relation, rather than an average over samples. | For the Fermat-grid blind cases, start from exact surviving chord intervals before adding speeds. The [prime-grid example](../fermat_dihedral_tournaments_and_lonely_chord_grids_20261127/README.md) shows why a modular optimum cannot replace that interval calculation. |
| [Thue–Morse valuation](https://github.com/1337isnot1337/solved-conjectures/blob/main/results/gawron-miska-ulas-thue-morse-valuation/PROOF.md), [CUT periodicity](https://github.com/1337isnot1337/solved-conjectures/blob/main/results/cut-arithmetic-periodicity/PROOF.md), and [pattern packing cutoffs](https://github.com/1337isnot1337/solved-conjectures/blob/main/results/blitvic-kammoun-steingrimsson-packing-cutoff/PROOF.md) | Exact recurrences can lift a defect, but finite initial agreement alone is not a uniform theorem. | Extend the [Pollock dyadic carry](../thue_morse_fermat_lonely_grid_20261127/README.md) by recording the carry state, not just finitely many coefficients of its centered series. For Chan–Pak a residue digit still has to come with a literal fixed-height lift. |

Run `node replay.js` for the exact simplex checks. The proved result here is the all-degree shifted-simplex theorem; the cross-problem tests are deliberately flagged separately.
