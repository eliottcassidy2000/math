---
title: "Endpoint tournaments, primitive Pythagorean norm fibers, and an elementary-symmetric Eisenstein sum"
date: 2026-10-19
status: "PROVED explicit fiberwise synthesis; INHERITED analytic identities; FINITE-EXACT controls"
script: 04-computation/experiments/endpoint_pythagorean_fibers_20261019.py
output: 05-knowledge/results/endpoint_pythagorean_fibers_20261019.json
---

# Endpoint tournaments and Pythagorean norm fibers

## Result

There is an exact connection between two of the repository's previously
separated E meanings. It passes through marked prime allocations:

\[
\boxed{\displaystyle
\operatorname{PPT}(C)\ \simeq\
\frac{\operatorname{Ham}(E_{r+2})\setminus\{I\}}
     {P_S\sim P_{S^c}},\qquad r=\omega(C),}                 \tag{1}
\]

provided `C>1` and every prime factor of `C` is `1 mod 4`. Here PPTs have
positive legs, identified under interchange; `I` is the increasing Hamiltonian
path. The bijection requires the ordered prime powers of `C` and selected
Gaussian prime factors. It is explicit in both directions below. It gives
`1,2,4` primitive shapes from the nonexceptional address quotients of
`E3,E4,E5`, respectively, at `C=5,65,1105`.

Summing the same fibers by norm yields

\[
\boxed{\displaystyle
\sum_{r\ge1}2^{r-1}e_r\!\left(\left\{\frac1{p^s-1}\right\}_{p\equiv1\ (4)}\right)
=P(s)=\frac{E(i,s)}{4(1+2^{-s})}-\frac12,
\qquad\Re s>1.}                                        \tag{2}
\]

The `e_r` are elementary symmetric functions of the **specified split-prime
sequence**. The Eisenstein family uses the normalization in the
[E3/E4/E5 atlas](e_atlas_20261019.md), and `P(s)` is the primitive-hypotenuse
series in [the odd-zeta note](ap_square_odd_zeta_20261019.md). The integer `r`
and complex variable `s` remain independent. Formula (2) is a sum over all
prime-support sizes, not a pointwise identification of E3 with ζ(5), E4 with
ζ(7), or E5 with ζ(9).

## 1. Inheritance, scope, and the live coordinates

The two closest proved mechanisms are the ordered Boolean path decoder in
[tournament path defects, §3](tournament_path_defects_20261014e.md) and the
Gaussian factor-choice torsor of
[THM-3334, §5](../../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md).
The latter already identifies a fixed-hypotenuse fiber with an affine
`F_2^r/<(1,...,1)>` torsor. The present contribution joins these two codecs,
then retains the prime-support grading in the established PPT Euler product.
There is no historical-priority claim for the classical counting formula.

The concept board is: ordered Hamiltonian paths; conjugate Gaussian prime
allocations; primitive Euclid parameters; signed square-sum labels;
elementary symmetric coefficients; Gaussian-point Eisenstein evaluation.
The anchor is the explicit path/PPT map; the niche is its support-size
Dirichlet refinement; the wildcard is its composition with signed AP
square-sum triangles.

The canonical hostile to dropping the prime restriction is `C=3`, which has
no primitive PPT. The corrected near miss is treating a permutation of
Boolean addresses as a graph automorphism. The least-used required data are
the **prime exponents and Gaussian phases**: counts ignore them, while the
actual sides and positivity do not.

All new unbounded statements below have proofs. The standard facts that
`Z[i]` is a unique factorization domain and that `p≡1 mod 4` splits are also
covered in Keith Conrad's
[*The Gaussian Integers*, §§6,8,9](https://kconrad.math.uconn.edu/blurbs/ugradnumthy/Zinotes.pdf).
The analytic normalization is inherited from the linked notes and rederived
in §5. Finite checks validate conventions, not the unbounded claims.

## 2. The tournament codec and its complement action

For `n≥3`, `E_n` has vertices `0,...,n−1` and arcs `i→j` for `i<j`, except
that `0→n−1` is replaced by `n−1→0`. Set `r=n−2`, with interior coordinate
set `[r]={1,...,r}`. For `S⊆[r]`, define

\[
P_S=(S\text{ increasing}),\ r+1,\ 0,\ (S^c\text{ increasing}).       \tag{3}
\]

Every displayed list is a directed Hamiltonian path. A path that avoids the
only backward edge must be increasing, so is `I=(0,1,...,r+1)`. Any other
path uses `r+1→0`; all entries before and after this edge must form increasing
runs. Therefore (3) is a bijection between subsets and nonexceptional paths,
and

\[
H(E_{r+2})=1+2^r.                                      \tag{4}
\]

The inverse reads the vertices before `r+1`. Thus address complement defines
an actual involution of the nonexceptional **path set**, by decoding `S^c`
again in increasing order. For `r≥1` it has no fixed path and has `2^{r−1}`
orbits.

This action is not induced by an automorphism of the tournament, including
at `E3`. Indeed, if a vertex automorphism induced it, its square would fix an
ordered Hamiltonian path, hence every vertex; the automorphism would be a
nonidentity involution. But a nonidentity involution of a finite vertex set
swaps some pair, reversing the oriented arc between them, which no tournament
automorphism can do. Equivalently, extend the action by fixing `I`: any
vertex automorphism fixing that ordered path is already the identity.
Literal reversal of a path is not this operation either. The ordered
address decoder is indispensable.

## 3. Explicit Gaussian decoder, inverse, and all symmetries

### Theorem 1: marked fiber bijection

Let

\[
C=\prod_{j=1}^r p_j^{e_j}>1,\qquad p_1<\cdots<p_r,\quad
p_j\equiv1\pmod4,\quad e_j\ge1.
\]

For each prime choose `π_j=a_j+ib_j` with norm `p_j`. One reproducible
choice is `a_j>b_j>0`. For `S⊆[r]`, form

\[
\alpha_S=\prod_{j\in S}\pi_j^{e_j}
          \prod_{j\notin S}\overline{\pi_j}^{e_j}=x+iy.                \tag{5}
\]

Set `n=max(|x|,|y|)`, `m=min(|x|,|y|)`, and return

\[
(a,b,C)=(n^2-m^2,\ 2nm,\ n^2+m^2).                   \tag{6}
\]

Here `(a,b)` is **parity ordered**: `a` is odd and `b` is even. To obtain the
size-ordered notation of the AP notes, take `(A,B)=(max(a,b),min(a,b))`.
Equations (3), (5), and (6) descend to the bijection (1).

**Proof.** The norm of (5) is `C`. Its coordinates are relatively prime:
a rational prime dividing both would divide the Gaussian integer `α_S`.
For an odd split prime this would require both conjugate Gaussian factors,
which (5) excludes. No inert prime occurs, and two is excluded by odd `C`.
The coordinates have opposite parity since their squared sum is odd.
Neither coordinate is zero: otherwise primitivity would force norm one.
Their absolute values cannot be equal since `C` is odd. Thus `n>m>0`, with
opposite parity and gcd one, and (6) is a primitive PPT.

Conversely, every primitive positive PPT has unique parameters `n>m>0`,
`gcd(n,m)=1`, `n+m` odd, once the odd leg is put first. Factor `n+im` in
`Z[i]`. An inert prime would divide both coordinates; the ramified prime
above two is excluded by odd norm. At each split prime power `p_j^{e_j}`,
all `e_j` copies go to exactly one of `π_j` and its conjugate, or again a
rational prime would divide both coordinates. Unique factorization gives
`n+im=u α_S` for one unit `u∈{1,−1,i,−i}` and one address `S`.
This also proves that no other positive hypotenuse can support a primitive
PPT.

Multiplication by a Gaussian unit and complex conjugation account for all
sign changes and coordinate interchanges of `(x,y)`. They do not change
(6). Two allocations `α_S,α_T` that are associates have `S=T`, by uniqueness
of the exponents at the pairwise nonassociate Gaussian primes. If one is
an associate of the other's conjugate, then `T=S^c`. Conversely,
`α_(S^c)=conj(α_S)` literally in (5). Therefore exactly each complement
pair, and no additional addresses, gives the same primitive triple.
Interchanging the legs creates no further identification: one is odd and
one is even, so their parity ordering is uniquely recoverable. QED.

The unit/conjugation action on these Gaussian representations has no
exceptional stabilizer: nonzero coordinates exclude the axes, and unequal
absolute coordinates exclude the diagonals. There are `4·2^r` primitive Gaussian
representations and eight per positive Euclid pair, yielding `2^{r−1}`.
At `C=1` there are only four axes, so that counting argument does not apply;
there is no nondegenerate PPT. This is why the constant term must be removed
before halving in the Dirichlet series.

### An executable inverse

Given a PPT with the marked norm data:

1. Recover its odd leg `a` by parity, then take the positive integer roots
   `n=sqrt((C+a)/2)` and `m=sqrt((C−a)/2)`.
2. Put `w=n+im`. For each `j`, determine whether `π_j` or `conj(π_j)` divides
   `w`. Exactly one has valuation `e_j`; place `j` in `S` precisely in the
   first case. Division is exact integer arithmetic: multiply by the
   conjugate and divide both coordinates by the norm.
3. Return the path pair `{P_S,P_(S^c)}`. For a unique representative, choose
   the member with `1∈S` and decode it by (3).

The possible unit relating `w` and `α_S` does not affect these valuations.
Changing a chosen Gaussian prime to an associate does not affect the
returned triple. Conjugating one selected prime toggles its coordinate;
reordering the primes permutes coordinates and requires a fresh ordered path
decoding. Thus the fiber correspondence is natural with its marks, rather
than an unmarked orientation on the set of triples. In particular it does
not construct a tournament whose vertices are the primitive triples.

## 4. E3, E4, E5 worked out

Use `π_5=2+i`, `π_13=3+2i`, and `π_17=4+i`, in that prime order. The following
rows choose `1∈S`; each row represents the listed path and its complementary
partner. Triple legs in this table are parity ordered.

| C | S | α_S | P_S | P_(S^c) | (a,b,C) |
|---:|---|---|---|---|---|
| 5 | {1} | 2+i | (1,2,0) | (2,0,1) | (3,4,5) |
| 65 | {1} | 8−i | (1,3,0,2) | (2,3,0,1) | (63,16,65) |
| 65 | {1,2} | 4+7i | (1,2,3,0) | (3,0,1,2) | (33,56,65) |
| 1105 | {1} | 31−12i | (1,4,0,2,3) | (2,3,4,0,1) | (817,744,1105) |
| 1105 | {1,2} | 23+24i | (1,2,4,0,3) | (3,4,0,1,2) | (47,1104,1105) |
| 1105 | {1,3} | 33+4i | (1,3,4,0,2) | (2,4,0,1,3) | (1073,264,1105) |
| 1105 | {1,2,3} | 9+32i | (1,2,3,4,0) | (4,0,1,2,3) | (943,576,1105) |

Thus the total Hamiltonian-path counts are `3,5,9`, whereas the relevant
quotient counts are `(3−1)/2=1`, `(5−1)/2=2`, and `(9−1)/2=4`.
The increasing path is excluded; it is not assigned an extra primitive
triple. The same `E3` address quotient also models the single triple
`(7,24,25)` from `(2+i)^2=3+4i`. Prime powers do not create additional
Boolean coordinates.

## 5. Prime support is the elementary-symmetric grading

Let `\mathcal C` contain one and every positive integer all of whose prime
factors are `1 mod 4`. For `Re(s)>1` and any complex marker `z`, define

\[
F(z,s)=\sum_{C\in\mathcal C}z^{\omega(C)}C^{-s}
=\prod_{p\equiv1\ (4)}\left(1+\frac{z p^{-s}}{1-p^{-s}}\right).
                                                                    \tag{7}
\]

The series at `C=1` is one. For fixed `Re(s)>1`, the sequence
`a_p(s)=p^{-s}/(1−p^{-s})=1/(p^s−1)` is absolutely summable. The product is
locally uniformly convergent in `z` and `s` on compact subsets of this
domain. Expanding by the primes selected, each used prime contributes any
positive exponent but exactly one factor of `z`. Therefore

\[
[z^r]F(z,s)=e_r(\{a_p(s)\})
=\sum_{\substack{C\in\mathcal C\\\omega(C)=r}}C^{-s}.                 \tag{8}
\]

In particular the exact norm-height contribution of the `E_(r+2)` path
quotients is

\[
P_r(s)=\sum_{\substack{\text{primitive PPTs}\\\omega(C)=r}}C^{-s}
=2^{r-1}e_r(\{a_p(s)\}),\qquad r\ge1.                               \tag{9}
\]

The first three stages are `e1`, `2e2`, and `4e3`, respectively. Each stage
has infinitely many arithmetically marked copies of its finite tournament.
A bare path has no norm weight until its prime powers have been supplied.

The bound `|e_r|≤(Σ_p|a_p|)^r/r!` justifies summing over `r`, so

\[
P(s)=\sum_{r\ge1}P_r(s)=\frac{F(2,s)-1}{2}.                           \tag{10}
\]

At `z=2` the local factors simplify:

\[
F(2,s)=\prod_{p\equiv1\ (4)}\frac{1+p^{-s}}{1-p^{-s}}
=\frac{\zeta(s)\beta(s)}{(1+2^{-s})\zeta(2s)}.                       \tag{11}
\]

Here `β(s)=Σ_{j≥0}(−1)^j/(2j+1)^s` is Dirichlet beta.
For clarity, at a split prime `ζβ/ζ(2s)` has factor
`(1+p^−s)/(1−p^−s)`; at an inert prime its factor is one; at two it has
factor `1+2^−s`, which must be removed. This proves the normalization
without overlooking the ramified prime.

With `τ=x+iy` and the atlas convention

\[
E(\tau,s)=\frac12\sum_{\gcd(c,d)=1}\frac{y^s}{|c\tau+d|^{2s}},
\]

the nonzero Gaussian lattice sum is `4ζ(s)β(s)`. Separating the positive
gcd of `(c,d)` contributes `ζ(2s)`, so

\[
E(i,s)=\frac{2\zeta(s)\beta(s)}{\zeta(2s)}
       =2(1+2^{-s})F(2,s).                            \tag{12}
\]

One may derive the lattice sum by summing the Gaussian ideal Euler product:
four units give four generators per ideal; split primes give two factors,
inert primes have norm `p²`, and the ramified prime has norm two. Equations
(10)–(12) prove (2). The Gaussian ideal series and the Eisenstein
normalization are used only in their absolutely convergent half-plane.

The atlas's constant Fourier coefficient of `E(τ,s)` contains `ζ(2s−1)`;
(12) instead evaluates the full Eisenstein series at `τ=i` and contains `ζ(s)`.
Hence the atlas's scattering parameters `3,4,5` and this note's norm-series
parameters `5,7,9` realize the same three odd zeta values by different
observations. None sets `r=s−2`, nor ties prime-support size to an odd-zeta
weight. The odd-zeta note supplies the exact rational factors expressing
`P(5),P(7),P(9)` in terms of `ζ(5)/π^5,ζ(7)/π^7,ζ(9)/π^9`.

## 6. The exact square-sum continuation, and its angle boundary

Compose (6) with the signed AP-square map from the odd-zeta note:

\[
[P_S]\longmapsto
2(C^2-4ab,\ C^2,\ C^2+4ab).                           \tag{13}
\]

This is a bijection from the marked path quotient to the base signed AP
square-sum triangles with middle label `2C²`. Positive integer square scale
`g²` gives all scales, with pair-sum roots
`2g(|a−b|,C,a+b)`. This continuation preserves all three square-sum edges
and the AP relation. It does not claim a Hamiltonian path through a whole
square-sum graph on a consecutive label interval.

Positive labels require `C²>4ab`. In the raw Gaussian phase

\[
t_S=\arg\alpha_S\equiv\sum_j\varepsilon_j e_j\arg\pi_j\pmod{2\pi},
\qquad\varepsilon_j=1\ (j\in S),\ -1\ (j\notin S),
\]

this condition is exactly

\[
2|\sin(4t_S)|<1.                                      \tag{14}
\]

Indeed `a=|x²−y²|`, `b=2|xy|`, and `4ab/C²=2|sin(4t_S)|`.
The predicate is unchanged by units and conjugation, so it is well defined
on the path quotient **with its arithmetic marks**. It is not determined
by the number of paths, by `r`, or by the address alone without the prime
powers.

At `C=65`, the `{1}` orbit gives the positive labels
`(386,8450,16514)`, whereas `{1,2}` gives
`(−6334,8450,23234)`. At `C=1105`, precisely the `{1,2}` and `{1,3}`
rows pass the positivity filter. Even fixing the prime support and the
single `E3` orbit is insufficient: `C=5` gives `(3,4,5)` and fails, whereas
`C=5^5=3125` gives `(237,3116,3125)` and passes. This is an explicit
exponent-sensitive hostile to forgetting phase.

The unrestricted series (10) consequently includes negative smallest
labels. Retaining only positive AP configurations requires the angular
filter from the odd-zeta note; it cannot be achieved by simply multiplying
(10) by one third or discarding a fixed fraction of every Boolean fiber.

### 6.1 Every split-prime power has an exact one-third exponent law

**Theorem 2.** Fix any prime `p≡1 mod 4`. For each exponent `e≥1`, there is
exactly one primitive PPT of hypotenuse `p^e`. The set of exponents for which
its associated AP square-sum triangle has positive labels has natural density
`1/3`:

\[
\lim_{N\to\infty}\frac1N
\#\{1\le e\le N:\ p^{2e}>4a_e b_e\}=\frac13.                         \tag{15}
\]

In particular there are infinitely many accepted and infinitely many rejected
prime powers although **every fiber uses the same single E3 address orbit**.

**Proof.** Choose `π_p=a+ib` with `a>b>0` and write `t=arg π_p`. The number
`t/π` is irrational. If it were rational, `π_p/conj(π_p)` would be a root
of unity in `Q(i)`, whose only roots of unity are `1,−1,i,−i`. The first two
force an axis, and the latter two a diagonal; neither is possible for this
Gaussian prime of odd norm.

Consequently `2t/π` is irrational. For every nonzero integer `k`, the
geometric-sum identity gives

\[
\frac1N\sum_{e=1}^N\exp(4ik e t)\longrightarrow0,
\]

because its denominator `1−exp(4ik t)` is nonzero while its numerator is
bounded. Thus the fractional parts of `2et/π` are uniformly distributed
on `[0,1)`: the displayed limit proves the assertion first for trigonometric
polynomials, then for intervals by continuous/trigonometric approximation.
Equation (14) accepts precisely the intervals where
`|sin(2πu)|<1/2`, whose total length is `1/3`. No boundary exponent occurs,
as equality would give the odd integer `p^(2e)=4a_e b_e`. This proves (15).
Uniqueness of the primitive triple follows from Theorem 1 with `r=1`. QED.

This is a density with respect to the **exponent index at one fixed prime**.
The global hypotenuse-cutoff density in the AP note and its very different
Berggren-depth fractions use different samples. None of these density
statements asserts that a convergent weighted Dirichlet series is one third
of another, or supplies an error rate in (15).

### 6.2 A weighted path signature retains exactly the missing phase

For the marked interior order, attach invertible variables `z_1,...,z_r` and
form the Laurent polynomial

\[
W_{E_{r+2}}(z_1,\ldots,z_r)
=\sum_{P_S\ne I}\prod_{j\in S}z_j\prod_{j\notin S}z_j^{-1}
=\prod_{j=1}^r(z_j+z_j^{-1}).                          \tag{16}
\]

The last equality follows by independently choosing one of the two terms
from each factor; equation (3) identifies the resulting choice with the
actual ordered path. The exceptional path contributes no term here.

For an integer `k≥0`, specialize

\[
z_j=\exp(8ik e_j\arg\pi_j).
\]

The weight of `P_S` is then `exp(8ik arg α_S)`, and address complement
conjugates that weight. Let `θ=arg(a+ib)∈(0,π/2)` for the parity-ordered
legs of a primitive triple. Since `a+ib` comes from squaring a normalized
Gaussian parameter, `cos(4kθ)=cos(8k arg α_S)`. Normalization changes only
sign and Gaussian unit phase, which this frequency removes. Consequently

\[
\boxed{\displaystyle
\frac12 W_{E_{r+2}}\bigl(e^{8ik e_1\arg\pi_1},\ldots,
                        e^{8ik e_r\arg\pi_r}\bigr)
=\sum_{(a,b,C)\in\operatorname{PPT}(C)}\cos(4k\theta).}              \tag{17}
\]

Thus (17) is exactly the finite-norm angular coefficient in the series
`S_k(s)` of [the odd-zeta note, §7](ap_square_odd_zeta_20261019.md).
Summing it over `C` gives the local factors

\[
1+2S_k(s)=\prod_{p\equiv1\ (4)}
\left(1+2\sum_{e\ge1}\cos(8ke\arg\pi_p)\,p^{-es}\right),
\qquad\Re s>1.                                        \tag{18}
\]

The prior note's Fourier/Fejér angular filter can therefore be applied to
these weighted path sums to recover the positive-AP selection. At `k=0`,
all variables equal one and `W(1,...,1)=H(E_(r+2))−1=2^r`;
(17) becomes the unweighted fiber count. That scalar specialization is
precisely where the phase information is erased.

This is a **marked-address weight** with supplied prime powers and Gaussian
phases. It is not claimed to be an intrinsic invariant of an unmarked
tournament. It does not assign Gaussian phases to its vertices without
those arithmetic data.

## 7. Connection contract and reproducibility

The lossless input is a nonexceptional ordered path together with
`(p_j,e_j,π_j)` and its complement identification. Equations (5)–(6) recover
the entire primitive triple; the inverse recovers the path orbit. Forgetting
the complement phase loses only the redundant global conjugation. Forgetting
the prime powers loses the norm, sides, angle, and square-sum labels.
Forgetting the path but retaining `C` preserves the middle label and
multiplicity, while losing the individual triple. Summing at `s` preserves
a weighted count, not an individual configuration or tournament orientation.

The script supplies finite, independent controls:

- Every `C≤100000`, including all excluded norms, is checked against a separate
  direct Euclid-parameter census. There are **9622** eligible norms,
  **15919** primitive triples, and **31838** raw subset addresses. The PPT
  counts at support sizes `1,2,3,4` are `4825,8142,2856,96`.
- Every Gaussian allocation in that universe checks its norm, its exact
  complement partner, all eight unit/conjugation transforms, leg parity,
  the explicit inverse, and agreement under leg interchange.
- Literal permutation enumeration checks all Hamiltonian paths for `E3`
  through `E8`; independent subset dynamic programming checks `E3` through
  `E10`. The counts are `3,5,9,17,33,65,129,257`.
- Literal vertex-permutation enumeration checks the automorphisms through
  `E7` and rejects every proposed automorphism-induced complement action.
- Independent finite Euler multiplication checks every formal coefficient
  through norm **10000**, with the prime-support variable retained.
  Exact rational elementary symmetric sums and all exponent vectors in
  `{0,1,2,3}^4` for primes `5,13,17,29` separately test the two gradings.
- Primitive integer-coordinate enumeration through Gaussian norm **1000**
  checks (12), including the norm-one axis and norm-two terms.
- At every eligible `C≤1000` and additionally `C=1105,32045,1185665`,
  modes `k=0,1,2,3` compare the literal weighted address sum, the factored
  Laurent polynomial, and the PPT angular sum in exact rational arithmetic.
- Exact Gaussian powers for each of `p=5,13,17,29,41` through exponent
  **1000** test the phase predicate, primitivity, and absence of boundary
  equality. These finite samples do not prove the exponent-density limit.
- Explicit display fibers additionally include `C=32045` and `C=1185665`,
  with four and five distinct primes. Invalid norms, nonprimitive triples,
  malformed paths, Boolean/float addresses, and wrong-fiber inverses are
  rejected.

From the repository root, run:

```bash
python3 04-computation/experiments/endpoint_pythagorean_fibers_20261019.py
python3 -O 04-computation/experiments/endpoint_pythagorean_fibers_20261019.py
```

Before integration, run the same file from its staging directory. Both modes
produce identical JSON with **604,921 passing checks**, and no check uses a
removable Python assertion. The JSON SHA256 is
`4b154006c1ea78c4a137d7922b97e7b8470ee46c701eb2f11d8b74e72ba1f033`.
An independent supplementary checker enumerates integer legs and square-sum
roots directly through `C≤500`: **2,400 checks**, **80** primitive triples,
and **80** base signed AP triangles. Run
`python3 04-computation/experiments/endpoint_pythagorean_fibers_20261019_audit.py`
(and the same command with `-O`); their outputs also agree.
This closes the specified exact synthesis. It does not prove a new odd-zeta
irrationality statement, a tournament path-spectrum classification, or a
global square-sum Hamiltonicity theorem.
