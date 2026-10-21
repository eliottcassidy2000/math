# Deficit-six centers: a quadratic cutoff, complete order ten, and literal ports

Research date: October 21, 2026. **PROVED current bound:35**, by combining the
incoming odd-shore/score-run theorem with the inherited closing-prefix theorem;
see [the combined proof](COMBINED_CUTOFF.md). The independent quadratic
bound248 in section2 remains valid but is superseded as the sharpest cutoff.
**FINITE-EXACT:** all128 necessary normalized order-ten presentations are
excluded. The known center polynomial is therefore

\[
R_6(z)=z^5+z^6+z^7+O(z^{12}),\qquad \deg R_6\le35.
\]

There remain **ten possible orders**, namely
`12,14,15,18,20,21,24,25,28,35`. Their coefficients are **OPEN**. Neither
the cutoff nor the order-ten computation establishes that the three known
centers are the only centers. The original91-order reduction and receipt below
are retained as valid weaker evidence; the combined receipt records the
stronger support. The safe recurrence onset is now48, improving section6's261.

## 1. Inheritance and precise objects

Let $T$ be a strong tournament, $A=\operatorname{Aut}(T)$, and $H$ its number of
directed Hamilton words. Write $h=H/|A|$. A word $P$ is blue when its positional
reversal $r_P$ is an anti-automorphism; let $b$ count blue $A$-orbits and put
$\Delta=h-b$. In a self-converse source, converse followed by word reversal is
an involution on path orbits whose fixed points are exactly the blue ones.
Consequently $\Delta$ is even. “Self-converse” is not an abbreviation for strong.

The retained inheritance is:

- **PROVED:** [the effective cubic bound](../effective_black_defect_order_bound_20261020/PROOF.txt)
  gives $q\le\Delta+2$ and the two connectivity bounds used below. Its general
  scope includes non-self-converse tournaments. The present strengthening uses
  self-converseness essentially.
- **PROVED, audited:** the [deficit-two packet](../deficit_two_centers_20261020/README.md)
  promotes its preserved manuscripts to audited status. In particular,
  [the prefix proof, §1](../deficit_two_centers_20261020/session15-deficit-centers/PREFIX_CONNECTIVITY_CUTOFF.md)
  supplies the stronger black-mate count $\Delta\ge2(t-1)$.
- **PROVED, audited:** [cycle-cut obstruction, §5](../paired_black_leaf_reduction_20261020/CYCLE_CUT_NEAR_PURITY_OBSTRUCTION.md)
  supplies a nonautomorphic Hamilton-cycle successor at every odd order at
  least seven. This uses its explicitly closing circulant reroutings; the
  deficit-two-only conclusions elsewhere in that note are not applied to six.
- **FINITE-EXACT inherited universes:** the complete ordinary class tables
  through eight and the complete order-nine table are named and hashed in
  [receipt.json](receipt.json). This checker filters them again and independently
  replays each of the three surviving sources. It does not newly enumerate all
  order-nine tournaments.
- [The deficiency-three center law](../deficiency_three_20261020/centers.md)
  gives the SCC factorization. The
  [session21 reconstruction](../session21_recovery_20261021/RECONSTRUCTED_METAGRAPH_PROOFS.md)
  explicitly has **RECONSTRUCTED** status: it records formerly checked formulas
  and missing certificate bytes, rather than being a fresh census. Its
  trivalent algebra is rederived below; its old 4095 cutoff is superseded here.

The classical input is **CITED**: Moon, *On Subtournaments of a Tournament*,
Theorem 1, p.297, gives a directed cycle of every length from three to the
order of a strong tournament. The primary theorem statement was checked in
[the published PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/5D7CCBE66AB07BB3F90F771D1CD3CD08/S0008439500052875a.pdf/on-subtournaments-of-a-tournament.pdf).
The repo's [Moon sidecar](../../reference/CORE-PAPERS-TOURNAMENT-STRONG-EARS-2026-08-25.md)
records the attribution correction and its limited role.

The live concepts are regular rotation orbits, exact endpoint attachment,
reflection-paired black orbits, and the unknown finite center alphabet. The
new operation averages inside one regular orbit, rather than bounding every
vertex by the global minimum degree. The canonical hostile is pure-blue
carousel $C_5$: two of its Hamilton words can differ by an even-order
permutation, so the reflected-support condition cannot be dropped. The least
used sidecar is the common suffix and its first vertex, which makes every
selected cycle cut an actual Hamilton word of the same source.

## 2. Orbit-local averaging gives an $O(\Delta^2)$ order bound

**Theorem.** Let $T$ be strong and self-converse, of order $n\ge3$, and assume
$\Delta(T)\le d$, where $d\ge0$ is even. Let $r$ be the least power of two at
least $d+4$. Then

\[
n\le\max\{3r-1,(d+2)(2r-1)\}<4(d+2)(d+4).                 \tag{1}
\]

In particular
$\Delta(T)>\sqrt{n/4+1}-3$. This inverse inequality is a lower bound, not a
classification or an equality assertion.

**The inherited frame.** Choose a Hamilton cycle, with successor permutation
$R$. Put $m=|A\cap\langle R\rangle|$ and $q=n/m$. Tournament automorphism groups
have odd order, so $m$ is odd. The subgroup is generated by $R^q$; its $q$
vertex orbits each have size $m$ and induce regular tournaments. The cycle
gives a directed perfect matching from each orbit to the next.

The cycle cuts occupy exactly $q$ path orbits, because their unique
coordinatewise comparison maps are powers of $R$. The number of blue cut
orbits is zero or $\gcd(2,q)$. Thus $q\le d+2$. The audited connectivity
estimates are

\[
\kappa(T)\ge\lceil m/3\rceil\quad(q=1),\qquad
\kappa(T)\ge\lceil m/2\rceil\quad(q\ge2).                 \tag{2}
\]

The second estimate uses the cyclic perfect matchings as well as regularity;
it is not asserted for an arbitrary regular tournament or for $q=1$.

**The new average.** Suppose first that $q\ge2$ and $m\ge2r+1$. In one regular
orbit $O$, choose a directed $r$-cycle with vertex set $S$. This is possible
because a nontrivial regular tournament is strong and Moon applies. Its
induced chords need not have any special pattern. Regularity gives exactly

\[
\#(S\longrightarrow O\setminus S)
 =r\frac{m-1}{2}-\binom r2=\frac{r(m-r)}2.                 \tag{3}
\]

Therefore some $v\in O\setminus S$ is beaten by at least $r/2$ vertices of
$S$. Equation (2) gives $\kappa(T)\ge r+1$, so $T-S$ is strong. Cut a Hamilton
cycle there to obtain a suffix $W$ beginning at this literal vertex $v$.
For each selected exit $x\to v$, cut the $S$-cycle immediately after $x$ and
append the same $W$.

These $t\ge r/2$ actual words have pairwise comparison permutations equal to
nontrivial powers of an $r$-cycle. Their orders are even, since $r$ is a power
of two. Their supports lie in the first $r$ positions, disjoint from the last
$r$ positions because $n\ge2r$.

For completeness, the inherited factor two has three obligations. Two such
words cannot be in one $A$-orbit: their unique comparison would be an
even-order automorphism. They cannot both be blue: composing their positional
reversals yields two disjoint conjugate factors of even order, hence an
even-order automorphism. Their black orbits cannot be converse mates either:
if the comparison is $u$, then $u r_P$ would be anti, and
$(u r_P)^2=u(r_Pur_P)$ gives the same contradiction. At most one selected word
is blue, and each other word requires a distinct black orbit and a distinct
unselected mate. Hence

\[
\Delta\ge2(t-1)\ge r-2>d.                                \tag{4}
\]

This contradiction proves $m\le2r-1$ for $q\ge2$. For $q=1$ the identical
argument works whenever $m\ge3r+1$, using the first bound in (2). Since $m$ is
odd and $3r$ is even, $m\le3r-1$. Multiplication by $q\le d+2$ proves (1).
Finally $r<2(d+4)$ gives the displayed quadratic bound. QED.

For $d=6$, $r=16$, so $q=1$ gives $n\le47$, and $q\ge2$ gives $m\le31$ and
$n\le8\cdot31=248$. An eight-cycle supplies only the average $t=4$ and
$\Delta\ge6$: this equality is not a contradiction to deficit six. The
checker retains that sharp logical boundary.

## 3. Honest residual orders and source-bound prefix receipts

This section records the original quadratic reduction. Its91-order support
is superseded by the ten-order intersection in [the combined proof](COMBINED_CUTOFF.md).

For even $n$, $q\ge2$ automatically because $m$ is odd. For odd $n\ge7$, the
inherited closing-prefix theorem ensures that some Hamilton-cycle successor
is not an automorphism, so one may also choose $q\ge2$. Applying §2 to that
cycle gives the necessary support

\[
n=qm,\qquad 2\le q\le8,\quad m\text{ odd},\quad1\le m\le31.              \tag{5}
\]

In particular $16\nmid n$, and order eleven is impossible. Before the new
order-ten computation there are 92 candidate orders above nine. Removing ten
leaves the 91 orders listed in the receipt. Equation (5) is a necessary
condition only. It does not assert an admissible tournament at every listed
order, or that every Hamilton cycle has a blue cut.

Two positive literal receipts illustrate both connectivity branches. The
first source is carousel $C_{49}$, with the cycle
$(0,3,6,\ldots,45)$ and suffix beginning at 1. The second has two interleaved
carousel-33 shores, cross rule $a_i\to b_j$ iff $i=j$, the cycle
$(0,4,8,\ldots,60)$ and suffix beginning at 2. Both retain eight exits and
prove $\Delta\ge14$. Their full bit rows, anti-isomorphism, prefix cycle,
suffix, eight Hamilton words, and eight complement targets are in the receipt.

`audit_prefix(source, receipt)` does not search for any word. It checks the
source identity, every literal join, the full Hamilton permutations, even
relative orders, disjoint reflected supports, and the targets obtained by
reversing precisely the non-path arcs. It rejects duplicated cycles, aliases,
repeated exits, substituted sources, wrong suffixes and wrong targets.
The bound is on the source's black orbit count, not on connectivity of its
metagraph or the identity of the neighboring target classes.

## 4. Complete normalized order-ten exclusion

Assume a deficit-six center has order ten. The cycle-cut bound says $q\le8$;
$m$ is an odd divisor of ten. Therefore $m=5,q=2$. Label one Hamilton cycle
$0,1,\ldots,9$. Rotation by two is an actual automorphism. Put
$a_i=2i,b_i=2i+1$ for $i\in\mathbb Z/5$.

Each shore is a labelled circulant tournament. Its connection set chooses one
element of each opposite pair, giving four choices per shore. Cross edges have
the form $a_i\to b_j$ iff $f(j-i)=1$. The marked Hamilton cycle forces
$f(0)=1$ and $f(-1)=0$, leaving three free bits. Thus all possible sources are
represented among exactly

\[
4\cdot4\cdot2^3=128                                      \tag{6}
\]

normalized frames. This is a complete necessary universe, not 128 distinct
isomorphism classes. There is no self-converse, reflection or score filter
silently imposed in its generation. Every generated source is strong because
its normalizing Hamilton cycle is checked literally.

For each source the checker computes $H$ by subset DP. It enumerates **all**
score-compatible permutations and tests every arc to obtain $\operatorname{Aut}$
and all anti-isomorphisms; it does not use only the rotation centralizer.
Every even-order-source anti-involution pairs the vertices without fixed
points. A blue word is uniquely obtained by ordering these pairs, selecting
one member of each for the first half, and reflecting the reversed half.
The first-half arcs and the middle arc suffice; the second-half arcs follow
from the checked anti-isomorphism. Every blue word has exactly one positional
reversal, so summing over anti-involutions counts it exactly once. Dividing
$H$ and the literal blue-word count by $|A|$ gives $h,b$.

Exactly 64 frames are self-converse. Their complete defect histogram is:

| $\Delta$ | Normalized frames |
|---:|---:|
| 1528 | 8 |
| 1628 | 8 |
| 2908 | 12 |
| 2926 | 12 |
| 3064 | 12 |
| 3114 | 12 |

None has deficit six. Each row retains its actual source, normalizing cycle,
automorphism list, anti-involutions, a blue word when one exists, a black word,
and that word's complement target. The separate eight-frame shore-size-three
control recovers the known order-six center in four frames, with
$(H,|A|,h,b)=(45,3,15,9)$; the other four have deficit twelve. Direct permutation
enumeration independently verifies its DP and blue-halfword counts.

## 5. The inherited small centers are retained as actual ports

Filtering the complete inherited tables through nine gives exactly these
three strong self-converse deficit-six sources:

| Order | Literal out-neighbor bit rows | $(H,|A|,h,b)$ |
|---:|---|---|
| 5 | $(4,1,18,7,11)$ | $(9,1,9,3)$ |
| 6 | $(36,17,10,19,37,14)$ | $(45,3,15,9)$ |
| 7 | $(112,73,35,21,70,26,44)$ | $(189,21,9,3)$ |

The six-vertex source has homogeneous pairs $(0,4),(1,3),(2,5)$ and quotient
$C_3$: it is the already present center $C_3[K_2]$, not a newly discovered
family. The seven-vertex source is Paley7; the receipt retains the exact
isomorphism to the connection set $\{1,2,4\}$.

For these three sources the present checker enumerates every permutation,
every Hamilton word, every automorphism and the complete path-orbit relation.
Two words belong to the same orbit exactly when they induce identical ordered
source rows; their unique coordinatewise map is then an automorphism. Each
of the six black orbits is retained with a literal word, its marked endpoints,
its converse mate and its actual complement target. Thus three mate pairs are
authenticated for each center. No equality of counts is used as a substitute
for port identity or target adjacency.

## 6. Enumeration payoff: order eighteen, onset at most 261

The denominator/minimality proof below remains current. The subsequent
combined cutoff35 improves its safe onset from261 to48; the original
degree248 synthetic control remains a test of the earlier uniform bound.

Let $p(z)=1-z-z^3$, $u(z)=p(z^2)$ and $A(z)=1/p(z)$. The inherited SCC laws,
using the now proved three-presentation alphabet, give the trivalent
paired-black vertex series

\[
G(z)=\frac{z^5}{2}\left(A(z)^2-A(z^2)\right)
       +A(z^2)R_6(z)+z^{10}A(z^2)^2(1+z+z^3).             \tag{7}
\]

This is the reconstructed formula's algebraic consequence, rederived here;
the unrecovered old checker is not described as having been rerun. The common
denominator is $p^2u^2$, of degree eighteen, and the numerator is

\[
N=\tfrac12z^5(u^2-p^2u)+p^2uR_6+z^{10}p^2(1+z+z^3).      \tag{8}
\]

The new cutoff gives $\deg N\le260$, so the resulting homogeneous recurrence
holds for every $n\ge261$. This improves the recovered safe onset 4108.

The denominator is exactly reduced independently of the unknown coefficients.
The polynomials $p,u$ are squarefree and coprime: at a root of $p$ one has
$u=2z^4\ne0$. There (8) reduces to $z^5u^2/2\ne0$. At a root of $u$, (8)
reduces to $z^{10}p^2(1+z+z^3)$; the last factor cannot vanish jointly with
$u$, since it would again force $2z^4=0$. Hence the minimal eventual
constant-coefficient recurrence order is eighteen. This is not a computed
list of the missing initial coefficients.

Exact polynomial controls include the synthetic choice
$R_6=z^5+z^6+z^7+z^{248}$. Its numerator has degree 260 and last coefficient
$-1$, so the homogeneous identity fails at 260 and holds from 261. This tests
the uniform algebraic onset; it does not assert a center of order 248.

## 7. Verification, limits and incoming connections

### Incoming fixed-degree rationality: a uniform center cutoff

Incoming commit `777c7240f7` proves
[rationality at every fixed fundamental degree](../fixed_degree_metagraph_rationality_20261021/README.md)
from the inherited finite-atom and center-divisor laws. It is a counting
theorem, not a wiring classification. Section2 above strengthens its finite
center input without changing its denominator bounds.

For degree k>=1, every strong self-converse center in its deficit-2k
decomposition has signature `(a,b)` and some odd half-count x>=1 satisfying
`x(x*a-b)=2k`. Thus

\[
0\le a-b=2k/x-(x-1)a\le2k.
\]

If r is the least power of two at least2k+4, **every** such strong center has
order at most `max(3r-1,(2k+2)(2r-1))<4(2k+2)(2k+4)`. This bounds the degree
of each of the incoming polynomials `R_(2k)` and `E_(x,2k)`; the empty center
has order zero and is kept separately. The quadratic bound also covers
zero-deficit strong centers, since the theorem assumes an upper bound on
deficit, not equality. It does not apply to ordinary non-self-converse
strong atoms `P_j`, and supplies none of the missing center coefficients.
The parent and an independent reader checked this direct substitution.

The later incoming [odd-shore proof](../odd_shore_deficit_centers_20261021/README.md)
sharpens both finite alphabets further. Its separately reviewed non-self-
converse argument gives `n<4j^2` for a strong atom with odd `h=j>=3`.
On the self-converse side, odd h forces b>=1 by the path-orbit involution;
substituting d=j-1 in that packet's section6 gives `n<2j(j+1)<=4j^2`.
Thus **every** strong atom with odd h=j>=3 has order `<4j^2`, by these
incoming mechanisms. The singleton and C3 at h=1 remain separate. This
stronger inherited corollary is credited to that packet, not to the earlier
averaging argument in section2. It does not determine the atom polynomials.

### Reproduction

Run from the repository root:

```text
python -B 05-knowledge/results/deficit_six_centers_20261021b/check.py --verify
python -B -O 05-knowledge/results/deficit_six_centers_20261021b/check.py --verify
python -B -O 05-knowledge/results/deficit_six_centers_20261021b/combined_check.py --verify
```

The saved receipt is deterministic and uses exact integers and rational
polynomial arithmetic. All checks remain active under `-O`. The declared
universe includes the 128 complete necessary order-ten frames, eight smaller
positive frames, three fully replayed center sources, the inherited census
filters, regular-orbit subset averages, literal prefix certificates and their
hostiles. It contains no search through the original91-order remainder or
the current ten-order remainder. The combined checker separately verifies
the necessary support and formal recurrence consequence.

The precise next task is to exclude or classify the normalized cyclic frames
at the ten remaining orders listed above, starting with twelve. A successful exclusion must retain
an actual source and a sufficiently large half-supported path clique, or
compute the full orbit deficit. Merely finding a removable cycle, a large
path count, or a small automorphism subgroup does not settle the coefficient.
There is no assertion here about active black-core connectivity or bridges.

Normal and optimized immutable verification both reproduce **7,041 checks**.
The LF receipt SHA-256 is
`8ea09c7747176477426a89184df29e5d385fe75e9fdfe11d8ccc70f6586d39af`.

**Independent proof/API audit: PASS.** The parent checked the complete cutoff,
normalization and generating-function proof. A separate peer checked the
source-bound prefix consumer, regular-orbit averaging, both connectivity
thresholds, distinct black mates, complete order-ten normalization and blue
halfword counting, and independently ran `-B -O --verify` against the saved
receipt. It matched exactly, with no requested correction. This audit covered
the original91-order reduction; the separately audited combined proof reduces
the unresolved support to ten orders.
