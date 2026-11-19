# Signed cancellation with two independently growing blocks

**Proved:** the six-arc patch for defect six extends to a two-parameter
family with unbounded canonical defect. Every admissible target with
canonical defect below 14 is now realized. At defect 14, the only
remaining family in this argument is $38\cdot2^j+1$; its first four
heights also have exact witnesses. This does not close Chan–Pak.

## One physical patch, two arbitrary widths

Use the rooted binary compiler and numerical defect $D=t-H$ from
[LOW_DEFECT_DESCENT](LOW_DEFECT_DESCENT.md). The canonical word
$0\,1^s\,0\,1^{m+4}$, with $s\ge0$ and $m\ge1$, has
$$
N=(12\cdot2^s+2)2^{m+3}+1,\qquad
D=\begin{cases}6 & s=0,\\ 2+3\cdot2^s & s\ge1.\end{cases}
$$
These follow by substituting in the old rooted recurrence: just before
the second zero its marked reading is $3\cdot2^s+1$. When $s\ge1$
the first zero contributes defect 2, and the second contributes
$3\cdot2^s$. For $s=0$ the two initial zeros are one run.

The quotient at $s=m=1$ is the word 01011111, with vertex labels in
birth order beginning at 0. Mark the first near-sink (vertex 3), and
the last near-sink (vertex 9). Reverse arcs on
$$
(0,1),(0,2),(0,6),(1,8),(4,5),(4,8).
$$
The output masks are 80, 765, 1017, 1009, 704, 977, 896, 769, 531,
1. Replace vertices 3 and 9 independently by transitive blocks of
widths $s$ and $m$, respectively (delete 3 for $s=0$). The unpatched
version is literally the canonical word above, since the two blocks
have the original uniform outside contacts. Neither marked vertex
is touched. The same reversal on the named quotient arcs is therefore
a physical patch at every width, without introducing new vertices.

For $s\ge1$, denote the compressed visit coefficients by $c_{r,q}$.
The inner marked block can have at most three runs: each run except
the first must be entered from one of its only two incoming exterior
vertices (1 and 2). The outer mark beats only 0 outside its own block,
so it can have at most two runs. Both statements hold before and
after the patch. Each compressed path word gives exactly
$r!{s\brace r}\,q!{m\brace q}$ paths by independently distributing
vertices among nonempty increasing runs.

Exact DP gives the following differences of patched and unpatched
counts; the columns are $m=1,2$.

| $s$ | $m=1$ | $m=2$ |
|---:|---:|---:|
| 1 | 8 | 8 |
| 2 | 14 | 14 |
| 3 | 26 | 26 |

The finite grid is a triangular change of basis, not a fit: the
matrices of $r!{s\brace r}$ and $q!{m\brace q}$ on these indices
have nonzero diagonals. It certifies exactly
$\Delta c_{1,1}=8$, $\Delta c_{2,1}=3$, all other differences
zero. Thus $\Delta H=8+3(2^s-2)=2+3\cdot2^s$ for **every**
$s,m\ge1$. In particular the cross terms involving two outer runs
cancel identically, while the inner mode intentionally changes.
For $s=0$ the quotient reduces to the
[defect-six certificate](MACNEILLE_DEFECT_REWRITE.md);
the two outer-width readings also verify it directly. This proves
exact placement at the original absolute target.

The cancellation is not linear in edge edits. At $s=m=1$, exact
Boolean Möbius expansion of the change in the two-outer-run
coefficient by interaction degree gives sums
$815,-976,153,8,0,0$, totaling zero. Only the empty and full edit
sets preserve that coefficient within these six edges. Looking
only at individual edits would predict a large drift of the
opposite kind. The replay computes all 64 evaluations.

## Low-defect descent

Write $j=r-1$ for a terminal run of $r\ge1$ ones. The positive
additions to defect are $2^{a+1}-2$ for an initial zero run $0^a$,
and $u(2^b-1)-1$ for a later zero run $0^b$ starting at reading $u$.
Up to defect 14 the exhaustive possibilities beyond defect six are:

| Defect | Word before the final ones | Target |
|---:|---|---|
| 8 | 010, 100, 1110 | $26\cdot2^j+1$, $22\cdot2^j+1$, $34\cdot2^j+1$ |
| 10 | none | |
| 12 | 1010 | $42\cdot2^j+1$ |
| 14 | 000, 0110, 1100 | $30\cdot2^j+1$, $50\cdot2^j+1$, $38\cdot2^j+1$ |

To see completeness, an initial zero run of length at least four
already costs 30. If its cost is 14, later runs are too expensive;
if it is 6, any later run adds at least 14. If its cost is 2, a later
single zero costs $3\cdot2^s$ after $s$ ones; only $s=1,2$ fit.
A longer later run or a third run exceeds the budget. With no
initial zeros, if the first zero run follows $s\ge1$ ones,
its cost is $2^s$ for length one, $3\cdot2^s+2$ for length two,
and at least 20 for length three. Two runs fit only when both are
singletons, with one initial and one intervening one: the cheapest
such costs are 2 and 10, and the next possibility costs at least
22. A third run then exceeds the bound.

The two-parameter patch supplies the 26 and 50 rays at $j\ge4$.
For other rays use the following marked tournaments, inflating their
mark by a transitive block of any positive width. The mark has
in- or outdegree one; this bounds visits by two in every inflation.
The columns $c_1,c_2$ are independently recovered from widths 1,2.

| Slope in the target | Starting $j$ | Outgoing masks | Mark | $(c_1,c_2)$ |
|---:|---:|---|---:|---|
| 22 | 0 | 46,60,24,48,33,4 | 5 | (23,11) |
| 34 | 1 | 126,252,248,48,96,192,8,89 | 6 | (69,34) |
| 42 | 3 | 8,9,19,52,67,343,399,319,31 | 7 | (337,168) |
| 30 | 2 | 126,124,248,48,224,192,8,75 | 6 | (121,60) |

The count is $c_1+c_2(2^w-2)=1+c_2\,2^w$ at width $w$.
Explicit small arc patches in the replay cover all earlier heights
on these six rays. The seventh ray, $38\cdot2^j+1$, has checked
patches for $j=0,1,2,3$ only in this packet. Even $j$ there also cannot
give the *least* counterexample: its target is divisible by 3 and is
larger than 63 for even $j\ge2$, so division by 3 cannot land at either
hole. Ordinal composition restores the factor 3.

## What the other 2–3 problems actually lend us

- **Disjointness slices.** In the [Diamond–Yehudayoff question](https://eccc.weizmann.ac.il/report/2021/148/revision/1/download/)
  the target coefficient is unchanged by reweighting forbidden
  words; giving 11 weight $-1/2$ makes cancellations away from the
  target and yields a uniform bound. On the tournament side, the
  roles differ: we are subtracting two *actual* counts, not giving
  paths negative multiplicity. The finite signed visit tensor tells
  us which growing terms vanish in their difference. The shared
  principle is to prove exact cancellation before passing from a
  local check to a uniform statement. An upper bound for one
  hyperplane does not imply a nonempty exact tournament fiber.
- **Catalan's equation.** If a proposed $+2$ mechanism is the
  difference of just two pure steps $2u,2v$ with $u,v$ 3-smooth
  (no prime factors outside 2 and 3), then $|u-v|=1$.
  Coprimality forces separate prime supports unless one term is 1.
  The only unordered possibilities are $\{1,2\}$, $\{2,3\}$,
  $\{3,4\}$, $\{8,9\}$, by the elementary 2–3 instance of
  [Mihăilescu's theorem](https://doi.org/10.1515/crll.2004.048)
  (the exponent-one cases are immediate). Large fresh adjacent
  pure steps cannot drive a universal construction. Mixed signed
  interactions, as in the patch above, are a different mechanism.
- **Collatz.** A block of accelerated odd steps composes to
  $x\mapsto(3^a x+B)/2^b$. The carry $B$ is fixed by the itinerary,
  not by the multiplier alone. This is the analogue of locking a
  growth coefficient without yet selecting the absolute count.
  [Tao's logarithmic-density result](https://arxiv.org/abs/1909.03562)
  makes the distinction between almost all and every starting
  point explicit; the same distinction matters for Chan–Pak.
- **Furstenberg.** Independent ×2 and ×3 motions in the
  [circle problem](https://math.ucsd.edu/~asalehig/F_Disjointness.pdf)
  yield a rigidity phenomenon about invariant sets. Its useful
  lesson here is to distinguish recurrence under growth from
  rigidity of the joint continuation profile. Orbit density in a
  completion is not a certificate of an exact integer witness.

A productive next generalization is to let the target choose a
multi-block quotient so that the *difference tensor* has no
unwanted mixed terms, rather than asking every single edit to
preserve every coefficient. The new family proves that this can
work with independently growing blocks and unbounded defect.
It does not prove that an arbitrary target has such a quotient.

Run the self-contained certificate with python3 signed_defect_replay.py.
