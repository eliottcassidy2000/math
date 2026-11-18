# Factorial coordinates and exact floors for tournament counts

**Proved below:** a single explicitly chosen strong self-converse tournament
count recovers the entire slot profile of an arbitrary source; exact finite
differences for transitive twin inflation; a lower bound on the core complexity
required to cover almost all odd targets. **Open:** selecting a physical core
and a compatible port or slot profile for every prescribed absolute count.

## Prior work used

- [Arithmetic braids II](../arithmetic_braids2_20260917_floor_reciprocity.md)
  uses floor reciprocity only after keeping the exceptional residue fibre.
  Here the discarded tail is bounded strictly before rounding, and the
  invisible structural coordinates are identified before counting.
- [Named boundary signatures](../session25_spectrum_wiring_20261021/session25-ports-review/AUDIT.md)
  uses falling factorials to label equality classes injectively. Here the
  equally important distinction is between ordered occupied runs and all
  available slots.
- [Weighted carrier](../tournament_future_boundary_20261020/README.md)
  supplies physical transitive twin modules in the metagraph. The floor below
  applies to any tournament quotient, not just a carrier.
- [Marked blue probes](../chanpak_blue_marked_probes_20261118/README.md)
  recovers the unmarked profile from several self-converse counts and
  distinguishes equal unmarked profiles by a port. The new single count
  recovers only the unmarked profile and does not remove the need for marks.

## 1. Run coefficients and a single blue reading

Write Q_j(X) for the number of spanning ordered covers by j **nonempty**
directed paths, and set it to zero outside 1..n for nonempty X of order n.
Then P_X(t)=sum_j Q_j binom(t,j) allows empty slots. With falling factorials
this is Q_j (t)_j/j!. The [Grinberg–Stanley odd-cycle packing identity](https://arxiv.org/abs/2307.05569)
gives P_X(-t)=(-1)^n P_X(t). In rising coordinates this gives the exact
interface constraint, for each 1<=j<=n,

    (-1)^n Q_j = sum_(k>=j) (-1)^k binom(k-1,j-1) Q_k.

Indeed (-t)_k=(-1)^k t^(overline k), and Vandermonde expands
binom(t+k-1,k) in binom(t,j). These conditions are necessary, not a test for
physical realizability.

Let Y_m=C3[X,X^op,TT_m] with cyclic block arcs. For m>=1 it is strong;
swapping the first two blocks by converse and reversing the transitive
chain is an anti-involution. It has a reversal-fixed Hamiltonian order:
choose a Hamilton path in the second block, then the increasing chain,
then the reflected reversed path in the first block. Thus this is an
actual blue frame. Put q_j=Q_j(X)=Q_j(X^op).

Compress a Hamilton path to its maximal block runs. The skeleton is forced
to follow the block cycle. If it uses j visits to the transitive block,
the three possible starting blocks and three possible ending blocks give

    H(Y_m) = sum_(j=1..n+1) A_j j! S(m,j),
    A_j = q_(j-1)^2 + 2q_(j-1)q_j + 3q_j^2 + 2q_j q_(j+1) + q_(j+1)^2.

In the skeleton, visits to the other two blocks are either both j (three
patterns), both j-1 or both j+1 (one each), or consecutive values adjacent
to j (two each). A pattern omitting a block contributes zero, as encoded
by q_0=0. Each transitive run has just one internal order, so its factor
counts onto maps from m vertices to j ordered runs.

By inclusion-exclusion put H(Y_m)=sum_(i=1..n+1) b_i i^m, where

    b_i = sum_(j>=i) (-1)^(j-i) binom(j,i) A_j.

For every source Q_j <= n! binom(n-1,j-1) <= K=n! 2^(n-1), by writing a
permutation and inserting j-1 cuts. Hence sum |b_i| <= M=9(n+1) K^2 2^(n+1).
Choose m=(n+1)(1+ceil(log_2(2M))). For every i>=2, the portion below i
divided by i^m has absolute value at most M((i-1)/i)^m < 1/2. This follows
from (n/(n+1))^(n+1)<1/2. Use integer bit lengths for the ceiling.

Given the one integer H(Y_m), peel from i=n+1 down to 1: divide the
remaining signed integer by i^m, add 1/2, take the floor, and subtract the
term found. The strict bound guarantees the correct signed digit without
a tie. Then invert the triangular binomial transform (diagonal 1) to get
A. Starting q_(n+1)=q_(n+2)=0, recover all q using

    q_(j-1) = isqrt(A_j - 2q_j^2 - 2q_j q_(j+1) - q_(j+1)^2) - q_j.

The radicand is exactly (q_(j-1)+q_j)^2. The height is O(n^2 log n).
This is storage in one witness, not a faster algorithm for evaluating it.
For a candidate target and a fixed n, it provides exact rejection tests:
bounded signed digits, integral nonnegative profile, q_n=n!, the
rising/falling reciprocity and a zero boundary q_0. Passing does not
supply an X with that profile.

A complementary low-digit fact: modulo p^a only visits j with
sum_(ell>=1) floor(j/p^ell)<a can contribute to the Stirling expansion,
because the other j! are divisible by p^a. High digits round away lower
bases; low prime-power digits kill high visit multiplicities. Neither
is a license to choose physical profile coordinates independently.

## 2. Physical forward differences and an exact twin floor

Let D be any d-vertex quotient. Replace each i by a nonempty transitive
chain of width m_i and orient all cross arcs according to D. More generally
than the cyclic formula, for each spanning directed quotient word w with
no consecutive repeated labels, let v_i(w) be its visit multiplicity.
Then the Hamilton count is

    sum_w product_i v_i(w)! S(m_i, v_i(w)).

There are finitely many nonzero terms at any width. This is a physical
statement, not a signed interpolation. If the quotient is self-converse
via a given anti-involution and opposite widths match, the blowup is
self-converse as well.

An equivalent expansion uses falling factorials in the widths themselves:
by the odd-cycle packing formula, group packings by the number k_i of
vertices used by cycles inside each chain, with singleton vertices left
over. Each choice of k_i vertices induces exactly the same ordered chain.
Thus H=sum_k c_k product_i binom(m_i,k_i) for nonnegative integers c_k
depending only on D and k, with c_0=1. Equivalently a kth forward
difference in m_i counts the weighted packings in which all k freshly
added marked twins occur in cycles: inclusion-exclusion cancels every
packing missing a mark. Both descriptions apply to same-chain or mixed
forward differences at arbitrary widths.

There is a clean physical dichotomy. If i lies on no directed triangle of
D, its in-neighbors all beat its out-neighbors. Then i is a singleton
strong component, its chain is a singleton-component interval in the
blowup, and its width is invisible to H. Otherwise choose a triangle
through i and one representative of each of its other vertices. Any odd
nonempty subset of the i-chain closes exactly one directed odd cycle on
those vertices and the two representatives. Keeping everything else
singleton proves

    H >= 1 + 2^(m_i),   so   m_i <= floor(log_2(H-1)).
    Delta_i^k H >= 2^(m_i)  for every k>=1 and m_i>=1.

For the latter use the same cycles containing the k marked twins and any
subset of the old chain of appropriate parity; there are exactly 2^(m_i-1)
such subsets, each cycle weighted 2. The former is sharp for a triangle
with the other widths 1, which gives H=1+2^m. Even a one-coordinate climb
can be far too steep to keep a +2 switch: beyond width one, no two distinct
widths along that climb differ by 2 in H. Width can be chosen from the
target by binary search and the floor bound, but hitting the target still
requires changing other physical data.

## 3. How big must the compressible core get?

In any quotient representation with H<=X, invisible widths can be
replaced by 1, and each active width lies in
1..L=max(1,floor(log_2(X-1))) for integer X>=2.
So the number of distinct counts obtainable from *all* quotients of
order at most d, not merely a fixed bank, is at most

    sum_(s=1..d) 2^binom(s,2) L^s.

In particular a fixed bank is polylogarithmically sparse. If
d <= (sqrt(2)-epsilon) sqrt(log_2 X), with 0<epsilon<sqrt(2),
the displayed bound is X^(1-delta+o(1)) for some delta>0. Thus almost
all odd integers up to X are outside the reach of every transitive blowup
with so small a quotient. This is a resource lower bound, not an
additional hole in the tournament spectrum. It holds even if the
selector is allowed to choose its quotient from the target and even if
one limits attention to the blue part.

## 4. A literal signed difference that does not pay the twin-width step

Replace a single quotient vertex by TT_3 or C_3; keep every cross arc
and every other block unchanged. The two hosts differ in one internal
arc. Their local ordered nonempty-cover profiles are (1,6,6) and
(3,6,6). In a full host slot cover, compress the replaced block's maximal
runs within each slot (the slot boundaries remain marked). Skeletons
using that block twice or more cancel between the hosts. Those using it
once are exactly the spanning slot covers of the contracted host, with
one ordinary visit at the contracted vertex. Therefore for every s>=1,

    P_s(cyclic triangle replacement) - P_s(transitive replacement)
        = 2 P_s(contracted host).

This works with arbitrary other blocks and all s simultaneously. For
k independent selected three-vertex modules, the mixed alternating
difference is 2^k P_s of the host with those k modules contracted:
apply the same cancellation in each coordinate. Only the one-coordinate
difference is itself a difference of two literal host counts; do not
identify the k-fold alternating sum with a single physical edge.

When the contracted host is strong and nontrivial, both inflated hosts
are strong. The switch is exact but the resulting Hamiltonian gap is
2H(contracted host), not 2 unless the contracted host is transitive.
For a continuation built using only genuine modules this identifies the
scale factor to cancel, rather than treating the new digit as independent.

**What to try next.** Use the decoder on a profile in which a carried
target digit is a *difference between two physical responses*, not an
increased width in one bag. The exact finite-difference bound explains
why one-sided inflation alone cannot close arbitrary residuals. A
successful selector has to coordinate changes between bags or increase
the quotient complexity, and still terminate at the absolute target.
Chan–Pak's [spectrum conjecture](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf)
is not settled here.

## Reproduction

Run `python3 verify.py` or `python3 -O verify.py`: 120 signed-profile
decodes and rising/falling checks, 330 direct subset-DP cyclic-host
controls, 3,420 finite-difference controls, 285 full-profile switch
controls. These are regression tests; the identities and the complexity
bound are proved above without extrapolating any finite census.
