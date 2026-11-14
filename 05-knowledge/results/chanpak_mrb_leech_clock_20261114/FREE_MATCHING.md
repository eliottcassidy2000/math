# The sparsity obstruction survives freely chosen port permutations

This strengthens [the parallel-glue theorem](PARALLEL_GLUE.md): the
matching of left and right ports may now be **any** permutation, and
both width and nesting depth are unrestricted. The resulting class
still has summable distinct \(H^{-119/125}\). It is not a theorem
about arbitrary cross-side reversal matrices or nonuniform connections
to the middle.

Let \(P_{h,\pi}(B)\) be the portal template with left terminals
\(a_1,\ldots,a_h\), middle tournament \(B\), right terminals
\(b_1,\ldots,b_h\), and exactly the reversed edges
\(b_{\pi(i)}\to a_i\) across the two terminal sides. The middle is
nonempty for the construction; the terminal-only version \(Q_\pi\)
has empty middle. The complement of the reversed matching consists
of the standard forward arcs at all inter-region pairs.

## 1. Factorial color gain from all used port edges

The number of Hamiltonian paths of \(P_{h,\pi}(B)\) that use
**every** reversed matching edge is exactly
\[
                  h!\,F_{h+1}(B).
\]
Contract each used reverse pair into an ordered block. The blocks
are ordered arbitrarily, and the middle uses the \(h+1\) slots before,
between and after them. The forward arcs between consecutive
matched pairs cannot themselves be reversed because the pairs are
distinct. A sequence using no middle vertices in a slot is allowed.
Conversely, a path using all the reversed edges necessarily has
this form: every terminal is already accounted for in one used pair.
This works for all matchings and for an empty middle. For \(h=0\)
the formula reads \(F_1(B)=H(B)\).

In the \(m\)-color observable, restrict to colorings in which both
ends of each matched pair share a color. If the numbers of pairs
in colors \(i\) are \(h_i\), restrict the color-\(i\) paths to
use all their \(h_i\) backward arcs. The previous formula replaces
the \(i\)-th color of the middle by \(h_i+1\) distinguished subcolors,
with factor \(h_i!\). Summing these refinements over all middle
colorings produces \(F_{m+h}(B)\). For fixed \(h_i\) the assignment
of pair labels to colors has multiplicity \(h!/\prod h_i!\).
Summing over all weak compositions gives the uniform lower bound
\[
 F_m(P_{h,\pi}(B))\ge(m)^{\overline h}F_{m+h}(B),
 \qquad (m)^{\overline h}=h!\binom{m+h-1}{m-1}.
\]

Among all matches, retain any two pairs together with the middle.
Their induced tournament is a crossing \(P_2(B)\) or a nested
double wrap. The preceding bound with \(h=2\) applies equally to
both. For arbitrary disjoint induced vertex groups, the odd-cycle
packing identity gives a lower bound on \(H\) by the product of
their individual \(H\); sum over colorings to get the analogous
bound on \(F_m\). With \(Q_{\pi'}\) the terminal-only tournament
on the remaining \(t=h-2\) pairs, this yields a second bound
\[
     F_m(P_{h,\pi}(B))
     \ge m(m+1)F_{m+2}(B)F_m(Q_{\pi'}).
\]

## 2. A uniform lower bound on terminal-only paths

Fix \(t\) and a permutation. Select \(k\ge1\) of the \(t\) matched
pairs whose reversed arcs will be used, in any order. Their ordered
block sequence has \(k+1\) connector slots, including the initial
and final slots. Send every unused left vertex to an eligible slot:
the initial slot or one after a selected left endpoint earlier
in left order. Send every unused right vertex to an eligible slot:
the final slot or one before a selected right endpoint later in
right order. In every slot list left then right, each in its
respective increasing order. The only potentially forbidden
crossings inside or at the boundary of slots involve an unused
matched pair with its ends in the same slot; guard against this
by giving such ends distinct slots. The resulting paths use
exactly the chosen \(k\) backward arcs, and are distinct by their
order and vertex assignments.

For an unused pair let \(L,R\) be its two eligibility counts.
There are at least \(LR\sqrt{q(L)q(R)}\) distinct-slot choices,
where \(q(1)=1/2\) and \(q(j)=1-1/j\) for \(j\ge2\).
In fact at most \(\min(L,R)\) same-slot pairs can be disallowed,
and unless \(L=R=1\) this discards at most \(1/\max(L,R)\) of
the choices. In the remaining case the unique eligible slots are
the initial and final ones, which are actually distinct for \(k\ge1\).

For a uniformly selected \(k\)-set of pairs, and an unused focal
pair sampled uniformly, its rank among the selected \(k\) pairs
plus itself is uniform on the \(k+1\) ranks in either of the
two independent linear orders (no claim of statistical independence
between the two ranks is needed). Therefore the average over the
chosen sets of the sum of the logarithmic lower bounds has mean
\[
 (t-k)\frac{2\log((k+1)!)-\log(2(k+1))}{k+1},
\]
because \(\prod_{j=1}^{k+1}q(j)=1/(2(k+1))\). Applying
arithmetic-geometric mean, and multiplying by the \(k!\) orders
of each chosen set, gives
\[
 H(Q_\pi)\ge\frac{t!}{(t-k)!}
       \left(\frac{((k+1)!)^2}{2(k+1)}\right)^{(t-k)/(k+1)}.
\]
For the replay let \(g(t)\) be the maximum of \(t!\) and these
bounds with the parenthesized \((k+1)\)-st root rounded down,
with \(g(0)=1\).

This gives certified lower bounds at moderate widths. For the
infinite tail, the following alternative bound is simpler. Select
the \(k\) leftmost matched pairs by left order, and among the
remaining pairs choose the \(k\) rightmost by right order. Alternate
the two groups, each in arbitrary internal order. All remaining
left and right vertices can go into any of the \(k\) first-group
to second-group connector slots; assigning the two endpoints of
each unused pair to different slots avoids every forbidden matched
forward crossing. This gives
\[
  H(Q_\pi)\ge(k!)^2[k(k-1)]^{t-2k}, \quad 2\le k\le t/2 .
\]
For \(h\ge256\), set \(t=h-2,\ k=\lfloor t/4\rfloor\).
Using \(k,k-1\ge h/5\), \(t-2k\ge t/2\), and
\(k!\ge(k/3)^k\) gives the lower bound
\[
   G_{\rm tail}(h-2)=(h/15)^{2h/5}(h/5)^{99h/100}.
\]
Also \(h!\le3h(2h/5)^h\), by the upper integral bound
on the sum of logarithms and \(5/2<e<3\).
For \(s=119/125\) rational exponent comparison at \(h=256\)
then gives
\[
       h! G_{\rm tail}(h-2)^{-s}\le3h(7/8)^h
       \quad(h\ge256).
\]
Indeed the requisite exponential base
\((2/5)15^{(2/5)s}5^{(99/100)s}/h^{(139/100)s-1}\)
decreases with \(h\), and at 256 is \(<7/8\), checked by
raising rational quantities to power 12500. The bound summed
over \(h\ge256\) is \(<10^{-9}\).

## 3. An exact small-color certificate for unbounded width and depth

For H and F_m we may count all descriptions. To avoid duplicating
an outer wrap as a larger portal width, if \(\pi(1)=h\) pull out
the pair \(a_1,b_h\) as a \(P_1\) and treat the rest as a smaller
portal tournament. Thus the number of necessary matching
symbols at width \(h\ge2\) is \(n_h=h!-(h-1)!\).
There is one symbol at \(h=1\). An atom is a base terrace
or singleton, or a portal operation on a nonempty ordinal
sequence of atoms. Compare sums of weights on descriptions with
bounded maximum nesting depth and then increase that depth.

Write \(S_m\) for the sum of \(F_m^{-s}\) over atom descriptions.
We prove \(S_3<.99,\ S_4<.49,\ S_5<.35\), and \(S_m<.30\)
for all \(m\ge6\). The strict finite certificate bounds
all depths inductively, so the limits also satisfy the displayed
bounds with nonstrict inequalities, which suffice for \(S_m<1\).
For a fixed width \(h\ge2\), improve the two color lower bounds
by writing
\[
 F_m(Q_\pi)\ge G_m(t):=
 \max\left\{
  m^{2t},\ 
  (m^4+4m^2)^{\lfloor t/2\rfloor}m^{2(t\bmod2)},\
  \max_{1\le j\le t}m^{2(t-j)+1}g(j)\right\}.
\]
The first is trivial; the second keeps disjoint doublets of
port pairs, each of which is a four-vertex tournament with
two directed triangles; the third restricts a chosen \(j\)-pair
subset to a single color and leaves every other terminal's color
arbitrary, using induced monotonicity. Choose the \(h\)-dependent
coefficient
\[
 c_{m,h}=\max\{(m)^{\overline h},\ m(m+1)G_m(h-2)\}.
\]
For the tail use only the \(G_{\rm tail}\) lower bound,
and interpret \(G_m\) there as at least that bound. In both cases
\(F_m(P_{h,\pi}(B))\ge c_{m,h}F_{m+2}(B)\), since \(F_m\)
is nondecreasing with the number of colors.

Define \(K_m=\sum_{h\ge2}n_h c_{m,h}^{-s}\).
It is nonincreasing in \(m\) and finite by the tail bound.
By exact rounding using only integers,
\[
 (K_3,K_4,K_5,K_6)<(.255,.083,.049,.034).
\]
In the upper-bound recurrence for \(S_m\) the base contribution
is at most \(A(m^{-s})\), with
\(A(z)=z+z^3+z^4/(1-2z)\).
At width one the preceding depth contributes at most
\(m^{-s}U_{m+1}/(1-U_{m+1})\); at widths at least two,
at most \(K_m U_{m+2}/(1-U_{m+2})\). Here
\[
 U_3=.99,\ U_4=.49,\ U_5=.35,\ U_m=.30\ (m\ge6).
\]
The exact replay checks these sums are strictly less than
\(U_3,U_4,U_5,U_6\) at \(m=3,4,5,6\) respectively.
Monotonicity of \(A(m^{-s})\), \(m^{-s}\) and \(K_m\) handles
all \(m\ge6\). This proves the induction for every finite
nesting depth. In particular \(F_3^{-s}\) and \(F_4^{-s}\)
sum over all ordered nonempty reservoirs by the geometric series.

To finish, the sum of \(F_2^{-s}\) over atom descriptions is finite:
on base classes use the [earlier one-wrap/terrace estimate](https://github.com/eliottcassidy2000/math/blob/edf094e6bc983f5e8c8599aecdd2e0c00783773c/05-knowledge/results/chanpak_frame_decoder_pressure_20261114/WRAP_CHANNELS.md);
at width one use \(F_2\ge2F_3(B)\); at \(h\ge2\) use
\(F_2\ge6F_4(B)H(Q_{\pi'})\). The width sum after this last
inequality converges by the proved tail despite the \(h!\)
many matchings. For ordinary widths at least two, use
\(H\ge2F_3(B)H(Q_{\pi'})\). Width-one ordinary outputs depend
only on the multiset of the constituent atom \((F_1,F_2)\) pairs
and are at least their \(F_2\) product. Their sum is therefore
bounded by a convergent multiset Euler product from the
\(F_2^{-s}\) atom sum. The base distinct \(H>1\) sum is finite at
\(s\) by the terrace estimate. A last Euler product over the distinct
atom totals \(H>1\) gives the claimed result for arbitrary
ordinal products. The unit factor has no effect and is omitted.

Thus neither changing the number of parallel ports arbitrarily nor
freely permuting their matching removes the quantitative obstruction.
The next candidate must go beyond disjoint reversed pairs between
two transitive boundaries: e.g. allow a terminal to participate in
several backwards choices, or let different terminals see different
subsets of the reservoir. This conclusion is conditional on none
of the speculative lattice analogies.

Run `python matching_certificate.py` for the exact inequalities.
