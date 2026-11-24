# Unit cycle seams, path leakage, and a 24-contact Golay restriction

**PROVED** for the exact seam and its path correction; **FINITE-EXACT**
for the binary code certificate. Inherits the complete global cycle
spectrum from [PROOF.md](PROOF.md), and the related
[binary comb](../cycle-bridge/PROOF.md). The path-spectrum problem is not
closed. We count cycles modulo rotation, not modulo reversal.

## Closing two ends selects a count; opening them leaks paths

The repo already proves that **every nonnegative integer is a Hamilton
cycle count**, with explicit positive realizations in strong tournaments
of feedback vertex number at most two. The operative difference from
paths is not just parity: insisting on a return can force selected
endpoints, and a path need not use that return. We make the leakage
exact on one fixed layout.

Take a transitive spine \(s,1,\ldots,m,t\), two additional vertices
\(a,b\), and contact bits \(u_j\in\{0,1\}\). Orient
\(a\to s,t\), \(s,t\to b\), \(b\to j\) for all interior \(j\), and
\(j\to a\) iff \(u_j=1\). The marked arc is \(b\to a\) when
\(\varepsilon=1\) and \(a\to b\) otherwise. Put
\[
 x=\sum_{j=1}^m u_j2^{j-1},\qquad r=\sum_{j=1}^m u_j2^{m-j}.
\]
Then the cycle count is \(C_\varepsilon=x+\varepsilon\).
Root at \(a\) and write \(a,P,b,Q,a\) with increasing spine words
\(P,Q\). Neither \(s\) nor \(t\) can lie in \(Q\): if present, \(s\)
would be first (but \(s\to b\)) and \(t\) last (but \(a\to t\)).
Thus \(P\) automatically fits, and \(Q=\varnothing\) contributes
\(\varepsilon\); otherwise its largest vertex \(j\) must have
\(u_j=1\), giving \(2^{j-1}\) subsets. For every positive \(x\) the
two realizations are strong and differ by exactly the one marked
arc. This is a literal unit cycle seam at every positive height
with at most \(\lfloor\log_2 x\rfloor+5\) vertices. The zero seam
has the same formula, without requiring strength at zero.

For path counts \(H_\varepsilon\), the identical edit instead satisfies
\[
 H_1-H_0=2^{m+1}+1-r-x+[u_1=0]+[u_m=1]\quad(m\ge1).
\]
Only paths using the marked arc matter. For \(P,a,b,Q\)
in the \(a\to b\) orientation, \(s\) must belong to \(P\) and \(t\)
to \(Q\). The interior vertices in \(P\) must form a nonempty
proper subset with bit 1 on its maximum, giving \(x-[u_m=1]\).
For \(P,b,a,Q\), first suppose \(t\in P\). If \(s\in Q\), all
\(2^m\) interior allocations work. If \(s\in P\), \(Q\) is empty
or has bit 0 on its minimum, giving \(2^m-r\) allocations.
If \(t\in Q\), \(P\) is empty or exactly \(\{s\}\); the latter
works iff \(u_1=0\). Subtraction proves the formula. It can be
rewritten
\[
 4+\sum_{j:u_j=0}(2^{j-1}+2^{m-j})+[u_1=0]-[u_m=0].
\]
In particular the gain is always at least 4, with equality exactly
on the all-one word. For \(m=0\) the gain is also 4 by direct count.
The extra paths are forced by the forward-only spine runs, not a
floating approximation to the cycles. This gives a precise reason
that even a universal literal \(+1\) cycle seam does not lift to a
universal \(+2\) Hamilton-path seam. It is not an obstruction to
path constructions outside this layout. In fact \(T-a\) in this
comb has \(D=1+2^m\) paths for \(m\ge1\). In a path \(P,b,Q\) the empty \(P\) is
impossible since \(b\) loses to \(s\). If \(Q\) is nonempty, \(P\)
must end at \(s\) or \(t\), the only vertices beating \(b\).
Ending at \(s\) gives one path (nonempty interior);
ending at \(t\) gives \(2^m-1\) paths (nonempty suffix of interiors).
There is also the one terminal path. Thus \(H_0\ge D\) by insertion
of \(a\), whereas the exact gain is at most \(2D\), with equality
only on the all-zero word. The ratio on every positive seam is
strictly less than three. This particular mark cannot reach \(7H_0\).

Other arcs need not have the same leakage, even with exactly the
same cycle gain. Flipping contact \(u_1\) gives a literal
\((C,H)=(6,81)\to(7,81)\) on one seven-vertex comb, and
\((3,31)\to(4,33)\) on one six-vertex comb, both checked below.
The marked arc and where it sits matter, not just the cycle gain.

## Leech coordinates as a physical restriction on the selector

Take \(m=24\), and restrict the 24 marked contacts to the extended
binary Golay code \(\mathcal G\), the actual binary code used in
[Leech–Sloane's 24-dimensional construction](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/6187C0A96E44CDF6BBDCD2DF41F615BC/S0008414X00053724a.pdf/sphere-packings-and-error-correcting-codes.pdf).
The orientation on \(a,b\) is **not** one of the 24 contact bits.
The result is a family of physical tournaments on 28 vertices, not
a copy of the lattice inside tournaments. Use the standard cyclic
23-bit generator \(1+z^2+z^4+z^5+z^6+z^{10}+z^{11}\) with an appended
parity bit. The replay verifies its dimension 12, self-duality
and enumerator \(1+759z^8+2576z^{12}+759z^{16}+z^{24}\).
These imply:

* Every prescribed pattern of at most seven selected contacts
  occurs with exactly the full-cube frequency. A failure of
  surjectivity on those coordinates would give a nonzero dual
  codeword of weight at most seven.
* At fixed \(\varepsilon\), the family still has just 4096 distinct
  physical cycle counts, not the \(2^{24}\) counts of the full
  fixed-\(\varepsilon\) register. Different base words differ
  in at least eight arcs on the given labels.
* No two allowed base cycle counts differ by a power of two. To
  add \(2^j\) to a binary word without overflowing, clear a run
  of \(t\) ones starting at \(j\), then set the next bit. The
  difference of the two Hamming weights is \(1-t\), whereas their
  XOR has weight \(1+t\). Both would have to be divisible by four,
  since \(\mathcal G\) is linear and doubly even. Their sum is two,
  impossible. In particular the 4096 unit seams using
  \(\varepsilon=0,1\) are disjoint and give 8192 distinct counts
  in all. Here **doubling** and a **mod-four** restraint manifest
  in the physical contacts; the restraint is not a congruence on
  ordinary Hamilton-path counts.

The local agreement is stronger than it first appears. With uniform
contact bits, any joint count moment
\(\mathbb E[C_\varepsilon^pH_\varepsilon^q]\) with \(p+2q\le7\)
is identical in the Golay and full-cube models. The cycle readout
has degree one; the path readout has degree at most two in the
bit variables, since a Hamilton path uses at most two incident
arcs at \(a\). Reducing monomials to squarefree support gives the
claim. The first difference in the pure cycle moments is explicit:
for weights \(w_j=2^{j-1}\),
\[
 \mathbb E_{\mathcal G}[C_\varepsilon^8]
  -\mathbb E_{\mathrm{cube}}[C_\varepsilon^8]
   =\frac{8!}{2^8}\sum_{O\ {\rm Golay\ octad}}\prod_{j\in O}w_j>0 .
\]
Indeed an eight-coordinate projection is either free or subject
only to even parity from the corresponding dual octad; on an
octad the all-one event has probability \(2^{-7}\) rather than
\(2^{-8}\). This is an operational example of local tests missing
large holes in an exact count spectrum. The code restriction, not
a tournament exclusion theorem, creates these holes. The Mathieu
coordinate symmetry preserves the code but does not preserve the
unequal weights or act by tournament isomorphisms on all selectors.

## Where the six labels and the Monster comparison fit

The inherited [Leech/MRB source audit](../../05-knowledge/results/chanpak_mrb_leech_clock_20261114/MRB_LEECH_SOURCES.md)
records the six labels in the order \(13\to2\to5\to3\to7\to4\),
coming from the old Monster/Baby valuation chain with the 11
branch dropped. The full pure-\(A\) rank-24 list also includes
9 and 25. The prior
[normalization check](../../05-knowledge/results/session32_square_monster_20261023/VALUATION_CYCLE_AND_QUOTIENT.md)
shows that \(4=13\) is not a legitimate modular quotient of the
chain. The [Monster/Baby partition](../../05-knowledge/results/session32_square_monster_20261023/monster/MONSTER_BABY_PRIME_PARTITION.md)
puts \(2,3,5,7,11,13\) in the same increased-valuation stratum;
one cannot substitute 4 for the prime 11 in that group-order
statement. The code supplies a different and exact bridge:
24 physical contacts, genuine Leech input, and a separation
between local tests and absolute coverage. Neither the Monster
nor the Baby Monster is acting on these tournaments.

The [sharp feedback-two maximum](FEEDBACK-TWO-MAXIMUM.md) is a
separate extension of the same spine model: the binary selector
gives a guaranteed interval, but the largest possible cycle count
for that feedback budget comes from a balanced pair of apex cuts,
not from the interval's upper endpoint.

Replay: python3 ports-and-golay-check.py.
