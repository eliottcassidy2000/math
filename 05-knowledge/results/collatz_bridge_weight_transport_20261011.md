# A refinement-stable source law on the new native Collatz bridges

2026-10-11. **PROVED, scoped and self-audited:** pointwise, setwise, and
conditional distortion bounds for the incoming all-height bridge family.
**INHERITED:** its native common-future words, ROOT223 certificate, and
integer-rank payment. **FINITE-EXACT:** the declared controls below.
**OPEN:** universal entry, positivity of an independently measured ROOT
score at every source, and Collatz. This is a quantitative extension of an
existing family, not a newly discovered ROOT basin.

## Inheritance and the actual increment

The live incoming commit `8300311a5b` supplies
[native run bridges, B2 and B4](collatz_atomic_run_bridge_20261011b.md).
It already proves the all-height family, grounds its t=0 section through
223, pays at least 399 integer-rank units per atomic G bridge, and bounds
the source-averaged number of such bridges. Those results are reused.
The [atomic-rank interface](collatz_atomic_rank_interface_20261011.md)
already controls episode lengths; it explicitly does not give the law
after conditioning on a guard. The present addition supplies that law
for this particular guard, uniformly under all further refinements.

The closest mechanism is the rational prior and source-floor compiler R6
of [renewal and source weights](collatz_renewal_source_20261011.md).
The canonical hostile is an episode ending at an ungrounded child.
The corrected near miss is conditioning an average over all sources on
a sparse guard without recomputing the measure. The least-used sidecar
here is the affine map on the actual source index, including its offset.

Board: **native word / source index / conditional law / retained ROOT
receipt / ramification / lost observer kernel**. The last two are examined
in the companion [source-frame and EIT note](source_frames_eit_20261011.md).

An older, stronger obstruction also remains live:
[THM-4483, forced-charges-two-place-ranks](../../01-canon/theorems/THM-4483-forced-charges-two-place-ranks.md)
already excludes height-summable valuation-bank ranks and identifies the
strict nonnegative-bank existence statement as equivalent to Collatz.
No claim below bypasses that obstruction.

## W1. The native map retains a disjoint branch address

Write U(n)=(3n+1)/2^{v_2(3n+1)} on positive odd integers. For j>=1,t>=0,
put K_j=2^{18j} and

\[
 N_j(t)=\frac{2264K_j-23}{27}+96K_jt,\qquad h(t)=223+256t.
\]

The inherited exact words are

\[
 N_j(t)\xrightarrow{(1,2,18j,1)}425+486t
       \xleftarrow{(1,1,1,1,3)}h(t).                 \tag{1}
\]

They are native words, including the exact valuations. Consequently
N_j(t) reaches ROOT iff h(t) reaches ROOT. For a positive certificate,
retain the child's literal ROOT word and splice its suffix at the common
endpoint; (1) alone does not supply that word.

The branch address can be recovered without ambiguity:

\[
 27N_j(t)+23=8K_j(283+324t),\qquad
 v_2(27N_j(t)+23)=18j+3.                              \tag{2}
\]

Thus all source progressions are disjoint. Within a fixed j the parameter
t is unique. Every source is 75 mod256, whereas every child is 223 mod256.
In particular the same bridge cannot simply be iterated on its own output
by ignoring its guard.

## W2. Pointwise distortion, retaining the affine offset

For any real r>0 define the probability on source indices m>=0 by

\[
 \nu_r(m)=(m+1)^{-r}-(m+2)^{-r}.
\]

Index m represents the odd integer 2m+1. For rational computations below
r is a positive integer. Let

\[
 b(t)=111+128t,\quad a_j(t)=(N_j(t)-1)/2,\quad
 s_j=3K_j/8,\quad d_j=2(K_j-1)/27,\quad
 \theta=9056/9153.
\]

These satisfy the exact index relation

\[
 a_j(t)+1=s_j(b(t)+1)-d_j,\qquad
 0<d_j/s_j<16/81.                                    \tag{3}
\]

**Proposition.** For every such r,j,t,

\[
 s_j^{-(r+1)}\nu_r(b(t))
 \le\nu_r(a_j(t))
 \le\theta^{-(r+1)}s_j^{-(r+1)}\nu_r(b(t)).           \tag{4}
\]

**Proof.** Write H=b(t)+1>=112 and use

\[
 \nu_r(m)=r\int_0^1(m+1+v)^{-r-1}\,dv.
\]

For 0<=v<=1, equation (3) gives

\[
 \theta s_j(H+v)\le s_jH-d_j+v\le s_j(H+v).           \tag{5}
\]

For the lower comparison in (5), the ratio
(s_jH-d_j+v)/(s_j(H+v)) decreases with v and increases
with H. Its minimum on H>=112,0<=v<=1 is therefore at
H=112,v=1. There it is at least

\[
 \frac{112-16/81}{113}=\frac{9056}{9153}=\theta.
\]

The upper comparison follows from d_j>=0 and s_j>=1. Raising (5) to
power -(r+1) and integrating proves (4). QED.

The constants are independent of j, t, any later residue modulus, and
the subset of children chosen. The exponent r+1 is the derivative decay
of the telescoping prior; it is not a heuristic dimension of the orbit.

## W3. Infinite branch sums and the conditional law

Define the target-ray mass and the all-branch pullback mass for any
subset B of the nonnegative integers:

\[
 \mu_r(B)=\sum_{t\in B}\nu_r(b(t)),\qquad
 \rho_r(B)=\sum_{j\ge1}\sum_{t\in B}\nu_r(a_j(t)).
\]

The disjointness in (2) means rho counts actual sources once. Tonelli's
theorem and (4) give, with an explicit geometric sum,

\[
 C_r\mu_r(B)\le\rho_r(B)\le
 \theta^{-(r+1)}C_r\mu_r(B),\qquad
 C_r=\frac{(8/3)^{r+1}}{2^{18(r+1)}-1}.                \tag{6}
\]

Both total masses are finite and positive. Normalize on the full child
ray and the union of all bridge-entry rays, respectively. Their induced
probabilities on t satisfy

\[
 \theta^{r+1}\frac{\mu_r(B)}{\mu_r(\mathbb N_0)}
 \le\frac{\rho_r(B)}{\rho_r(\mathbb N_0)}
 \le\theta^{-(r+1)}\frac{\mu_r(B)}{\mu_r(\mathbb N_0)}. \tag{7}
\]

The same bounds hold when conditioning on any single fixed bridge j.
This follows by normalizing (4), not by treating its source guard as
independent of the prior. Equations (6)--(7) hold for singletons, residue
cells at arbitrary depth, and every subset; thus the comparison survives
refinement. More generally (6) holds with a nonnegative score V(t) in
place of an indicator whenever the right-hand sum is finite.

For r=6 the conditional distortion factor is

\[
 \theta^{-7}=(9153/9056)^7<1.078.                      \tag{8}
\]

This gives a precise version of a reusable move: an authenticated
arithmetic bridge can transfer an independently supplied source-specific
score without the comparison degrading at deeper digit resolutions.
The lower bound remains relative to a genuine target score. It does not
manufacture positivity where that score is still unknown.

## W4. What the zeta(7) connection pays, and what it does not

R6 of the inherited renewal note already proves, for q>m,

\[
 E_q(m)=\sum_{k\ge1}\nu_r(m+kq)
       \le(r+1)q^{-r-1},\qquad
 q^{r+1}E_q(m)\longrightarrow r\zeta(r+1).
\]

Setting r=6 gives 7/q^7 and 6*zeta(7). This is an existing theorem's
specialization, not a new irrationality result. It has a concrete common
mechanism with (4): differentiating a degree-six reciprocal tail gives
power seven. No arithmetic information about zeta(7) is needed for either
inequality. In particular there is no unproved transcendence premise in
the bridge law.

The apparent mass-expansion shortcut fails for a quantitative reason:

\[
 0<C_6\le\theta^{-7}C_6<10^{-34}.                     \tag{9}
\]

The infinitely many very large sources in this family import only a tiny
amount of the fixed probability mass relative to their child ray. Their
large rank payment does not turn that mass into a supercritical gain.
Moreover their guard and their child guard differ. Both losses must be
paid before using repeated pullback to exclude a hypothetical bad set.

The prior [adaptive mixture flow, P5 and Section7](collatz_adaptive_mixture_flow_20261005.md)
already supplies exact ROOT-certificate weight transport and a normalized
injection measure on rooted multiples of three. Its missing full support
cannot be inferred from normalization. Our fixed source prior gives a
different measure, with explicitly controlled conditional distortion;
these two measures have not been identified.

The older [credit/zeta bridge](collatz_credit_zeta_bridge_20261005.md)
also already isolates the missing integral-readout hypothesis in an
Apery-style argument. A rational-valued atomic weight or a sharp tail
bound does not supply integrality of an infinite expectation.

## Next helper questions

1. Which other authenticated bridge banks admit a uniform analogue of
   (5), and can their target guards cover the actual pending-child set?
2. Can a multi-guard return network retain all offsets, overlap corrections,
   and ROOT receipts while giving a source-specific lower score before
   a complete ROOT word is known? Test its signed/negative-cycle analogue
   before treating a generic contraction as evidence of ROOT entry.
3. Is there a computable positive comparison from the fixed nu_r prior to
   the rooted injection measure at every multiple of three? A proof would
   supply the missing full support; finite averages cannot replace it.

These questions strengthen a universally defined prior on a proved
entry substructure. No rule covering every integer has been established.

## Exact reproduction

[Checker](../../04-computation/experiments/source_frames_bridge_weights_20261011.py)
and [saved output](source_frames_bridge_weights_20261011.json) share the
planar controls. Bridge universe: j=1..6; t=0..20,127,1024,10^9;
r=1,2,6,10. The native words are replayed independently of the incoming
compiler. Every subset of t=0..7 is compared after actual finite-box
normalization for j=1,2 and r=1,2,6. Wrong valuation controls are rejected;
the infinite geometric tail is an exact rational identity. General scope
comes from (3)--(7), not from the finite enumeration.

```text
python3 -B 04-computation/experiments/source_frames_bridge_weights_20261011.py
python3 -B -O 04-computation/experiments/source_frames_bridge_weights_20261011.py
```

Both modes byte-agree on **3,036 always-active exact checks** across the
two notes. The local artifact links and agent-facing documentation gates
also pass.
