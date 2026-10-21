# Linked tail geometry and exact two-block terminal constraints

Status: **PROVED** source-state conditional bounds and their comparison rules;
**FINITE-EXACT** consumers of the four frozen stopped-prefix states below.
Uniform native `153N` coverage and arbitrary Collatz entry remain **OPEN**.
No source is advanced here, and no pending suffix is read.

## 1. Inheritance and the retained information

The nearest proved mechanism is the synchronous four-vector envelope in
[future tail blocks](../future_tail_blocks_20261021b/README.md). Its pending
blocks are fibres of the **same current native P node at the same depth**,
not arbitrary zero-sum selections. The independent incoming
[route-free suffix gate](../routefree_suffix_gate_20261021/README.md) and its
[mixed zonotope refinement](../routefree_suffix_gate_20261021/session22-residual/MIXED-COROLLARY.md)
already prove signed-sum domination with arbitrary future destinations and
nonarrival. Here that positive-norm mechanism is applied to *linked,
coalesced old/new blocks*, through an explicit rational change of quadratic
coordinates. The all-future stopping idea is inherited.

The canonical hostile is N7 at cut8. The corrected near miss is to erase
opposite charges at different pending nodes, or to maximize an indefinite
quadratic form only at the vertices of a zonotope. The least-used sidecar
is the common continuation sign of a block's **two** channels. The live
concept board is: literal source / common current node / linked old-new
channels / signed angular geometry / terminal assignment / running selector.

The [cross-scale obstruction](../tail_scale_certificate_obstruction_20261021/README.md)
still applies: these estimates prove neither a fixed `D+6` scale law nor
old-source completion on every scale. Old completion along all scale levels
would itself imply universal native entry. No such premise is used here.

Fix a finite stopped-prefix state at depth D. Every terminal t in the fixed
627-point frontier is retained, including empty ports, with
`c_t=(W_t,H_t)`. Every distinct pending current node has a source-disjoint
block `b_i=(w_i,h_i)`. The old/new coordinates include actual past signs.
Write

\[
 Q(a,b)=a^2-ab+b^2,\qquad
 \beta(v,z)=Q(v+z)-Q(v)-Q(z),
\]
\[
 \delta(w,h)=Q(h)+\beta(w,h)=Q(w+h)-Q(w),
\]
\[
 \kappa((w,h),(w',h'))
 =\beta(h,h')+\beta(w,h')+\beta(w',h).
\]

The present physical increment is `Delta_D=sum_t delta(c_t)`. At any later
cut, each pending block is absent or contributes its entire four-vector
with one sign to one terminal. These choices may be relaxed independently
for a sufficient upper bound; actual dynamics is a subset. In particular
nonarrival is allowed forever. The complete source marks, native P clocks,
and joint mod15/mod35 registers remain in the inherited carrier `(N,D,s)`.
No change of charge coordinates below is a conjugacy of P or a new source.

## 2. A linked positive majorant for the indefinite form

For any positive rational lambda set

\[
 g_\lambda(w,h)=\lambda w+\frac{\lambda+\lambda^{-1}}2h,
 \quad
 f_\lambda(w,h)=\lambda w+\frac{\lambda-\lambda^{-1}}2h.
\]

Direct expansion gives the two exact identities

\[
 \delta(b)=Q(g_\lambda(b))-Q(f_\lambda(b)),\qquad
 \kappa(b,c)=\beta(g_\lambda(b),g_\lambda(c))
             -\beta(f_\lambda(b),f_\lambda(c)).            \tag{1}
\]

This is a rational hyperbolic change of the positive-minus-positive
quadratic coordinates. Its use here needs only the displayed algebra.
Both channels of each source fibre stay linked. If lambda=p/q, the checker
uses integer numerators

\[
 G_b=2p^2w+(p^2+q^2)h,\quad F_b=2p^2w+(p^2-q^2)h,
 \quad g=G_b/(2pq),\ f=F_b/(2pq).                         \tag{2}
\]

For planar generators z_i define the signed-sum maximum

\[
 \mathcal D(z)=\max_{\epsilon_i\in\{-1,1\}}
                    Q\!\left(\sum_i\epsilon_i z_i\right).
\]

**Positive-norm destination lemma.** Arbitrary disjoint subsets assigned to
arbitrary terminals, with arbitrary signs and omitted members, satisfy

\[
 \sum_t Q\!\left(\sum_{i\text{ assigned }t}\epsilon_i z_i\right)
 \le\mathcal D(z).                                      \tag{3}
\]

To prove this, multiply each entire terminal sum by an independent symmetric
sign. The expectation of the norm of their sum is the left side. Add each
omitted generator with an independent symmetric sign; its addition increases
that expectation by Q(z_i). Every resulting complete signed sum is bounded
by D. These are proof signs, not altered native clocks. This is the positive
norm argument inherited from the incoming gate.

The zonotope `sum_i[-z_i,z_i]` is planar. Convexity of Q places a maximum
at a vertex. Orient all nonzero generators into one half-plane, sort by
exact determinant signs, start at minus their sum, and successively add
twice each generator; also retain the opposite boundary. This lists every
vertex. Collinear intermediate points may be retained harmlessly. At most
`2m+2` points are evaluated for m nonzero generators, regardless of their
directions. **No twelve-direction assumption is made for coalesced blocks.**

Set `m_i=max_t |kappa(c_t,b_i)|` over all marked terminals and `L=sum_i m_i`.
By (1), (3), and the exact current/pending polarization, every later cut obeys

\[
 \Delta_h\le \Delta_D+\mathcal D(g_\lambda(b_i))+L,
 \qquad h\ge D.                                        \tag{4}
\]

Indeed the pending quadratic part is
`sum_t [Q(sum g)-Q(sum f)] <= sum_t Q(sum g)`; its current-port interaction
is at most L. This proof does not assume nonnegative delta_i, an eventual
arrival, or completion of the old channel. Integer arithmetic may floor
the rational right side because the physical increment is an integer.

When every **aggregate pending old vector** w_i is zero, lambda=1 gives
g_i=h_i and

\[
 \mathcal D(h_i)
 =\sum_iQ(h_i)+\max_\epsilon\sum_{i<j}\epsilon_i\epsilon_j\beta(h_i,h_j)
 \le\sum_iQ(h_i)+\sum_{i<j}|\beta(h_i,h_j)|.
\]

Thus (4) is no worse than the prior U_D in that case. This is a vector
condition, not an assertion that every individual old source has arrived.
For general mixed channels neither bound dominates the other universally.
In the zero-pending-old case, lambda1 is globally best within this particular
projection family: `D(g_lambda)=((lambda+1/lambda)/2)^2 D(h)` and L is unchanged.
This is not an optimality claim among all future bounds.

## 3. A computable sufficient coalescence condition

Define

\[
 P=\sum_i\sqrt{Q(g_i)},\ R=\sum_i\sqrt{Q(f_i)},\quad
 G=\max_t\sqrt{Q(g(c_t))},\ F=\max_t\sqrt{Q(f(c_t))}.
\]

Triangle inequality and Cauchy--Schwarz applied to (1) give

\[
 \mathcal D(g_i)\le P^2,\quad L\le2GP+2FR,
 \quad\Delta_h\le\Delta_D+P^2+2GP+2FR.                 \tag{5}
\]

This supplies an explicit conditional uniform route. At any N, suppose a
finite D and a positive rational lambda satisfy

\[
 \Delta_D\le(153-\varepsilon)N,\quad
 P\le a\sqrt N,\ R\le b\sqrt N,\ G\le c\sqrt N,\ F\le e\sqrt N,
 \quad a^2+2ca+2eb\le\varepsilon.                       \tag{6}
\]

Then every h>=D satisfies the original153N bound. To conclude all horizons
uniformly for all N>=N0 one still must prove these hypotheses for every
such N, and cover the earlier cuts. Neither assertion is supplied by four
examples. D and lambda may depend on N; no old-completion premise is hidden
in (6). The quantities are **norm-weighted coalesced perimeters**, not just
numbers of nodes or numbers of source addresses.

`perimeter_majorant` gives a conservative exact rational version of (5):
replace each numerator norm by its integer ceiling square root and retain
the common denominator2pq. This is independently computable from the
finite state. It is weaker than the angular bound, but its four radius
quantities state what a uniform coalescence proof would need to control.

## 4. Exact cluster bounds retain destination conflicts

For a set C of pending node IDs let H(C) be the exact maximum of

\[
 \sum_{i\text{ active}}[\delta_i+\sigma_i\kappa(c_{t_i},b_i)]
 +\sum_{i<j\text{ active},\ t_i=t_j}\sigma_i\sigma_j\kappa(b_i,b_j)
                                                               \tag{7}
\]

over absence, signs and marked terminal choices. For any partition of the
pending nodes, a valid all-future bound is

\[
 \Delta_D+\sum_C H(C)
 +\sum_{i<j\text{ in different clusters}}|\kappa(b_i,b_j)|.       \tag{8}
\]

Each original interaction occurs exactly once. Singleton clusters recover
U_D. Merging two clusters cannot increase (8), because the new joint maximum
is at most the sum of the old maxima plus their intervening absolute pair
bills. This monotonicity concerns **partition refinement at one fixed
state**, not the independently chosen partitions at different horizons.

For two blocks (7) has a cheap exact solver: compare neither, either alone,
both at the same terminal with all four sign choices, and both at different
terminals. The last case uses the two best absolute current-port supports
for each block; if their best destinations agree, one uses its second best.
The maximum includes absence even when both nonzero delta values are
negative. The checker retains a maximizing pair of marked choices.

Selecting disjoint positive-gain pairs yields a certificate with exact
savings against U_D. The producer uses a deterministic greedy selection;
**it does not claim an optimal matching or an exact global future maximum**.
The consumer accepts any disjoint actual-node pairing and recomputes every
local maximum. Reusing a node in two selected pairs is rejected.

## 5. Hostiles, selector and unchanged ownership

Let `b=((-5,0),(8,0))` and take initially empty current ports.
Here delta(b)=-16. Maximizing delta only at the two vertices ±b gives -16,
whereas optional nonarrival gives0. Convexity of Q cannot be transferred to
the indefinite delta without the majorant in section2.

With two distinct pending copies of b, lambda2 makes g(b)=0, hence (4)
gives0; lambda1 instead gives36, and U gives32. After one copy arrives,
the current value is-16 and the remaining current-port support is32, so
the lambda2 bound becomes16. **The new gate is sound but not inherently
horizon-monotone.** Do not erase the previous certificate when advancing.

Use the minimum of U_D and every source-authenticated all-future bound
previously computed at a cut r<=D, including (4), (8), and the incoming
mixed gate with its own source ledger. This selector bounds every h>=D.
It is nonincreasing when certificates are retained and U_D is advanced.
A minimum is a choice of complete upper certificates, not added credit.

The strict two-block example with current ports `(0,0,1,0)` and `(0,0,0,1)`,
and pending vectors `(-1,0,-1,-1)` and `(-1,1,1,-1)`, has an old local bill8
but exact local maximum5. One maximizing receipt omits a block. In this particular
two-block example an all-active optimizer also attains5; the single-block
negative-delta example above proves why absence cannot generally be removed.
Separately, opposite new
unit charges at different pending nodes cannot be erased, because just one
may arrive at the observed cut.

For an authenticated ownership ledger `Delta_h=A_h+R'_h+J'_h` with
`A_h+M_h=153N`, any complete bound B supplies

\[
 R'_h\le B-A_h-J'_h.
\]

When B<=153N this is precisely the compensated target `R'<=M-J'`.
No pair interaction, source norm, or negative credit has been removed or
claimed a second time. Neither the physical increment nor R' is improved
by changing which complete certificate is selected.

## 6. Exact frozen-prefix comparison

The sole native input is the independently audited finite-prefix receipt
from `future_tail_blocks_20261021b`, LF SHA256
`7731a9a9ddf8e5f0478b8d6455e53c6130870f2004353db2e0fefa345bb509c0`.
Its fixed frontier hash is checked separately. The consumer reconstructs
all627 marked current ports, retains all distinct pending node IDs, checks
the same source hashes and witness clocks/registers, and independently
recovers U_D. This hash-pinned handoff relies on that packet's source/edge
replay; a hash alone is not a replacement for its mathematical audit.

| N | D | Pending blocks | U_D | lambda1 angular bound | Pair gain | Best selected bound |
|---:|---:|---:|---:|---:|---:|---:|
|7|12|5|1153|1153|1|1152|
|4096|105|195|557164|551844|180|551844|
|8192|110|348|1138407|1118881|201|1118881|
|32768|130|760|4241482|4142214|1269|4142214|

The angular gains at the positive scales are5320,19526,99268. Their pending
aggregate old vectors are zero, so the domination proved in section2
applies. The selected two-block pairs number45,75,319 respectively. N7
still exceeds its budget1071 and has an earlier actual violation. The three
positive cases retain their inherited earlier-cut certificates and pass
without completing their pending trajectories. No new N or new native cut
was sampled. No claim that the best listed bound is the exact future
maximum is made.

The coarser ceiling-root perimeter bounds (5) are609156,1276111,4759272
at the three positive scales. They independently pass at4096 and32768 but
miss the8192 budget1253376. The exact angular bound passes there; source
counts or radius sums alone can lose relevant angular information.

## 7. Verification and boundary

```text
python -B 05-knowledge/results/tail_geometry_20261021c/check.py --write
python -B -O 05-knowledge/results/tail_geometry_20261021c/check.py
```

Controls independently enumerate all full sign assignments for every
multiset of at most four vectors from `{-1,0,1}^2`. For each ordered pair
of four-vectors in `{-1,0,1}^4`, they enumerate every optional sign/two-port
assignment, compare the exact pair optimizer, and check three rational
projections against the independent physical energy formula. Both
hyperbolic identities are checked. Literal bool/float, wrong-source,
pair-owner reuse and forged-certificate controls are retained. The native
receipt carries each projection and selected pair certificate, not merely
its final scalar.

The next unproved obligation is source-specific control of the coalesced
geometry in (6), or a stronger computable certificate that eventually
passes at every N. More angular freedom or a better local matching alone
does not prove that global stopping statement.

Normal generated receipt: **76,851 checks**. Optimized replay: **76,852
checks**, including the final saved-equality gate. Both reproduce LF SHA256
`fe7d8190a61e3e883f1b82ad38d6d11a9abeebc9da0833aafa2eb5fa9c16aae9`.
Parent proof review of sections2–5 passed; independent peer audit is pending.
