# Weighted destination-coupled tail bounds

Status: **PROVED** finite-state upper-bound and allocation theorems;
**FINITE-EXACT** consumers of four frozen native prefixes; **OPEN** uniform
stopping/153N coverage. No unresolved future suffix, new source census, or
ROOT completion is used. Supplied finite linking paths in §4 are replayed
exactly. This packet changes an upper relaxation, not physical
energy or ownership of a source interaction.

## 1. Inheritance, target, and retained state

The closest proved mechanism is the linked old/new Lorentz gate in
[tail_geometry_20261021c](../tail_geometry_20261021c/README.md), fed by
[future_tail_blocks_20261021b](../future_tail_blocks_20261021b/README.md).
The current incoming review
[session22-latest-tail-review](../session22_tournament_tail_synthesis_20261021/session22-latest-tail-review/REVIEW.md)
does not close its uniform source-entry obligation. Its independent angular
replay is distinct from replaying the native source trajectories.

Our anchor is a finite-prefix certificate for all later horizons, toward
the exact annular increment budget153N. The niche is a fractional allocation
of unary destination costs. The exploratory matrix direction has a sharp
boundary: on a pure-new pending subspace, every positive quadratic majorant
restricts to a form at least as large as the existing Q. Changing only that
metric cannot improve its exact signed-sum maximum. The present improvement
instead retains **which current terminal supplies each field value**.

The live objects are linked channels, literal terminals, optional nonarrival,
positive angular norms, rational factor allocations, and actual source
clocks. The canonical hostiles remain: negative mixed-channel increment does
not force arrival; equal vectors at different current nodes cannot be merged
without an actual common future; and N7 is not a positive153N control.
The least-used sidecar here is the joint two-coordinate field at a common
terminal, rather than the two separate marginal support values.

Write Q(a,b)=a²−ab+b², with polarization β. A linked block is
b=(w,h), where w is its old signed charge and h its new signed charge. Set

\[
 g(b)=w+h,\qquad
 \delta(b)=Q(h)+\beta(w,h)=Q(g(b))-Q(w),
\]
\[
 \kappa(b,b')=\beta(h,h')+\beta(w,h')+\beta(w',h)
             =\beta(g(b),g(b'))-\beta(w,w').                 \tag{1}
\]

A state has current terminal rows c_t, with every actual frontier terminal
retained, and n separately marked pending nodes b_i. The prefix increment is
Δ=Σ_t δ(c_t). A future pending choice x_i is either absent or (t,σ), σ∈{−1,1}.
Actual deterministic evolution chooses some subset of these assignments, so
maximizing over the entire finite set is a sound relaxation even when a
source never reaches the frontier. Define f_i(t)=κ(c_t,b_i).

The exact assignment objective is

\[
 \Delta+\sum_{i\ active}[\delta_i+\sigma_i f_i(t_i)]
 +\sum_{i<j\ active,\ t_i=t_j}\sigma_i\sigma_j\kappa_{ij}.   \tag{2}
\]

This is a map from a **fixed source-prefix state** to an upper certificate.
It preserves the full linked vectors, node/terminal IDs, first-arrival
boundary, and choice of future common sign. It relaxes actual reachability,
future synchronization, and sign constraints. Those lost facts cannot be
recovered from mod15/mod35 alone. The frozen source receipt retains those
registers, each witness's full integer primitive-P clock and edge-clock
word, and hashes of all source rows and prefix edges. P clocks are not
silently replaced by literal fully-odd-U clocks.

## 2. A weighted pair-factor theorem

Let

\[
 D_g=\max_{\epsilon_i\in\{-1,1\}}Q\!\left(\sum_i\epsilon_i g_i\right).
\]

The inherited positive-norm destination lemma gives
Σ_t Q(Σ_{i→t}σ_i g_i)≤D_g, including omitted blocks. One proof randomizes
an independent sign on each occupied destination and on each omitted
generator: the expectation adds the destination norms and nonnegative
omitted diagonal terms, while every complete sign word is bounded by D_g.

Fix rational α∈[0,1]. Decompose (2) into α times that positive destination
norm and a residual objective. Its unary and pair terms are

\[
 l_i(\varnothing)=0,\qquad
 l_i(t,\sigma)=\delta_i-\alpha Q(g_i)+\sigma f_i(t),
\]
\[
 q_{ij}(x_i,x_j)=
 \begin{cases}
 \sigma_i\sigma_j[\kappa_{ij}-\alpha\beta(g_i,g_j)],
      &x_i,x_j\text{ active at the same }t,\\
 0,&\text{otherwise}.
 \end{cases}                                                \tag{3}
\]

For n≥2 assign the unary l_i equally to its n−1 incident pairs. Then every
actual future assignment satisfies

\[
 \boxed{\ B_\alpha=\Delta+\alpha D_g+
 \frac1{n-1}\sum_{i<j}m_{ij}(\alpha)\ },                    \tag{4}
\]
where

\[
 m_{ij}(\alpha)=\max_{x_i,x_j}
 [l_i(x_i)+l_j(x_j)+(n-1)q_{ij}(x_i,x_j)].                   \tag{5}
\]

This is proved by exact objective decomposition followed by independent
factor maximization. A node appears in several factors, but the total
weight of its unary is exactly1; each original pair coefficient is exactly1.
There is no withdrawal of source norm, repeated payment, or mutation of R′.
Physical increments are integers, so floor(B_α) is still an upper bound.
For n0 or1 the checker uses the exact zero/single-block formula separately.

Each factor has at most eight affine candidates in α:

* both absent;
* only i, or only j, with its maximum absolute terminal field;
* both at the same literal terminal, for each of four sign pairs;
* both at different terminals, maximizing the two absolute fields subject
  to the terminal IDs being different.

The last case needs only the top two terminal entries of each support list.
The same-terminal case must maximize σf_i(t)+τf_j(t) at a common t; separate
maxima are generally incompatible. Hence B_α is rational, convex, piecewise
linear. Its exact minimum on[0,1] is attained at an endpoint or an envelope
breakpoint. `optimize` sweeps the exact slope jumps and retains left/right
subgradients; it does not use floating-point fitting.

For mixed blocks (3) can be indefinite and have negative diagonal terms.
Absence is retained in every local factor. For example b=(−5,0,8,0) with
zero current fields has δ=−16; its correct optional maximum is0, not−16.

## 3. Exact native improvements and what produces them

These are the **same** four pinned prefix states used by the previous gate.
There is no new route discovery or larger stopping census.

|N|P-depth D|Pending nodes|Prior U|Prior angular|Exact α|New bound|Gain over angular|
|---:|---:|---:|---:|---:|---:|---:|---:|
|7|12|5|1153|1153|0|1149|4|
|4096|105|195|557164|551844|149/194|551169|675|
|8192|110|348|1138407|1118881|310/347|1118290|591|
|32768|130|760|4241482|4142214|641/759|4133443|8771|

The positive budgets are626688,1253376,5013504. N7's1149 still exceeds1071;
this neither certifies N7 nor refutes the original target restricted to
N≥4096. Earlier cuts of each positive case remain covered by the inherited
source-prefix receipt. Native provenance is inherited from its pinned
replay; this consumer independently checks the exact algebraic state and
its retained source/clock witnesses.

All four finite states happen to have w_i=0 in each pending aggregate. In
that case α1 recovers the old angular gate exactly. Here is an explicit
accounting of the improvement. Let L_i=max_t|f_i(t)|, L=ΣL_i, and

\[
 G_0=\sum_iQ(h_i)+\sum_{i<j}|\beta(h_i,h_j)|-D_g\ge0.
\]

Define the nonnegative local relaxation deficit

\[
 d_{ij}(\alpha)=(1-\alpha)
 [Q(h_i)+Q(h_j)+(n-1)|\beta(h_i,h_j)|]+L_i+L_j-m_{ij}(\alpha).
\]

Then (4) is exactly

\[
 B_\alpha=\Delta+D_g+L-\Gamma_\alpha,\qquad
 \Gamma_\alpha=\frac{\sum_{i<j}d_{ij}(\alpha)}{n-1}
                        -(1-\alpha)G_0.                    \tag{6}
\]

The checker independently sums these per-type deficits with the actual
pair multiplicities. Identical full four-vectors may share an algebraic
calculation; their marked nodes are not identified. Equation(6) measures
lost slack in the upper relaxation, **not physical negative pair credit**.
Optimizing α retains α1, so it cannot worsen these pure-new angular bounds.

A sharp cheap joint-field hostile has pending h_i=(1,0), h_j=(0,1), zero old
channels and κ_ij=−1. Current h ports (3,3),(0,0) yield local maximum7;
ports (1,−1),(−1,1) yield9. Both have exactly the same separate supports
L_i=L_j=3. Marginal support values therefore do not determine the pair
factor. The literal joint destination coordinate is indispensable.

## 4. Composition with authenticated optional pulses

The separate [lagged pulse packet](../lagged_tail_pulses_20261022/README.md)
supplies source-disjoint pending-node pairs u,v, a literal native path
P^ell(u)=v avoiding the frontier, its full integer clock a, and
b_v+(−1)^a b_u=0 in **both** linked channels. If v first arrives at native
delay τ with subsequent sign σ, their combined contribution is0 beforeτ,
σb_v on[τ,τ+ell), and0 afterward. If v never arrives, neither does u and
the contribution stays0. Thus one typed optional pulse replaces the two
pending owners for any fixed-cut choice upper bound. It is not a
synchronous native node and may not be erased before its pulse ends.

This consumer independently replays all supplied finite path triples,
checks exact integer types, frontier avoidance, clocks and distinct owner
IDs, and reconstructs the exported pulse relaxation. It does not search for
new paths or examine any unresolved suffix. Every retained source hash and
prefix-clock handoff remains attached to the original state.

Applying (4) to those typed optional pulses gives:

|N|Pulses|Remaining optional objects|Pulse angular|Exact α|Weighted pulse bound|Additional gain|
|---:|---:|---:|---:|---:|---:|---:|
|7|0|5|1153|0|1149|4|
|4096|39|156|521233|22/31|520738|495|
|8192|63|285|1039528|247/284|1039055|473|
|32768|144|616|3856740|497/615|3848613|8127|

These are **compositions of certified relaxations**, not sums of gains
from bounds on incompatible states. The original sources/energy are
unchanged. The use of an early-node integer as the pulse's serialization
key does not change its explicit optional-pulse type.

## 5. Adaptive weights are a genuine finite LP

For a **fixed** rational α, let w_{i,ij}≥0 with Σ_{j≠i}w_{i,ij}=1.
Replace (4)'s equal split by

\[
 B(\alpha,w)=\Delta+\alpha D_g+
 \sum_{i<j}\max_{x_i,x_j}
 [w_{i,ij}l_i(x_i)+w_{j,ij}l_j(x_j)+q_{ij}(x_i,x_j)].        \tag{7}
\]

The same exact decomposition proves the bound. Minimizing the final sum
is a finite rational linear program: introduce y_ij and impose one
epigraph inequality per marked pair choice, together with the allocation
simplexes. An exact maximizing-pair oracle supplies all the constraints
using the same shared-terminal/top-two decoder. The feasible weight domain
is compact; the convex finite maximum is continuous, and an optimum can be
chosen rational. Native weights in §3 are **not** claimed to solve this LP.
Joint optimization of α and these particular weights is not asserted to be
the same LP: the displayed parameterization contains products αw.

A useful dual certificate has a probability distribution π_ij over each
factor's marked choices and numbers z_i satisfying

\[
 z_i\le\mathbb E_{\pi_{ij}}l_i(x_i)\quad\text{for every }j\ne i.
\]

Its lower bound for the residual LP objective is
Σ_ij Eπ q_ij+Σ_i z_i. This follows directly by multiplying the epigraph
inequalities and using the allocation sums. A single complete assignment
gives a feasible dual by placing every π on its corresponding two choices
and taking z_i=l_i(x_i). This is also an independent lower-bound mechanism.

There is a small exact positive control. With current ports
(0,0,1,0),(0,0,0,1), take pending

\[
 b_1=(0,0,1,0),\quad b_2=b_3=(1,0,-1,0).
\]

The best equal-split bound is6, at α0. Allocate b1's unary to pair12, b2's
to pair23 and b3's to pair23. All three interactions remain, including
pair13 whose two unary weights are0. Equation(7) gives5. A fully enumerated
actual assignment attains5, proving the adaptive LP optimum **for this
specific abstract state**. It is not claimed to be a native Collatz state.

`allocation_bound` authenticates every unordered pair ID exactly once,
nonnegative exact rational weights, and each unary's total mass1. Missing,
repeated or underallocated factors are rejected. This is a lossless
allocation interface; it does not claim that arbitrary supplied weights
are optimal.

## 6. The remaining uniform target and stopping logic

For each source scale N, construct a genuine finite stopped-prefix state
with complete old/new source ownership, literal frontier and actual clocks.
Any feasible α,w with B(α,w)≤153N certifies **all** later horizons, including
nonarrival. Earlier horizons need their retained direct bounds or another
complete certificate. Thus the quantified sufficient condition is:

> For every N in the claimed universe, there is such a finite prefix,
> together with a checked allocation whose bound is at most153N, and a
> certificate covering every earlier cut.

This existence condition is **OPEN**. It is not supplied by the finite
table. A more concrete pure-new sufficient inequality, if that channel
hypothesis is actually authenticated at the chosen cut, is

\[
 \frac{\sum d_{ij}(\alpha)}{n-1}-(1-\alpha)G_0
 \ \ge\ \Delta+D_g+L-153N.                                 \tag{8}
\]

It asks for enough **joint terminal-field incompatibility** to cover the
remaining angular excess. No cancellation of physical source norm or R′
is presumed. Pure-new aggregate conditions are not silently imposed at
all scales; the mixed theorem(4)/(7) remains available instead.

Keep the minimum of all already authenticated complete all-future bounds
at earlier cuts, together with the inherited monotone U gate. No horizon
monotonicity of newly chosen allocations is needed or claimed. In
particular the inherited λ2 gate's nonmonotonicity hostile remains relevant
to the stopping interface. One may take minima of whole certificates;
one may not add their apparent gains.

For an existing ownership decomposition Δ_future=A+R′+J′ with A+M=153N,
this packet only gives R′≤B−A−J′. It does not alter A, J′, R′, any pair
address, or the original moving credit. A paying complete gate yields the
compensated inequality R′≤M−J′. It yields R′≤M only with the additional
sign premise J′≥0; that premise is not assumed here. Uniform existence of
the paying gate, and hence of that compensated bound, remains OPEN.

## 7. Reproduction and bounded validation

Run from the repository root:

```text
python -B 05-knowledge/results/weighted_tail_20261022/check.py
python -B -O 05-knowledge/results/weighted_tail_20261022/check.py
```

`--write` regenerates only this directory's receipt. The consumer pins the
LF-normalized old geometry code, prefix receipt, frontier and pulse receipt;
the latter authenticates the supplied finite path bank. The bounded
controls enumerate the complete optional signed pair assignments for six
specified linked vectors, one/two current ports, degrees1/2/5 and four
rational parameters; all64 ordered three-block states over the declared
four-vector bank; exact native envelope breakpoints; and an adaptive
allocation with a matching actual lower witness. Boolean/float rationals,
wrong-source certificates, omitted/repeated pair IDs and failed unary mass
are explicit hostiles. No assertion is required for production guards.

The signed positive-norm bound, the factor maximizer and the global
one-parameter optimizer have separate control paths. The native source
trajectories remain the authenticated inherited input; the additionally
supplied finite pulse paths are independently replayed, not discovered by
this consumer. Neither is a new empirical convergence result.

The normal frozen receipt has19,434 body checks and LF SHA256
`2d17a1124daa2c4859c2a3bee039da8ad728e667cac918b5e93b167113fdac2d`.
The verification command adds one saved-equality gate. Final independent
proof/API and optimized comparison status will be recorded below.
