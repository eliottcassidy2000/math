# Receipt-policy geometry, adaptive termination, and refinement

**Status:** PROVED elementary criteria; FINITE-EXACT authenticated controls.
**Scope:** finite stored proof graphs and precisely defined refinements.
**Date:** started 2026-10-09; completed 2026-10-11.

The simplex/cube/sphere representations describe how to choose among already
authenticated proof steps. They do not produce a missing step or discharge an
open terminal. The additional result here is a common-rank certificate that
remains valid under adaptive choices and arbitrary hidden-state refinement,
provided an exact projected-transition condition is retained.

## 1. Inheritance and the connection contract

The closest proved mechanism is
[finite receipt-network grounding](paper_local_global_grounding_20261007.md):
an authenticated common-future graph has a finite ROOT path exactly when its
source admits the corresponding finite signed flow or positive effective
capacity. That package already exports actual ROOT words. We reuse its exact
receipt verifier, rational linear solver, and path exporter.

[The Poisson source dual](collatz_poisson_source_dual_20261005.md) and
[Poisson patch calculus](collatz_poisson_patch_20261005.md) already distinguish
resolvent existence, positive support, and a source-specific certificate;
they also retain boundary obligations under elimination. The
[atomic prefix ledger](atomic_prefix_certificate_20261005.md) counts only
independently verified source certificates. None of these is replaced or
promoted to universal coverage here.

The live objects are authenticated receipt edges, row probability simplexes,
a common rank, the terminal fiber, and hidden-state refinement. The map from
receipts to a matrix preserves their labelled availability but forgets actual
valuation words and seam clocks. Those receipts remain mandatory sidecars.
The map from probabilities to square roots preserves probabilities exactly;
the shape alone proves neither a transition nor a ROOT claim.

Two cheap hostiles organize the result: a finite closed component containing
27 and 31, and a time-dependent positive exit probability whose total chance
of eventual exit is only one half. Infinite escape is a separate boundary.

## 2. Which geometric object?

At a proof state with \(d\) available outgoing receipts, a randomized choice
is \(p\in\Delta_{d-1}\). Independent choices at all states give
\(\prod_i\Delta_{d_i-1}\). This is a cube when every nontrivial row has exactly
two options; a row with three options contributes a triangle instead.

The lossless substitution \(z_{ij}=\sqrt{p_{ij}}\) gives
\[
z_{ij}\ge0,\qquad \sum_j z_{ij}^2=1.
\]
Thus the same policy is represented by a product of positive orthants of
spheres, not generally by one sphere. If coordinates are rescaled by positive
weights \(w_{ij}\), the corresponding equation becomes
\(\sum_j z_{ij}^2/w_{ij}=1\), an ellipsoid in each row.

A uniform floor \(p_{ij}\ge\epsilon>0\) on required edges keeps policies
away from the trapping boundary. In the binary case, independently restricting
the two probabilities to \([\epsilon,1-\epsilon]\) gives a smaller product
box. Merely lying in the open simplex at every time does not supply one
uniform floor.

## 3. Fixed finite policy: resolvent and deterministic extraction

Let \(V\) be finite, with ROOT denoted by 1. A nonroot row sums to one;
ROOT is killed. Every positive matrix arc must have an authenticated receipt
in its indicated direction. A common-future receipt can be used in either
direction because it is a proof implication; this does not assert that
backward proof motion is a physical Collatz step. An identity self-loop is a
tautological proof step and adds no grounding.

Write \(Q\) for the nonroot submatrix. For a fixed policy the following are
equivalent:

1. Every vertex has a positive-support path to ROOT.
2. \(\rho(Q)<1\).
3. \(h=(I-Q)^{-1}{\bf1}\) exists, is nonnegative, and satisfies
   \(h-Qh={\bf1}\).

For 1 implies 2, select one simple ROOT path from every vertex. There are
finite \(d\) and \(\delta>0\) such that ROOT is reached within \(d\) steps
with probability at least \(\delta\), from every vertex. Therefore
\(\|Q^d\|_\infty\le1-\delta\). The convergent Neumann sum proves 3.
Conversely, if ROOT is unreachable from a vertex, its reachable nonroot
support is closed. At a minimum of \(h\) on that finite set, \(Qh\ge h\),
contradicting unit decrease. This also proves uniqueness and positivity.

More generally, any nonnegative vector with \(h_1=0\) and
\[
h_i-\sum_j p_{ij}h_j\ge1 \quad(i\ne1)                         \tag{1}
\]
is a sufficient rational certificate. Each nonroot \(h_i\ge1\), and
\[
\theta=\max_{i\ne1}\frac{\sum_jp_{ij}h_j}{h_i}<1.
\]
This is an exact spectral-radius upper certificate, with no numerical
eigenvalue computation.

Some supported successor has strictly smaller rank. Repeatedly choose one,
retaining its receipt: at most \(|V|-1\) choices reach ROOT. The inherited
seam-clock exporter then performs bounded replay and returns the strict
first-hit valuation word. It stops at ROOT rather than padding a root cycle.
Rank counts proof decisions, not odd Collatz steps, bit operations, or source
size. Increasing integer labels along the selected proof path is allowed.

For a source-specific query, only its supported reachable subgraph is needed.
A disconnected unresolved component elsewhere need not block that query.
Missing row mass is not silently treated as a terminal proof.

## 4. Actual Collatz controls and the terminal obligation

The stored square example has common-future receipts
\[
7\xrightarrow{(1,1,2,3)}5\xleftarrow{(1)}3,\qquad
3\xrightarrow{(1)}5,\qquad 5\xrightarrow{(4)}1.
\]
Choose 7 to 3 with probability 1, and each of the two choices at 3 and 5
with probability one half. In order \((1,3,5,7)\),
\[
h=(0,8,5,9).
\]
The rank compiler selects \(7\to3\to5\to1\) and authenticates
the actual ROOT word \((1,1,2,3,4)\).

The parent-suggested control is
\[
31\xrightarrow{(1)}47\xleftarrow{(1,2,1)}27.
\]
With only these two proof directions, \(\{27,31\}\) is closed and no choice
of internal positive weights supplies a ROOT exit. The script separately
stores the explicit 41-letter, cost-70 first-hit ROOT word of 27. That finite
data was extracted by a declared 128-step bounded check during construction;
production only replays the frozen tuple. Adding this supplied receipt gives
the policy \(31\to27\), with 27 choosing 31 or ROOT equally. Its rank is
\((h_1,h_{27},h_{31})=(0,3,4)\). The exported ROOT word of 31 has 39 odd
steps and cost 67. The theorem does not assume the terminal word from its
mere label.

This is also the relevant boundary for a closed zero-defect Mersenne fan:
changing weights inside the stored component cannot create a missing exit.
The matrix quotient still owes an actual authenticated edge to a grounded
terminal. No giant fan is materialized or newly grounded here.

## 5. Adaptive policies: positivity fails, common rank repairs it

Consider one nonroot state with a self-loop and an already authenticated
ROOT exit. At step \(t\ge0\), use the exit with probability
\(p_t=1/(t+2)^2\). Every chosen row is strictly positive; freezing any one
of them gives a proper finite policy with expected time \(1/p_t\). Yet
\[
\Pr(\tau>N)=\prod_{t=0}^{N-1}(1-p_t)
           =\frac{N+2}{2(N+1)}\longrightarrow\frac12.           \tag{2}
\]
The first failed implication is replacing a fixed policy by an adaptive
sequence without a shared progress certificate.

There are two useful sufficient repairs.

**Uniform paths.** On a fixed finite support graph with ROOT paths of length
at most \(d\), suppose every edge of each selected path has conditional
probability at least \(\epsilon>0\), after every history. Then
\[
\Pr(\tau>t)\le(1-\epsilon^d)^{\lfloor t/d\rfloor}.              \tag{3}
\]
The conditional probability of following a selected path is at least
\(\epsilon^d\); iterate this bound on disjoint blocks. Independence of
successive policies is not required.

**Common rank.** Let \(\mathcal P_i\) be a finite row polytope. Check
\[
h_i\ge1+\max_{p\in\mathcal P_i}\sum_jp_jh_j,\qquad h_1=0.       \tag{4}
\]
Only the listed vertices need checking. Every history-dependent choice
inside these polytopes then satisfies the same conditional drift.
With \(X_t\) stopped at ROOT, telescoping (4) gives
\[
\sum_{t<T}\Pr(\tau>t)\le h(X_0),\qquad
\mathbb E\tau\le h(X_0).
\]
Taking the maximum ratio over all vertices also gives \(\theta<1\) and
\[
\Pr(\tau>t)\le \mathbb E h(X_t)\le h(X_0)\theta^t,              \tag{5}
\]
because rank is at least one outside ROOT. The bound may be capped at one.
Almost-sure termination does not mean every probability-zero infinite
realization is forbidden.

For the square graph with binary probabilities in \([1/4,3/4]\), a common
rank is
\[
(h_1,h_3,h_5,h_7)=(0,32,25,33),\qquad\theta=32/33.
\]
There are four product vertices. Formula (3) separately gives
\((63/64)^{\lfloor t/3\rfloor}\). These are policy-level guarantees on
authenticated edges, not a new integer-coverage theorem.

## 6. A source-specific condition that survives refinement

Here is the additional reusable transfer. Keep a finite coarse policy
system satisfying (4). Let \(Z\) be a finite or countable collection of
refined proof states, with projection \(\pi:Z\to V\). Require:

- the fiber of ROOT consists exactly of genuine completed proof states;
- at each nonterminal state and after each history, the conditional
  distribution of \(\pi(Z_{t+1})\) belongs to
  \(\mathcal P_{\pi(Z_t)}\);
- the concrete transition still retains its authenticated receipt or
  continuation sidecar.

Then \(H(z)=h(\pi(z))\) satisfies the same drift and contraction estimates.
Consequently (5) and \(\mathbb E\tau\le h(\pi(z_0))\) hold, even if forgetting
the hidden state produces a non-Markov process. No finite bound on the
number of hidden states is needed: \(H\) is bounded by the fixed finite
coarse maximum. The proof is conditional expectation followed by the same
telescoping and geometric estimates. It applies just to states reachable
from the requested source if that is the declared system.

The API verifies finite instances using supplied exact convex coefficients
for every pushed-forward row. It refuses a nonterminal falsely projected
to ROOT or a new trapping transition outside the old polytope. The abstract
refinement validator does not itself authenticate Collatz labels; the
separate ROOT compiler requires actual receipts.
For countable systems, the condition at every reachable state is a theorem
premise, not something certified by these finite examples. Effective proof
export also needs effective access to the authenticated transitions.

This condition is stronger than preserving an abstract graph shape. Inserting
arbitrarily slow hidden waits can destroy unit drift. Formula (2) is the
smallest temporal warning. Similarly, adding alternatives enlarges a
polytope and requires rechecking its new vertices. Restricting alternatives,
or refining only their hidden descriptions while preserving the conditional
pushforward, automatically preserves the certificate.

The source-specific rank is a checkable sufficient obligation. It is not an
independently obtained positive floor for an unresolved source. Compared with
the inherited Poisson dual, it controls an undiscounted adaptive proof
procedure; compared with the atomic ledger, it certifies the procedure
that produces a stored ROOT word rather than assigning mass before that
word is authenticated.

## 7. Why finiteness matters, even with positive transitions

On the abstract states \(0,1,2,\ldots\), let 0 be ROOT; from each \(n>0\),
move to \(n-1\) with probability \(1/3\) and to \(n+1\) with probability
\(2/3\). Every state has a finite ROOT path and all permitted transitions
are positive. Nevertheless
\[
\Pr_n(\tau_0<\infty)=2^{-n}.
\]
Indeed, on the finite interval \([0,N]\) the unique harmonic boundary
solution is
\[
\Pr_n(\tau_0<\tau_N)=
\frac{2^{-n}-2^{-N}}{1-2^{-N}}.
\]
Let \(N\to\infty\). Every finite path hitting 0 has a finite maximum, so the
increasing union of these events is exactly eventual absorption. The killed
operator on bounded functions has norm one for every power: start above
the time horizon. Thus its spectral radius is one.

This does not deny finite proof-path existence; it refutes the attempted
infinite-state upgrade from path existence to almost-sure random termination.
It is not a Collatz counterexample. A valid common rank of the refinement
form above prevents precisely this escape.

## 8. Reproducible scope

The companion script uses standard-library exact fractions and explicit
checks, including under optimization. It checks 512 directed support graphs
on three nonroot nodes, 81 strict-interior square policies, all four robust
vertices, 49 adaptive survival horizons, 192 hidden-mode refinements, the
actual 7 and 27/31 receipts, 256 temporal hostile prefixes, and 31 finite
gambler's-ruin intervals. Malformed sources, float probabilities, missing
terminal receipts, ROOT padding, false terminal fibers, missing row mass,
and false rank claims are rejected.

The only declared orbit discovery was the bounded extraction of the frozen
27 word. Production consists of supplied-word replay, exact matrix/rank
checks, path extraction, and clock-bounded export.

Run from the repository root:

    python -B 04-computation/experiments/collatz_receipt_flow_geometry_20261009.py
    python -B -O 04-computation/experiments/collatz_receipt_flow_geometry_20261009.py

Normal and optimized runs agree with the saved output: **2,206 checks**.
LF-normalized SHA256:
`0952123f84717e653543ddd273eea4e437ea05447ee65c4bcd466d44bfe36845`.

Next concrete target: furnish a coarse rank and exact projected-transition
contract for one currently open child port, retaining a genuinely grounded
terminal. Reweighting its closed old component or drawing it as a sphere
cannot substitute for that missing receipt.
