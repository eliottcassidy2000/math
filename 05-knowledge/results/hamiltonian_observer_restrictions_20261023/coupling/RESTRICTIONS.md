# Hamiltonian observers: native obstruction and auxiliary carrier

Status: **PROVED** for the finite-partition criteria, the all-size obstruction,
and the auxiliary carrier below; **FINITE-EXACT** for the enumerated controls.
The proofs are included. No claim here identifies an auxiliary move with a
native tournament metagraph edge or solves its all-size wiring problem.

## Inheritance and live concepts

Closest proved mechanism: the full-word-cube graph/tournament coupling has
gap exactly `8/9` in order at least three; see
[`session42.../coupling/COUPLING_OPERATOR.md`](../../session42_shared_lattices_square_colors_20261023/coupling/COUPLING_OPERATOR.md)
and its signed-Walsh input
[`session42.../sign-moments/THEOREM.md`](../../session42_shared_lattices_square_colors_20261023/sign-moments/THEOREM.md).
The inherited theorem does not assert a native fixed-path restriction.

Canonical hostile: a marked first vertex makes graph degree and tournament
indegree the same shared observable. Corrected near miss: discarding a forced
Hamiltonian path's `n-1` bits leaves a full cube, but its *native* tournament
quotient is not the whole-cube tournament quotient on the remaining bits.
Least-used sidecar: the actual observer word, not just its two endpoints.
Compare the exact path/converse/integer-guard scope in
[`odd_cell_completion.../seven/CANONICAL.md`](../../odd_cell_completion_20261023/seven/CANONICAL.md).

The concept board is: conditional expectation; supplied Hamiltonian word;
native versus auxiliary quotient; retained marks; positive weighting; actual
return clocks. The last belongs to the separate
[`TRACE_WITH_CLOCKS.md`](../trace/TRACE_WITH_CLOCKS.md).

## 1. A criterion that remains valid after restriction

Let `Omega` be a finite positively weighted space, and let `A` and `B` be two
observations, equivalently two partitions. Write

\[
 W_{ab}=w\{x:A(x)=a, B(x)=b\},\quad
 r=W\mathbf1,\quad c=W^t\mathbf1,\quad Z=\sum_a r_a.
\]

Delete zero-mass rows and columns. The alternating chain on `A` is

\[
 P=D_r^{-1}WD_c^{-1}W^t.                                      \tag{1}
\]

It is stochastic, reversible with stationary probability `r/Z`, and positive
semidefinite. Its Dirichlet form has the useful conditional-variance identity

\[
 \langle f,(I-P)f\rangle_{r/Z}
       =\mathbb E[\operatorname{Var}(f(A)\mid B)].              \tag{2}
\]

**Exact restricted criterion.** A Poincare gap at least `gamma` holds if and
only if

\[
 (1-\gamma)D_r+\frac\gamma Zrr^t-WD_c^{-1}W^t\succeq0.         \tag{3}
\]

In particular the `8/9` bound is equivalent to the centered normalized
incidence matrix

\[
 D_r^{-1/2}WD_c^{-1/2}-\frac{\sqrt r\sqrt c^{\,t}}Z
\]

having operator norm at most `1/3`. These are finite linear-algebra
equivalences, with no transitivity, cube, or group-action hypothesis.
Equation (1) is the composition of two conditional expectations; expanding
its quadratic form proves (2), and subtracting `gamma Var(f)` proves (3).

**Shared-observable obstruction.** If a nonconstant `J` is measurable from
both `A` and `B`, then `PJ=J`, so the global Poincare gap is zero. Equivalently,
the bipartite support of `W` is disconnected. Each connected component can
be studied with its own constants and mass; deleting the shared observable
from notation does not make it disappear.

**Coarsening is safe on an unchanged probability space.** Replacing either
observation by a function of it cannot increase maximal correlation, because
the normalized centered test-function spaces shrink. Thus any coarsening of
both whole-cube observations retains the `8/9` bound. This applies to the
inherited complement/converse coarsening. It does not justify conditioning,
refining observations, or substituting a different target map.

**Strictly positive weights have a quantitative repair.** Suppose the
unweighted/base version has gap `gamma`, and multiply its point masses by
`w(x)` with `0<a<=w(x)<=b`. With the same observations, the new coupling has
gap at least

\[
                 \gamma\frac ab.                            \tag{4}
\]

For a direct proof use
`E_w Var_w(f|B) = inf_h E_w[(f-h(B))^2]`. If `Z_w=E_0 w`, this is at least
`(a/Z_w) E_0 Var_0(f|B)`. Also
`Var_w(f)<= (b/Z_w) Var_0(f)`. Combining these with the base Poincare
inequality proves (4). The `h` infimum includes every function on the actual
`B`-cells, so no invariance of `w` is required. Hard deletion has `a=0` and
is deliberately outside this positive conclusion.

## 2. Supplied Hamiltonian-path compression and the native failure

Fix the actual path `P=(v_0,...,v_{n-1})` and force all its consecutive arcs
forward. A word `x` has one free bit for every pair `(i,j)` with `j>=i+2`;
bit one means `v_j -> v_i`. There are

\[
         2^{\binom{n-1}{2}}
\]

words. Encode the free bit `(i,j)` as the undirected edge `{i,j-1}` in a graph
`G_x` on `0,...,n-2`. Let `T_x` be the *native* `n`-vertex tournament obtained
by restoring the forced path arcs. The two observations being tested are

\[
             A(x)=[G_x],\qquad B_{\rm native}(x)=[T_x].        \tag{5}
\]

**All-size native obstruction.** For every `n>=3`, (5) has gap zero.

Proof. `G_x` is empty exactly for the zero word. A transitive tournament has
a unique Hamiltonian vertex order. Since our order `P` is Hamiltonian, a
transitive `T_x` must have every arc forward in that very order; hence `x=0`.
Thus the empty graph and transitive native tournament form an isolated
incidence component, and the domain contains another word. The indicator of
this component is a nonconstant shared observable. Converse-merging native
tournament classes leaves this component isolated. Refining either
observation by more marks also leaves it isolated. QED.

This argument is about ordinary graph isomorphism at the source. Pairing
source graphs by complementation changes the observation and can join the
singleton to another component. That change is tested separately below; it
does not recover `8/9` even in order four.

### The four-vertex hostile survives the obvious repairs

For `n=4`, the four source graph types on three vertices are empty,
one-edge, two-edge path, triangle. On the seven nonzero words their masses
are `(3,3,1)`. The native coupling on these three graph types is exactly

\[
 P_* =
 \begin{pmatrix}
 11/15&1/5&1/15\\
 1/5&3/5&1/5\\
 1/5&3/5&1/5
 \end{pmatrix}.                                               \tag{6}
\]

Its spectrum is `{1,8/15,0}`: the centered eigenvector `(4,-3,-3)` has
eigenvalue `8/15`; trace and rank supply the remaining value. The exact
nontransitive-component gap is therefore `7/15`, rather than `8/9`.
Converse-merging the two diamond targets does not change (6).

Even combining compressed source complement-pairs and native target
converse-pairs yields a connected two-source-state chain whose exact gap is
only `8/15`. Thus connectivity, deleting transitivity, and these paired
symmetries are each insufficient to recover the whole-cube bound.

The finite counts of source components for `n=3,4,5,6` are respectively
`(1,1)`, `(1,3)`, `(1,10)`, `(1,33)`. The subsequent constructive proof in
[`NATIVE_CONNECTEDNESS.md`](NATIVE_CONNECTEDNESS.md) establishes exactly two
components at every order: a guarded strong-four descent, with separate
matching/star/complete repairs, connects every nonempty source to the
one-edge class. It also supplies an explicit native gap bound depending on
`n`, without recovering `8/9`.

## 3. The positive carrier: use the right tournament quotient

There is another tournament attached to exactly the same free word: on the
`m=n-1` compressed positions, orient `{i,k}` backward precisely when it is an
edge of `G_x`. Call it `U_x`, and use

\[
          A(x)=[G_x],\qquad B_{\rm aux}(x)=[U_x].              \tag{7}
\]

**Auxiliary Hamiltonian carrier theorem.** The coupling (7) has gap exactly
`8/9` for `n>=4` (and gap one for `n=3`). Its self-adjoint word-space
realization is `Q=E_A E_B E_A`, which has the same nonzero spectrum as the
source chain. In particular it has the same gap and is PSD.

Proof. The compression is a bijection of our supplied-path words with the
entire `m`-vertex edge cube, carrying uniform measure to uniform measure.
The observations in (7) are exactly the two observations of the inherited
whole-cube theorem. Its bound and sharpness transfer through this bijection.
The word-space statement follows by restriction to the range of `E_A`;
its orthogonal complement is killed by the rightmost `E_A`. QED.

This is a real positive statement about actual path-marked native states:
every intermediate cube word still reconstructs a native tournament
containing the same supplied path and the same vertex positions. But the
intermediate quotient `U_x` has `n-1` vertices, and must not be relabelled as
the native tournament `T_x`. The operator may change native tournament
isomorphism class, other marked port incidences, and electrical mod-8 data.
It is not asserted to be a walk along legal `K_P` metagraph moves.

A useful local variant is to fix a transitive prefix, orient every arc from
that prefix to a free suffix forward, and run the same construction on the
suffix. The fixed prefix and every one of its incidences are then retained
exactly. A suffix of at least four vertices carries the same `8/9` operator.
This is a sufficient free-room condition, not a theorem that arbitrary
attachments admit such a suffix.

## 4. Actual marks can carry a conserved quantity

There is a separate, very small full-cube marking obstruction. Keep the
first observer vertex `0` fixed under both ordinary graph and affine
tournament relabellings. Its graph degree is precisely its tournament
indegree, because every incident backward bit points toward vertex `0`.
Thus the degree `0,...,n-1` is a common observable and the globally marked
coupling has at least `n` components.

The source/sink boundary strata `degree=0` and `degree=n-1` are positive
controls: deleting the marked root identifies each with the unrestricted
cube on `n-1` vertices, so each has gap `8/9` for `n>=4` and gap one for
`n=3`. This illustrates why a retained owner mark must be examined before
invoking a mixing theorem, and why a component-specific theorem can survive.

An actual electrical mod-8 receipt is a separate retained coordinate. No
test in this packet substitutes a bit-count residue or Hamiltonian-path
count residue for it. Arbitrary hard guards are handled by the return-chain
construction in [`TRACE_WITH_CLOCKS.md`](../trace/TRACE_WITH_CLOCKS.md), with
all the excursion clocks and carrier restrictions stated there.

## 5. Connection record and reproduction

| Item | Exact meaning |
|---|---|
| Source | Uniform free bits of a supplied Hamiltonian `n`-vertex word |
| Target | Graph/tournament observations on the compressed `n-1`-vertex cube |
| Map | Native free chord `(i,j)` becomes compressed pair `(i,j-1)` |
| Preserved | Every bit, supplied path, named positions, measure, auxiliary operator spectrum |
| Destroyed by a quotient | Particular word, native source ownership, arbitrary endpoint incidence |
| Necessary sidecar | Actual word/owner lift, port records, exact guards, and return clocks when tracing |
| Minimal decisive hostile | Native `n=3` is already disconnected; native nontransitive `n=4` has gap `7/15` |
| Positive repair | Auxiliary quotient; safe coarsening; positive weights; free suffix; authenticated return chain |

`check.py` imports only the earlier exact canonization routines, enumerates
every supplied-path word for native `n=3..6`, and uses rational LDL for all
PSD and gap claims. It verifies the exact eigenvector of (6), source pairing,
the auxiliary sharp gap, complement/converse coarsening, and non-invariant
positive weights `w(x)=1+(popcount(x) mod 3)` with the guaranteed `8/27` gap.
Separate controls enumerate all first-root-marked full-cube words for
`n=3..5`. There are no `assert` statements and no floating spectral tests.

```text
python -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/coupling/check.py
python -O -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/coupling/check.py
```

The nontransitive native incidence component is now connected by the
separate constructive proof. Sharper bounds, and a descent that retains a
fixed native source and its actual indexed ports throughout, remain open.
Neither whole-cube Fourier diagonalization nor an auxiliary return-chain
gap supplies those stronger native statements for free.
