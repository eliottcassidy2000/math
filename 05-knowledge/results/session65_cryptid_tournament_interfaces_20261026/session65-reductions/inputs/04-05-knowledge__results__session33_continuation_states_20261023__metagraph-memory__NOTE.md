# Exact memory of fixed-path complement under node observation

Status: PROVED finite-involution statements, specialized to the actual tournament tilings; FINITE-EXACT independent replay n=3,4,5 and color witness n=6. No repository mutation or publication. No claim about LRC closure or deletion/extension operations.

## Inheritance and novelty boundary

Inherited: THM-796-three-sorted-recursive-tiling-line-node-incidence.md defines the tiling/node incidence and already diagnoses non-Markov node compression. THM-840-hamming-five-continuation-congruence-boundary.md supplies the operation-kernel criterion and minimum linear sidecar rank. Section 2 / Proposition 0 of 05-knowledge/results/boundary_support_20261020_metagraph.md already proves the Gram identity and rank r-c-b. Those are not new here. The Halmos block construction below is classical scaffolding, not new operator theory; see P. L. Robinson, [Julia operators and Halmos dilations](https://arxiv.org/abs/1803.09329), 2018. Targeted searches of theorem/results Markdown found no existing Halmos, 2r-c-b, involution-dilation, or BD^k formula in this snapshot. This is not a claim of literature priority.

Closest mechanism: continuation congruence. Canonical hostile: order-four node compression. Corrected near miss: a linear dilation is not a positive labelled tournament realization. Least-used sidecar: the complement target of an actual tiling. Live board: continuation congruence; Gram defect; deterministic states; positive realizability; color; literal involution versus rerandomization. Anchor: exact lost memory; niche: minimum observer dimension; wildcard: color obstruction.

## 1. Actual source and deterministic quotient

For n>=3, X is all fixed-path tournament tilings: vertices 0,...,n-1, forced arcs i->i+1, and one bit for every lexicographically ordered free edge (i,j), j>i+1. In the checker a one means i->j, a zero j->i. The actual operation tau toggles every free bit and preserves forced path arcs. It is NOT full tournament converse and NOT a single-arc reversal. Let o(x) be the tournament isomorphism class after merging converse classes. Let sigma reflect the free-edge array through (i,j)->(n-1-j,n-1-i). A line is blue iff sigma(x)=x, otherwise black; this predicate is tau-invariant.

For the operation alphabet consisting solely of tau and the output o, the coarsest deterministic continuation-complete quotient is

    q(x)=(o(x),o(tau x)),       q(tau x)=swap(q(x)).

Proof: its update is well defined. Any adequate state must recover both current output and the output after one tau, so equality of adequate states implies equality of q. All later outputs repeat these two because tau²=id. Thus q is minimal in the refinement sense, not merely minimal cardinality. Its number of states is |{(i,j): A_ij>0}|. Parallel lines and distinct presentations with the same endpoint pair disappear; even diagonal pairs may represent two or more tilings because tau has no fixed points.

If color is also an output, the minimal state is (o(x),o(tau x),color(x)), with update swapping endpoints and preserving color. It is strictly finer exactly at endpoint pairs supporting both colors. This statement concerns tau only; adding endpoint deletion, extension, or other operations can require further refinement, as THM-796 already demonstrates.

## 2. Linear observer and its exact minimum dimension

Let H=R^X, P the symmetric permutation of tau, J the fibre indicator matrix, F=J*J=diag(f_i), U=JF^(-1/2), V=im U and K=U*PU. Let r be the number of nonempty fibres. For a source vector h, observe y=U*h. Its entire future is determined by y and U*Ph. Hence the unobservable space is

    ker U* intersect ker U*P = (V+PV)^perp.

The minimal linear realization, for all signed source vectors, therefore has dimension

    dim(V+PV)=r+rank((I-UU*)PU)
             =r+rank(I-K²)=2r-c-b.

Here c is the number of support components and b the number of bipartite components, with loops excluding bipartiteness. The final equality is the inherited defect-rank theorem, not a fresh graph-rank result. Minimality follows because any linear state from which both outputs can be recovered has kernel contained in the displayed unobservable space. Sufficiency follows because V+PV is P-invariant. Equivalently the minimum additional linear sidecar rank is r-c-b, exactly THM-840 applied to current and next node observations.

This minimum is not the cardinality of the deterministic quotient. At n=5 it is 19, whereas the deterministic quotient has 38 states. Allowing signed linear superposition compresses further than a deterministic state partition. The image of nonnegative tiling masses in this minimum realization is a special polyhedral cone, not automatically a coordinate orthant.

### 2a. Positive/stochastic realizations have the larger deterministic minimum

This strengthens the positive-realizability boundary under a precise initialization requirement. Suppose a fixed h-state stochastic system has column-stochastic transition T and node-output matrix O, and every pure tiling x may be encoded by some probability vector p_x, with

    O p_x=e_(o(x)),       O T p_x=e_(o(tau x)).

All-times reproduction may also be required, but these first two times already suffice for the lower bound. For q(x)=(i,j), both displayed outputs are extreme points of their probability simplex. Thus for every hidden state a with (p_x)_a>0, necessarily

    O e_a=e_i,            O T e_a=e_j.

One hidden state cannot satisfy these identities for two distinct ordered pairs. Each realized pair needs a nonempty set of hidden states, and these sets are disjoint. Consequently h is at least the number of realized ordered pairs. This remains true even if x->p_x is nonlinear and even if no full-state intertwining T p_x=p_(tau x) is required.

The pair-state system attains the bound: encode x by the point mass at q(x), transition by swapping the pair, and output its first endpoint. It reproduces every time exactly. Hence the minimum number of positive stochastic hidden states is exactly |supp A|, the deterministic minimum, not the smaller signed linear dimension. At n=5 this proves an exact 38-versus-19 separation. Adding color as an output gives the analogous colored-pair minimum by extremality at times zero and one.

Critical scope: every pure tiling initial condition must be represented. If only uniform node-fibre initial conditions are required, the time-one output is generally mixed; the extremality step fails and this lower bound is not established. The result is about a positive hidden-state model, not construction of actual tournaments from an arbitrary matrix.

## 3. Minimum selfadjoint involution compression realization

Write Delta=(I-K²)^(1/2) and E=im Delta. Functional calculus makes E invariant under K and Delta commute with K. On R^r direct-sum E define

    L = [[K, Delta|E], [Delta, -K|E]].

Then L*=L and L²=I; its compression to the first summand is K. Its dimension is r+rank Delta. Conversely every selfadjoint involution with upper-left block K has off-diagonal block B satisfying BB*=I-K², so its auxiliary space has dimension at least rank Delta. This is the defect-space restriction of the classical Halmos construction.

The actual restriction P|(V+PV) is unitarily equivalent to this realization by an equivalence fixing the observed copy R^r: polar decomposition of (I-UU*)PU identifies its range with E, and the involution block equations identify the lower-right restriction with -K|E. No arbitrary invisible complement is part of this minimum realization.

Terminology warning: this is a one-step compression realization by an involution, NOT a power dilation with compression L^t=K^t. Indeed its even compressed powers are I. Nor does this construction produce a positive labelled fixed-point-free permutation, integer fibre sizes, blue/black structure, or actual tournament representatives. Those are separate realization constraints.

## 4. Exact memory kernel and initial-condition boundary

Relative to H=V direct-sum V^perp, write

    P=[[K,B],[B*,D]].

The involution equations give BB*=I-K² and BD=-KB. Therefore for every integer k>=0,

    BD^kB*=(-1)^k K^k(I-K²).

For source evolution h_(t+1)=Ph_t, with coordinates (y_t,z_t), elimination of the hidden variable gives, for t>=0,

    y_(t+1)=K y_t + sum_(s=0)^(t-1) BD^(t-1-s)B* y_s + BD^t z_0.

The sum is empty at t=0. Every convolution kernel is therefore determined by K. The initial forcing is

    BD^t z_0=(-K)^t eta,   eta=Bz_0=y_1-Ky_0.

K alone does NOT determine eta for arbitrary hidden initial data. The safe claim is: K determines the full projected impulse response for fibre-uniform initial data z_0=0, all memory-kernel coefficients, and the forcing propagation rule; a nonuniform initial state still needs its initial sidecar eta. This precisely marks what "memory completely determined by K" can mean.

Directly, U*P^tU alternates I,K. In contrast, K^t is the process with an intermediate projection UU* between successive literal complements: each step redistributes mass uniformly in the current node fibre. For probability row kernels this is T=F^(-1)A, similar to K. Actual two-step node evolution from a uniform fibre is I, whereas the rerandomized process is T². Equality for all times holds iff I-K²=0, equivalently every fibre maps wholesale to one fibre.

## 5. Exact actual-source replay

Reproduce from the enclosing workspace:

    python3 metagraph-memory/check.py

The script imports no repository producer or data. It exhausts all 2^binom(n-1,2) masks, canonicalizes by all degree-compatible vertex permutations for each tournament and its converse, and constructs J,P,A directly. Sorting outdegrees first is a canonical relabelling invariant; within each equal-degree block all permutations are exhausted, so no isomorphisms are discarded. Integer labels are sorted canonical codes in this checker only, not the repository's node IDs.

All ranks and matrix identities use exact SymPy rationals, no floating thresholds. Checks use explicit require/raise, so remain active under python3 -O; normal and optimized runs are compared. For n=3,4,5 it checks Gram equality, direct rank [J,PJ], support-component rank formula, memory kernels k=0,...,4, projected times t=0,...,6, and independently iterated deterministic partition refinement. The n=3 zero-defect case is a positive control; n=4 and n=5 are strict Markov-failure controls. At n=6 it does only the inexpensive actual-source pair/color census, not matrix ranks.

| n | tilings | nodes r | pair states | colored pair states | linear dimension | defect rank |
|---|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 2 | 2 | 2 | 0 |
| 4 | 8 | 3 | 5 | 5 | 5 | 2 |
| 5 | 64 | 10 | 38 | 38 | 19 | 9 |
| 6 | 1024 | 34 | 356 | 364 | not replayed | not replayed |

At n=6 four unordered endpoint supports carry both colors (eight ordered pairs). An explicit witness in the stated bit gauge is x=521 (blue) versus x=27 (black), both with q(x)=(12,14). Their complements are 502 and 996. The node IDs are local, so the masks, gauge, and canonicalizer are the reproducible witness. All witnesses are in results.json. This matches the inherited THM-796 count of 187 colored versus 183 plain unordered supports, but was independently computed.

## 6. Preserved/lost-coordinate ledger

- Actual tiling -> ordered endpoint pair: preserves every tau-word node output; loses presentation, parallel-line identity, color at mixed supports, and recursive faces. Sidecar for colored tau output: one color bit only where needed.
- Signed tiling vectors -> minimum linear observer state: preserves all linear node outputs under tau; loses the full source vector and positive labelled realization. Sidecar over current node output: r-c-b independent linear coordinates.
- K -> minimum involution block model: preserves node impulse response and all memory kernels; loses nonuniform initial forcing unless supplied, literal source labels, fibre integers, colors and tournament realizability.
- K -> K^t: changes the source operation by inserting fibre rerandomization. It is not an exact literal-complement continuation, except at zero defect.

The new scoped consequence is an exact three-way distinction between minimal deterministic continuation, minimal linear memory, and positive labelled realization for the actual fixed-path complement. No all-order metagraph connectivity, blue classification, or LRC theorem is asserted.
