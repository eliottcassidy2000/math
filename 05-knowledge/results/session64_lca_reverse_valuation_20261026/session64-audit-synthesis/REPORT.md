# Independent synthesis review: exact guards, local threshold, and separate carriers

2026-10-26. **ACCEPTED at the pinned final bytes.** The Session64 root
`SYNTHESIS.md`, `CONNECTIONS.md`, and `newton-lca/` proof packet accurately
combine the separately audited results. No unresolved mathematical correction
is required. This review checks the synthesis's implications and source scope;
it does not claim another independent execution of every long native checker.
The small root LCA checker was replayed normally and under optimized Python.

## 1. The sharp endpoint-four claim is local and fully quantified

The accepted statement is: for every n>=8, every directed, closed, blue marked
Hamiltonian path P with endpoint outdegrees (n-5,4) has three successive
directed black K moves to its exact labelled K_P target, without a clone,
module, strongness, port-disjointness, or isomorphism filter.

This conclusion combines separate, compatible parts:

| Orders | Full marked-source coverage | New exceptional-source repairs |
|---|---|---:|
| 8 | 40 nonempty cells, dimension 4 each | 16 |
| 9 | 80 nonempty cells, dimension 7 each | 5 |
| 10 | 140 nonempty cells, dimension 10 each | 9 |
| 11 | 224 nonempty cells, dimension 14 each | 7 |
| 12 | 336 nonempty cells, dimension 18 each | 16 |
| 13 | 480 nonempty cells, dimension 23 each | 0 |
| 14 | 660 nonempty cells, dimension 28 each | 4 |
| >=15 | Prior audited unfiltered theorem | Inherited structural ledgers |

The lower-boundary audit independently reconstructs every source cell and its
complete guarded leaf partition. At n=8 the two exceptions are **illegal
old leaves**, not blue subloci: both entire dimension-three leaves are
exhausted. The n=9/n=10 exceptions are exactly five/nine affine singletons.
At n=11 all seven actual blue-locus sources are included. At n=12/n=14 the
six old failing cells have exactly 8+8+1+1+1+1 offending-blue completions.
Uniform old ledgers handle each complementary source family. Therefore the
root does not promote a successful finite search menu to a universal theorem.

The inherited `PORT_FRONTIER.md` was inspected at its exact quantifiers:
its n>=15 base partition plus complementary h-families has the same unfiltered
source and literal three-move target. Finite audits explicitly recheck the
lower orders rather than extending those formulas below their proof bounds.

The n=7 source `(6,4,8,19,103,79,15)` has target
`(122,124,121,112,32,64,16)`. Exhausting all 5040 target permutations gives
nine Hamiltonian paths, all blue. Own-K invariance of a word's color means
that no black move can enter this target. Thus 8 is the optimal uniform
integer threshold for **this closed marked endpoint-four replacement
statement**, not for whole black-core connectivity. The root explicitly
retains the latter as open, together with open marked paths and other layers.

Literal clone claims are also correctly separated: the original twenty chains
are clone-free at all eighty states, but that is an observation, not a source
filter. The n=11 repair forms the ordered pair (5,4) after its second move and
splits it at the target. The lower audit additionally checks all clone blocks
by exhausting vertex subsets at every literal state; no quotient lifting is
used to prove any new K move.

## 2. The bridge gate is sufficient, not an existence theorem

Write P=(0,...,z), B=K_P(T), k=d_T(z), I=N_T^-(0) minus {z}, and choose a
1-to-z Hamiltonian word W of B minus vertex 0. The cyclic words
`(1,...,z,0), (W,0), (0,W)` have the correct joins and exact path-edge XOR.
The degree condition d_T(1)!=k certifies the first two as black. For the
third, `g=r_(0,W) r_P` is an automorphism of B if and only if that word is
blue. Such g fixes 0,z and must preserve I union {1}; tournament
automorphisms have odd order. Consequently nonpreservation of the actual
labelled set or even order of g certifies the third word black.

The independent bridge proof supports all these implications for every size.
The root does not infer that a suitable W always exists, that the degree
guard is necessary, or that merely preserving the port-set size is enough.
Its equality hostile is correct: the n=11 transient-clone source has
`d_T(1)=k=4` and nevertheless has valid three-black cyclic bridges.

## 3. Counting: structural contraction and arithmetic filtering remain distinct

The finite-DAG potential V is a minimum suffix valuation, B=V(initial), and
`J'=J+V(w)-V(v)>=0`. The gauge `Z(v)=2^(V(v)-B)Y(v)` is integral on every
reachable prefix contribution, even where the displayed coordinate exponent
is negative. Being state-defined, it commutes with actual weighted mergers.
This is not integrality on arbitrary unrelated integer assignments.

For forward reduced minimum P, a minimum suffix has reduced cost zero.
Hence an edge belongs to some complete path of reduced cost <K exactly when
`P(v)+J'<K`. This is the exact support of potentially visible path
contributions, not a promise that their sum is nonzero. It removes modular
updates, not graph preprocessing or structural states.

The separate contraction deletes only unique unit-factor, zero-exponent
unaffordable-stage/early-completion skips. Deleting and uniquely reinserting
those edges gives a complete weighted-path bijection. V descends pointwise;
P and accumulated Y/Z agree only at surviving aggregates, not at every removed
history. The N=8 active `(3,3,(0,2))` has two choices and k=0 cost 3, so it is
correctly excluded from the contraction rule. Active stage and defect data
are retained.

The N=100 figures agree in producer and independent receipts: input states
13,444 -> 10,547, edges 50,044 -> 47,147, and K=2 modular updates 568 -> 64.
Input-state totals exclude the common terminal. Root wording preserves the
cost boundary: no measured runtime, peak-memory or asymptotic speedup, and
no newly proved equivalence of active defect states.

## 4. Newton/LCA: exact full image, scoped operations, signed core lane

For nodes p,p^2,... and q_j(X)=product_(i=1)^j(X-p^i), evaluations vanish
at earlier nodes and have exact valuation tau_j=j(j+1)/2 at every later
node. This is the sum of the relevant low-digit-trie LCA depths. The monic
Newton basis is integral-unimodular, and rowwise induction gives the full
kernel, not merely sufficient kernel elements. Thus the full polynomial
observation image has cyclic factors `Z/(p^(K-tau_j))` for j<S and tau_j<K.
For degree <h additionally require j<h. A sparse exponent lattice or one
reachable counting fibre need not attain this full image.

At K=64 the positive precisions are
`64,63,61,58,54,49,43,36,28,19,9`: eleven coordinates, 484 residue bits,
versus the previous 64 full 64-bit coordinates (4096 bits). These are
mathematical payload sizes, not Python-memory or runtime figures.

The native signed clone-excess transition reads only retained positive source
nodes, multiplies by integers and nonnegative powers, and therefore preserves
the exact kernel. The separate terminal observation is covered by the native
proof; no unavailable P(1), division by 2, or division by X is inferred from
the quotient. Semantic modular-zero support remains present. Ordinary
A000568 was already scalar; the Newton gain concerns the signed core bank and
absolute residues, not normalized digits after truncation.

The fresh independent literal census supports the newly added hostile:
naive ordinary-series substitution A(z/(1+z)) gives 0 versus actual 2 core
classes at n=4 and 20 versus 26 at n=6. It cannot replace the cycle inventory.
The root does not conflate ordinary tournament counts with pair-free counts.

The root LCA checker was replayed in both Python modes and exactly matches
its sealed receipts: 588 diagonal/LCA checks, 58,968 kernel-transition
evaluations, 72 exhaustive coefficient-image cases. This finite check is
separate from the independent native Newton audit through n=48.

## 5. Dynamics: coding metric, integer ancestry, resource, and clock

The valuation shift is on the positive-integer alphabet over the two correctly
punctured odd 2-adic domains. The same-code map preserves the metric formula
`v2(x-y)=A+min(p,q)` but not positive integrality. Constant words identify
ordinary 1 with twisted 3/5 and ordinary -1 with twisted 1. Collatz integer
sibling mergers such as 3,13 -> 5 therefore do not transport into the twisted
integer locus, whose nontrivial components are unbranched escaping rays.
Common low-bit prefixes and common forward futures remain different predicates.

The debt monoid retains both total increment and minimum prefix. Its
`R(w)=1-min D_j` uses empty and final prefixes, expressing positivity at the
returned odd section; execution of a fixed number of blocks omits that final
continuation guard explicitly. The (1,2)/(2,1) hostile correctly shows why
net drift alone is insufficient.

The finite control has 10,000 odd returns, 19,939 parameter macros, final
resource 9,884, minimum 6 at odd sections and minimum 5 over macros from the
first odd section. Blank setup starts at 0. These observation cuts are
consistent, and no macro count is presented as a native Turing-machine clock.
The native bridge is cited rather than rebuilt; the invalid printed boundary
clock is excluded. Haar ruin probability q^a is expressly a random-code
statement, not a probability or survival proof for fixed integer 27.

## 6. Corrections, provenance, and final disposition

Two small clarity changes requested during this review are present in the
pinned root bytes: the endpoint theorem now spells out its complete marked
path/degree quantifier; the normalized triangular matrix now has **units on
its diagonal**, not an ambiguous “unit diagonal” that could mean all ones.
For example, q_2(8)/2^3=3. Neither changes the theorem or implementation.
The matching producer Newton wording was clarified by its independent auditor.

`SOURCE_PINS.sha256` records the exact reviewed sources and sealed proof
packets. `verification.json` records checksum-manifest validation, matched
normal/optimized receipts and the local LCA replay. `SHA256SUMS` seals this
review's own artifacts. Workspace source prefixes in the pins are provenance
locators; publication may relocate files without changing the pinned bytes.

**Final disposition: PASS.** No outstanding mathematical correction,
permission request, denied action, uncertain execution, or audit blocker.
Only this scratch audit directory was written; no repository edit or push was
performed by this review.
