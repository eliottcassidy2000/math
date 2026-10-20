# Canonical chain growth: exact transition law and its counting boundary

**20 October 2026. Status: PROVED elementary identities below, using the
inherited canonical-chain theorem; FINITE-EXACT controls explicitly bounded.
No new closed recurrence for unlabelled tournament classes is claimed.**

This note is a targeted inheritance audit and an independent proof of the
vertex-extension mechanism. Its central positive result is an exact transition
kernel on `(number of vertices, number of canonical chains)`. That observable
really is sufficient for its stated labelled counting task. It does not
determine the unlabelled child, its automorphism group, or its line incidence.

## 1. Inheritance, hostile examples, and retained coordinates

The closest proved mechanism is
[THM-4463, homogeneous tournament contacts](../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md),
section 2. For a tournament `T`, form the graph `J(T)` whose edges are
homogeneous pairs: `{u,v}` is an edge exactly when

    T(u,w)=T(v,w) for every w outside {u,v}.

Orient those edges by the existing arcs of `T`. THM-4463 proves that the
components are consistently directed paths and that their vertex sets are
exactly the maximal transitive modules. Let their lengths be
`l_1,...,l_b`, with sum `n`. Their quotient is homogeneous-pair-free.
An isolated vertex is a chain of length one.

The proof matters: a vertex cannot have two homogeneous outgoing neighbors
or two homogeneous incoming neighbors; the third vertex would distinguish
one alleged pair. Directed cycles are impossible by propagation of
homogeneity. A path is transitive and every external vertex sees it
uniformly. Conversely, consecutive vertices of any transitive module are
homogeneous in the whole tournament. This proves canonical maximality.

Other directly relevant prior work is:

- [THM-4169, prime-parent augmentation](../../../01-canon/theorems/THM-4169-prime-parent-one-vertex-augmentation-and-quartic-johnson-transfer.md),
  sections 1 and 2: a prime parent of order `q>=3` has exactly `2q` clone
  attachments, two uniform attachments, and `2^q-2-2q` prime attachments.
  Its root-preserving orbit count explicitly does not count unrooted children.
  The clone count below extends the homogeneous-pair part to every parent,
  with `n+b` replacing `2n`.
- [THM-4131, independently audited augmentation through order eight](../../../01-canon/theorems/THM-4131-strong-tournament-centrality-through-order-eight.md).
  Its independent script
  `04-computation/tournament_strong_centrality_through_order_eight_thm4131_independent_audit.py`
  already enumerates every one-vertex attachment of every parent class and
  deduplicates canonical child codes. This is exact generation without
  Burnside, not a newly invented algorithm here.
- [THM-3235, transitive blowup and observer twins](../../../01-canon/theorems/THM-3235-transitive-tournament-blowup-saddle-and-hypergeometric-decimation.md)
  and [THM-3121, the walk-content kernel](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md)
  supply the existing substitution
  response theory. The preceding
  [growth note](../metagraph_growth_20261020/growth.md) gives the exact
  clone-run multinomial and a hostile showing that block words do not
  determine ordinary complement targets.
- [THM-796, three-sorted recursion](../../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md)
  proves why a merged tournament node is not generally a deletion or line
  update state. The labelled scalar kernel below is a different, precisely
  specified observable; it does not contradict that obstruction.

The canonical hostile to naive core compression is the transitive parent:
alternating the new vertex's relations along its one chain can create a
pair-free child of arbitrarily large order. Thus one old core vertex can
give birth to an arbitrarily large new core in one extension. The corrected
near miss is the suggestion that the new vertex may join two distinct old
chains: this is impossible, as proved next.

The concept board is canonical chains, attachment words, switch positions,
clone insertion, weighted quotient, rooted augmentation, and lost symmetry.
The anchor is exact extension; the niche is labelled counting; the wildcard
is alternating generation of the pair-free core bank. These are related by
literal maps below, rather than by an analogy to another open problem.

## 2. The edge update and the stronger clone dichotomy

Adjoin a new vertex `x` to `T`. Write `z_v=1` when `x->v`, and zero
otherwise. There are `2^n` attachment words. For each canonical chain use
its intrinsic order `v_1->...->v_l`. Let `s` be the total number of switches
`z_(v_i)!=z_(v_(i+1))` along all chains. Let `d` be the degree of `x` in
`J(T+x)`.

**Proposition 1.** The old homogeneous edges retained in `J(T+x)` are
exactly the old edges with equal endpoint bits. No old nonedge becomes a
homogeneous edge. The only new edges touch `x`, and `0<=d<=2`. Consequently,
if `b'` is the number of canonical chains of `T+x`, then

    b' = b+1+s-d.                                           (1)

**Proof.** Every old witness distinguishing a pair remains present after
adjoining `x`. An old homogeneous pair acquires exactly one additional test:
whether `x` sees both endpoints the same way. This proves the edge statement.
THM-4463 makes both homogeneous graphs forests of paths. The old forest has
`n-b` edges; the new forest has `n-b-s+d` edges on `n+1` vertices.
Components equal vertices minus edges, proving (1). QED.

There is a stronger restriction than just `d<=2`.

**Proposition 2, clone dichotomy.** If `d>0`, the new vertex is inserted
into a single existing canonical chain, every other chain is unchanged,
and the core quotient is unchanged. Exactly two cases occur:

    insertion at a chain endpoint:  (s,d)=(0,1);
    insertion inside a chain:       (s,d)=(1,2).              (2)

In both cases `b'=b`. If `d=0`, then `b'=b+1+s`. In particular canonical
core order never decreases under vertex extension.

**Proof.** Suppose `{x,u}` is homogeneous and `u=v_j` belongs to one old
chain. For every old vertex `w!=u`, homogeneity forces
`x->w` exactly when `u->w`. On every other old chain the new bits are
therefore constant and agree with the quotient relation of `u`. On the
chain of `u` they are zero before `u`, one after `u`, with the mutual bit
at `u` free. Thus the whole chain word is `0^a 1^(l-a)` for some
`0<=a<=l`, and `x` is inserted at that cut of the chain.

At an endpoint it has one homogeneous neighbor. At an interior cut it
replaces the deleted homogeneous edge by two homogeneous edges through
`x`. No other old chain changes. This proves (2).

Equivalently, if `u->x->v` are two new homogeneous edges, then every
other old vertex sees `u` and `v` identically, and homogeneity forces
`u->v`. Hence `{u,v}` was an old homogeneous pair: the two pieces joined
through `x` belong to the same old chain. Joining distinct old chains
cannot occur. QED.

This proves the proposed split-and-recompress algorithm, with a sharper
recompression rule. Either all changes amount to lengthening one old
chain, or there is no recompression at all after splitting the old chains
at their binary switches and adjoining isolated `x`.

The weighted quotient and all within-chain binary words determine the
extended tournament exactly. Keeping only the number of switches does
not determine that tournament; its sufficiency below concerns only the
number of child chains after summing all attachments.

## 3. Exact labelled transition kernel

**Proposition 3.** A fixed labelled parent of order `n` with `b` canonical
chains has exactly `n+b` clone attachments. Among them `2b` are endpoint
insertions and `n-b` are interior insertions.

Each chain of length `l_i` has `l_i+1` insertion cuts. Its old transitive
order fixes the attachment at every other vertex, so each cut yields one
word. Distinct cuts give distinct words. Cuts belonging to different
chains cannot give the same word: that would make `x` homogeneous with
vertices of distinct old chains, contradicting Proposition 2. Summing
`l_i+1` proves the count.

For arbitrary attachments, each chain permits a free initial bit and a
free decision to switch at each of its `l_i-1` edges. Hence the number
with exactly `s` switches is

    2^b * binom(n-b,s).                                     (3)

Combining (2) and (3) gives an exact polynomial in the child core size:

    K_(n,b)(y)
      = (n+b)y^b
        + y^(b+1) [2^b(1+y)^(n-b) - 2b - (n-b)y].           (4)

Explicitly, the counts of attachments are:

| Child core order | Number of attachment words |
|---|---:|
| `b` | `n+b` |
| `b+1` | `2^b-2b` |
| `b+2` | `(2^b-1)(n-b)` |
| `b+1+s`, `2<=s<=n-b` | `2^b binom(n-b,s)` |

Terms outside the indicated range are zero. The row sum is `2^n`.
All these coefficients are nonnegative on actual parent states. There
is no tournament whose canonical pair-free quotient has order two.

The distribution depends only on `(n,b)`, not on the quotient tournament
or the individual chain lengths. This is a genuine exact aggregation:
under uniformly random independent attachment bits, core order has the
transition probabilities obtained by dividing (4) by `2^n`.

Let `L_(n,b)` count tournaments on a fixed labelled vertex set of size `n`
whose core order is `b`. Deleting the designated last label is a unique
parent operation, so

    L_(1,1)=1,
    L_(n+1,c)=sum_b L_(n,b) [y^c] K_(n,b)(y).                (5)

There is no isomorphism quotient or parent multiplicity in (5). It is
an exact recurrence for labelled strata, and its row sums recover
`2^binom(n,2)`.

As an independent first-moment check, a given pair is homogeneous with
probability `2^(-(n-2))`: each of the other `n-2` vertices must see both
endpoints alike. Since `b=n-|E(J)|`, for a uniform labelled tournament

    E[b] = n - binom(n,2) / 2^(n-2).                        (6)

The mean of (4) gives the same identity recursively. This checks a
consequence using individual pair events instead of chain-run counting.

## 4. Alternating generation of pair-free cores

**Corollary.** The child `T+x` is pair-free exactly when every old
homogeneous edge is broken and `d=0`. Equivalently:

1. the attachment word strictly alternates along each old chain;
2. `x` has no homogeneous neighbor.

There are only `2^b` candidate alternating words, one starting bit per
old chain, instead of `2^n` arbitrary attachments.

This gives a complete recursive construction of pair-free tournaments:
take every parent of order `n`, express it by its smaller pair-free
weighted core, try its alternating words, reject clone cases, and
canonicalize the resulting pair-free child. Completeness follows by
deleting any vertex of any pair-free child and applying the criterion.
Duplicate removal remains necessary.

The number of successful candidates can be read exactly:

    n=b:     2^b-2b;
    n=b+1:   2^b-1;
    n>=b+2:  2^b.                                          (7)

For `n=b+1`, the single old homogeneous edge can be restored only by its
unique interior clone attachment. For `n>=b+2`, all-edge breaking has
at least two switches, so Proposition 2 excludes every clone attachment.

The hostile is also constructive. For a transitive parent of order
`n>=3`, both alternating attachment words produce pair-free children
of order `n+1`, despite the parent having a one-vertex core. Thus the
grammar has a concise production rule but no bound on its child core
order in terms of parent core order alone.

## 5. Labelled core EGF and the unlabelled boundary

Let `t_n=2^binom(n,2)` count labelled tournaments, with `t_0=1`, and let
`c_n` count labelled pair-free tournaments, with `c_0=1`. Their exponential
generating functions are formal power series

    T(z)=sum_(n>=0) t_n z^n/n!,
    C(z)=sum_(n>=0) c_n z^n/n!.

No analytic convergence is being asserted. The unique canonical-chain
decomposition gives

    T(z)=C(z/(1-z)),
    C(z)=T(z/(1+z)).                                       (8)

For a direct count, fix core order `b`. Partition the `n` labels into
`b` nonempty ordered chains and put a labelled pair-free tournament on
their block set. There are

    L_(n,b) = c_b * n!/b! * binom(n-1,b-1),                 (9)

because the `binom(n-1,b-1)` positive compositions give lengths, a
permutation of all `n` labels orders and fills the chains, and the
temporary numbering of the `b` blocks is counted `b!` times. The quotient
is pair-free, so every such expansion has exactly those canonical chains.
This proves (9) without a generating-function assumption, and summing
over `b` proves (8).

For every `n>=1`, coefficient extraction gives

    c_n = n! sum_(k=1)^n
          (-1)^(n-k) binom(n-1,k-1) t_k/k!.                 (10)

Equations (4), (9), and (10) are independent descriptions of the same
labelled strata. The root/core quotient, chain order, and labels are
load-bearing. Removing the labels does not turn (8) into an identity
for the ordinary generating function of A000568.

For example, consider the pair-free directed triangle as quotient.
At total order six the ten ordered positive weight triples summing to
six fall into four cyclic-rotation orbits:

    (4,1,1), (3,2,1), (3,1,2), (2,2,2).

The first three weighted tournaments have trivial automorphism group;
the last has automorphism group of order three. Their labelled masses
are `720,720,720,240`, summing to `2400`, exactly (9) with `c_3=2`.
Counting unlabelled classes gives four, and requires the cyclic action
on the weight positions. A weight multiset would also wrongly merge
the middle two classes. This is the minimal relevant distinction:
weighted quotient up to its actual automorphisms, not just a list of
chain lengths.

There are exact Burnside-free class recurrences, but they still need
isomorphism work. One can enumerate all attachment words to every parent
class, canonicalize children, and deduplicate, as the inherited code does.
Alternatively, choose an isomorphism-invariant canonical vertex orbit
of each child. Keep precisely root-preserving extension orbits whose
new vertex lies in that chosen orbit. Every unrooted child then appears
once. Parent attachment orbits and the child's canonical deletion orbit
must both be retained; a rooted augmentation count alone overcounts.

For a useful weighted identity, let `r(S)` be the number of vertex orbits
of `Aut(S)`. Summing `1/r(S)` over every rooted extension orbit of every
parent class counts each unrooted child exactly once. The equivalent
all-attachment identity is

    a_(n+1) = 1/(n+1) *
              sum_[T], |T|=n  1/|Aut(T)| *
                 sum_(z in {0,1}^n) |Aut(T+x_z)|.            (11)

Indeed, a rooted child `(S,x)` occurs in
`|Aut(T)|/|Aut(S)_x|` attachment words over its parent class.
Its contribution to (11) is
`|Aut(S)|/((n+1)|Aut(S)_x|)`, the size of that vertex orbit divided
by `n+1`. Summing over the child vertex orbits gives one.

Equation (11) is exact bookkeeping, not a compressed algorithm:
the child automorphism factor is precisely information absent from
the scalar labelled kernel. No claim is made that it can be inferred
from `(n,b)`.

## 6. Controls, scope, and next audit

A separate read-only Python probe in this inheritance pass generated
every labelled parent through order five and every new-vertex attachment:
`33,866` extensions in total. It computed homogeneous pairs directly
from adjacency rows, rather than calling the canonical-chain
implementation. It verified:

- retention of old pairs exactly at equal attachment bits;
- absence of newly created old-old pairs;
- equation (1), `b'>=b`, and both cases of (2);
- that two new homogeneous neighbors were an old homogeneous pair;
- the alternating criterion for pair-free children.

The `d>0` cases split as `8,030` endpoint insertions and `1,390`
interior insertions. Transitive parents of orders three through seven
with both alternating words gave pair-free children; at order two one
alternating word is the interior clone and the other gives a triangle.
These finite checks support but do not replace the proofs.

A separate exact-rational coefficient check agreed between (5), (9) and
(10) through order eight. The labelled pair-free counts `c_1,...,c_8`
were `1,0,2,16,344,14528,1172816,180364032`. This second check compares
algebraic formulas and does not constitute another exhaustive class census.

The parent session's independent reproducible checker supplies the saved
kernel/EGF/core-bank controls. This note does not promote a finite
absence of collisions to an all-size quotient theorem. The strongest
remaining target is a smaller operation-stable carrier for weighted
pair-free cores that preserves child class, clone allocations, and
ordinary line transport. The exact scalar labelled kernel is a positive
example of safe compression; equations (10) and (11) pinpoint why its
isomorphism-class analogue needs more information.

## 7. Full cycle-index transport and the undetermined boundary

The scalar labelled series in (8) is the identity-permutation part of a
stronger exact identity. Define the full cycle indices

    Z_T(p)=sum_(n>=0) (1/n!) sum_(g in S_n)
             Fix_T(g) product_i p_i^(m_i(g)),
    Z_C(p)=sum_(n>=0) (1/n!) sum_(g in S_n)
             Fix_C(g) product_i p_i^(m_i(g)),               (12)

where `Fix_T` counts labelled tournaments fixed by the permutation and
`Fix_C` counts pair-free ones. Include the empty object in both series.
Work over rational formal series with weighted degree `deg(p_i)=i`.
Every stated coefficient has a finite calculation; analytic convergence
is neither needed nor asserted.

**Proposition 4, chain substitution.**

    Z_T(p)=Z_C(p_i/(1-p_i)),
    Z_C(p)=Z_T(p_i/(1+p_i)).                              (13)

Here the notation means simultaneous substitution for every variable.
The ordered nonempty-chain species has cycle index `p_1/(1-p_1)`:
only the identity permutation fixes a linear order. Substituting these
chains at the vertices of the canonical pair-free quotient gives the
first identity. Explicitly, a quotient-automorphism cycle of length `i`
requires equal chain lengths around that cycle. Length `k>=1` creates
`k` vertex cycles of length `i`, contributing `p_i^k`. Summing over `k`
gives `p_i/(1-p_i)`. This is precisely the automorphism action induced
by the unique canonical decomposition, not an arbitrary representation.
The substitutions are formal inverses, proving the second identity.

**Proposition 5, fixed-root extension.**

    partial_(p_1) Z_T(p) = Z_T(2p_1,2p_2,...).             (14)

To prove this, fix an old permutation `g` and add one distinguished
fixed vertex. An invariant tournament on the old vertices has exactly
two choices of its relation to the new vertex on each vertex cycle of
`g`. Hence the number of invariant extensions is
`2^(sum_i m_i(g)) Fix_T(g)`. Differentiating the cycle index at `p_1`
marks and removes one fixed point, including the correct factorial
normalization. This proves (14).

Combining (13) and (14) by the formal chain rule gives the requested
core identity:

    (1+p_1)^2 partial_(p_1) Z_C(p)
       = Z_C(2p_i/(1-p_i)).                              (15)

For clarity, first write `q_i=p_i/(1-p_i)`. The left derivative in
(14) is `(1-p_1)^(-2) partial_(q_1) Z_C(q)`. On the right,
`Z_T(2p)=Z_C(2p_i/(1-2p_i))`. Now put `p_i=q_i/(1+q_i)` and
rename `q` as `p`, obtaining (15).

This identity has an exact coefficient recurrence. Write

    Z_C(p)=sum_(r>=0) c_r(bar p) p_1^r,
    bar p=(p_2,p_3,...),
    Phi(bar p)_i=2p_i/(1-p_i).

Then

    c_1 = c_0(Phi(bar p)),                               (16)

and, for `r>=1`,

    (r+1)c_(r+1)
      = sum_(j=1)^r 2^j binom(r-1,j-1) c_j(Phi(bar p))
        - 2r c_r - (r-1)c_(r-1).                        (17)

Indeed, the coefficient of `p_1^r` on the left of (15) is
`(r+1)c_(r+1)+2r c_r+(r-1)c_(r-1)`. On the right, for
`r>=1`, expand
`(2p_1/(1-p_1))^j`; its coefficient is
`2^j binom(r-1,j-1)`. The constant coefficient gives (16).

**Boundary.** Equations (16)--(17) uniquely reconstruct every coefficient
containing `p_1` from the formal boundary

    c_0(bar p)=Z_C(0,p_2,p_3,...).                        (18)

They do not determine that boundary. Any prescribed rational formal
boundary yields a unique formal solution of (15) by the same recurrence;
being a genuine cycle index of a combinatorial class is an additional
constraint. In particular the labelled specialization, even together
with the prohibition of even permutation cycles, does not determine the
coefficient of `p_3`.

Precisely, (18) records **pair-free cores with fixed-point-free
automorphisms**, including the empty term. It does not mean a new
subclass of tournaments defined by having no fixed vertices. Ordinary
and pair-free tournament automorphism groups have odd order, so terms
containing any even-length permutation cycle vanish. The surviving
nonconstant boundary begins

    1 + (2/3)p_3 + (4/5)p_5 + terms of higher weighted degree. (19)

The first independent coefficient audit enumerated every labelled
tournament through order five and every permutation of its vertices,
testing invariance and pair-freeness directly from adjacency. It verified
both substitutions in (13) through weighted degree five and both
derivative identities (14)--(15) through degree four. The latter cutoff
is necessary: differentiating a degree-five input only verifies degree
four output. This test imported no existing enumeration or cycle-index
implementation.

## 8. What the functional identity adds beyond the unit-tail kernel

[THM-514, the unit-tail core kernel](../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md)
already peels the fixed-point part of a permutation cycle type. It uses
the Davis/Burnside fixed-point formula, with the remaining permutation
cycle lengths at least three, and proves the universal tail law
`F'(x)=F(2x)`. Its "core" is the residual **permutation cycle type**.
The core in this note is the canonical **pair-free tournament quotient**.
These are different objects.

Equation (15) transports the unit-tail mechanism through an explicit
canonical-chain substitution and gives a recursive law for the actual
pair-free tournament cycle index. Its proof uses invariant attachment
choices rather than importing the whole Burnside formula. Nevertheless
it does not remove the fixed-point-free symmetry boundary: the
unresolved data of (18) are the counterpart of THM-514's hard kernel.
The constructive advance is the precise change of coordinates and
operation law, not an autonomous computation of A000568.

There is a complete hierarchy that exposes the remaining cost. For any
odd positive integer `ell`, adjoining one permutation cycle of length
`ell` gives

    partial_(p_ell) Z_T(p)
      = kappa_ell Z_T(alpha_j p_j),
    kappa_ell=2^((ell-1)/2)/ell,
    alpha_j=2^gcd(ell,j).                               (20)

To verify (20), within the new odd cycle the unordered pairs have
`(ell-1)/2` orientation orbits. Between it and an old cycle of length
`j`, there are `gcd(ell,j)` pair orbits, each with two orientation
choices. The cycle-index differentiation contributes `1/ell`.
Old types with an even cycle have no invariant tournament and contribute
zero on both sides. Thus the formula also holds with all variables
present. Applying the same chain substitution gives

    (1+p_ell)^2 partial_(p_ell) Z_C(p)
      = kappa_ell Z_C(
          alpha_j p_j / (1-(alpha_j-1)p_j)).              (21)

For `ell=1` this is (15). The hierarchy for all odd `ell`, together
with the vanishing of even-cycle terms and constant term one, does
determine the full cycle index recursively. It does so by recovering
the same unbounded gcd interactions that appear in the Davis formula
and THM-514's prime-tail ladder. It therefore supplies no demonstrated
complexity improvement over that formula. Using only finitely many
equations leaves a corresponding new-cycle boundary.

The nearby
[sign-aware augmentation compiler](../sign_aware_augmentation_20261020/signed_root_augmentation.md)
retains an actual parent stabilizer and its character while adding
distinguished roots; it explicitly audits the defect when those roots
are forgotten. This agrees with the distinction here. Equations
(14)--(21) are fixed-permutation inventories, hence suitable for
counting. They do not give a transporter or recover which unrooted child
receives an attachment. A full cycle index can determine many orbit
counts while still failing to determine an operation on individual
objects. The canonical root-orbit rule or the explicit weights in
(11) remain necessary for unlabelled class generation.

## 9. Independent source audit of the growth checker

The root session's `growth_check.py` was read independently after this
note's proofs and read-only probes. Its local universe consists of all
labelled parents through order five and all attachments; its mask
convention is the complement of `z_v` used here, which preserves
switches, clone cases, and the alternating gate.

The alternating bank correctly starts with the singleton core, obtains
new pair-free cores from every prior weighted class, and reconstructs
all classes from positive weights on the accumulated core bank. The
per-child mass `Aut(child)/(n*Aut(parent))` is equation (11) restricted
to pair-free children, with `n` denoting child order in that script.
Its explicit check that each child receives mass one tests the
root-forgetting correction, not just a total count.

The canonicalizer, chain routine, and fixed-path comparison are imported
from the frozen preceding packet, so those particular controls share
those implementations. The inheritance pass's direct homogeneous-pair
probe and labelled/permutation cycle-index audit do not. The source
correctly uses explicit exceptions rather than optimization-disabled
assertions. No mathematical correctness issue was found in the advertised
standalone execution path. This is a source audit; it does not claim an
additional complete order-eight execution by this reviewer.

## 10. Incoming pure-blue alphabet claim: a proof-test gap, not a refutation

During this audit, incoming commit `70ca1d8210` promoted the previously
reserved namespace to
`01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md`,
explicitly marked pre-audit. Its odd-order argument says that a Hamiltonian
path is non-palindromic when the sums of its reversal partners are not
constant in cyclic coordinates. That tests whether positional reversal
is an affine reflection `v -> c-v`; it does not test all
anti-automorphisms of a composite-order circulant.

There is an explicit witness to the failed implication. In the circulant
on `Z_9` with connection set `{1,3,4,7}`, the Hamiltonian path

    (0,1,2,3,7,5,6,4,8)

has positional reversal, written as images of vertices `0,...,8`,

    (8,4,6,5,1,3,2,7,0).

Direct inspection of every unordered pair verifies that this reversal is
an anti-automorphism. But the partner sums modulo nine are both five and
eight. Thus this is a palindromic path with a nonaffine reversal. The
incoming companion `atom_alphabet_proofcheck.py`, function `palindromic`,
uses only the constant-sum test in its odd-circulant section and would
misclassify this path.

This does **not** refute that theorem's special three-interval path
claim, or its proposed alphabet classification. An independent direct
anti-automorphism test checked all `103,423` eligible three-interval
paths across every connection set containing one at odd orders
`7,9,...,25`, and found no counterexample to that specific path claim.
The distinction is essential: the displayed general inference is false,
while the specialized claim may admit another proof. An actual arc
witness or an additional rigidity argument is needed.

The new weighted-core constructions in this packet do not depend on
THM-4730, so this proof gap does not affect them. The safe inheritance is
the previously proved stack center law, together with explicitly bounded
pure-blue census results, until the new odd-order argument is repaired
and independently audited.

**Reproducibility and repair update.** The witness and full finite probe are
now retained in [incoming_audit.py](incoming_audit.py) and its matching
[output](incoming_audit.out). The incoming checker was repaired to test
actual arc reversal; the all-order theorem and dependent routers carry an
explicit proof-gap override. Its selected-path normal run still passes.
The corrected universe totals are 2,044 odd-circulant connection sets
(2,035 three-interval paths and nine descending paths) and 1,198,372
bicirculant mixing cycles. Neither finite success repairs the general
argument automatically.
