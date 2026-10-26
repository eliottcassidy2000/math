# Rooted incidence, exact augmentation, and an infinite family-lumping obstruction

**2026-10-26. Status:** PROVED elementary statements below; FINITE-EXACT
verification explicitly bounded. This packet independently verifies existing
all-order enumeration mechanisms and sharpens their family-state obstruction.
It does not claim a new formula for A000568, a polynomial algorithm, or a
general impossibility theorem for scalar recurrences.

## Results first

1. There is an infinite, explicit failure of one-step **unlabelled** family
   lumping even when the actual pair-free core, total order, SCC-size sequence,
   and full automorphism cycle index agree. For every `k>=2`, the two rigid
   strong tournaments

       A_k = C3[TT_(3k-3), TT_1, TT_1],
       B_k = C3[TT_(k-1), TT_k, TT_k]

   have order `3k-1` and core `C3`, but their symmetry-corrected transition
   mass to core order three differs by exactly `2/3`. The lost observable is
   alignment of weights with **future** symmetry, not current symmetry.
2. Canonical selection must choose a **root orbit**. Selecting a minimum
   unrooted deletion card fails already at order two. Dividing the number of
   attachment orbits by the child order first fails at order three: `4/3 != 2`.
3. The fresh standalone generator accepts each unrooted child exactly once,
   without relying on global duplicate removal. It independently recovers
   all `6880` classes through order eight. Every child separately has mass
   one in three corrected rooted-incidence formulas.
4. Insertion merges exactly one interval of the old ordered SCCs, or inserts
   a singleton at one cut. This gives a product formula for strong **rooted**
   extension orbits, not strong unrooted classes.

## Inheritance, scope and priority

Read against these existing results, rather than treating this as a fresh
solution of the enumeration problem:

- `01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md`:
  unique maximal transitive chains and their pair-free quotient.
- `05-knowledge/results/weighted_core_transport_20261020/README.md`,
  `inheritance.md`, `core_count.md`: all-order structural recurrence,
  weighted-core orbit counts, exact attachment symmetry weights, labelled
  `(n,b)` transition kernel, and equal-cycle-index core-growth hostiles.
- `05-knowledge/results/tournament_future_boundary_20261020/README.md`:
  closed finite-horizon count recurrence and its wordwise-minimal future
  response vector. It already computes ordinary counts through order 200.
- `05-knowledge/results/defect_boundary_20261020/boundary.md`:
  live divisor elimination and the full fixed-point-free operator boundary.
- `01-canon/theorems/THM-4169-prime-parent-one-vertex-augmentation-and-quartic-johnson-transfer.md`,
  section 6: exact rooted-card incidence in the restricted prime bank.
- `01-canon/theorems/THM-260-rooted-tournament-layer-decomposition.md`:
  rooted counts and their automorphism deficit; retain its denominator repair.

The historical `04-computation/a000568_fiber_recursion.py` begins with the
invalid divide-by-order formula and immediately acknowledges a counting
problem. Its SHA-256 at inherited commit `b4234f728c05` is
`77e2554bb2857246a7a92287a99a0d984096d1637aea009ef82582057ad9ef94`.
The correct raw-attachment formula was already present in the October 20
packet. This note clarifies the old failure; it does not pretend the repair
was missing from current canon. The historical wording "overcounts" should
not be attached to the divided expression: that expression **undercounts**.

Canonical augmentation is established prior art: Brendan McKay,
[Isomorph-free exhaustive generation](https://users.cecs.anu.edu.au/~bdm/papers/orderly.pdf),
J. Algorithms 26 (1998), 306–324, author version with errata (checked
2026-10-26). The construction below is a self-contained specialization, not
a priority claim. The sequence comparison is [A000568](https://oeis.org/A000568).

Concept board: weighted cores; attachment bits; root orbits; prospective
stabilizers; SCC intervals; exact family lumpability. Anchor: complete
enumeration. Niche: a smaller strong-rooted counting observable. Wildcard:
symmetry appearing after one clone insertion.

## 1. The exact rooted-incidence correspondence

Let `U_n` be the set of unlabelled tournament classes of order `n`, including
one empty class. For `T in U_n`, let `X(T)={0,1}^{V(T)}` be attachment masks,
where `z_v=1` means the new vertex `x` beats `v`. Let `G_T=Aut(T)`.

There is a bijection

    disjoint_union_[T in U_n] X(T)/G_T
       <--> {rooted tournament classes [(U,x)] : |U|=n+1}.       (1)

Deleting the root constructs the inverse: an identification `U-x -> T`
determines a mask, and changing the identification is exactly the `G_T`
action. Root-preserving isomorphisms are exactly the same equivalence.

Write `r(U)=|V(U)/Aut(U)|`. Forgetting the root has fibre size `r(U)`, not
`n+1`. Consequently, with one representative of each attachment orbit,

    a_(n+1) = sum_[T] sum_[z] 1/r(U_z)
             = 1/(n+1) sum_[T] sum_[z] |Orb_Aut(U_z)(x)|.       (2)

The second identity uses the fact that the sizes of all vertex orbits sum
to `n+1`. Since `Aut(U_z,x)` is the attachment stabilizer in `G_T`, lifting
the orbit sum to all masks gives the inherited formula

    a_(n+1) = 1/(n+1) sum_[T] 1/|G_T|
                          sum_(z in X(T)) |Aut(U_z)|.           (3)

These identities hold separately for every isomorphism-invariant target
family, by restricting the children. They remain sums over actual objects;
they do not magically evaluate the child automorphisms.

At order three the transitive tournament has three vertex orbits and `C3`
has one. Thus there are four rooted classes but two unrooted classes.
The bad formula gives `4/3`, whereas (2) gives `3*(1/3)+1=2`.
Orders one and two have only rigid tournaments, so three is the first failure.
Indeed the divided expression strictly undercounts at every child order
`N>=3`: each `r(U)<=N`, and `C3 -> TT_(N-3)` has a nontrivial automorphism
(use `C3` itself at `N=3`), hence strictly fewer than `N` vertex orbits.

## 2. Exactly-once canonical augmentation

Choose any complete rooted canonical code `c(U,v)`. Its minimum over roots
is attained on exactly one automorphism orbit, because equality of rooted
codes is equivalent to rooted isomorphism. Call this orbit `D(U)`.

Generate one attachment representative per `G_T`-orbit for every parent
class `T`. Accept `(U,x)` exactly when `x in D(U)`.

**Proof of completeness and uniqueness.** Every child has a nonempty selected
root orbit. Delete a selected root; its parent is in the complete preceding
level, and (1) generates its rooted class. Conversely, all accepted roots in
one child belong to a single orbit, so (1) supplies exactly one accepted
presentation. Induction from the empty tournament proves an all-order
enumerator. No global child-set deduplication is needed for correctness.

Using only `min_v canon(U-v)` is not enough: in the unique two-vertex
tournament both deletion cards are singletons, but the source and sink are
different root orbits. The minimum-card test accepts that child twice.

The verifier does store the resulting class set for auditing, but checks
that every accepted code occurs **exactly once before** making that set.
This distinguishes exact-once generation from a generator repaired by
post-hoc deduplication. It makes no speed advantage claim.

## 3. What a scalar family recurrence actually requires

For a child family `G`, define the exact rational transition

    K_n(T,G) = 1/((n+1)|Aut(T)|)
                 sum_(z : T+x_z in G) |Aut(T+x_z)|.            (4)

Then `a_(n+1)(G)=sum_[T] K_n(T,G)`. If a proposed parent partition `F`
has a common value of (4) for every member of each block and every target
block, its family counts satisfy the corresponding matrix recurrence.
Conversely, this constancy is necessary for a transfer that works for
**arbitrary weights assigned to parent classes**: compare unit mass at two
members of the same block. This is the precise strong-lumpability criterion.

Nonconstancy does not forbid an accidental identity for the one actual
global sequence, analytic averaging over a whole family, or a different
symbolic count method. No such general no-go is asserted.

### Same-core transfer and prospective stabilizers

Write a tournament in its unique canonical form `T=Q[TT_(l_i)]`, with
`G=Aut(Q)`, `H=Stab_G(l)`, and total size `n`. There are exactly `l_i+1`
clone insertion masks at chain `i`. Each gives the weight vector `l+e_i`
and has automorphism group `H_i=Stab_G(l+e_i)`. All other masks strictly
increase core order, by the inherited clone-growth dichotomy. Therefore

    K_n(T, core order |Q|)
       = sum_i (l_i+1)|H_i| / ((n+1)|H|).                     (5)

To refine to a particular weight orbit `[m]`, restrict the sum to indices
with `[l+e_i]=[m]`. This is a partial flow: other parent cores may also
produce a given child. The present stabilizer `H` is not enough; (5) sees
the stabilizers after incrementing the weights.

### Infinite hostile family

Let `k>=2`, `n=3k-1`, and use `A_k,B_k` from the opening. Their core is the
oriented triangle, whose automorphisms are its three rotations. A positive
weight vector has nontrivial stabilizer iff all three weights are equal.
Thus both parents are rigid, strong, have the same actual core `C3`, and
have cycle index `p_1^n`.

For `A_k`, no one-coordinate increment makes the three weights equal.
Every clone child is rigid. The sum of insertion-cut counts is `n+3=3k+2`.
For `B_k`, the unique deficient chain has length `k-1`, hence has exactly
`k` insertion cuts. Those cuts produce the uniform child `(k,k,k)`, with
automorphism order three; every other clone child is rigid. Equation (5)
gives

    K_n(A_k, core order 3) = (3k+2)/(3k),
    K_n(B_k, core order 3) = (3k+2+2k)/(3k),
    difference = 2/3.                                      (6)

This is an all-order obstruction to the specified universal family
transfer. It does not contradict the exact labelled `(n,b)` kernel:
the two parents have identical raw mask counts but different symmetry mass.

At `k=2`, independent direct permutation enumeration gives:

| Parent | Raw counts at child core orders 3,4,5,6 | Corrected masses at 3,4,5,6 | Distinct children |
|---|---|---|---:|
| `C3[3,1,1]` | `8,2,14,8` | `4/3,1/3,7/3,2` | 24 |
| `C3[1,2,2]` | `8,2,14,8` | `2,1/3,3,4/3` | 23 |

Both have 32 rooted attachment classes. Their child-automorphism histograms
are respectively `30*1 + 2*3` and `28*1 + 4*3`. The explicit row masks are
`(16,1,3,7,14)` and `(2,24,3,5,13)`; bit `j` in row `i` means `i -> j`.
The second is isomorphic to the stated weighted triangle; clone-size
multisets are `(1,1,3)` and `(1,2,2)`.

Order five is the first hostile in the exact tested family key
`(order, core order, ordered SCC sizes, full Aut(T) cycle index)`.
All classes through order four have distinct such keys. It is also the
first possible split between positive weight orbits on a fixed `C3` core:
there is only one such orbit at each total size three and four.

## 4. The SCC interval law and a genuinely smaller rooted observer

Let the old SCCs be `C_1 -> ... -> C_s`. Set

    L = min{i : x beats some vertex in C_i}, default s+1;
    R = max{i : some vertex in C_i beats x}, default 0.

If `L>R`, then necessarily `L=R+1` and `x` is a singleton SCC inserted at
that cut. If `L<=R`, exactly `C_L,...,C_R` merge with `x`; all other SCCs
remain unchanged and in the same order.

**Proof.** If `L<R`, the directed path from `C_L` through successive
components to `C_R`, together with `x -> C_L` and `C_R -> x`, makes every
vertex in that interval mutually reachable with `x`. If `L=R`, the two
opposite attachment directions occur within one old strong component;
its internal paths give the same conclusion. Earlier components all beat
`x` and later components all lose to `x`, so none can join. If `L>R`, no
component can lie strictly between them: any vertex there must either
beat or lose to `x`, contradicting one extremal definition. QED.

Let `E(C)` be the number of subset orbits of `Aut(C)`, computable by
`|Aut(C)|^-1 sum_g 2^(number of vertex cycles of g)`. Ordered SCCs cannot
be permuted by an automorphism, so `Aut(T)=product_i Aut(C_i)`.
The number of attachment orbits producing a strong child is therefore

    1                                       if T is empty;
    E(C_1)-2                                if s=1;
    (E(C_1)-1)(E(C_s)-1) product_(1<i<s) E(C_i)  if s>=2.     (7)

For `s>=2` the new vertex must beat at least one vertex of the first SCC
and lose to at least one in the last; each forbidden endpoint condition
removes one subset orbit. For one SCC, remove the empty and full masks.
This is a valid smaller statistic for this single rooted counting query.
It cannot replace the child-orbit correction in (2)–(4).

Writing `A(x)` for all unrooted classes, `S(x)` for nonempty strong
unrooted classes, `B(x)` for rooted classes, and `Q(x)` for strong rooted
classes, the unique ordered SCC decomposition gives

    A = 1/(1-S),              B = A^2 Q.                      (8)

These are formal series identities. They do not evaluate the unknown
strong counts on their own and are not claimed as new prior art.

## 5. Verification, loss ledger and remaining boundary

`verify_recursion.py` is standalone standard-library source. It uses exact
integers/rationals and explicit exceptions, never disabled `assert`s.

| n | Unrooted | Rooted | Strong unrooted | Strong rooted |
|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 1 |
| 2 | 1 | 2 | 0 | 0 |
| 3 | 2 | 4 | 1 | 1 |
| 4 | 4 | 12 | 1 | 4 |
| 5 | 12 | 48 | 6 | 24 |
| 6 | 56 | 296 | 35 | 194 |
| 7 | 456 | 3040 | 353 | 2395 |
| 8 | 6880 | 54256 | 6008 | 47588 |

It checks: exact-once canonical augmentation through eight; each child
separately in (1)–(3); independent odd-partition Burnside counts; total
labelled orbit mass; full labelled class universes through five; every
labelled parent through five and every attachment (`33,867` trials) in the
SCC interval law; (7) for every generated parent through seven; and (8).

`audit_hostile.py` independently uses matrix entries and all `6!`
permutations per child, without degree pruning or imports from the other
verifier. It checks all 64 attachments, `46,080` child permutations,
the full raw and weighted histograms, and distinct-child counts. It emits
every `[attachment mask, child core order, child automorphism order]`.

Reproduce:

```sh
python3 verify_recursion.py --max-order 8 > receipt-order8.json
python3 -O verify_recursion.py --max-order 8 > replay-order8.json
cmp receipt-order8.json replay-order8.json
python3 audit_hostile.py > hostile-receipt.json
python3 -O audit_hostile.py > hostile-replay.json
cmp hostile-receipt.json hostile-replay.json
```

Both normal/optimized comparisons passed. No inherited executable or
repository script was executed as a mathematical verifier.

**Connection contract.** Full tournament -> weighted core preserves the
unmarked tournament class; forgetting the weights or their prospective
stabilizers loses (5), as (6) shows. Pointed extension -> unrooted child
forgets the root orbit; (2) supplies the correction. SCC tuple -> endpoint
subset-orbit counts preserves only the query (7), not child identities.
None of these claims turns a current count into a future-complete state.

The complete weighted core remains a sufficient structural state. The
previous future-response vector remains sufficient, and wordwise-minimal
for its specified finite-horizon Burnside continuation task. A smaller
all-order structural family state is still OPEN. The exact obstruction
proved here is narrower and useful: present core, size, SCC sizes and
present automorphism cycle index cannot by themselves supply the universal
unlabelled transfer. Prospective symmetry is already visible at order five.
