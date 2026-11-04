# Integral forest-polynomial transport for a growing port grammar

2026-11-03. **PROVED** typed birth, flip, converse and boundary-forget rules;
**FINITE-EXACT** full-profile controls and the literal inner-defect transport.
Two-branch return closure and global Hamiltonian-count completeness remain
**OPEN**. No global theorem identifier is reserved here.

The new carrier has finitely many polynomial coordinates at each fixed active
boundary size. It retains unbounded component depth. It is substantially
smaller than a full named-vertex squarefree matrix when the active boundary
is small, but its compression has an explicit operational restriction: a new
vertex must see every forgotten vertex in the same direction.

## 1. Inheritance and the actual new extension

The closest proved mechanism is
[`ACTIVE_PORTS.md`](../../adaptive_metagraph_reformulation_20261103/ACTIVE_PORTS.md),
the endpoint-pairing carrier with uniform future incidence to forgotten
vertices. That carrier prunes states with two permanently forgotten ends
under its one-direction birth language. Here future uniform directions may
change, so those states are retained as polynomial degree. The new extension
adds exact **edge reversals** and their mixed terms, arbitrary boundary
forgetting, and explicit transport of the rooted defect.

The August 2 mechanism behind component decomposition is **THM-3121**,
`01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`.
The need for an unbounded path-cover jet is inherited from **THM-3134**,
`01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md`.
The mixed-insertion warning is **THM-4099**,
`01-canon/theorems/THM-4099-squarefree-gap-transfer-and-mixed-insertion-boundary.md`.
No priority claim is made for their general path-cover principles.

The canonical hostile is the previous 6/8 port switch: the entire rooted
isomorphism type at31 agrees, but a newly distinguished physical port6 changes
the return. The corrected near miss is dropping component pairing, or
dropping higher component layers because the present observer only reads one
or two. The least-used sidecar here is the **number of fully anonymous path
components**, rather than another scalar clock. The live board consists of
endpoint pairings, a polynomial component counter, sparse flip supports,
positive root defects, and actual source words.

## 2. The state and its size

Let T be an actual tournament, with active physical ports P. A spanning
directed path cover is an **unordered** collection of nonempty disjoint
directed paths covering every vertex. Singleton paths are allowed.

For each component record its paired `(start,end)`. Retain an endpoint's
physical name if it lies in P, and otherwise replace it by `*`. Write

* `k` for the number of components of type `(*,*)`;
* `kappa` for the unordered multiset of the remaining endpoint pairs.

Let `c_T(k,kappa)` be the number of actual covers with this signature and set

```
C_T[kappa](z) = sum_k c_T(k,kappa) z^k.                 (T1)
```

An active vertex occurs at most once as an endpoint, except for `(p,p)`,
which represents its singleton component. It cannot be the endpoints of
two different components. Ports not appearing in kappa are internal to
their components. Their component positions need not be stored for the
operations below: the cutting identities supply them when an edge changes.

If `g_p` counts formal anchored signatures on p available names, then

```
g_0=1, g_1=4,
g_p=4g_(p-1)+2(p-1)g_(p-2).                           (T2)
```

Classify the least name as absent from endpoints, a singleton, `(*,p)`,
`(p,*)`, or paired with another name in one of two directions. Thus
`g_p=1,4,18,88,460,2544,14776,...`. These are upper bounds on indexed
coordinates, not ranks or claims of reachability. At order n the state has
at most `g_p(n+1)` integer coefficients. Initialization for an arbitrary
given graph can cost more; this is not a polynomial-time arbitrary-tournament
counting theorem. If every vertex must stay active, p itself grows with n.

The state is an observer, not a complete graph encoding. Alongside T1 retain
the active name list and its induced orientation as finite control. This
control supplies the actual old direction for a named reversal; birth,
forgetting and converse update it directly. Keep the actual construction
word and its cut instructions as the full source sidecar. Strongness and
native Hamilton-path guard certificates remain separate obligations.

## 3. Named joins, edge reversals and mixed terms

For distinct active ports a,b, define `J_ab` on a signature as follows. If
there is a component ending at a and a **different** component starting at
b, join them by a→b. Its new endpoints are the first component's start and
the second component's end. If both become anonymous, increment k. Otherwise
keep the new anchored pair. If either endpoint is unavailable or the two
components are the same, the result is zero.

No multiplicity or division is needed for this operation: each named endpoint
belongs to at most one component. In polynomial coordinates the matrix entries
of J are 0, 1 or z.

If the actual edge a→b is reversed, the exact update is

```
C_(T with b→a) = (I + J_ba - J_ab) C_T.                (T3)
```

Proof. Every old cover using a→b has a unique preimage obtained by cutting
that edge. It is subtracted by J_ab. Every newly allowed cover using b→a
has the unique cut preimage added by J_ba. That preimage is an old cover:
it cannot contain a→b, which would either connect the two components already
or produce a cycle after the new join. Covers using neither edge are unchanged.
The proof is at the level of actual covers and then pushes forward to T1.

All named joins commute. Adding a specified finite set of arcs to a forest
either respects all indegree/outdegree bounds and creates no directed cycle,
in which case every addition order gives the same final forest, or every
order is rejected at some step. In particular

```
J_ab^2=J_ab J_ba=J_ba J_ab=0.
```

For `R=J_ba-J_ab`, `R^2=0`; the reverse-oriented update is `I-R`, and
`(I-R)(I+R)=I`. The fixed operator `I+R` is **not** an involution over Z.

For a batch of distinct edge reversals on a fixed vertex set, use the product
of the T3 operators. They commute, but their mixed terms must be retained;
adding only first-order changes is incorrect. Births need not commute with
these operators.

Let lambda be the maximum number of edges in a linear forest contained in
the undirected flip-support graph. A nonzero product of selected joins must
itself form a linear forest, so the batch has join degree at most lambda.
Each join lowers the number of path components by one. Consequently the
post-batch one-component observer needs only initial layers through
`1+lambda`; its two-component observers need layers through `2+lambda`.
An arbitrarily large star has lambda≤2. The literal support `P3 disjoint K2`
has lambda=3. This is a local response bound, not a fixed truncation valid
under indefinitely many later operations.

## 4. Positive birth, exact forgetting and converse

For a new active vertex x, specify any incoming subset of active ports and
one bit epsilon giving the orientation from **every forgotten vertex** to x.
In each old cover choose at most one component entering x and at most one
component leaving x. Test its endpoint against the specified active cut or
the uniform anonymous bit. Choosing the same component twice is forbidden,
because that would create a cycle. Join the selected components through x,
or add x as a singleton if neither is chosen.

Choosing an anonymous component contributes its multiplicity k. In polynomial
coordinates this is differentiation: `d(z^k)/dz=k z^(k-1)`. A single uniform
birth cannot choose anonymous components in both roles, because their common
orientation would have to be both directions. Thus one birth uses at most
first derivatives. Compositions may use higher derivatives. The update is
a nonnegative integral matrix of polynomial differential operators.

Every new cover determines its unique old cover by deleting x and its incident
chosen arcs, keeping two old neighbors in separate components when x was
internal. This proves the birth rule without adding a nonexistent shortcut.

Forgetting a port replaces that name by `*` in every endpoint pair; a pair
becoming `(*,*)` multiplies the coordinate by z. This is a nonnegative exact
pushforward. Subsequent births must obey the uniform-anonymous rule, and
subsequent named flips and root selections must use active ports. A forgotten
name cannot later be distinguished without replaying or recomputing an
appropriate richer source interface.

That restriction is necessary even when the **whole** polynomial is retained.
The two labeled transitive pairs `0→1` and `1→0` have the same completely
forgotten profile `z+z^2`. Attempt to reactivate0 and add x with `0→x→1`:
the outputs have H=1 and H=3. The earlier profile cannot determine that cut.

Converse is the exact permutation `kappa:(a,b)↦(b,a)`, leaving k unchanged.
It exchanges start/end observers with no extra count coordinate. This is a
raw tournament operation; a legal native complement word is not inferred.

Component pairing is also essential. The two signatures
`{(a,b),(c,d)}` and `{(a,d),(c,b)}` have the same start/end marginals.
The b→c join is legal for the first and creates a cycle for the second.

## 5. Defect observers and the next layer

Write `d=k+|kappa|`. For a retained root r,

```
H = total weight with d=1,
E_r(2) = total weight with d=1 or2 and a component ending at r,
b_r = E_r(2)-H.                                      (T4)
```

For d=2 the r-ending component uniquely specifies the first ordered slot.
For d=1 the second slot is empty. A joint response `J_T(r,u)` is the total
weight with d=2, one component ending at r, and the other starting at u.
Splitting Hamilton paths immediately after r gives the positive defect
partition

```
b_r(T) = sum_(u→r) J_T(r,u).                          (T5)
```

T5 uses named incoming ports; if some are forgotten, T4 still evaluates the
defect exactly, but a positive name-by-name partition needs their retained
names or an appropriate endpoint color refinement.

An edge flip needs two-component joint data for H and three-component data
for E_r(2). Explicitly, let `e_a`, `s_b`, `G_br` and P denote respectively
the simple-path sums ending at a, starting at b, from b to r, and all nonempty
simple paths, with squarefree vertex weights. This notation is only an
independent readout evaluator; it is not the stored T1 state. Then

```
H(J_ab C_T)   = [V] e_a s_b,
E_r(2)(J_ab C_T) = [V] e_a G_br (1+P) + e_a s_b e_r.   (T6)
```

The first term has r at the end of the b-starting joined component, with an
optional unjoined component. The second has r at the end of the third
component. The squarefree product enforces disjointness. T6 also covers
coincident r=a or r=b: overlapping factors vanish and a root-to-itself path
is its singleton.

For a root-incident reversal `r→v` to `v→r`, an especially small identity is

```
E'_r(2)-E_r(2) = E_v(T-r,2),
H'-H = J_T(v,r)-J_T(r,v).                             (T7)
```

The removed r→v edge cannot occur in an r-ending first slot. A new cover
using v→r is exactly a v-ending first-slot cover of T-r, followed by r.
The inverse orientation negates the first increment, with the same deletion.

## 6. The exact 0→2960 collar transport

The carrier is the actual induced B14 from the saved
512425→1537279 return, not either full outer tournament. Preserve physical
names from `zero_ready_strong_switch_20261103/history_phase/collar.json`.
The old root32 has sole in-port33; the path of actual suffix flips is:

| Old arc reversed | H after | E32(2) after | b32 after |
|---|---:|---:|---:|
| initial | 46073 | 46073 | 0 |
| 4→33 | 40243 | 40243 | 0 |
| 32→22 | 73259 | 75961 | 2702 |
| 33→31 | 80319 | 83279 | 2960 |

The corresponding plus/minus join contributions are

| Reversal | H: added, removed | E32(2): added, removed |
|---|---|---|
| 4→33 | 1488,7318 | 1488,7318 |
| 32→22 | 34827,1811 | 35718,0 |
| 33→31 | 7932,872 | 8191,873 |

The middle 35718 is independently `E22(B14-32,2)` at that stage, checking T7.
At the end, T5 is exactly `2960=2070+890`, from incoming names22 and33.

The eight actual edge-flip corners provide a mixed-response control. In the
column order `(H,E32(2),end32,two-component covers ending32,b32)`, the third
mixed Boolean difference is

```
(396,396,0,396,0).
```

The nonzero two-component term confirms an actual contribution from the
initial **five-component** layer under the three-flip batch. Simply discarding
that layer loses 396. This does not prove that every conceivable alternative
coordinate compression has the same dimension. The defect's third mixed
term cancels; the next section explains why structurally.

## 7. A genuinely closed affine subgrammar for the inner defect

The inherited rooted cut identity, recorded in
[`ADAPTIVE.md` §3](../../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md),
gives, on D=T-r with root in-port set I,

```
b_r = sum_(i in I) start_i(D)
      + sum_({i,j} subset I) kappa_ij(D),
kappa_ij(D)=J_D(i,j)+J_D(j,i).                        (T8)
```

Now protect a vertex a dominating D-a, retain exactly two in-ports `{u,w}`
different from a, and vary only the w-star inside D while keeping a→w.
Both start terms in T8 are zero. Let Q=D-w, which is unchanged by every
allowed star reversal. Peeling w off its endpoint in each ordered cover gives

```
kappa_uw(D) = end_u(Q)+start_u(Q)
             + sum_(w→x) J_Q(u,x)
             + sum_(x→w) J_Q(x,u).                   (T9)
```

The singleton-w cases give the first two terms. All other w-starting or
w-ending paths have a unique adjacent x, giving the two sums. The x=u
terms are zero, so the orientation of u/w does not affect kappa.

Thus an arbitrary sequence of allowed star reversals transports the defect
with one fixed integer arm weight per name. Reversing w→x adds
`J_Q(x,u)-J_Q(u,x)`; reversing it back subtracts that number. Keep the actual
orientation bits. This is a closed affine defect grammar, with unbounded
numbers of possible arms as Q varies. It does not close the full H observer:
Hamilton paths can use two w-incident edges, whereas this endpoint-pair kernel
can use at most one. That proves the vanishing third defect term in §6.
In that cube the root-inport22 bit multiplies kappa; its coefficient is
affine in the other two33-star bits, so their joint second difference is zero.

For the actual collar, `(r,u,w,a)=(32,22,33,40)` and Q has12 vertices.
Its relevant exact values are

```
H(Q)=781, end22(Q)=3, start22(Q)=0,
(J_Q(22,4), J_Q(4,22))=(115,313),
(J_Q(22,31),J_Q(31,22))=(1,259).
```

Before either33-star edit the pair kernel is2900. The actual arm changes
give `2900-198+258=2960`, with2702 when the added in-port22 first appears.

T9 also gives a sharp all-pattern bound. With40→33 protected, independently
choose the smaller or larger of the two values at each remaining effective
arm. Across all `2^10=1024` effective orientation patterns,

```
614 <= b32 <= 7838.                                  (T10)
```

The checker constructs literal minimizers and maximizers; both are strong.
Their `(H,E32(2))` pairs are `(8281,8895)` and `(198959,206797)`.
The 1024 effective patterns give 984 distinct defect values. This is an exact
obstruction to paying the defect
down to zero by any number of these star edits with `{22,33}` still fixed.
It imposes no fixed-H constraint and does not rule out changing the in-port
set, changing Q, changing the protected prefix or changing the root.

## 8. A reusable growing-depth macro and its exact missing carry

With no retained ports, T1 is the ordinary unordered path-cover polynomial
`p_T(z)=sum_k pc_T(k)z^k`. Add x with `T→x`, then y with `x→y→T`, and forget
the two new names. This is the literal strong cyclic extension `C3(T,1,1)`.
Its exact positive operator is

```
C z^k = z^(k+2)+(2k+1)z^(k+1)+k(k+2)z^k+k(k-1)z^(k-1),
C = z^2+z+(2z^2+3z)∂+(z^2+z)∂^2.                   (T11)
```

It follows either by the two birth rules or by the older THM-3121 C3 kernel.
In particular `H(C3(T,1,1))=3H(T)+2pc_T(2)`. This is a reusable mechanism,
but the carry is a growing path-cover count rather than the desired 2 or4.

The coefficient of `pc_T(m+1)` in H after m such macros is
`(m+1)!m!`: its unique degree descent takes the `k(k-1)` term every time.
This proves why a fixed truncation of component degree is unsafe for arbitrary
future words, even though a prescribed finite batch has the sharp finite
lookahead bound of §3. It does not prohibit other nonlinear encodings.

There is also an exact factorization of the inherited adaptive reset. Retain
an old physical root w; after the cyclic macro, reverse x against every old
vertex except w. Now x has sole old in-port w and y has sole in-port x.
The resulting two-birth construction is the old N-then-R reset from
[`ADAPTIVE.md` (2.3)](../../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md):

```
H_final = H+2E_w(2) = 3H+2b_w,
E_x(2)=H_final.                                     (T12)
```

The x-star flip operator cancels the uncontrolled carry from T11 and replaces
it with the precisely transported rooted defect. If an actual source has
`b_w=1` or2, T12 gives its desired +2 or+4 child and a zero-ready marked
output. For a strong source this is again strong. The complete numerical
profile is checked, not just H. This inherits the earlier construction; the
new contribution is its explicit integral operator factorization and its
compatibility with the growing-port state.

**Still open:** a source-preserving rule that finds the required defect1 and2
ports, or an equivalent richer return condition, after every branch. Zero-ready
output alone does not supply that rule. A small formal vector satisfying the
carry equation is not automatically realized by a tournament.

## 9. Relation to the incoming multi-in-port construction

The incoming `session181_multiport_defect_20261104/multiport/PROOF.md` at commit
`b38b531963` evaluates an attachment with in-port set W by
`H=H(A)(h-J_W)+P_A(2)F_W` and root defect `H(A)J_W`. T1 retains its fields:
`F_W=sum_(w in W)E_w(2)` and `J_W=sum_(u,v in W)J(u,v)`, with the joint
slot disjointness already enforced by component pairing. A nonempty transitive
prefix can be grown by uniformly dominating births, forgetting a prefix name
only when it will no longer be distinguished. The final root birth has the
named incoming cut W and dominates the anonymous remainder.

Thus the formula is a readout of a compositional state, rather than a scalar
state reset. Its literal512425 graph and the earlier zero-ready512425 graph
remain different physical sources. Equality of their H values does not
identify their port states, charges, or native guards.

## Reproduction and status boundary

```
python 05-knowledge/results/multiport_defect_grammar_20261103/transport/check_transport.py
python -O 05-knowledge/results/multiport_defect_grammar_20261103/transport/check_transport.py --output 05-knowledge/results/multiport_defect_grammar_20261103/transport/verification_optimized.json
```

Every check uses explicit failure tests. The small universe is every labeled
tournament through order4, every active subset, every permitted uniform bulk
birth, every active flip, every port forget and converse. Independent
permutation-and-cut enumeration compares entire endpoint-profile vectors.
The receipt also checks commuting joins, all two/three-slot readouts,75 cyclic
macros and285 full adaptive reset factorizations. The B14 control uses an
independent endpoint subset DP and explicit two/three-factor convolutions,
with the actual source bytes pinned. It checks every flip-cube corner and
all1024 effective star patterns. No conjectural spectrum exclusion or
isomorphism filter is used.

The count state transports exact integers and actual joint ports. Positivity
is proved for births/forgetting and for edge-flip results on actual source
states. Arbitrary signed formal vectors are not declared realizable. Strongness,
native complement guards and finite-field phase conversions must still be
checked against the retained source word or literal adjacency.
