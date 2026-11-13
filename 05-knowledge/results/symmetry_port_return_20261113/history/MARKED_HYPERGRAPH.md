# A marked cyclic-triangle congruence for arbitrary root-child programs

**Status: PROVED**, using the inherited cyclic-triangle path-cover invariant
in [THM-4803](../../../../01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md).
The finite receipt is a control, not the proof. This is an operational
extension of that invariant; it is not a new claim about the classification
of cyclic-triangle fibres.

## 1. Object, allowed operations, and exact theorem

For a tournament `T` on a finite named set `V`, let `C(T)` be the
3-uniform hypergraph whose edges are its directed cyclic triangles. For a
named vertex `r`, its **star** is the full list of orientations between
`r` and `V\{r}`. A root-child birth at `r` adds a new vertex `z` with

```
r -> z -> V\{r}.                                      (1)
```

Fix a set `P` of initially permitted ports. At every subsequent step, the
chosen root may be any member of `P` or any vertex born in an earlier step.
We retain names and may keep any number of these available ports. This
contains the usual `R/N` programs, but also arbitrary interleaving among
several retained ports.

**Theorem.** Suppose two tournaments `T,T'` on the same named set satisfy

```
C(T)=C(T'),
star_T(r)=star_T'(r) for every r in P.                  (2)
```

Apply any common finite root-child program of the specified kind. At every
stage:

1. the resulting labelled cyclic-triangle hypergraphs agree;
2. the full stars of every available port agree;
3. the entire ordinary path-cover profiles, in particular Hamiltonian-path
   counts, agree.

The same conclusion holds for an adaptive deterministic program whose next
root/stop decision is a function only of the common names, its previous
choices, and previously observed path-cover profiles. Randomized programs
can be coupled by using the same random choices.

**Proof.** In a birth (1), an old triple is unchanged. A triple containing
`z` but not `r` is transitive, since `z` beats its other two vertices. The
remaining new triples are exactly

```
{r,z,v} is cyclic  iff  v -> r.                        (3)
```

Thus the new hypergraph is determined by the old hypergraph and the chosen
root's star. Stars of old available ports acquire the same one new arc in
both tournaments; the new vertex's entire star is prescribed by (1).
This proves (1) and (2) of the theorem by induction. The inherited theorem
THM-4803 says that the full ordered path-cover profile is constant on each
labelled cyclic-triangle fibre, giving (3). An adaptive program receives
the same observations, so makes the same next choice and remains inside
the induction. QED.

No assertion is made for a later choice of a previously unmarked old vertex:
its star was not part of (2). Nor does the theorem cover arbitrary newborn
cuts. These two ways to release the protected interface are material.

## 2. The complete four-coordinate response is already fixed

For a root `r`, use the established conventions

```
h = H(T),                       q = H(T-r),
s = number of spanning paths beginning at r,
F = ordered two-slot covers whose first slot is nonempty and ends at r;
    the second slot may be empty.
```

The inherited root-retaining birth map is

```
R(h,s,q,F)=(h+F-s, q, q, 2F-q).                       (4)
```

Let `h_R,h_RR` be the Hamiltonian counts after one and two births at `r`.
Then

```
h_R  = h+F-s,
h_RR = h+3F-s-2q,
F = q+(h_RR-h_R)/2,       s = h+F-h_R.                (5)
```

Under (2), `h` agrees by THM-4803 and `q` agrees by applying the same
invariant to the induced deletion. The theorem gives common `h_R,h_RR`,
so (5) gives common `s,F`. Hence the marked state, every continuation gcd
derived from that state, and the readiness defect `q-s` are fixed. This
is a finite observation proof of the state consequence; it does not assume
a new theorem about marked modular decomposition.

## 3. Indegree one removes even the star hypothesis

**Corollary.** Suppose `C(T)=C(T')` on the same named set, and the same
named `r` has indegree one in both tournaments. Then the full marked states
at `r` agree, even if its unique predecessor changes. Consequently their
unit-readiness defects and every `R/N` continuation agree.

**Proof.** Let `p` and `q` be the unique predecessors. If `p=q`, the stars
agree and Section 2 applies. Suppose `p!=q`. Every cyclic triangle through
`r` contains `p` in `T` and `q` in `T'`, so the common hypergraph contains
at most one such triangle, namely `{r,p,q}`.

If it contains none, `p` is a global source in `T`: otherwise `v->p`
would give `p->r->v->p`. Likewise `q` is a global source in `T'`. Thus the
common hypergraph has no edge through either `p` or `q`.

If it contains `{r,p,q}`, the same argument shows that `p` dominates every
vertex other than `q` in `T`; hence every cyclic triangle through `p`
contains `q`. In `T'`, every cyclic triangle through `q` contains `p`.
The common hypergraph therefore has no edge containing exactly one of
`p,q`. In both cases the transposition `(p q)` preserves that hypergraph.
Apply it to the names of `T'`. It fixes `r`, changes its unique predecessor
from `q` to `p`, and keeps the common hypergraph. Section 2 now proves
equality of marked states. Renaming vertices other than the root does not
change the state. QED.

Thus a same-count unit-ready rearm which retains the named root and
indegree one must change its cyclic-triangle hypergraph. This statement
allows arbitrary neutral changes of the old graph and predecessor; it is
stronger than the fixed-star special case. It does not prohibit neutral
repairs such as Section 5, which make precisely that hypergraph change.

## 4. What symmetry breaking can and cannot accomplish

The inherited October 26 example

```
A=C3 -> 1,       B=1 -> C3,
U=C3[A,A,A],     U'=C3[B,A,A]
```

has the same labelled cyclic-triangle hypergraph and full cover profile,
but automorphism orders `81` and `27`. It is recorded in
[the corrected fibre result](../../tournament_composition_classification_20261026/fibre/RESULT.md).
Put the marked root in one of the two unchanged leaves. Its entire star is
unchanged, so our theorem proves that **every allowed future port program
has the same response**, despite the global automorphism change.

This is compatible with the current ports lane's positive construction:
one can break a global leaf permutation while preserving an already
primitive fixed port. That lane uses equality of the substituted module's
cover profile, which is a different sufficient congruence and need not
preserve the labelled cyclic-triangle hypergraph. Neither congruence can
improve the response at the protected port. Selecting another old port or
changing its physical star is a distinct operation and requires a separate
sidecar.

## 5. A sharp physical boundary from the five cores

The established same-root-star one-arc repair has rows

```
T =(22,8,2, 5,14),       (h,s,q,F)=(9,3,3, 9),
T'=(22,8,2,21, 6),       (h,s,q,F)=(9,4,5,10),
root=0.
```

Only the arc between vertices `3,4` changes. The root star and `h` are
unchanged, but the cyclic triangles change from

```
{013,034,123}  to  {013,123,134}.                      (6)
```

Its readiness defect changes `0 -> 1`, and the common first `N` count is
`15`, followed by `NR` counts `27` and `29`. This is inherited from
[ONE_ARC_H9](../../h17_exact_repairs_20261112/history/ONE_ARC_H9.md), not a
new example. It shows concretely how a successful port repair leaves the
congruence fibre: a root-star-preserving hypergraph trade releases the
necessary information. Equal scalar count, or even equal unmarked cover
polynomial, does not replace (2).

## 6. Historical transfer and exact controls

The source mechanism was not a claimed identification of LRC and
tournaments. July 25's
[THM-2290 context-selected kernel](../../../../01-canon/theorems/THM-2290-context-selected-colored-pair-kernel-is-hafnian-complete.md)
asks whether a quotient remains exact under every allowed context. July
29's
[THM-2894 residual semilattice](../../../../01-canon/theorems/THM-2894-unmarked-residual-semilattice-order-and-group-clutch-no-go.md)
proves that information lost by a context congruence cannot be restored by
more operations of the same type. Here the transported test is: retain
`C(T)` and port stars, then derive the exact update (3). The contraction is
the tournament path-cover invariant, not a hafnian or an owner-flow kernel.

`hypergraph_probe.py` enumerates every labelled tournament of orders `2..5`
at root `0`; every common `(C(T),star(0))` fibre has a single marked state.
It also tests the indegree-one corollary after discarding the star from the
fibre key.
It also checks the birth hypergraph identity directly, a small complete
multiport program universe, the inherited automorphism-changing example,
and the H9 escape boundary. Reproduce with

```
python3 -B 05-knowledge/results/symmetry_port_return_20261113/history/hypergraph_probe.py
```

The theorem's all-order quantifier comes from (3), not from this finite
universe. There is no claimed converse: the marked state may agree for
different labelled cyclic-triangle hypergraphs.
