# Signed deletion cards, exact root forgetting, and A000568

**Status: PROVED (E1–E4; independent audit PASS).** The results E1–E4 have self-contained proofs below. [Independent proof audit](../audit_even/PROOF_AUDIT.md). Classical equinumerosity and Burnside enumeration are CITED; the order-six census is FINITE-EXACT. This is an unrestricted, state-rich augmentation recurrence, not a closed scalar recurrence in earlier A000568 values, a new natural bijection, or a claim of literature priority.

## 1. Inheritance and the exact meaning of “even”

An even graph here is a finite simple graph for which **every automorphism reverses an even number of reference-oriented edges**. It need not have even vertex degrees. Let

\[
 \epsilon_X(g)=(-1)^{\#\{ij\in E(X):i<j,\ g(i)>g(j)\}},
 \qquad g\in\operatorname{Aut}(X).
\]

This is a character, independent of the chosen reference orientation. A graph is sign-even when this character is trivial, and sign-odd otherwise. All graph and tournament counts below mean **isomorphism types on exactly n vertices**; tournaments are not identified with their converse unless they are isomorphic.

Royle, Praeger, Glasby, Freedman and Devillers, *Tournaments and even graphs are equinumerous*, Theorem 1.1, prove equality of these unlabelled counts; their Lemma 4.1 proves the character property. Their Theorem 3.3 gives the classical tournament Burnside sum. We use [arXiv:2204.01947v2](https://arxiv.org/pdf/2204.01947v2), §§1, 3, 4; the tiny orders 0, 1, 2 are also checked directly. Write their common sequence as `a_n=e_n=A000568(n)`, with `a_0=1`.

The distinctions have immediate hostile examples:

* At n=2 there are two labelled tournaments but only one labelled sign-even graph, the empty graph. Both have one isomorphism type. Thus count equality does not supply a labelled or equivariant bijection.
* At n=3, `P3` is sign-even although its endpoints have odd degree. `K3` has even degrees but is sign-odd: a reflection reverses an odd number of its three edges. Sign-even and Eulerian are different predicates.
* The sign-even graph `P3` has a sign-odd deletion card `K2`. Restricting augmentation to even parents loses valid even children.

The immediate inherited proof is [signed_root_augmentation.md](../../sign_aware_augmentation_20261020/signed_root_augmentation.md), particularly Theorems 1–3 and §6, with [ROOT_REVIEW.md](../../sign_aware_augmentation_20261020/ROOT_REVIEW.md). It already gives the relative-even root compiler, the first attachment-fibre obstruction, and the safe recurrence

\[
 e_N=\sum_{[G]}\sum_{[S]\text{ relatively even}}
       \frac{\mathbf 1_{G+v(S)\text{ globally even}}}{r(G+v(S))},
 \tag{I}
\]

where `r(X)` is the number of vertex orbits. Its explicit remaining sidecar is the sign of **root-moving transporters**. We supply that sidecar using marked deletion cards and prove that its cancellations localize under disjoint union. This is a refinement of (I), not a first unrestricted counting formula.

Other recovered boundaries are [THM-514, unit-tail core kernel](../../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md), which already supplies a cycle-profile all-order formula; [THM-4169, prime-parent augmentation](../../../../01-canon/theorems/THM-4169-prime-parent-one-vertex-augmentation-and-quartic-johnson-transfer.md), which retains root/card incidence weights; and the [primitive Hopf/Lie packet](../../primitive_hopf_lie_bijection_20261020/PRIMITIVE_HOPF_LIE_BIJECTION.md), which already relates connected even graphs to Lyndon words in strong tournament types. Its global matching is not a local attachment bijection. The October 20 [symmetry correction](../../c15_symmetry_counterexample_20261020/C15_ORDER_TEN_COUNTEREXAMPLE.md) prevents promotion of count equality to arbitrary symmetry-preserving coupling.

The concept board is: sign character / attachment orbit / marked deletion card / root-moving transport / connected-component composition. The anchor is exact unrooted enumeration; the niche is localization of the odd-root defect; the wildcard is the rank-one signed transport matrix in E2. Their common missing coordinate is the sign on a root-moving map.

## 2. E1: the order-n marked-card recurrence

Fix `N=n+1>=1`. For each unlabelled simple graph `G` of order n choose one representative on `[n]={0,...,n-1}`, and put `A=Aut(G)`. Choose one representative `S` of each A-orbit on subsets of `[n]`. Adjoin a distinguished last vertex `v=n` adjacent exactly to S, obtaining X. Put

\[
 K=A_S=\{k\in A:k(S)=S\}.
\]

For each `w∈V(X)`, let `F_w(G,S)` be the set of graph isomorphisms

\[
 f:G\longrightarrow X-w,\qquad f(S)=N_X(w).
 \tag{1}
\]

Extend f to the bijection `f̂:V(X)→V(X)` by `f̂(v)=w`. The marked-neighbourhood condition in (1) makes f̂ an automorphism of X. Define

\[
 Q(G,S)=\frac1{|K|}\sum_{w\in V(X)}\ \sum_{f\in F_w(G,S)}
            \epsilon_X(\widehat f).
 \tag{2}
\]

**E1.** For every N≥1,

\[
 \boxed{\quad a_N=e_N=\frac1N
       \sum_{[G]:|G|=N-1}\ \sum_{[S]\in\operatorname{Aut}(G)\backslash2^{V(G)}}
                         Q(G,S).\quad}
 \tag{3}
\]

For each attachment, Q is exactly the size of the root's automorphism orbit when X is sign-even, and is zero when X is sign-odd. In particular Q is a nonnegative integer despite its signed definition.

**Proof.** Extension by `v↦w` identifies `F_w` with all `g∈Aut(X)` satisfying `g(v)=w`. These fibres partition `Aut(X)`. The stabilizer of v is exactly K, extended by fixing v. Since v is last in the reference order, no edge incident to v reverses under this stabilizer, so its edge character is the restriction of `ε_G`.

Consequently

\[
 Q(G,S)=\frac{1}{|\operatorname{Aut}(X)_v|}
          \sum_{g\in\operatorname{Aut}(X)}\epsilon_X(g).
 \tag{4}
\]

If the character is trivial this equals `|Aut(X)|/|Aut(X)_v|=|Orb_X(v)|`. Otherwise multiplication by any negative-sign automorphism pairs positive and negative group elements, and the sum vanishes.

Pairs `([G],[S])` classify rooted graph isomorphism types: delete the distinguished root for the forward map; root-preserving isomorphisms are precisely parent isomorphisms carrying the neighbourhood. For a fixed unrooted X, those rooted types correspond to its vertex orbits. If X is sign-even, their orbit sizes sum to N. If X is sign-odd all their charges are zero. Sum (4) and divide by N. Finally use the cited equality `a_N=e_N`. ∎

This is an actual one-order computation: (1) involves graphs of order n and an explicit neighbourhood mark. Its implementation may enumerate card isomorphisms directly or use an order-n graph catalogue with isomorphism transporters. It needs no canonical representative or automorphism oracle for an order-N child. It still has full isomorphism-search cost; a card implementation in effect recovers the relevant child automorphisms. There is no claimed complexity improvement or recurrence using only the numbers `a_0,...,a_(N-1)`.

## 3. E2: choice-independent signs and composition

Call the rooted attachment relatively even if `ε_G|K=1`. If it is not, every nonempty `F_w` already has zero sign sum: it is a coset of K and K has as many positive as negative elements. Such attachments may be omitted from (3).

For a relatively-even attachment and `w∈Orb_X(v)`, choose any one map in `F_w` and set

\[
 \eta_v(w)=\epsilon_X(\widehat f)\in\{+1,-1\}.
\]

**E2.** This does not depend on the chosen marked-card isomorphism, and

\[
 Q(G,S)=\sum_{w\in\operatorname{Orb}_X(v)}\eta_v(w).
 \tag{5}
\]

If X is sign-even, every term is +1. If X is sign-odd, exactly half of the reachable roots have each sign. Thus the filtered form of (3) needs one transporter per reachable root, rather than every transporter.

**Proof.** Two automorphisms taking v to w differ by an element of the v-stabilizer, whose sign is +1. Every nonempty fibre has |K| maps of that same sign, giving (5). If X is odd, its index-two sign kernel has exactly two orbits of equal size inside `Orb_X(v)`: the stabilizer is contained in that kernel. These are the positive and negative transporter classes. Equivalently the balanced assertion follows directly from (4). ∎

There is also exact composition, not just a root count. On any relatively-even vertex orbit O define `η_u(w)` using any automorphism taking u to w. Then

\[
 \eta_u(w)\eta_w(t)=\eta_u(t),\qquad \eta_u(u)=1.
 \tag{6}
\]

Indeed signs multiply under composition and the choices are immaterial. Fix `b∈O` and put `s_w=η_b(w)`. The matrix `B_(u,w)=η_u(w)` on O is `s s^T`, so it has rank one over Q and satisfies `B²=|O|B`. Changing b changes s by a common sign. For an odd X, `sum s_w=0` and therefore `B·1=0`; for an even X, B is the all-ones matrix. This is a signed card-transport identity; it is not a tournament adjacency matrix or a label-free bijection.

## 4. E3: the exact disjoint-union localization

We first recall the inherited component-even lemma. A graph is sign-even iff all its connected components are sign-even. To see the nontrivial direction, reference-orient corresponding copies of an isomorphic connected component identically. Permuting whole component copies then has sign +1. Every remaining automorphism acts within components, so its sign is the product of the component signs. Reference independence removes the orientation choice.

**E3.** A rooted graph `(X,v)` is relatively even iff every component not containing v is sign-even and the root stabilizer inside the component containing v is sign-even. Consequently, if X itself is sign-odd but `(X,v)` is relatively even, then:

1. X has exactly one sign-odd connected component C;
2. C occurs with multiplicity one and contains v;
3. every other component is sign-even;
4. `Orb_X(v)=Orb_C(v)`, and all transporter signs in (5) are exactly the signs of their restrictions to C.

**Proof.** Any odd automorphism of a component not containing v extends by identity to an odd automorphism fixing v. An odd root-fixing automorphism of its own component does the same. These prove necessity. Conversely an automorphism fixing v must preserve its component. On the rest it consists of permutations of isomorphic components and internal component automorphisms, all with sign +1 by the preceding lemma. The sign is therefore just its root-component sign, proving sufficiency.

If X is globally odd, at least one component is odd. Necessity places every such component at v, so there is only one. A second copy would be another odd component not containing v, which is forbidden. Thus no automorphism can move the root component to a different component. Every root transporter restricts to C, the rest has sign +1, and every automorphism of C extends by identity. This gives the orbit and sign statements. ∎

Hence cancellation of the negative branches in E2 is confined to one odd connected component. Adjoining any multiset of sign-even components leaves that local cancellation unchanged. Arbitrary graph parents still have to be admitted: attachment may merge several previously separate components, and evenness is not hereditary.

## 5. E4: a factorization of the rooted defect

Let `o_n` count relatively-even rooted isomorphism types whose underlying graph is **connected and sign-odd**. Let `d_n` count relatively-even rooted types of arbitrary sign-odd graphs. Set `o_0=d_0=0`, `E(z)=sum_(n>=0)e_n z^n`, `O(z)=sum_(n>=1)o_n z^n`, and `D(z)=sum_(n>=1)d_n z^n`.

**E4.** There is a size-preserving bijection

\[
 \{\text{relatively-even rooted odd graphs}\}
 \longleftrightarrow
 \{\text{relatively-even rooted connected odd graphs}\}
 \times\{\text{unrooted even graphs}\},
\]

given by taking the root component and its complementary components. Consequently

\[
 \boxed{D(z)=O(z)E(z),\qquad d_n=\sum_{j=0}^n e_j o_{n-j}.}
 \tag{7}
\]

**Proof.** E3 gives the unique decomposition and its inverse, disjoint union. An even complementary graph contains no component isomorphic to the odd root component, so no hidden exchange changes the rooted type. ∎

The inherited relative-even theorem says pointed tournament types satisfy `t_n^bullet=e_n^bullet+d_n`, where `e_n^bullet` counts vertex orbits of even graphs, not labelled vertices. Formula (7) therefore identifies exactly the connected source of the inherited root-forgetting defect. It is not the incorrect identity `t_n^bullet=N a_n`.

The finite-exact sequences through n=6, indexed from n=0, are

```
e_n = 1,1,1,2,4,12,56
o_n = 0,0,1,0,4,9,58
d_n = 0,0,1,1,5,15,75
t_n^bullet = 0,1,2,4,12,48,296.
```

These values do not assert a closed scalar formula for `o_n`.

## 6. Smallest failure boundaries

For N=2, the attachment producing `K2` has trivial root stabilizer, but the two possible roots have transport signs `+1,-1`. Its charge is zero. The empty graph `I2` contributes charge 2. Thus (3) gives `a_2=1` rather than confusing two relatively-even rooted objects with two globally even children.

For N=3, `P3` has a center orbit of size 1 and a leaf orbit of size 2. The center attachment has parent `I2`; the leaf attachment has the odd parent `K2`. Their normalized contributions are `1/3` and `2/3`. Including `I3`, the even-parent-only sum would be `4/3`, while the correct sum is 2. This is the first failed even-parent-only recurrence.

The neighbourhood mark in (1) is indispensable even in that example: when v is a leaf of `P3`, its deletion card is `K2`. The two unmarked automorphisms of that card do not both extend by fixing v; only the one preserving v's singleton neighbourhood does. Dropping the mark includes a map that is not a child automorphism, so its supposed character/cancellation is mistyped.

No count equality supplies a particular graph attached to each tournament, preserves all automorphism actions, or intertwines parent deletion. The inherited rooted-fibre obstruction is already `4≠3` at the order-two tournament parent and the sole even graph parent. Our signed all-parent operator repairs enumeration, not that nonexistent fibre-preserving pairing.

## 7. Classical formula and verification

For comparison, the cited Theorem 3.3, grouped by odd cycle lengths, gives

\[
 a_n=\sum_{\lambda\vdash n\text{ with odd parts}}
 \frac{2^{\sum_i(\lambda_i-1)/2+\sum_{i<j}\gcd(\lambda_i,\lambda_j)}}
 {\prod_j j^{m_j}m_j!}.
 \tag{8}
\]

This is classical Burnside enumeration; the repository already refines it in THM-514. Formula (3) is a different, card-based implementation of unrooted counting with explicit cancellation. It does not make (8) new or remove its state requirements.

[check.py](check.py) uses only the Python standard library and imports no repository producer. Its universe is all 209 graph types and 77 tournament types of orders 0 through 6, all their automorphisms, and all 663 rooted attachments to graph parents of orders 0 through 5. The card routine enumerates order-n marked isomorphisms directly. An independent full-child automorphism route checks each individual transport fibre, the stabilizer restriction, its charge, and the sum per unrooted child. Edge sign is checked again using self-paired edge orbits. Direct tournament relabelling, direct sign-even graph enumeration, and the partition formula (8) give three counting paths. Component localization and (7) are checked on the complete small universe; classical (8) is also evaluated through order 12, without claiming an exhaustive graph census there.

Reproduce from this directory:

```
python3 check.py > receipt.json
python3 -O check.py > receipt_optimized.json
cmp receipt.json receipt_optimized.json
```

The two receipts are byte-identical: **17,990 explicit checks, PASS**. `-O` does not remove the checks. See [SOURCES.md](SOURCES.md) and [SOURCE_MANIFEST.json](SOURCE_MANIFEST.json) for inherited source pins and the primary-source reading boundary.

The source/target contract is precise: all parent graph types plus neighbourhood orbits map to rooted child graphs; marked-card maps carry root-transport signs; weighted forgetting preserves the unrooted even-graph count and hence A000568. Dropping the neighbourhood destroys extension validity; dropping transporter signs destroys cancellation; dropping odd parents loses `P3`; retaining only prior scalar counts loses the physical state needed to run (3). The remaining target is a smaller sufficient state or a constructive correspondence that also intertwines tournament composition. Neither is established here.
