# Flat blue presentations above nonflat reversal markings

Session93, 2026-10-29. **PROVED elementary synthesis; FINITE-EXACT controls.**
Source pin `59cb22d0e93b6c364c4ec5d2954d9fe14ed468d8`.
No external-priority claim or black-core connectivity claim.

## Inheritance / concept board

Anchor: actual black/blue native lines. Niche: coprime fixed-orbit descent.
Wildcard: a torsor makes relative port position explicit.
Live board: marking transport; centralizer holonomy; blue path presentations;
native partner operation; joint named-port alignment.

The closest proved mechanisms are Session90's
`carriers/TRANSPORT_AND_HOLONOMY.md` and `oddness/FIXED_ORBIT_DESCENT.md` in
`05-knowledge/results/session90_even_graphs_odd_descent_20261029/`.
The canonical hostile pair for *this* question is W9 versus Paley7 below.
The corrected near miss is calling a marking triangle a black native move.
The least-used sidecar is the actual path, as opposed to its fixed-path count.

Exact native definitions come from:
- `01-canon/theorems/THM-643-blue-black-line-parity-structure.md`, with its
  October20 corrections overriding the historical pure-blue/leaf extrapolation;
- `01-canon/theorems/THM-644-gridsym-fiber-law-and-anti-redei.md`; its old
  conjectural parity status is superseded by `THM-647-anti-redei-theorem.md`;
- `01-canon/theorems/THM-4819-regular-tournaments-under-native-moves-the-transitive-slab-exact-distance-and-the-prescribed-black-core.md`.

## 1. Types before analogies

Let T be an odd-order self-converse tournament, A=Aut(T), and I its set of
anti-involutions. Let X be its ordered Hamilton paths. Automorphisms act freely
on X because a permutation fixing an ordered full path fixes every vertex.
For s in I define

    X_s = {P : s(P)=reverse(P)},     C_s=C_A(s).

Let rho_P be the positional reversal permutation of P. The ordinary
presentation [P] in X/A is **blue** precisely when rho_P is an anti-automorphism;
otherwise it is **black**. This is the native grid-symmetry definition, not an
arbitrary colouring. All paths belonging to X_s are blue, but a blue ordinary
orbit contains paths with several conjugate reversal marks.

The inherited odd coprime descent identifies

    B(T) = (X/A)^(converse)  ≅  X_s/C_s.

Consequently H=|X|, H/|A|=|X/A|, h_s=|X_s|, and b=|B(T)|=h_s/|C_s|
are four distinct statistics. The global blue labelled-path count is |I|h_s,
not h_s. Different paths have unique positional reversal, so the X_s are disjoint.

For a directed Hamilton path P, the **native partner** K_P(T) reverses every
arc not on P and preserves all path arcs. P is still a Hamilton path of its
partner. The pair of presentations is the native line; its colour is blue iff
rho_P is anti for T. If blue, rho_P is also anti for K_P(T), since the set of
unordered path edges is invariant under rho_P. Thus both line endpoints are blue.

## 2. Quotient-flat transport theorem

**Theorem.** If g in A satisfies gsg^-1=t, its map P↦gP induces a bijection

    F_ts : X_s/C_s → X_t/C_t

independent of the chosen g. Moreover F_ut F_ts=F_us.

**Proof.** The intertwining equation sends X_s to X_t and conjugates C_s onto
C_t, proving well-definedness. If g' is another transporter, g'=gc for c in C_s;
therefore g'P and gP are in the same C_t-orbit. The product of two transporters
is a transporter, giving strict composition. Equivalently every map identifies
its input with exactly the same blue ordinary orbit in X/A. ∎

For the inherited pair-local convention g(t,s)=sqrt(ts), the marking triangle

    h(u,t,s)=g(u,s)^-1 g(u,t)g(t,s)

lies in C_s. Therefore it is invisible on B(T). But if h≠1, it moves **every**
P in X_s, by freeness of the path action. Thus nonflat named-path transport and
flat blue-presentation transport coexist without contradiction.

This does not say a particular marking transport convention is forced: a
chosen global base gives flat transport even on named paths, as Session90
already explains. The obstruction concerns the specified pair-local convention.

## 3. Native-operation consequence, and its precise limit

For every g in A and every Hamilton path P (blue or black),

    K_{gP}(T) = g · K_P(T),

where g· means relabelling the tournament. This follows because g carries the
set of path edges of P to that of gP, and relabelling commutes with reversing
the complement. Positional reversal conjugates as rho_(gP)=g rho_P g^-1, so
colour is preserved as well.

**Corollary.** Every marking triangle preserves the source presentation and
native partner **isomorphism class**. It cannot turn a blue line black or supply
a black move from a new class. Literal labelled partners can differ.

Even an automorphism transporting a black path merely changes the representative
of that same black presentation; it does not perform the native move. The
partner may have a different automorphism group or mark centralizer from T.
Nothing here identifies those groups.

Hence the prescribed black-core theorem is complementary, not explained by
marking holonomy. THM-4819(D) supplies an actual word of m black native moves
from every regular T on 2m+1≥7 vertices to either prescribed singleton role plus
a strong core. Its named-start/cycle construction changes adjacency. Our
marking loop leaves the tournament isomorphism class unchanged. A claim that
holonomy generation alone proves black-core connectivity fails at this type
boundary. The general connectivity question is not settled here.

### A positive black-core certificate bridge

There is nevertheless an exact way to reuse a black-core certificate. Suppose
T_0=T and T_(i+1)=K_(P_i)(T_i) form a legal native word. For any permutation g,
relabel **every** T_i and P_i by g. The preceding identity gives a legal word
of the same length and colour sequence, preserving strongness, and carrying a
prescribed singleton v and its role to gv with the same role. For a marking
holonomy g in C_s the source is again T.

In odd order, s has a unique fixed vertex r, and every g in C_s fixes r:
s(gr)=g(sr)=gr. Therefore all holonomy translates of a prescribed-r black-core
word retain **that same named singleton**, not merely its vertex orbit. For a
noncentral named vertex v one must instead retain gv or restrict to its
stabilizer C_s∩A_v. Distinct g give distinct full path-word certificates when
the word is nonempty, since the first path has trivial stabilizer.

Thus holonomy does transport actual black-core certificates without changing
the central pole, but it cannot replace finding the certificate or change its
length. The first path-start guard and all later named starts must travel with
the same g; individually rematching them is not justified. This is the precise
positive connection to THM-4819's prescribed-singleton construction.

## 4. Exact holonomy sheet count

Let L_s≤C_s be the subgroup generated by all fixed-base triangle errors
h(u,t,s). It need not be assumed normal. Then

    X_s/L_s → X_s/C_s

has exactly [C_s:L_s] points over every blue presentation. Each blue fibre in
X_s is a free C_s-torsor of cardinality |C_s|; its L_s-orbits all have cardinality
|L_s|. This proves the assertion without normality.

If L_s=C_s, the holonomy orbits on fixed paths are exactly blue presentations.
If L_s=1, holonomy preserves every named fixed path, but the C_s quotient may
still forget information. Thus trivial holonomy does not mean trivial centralizer.
The recursive-family assertion L_s=C_s is handled separately in Session93;
this note uses only the verified W9 instance and states the general consequence
conditionally.

## 5. A port coordinate cannot be discarded merely because blue transport is flat

Fix one blue presentation and its C_s-torsor Y⊂X_s. Let D be any C_s-set of
additional data, for example a vertex or tuple of named ports. The joint
quotient (Y×D)/C_s is, after choosing p0 in Y, in bijection with **D**, not
usually D/C_s: write uniquely P=c p0 and send [(P,d)] to c^-1 d.

For D=V(T) the invariant is the position P^-1(v). Moving both P and v preserves
it. Moving P around a holonomy loop while keeping a named vertex v frozen can
change it. Quotienting the path and the port separately loses this alignment.
This elementary torsor calculation explains exactly when a bit is free: only
when the later operation is invariant under the simultaneous residual action.

In W9 the known triangle h=(0,1,2,4,5,3,8,6,7) moves port3 to4, while fixing
the anti-involution. For example the fixed path (1,3,4,7,0,5,8,6,2) becomes
(1,4,5,6,0,3,7,8,2); frozen port3 moves from position2 to position6
(one-based). Their literal labelled native partners differ but are isomorphic
by h. Every full fixed path detects a nonidentity holonomy. The test
“current vertex is3” in Session92's folded observer is an explicit operational
version of this issue (16 minimal weighted-partition states become35), though
that partition size is not deduced from the torsor calculation.

Session92's actual strong-child ordered-run control is sharper in a different
direction: all individual cursor responses agree, but the two path-dependent
native partners have H=2,483,993 and2,491,877 (or2,899,675 and2,906,091 with all
three children strong). Those tuples are **not related by one common child
automorphism**. Our equivariance theorem therefore does not identify their
partners. The exact distinction is simultaneous transport versus separately
chosen pointwise transports, not a failure of equivariance.
Sources: `session92_joint_order_carriers_20261029/runs/RECOVERED_RUN_WITNESS.md`
and `observers/FOLDED_OBSERVER_COSTS.md` in the results directory.

## 6. Hostile controls and reproducibility

`check_blue_holonomy.py` is standalone standard-library Python, with explicit
checks active under optimization. Normal and `python3 -O` outputs agree, with
64,018 checks. It constructs W9=C3[C3,C3,C3], uses its inherited full wreath
Aut description (and verifies each supplied permutation); independently searches
all7! permutations for Paley7; enumerates all Hamilton paths; checks every
fixed-base marking triangle and its action on every fixed path; checks partner
relabel equivariance and colour for all paths under one nonidentity automorphism. It also replays a four-step actual native word
under all centralizer elements, checking each intermediate colour, the native
involution, and the fixed named pole (these test words are not claimed to be
black-core geodesics).
The structural proof of fullness of W9's wreath automorphism group is inherited
from Session90 `carriers/ITERATED_TRIANGLE_FAMILY.md`, not re-proved by this test.

|T|H|A|C_s|I|h_s|blue presentations|black presentations|L_s|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|W9|3159|81|3|27|33|11|28|3|
|Paley7|189|21|3|7|9|3|6|1|

For W9,432 of729 fixed-base triangles are nonidentity. For Paley7, all49 are
identity. In both examples the set of triangle values happens already to be
a subgroup, verified explicitly, but no universal closure assertion is made.

The two controls refute both tempting shortcuts:
- flat holonomy does not imply a pure-blue tournament (Paley7 has6 black presentations);
- full centralizer holonomy does not imply a pure-black tournament (W9 has11 blue presentations).

No output counts target isomorphism classes or merges parallel native lines.
No general improvement to A000568 enumeration follows solely from these counts.
What is reusable is a safe implementation rule: memoize native partner class
by the **joint ordinary path orbit**, and transport all named decorations with
its witness. Do not replace that joint orbit by marginal port responses.
