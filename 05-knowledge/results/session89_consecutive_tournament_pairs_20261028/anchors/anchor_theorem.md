# Exact orientation-anchor boundary for matrix / tournament-pair codecs

**Status:** PROVED (elementary theorem and formulas); FINITE-EXACT (reported tests).
**Session:** 89, 2026-10-28.
**Read-only inheritance pin:** `03408f84ac43fcdcfb0887e5d6bc1fc8792415cf`.
**Scope:** all binary square matrices, with no matrix filter; ordered natural
consecutive tournament pairs with a common labelled core and a fixed extra root.
No claim of a new literature theorem, a canonical anchor, or a GI lower bound.

## 1. Objects and exact theorem

Let `V` be a finite `n`-element set, `E=binom(V,2)`, and let `H≤Sym(V)` be a
**faithful permutation subgroup**. Define

- `X(V)={0,1}^{V×V}`, acted on by simultaneous row/column relabelling:
  `(hM)_(hi,hj)=M_ij`;
- `Y(V)=Tour(V) × Tour(V∪{*})`, with the second tournament's extra vertex `*`
  distinguished and fixed. Both tournaments are relabelled by the same `h` on
  their shared core `V`. The pair is ordered.

Both labelled sets have size `2^(n²)`, because
`binom(n,2)+binom(n+1,2)=n²`. Equality of labelled sizes alone says nothing about
these actions.

**Theorem.** The following are equivalent:

1. An `H`-equivariant bijection `X(V)→Y(V)` exists.
2. An `H`-equivariant map `X(V)→Y(V)` exists, even without injectivity.
3. There is an `H`-invariant tournament on `V`.
4. No element of `H` interchanges the two endpoints of an unordered pair.
5. `|H|` is odd.

These include `n=0,1`, where `H` is trivial.

**Proof.** (1) implies (2). The zero matrix is fixed by all of `H`. Its image
under (2) is an `H`-fixed pair, so its first tournament proves (3). An invariant
tournament cannot have an element reverse an edge, giving (4). If `|H|` were
even, Cauchy's theorem would give a nonidentity involution. In a faithful
permutation action that involution contains a 2-cycle, contradicting (4).
Thus (4) implies (5).

Conversely, suppose `|H|` is odd. For each `H`-orbit of unordered pairs, select
one pair `{u,v}` and orient it `u→v`, then propagate this orientation by `H`.
A conflict would mean some element sends `(u,v)` to `(v,u)`; its order would be
even. Thus the propagation is well-defined and yields an invariant tournament
`A`. This proves (5) implies (3). The next construction proves (3) implies (1).

### Codec from a fixed invariant anchor

For every anchor arc `u→v` in `A`, set

```
T(u,v) = M_uv,
U(u,v) = M_vu,
U(*,i) = M_ii.
```

Here `T(u,v)=1` means that `T` has arc `u→v`; bit zero means the reverse arc.
Thus bit one means “follows the anchor,” consistently in both core tournaments.
Every matrix coordinate occurs exactly once, giving the inverse by reading the
same arcs. Since every `h∈H` preserves `A`, it preserves which matrix coordinate
is assigned to which tournament, and the codec commutes with `h`.

This is an explicit affine binary-coordinate bijection. It is not presented as
a canonically chosen codec. The anchor `A` is the required choice/sidecar.

### Nonfaithful actions

For an abstract finite group `G` acting through `ρ:G→S_n`, the correct condition
is that `|ρ(G)|` is odd, not that `|G|` is odd. An even kernel is invisible.
For example, an abstract `C2` acting trivially on two vertices admits both
anchors and an equivariant codec. The theorem's `H≤S_n` wording avoids this trap.

## 2. Exact number of anchors and fixed objects

Put

```
v_H = |V/H|,
e_H = |E/H|,
o_H = |(V×V)/H|.
```

Every unordered-edge orbit permits exactly two propagated orientations when
there is no reversal; choices on distinct orbits are independent. Consequently

```
|Tour(V)^H| = 2^e_H  if |H| is odd,
             0      if |H| is even.
```

More generally, the same formula works for any specified graph edge set if no
edge orbit contains a reversal; for a graph other than the complete graph,
this no-reversal condition need not force the whole faithful group to be odd.

For the complete graph and an odd group, each unordered-edge orbit splits into
exactly two directed-edge orbits. The diagonal orbits are the vertex orbits, so

```
o_H = 2e_H + v_H,
|X(V)^H| = 2^o_H,
|Y(V)^H| = 2^(2e_H+v_H).
```

For even `H`, the first count is still positive and the second is zero. This
is a fixed-set obstruction, not merely a mismatch of total orbit counts.
For odd `H`, the same formulas hold for every subgroup of `H`, in agreement
with the explicit `H`-set isomorphism.

### Efficient construction from generators

An explicitly enumerated group is unnecessary. Given `s` permutation
generators, impose the parity equations on unordered pairs (in any fixed
reference orientation):

```
x_(ge)=x_e XOR c(g,e).
```

Propagate these constraints on the graph with one vertex per unordered pair
and one generator-transport edge per pair per generator. A contradiction
means no invariant anchor. Otherwise each connected component contributes one
free bit; its two solutions give the two choices for that unordered-edge
orbit. A representative anchor and all free components are obtained in
`O((s+1)n²)` elementary visits, including initialization when `s=0`.
This is a **compact affine solution family**: one particular anchor and the
free orbit components. The code exposes this as `solve_anchor_constraints`;
its separate `anchors_by_generators` wrapper explicitly enumerates all `2^e`
anchors and has exponential output cost in general. Neither operation computes
`Aut(A)`, an arbitrary matrix automorphism group, or a GI canonical form.

A useful hostile is `A4` on four vertices, generated by `(012)` and `(013)`.
Both generators have odd order three, but their generated group has order
12 and permits no anchor. Testing individual generator orders is insufficient;
parity propagation catches the composite-word obstruction.

**Short certificates.** A positive certificate is the returned tournament:
a verifier checks that each supplied generator preserves every oriented edge.
For a negative certificate, retain predecessor edges and signed generator
labels during propagation. A parity conflict at a constraint closes a walk
from its component's root edge back to that same unordered edge, with odd
orientation reversal. Its signed generator word therefore swaps the root
edge's endpoints. Evaluating the word directly certifies this without listing
`H`; its resulting permutation has even order. The word has length at most
`2|E|−1` because it uses two tree paths and one conflict edge. This is a
polynomial-length certificate, not a claim of shortest-word optimality.

The implementation returns that word and independently verifies its action.
For `a=(012)`, `b=(013)`, the `A4` output is the three-letter word
`b^(-1), a, b^(-1)`, applied left to right. It evaluates to `(01)(23)` and
reverses the edge `{0,1}`. This ties the anchor test to the previous session's
short-certificate theme without importing a general isomorphism algorithm.
The odd-order conclusion is about the faithful permutation image, not an
abstract acting group's potentially invisible kernel.

## 3. Orientation torsor, reversal cocycle, and gauge

Fix any reference anchor `A`, not necessarily invariant. For `e={u,v}`, write
its anchor direction as `a_e→b_e`. An ordinary matrix has local coordinates

```
z_e=(x_e,y_e)=(M_(a_e,b_e), M_(b_e,a_e)).
```

A natural pair has coordinates `(x_e,y_e)` indicating whether each core
orientation agrees with `A`; its root coordinates are `d_i=U(*,i)`.

Define the reversal indicator `c_A(g,e)` to be one if the image under `g` of
the anchor arc on `e` is opposite the anchor arc on `ge`, and zero otherwise.
With `gh` meaning apply `h` first, then `g`, it satisfies

```
c_A(gh,e)=c_A(h,e) XOR c_A(g,he).
```

The two actions on the same local bit coordinates differ:

```
ordinary matrices: z'_(ge) = swap^(c_A(g,e))(z_e),
natural pairs:     z'_(ge) = z_e XOR (c_A(g,e),c_A(g,e)).
```

Both root/diagonal coordinates are simply permuted. Invariance of `A` means
`c_A(h,e)=0` for all `h∈H`, removing exactly this discrepancy.

### Coboundary and odd-order averaging

If anchor `B` is obtained by reversing `A` on a set encoded by `f∈F₂^E`, then

```
c_B(g,e)=c_A(g,e) XOR f_e XOR f_(ge).
```

Thus finding an invariant anchor is exactly trivializing this affine cocycle.
In destination-index convention the orientation torsor action is
`g⋆x=P_g x+t_g`, with `t_gh=P_g t_h+t_g`. For odd `H`, put
`b=sum_(h∈H) t_h` over `F₂`. Then

```
g⋆b = sum_h (P_g t_h+t_g) = sum_h t_gh = b,
```

where the first equality uses `|H|=1` in `F₂`. This proves directly that the
affine torsor has a fixed point; it is the elementary odd-order averaging
mechanism behind the coboundary statement. No Feit–Thompson theorem is needed.

### Changing the codec, not just its coordinate chart

Let `F_A` and `F_B` be the two codecs. On an edge where the anchor is reversed,
the physical output of `F_B(M)`, expressed in the old `A` chart, is obtained
from the physical output of `F_A(M)` by

```
J(x,y)=(1 XOR y, 1 XOR x).
```

On other edges and root bits it is unchanged. `J²=id`. Equivalently
`F_B F_A^(-1)` is this edgewise complement-plus-swap gauge.

Do not confuse this with merely re-expressing one fixed tournament pair in the
new chart: that chart change is simultaneous complement, with no swap.

## 4. Why the anchor is exactly a symmetric-tie sidecar

There is a canonical, relabelling-compatible decomposition of every binary
relation into:

- diagonal loops;
- absent pairs `00`;
- bidirected pairs `11`;
- singly directed pairs `10` or `01`, retaining their intrinsic direction.

Equivalently, the matrix determines an oriented graph on the asymmetric
support plus disjoint undirected pair types and loops. It does **not**
canonically orient its tied pairs. This is safer over `F₂` than attempting a
half-sum/half-difference symmetric/skew decomposition (division by two is
unavailable).

The gauge `J` fixes `10` and `01` and exchanges `00` with `11`. Therefore the
codec's physical output on the asymmetric support is independent of the
anchor: `T` follows the source's unique arc and `U` opposes it. Only symmetric
ties require a chosen orientation.

In coordinates `q=x XOR y`, an edge reversal acts as

```
ordinary source: (x,q) -> (x XOR q,q),
natural target:  (x,q) -> (x XOR 1,q).
```

These agree on `q=1` and disagree on `q=0`. The zero matrix is the extreme
hostile: every off-diagonal pair is a tie, yet every output tournament must
orient every pair.

## 5. Names, roots, class representatives, and global accounting

1. For `n≥2`, `S_n` has even order, so no global `S_n`-equivariant map exists.
   This is stronger than saying that the familiar upper/lower-triangle codec
   happens not to be equivariant.
2. A single named vertex **inside the common core** leaves `S_(n−1)`, so by
   itself it does not repair the obstruction for `n≥3`.
3. A distinguished extra root `*` is already part of `Y(V)` and leaves the
   entire core `S_n`; it alone does not repair the obstruction for `n≥2`.
4. If the only sidecar is `k` distinct, individually named core vertices,
   the residual group is `S_(n−k)`. For `n≥1` the minimum is exactly `n−1`
   names. For `n=0`, no names are needed.
5. An already oriented tournament representative `A`, with or without an
   additional named vertex, has odd automorphism group by the involution
   argument. Its stabilizer therefore supports the local codec. Rooting may
   reduce the stabilizer further but is not needed for parity.
6. This does **not** say that separately chosen unlabelled tournament classes
   can be multiplied freely. Fixing `A` replaces the symmetry group by
   `Aut(A)`; varying `A` adds a genuine object. Relative core alignment and
   automorphism actions still have to be accounted for by the paired-class
   generation procedure.

There is a fully natural enlarged correspondence

```
(A,M) <-> (A,F_A(M)),
```

when `A` is carried and relabelled along with everything else. On the left is
a matrix with an anchor; on the right is a tournament pair with an anchor.
Forgetting the anchor is not free: each side has `2^binom(n,2)` possible
labelled anchors over each underlying matrix or underlying pair.

### Odd subgroup does not mean canonical choice or exact full automorphism group

- `C3≤S3` has exactly two invariant anchors, the opposite directed 3-cycles.
  Its normalizer `S3` interchanges them, so one cannot choose an anchor from
  this subgroup in a conjugation-natural way under all ambient relabellings.
- The construction only guarantees `H≤Aut(A)`, never `H=Aut(A)`. For the
  regular `C3×C3` action on nine vertices there are sixteen invariant anchors;
  the checker finds that **all sixteen have full automorphism order 81**, not
  nine. Nevertheless every one gives the required equivariant codec for `H`.

## 6. Relation to Feit–Thompson and isomorphism algorithms

**CITED connection, not a proof dependency.** Tournament automorphism groups
have odd order by the elementary involution argument above. Feit–Thompson then
implies they are solvable. The primary paper states the odd-order theorem in
its opening theorem:
[Feit and Thompson, *Solvability of groups of odd order*, Chapter I (1963)](https://msp.org/pjm/1963/13-3/pjm-v13-n3-p01-p.pdf).

Solvable automorphism groups are used in the tournament isomorphism and
canonization literature; see [Arvind, Das and Mukhopadhyay, *Isomorphism and
canonization of tournaments and hypertournaments* (2010)](https://www.sciencedirect.com/science/article/pii/S0022000009000774).
The useful bridge here is a structural one: an oriented anchor moves a local
problem into an odd-order stabilizer regime. It is not an automatic complexity
bound for the full, unanchored enumeration problem.

The distinction between containing a prescribed odd group and having exactly
that full automorphism group is substantive. The recent primary paper
[Du, *A note on tournament m-semiregular representations of finite groups*
(2025)](https://www.sciencedirect.com/science/article/pii/S0097316524000918)
records the regular-representation exceptions `C3²,C3³`. Our `C3²` computation
is independent finite-exact evidence for the relevant nine-vertex instance.

**Boundary:** this note rules out the stated same-core, faithful equivariant
codec (indeed any equivariant map from the whole matrix space) when an even
symmetry remains. It does **not** prove that no arbitrary many-one or Turing
reduction from graph isomorphism to a tournament problem can exist. Nor does
it rule out gadget encodings or external choices unless they are required to
retain the forbidden equivariance/faithful automorphism structure.

## 7. Reproduction and exact evidence

Run from this directory:

```
python3 check_anchors.py > check_anchors.out.json
python3 -O check_anchors.py > check_anchors.optimized.out.json
cmp check_anchors.out.json check_anchors.optimized.out.json
```

The checker uses only the Python standard library. All validations use
explicit `require` calls rather than removable Python `assert` statements;
normal and optimized (`-O`) runs produced byte-identical PASS output. It independently represents
tournaments by adjacency rows and compares direct adjacency relabelling with
the codec. Its universe has no graph or matrix restrictions.

- Enumerates all subgroups of `S_n` for `n=0,1,2,3,4` (counts `1,1,2,6,30`).
- For each subgroup, enumerates all invariant anchors two independent ways:
  directed-orbit propagation, and brute-force tournament enumeration. A third
  generator-only parity-constraint path agrees in every case. Its compact
  positive certificates and negative reversal words are independently checked.
- Verifies the zero/even and `2^e_H`/odd anchor counts.
- Exhausts matrices, anchors, and group elements through `n=3`; checks the
  affine basis at `n=4` (equality of affine maps is determined by zero and a
  coordinate basis).
- Checks the reversal cocycle for every anchor and pair of permutations through
  `n=3`; verifies the physical anchor-change gauge independently.
- Uses all sixteen regular `C3×C3` anchors at `n=9`, checks the codec on the
  affine basis, and counts their full automorphisms by enumerating every
  zero-fixing permutation. The verified regular translation subgroup is
  transitive, so orbit–stabilizer recovers the full order.
- Includes a nonfaithful even-group/trivial-image control and the `A4`
  odd-order-generators/even-generated-group control.

| Acting subgroup | Anchor count | Fixed matrices | Fixed natural pairs |
|---|---:|---:|---:|
| `C2` on 2 vertices | 0 | 4 | 0 |
| `C3` on 3 vertices | 2 | 8 | 8 |
| `A4` on 4 vertices | 0 | 4 | 0 |
| regular `C3×C3` on 9 vertices | 16 | 512 | 512 |

The output is an audit of the formulas and implementation, not a replacement
for the all-order proof. No repository or global ledger writes are performed
by the checker.
