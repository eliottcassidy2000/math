# A three-coordinate transfer with distinct entry and exit

Session 197. This is an all-order physical extension of the one-seam rooted interface, together with a bounded-span example having a genuine constant spectral component. It is not an arbitrary-target construction.

## 1. Definitions

Let `(A,a)` be a nonempty rooted tournament. Write `h=H(A)`, `e=e_a(A)` for Hamiltonian paths ending at `a`, and `R=R_a(A)` for ordered spanning path covers `(P,Q)` where `P` may be empty and `Q` is nonempty and starts at `a`. Each nonempty slot is a directed path. The empty path has one choice in these cover sums.

Let `G` be another nonempty tournament, with designated entry `b` and exit `c`, which may coincide. Form `T` by directing all cross edges from `A` to `G`, except for `b->a`. The returned root is `c`, not necessarily `b`.

Define the following constants of the gadget:

- `g=H(G)` and `t=s_b(G)`, its number of Hamiltonian paths starting at `b`;
- `L=L_b(G)`, the number of ordered covers `(U,V)` where `U` is nonempty and ends at `b`, and `V` may be empty;
- `E=e_c(G)`;
- `F=F_(b,c)(G)`, the number of Hamiltonian paths starting at `b` and ending at `c`;
- `D=D_(b,c)(G)`, the number of ordered covers `(U,V)` with both paths nonempty, `U` ending at `b`, and `V` ending at `c`.

Thus `D=0` if `b=c`. Also `F=0` if `b=c` and `|G|>1`, while `F=1` for a singleton gadget.

Make an augmented gadget `G^+` by adding `z` with `z->c` and `x->z` for every `x` in `G-{c}`. Keep the entry at the old vertex `b`, and put

`g^+=H(G^+)`, `t^+=s_b(G^+)`, `L^+=L_b(G^+)`.

## 2. The exact returning matrix

**Theorem.** The output `(h',e',R')=(H(T),e_c(T),R_c(T))` depends only on `(h,e,R)` and is

```
h' = g h - t e + L R,
e' = E h - F e + D R,
R' = (g^+-g+E) h + (t-t^+-F) e + (L^+-L+D) R.
```

The entries are integers computed from actual gadget path counts. This statement holds for all nonempty sources and gadgets, including coincident terminals and the singleton boundary cases.

### Proof of the first row

A Hamiltonian path avoiding the only backward cross edge consists of an `A` path followed by a `G` path. Of the `g h` pairs, precisely `t e` try to use the absent edge `a->b`.

If the path uses `b->a`, it has the unique run pattern

`A_before, G_before, A_after, G_after`,

with `A_before` and `G_after` allowed to be empty, `G_before` ending at `b`, and `A_after` starting at `a`. The two `A` pieces form a cover counted by `R`, and the two `G` pieces one counted by `L`. Other joins are automatically legal: `a` is already in the second `A` piece and `b` in the first `G` piece. This gives `L R`.

### Proof of the second row

Require the entire path to end at `c`. In the avoiding case there are `E h` joins, and exactly `F e` are forbidden. In the using case the last `G` run must be nonempty and end at `c`, so its two `G` pieces are counted by `D`. This gives `D R`. If `b=c`, the using case is impossible because the path already visited `c` before crossing to `A`; correspondingly `D=0`.

### Proof of the third row by regrouping

Attach one more vertex `z` to `T`, directing `z->c` and all other edges toward `z`; call the result `T^+`. The first-row formula applied with singleton gadget and old root `c` gives

`H(T^+)=H(T)-e_c(T)+R_c(T)`.

Regroup the same physical graph as source `A` and gadget `G^+`. Its only reverse cross edge is still `b->a`, since `z` loses to every vertex of `A`. Thus

`H(T^+)=g^+ h - t^+ e + L^+ R`.

Subtract the first-row expression for `H(T)` and add the second-row expression for `e_c(T)`. This proves the third row without any unproved assumption that a fitted three-state model is sufficient.

## 3. A literal bounded-span compiler

If the source has an ordering of backward span at most `w` with root `a` last, and the gadget has such an ordering with exit `c` last and entry `b` in one of its first `w` positions, concatenate these orders. All old cross edges point forward except `b->a`, whose distance is the position of `b` in the gadget, hence at most `w`.

Consequently an arbitrary sequence of gadgets meeting these conditions gives actual span-`w` tournaments and is tracked exactly by the three-coordinate matrix. The gadgets may differ at every step. This is a sufficient construction class, not a classification of all span-`w` tournaments, and it does not imply that its scalar output is complete.

## 4. A second-order affine clock at width three

Take the ordered four-vertex gadget with adjacency bit rows

`G=[6,4,8,3]`, entry `b=2`, exit `c=3` (zero-based indices).

Its edges are `0->1, 0->2, 1->2, 2->3, 3->0, 3->1`. It has backward span three. Its constants are

`(g,t,L,E,F,D)=(5,1,7,1,0,1)`,

`(g^+,t^+,L^+)=(11,2,20)`.

The theorem therefore gives

```
M = [[5,-1,7],
     [1, 0,1],
     [7,-1,14]].
```

Start with `T_1=G` rooted at `3`, whose profile is `(5,1,7)`, and repeatedly attach the same gadget. The resulting ordered graph is exactly the width-three mask word `(0013)^k`: bit zero addresses the oldest retained vertex, so the third new vertex points back to the previous exit, and no other reverse cross edge occurs.

The row `ell=(-2,1,1)` obeys `ell M=ell`. Hence

`-2h_k+e_k+R_k=-2` for every `k>=1`.

Let `f=(1,0,0)`. Direct multiplication gives

`f(M^2-18M+5I)=6 ell`.

Applying this identity to the physical state at step `k` proves the exact scalar recurrence

`h_(k+2)=18h_(k+1)-5h_k-12`, for all `k>=1`,

with `h_1=5`, `h_2=73`. Its first values are

`5,73,1277,22609,400565,7097113,125745197`.

Equivalently, `h_k-1` satisfies the homogeneous recurrence with characteristic polynomial `x^2-18x+5`. The full matrix characteristic polynomial is

`(x-1)(x^2-18x+5)`.

The constant component is physically observable: the second-order recurrence has the nonzero constant term `-12`. In particular, the tempting blanket assertion that every nontransitive periodic ribbon has no observable eigenvalue one is false. This example does not contradict the separately proved exclusions of fixed first-order affine clocks with multipliers two, three, and four. It is a second-order law and supplies no inverse compiler for arbitrary target counts.

## 5. Checks and source scope

Run `python session197/ribbon-terminals/verify.py`. The standalone standard-library checker computes all constants from literal adjacency graphs, checks every source of orders one through three against every gadget of orders one through four with every choice of entry and exit, and verifies all three output coordinates by independent subset dynamic programming. It separately checks the repeated literal word graphs, the matrix identities, and the ordered-span compiler's premises.

The inherited first row is also recorded in the specifically named `../../session197-lowcycle/ROOTED.md`, section 1, which identifies its S195 rooted-interface source. The first-row proof is reproduced here. The distinct-terminal ending row and augmentation/regrouping argument are independently developed consequences; no priority or novelty claim is made.
