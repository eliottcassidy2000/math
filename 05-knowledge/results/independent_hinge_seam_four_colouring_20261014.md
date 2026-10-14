# Independent hinges decouple: a global four-colouring and the next affine operation

**Date:** 2026-10-14. **Status: PROVED; FINITE-EXACT certificates. Novelty unconfirmed.**

**Correction / supersession:** Hypothesis H-G in [Euclid rigid hinges and clutch growth](euclid_rigid_hinges_and_clutch_growth_20261014.md), published at `be3fdf3cb11c895779bee21288b5584c2ee717f9`, is **REFUTED**. Every graph in that proposed family has an explicit four-colouring. The finite tests led to a theorem for the entire infinite union, not merely two failed search stages. All other claims of the earlier note, including its independent H-J, retain their stated scopes.

## 1. Useful result first

Independent quadratic hinge rotations about one common centre cannot make this triangular-lattice construction exceed four colours. Algebraic independence forces every cross-layer unit edge to change just one rotation index, and an elementary denominator argument places both ends in the same one of the triangular lattice's three colour classes. Split that class by a parity rule and four colours suffice globally.

The obstruction survives overlapping seam norms. Multiplicative dependence among the rotation denominators alone does not fix the construction: the parity rule simply sums over all eligible rotation indices. A productive next operation must break the common-centre decomposition or introduce genuine radical dependencies that permit mixed-index seams.

## 2. Exact theorem

Let `d_1,...,d_r` be integers at least 2, each divisible by 3. Set

`D_i=4d_i-1`, `cos(theta_i)=1-1/(2d_i)`, `sin(theta_i)=sqrt(D_i)/(2d_i)`.

Assume the squareclasses `3,D_1,...,D_r` are linearly independent in `Q*/Q*2`. Let `R_i` be the positive-angle rotation and let

`Lambda=Z e_1+Z e_2`, `e_1=(1,0)`, `e_2=(1/2,sqrt(3)/2)`.

Consider the complete unit-distance graph on

`X = union_(a in Z^r) R_1^(a_1)...R_r^(a_r) Lambda`.

### Theorem

1. Every nonzero point of `X` has a unique label `(a,u)`, with `a in Z^r`, `0!=u in Lambda`. The origin is shared by all layers.
2. If two nonzero points in different layers are unit-distance, their layer indices differ in exactly one coordinate, say `i`, by a nonzero integer `k`. Their lattice vectors are equal or opposite, and their common squared norm is exactly `d_i^|k|`.
3. The graph on `X` is four-colourable by the explicit formula below.

No assertion that every candidate norm actually produces an edge is needed. The norm condition is necessary, not sufficient.

### 2.1. Radical separation of layers

Write `u=(p+q/2,sqrt(3)q/2)` and similarly `v`. Their dot product is rational, and their determinant is `sqrt(3)` times a rational number.

For each nonzero integer `k`, both `cos(k theta_i)` and the rational coefficient of `sqrt(D_i)` in `sin(k theta_i)` are nonzero. Otherwise the unit complex number `exp(i theta_i)` would be a root of unity. Its rational trace `2-1/d_i` would then be an algebraic integer, hence an integer, which it is not.

For a relative rotation involving at least two nonzero indices `i,j`, the squared-distance expansion contains:

- the basis term `sqrt(D_i D_j)` with nonzero multiplier times `u dot v`;
- the basis term `sqrt(3D_i)` with nonzero multiplier times `det(u,v)/sqrt(3)`.

All other rotation factors contributing to these coefficients are nonzero cosines. Independence of the squareclasses makes these distinct irrational basis coordinates; neither can be cancelled by the rational target distance or by another subset of radicals. Thus a rational squared distance forces both dot product and determinant to vanish. For nonzero real planar vectors that is impossible.

With just one nonzero relative rotation index, the irrational coefficient instead forces `det(u,v)=0`.

The same argument applies to a putative coincidence, whose squared distance is zero. In the single-index case, coincidence between collinear nonzero vectors would require the rotation to be `+1` or `-1`, already excluded. Hence nonzero labels are unique. This resolves representation and negative-vector ambiguities before defining the colouring.

### 2.2. Collinearity and the exact denominator obstruction

For a single relative rotation, write `u=s w`, `v=t w`, where `w` is a primitive triangular-lattice vector and `s,t` are nonzero signed integers. Put `Q=|w|^2`, a positive integer. If `|s|!=|t|`, the squared distance is

`Q (|s|-|t|)^2 + 2Q |st| (1-sign(st) cos(k theta_i)) > 1`.

Therefore a unit edge requires `s=+/-t`, and the common norm `q=|u|^2` satisfies

`2q (1-sigma cos(k theta_i))=1`, where `sigma=+1` or `-1`.

For `h=|k|`, define the integer sequence

`C_0=2`, `C_1=2d_i-1`,

`C_h=(2d_i-1)C_(h-1)-d_i^2 C_(h-2)`.

The trace identity gives `C_h=2d_i^h cos(h theta_i)`. Reduction modulo each prime dividing `d_i` shows `gcd(C_h,d_i)=1` for `h>=1`. The unit equation becomes

`q = d_i^h / (2d_i^h-sigma C_h)`.

The positive integer denominator is coprime to `d_i`. Since `q` is an integer, that denominator must be 1. Thus `q=d_i^h`, proving the seam classification.

This trace recurrence is the actual prime-power connection: denominators are restricted to powers of one base, whereas the numerator is coprime to that base. It is not a heuristic appeal to prime density.

### 2.3. Explicit global colouring

Give the triangular vector `u=p e_1+q e_2` its usual colour

`c(u)=p-q mod 3`.

If `c(u)=1` or `2`, keep that colour. For `c(u)=0`, put `N=|u|^2` and let

`I(N)={i: N=d_i^(h_i) for some positive integer h_i}`.

For each `i`, the exponent `h_i` is unique because `d_i>=2`. Colour the point with label `(a,u)` by

`0` if `sum_(i in I(N)) floor(a_i/h_i)` is even, and `3` otherwise.

Colour the origin 0. If `I(N)` is empty the sum is zero.

Within any one layer, every unit edge changes the triangular colour. Replacing its zero-colour class by colours 0 and 3 cannot create a conflict. For a cross-layer unit edge, the seam theorem gives `N=d_i^h`, divisible by 3. The identity

`p^2+pq+q^2 = (p-q)^2 mod 3`

then puts both endpoints in the zero-colour class. Their indices differ only in coordinate `i`, by exactly `+h` or `-h`; its floor term changes by one and every other term stays fixed. Hence their colours differ. Changing the lattice vector to its negative preserves this zero class and its norm, so it causes no ambiguity. Edges incident with the common origin have norm one and triangular colour 1 or 2.

This is a proper four-colouring of the entire infinite graph. It remains valid when different bases have common positive powers, because every eligible index is included in the sum. ∎

## 3. The original mixed-parity hypothesis is exactly settled

Take `(d_1,d_2)=(3,12)`, giving `(D_1,D_2)=(11,47)`. The squareclasses `3,11,47` are independent. The original `G_m` is the induced subgraph obtained by `|a|,|b|,|p|,|q|<=m`. It contains the first rigid hinge for `m>=1`, so

`chi(G_m)=4 for every m>=1`,

and the whole infinite union also has chromatic number four.

The first exact searches found:

- `m=1`: 73 distinct vertices, 156 unit edges, four-colourable;
- `m=2`: 601 distinct vertices, 1,560 unit edges, four-colourable.

Glucose initially produced colourings during discovery. The retained reproduction script has **no SAT dependency**: it enumerates all pairs in the exact field basis, constructs the theorem's formula colouring directly, verifies every edge, checks uniqueness of nonzero labels, and checks every nonzero cross-layer edge against the seam theorem. The stored certificates contain the labels, complete edge lists and colours.

This is stronger than the earlier odd-tower local four-colouring. The present result also covers mixed local squareclasses and even hinges; it identifies a geometric/algebraic separation specific to common-centre unions of rotated lattices. It does not four-colour the entire number-field plane `Q(sqrt(3),sqrt(11),sqrt(47))^2`.

## 4. What failed, and what remains worth trying

The first failed implication in H-G was:

`loss of the old 2-adic colouring certificate => mixed rotations couple their colour constraints`.

The global seam theorem supplies a different colouring certificate. Fresh global fields and a changed local squareclass do not overcome it.

The loss ledger is explicit:

- **Source:** actual exact-coordinate unit-distance configurations.
- **Map:** express each point as a lattice vector plus its rotation-layer label.
- **Preserved:** every edge; base colour; squared norm; the unique active rotation index of each seam.
- **Restored sidecar:** signed lattice vector, all layer indices, and the denominator-power exponent.
- **Decisive obstruction:** no seam can involve two active independent radicals at once.

Norm overlap by itself is not a repair. Even if `d_i^h=d_j^k`, the floor-parity sum handles both seam directions. Increasing the number of independent common-centre hinges is likewise blocked by the theorem.

### Concrete next operation: affine attachment with a shared unit edge

Use translations to change the field-coefficient equations rather than merely adding a new independent angle. Let `t=e_1` and define

`Y_m=G_m union (t+G_m)`.

Both copies contain a unit edge from 0 to `e_1`; they share at least that edge, and their base triangular patches overlap substantially. Retain all point identifications and all unit edges. This is an explicitly realizable graph, with at most twice the old vertex count, over the same exact field.

For a general translated cross-distance,

`|u-(t+Rv)|^2=|u-t|^2+|v|^2-2 (u-t) dot (Rv)`.

The previous proof no longer compares two lattice vectors about a common origin with one fixed squared norm. This allows different interfaces to impose compatibility on the same colour classes. It is not a guarantee of higher chromatic number.

**Minimal next test:** compute `Y_1` and `Y_2`, retaining the shared edge and the complete 4-colour boundary relation on the overlap. The important measurement is whether every colouring of one copy has an incompatible overlap restriction in the other, not just whether each copy is four-chromatic. If both remain four-colourable, return the actual overlap relation or a coherent extension rule. These affine tests are proposed here, not claimed to have been run.

A second lawful route is a dependent-radical design in which the basis separation genuinely fails. Repeated squareclasses alone may still preserve a colouring; they require a fresh seam analysis. For example `12*378^2-1=47*191^2` reuses the radical of the `n=2` hinge. This arithmetic identity is an available design parameter, not evidence for a five-chromatic graph.

## 5. Reproduction and status

Run from repository root:

`python3 04-computation/experiments/independent_hinge_seam_four_colouring_20261014.py`

`python3 -O 04-computation/experiments/independent_hinge_seam_four_colouring_20261014.py`

Both outputs match `05-knowledge/results/independent_hinge_seam_four_colouring_20261014.json` byte-for-byte. The script uses only Python's standard library. Every pair is tested in the basis

`1,sqrt(3),sqrt(11),sqrt(33),sqrt(47),sqrt(141),sqrt(517),sqrt(1551)`.

It first uses the exact rational component as a rejection test, then all other field coefficients; no approximate-distance cutoff occurs. The two tests comprise 2,628 and 180,300 unordered pairs. The finite tests validate the retained certificates; the proof above establishes all stages and the generalized theorem.

No claim is made about a new lower bound for the plane, all graphs over these fields, or literature priority. H-G is refuted for a structural reason; the successor target changes the geometric operation instead of repeating the same search at larger size.
