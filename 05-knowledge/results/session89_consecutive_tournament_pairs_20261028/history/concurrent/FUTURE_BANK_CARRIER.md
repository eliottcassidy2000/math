# Exact smallest future-decision carrier for the rectangle code

**Status:** elementary all-size theorem, independently checked in finite universes. This refines the named-bank obstruction in [BANK_MATRIX_BOUNDARY.md](BANK_MATRIX_BOUNDARY.md); it does not make the square-sigma predicate an invariant of natural pair isomorphism classes.

All arithmetic below is over `F_2`. Let `V_n=F_2^(n×n)` and

```
C_n = { M : M_ij = a_i + b_j }.
```

The decision is `M∈C_n`. For the explicitly ordered nine-prime bank at `n=3`, this is exactly whether the selected-prime product has square sigma. Allowed operations are **each named single-cell toggle**, and **natural simultaneous paired-tournament relabelings** `alpha_g`, `g∈S_n`. No arbitrary row/column action, unknown prime dictionary, or forgotten label convention is substituted.

## The action and the operational quotient

Use `x_ij=M_ij`, `y_ij=M_ji` for `i<j`, and `d_i=M_ii`. Write
`alpha_g(M)=P_g(M)+k_g`. With `u=g^(-1)(i)`, `v=g^(-1)(j)`, `l=min(u,v)`, `h=max(u,v)`, its coordinates for `i<j` are

```
(alpha_g M)_ii = M_uu;
(alpha_g M)_ij = M_lh + [u>v];
(alpha_g M)_ji = M_hl + [u>v].
```

Thus `P_g` is a coordinate permutation, but generally not ordinary simultaneous matrix conjugation. Every allowed operation sequence has the form `M -> alpha_g(M)+t`, and every such form is available: named toggles generate every translation `t∈V_n`.

**Theorem (exact operational carrier).** Two inputs have the same eventual decision after **every** finite allowed operation sequence if and only if their difference belongs to

```
D_n = intersection_(g∈S_n) P_g^(-1)(C_n).
```

Consequently the coarsest exact deterministic state space retaining those decisions and transitions is `V_n/D_n`. For `n>=3`,

```
D_n = { M : M_ij = a_i + a_j + c },
dim D_n = n,
number of necessary states = 2^(n²−n).
```

It requires exactly `n²−n` fixed-width binary state bits, without a linearity assumption on the competing exact deterministic carrier. The exceptions are `D_2=C_2`, dimension 3 (two states, one bit), and `D_1=V_1` (one state, zero bits).

### Proof of the operational statement and minimality

For a linear code `C`, its indicator has the same values on every translate of `x` and `y` exactly when `x+y∈C`: if the difference is in `C`, all memberships agree; otherwise translation by `x` sends `x` to `0∈C` and `y` outside `C`. Apply this after every `alpha_g`. Differences remove the affine offsets, giving the stated intersection.

More explicitly, if `delta=x+y` is outside `D_n`, choose `g` with `P_g delta` outside `C_n`. Relabel by `g` and then toggle the cells of `alpha_g(x)`. The resulting decisions are true for `x` and false for `y`. Thus distinct cosets really require distinct deterministic states. Conversely, `D_n⊂C_n` and `P_g D_n=D_n`; each toggle and affine relabeling descends to the quotient and its decision is well-defined. This proves both sufficiency and minimality for the declared operations. Only the linear parts preserve the kernel: `P_g D_n=D_n`; generally `alpha_g(D_n)=k_g+D_n`, not `D_n`. The induced quotient action is affine, not an assertion that the zero coset is fixed.

### Proof of the common kernel

The identity permutation forces `M∈D_n` to lie in `C_n`; write `M_ij=a_i+b_j` and set `d_i=a_i+b_i=M_ii`. Let `s=(i,i+1)` be an adjacent transposition. Its linear action differs from ordinary simultaneous conjugation `Q_s` only at the two cells of this reversed pair:

```
P_s M + Q_s M = (d_i+d_(i+1)) (E_(i,i+1)+E_(i+1,i)).
```

Ordinary conjugation preserves `C_n`. For `n>=3`, `E_ij+E_ji` is not in `C_n`: choose `k` different from `i,j`; the rectangle with rows `i,k` and columns `j,k` has parity one. Therefore `P_s M∈C_n` forces `d_i=d_(i+1)`. All adjacent transpositions force a common value `c`, whence `b_i=a_i+c` and `M_ij=a_i+a_j+c`.

Conversely such matrices are symmetric, so on them every `P_g` equals ordinary conjugation, which preserves the displayed form. They lie in `D_n`. Their parameters `(a,c)` have only the kernel `(a,c)=(constant,0)`, so the dimension is `n`. For `n=2`, `C_2` is the single even-total-parity constraint, preserved by the coordinate permutation `P_g`; for `n=1` there is no constraint.

## Explicit minimal coordinates, not just a state count

For `n>=3` define the zero-first-row representative

```
Q(M)_ij = M_ij + M_0i + M_0j + M_00.
```

Its first row vanishes. The difference `M+Q(M)` lies in `D_n`, and restriction of `D_n` to the first row is a bijection onto `F_2^n`. Hence each coset has exactly one such representative. Store its `n(n−1)` remaining cells. To update, apply the declared toggle or affine action to the representative and normalize again with `Q`; to decide, test its rectangle syndrome. No original first-row information is needed for these declared future decisions.

For a matrix in `C_n`, this representative satisfies

```
Q(M)_ij = (a_i+b_i) + (a_0+b_0),
```

independent of `j`. Thus the current-decision subspace `C_n/D_n` consists precisely of rectangles whose every row is constant, and has dimension `n−1`. Its full syndrome can be stored as the differences between each row's first entry and its other entries.

The checker exposes importable `canonicalize`, `carrier_toggle`, `carrier_relabel`, `carrier_predicate`, and `carrier_syndrome` functions for `n>=3`. It tests those exact functions, not merely a parallel formula.

The complete **current rectangle syndrome** has `(n−1)²` bits. It is already the minimal carrier for all future named toggles *without* relabeling, namely `V_n/C_n`. Adding natural relabelings costs exactly `n−1` more bits when `n>=3`. The isolated present yes/no answer itself is only one bit; it must not be confused with the complete syndrome or an operation-closed carrier.

**Source / target / map / preserved / lost:** source is the labelled matrix with its named arithmetic-bank interpretation; target is `V_n/D_n`; map is the quotient, explicitly `Q`; preserved are all future membership decisions under the stated toggles and natural relabelings, with deterministic updates. Lost are the full matrix and prime subset, exact integer product, and any observable not proved constant on these cosets. This is not a quotient by natural pair isomorphism.

## Exact nine-port consequence and boundary with F64

At `n=3`, the present four-bit syndrome refines to a **six-bit** future carrier, not the full nine-bit matrix. Its 64 states have eight preimages each, and its kernel has weight enumerator

```
1 + 3 z^4 + 3 z^5 + z^9.
```

The six-bit representative is a `2×3` binary rectangle. Regard its three columns as vectors `u_0,u_1,u_2` in a chosen `F_2^2` basis. The current-decision kernel is exactly

```
{(u,u,u): u∈F_2^2},
```

and a four-bit syndrome is `(u_0+u_1, u_0+u_2)`. In particular the nine states with one nonzero column and the other two zero have nine distinct nonzero syndromes. This is the `r=3, k=2` representative `L(x,y,z)=(x+y,x+z)` already given in [Session 87 BLOCK_COMPRESSION.md, Section 2](https://github.com/eliottcassidy2000/math/blob/7a0bab9bfa7bb87c5a10c7ba239592d90cfafd51/05-knowledge/results/session87_four_vertices_enumeration_20261028/BLOCK_COMPRESSION.md). The explicit isomorphism to that three-block/two-bit linear decision model is: source is the chosen `2×3` additive-coordinate model; target is three named `F_2^2` blocks; map reads the columns; preserved are addition, the diagonal decision kernel, and all six state bits. No information is lost in this column regrouping. It is **not** an identification of these nine selected states with the original nine prime-coordinate lifts. Any further interpretation as an `F64` coordinate system requires a declared additive basis and grouping; no multiplication or Frobenius compatibility is asserted.

The earlier named nine-port map uses six-bit columns, in this same cell order,

```
[ 1,  4,  5]
[34,  2, 32]
[16, 24,  8],
```

followed by an invertible basis map into `F64`. Its kernel instead consists of unions of whole rows, with enumerator

```
1 + 3 z^3 + 3 z^6 + z^9.
```

Both quotients have abstract binary dimension six, but the kernels differ even up to a permutation of the nine cells, since cell permutations preserve Hamming weights. Thus this theorem supplies **no prime-index permutation or canonical field identification** with the earlier F64 carrier. In particular the older port map is not automatically a sufficient statistic for these newly declared future contexts. The field's dimension alone is not evidence of equivalence.

## Verification scope

`check_future_bank_carrier.py` is standalone and imports no producer code. It checks transported rectangle-row ranks over **every permutation** for `n=1,...,8`; the ranks are `0,1,6,12,20,30,42,56`. At `n=3` it enumerates all 512 matrices and all `6×512=3072` relabel-then-translate decision contexts per input. The resulting full truth profiles, six-coordinate linear profiles, and affine rectangle profiles induce exactly the same 64 eight-element classes as `D_3`. It also checks explicit normalization, all named toggle and natural-relabel transitions, and both kernel weight enumerators. Normal and optimized Python runs agree at **1,113,285 always-active checks**.

This is an exact minimality result for a specified finite-state decision problem, not a tournament-isomorphism algorithm, an arithmetic classification of arbitrary prime banks, or a claim that equal-dimensional encodings preserve one another's operations.
