# THM-2825 Schur/scale bridge: norm-one tripotent, exact-offset no-go,
# and the two-layer linking algebra

**Status:** `FINITE-EXACT SCRATCH CANDIDATE`, 2026-07-28.  This note makes no
canonical theorem claim.  At the time of the audit, the canonical THM-2825
file was still a `RESERVED / UNPROVED EMPTY STUB`; the computations below
audit the current primary proof companion independently of that status.

## Inheritance pass

- **Closest proved/imported mechanism:** Beke--Goh--Hatami--Jaffe--Naylor,
  arXiv:2607.14316v2, Theorem 1.1, characterizes Boolean Schur multipliers of
  bounded norm as finite signed sums of blocky (contractive idempotent)
  multipliers.  The paper recalls that a nonzero blocky Boolean matrix has
  Schur norm exactly one.
- **Corrected near miss:** the current THM-2825 candidate records the semantic
  sign as the difference of two norm-one Boolean idempotents.  That proves
  norm at most two, but misses a stronger factorization.
- **Canonical hostile:** forgetting half-step scale turns each cell into a
  complete bipartite relation with 189--526 common choices per right piece.
  The relation contains no unique unscaled mate.
- **Least-used sidecar:** ordinary matrix composition of the rectangular
  collar matchings.  Schur composition and ordinary composition behave
  oppositely here; the first destroys offset addition while the second
  exposes a two-state linking algebra.

The live concept board was:

1. blocky Schur idempotents;
2. the semantic `C_2` character;
3. exact half-step offset matchings;
4. the metric collars `J_1,J_2`;
5. ordinary linking/TRO composition; and
6. physical source-carrier mismatch.

## 1. The signed semantic symbol is itself norm one

Let the row and column sets be the labelled right and common pieces,
respectively.  Write

```text
E = union_c R_c x M_c
```

for the 193-cell complete-bipartite relation.  Encode the two semantic values
`(C,C)` and `(0,0)` by signs `epsilon=+1,-1`, and put

```text
S(r,m) = 1_E(r,m) epsilon(r) epsilon(m).
```

Our convention is

```text
||A||_m = sup_(0<||X||<infinity) ||A o X||/||X||,
```

for matrices `X:l2(M)->l2(R)`, exactly the convention in
arXiv:2607.14316v2.  Rows and columns are **labelled by their cell**, so each
belongs to one and only one of the 193 blocks; no physical-interval collision
across labels is quotiented out.

The THM-2825 parity law says this is also `(-1)^delta(r,m)` on `E`.

### Proposition 1 (norm-one contractive tripotent)

The Schur multiplier norm and the `gamma_2` factorization norm of `S` are
both exactly one.  Moreover,

```text
S o S = E,                  S o S o S = S,
P_equal    = (E+S)/2,
P_opposite = (E-S)/2.
```

Here `o` is entrywise/Schur product, and the last two symbols are orthogonal
contractive Boolean idempotents.

### Proof

Let `H=l2({nonempty cells})`.  For `r in R_c` and `m in M_c`, take

```text
u_r = epsilon(r) e_c,       v_m = epsilon(m) e_c.
```

Then `S(r,m)=<u_r,v_m>`, every row and column factor has norm one, and hence
`gamma_2(S)<=1`.  Since `S` has a unit entry, the max-entry lower bound gives
`gamma_2(S)>=1`.  Equivalently, if `D_R,D_M` are the two diagonal sign
unitaries, then

```text
M_S(X) = D_R M_E(X) D_M,
```

so `S` is a unitary row/column gauge of the blocky contraction `E`.
The polynomial identities are entrywise.  This proof is elementary and does
not require the decomposition theorem.

**Strengthening over the current two-term wording:** signs cost *zero*, not a
triangle-inequality factor of two, because they separate into one row phase
and one column phase.

**Type guard:** `S` is generally `{-1,0,+1}`-valued and tripotent, not a
Boolean idempotent.  The external decomposition theorem applies directly to
`E`, `P_equal`, and `P_opposite`; the norm-one statement for `S` comes from
the explicit phase gauge above, not by misapplying that theorem.

## 2. Exact offset resolution and the sharp scale boundary

For every centered half-step offset `k`, let

```text
J_k(r,m)=1  iff  left(m)-left(r)=k h.
```

Translation is injective in each coordinate, so every nonempty `J_k` is a
partial matching and therefore a contractive idempotent Schur multiplier.
The exact audit proves

```text
E = sum_k J_k,                S = sum_k (-1)^k J_k,
J_k o J_l = 0  for k != l.
```

There are exactly `1043` nonempty offsets, ranging from `-560` to `582`.
The edge totals are `97661` on even offsets and `97856` on odd offsets.
Every nonzero `J_k` has norm one, regardless of its distance.

The decisive new census is

```text
|J_1|=|J_2|=587,
|J_k|<587 for every k not in {1,2}.
```

More precisely,

```text
intersection_(r in R) {k : r+k h lies in M_(cell(r))} = {1,2},
max_(k notin {1,2}) |J_k| = 527.
```

This set intersection is how the script checks saturation; it is not inferred
from the nearest-neighbour result.  The complete sorted `(k,|J_k|)` census is
hashed in the output.  All positive offsets `1,...,582` occur.  There are no
offsets from `-49` through `0` (the nearest negative offset is `-50`), while
offsets `1,2` have size 587 and offsets `3,...,14` have size 527.  The
direction is therefore a genuine one-sided boundary feature.

Thus:

- `J_1` is the unique **odd** exact-offset matching with full right domain;
- `J_2` is the unique **even** exact-offset matching with full right domain;
- the geometric nearest/opposite and nearest/same collars can be
  recharacterized as `parity + domain saturation`, once the exact-offset
  grading has been supplied;
- the common bank has a uniform initial collar of depth exactly two.  Sixty
  rows already lose offset `3`; no third exact-offset matching is uniform.

This is a positive metric-free characterization only **after** the offset
grading is present.  Schur norm alone cannot supply that grading: all 1043
nonempty offset matchings have the same norm one.

## 3. Why Schur composition cannot recover the collar

Even enrich the coarse relation by:

- all 193 cell projections;
- separate row and column semantic-colour projections; and
- arbitrary sums, complements, and Schur products of those masks.

The nonzero atoms of this Boolean/Schur algebra **inside `E`** are the
semantic rectangles

```text
R_(c,a) x M_(c,b),       a,b in {+1,-1}.
```

There are exactly `414` such atoms, `207` of each relative sign.  Every
column side has size at least `94`.  The matching `J_1` cuts each of the 207
opposition atoms in a graph of size `|R_(c,a)|`, strictly smaller than the
rectangle size `|R_(c,a)| |M_(c,-a)|`; `J_2` does the same to all 207 equality
atoms.  Hence neither collar lies in this algebra.

This gives two equivalent no-go mechanisms.

1. **Atom obstruction:** an element generated from cell and colour masks is
   constant on every semantic rectangle, whereas a collar selects one column
   per row.
2. **Automorphism obstruction:** independent permutations inside each row-
   and column-colour class preserve all generators but move the nearest
   graph.  Any invariant Boolean selector is a union of whole rectangles,
   never `J_1` or `J_2`.

Composition does not help because composition of Schur multipliers is
entrywise multiplication.  Distinct exact offsets annihilate:

```text
M_(J_k) M_(J_l) = 0,       k != l.
```

The additive law `k+l` of translations is lost.  In particular the algebra
generated only by `E` and `S` collapses all offsets to their parity classes.

## 4. Ordinary composition exposes the missing two-state carrier

Now regard `J_i` as a rectangular matrix

```text
J_i : l2(M) -> l2(R),       i=1,2.
```

Both matchings have every right row and injective common columns.  The audit
finds that their two 587-column images are disjoint.  Therefore

```text
J_i J_j^* = delta_(i,j) I_R.
```

Put

```text
V_ij = J_i^* J_j.
```

Then

```text
V_ij V_kl = delta_(j,k) V_il.
```

So the two uniform collar layers generate an exact

```text
M_2 tensor I_587
```

linking algebra on their 1174-dimensional common-column subspace.  In
particular,

```text
X=V_12+V_21,       Z=V_22-V_11,       Q=V_11+V_22

X^2=Z^2=Q,         XZ=-ZX.
```

`X` swaps the `+h` and `+2h` common collars, hence is a partial translation
by one half-step; `Z` is the relative semantic phase.  This is the
noncommutative two-state/Pauli refinement of the commutative parity symbol.

At the set level, the same operator is the fixed-point-free involution

```text
x(j_1(r))=j_2(r),       x(j_2(r))=j_1(r).
```

It partitions the 1174 common pieces into 587 adjacent pairs and reverses
the delayed semantic value on every pair.  Thus the *paired* collars do
produce an abstract free `C_2`-torsor over the labelled right set, even
though `j_1` alone does not act on its right-domain.

The direction matters:

```text
rank(J_1^*J_2)=587,        rank(J_1J_2^*)=0.
```

Thus the useful shift lives on the **common-column** side.  It does not turn
the right cofiber into a closed two-state fibre.  The torsor also depends on
the already metric-selected correspondence through the right labels, so it
is not the missing intrinsic physical bibundle.

## 5. Consequence for HYP-9046

### Positive import

Whenever a signed wall symbol has the form

```text
signed_support(x,y)=alpha(x) 1_A(x,y) beta(y),
|alpha|=|beta|=1,
```

its Schur norm is exactly the norm of the unsigned support.  Factorable signs
should therefore be gauged away before invoking the Beke--Goh--Hatami--
Jaffe--Naylor decomposition.  The THM-2825 parity sign is a complete finite
positive control: it has norm one directly.

This separates two obligations that HYP-9046 currently places too close
together:

1. bound the norm of the Boolean wall incidence pattern; and
2. preserve the analytic signs in the pairing.

Row/column-separable signs do not affect obligation 1 at all.

### Exact obstruction

The external theorem is an existence theorem for a signed sum of blocky
Boolean masks.  It does not make the decomposition canonical, retain an
offset label, or turn Schur composition into translation composition.
THM-2825 shows the loss in the cleanest possible form:

- the coarse signed multiplier is already norm one;
- every exact scale slice is also norm one;
- only an external offset/order sidecar says which slice is nearest;
- only ordinary linking composition recovers the half-step shift; and
- that shift appears only after `J_1,J_2` were selected.

Therefore contractivity cannot by itself recover metric or chronology.
For HYP-9046 the true load-bearing test remains the uniform Schur norm of the
actual large-diameter relation/window incidence matrix.  This 193-cell
positive control neither proves nor meaningfully estimates that norm.

### Physical-support boundary

The current THM-2825 physical audit gives

```text
R source carrier  = empty,
J_1(R),J_2(R) source carrier = delta_0
```

on all 587 collar pairs.  Restricting or reweighting their pair mask cannot
alter either vertex sidecar.  The `M_2` linking algebra therefore repairs the
**common-side semantic/metric grammar**, not the absent physical source
support.  A lawful LRC application still needs a co-supported physical
translation or an explicitly broader two-object carrier.

## 6. Suggested compact theorem addition

If THM-2825 is promoted, a concise safe strengthening is:

> The signed semantic symbol is a norm-one Schur tripotent, not merely a
> two-term signed sum of contractive idempotents.  Its exact-offset
> resolution consists of 1043 pairwise Schur-orthogonal contractive matching
> idempotents.  Exactly `J_1` and `J_2` have the full 587-row domain, uniquely
> selected by odd/even parity plus saturation.  Their column images are
> disjoint, and `V_ij=J_i^*J_j` form matrix units for
> `M_2 tensor I_587`; nevertheless the cell/semantic Schur algebra cannot
> generate either collar, and the linking shift remains common-side rather
> than physically co-supported.

## 7. Reproduction

```bash
python3 .scratch/lrc_schur_bridge_20260728/schur_scale_algebra_audit.py \
  > /tmp/lrc_schur_bridge.out
python3 -O .scratch/lrc_schur_bridge_20260728/schur_scale_algebra_audit.py \
  > /tmp/lrc_schur_bridge.opt.out
cmp /tmp/lrc_schur_bridge.out /tmp/lrc_schur_bridge.opt.out
```

Evidence hashes (LF-normalized bytes):

```text
script: 6d9a96d6c2e1c0e631bf0160a4113705a973baefe640683827f33c11ced67f26
output: 6f90f04c1efa0a6be5f19bdb5e9a353dcda5e7ebb70cbabb8a34f9d62fc5b1d4
```

The script pins the current primary THM-2825 companion hash and contains no
executable Python `assert`.
