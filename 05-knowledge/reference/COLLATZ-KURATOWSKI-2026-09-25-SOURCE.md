> **HISTORICAL INPUT / CURRENT-TRUTH WARNING.** This is the user's pasted
> 2026-09-25 Kuratowski--Tutte note, preserved for provenance. Its status
> claims are source claims, not repository endorsements. The orbit corollary
> in section 3.3, the completeness claims in sections 7.2/7.4, and the broad
> proof-method exclusion in section 7.1 need the repairs in
> [the current reframe and audit](../results/kuratowski_reframe_20260925.md).
> Original attachment: ec4421e9-d293-4596-b9d3-aa8e3562b85b/Pasted text.txt.
# The Kuratowski–Tutte triple, the (3,4,5) root, and where Collatz edges live in the Berggren tree

**Date:** 2026-09-25. **Session:** cowork-20260925-kuratowski (machine `kind-pasteur`).
**Status:** PROVED elementary statements (sections 1–5 and 7.1, each with
its own proof and exhaustive finite control); FINITE-EXACT controls;
CITED classical graph theory (section 6); ANALOGY / NUMEROLOGY typed
explicitly where no map exists; **Collatz remains OPEN**. No novelty claim
beyond the repository; no Lean claim. Two independent adversarial audits
(recorded in section 9) found and repaired eight wording/sign/quantifier
defects in the first version; the corrected statements are the ones below. The owner's brief asked for an
isomorphism between the triple {Petersen, K_5, K_{3,3}} and the root (3,4,5)
of the primitive Pythagorean triples, with "complete graphs ↔ difference-1
leg, bipartite ↔ doubled leg, snarks ↔ the family approaching the isosceles
right triangle", and for a creative pursuit of *universal coverage*.

## Inheritance and concept board

Closest proved mechanisms: the odd-root chart and the Berggren branches on
root pairs `B1(s,t)=(s+2t,t)`, `B2=(2s+t,s)`, `B3=(2s-t,s)` with the unique
parent map on three cones ([THM-3756](../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md);
[ternary_berggren §1](ternary_berggren_20260925.md)); the marked triple
`Psi(n)=(2n,n^2-1,n^2+1)/d(n)` and the spine `Phi(n)` with the middle-child
reading of the plus shortcut ([zenodo_triples §2, §4](zenodo_triples_20260925.md));
the sibling ray `(x_j,y)=B1^{H(j)}(x_0,y)` ([edge transport §6](collatz_mod6_20260921_berggren_edge_transport.md),
[ternary_berggren §3](ternary_berggren_20260925.md)); the consecutive-edge
incomparability theorem ([ternary_berggren §5](ternary_berggren_20260925.md)).
Canonical hostile: the legal edge `7->11`, whose Berggren parent `(7,3)` is
not a plus edge. Corrected near miss (this note): reading the minus-sheet
comparability `27->5->7` as a structural sign difference; it is the 2-cycle
`{5,7}` at the degenerate target 5. Least-used sidecar: the **target** of an
edge and its residue modulo 3, which turn out to fix the whole spine part of
the Berggren address.

| Live concept | Retained structure | Exact boundary / hostile |
|---|---|---|
| Three Berggren rays from (3,4,5) | A: odd states, C: even states, B: Pell wall | pairwise meet only at the root |
| Three children of a spine point | m-coordinates `n+2`, `(3n+1)/2`, `(3n-1)/2` | second parameters `(n±1)/2` must be kept |
| Address of a Collatz edge | `A^{r(y)}` + letter for `y mod 3` + tail for `(sign,k)` | targets 1 and 5 degenerate |
| Edge slope bands | `k=1: [4,5)`, `k=2: (7,8]`, `k>=3: <=9/4` | the Pell wall `1+sqrt2` separates `k<=2` from `k>=3` |
| Graph triple | girth `(3,4,5)`, chromatic index `(3,4,5)`, three Moore cages | no map to any Collatz object |
| Universal coverage | obstruction set = cycle minima of the rational sheets | no finite 2-adic obstruction set |

Anchor: the exact address of every Collatz edge in the Berggren tree.
Niche: the three-letter graph on odd integers and its Kuratowski thresholds.
Wildcard: the Kuratowski/Tutte "finite obstruction set" template for coverage.

Notation. `U_σ(x)=(3x+σ)/2^{v_2(3x+σ)}` on odd `x`, `σ=±1`; `k=v_2(3x+σ)`
is the halving count of the edge `x->y`. Root pair of an edge:
`(s,t)=(max(x,y),min(x,y))`, odd, coprime. Euclid parameters `(m,n)` with
`(a,b,c)=(m^2-n^2,2mn,m^2+n^2)`, `s=m+n`, `t=m-n`. Berggren letters on Euclid
pairs: `A(m,n)=(2m-n,m)`, `B(m,n)=(2m+n,m)`, `C(m,n)=(m+2n,n)`; these are
`B1,B2,B3` on root pairs. The address of a pair is its word from the root
`(2,1)<->(3,1)<->(3,4,5)`.

## 1. The three rays are the odd states, the even states, and the Pell wall

**1.1 (PROVED).** Iterating one letter from the root gives three rays:

```text
A^j(3,4,5) : c-b = 1,   = Phi(n) = (n,(n^2-1)/2,(n^2+1)/2),   n = 2j+3 odd;
C^j(3,4,5) : c-a = 2,   = (4k^2-1, 4k, 4k^2+1) = Psi(2k),       2k = 2j+2 even;
B^j(3,4,5) : |a-b| = 1, hypotenuses 5, 29, 169, 985, ... (Pell).
```

Pairwise intersections are exactly `{(3,4,5)}`: `c-b=1, c-a=2` gives
`(c-2)^2+(c-1)^2=c^2`, i.e. `c^2-6c+5=0`, `c=5`; `c-b=1, |a-b|=1` forces
`a=c-2` (since `a=c` is impossible) and the same equation; `c-a=2, |a-b|=1`
gives either `b=c-1` (same equation) or `b=c-3`, and `(c-2)^2+(c-3)^2=c^2`
has no integer root. Script `rays_dictionary.py`, D1, checks this
symbolically and on the first 60 members of each ray.

**1.2 (PROVED).** For every `n>=2`, `Psi(n)` is primitive and

```text
n odd  : Psi(n) = A^{(n-3)/2}(3,4,5)      (spine depth (n-3)/2),
n even : Psi(n) = C^{(n-2)/2}(3,4,5)      (mark on the even leg 2n).
```

`Psi(n)` lies on the B ray only for `n in {2,3}`: the root is the marked
triple of **both** 2 (`Psi(2)=(4,3,5)`, mark on the even leg) and 3
(`Psi(3)=(3,4,5)`, mark on the odd leg). Uniformly, `Psi(n)` is the
primitive part of the Gaussian square `(n+i)^2`, and its hypotenuse is
`N(n+i)/d(n)`. The root generator `3+i=(1+i)(2-i)` is the product of the
primes above 2 and 5. Checked for `n<=20000` (D2).

So the owner's three families are exactly: **difference-1 leg = the odd
integers**, **doubled leg = the even integers**, **the family approaching the
isosceles right triangle = the Pell ray, which contains no integer state
except the root**. The Collatz orbit of the root's odd state is
`(3,4,5) -> (5,12,13) -> (1,0,1)`: one A step, then a jump to the degenerate
root (the repository's spine-jump reading, zenodo_triples §3).

## 2. The three children of a spine point are the successor, the 3n+1 shortcut, and the 3n−1 shortcut

**2.1 (PROVED).** Let `n` be odd, so `Phi(n)` has Euclid pair
`(m,n_E)=((n+1)/2,(n-1)/2)`. Its three Berggren children are

```text
A: ((n+3)/2, (n+1)/2)  = Phi(n+2),
B: ((3n+1)/2, (n+1)/2) = (n(2n+1), (3n+1)(n+1)/2, (5n^2+4n+1)/2),   root pair (2n+1, n),
C: ((3n-1)/2, (n-1)/2) = (n(2n-1), (3n-1)(n-1)/2, (5n^2-4n+1)/2),   root pair (2n-1, n).
```

The B child's first parameter is the plus shortcut `T_+(n)=(3n+1)/2`
(zenodo_triples §4); the C child's first parameter is the **minus shortcut**
`T_-(n)=(3n-1)/2`. Direct expansion; checked for all odd `n<=200001` (D3).

**2.2 (PROVED, elementary).** Iterating one letter from `Phi(n)` reads three
affine families in the first parameter:

```text
A^j : n -> n+2j,
B^j : n -> (q_{j+1} n + q_j)/2,   q = 1,1,3,7,17,41,99,...  (q_{j+1}=2q_j+q_{j-1}),
C^j : n -> ((2j+1) n - (2j-1))/2.
```

(Induction: `B` acts as `(m,n')->(2m+n',m)`, so the pair of first parameters
`(m_j,m_{j-1})` satisfies the Pell recurrence with `m_0=(n+1)/2`,
`m_1=(3n+1)/2`; `C^j` adds `j(n-1)` to `m`.) The Collatz shortcuts are the
`j=1` members of the B and C families, and `3` is the only multiplier
`q_{j+1}` or `2j+1` below the halving scale `4`, i.e. the only member of
either family that contracts on average; the next ones are `7n+3`, `17n+7`
(Pell–Collatz) and `5n-3`, `7n-5`. The multiplier ratio `q_{j+1}/q_j`
tends to the silver ratio `1+sqrt2`, the slope of the isosceles limit.
This is a reading of the tree, not a dynamical theorem: the B and C rays
from `Phi(n)` are not orbits.

## 3. The Address Theorem: every Collatz edge is `A^{r(y)}` · (a letter for `y mod 3`) · (a tail for the sign and the halving count)

Define the **spine attachment** of an odd target `y` not divisible by 3 as
the odd integer `n*(y)` among `(y-2)/3, (y+2)/3`:

```text
n*(y) = (y-2)/3  if y = 2 (mod 3),      n*(y) = (y+2)/3  if y = 1 (mod 3),
r(y) = (n*(y)-3)/2  =  (y-11)/6  resp.  (y-7)/6.
```

**Theorem 3.1 (PROVED).** Let `x->y` be a legal edge of `U_σ` with halving
count `k`, `y` odd, `y not in {1,5}`. Then the Berggren address of its root
pair is

```text
plus  (σ=+1):  y=2 mod 3 (k odd) : A^{r(y)} BC          (k=1),
                                    A^{r(y)} BCB A^{4(4^j-1)/3}   (k=3+2j),
               y=1 mod 3 (k even): A^{r(y)} CCC A^{2(4^j-1)/3}   (k=2+2j);
minus (σ=-1):  y=1 mod 3 (k odd) : A^{r(y)} CC          (k=1),
                                    A^{r(y)} CCB A^{4(4^j-1)/3}   (k=3+2j),
               y=2 mod 3 (k even): A^{r(y)} BCC A^{2(4^j-1)/3}   (k=2+2j).
```

In words: climb the spine to `Phi(n*(y))`, i.e. to one third of the target;
take the B child if `y=2 mod 3` and the C child if `y=1 mod 3` (this
letter depends only on `y`); then a C; then the remaining letters record
the sign and the halving count. For a fixed target and sign the parity of
`k` is forced by `(-1)^k y = σ (mod 3)`, so the inverse fibre of `y` is a
single sibling ladder `x -> x + 2^k y = 4x + σ`: it is the B1 ray with base
`((4y-σ)/3, y)` (even `k`) or `((8y-σ)/3, y)` (odd `k>=3`), and in the odd
case the `k=1` edge is the B2 parent of the ladder base.

*Proof.* Write the parent map on root pairs `(s,t)`, `s>t`:
`(s-2t,t)` if `s>3t`, `(t,s-2t)` if `2t<s<3t`, `(t,2t-s)` if `t<s<2t`;
these undo `B1=A`, `B2=B`, `B3=C` respectively. Take `σ=+1`.

*k=1.* `x=(2y-1)/3<y`, pair `(y,x)`, and `x<y<2x` for `y>2`, so the parent
is `(x,2x-y)=(x,(y-2)/3)`. Now `x/((y-2)/3)=(2y-1)/(y-2)=2+3/(y-2)` lies in
`(2,3)` for `y>5`, so the next parent is `((y-2)/3, x-2(y-2)/3)=((y-2)/3,1)`,
a spine point of depth `r(y)=(y-11)/6` provided `(y-2)/3>=3`, i.e. `y>=11`.
Read upward: `B` then `C`. (For `y=5` the pair `(5,3)` is itself the C child
of the root, so the chain has one step and the address is `C`; that is the
exception `3->5`.)

*k=2.* `x=(4y-1)/3`, `y<x<2y`: parent `(y,2y-x)=(y,(2y+1)/3)`; ratio
`3y/(2y+1) in (1,2)`: parent `((2y+1)/3,(y+2)/3)`; ratio
`(2y+1)/(y+2)=2-3/(y+2) in (1,2)`: parent `((y+2)/3,1)`, depth `(y-7)/6`,
valid for `y>=7`. Tail `CCC`.

*k=3.* `x=(8y-1)/3`, `2y<x<3y`: parent `(y,(2y-1)/3)`; ratio
`3y/(2y-1) in (1,2)`: parent `((2y-1)/3,(y-2)/3)`; ratio
`(2y-1)/(y-2)=2+3/(y-2) in (2,3)` for `y>5`: parent `((y-2)/3,1)`. Tail
`BCB`, depth `(y-11)/6`, valid for `y>=11`; at `y=5` the chain reaches the
root after `CB` (the exception `13->5`).

*k>=4, ladder.* `x_{k+2}=(2^{k+2}y-σ)/3=4x_k+σ=x_k+2^k y`, so
`(x_{k+2},y)=B1^{2^{k-1}}(x_k,y)`; the parent chain from `(x_{k+2},y)`
subtracts `2y` while `s>3y`, and since `x_k>y` for `k>=2` (and `y>=3`) it
passes through `(x_k,y)` after exactly `2^{k-1}` A steps. Summing `2^{k-1}`
over the ladder gives the exponents `2(4^j-1)/3` (even) and `4(4^j-1)/3`
(odd). For the odd ladder the base is `(x_3,y)=(2y+x_1,y)=B2(y,x_1)`, which
links the `k=1` pair to the `k=3` pair by one B.

The minus case is identical with `x=(2^k y+1)/3`; the three reductions give
`CC` (via `(x,(y+2)/3)`, `((y+2)/3,1)`), `BCC` (via `(y,(2y-1)/3)`,
`((2y-1)/3,(y-2)/3)`, `((y-2)/3,1)`) and `CCB` (via `(y,(2y+1)/3)`,
`((2y+1)/3,(y+2)/3)`, `((y+2)/3,1)`). The residue conditions `y=1` or
`2 (mod 3)` are exactly the integrality conditions on `x`. ∎

*Exceptions.* Target `1`: the pairs `((4^j-1)/3,1)` (plus) and
`((2^{2j+1}+1)/3,1)` (minus) lie on the A ray itself, addresses `A^{(x-3)/2}`.
Target `5`: plus odd ladder `3->5, 13->5, 53->5, ...` has addresses
`C, CB, CBA^4, CBA^{20}, ...`; minus even ladder `7->5, 27->5, 107->5, ...`
has `CC, CCAA, CCA^{10}, ...`.

*Controls.* `edge_addresses.py`: all 7000 edges with target `<=3001`,
`k<=14`, for each sign: **0 exceptions with `y>=7`**; the 7 exceptions per
sign are the target-5 ladders above. The audit re-derived all 3990 edges
per sign with `y<=2000`, `k<=12`, `y not in {1,5}` from the parent map alone
with 0 mismatches. The orbit script `orbit_comparability.py` uses the parent
map only.

**Corollary 3.2 (comparability, PROVED).** Two distinct edges of the same
sign whose targets are not in `{1,5}` are Berggren-comparable (one an
ancestor of the other) **iff they have the same target**; they then lie on
the single sibling ladder of that target (with its `k=1` edge). Proof:
addresses with different `r(y)` or different first letter are incomparable;
with the same target the tails `BC ≺ BCB ≺ BCBA^4 ≺ ...` (resp.
`CCC ≺ CCCAA ≺ ...`, and the minus analogues) are pairwise prefix-related.
Exhaustive check: all pairs among the 1194 edges with target in `[3,600]`,
`k<=12`, per sign: 2991 (plus) / 2996 (minus) comparable pairs, of which
exactly 6 / 11 are not same-target, and every one of those involves the
target-5 ladders listed above (`edge_addresses.py`, second block).

**Corollary 3.3 (orbits, PROVED).** Two distinct edges of one orbit whose
targets are both outside `{1,5}` are incomparable (an orbit visits each odd
value once, so the targets differ, and 3.2 applies). The exceptions are
exactly of two kinds. Edges into 1, of either sign, lie on the A ray
(`(x,1)` with `x=(4^j-1)/3`, resp. `(2^{2j+1}+1)/3`, address `A^{(x-3)/2}`)
and are comparable with every edge whose spine exponent `r(y)` is at least
`(x-3)/2`; consecutive edges `x->y->1` with `y>=7` are still incomparable,
because `r(y)<=(y-7)/6<(y-3)/2` puts a B or C at a position where the trunk
pair has an A. On the minus sheet the 2-cycle edges
`5->7` and `7->5` share the pair `(7,5)` with address `CC`, which is a prefix
of the even ladder into 5 (`27->5, 107->5, ...`, addresses `CCA^{2(4^j-1)/3}`)
and of the odd ladder into 7 (`19->7, 75->7, ...`, addresses
`CCB A^{4(4^j-1)/3}`); so the comparable consecutive minus pairs are exactly
those of the two types `x->5->7` (`x` in the even ladder of 5) and `x->7->5`
(`x` in the odd ladder of 7), e.g. `27->5->7` and `19->7->5` (audited
exhaustively for `x<=200000`). Hence the consecutive-edge theorem of
ternary_berggren §5 holds for the plus sheet with no exception, and for the
minus sheet exactly outside the 2-cycle `{5,7}`; the sign-decisive example
is the 2-cycle at the degenerate target 5, not a structural difference
between the sheets. Finite control (`orbit_comparability.py`): all orbits
from odd starts `3..20000`, plus to 1 and minus to first cycle entry
(cycle edges outside this universe), 6.90 million resp. 2.82 million edge
pairs; every comparable pair found involves an edge into 1, and a rerun
excluding those edges finds **zero** comparable pairs for both signs.

**3.4 Matrix form: why the sign is the first letter and the halving count
is a tree address (PROVED).** In root coordinates `A=[[1,2],[0,1]]`,
`B=[[2,1],[1,0]]`, `C=[[2,-1],[1,0]]` (determinants `1,-1,1`), and the node
at address `w=w_1...w_l` is `M_w(3,1)^T` with `M_w=M_{w_l}...M_{w_1}`.
Applied to a spine point `(n*,1)` a word `w=Lw'` with first letter
`L in {B,C}` gives, because `A(1,0)=(1,0)`, `B(1,0)=C(1,0)=(2,1)` and
`B(0,1)=-C(0,1)=(1,0)`,

```text
pair = n* . N(w') + eps_L . F(w'),     N(w') = M_{w'}(2,1)^T,  F(w') = M_{w'}(1,0)^T,
eps_B = +1, eps_C = -1,                gamma s' - alpha t' = -det(M_w),  (alpha,gamma)=N(w').
```

`N(w')` is the node at address `w'` of the tree rooted at `(2,1)`. Hence a
spine-attached word gives Collatz edges for **all** `n*` iff `N(w')=(2^k,3)`
(a fall with halving count `k`) or `N(w')=(3,2)` (the rise), and **the sign
`σ` is the choice of the first letter**, since `L` only flips the constant
term. (Individual pairs at the degenerate targets are not of this form:
`(5,3)=C(3,1)` is `3->5` with `N=(2,1)`.) Tree
addresses are unique, so each `k` has exactly one `w'`: the parent-map
address of `(2^k,3)` in the `(2,1)`-rooted tree, `CC A^{2(4^j-1)/3}` or
`CB A^{4(4^j-1)/3}`, and `C` for `(3,2)`. This reproves the tails of 3.1 at
once, explains the plus/minus pairing (`BC|CC`, `CCC|BCC`, `BCB|CCB`, ...),
and shows that every word from a spine point defines some family
`gamma x -/+ 1 = alpha y`; the Collatz families are exactly those with
`(alpha,gamma)=(2^k,3)` or `(3,2)`. Script `spine_words_matrix.py` checks the
identities for all words of length `<=6` and the addresses for `k<=20`.

**Connection contract.** Source: the Collatz edge `x->y` with its sign and
halving count. Target: a Berggren address. Map: odd-Euclid reduction of the
root pair. Preserved: target (as `r(y)` and the first letter), sign and `k`
(as the tail). Destroyed: nothing for targets outside `{1,5}` (the map is
injective there; `5->7` and `7->5` share `CC`). What it
shows: **Berggren ancestry is target-indexed and dynamics-blind**: the tree
sorts the inverse fibres of each `y` next to the spine point `Phi(n*(y))`
at one third of `y`, and the forward orbit `y -> U(y)` appears only as the
change of attachment point `n*(y) -> n*(U(y))`, i.e. as the orbit itself
scaled by `1/3`. (The attachment point is itself a shortcut preimage of a
neighbour of `y/2`: `3n*(y)+1 = y-1` or `y+3`, so `T(n*(y)) = (y-1)/2` or
`(y+3)/2`, whichever is `2 mod 3`; this is exact but is again only the
statement `n* ≈ y/3`.) No new well-founded quantity is produced; this is the
sharp reason (complementing the J-invariant obstruction of ternary_berggren
§6) why Berggren's universal coverage does not transport.

## 4. Slope bands and the Pell wall

For an edge with root pair `(s,t)` write the Euclid slope `σ_E=m/n=(s+t)/(s-t)`.
The near-isosceles triangles have `σ_E -> 1+sqrt2`.

**4.1 (PROVED).** For plus edges, `σ_E=(5x+1)/(x+1) in [4,5)` when `k=1`,
`(7x+1)/(x-1) in (7,8]` when `k=2`, and
`((2^k+3)x+1)/((2^k-3)x-1) <= 9/4` when `k>=3` (equality only at `13->5`).
For minus edges: `k=1: (5,6]`, `k=2: [6,7)`, `k>=3: < 11/5` (the `k=3`
band is `[2, 11/5)`, with `2` at `3->1`). Hence:

- the Pell wall `1+sqrt2 = 2.414...` separates the edges with `k>=3` from
  those with `k<=2`, for both signs; no edge is near-isosceles (this recovers
  the angle bound `tan θ <= 65/72` of odd_square_20260921_edges);
- a 2-cycle `x->y->x` consists of a rising edge (`k=1`) and a falling edge
  (`k>=2`) with the same triangle, hence the same `σ_E`. For plus the rising
  band `[4,5)` is disjoint from both falling bands `(7,8]` and `(0,9/4]`:
  **the plus map has no 2-cycle** (a one-line geometric proof of a classical
  fact);
- for minus the rising band `(5,6]` meets the falling bands only at
  `σ_E=6`, attained only by `5->7` (`(5x-1)/(x-1)=6` forces `x=5`) and
  `7->5` (`(7x-1)/(x+1)=6` forces `x=7`), triangle `(35,12,37)`: **`{5,7}` is
  the only positive 2-cycle of the minus map**.

Verified for all odd `x<10^6`, both signs (D4). The band formulas follow
from `y=(3x+σ)/2^k` by substitution.

## 5. The three edge families on odd integers meet only on the spine segment 1–3–5–7–9

Let `A={x,x+2}`, `B={x,U_+(x)}`, `C={x,U_-(x)}` (edges on odd integers,
loops excluded).

**5.1 (PROVED).** `A∩B={{3,5},{7,9}}`, `A∩C={{1,3},{5,7}}`, `B∩C=∅`.
Proof: `3x+σ=2^k(x+2)` gives `x(3-2^k)=2^{k+1}-σ`, positive only for
`k=1`: `x=4-σ`, i.e. `3->5` (plus) and `5->7` (minus). `3x+σ=2^k(x-2)` gives,
for `k>=2`, `x=2+(6+σ)/(2^k-3)`, integral only when `2^k-3` divides 7
(`σ=+1`: `k=2`, `x=9`, the edge `9->7`) or 5 (`σ=-1`: `k=2`, `x=7`, the edge
`7->5`; `k=3`, `x=3`, the edge `3->1`); `k=1` gives `x<0`. A common B/C edge
`{x,y}` has either the same source (`U_+(x)=U_-(x)`, which forces
`2^b(3x+1)=2^a(3x-1)`, hence `x=1`, the loop) or is a mixed 2-cycle
`3x+1=2^a y`, `3y-1=2^b x`, i.e.
`x(2^{a+b}-9)=3-2^a`; for `a>=3` the right side is negative, forcing
`a+b<=3`, impossible, and `a<=2` gives only `x=y=1`. Exhaustive check
`x<10^6`, including a direct mixed-2-cycle search
(`edge_families_intersection.py`).

So the spine `1-3-5-7-9` carries, alternately, the minus edge `3->1`, the plus
edge `3->5`, the minus 2-cycle `5<->7`, and the plus edge `9->7`; this is the
repository's seven-pair root cluster read on the spine. The three families
meet only at this root cluster, exactly as the three rays meet only at
`(3,4,5)=Phi(3)`.

## 6. The graph triple: what is exact, what is analogy, what is numerology

Exact invariants (`graph_triple_invariants.py`; all CITED classical):

| | V | E | degree | girth | χ | χ′ | α | τ (spanning trees) | Aut | spectrum |
|---|---|---|---|---|---|---|---|---|---|---|
| K_5 | 5 | 10 | 4 | **3** | 5 | **5** | 1 | 5^3 | 120 = S_5 | 4, (−1)^4 |
| K_{3,3} | 6 | 9 | 3 | **4** | 2 | **3** | 3 | 3^4 | 72 | 3, 0^4, −3 |
| Petersen | 10 | 15 | 3 | **5** | 3 | **4** | 4 | 2^4·5^3 | 120 = S_5 | 3, 1^5, (−2)^4 |

**6.1 (CITED, exact).** Among the seven invariants tested (girth, χ, χ′,
degree, α, diameter, cycle rank), exactly two realise a bijection of
the triple onto `{3,4,5}`: the **girth** (K_5, K_{3,3}, Petersen) = (3,4,5),
and the **chromatic index** (K_{3,3}, Petersen, K_5) = (3,4,5). The girth
bijection has structure behind it: the three graphs are the unique Moore
graphs (cages) of parameters (degree, girth) = (4,3), (3,4), (3,5), attaining
the Moore bounds 5, 6, 10. The chromatic-index bijection is the Tutte/Tait
side: K_{3,3} is class 1, Petersen is the smallest snark (class 2, the
obstruction in Tutte's 4-flow conjecture, proved for cubic graphs), K_5 is
class 2 as an odd complete graph. The cubic cage sequence continues
K_4, K_{3,3}, Petersen, Heawood (the Fano incidence graph), McGee,
Tutte–Coxeter with orders 4, 6, 10, 14, 24, 30; the repository's Fano lane
([paley_fano_octonion_design](collatz_mod6_20260922_paley_fano_octonion_design.md))
is the girth-6 member.

**6.2 (structural coincidence, typed ANALOGY).** `Aut(K_5)=Aut(Petersen)=S_5
= PGL(2,F_5)`. The group `PGL(2,Z)` (determinants `±1`), whose positive
monoid is the Berggren tree, reduces modulo the root hypotenuse 5 onto
`PSL(2,5) ≅ A_5` (since `-1=2^2` is a square modulo 5), the index-2
rotation subgroup of `S_5`; the conic `x^2+y^2=z^2` over `F_5` has 6 points,
on which `S_5` acts through its exotic transitive action.
Both K_{3,3} and Petersen have exactly six perfect matchings (K_5 has none),
again `S_6`-sets. These are true statements about the two sides separately;
**no map** carrying a Collatz predicate has been specified, so they are
analogies, not connections.

**6.3 (NUMEROLOGY, recorded as hostile).** Petersen's spectral multiplicities
(1,5,4) with valency 3; `τ(K_{3,3})=3^4`, `τ(K_5)=5^3`,
`τ(Petersen)=4^2·5^3`; the 1-cycle denominators `2^k-3 = 5,13,29,61,125`
are hypotenuses of `Phi(3), Phi(5), B(Phi(3)), Phi(11), ...` for `k<=7` and
fail at `k=8` (`253=11·23`). None of these has a map; they are listed so that
they are not rediscovered as evidence.

**Contract for a future map.** A genuine isomorphism would have to send a
PPT to a graph so that the three Berggren letters become graph operations
with the odd states, even states and Pell ray landing in the complete,
complete-bipartite and snark families, and the Collatz predicate
`3x+1=2^k y` becoming a graph predicate. Section 3 shows what any such map
must reproduce: the target-indexed fan structure of the inverse fibres.
Nothing in sections 1–5 supplies it.

## 7. Universal coverage through the Kuratowski–Tutte template

Kuratowski/Wagner: planarity is a minor-closed property with a finite
obstruction set `{K_5, K_{3,3}}`; Robertson–Seymour: every minor-closed
property has one. Tutte's generation theorems (wheels for 3-connected
graphs; K_4 for cubic 3-connected graphs) and Berggren's theorem (every PPT
from `(3,4,5)`) are *universal coverage* theorems proved by a well-founded
descent. The Collatz conjecture is a universal coverage statement for the
inverse tree from 1. What the template gives and does not give:

**7.1 (PROVED, elementary).** There is no finite 2-adic obstruction set: for
every `k`, the class `n = -1 (mod 2^{k+1})` has `k` consecutive rises
(`n = -1 mod 2^{j+1}` gives `n = 3 mod 4` and `U(n)=(3n+1)/2 = -1 mod 2^j`),
so its accelerated stopping time exceeds `k`, and the witness is sharp:
`n=2^k-1` rises only `k-1` times. The proportion of such classes is `2^{-k}`
of the odd classes (exactly 2 of `2^{k+1}` at every `k<=14`,
`obstruction_sheets.py`; CITED Terras 1976 for density zero of infinite
stopping time). Hence no certificate that inspects
a fixed number of low bits can be universal; the repository's SHEET no-go
(any proof must use a hypothesis true for `3n+1` and false for `3n-1`) is
the sign-law form of the same fact.

**7.2 (framing, with FINITE-EXACT census).** On the odd rationals the map
`U_+` preserves the prime-to-3 part `q` of the reduced denominator
(ternary_triples §3), and on the sheet of denominator `q` it is the `3n+q`
map on odd integers. The obstructions to coverage on sheet `q` are its
cycles; the 2-adic obstruction set intersected with the rationals is the set
of cycle minima over all sheets. Census of cycle minima (odd seeds `<=10^5`,
no overflow): `q=1: {1}`; `q=5: {1,5,19,23,187,347}`; `q=7: {5,7}`;
`q=11: {1,11,13}`; `q=13:` ten cycles; every `q<=65` coprime to 6 has at
least one proper cycle (minimum coprime to `q`). The 1-cycles
`x=1/(2^k-3)` live on the sheets `q=2^k-3=1,5,13,29,61,...`, all sharing
the shape of the root; their unoriented triangles are the reflected spine
points `Phi(2^k-3)` (D5), and the two integer members are `k=1` (the
negative root `-1`) and `k=2` (the root `1`). In the Kuratowski–Tutte
reading, the "Petersen" of the rational system is the first proper
obstruction `1/5` with triangle `(5,-12,13)`, the reflection of the A child
of the root. Background CITED (bibliographic data checked, texts not
re-read this session): Lagarias, [*The set of rational cycles for the 3x+1
problem*, Acta Arith. 56 (1990), 33–53](https://eudml.org/doc/206298), and
Belaga–Mignotte, [*Embedding the 3x+1 conjecture in a 3x+d context*,
Experimental Math. 7:2 (1998)](https://projecteuclid.org/journals/experimental-mathematics/volume-7/issue-2/Embedding-the-3x1-conjecture-in-a-3xd-context/em/1048515662.full);
neither is used for a proof step here.

**7.3 (REFUTED as a mechanism; FINITE-EXACT).** Planarity of the three-letter
graph. Let `Γ_N` be the graph on odd `1..N` with the A, B, C edges of
section 5. `A+B` is first nonplanar at `N=29` (K_{3,3} subdivision with branch
vertices 13,17,19,21,23,25), `A+C` at `N=35`, `B+C` at `N=37`, `A+B+C` at
`N=19`; the pure B graph (a forest) and the pure C graph (a forest plus one
7-cycle through 17; the 2-cycle `{5,7}` is a single edge and the fixed point
1 a loop, excluded) are planar for all odd `N<4001`. Since A is a
Hamiltonian path, `A+B` together with the closing edge `{1,N}` is planar iff
the Collatz chords admit a two-page book embedding in the natural cyclic
order (without the closing edge the implication is one-way: at `N=25` and
`N=27`, `A+B` is planar although the chords `(1,21),(15,23),(19,25)` pairwise
cross). The chromatic number of the crossing graph of the `3n+1` chords is
3, 4, 8 at `N=31, 63, 127` and at `N=255` the crossing graph contains a
15-clique, so the book thickness of the Collatz chord system in the natural
order grows without any sign of a bound. Planarity thresholds near `15`/`16`
(decoder_minors) and `19`/`29` here are numerical juxtapositions with no
transported predicate.

**7.4 What survives.** The Address Theorem gives a complete, sign-symmetric
description of the Collatz edge set inside the Berggren tree; it explains
the previous incomparability results and shows the tree cannot see the
orbit. The Pell wall separates the halving-count bands but is never touched.
The three families meet only at the root, on both carriers. The obstruction
set of the rational extension is infinite and parametrised by periodic
parity words; coverage of the integer sheet is exactly the statement that
this set meets the positive integers only at 1. A proof therefore cannot be a
finite-obstruction argument at any 2-adic resolution, and cannot use
Berggren ancestry as a rank. Open obligations, unchanged: a sign-sensitive
well-founded rank on the integer boundary, or a common-future certificate
with a smaller partner for every `n = 3 (mod 4)`.

## 8. Reproduction

From the repository root:

```text
python3 04-computation/experiments/cowork_20260925_kuratowski/rays_dictionary.py
python3 04-computation/experiments/cowork_20260925_kuratowski/edge_addresses.py 3001 14
python3 04-computation/experiments/cowork_20260925_kuratowski/orbit_comparability.py 20000
python3 04-computation/experiments/cowork_20260925_kuratowski/edge_families_intersection.py
python3 04-computation/experiments/cowork_20260925_kuratowski/graph_triple_invariants.py
python3 04-computation/experiments/cowork_20260925_kuratowski/obstruction_sheets.py
python3 04-computation/experiments/cowork_20260925_kuratowski/spine_chords_planarity.py
python3 04-computation/experiments/cowork_20260925_kuratowski/spine_words_matrix.py
```

Frozen outputs: `*.out` next to the scripts. Universes are explicit in each
script; no sampling. `rays_dictionary.py` uses sympy for the two quadratic
equations of 1.1 only; `graph_triple_invariants.py` and
`spine_chords_planarity.py` use networkx (planarity test with counterexample,
VF2 automorphism count); everything else is standard-library integer
arithmetic.

## 9. Independent audit and correction lineage

A read-only adversarial audit (separate agent, own scripts under `/tmp`,
no imports from this lane) re-derived the Address Theorem for all 3990
edges per sign with `y<=2000`, `k<=12`, `y not in {1,5}` (0 mismatches),
checked every reduction inequality of the proof of 3.1, the ladder lemma,
the slope bands with their boundary values, the family intersections to
`10^5`, the children formulas to `n<=1001`, the graph table, the Moore
bounds, the `2^k-3` hypotenuse list, and the sheet census. It found the
following defects in the first draft, all repaired above:

- the sibling ladder was written `x->4x+1` for both signs; the sign-free
  form is `x->x+2^k y = 4x+σ`;
- "for `y=5` the chain stops at `(1,1)`" was wrong: `(5,3)` is the C child of
  the root;
- 3.2's "same parity of halving count" is vacuous (the parity is forced by
  the target and the sign) and its control numbers were misreported;
- 3.3 omitted that edges into 1 are exceptions on the minus sheet as well,
  and that the odd ladder into 7 (`19->7->5`) is comparable with the cycle
  edge `7->5`;
- 5.1 had the mixed-2-cycle equation with the wrong sign (`2^a-3` for
  `3-2^a`) and a wrong closed form; the script enumerated the wrong equation
  and was corrected to the right one plus a direct search;
- 6.2 claimed `PGL(2,Z)` surjects onto `PGL(2,5)`; the image is `PSL(2,5)`;
- 7.1 was off by one: `k` rises need `n = -1 mod 2^{k+1}`, not `2^k`;
- 7.3 stated a false "iff" (a Hamiltonian path plus chords can be planar
  with a non-bipartite crossing graph; the equivalence needs the closing
  edge) and used "unbounded in this range".

A second audit re-verified the repairs and the new matrix form (3.4): the
identities for all words of length `<=8`, the `(2,1)`-tree addresses for
`k<=22`, the exhaustive list of comparable consecutive minus pairs to
`x<=200000`, the band and 2-cycle statements, the `PSL(2,5)` image, and the
planarity equivalence with the closing edge for all odd `N<=401`. It added
the same-source case to 5.1 and the "for all `n*`" qualifier to 3.4.

Repository correction lineage: none of the inherited theorems is changed.
The reading "the sign is decisive" in ternary_berggren §5 stands as a true
statement; this note supplies its mechanism (the 2-cycle at target 5) and
the sign-symmetric general structure.