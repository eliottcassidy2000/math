# Exact boundary of the prime-circulant unit-power-sum detector

Status: **PROVED, independently audited PASS, with exact certificates.**
The independent review is in [independent_audit.md](independent_audit.md).

Scope: labelled tournament connection sets on the additive group of an odd
prime field, and their multiplicative stabilizers. “Rigid” below means
**multiplier-rigid**. Every circulant tournament still has its translation
automorphisms; this note does not call its entire automorphism group trivial.

This extends Proposition D of
[THM-4724, Proposition D](../../../01-canon/theorems/THM-4724-affine-fixed-point-criterion-restricted-twisted-burnside-for-every-submodule-even-n-switching-half-and-the-dc-obstruction-to-uniform-circulant-rules.md).
Its two positive regimes are retained, the remaining regimes are decided,
and the two-odd-prime census is computed in closed form. No assertion here
depends on Theorem W's arithmetic argument, a general local-rule no-go
theorem, or complex Fourier analysis.

## 1. Statement

Let `p` be an odd prime, `n=p−1=2^e m`, with `m` odd and `e≥1`. A tournament
connection set is a subset `S⊂F_p^×` with `S⊔(−S)=F_p^×`. Put

`σ_k(S)=Σ_{s∈S}s^k ∈ F_p`.

Call `S` traceless if `σ_k(S)=0` whenever `gcd(k,n)=1`, and put
`K(S)={a∈F_p^×:aS=S}`.

**Theorem 1 (exact boundary).** The implication

`S traceless ⇒ K(S)≠{1}`

holds for every tournament connection set `S` if and only if one of the
following mutually exclusive conditions holds:

1. `m=1`;
2. `m=q^a` for one odd prime `q` and `a≥1`;
3. `e=1` and `m=qr` for two distinct odd primes `q,r`.

The reverse implication `K(S)≠{1} ⇒ S traceless` holds at every odd prime.
Thus these are exactly the primes where tracelessness detects nontrivial
multiplier symmetry.

In case 1 no traceless tournament set exists and every multiplier
stabilizer is trivial. In case 2 tracelessness is equivalent to invariance
under the unique subgroup of order `q`.

**Theorem 2 (all two-odd-prime counts).** If

`p−1=2^e q^a r^b`, where `q,r` are distinct odd primes and `a,b≥1`, define

`B=2^(e−1) q^(a−1) r^(b−1)` and `C=2^q+2^r−2`.

The number of labelled traceless tournament connection sets is exactly

`N_tr=C^B`.

The number of labelled multiplier-rigid traceless sets is exactly

`N_rigid=C^B−2^(rB)−2^(qB)+2^B`.

There is no quotient by multipliers or tournament isomorphisms in these
counts. In particular, rigidity is impossible exactly when `B=1`.

## 2. Inheritance and concept board

- Closest inherited mechanism: the `F_p`-Fourier description in
  THM-4724 Proposition D, and its two-coordinate Boolean-array lemma.
- Canonical hostile: `p=61`, where different antipodal pairs can select
  different prime-coordinate directions.
- Corrected near miss: “at most two odd primes” is safe only when the
  repeated exponents and the two-adic slices have also been controlled.
- Least-used sidecar: the lower digits of a primary cyclic coordinate
  survive the prime-order projector and supply extra independent blocks.
- Live comparison board: unit Fourier support; prime-order averaging
  projectors; antipodal/lower-digit blocks; multiplier rigidity;
  square-class-sensitive normalization.

The central connection is:

`tournament sign function → independent primary top-digit blocks`.

It preserves the exact unit-power-sum condition and prime-order
stabilizer tests. Merely counting distinct odd primes loses the
two-adic slice and lower-digit block labels. Keeping those labels supplies
both the corrected boundary and the exact census.

## 3. The projector identity over F_p

Choose a primitive root `g` of `F_p^×`, identify `F_p^×` with `C_n` by
`i↦g^i`, and write `f(i)=1` for `g^i∈S`, and `f(i)=−1` otherwise. Then

`f(i+n/2)=−f(i)`.

For `k mod n`, use the Fourier sum

`fhat(k)=Σ_{i mod n} f(i)g^(ki)`.

For odd `k`, `fhat(k)=2σ_k(S)`; for even `k`, `fhat(k)=0` by antipodality.
The Fourier matrix is invertible over `F_p` because its roots are distinct
and `p∤n`.

For each odd prime `q|m`, let `Q_q≤C_n` be the subgroup of order `q` and
let

`A_q f(i)=q^(−1) Σ_{u∈Q_q} f(i+u)`.

On Fourier mode `k`, this projector is the identity if `q|k` and zero
otherwise. The projectors commute. Therefore, on odd functions,

`f is traceless ⇔ [∏_{q|m}(I−A_q)]f=0`.                     (3.1)

Equivalently, a traceless odd function belongs to the sum of the
`Q_q`-invariant subspaces. One can see the latter without a dimension
argument by expanding `I−∏(I−A_q)`; each resulting term has image in at
least one `Q_q`-invariant subspace.

The product over an empty prime set is the identity. Thus for `m=1`,
tracelessness would force `f=0`, which contradicts `f∈{−1,1}`.

If `m=q^a`, (3.1) says `A_q f=f`, exactly invariance under `Q_q`.

For the universal easy direction, if `aS=S` with `a≠1` and `gcd(k,n)=1`,
then `σ_k(S)=a^kσ_k(S)`, while `a^k≠1`; hence `σ_k(S)=0`.

## 4. Exact primary blocks and the two-prime count

Use the CRT coordinates

`C_n = C_(2^e) × C_(q^a) × C_(r^b)`.

Write the last two coordinates uniquely as

`x=x_0+q^(a−1)u`, with `0≤x_0<q^(a−1)`, `u∈C_q`,

`y=y_0+r^(b−1)v`, with `0≤y_0<r^(b−1)`, `v∈C_r`.

The group `Q_q` changes `u` and fixes `t,x_0,y_0,v`; similarly `Q_r`
changes `v` and fixes the other coordinates. Consequently (3.1) acts
separately on each top-digit `q×r` rectangle with fixed `t,x_0,y_0`.
Antipodality pairs the `t` and `t+2^(e−1)` rectangles as negatives.
There are exactly `B=2^(e−1)q^(a−1)r^(b−1)` independent rectangles.

On one rectangle, `(I−A_q)(I−A_r)F=0` says

`F(u,v)=U(u)+V(v)`,

by subtracting the two coordinate averages and adding the grand average.
Conversely every such additive rectangle lies in this kernel.

**Boolean rectangle lemma.** Over any field of characteristic not two,
if `U(u)+V(v)∈{−1,1}` at every entry, then `U` or `V` is constant.

Indeed, `F(u,v)^2=1`. Its mixed difference at `u,u',v,v'` is

`2[U(u)−U(u')][V(v)−V(v')]=0`.

If `U` has two distinct values, this forces every difference of `V` to
vanish. This proof does not require an ordering, positivity, or an
unspoken characteristic-zero hypothesis.

Thus each rectangle is an arbitrary sign function of `u` alone or of `v`
alone. There are `2^q+2^r−2=C` choices: only the two constant rectangles
are counted in both families. Independent blocks give `N_tr=C^B`.

A tournament multiplier stabilizer has odd order, since a subgroup of
even order in the cyclic group `F_p^×` contains `−1`, which cannot fix
`S`. Every nontrivial subgroup of odd order contains a prime-order
subgroup. With only the odd primes `q,r`, a set is nonrigid precisely
when it is `Q_q`-invariant or `Q_r`-invariant.

The `Q_q`-invariant tournament sets number `2^(n/(2q))=2^(rB)`.
The `Q_r`-invariant ones number `2^(qB)`. Their intersection consists of
the sets invariant under the subgroup of order `qr`, and numbers `2^B`.
All these sets are traceless by the universal easy direction. Subtracting
this union from the traceless population proves Theorem 2.

For `B=1`, the count of rigid sets is zero. For `B≥2`, choose one block
to be a nonconstant function of `u` and a second block to be a
nonconstant function of `v`; fill every other block with a constant.
The first chosen block breaks `Q_q`, and the second breaks `Q_r`.
Thus the rigid count is strictly positive, independently of merely
inspecting its closed formula.

This proves both the positive `e=a=b=1` case and every negative
two-odd-prime case, retaining all prime-power factors.

## 5. Explicit formulas and prime-order stabilizer witnesses

For a tuple `z`, define the point-sign

`η(z)=+1 if every coordinate of z is zero, and −1 otherwise`.

The empty tuple has point-sign `+1`. Each sign function below gives
`S={g^i:f(i)=+1}`.

Write

`C_n=C_(2^e)×∏_{j=1}^s C_(q_j^(a_j))`,

with coordinates `(t,x_1,...,x_s)`, and let `h=2^(e−1)`.
The unique subgroup of order `q_j` shifts only `x_j`; it is generated
there by the nonzero increment `d_j=q_j^(a_j−1)`.

### 5A. Extra antipodal slices: e≥2 and s≥2

For `0≤t<h`, set

`F(t,x)=η(x_2,...,x_s)` if `t=0`,

`F(t,x)=η(x_1,x_3,...,x_s)` if `1≤t<h`,

and define `f(t+h,x)=−F(t,x)`, `f(t,x)=F(t,x)`.

The first kind of slice is invariant under `Q_(q_1)` and the second
under `Q_(q_2)`. Slice indicators depend only on `t`, so the complete
function is a sum of a `Q_(q_1)`-invariant function and a
`Q_(q_2)`-invariant function. It is therefore traceless.

At `t=1`, shifting `x_1` from zero to `d_1` changes its sign. At `t=0`,
shifting any `x_j`, `j≥2`, from zero to `d_j` changes its sign. Every
odd prime-order subgroup is broken; antipodality excludes even-order
stabilizers. Hence the multiplier stabilizer is trivial.

### 5B. Lower-digit selector: e=1, s=2, one repeated odd prime

Rename the factors as `C_(q^a)×C_(r^b)` with `b≥2`; put
`d=r^(b−1)`. Define

`F(x,y)=η(y/d)` if `y≡0 mod d`,

`F(x,y)=η(x)` otherwise,

and `f(t,x,y)=(−1)^t F(x,y)`.

The first branch, including its support indicator, is independent of
`x`, hence `Q_q`-invariant. The second branch is constant on the
order-`r` orbits in `y`, because its selector depends only on `y mod d`;
hence it is `Q_r`-invariant. Their sum is traceless.

The shift `y:0→d` changes the sign in the first branch and breaks
`Q_r`. On the slice `y=1` (which exists outside the selected fiber
because `d≥r≥3`), the shift `x:0→q^(a−1)` changes the sign and breaks
`Q_q`. This proves rigidity for all `a≥1,b≥2`, without inflating a
smaller witness through a quotient.

### 5C. Third primary-coordinate selector: e=1 and s≥3

Define

`F(x)=η(x_2)η(x_4,...,x_s)` if `x_3=0`,

`F(x)=η(x_1)η(x_4,...,x_s)` if `x_3≠0`,

and `f(t,x)=(−1)^t F(x)`.

The first term is independent of `x_1`; the second is independent of
`x_2`. Their support indicators are also independent of those
coordinates. Thus `f` lies in the sum of the `Q_(q_1)`- and
`Q_(q_2)`-invariant subspaces and is traceless.

The following explicit witnesses break all prime-order subgroups:

- `Q_(q_1)`: choose `x_3=1`, all other coordinates zero; change
  `x_1:0→d_1`.
- `Q_(q_2)`: choose all coordinates zero; change `x_2:0→d_2`.
- `Q_(q_3)`: choose `x_1=0`, `x_2=d_2`, and the other coordinates zero;
  change `x_3:0→d_3`. The two branch values are `−1` and `+1`.
- For `j≥4`: choose all coordinates zero and change `x_j:0→d_j`;
  the extra point-sign factor changes from `+1` to `−1`.

All exponents `a_j` can be arbitrarily large. No unused coordinate or
quotient-lift kernel is hidden in this formula.

Theorem 1 now follows: the three positive regimes were proved in
sections 3–4; all remaining factorizations fall under 5A, 5B, or 5C.

## 6. Representative primes and exact counts

The full explicit connection sets and every unit-index power sum are in
[verification.json](verification.json), produced independently of the proof by
[verify_boundary.py](verify_boundary.py).

| p | p−1 | mechanism | labelled traceless count | labelled rigid traceless count |
|---|---|---|---:|---:|
| 61 | `2²·3·5` | antipodal slices | 1,444 | 360 |
| 127 | `2·3²·7` | lower-digit selector | 2,406,104 | 308,448 |
| 151 | `2·3·5²` | lower-digit selector | 79,235,168 | 45,648,000 |
| 271 | `2·3³·5` | lower-digit selector | 165,216,101,262,848 | 130,031,594,956,800 |
| 211 | `2·3·5·7` | third primary coordinate | not enumerated | at least one explicit witness |
| 331 | `2·3·5·11` | third primary coordinate | not enumerated | at least one explicit witness |

The counts for 127, 151, and 271 are **PROVED formula evaluations**, not
claims of exhaustive enumeration. The count at 61 is both proved and
independently **FINITE-EXACT** by a complete `15+15` meet-in-the-middle
census of all `2^30` connection sets.

The independently supplied original 151 witness uses primitive root 6
and agrees with formula 5B for `C_2×C_3×C_25`.

The smallest prime where the detector fails is **61**. The smallest
prime congruent to 3 modulo 4 where it fails is **127**; 211 remains a
valid later example.
These minimality statements follow from the complete theorem and the
finite factorizations of smaller primes, not a search through all
tournaments of order below 127.

The incoming note's “p=151,331 ... with a square odd prime factor”
parenthetical is arithmetically inaccurate: `331−1=330=2·3·5·11` is
squarefree. It belongs to the third-distinct-prime regime.

## 7. Exact computation and audit envelope

Reproduce from this directory:

`python3 verify_boundary.py --bound 1000 --output verification.json`

The script uses only Python's standard library and exact integers
modulo `p`. It performs:

1. A complete Gray-code census of all `2^((p−1)/2)` tournament
   connection sets for every odd prime `3≤p≤31`. Traceless counts are
   respectively `0,0,2,2,4,0,8,2,4,38` at
   `p=3,5,7,11,13,17,19,23,29,31`; all rigid-traceless counts are zero.
2. A complete exact `p=61` meet-in-the-middle census. It gives 1,444
   traceless sets: 360 of stabilizer order 1, 1,020 of order 3,
   60 of order 5, and 4 of order 15.
3. One explicit construction for **every one of the 62 bad-regime
   primes at most 1,000**. On each constructed set it computes all
   unit-index power sums directly from actual residues, checks the
   tournament property, and tests every multiplicative stabilizer.
   It separately checks every additive period of the exponent sign
   sequence and compares the two stabilizer calculations.

Step 3 is a **complete prime census of the displayed construction
within its bound**, not a census of all connection sets at those
primes. The structural proofs, rather than this finite check, cover
arbitrary prime powers and all primes.

## 8. Scope corrections and non-consequences

The last sentence of incoming THM-4724 section 4 says exactly:

> Theorem W shows that no local rule does this; the normalised matching problem is the honest residue of HYP-9341 (c).

Section 3 of the same file expressly leaves general dilation-rule
non-existence open. The quoted sentence therefore broadens the stated
result without justification. The narrow survivor is:

“Theorem W addresses the derivative family `D_c`; the forcing analysis
constrains general dilation rules, whose non-existence remains OPEN.”

This note records that scope slip rather than silently altering the
incoming source. Its own boundary theorem proves neither the
derivative-family arithmetic claim nor a general local-rule no-go
theorem, and supplies no bijection between tournaments and even graphs.

For rigid tournament sets in the good regimes, some unit-index
coefficient is nonzero. The least such index `k` in `1≤k<p−1` and the unique solution
of `a^kσ_k(S)=1` do canonically normalize each multiplier orbit: power by
`k` is a permutation of `F_p^×`. The theorem identifies exactly when this
particular detector is guaranteed to offer such an index for every
rigid tournament set.

For symmetric graph connection sets, coefficients with
`gcd(k,p−1)=2` scale only through squares. A corresponding normalization
needs the coefficient's square class and is unique only modulo `±1`;
one cannot silently normalize every nonzero coefficient to 1. For
example, at `p=7`, the even circulant connection set
`E={±1,±2}` has multiplier stabilizer `{±1}` and `σ_2(E)=3`, a nonsquare.
Its complement is a 7-cycle, so its automorphism group is dihedral;
each reflection reverses exactly two fixed edges of `E`, which verifies
the required even-graph condition.
Thus `a²σ_2(E)=1` has no solution. Choosing a fixed target in each
square class repairs this normalization step when such a coefficient
exists; this note does not prove that it always exists for rigid
symmetric graph sets.

## 9. Mechanism, failure boundary, and next question

The failure mechanism is not mysterious modular cancellation. With
exactly two odd primes, primitive Fourier annihilation is an additive
top-digit rectangle condition. Booleanity forces one direction on each
rectangle, but it does not force the same direction on different
rectangles. The missing label was the independent block address.

Extra antipodal pairs (`e≥2`) and lower primary digits (`a≥2` or `b≥2`)
both create such extra block addresses. A third distinct odd prime
supplies a selector even with one antipodal pair and no repeated
prime. These are precisely the three ways to escape the good regimes.

The strongest survivor of the original detector is its exact boundary
and the two-prime census above. A useful subsequent question is a
closed census of Boolean arrays in the sum of coordinate-omitting
subspaces for three or more distinct odd primes; this is separate
from, and not needed for, the existence/rigidity classification.
