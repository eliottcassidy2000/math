# The fixed sporadic map has a global level-three pole-seven norm and one new prime image

**Status: SUPERSEDED AS THE CURRENT FRONTIER BY PROVED THM-3495; retained as
the independent pole/residue derivation.**  The proof uses the fixed-map geometry in
[THM-2473](../01-canon/theorems/THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy.md),
[THM-2576](../01-canon/theorems/THM-2576-composite-jelonek-image-divisor-and-two-component-nonproperness-law.md),
and [THM-2582](../01-canon/theorems/THM-2582-odd-block-discriminant-tower-and-composite-jelonek-square-class.md),
plus the pinned finite-exact slice computation.  THM-3495 subsequently froze
the global primitive normalization and closed the degree-27 separability gate.

## Typing correction and inheritance

The three pinned `N(H)` slices belong to the fixed sporadic generic-degree-three
map `F` of
[THM-1300](../01-canon/theorems/THM-1300-jacobian-counterexample-dixmier-A3-explicit.md),
not to the degree-four weighted-lift `G` called the quartic G1 witness in
[THM-3438](../01-canon/theorems/THM-3438-weighted-lift-keller-degree-spectrum.md).
The latter already resolves the G1 existence gate but has different inverse
and Jelonek data.  Nothing below is transported to that quartic map.

The closest proved mechanism is THM-2582's discriminant-of-a-norm law.  Its
canonical hostile is the false level-two guess that the old factor `L` remains
odd: the exact identity `N(L)=H/(64L)` cancels it.  The corrected near miss at
level three was to promote three irreducible slices to a global divisor law.
The least-used sidecar is the valuation at the **two escaping inverse sheets**
over the generic point of `V(L)`.  That sidecar is enough to globalize the
denominator without expanding the new numerator.

The live board was:

1. the function-field norm `N(H)`;
2. the escaping pair of the depressed cubic;
3. its unique finite survivor;
4. the prime image of `V(H)` and its generic multiplicity; and
5. the reduced nonproperness-set tower.

## Exact global statement

Put `R=Q[a,b,c]`, `K=Frac(R)`, and

```text
L = 27a^2c^2 - 18abc + 16a + b^3c - b^2,
T = 4 - 3bc,
E(w) = Lw^3 + Tw - 2c.
```

Let `H` be THM-2576's primitive irreducible 361-term equation of
`closure(F(V(L)))`, and let `N` be the degree-three norm from the generic
inverse algebra `K[w]/(E)`.  Define

```text
J := L^7 N(H).                                             (1)
```

This reflection's local `J` is the rationally normalized `J_res`.  In
THM-3495's primitive integral notation,

```text
J_prim=2^35 J_res,      J_res=J_prim/2^35.                (1a)
```

All hypersurface statements below are scale-invariant.  Square classes must
use (1a): `[-J_res]=[-2J_prim]`.

Then the proof below gives:

```text
J lies in Q[a,b,c] and L does not divide J;                (2)
J is irreducible over Q;                                   (3)
closure(F(V(H))) = V(J);                                   (4)
the generic divisor/image multiplicity is one;             (5)
J is distinct from both L and H.                            (6)
```

Thus the exact global norm shape is

```text
N(H) = J/L^7,                                              (7)
```

where the normalization in (1) absorbs the rational scalar.  This obtains the
global factor structure without printing the potentially large coefficient
ledger of `J`.

Combining (4)--(6) with THM-2576's all-iterate set law gives the reduced
third-iterate nonproperness set

```text
S_(F^3) = V(L) union V(H) union V(J) = V(L H J),           (8)
```

with exactly three irreducible components.  This is a statement about one
fixed map and one iterate grade, not a Keller-family classification.

## Why the pole is exactly seven

THM-2473 proves that `F` is finite etale of degree three away from `V(L)`.
Consequently `N(H)` belongs to `R[1/L]`: no irreducible divisor other than
`L` can occur in its denominator.  It remains to compute the single valuation
`v_L`.

At the generic point of `L=0`, with `A=12a-b^2` and `AT != 0`, the cubic
has one finite root and two escaping roots:

```text
w_0 -> 2c/T,
w_1,w_2 -> infinity,
L w_1 w_2 -> T.                                           (9)
```

The exact inverse graph is

```text
q_x = w,
q_y = b - 3aw((9ac-b)w+2)/(Aw^2+bw+2),
q_z = (2w-c-3w^2q_y)/w^3.                                (10)
```

The companion verifies all three rows `F(q(w))=(a,b,c)` modulo `E`.  Along an
escaping root,

```text
q_y -> eta := (15ab-b^3-27a^2c)/A,
q_z = -3eta/w + O(w^-2).                                 (11)
```

Among the 361 monomials `a^i b^j c^k` of `H`, the maximum of `i-k` is
exactly seven.  Its four-term face evaluates exactly to

```text
H(w,eta,-3eta/w) / w^7
  -> C eta^3,
C = 83958031872 = 2^9 3^6 11^3 13^2.                    (12)
```

The finite survivor contributes `H(q(w_0))`, which is nonzero in the generic
function field of `V(L)`.  Equations (9)--(12) therefore give

```text
L^7 product_(i=0)^2 H(q(w_i))
 -> C^2 T^7 eta^6 H(q(w_0)) != 0.                        (13)
```

Hence

```text
v_L(N(H)) = -7.                                           (14)
```

Together with `N(H) in R[1/L]`, (14) proves (2) and (7).  The mechanism is
not “three slices agree”: it is the even escaping pair, each carrying the
order-seven face with half-integral `L`-valuation, plus one nonvanishing
finite sheet.

## Exact boundary residue on the normalization

On THM-2576's normalization

```text
nu(tau,lambda) = (
  lambda^2(3-tau lambda)/27,
  lambda(4-tau lambda)/3,
  tau
),                                                        (15)
```

the quantities controlling (13) simplify to

```text
A o nu   = -lambda^2(lambda tau-2)^2/9,
T o nu   =  (lambda tau-2)^2,
eta o nu =  lambda/3.                                     (16)
```

The finite inverse point is

```text
q_0 = (
  2tau/(lambda tau-2)^2,
  -lambda(3lambda tau-8)/6,
  -lambda^2(lambda tau-2)^2
    ((lambda tau)^2-8lambda tau+14)/8
).                                                        (17)
```

Exact substitution into the transported 361-term `H` gives

```text
H(q_0) = lambda^2 P(tau,lambda)
         / [2^35 3^21 (lambda tau-2)^14],                 (18)
```

where `P` is primitive with

```text
deg_tau P = 57,  deg_lambda P = 84,
total degree = 141,  terms = 527,
coefficient-ledger SHA-256
  5f15ad89e0788eafcdffecf6fa0f0204224d27fb6ecf94cd46ed1982da678333.
```

FLINT factors `P` over `Z` as one exponent-one factor, so `P` is irreducible
over `Q`.  Substituting (16)--(18) into (13) cancels the entire exceptional
seam power and yields the exact residue identity

```text
nu^*(J restricted to L=0)
 = [11^6 13^4/(2^17 3^15)] lambda^8 P(tau,lambda).        (19)
```

The root split itself excludes the hostile seam `lambda tau=2`, where both
`A` and `T` vanish.  Equation (19), derived on the dense lawful chart, extends
regularly across it; the exact control `P(1,2) != 0` shows that the canceled
formula does not silently turn that seam into the whole residue.

## Why the numerator is one prime, not a hidden power

Over `U=Spec R minus V(L)`, the norm vanishes exactly when one point of the
finite etale inverse fibre lies on `V(H)`.  Since `H` is irreducible and its
intersection with `F^{-1}(U)` is nonempty, its finite image is one irreducible
prime divisor, say `V(J_0)`.  Therefore

```text
J = u J_0^e                                                (20)
```

for a rational unit `u` and an integer image multiplicity `e>=1`.

The pinned exact slice `(b,c)=(1,2)` gives

```text
J(A,1,2) = K_(1,2)(A)/2^21,
deg K_(1,2)=86,
K_(1,2) irreducible.                                      (21)
```

If `e>1`, the nonconstant specialization in (21) would remain an `e`th power,
contradicting its exponent-one irreducibility.  Thus `e=1`, proving (3)--(5).
This use of the slice is lawful: globality and unique support come first from
the finite-map/divisor argument; the slice supplies only the remaining
integer multiplicity.

A small exact point gives independent rank and distinctness controls:

```text
(0,-1,-1) in V(L),
F(0,-1,-1) = (3,-1,0) in V(H),
F(3,-1,0) = (10,-46,33) in V(J),

L(10,-46,33) = -504,
H(10,-46,33) = -1402696598666966597632.                  (22)
```

The gradient of `H` at `(3,-1,0)` is nonzero, and `det J_F=-2`; hence the
restriction has rank two there.  Equation (22) separates the new surface from
both older components without relying on degree numerology.

## Discriminant consequence and the gate closed by THM-3495

The norm part of the level-three recursion is now exact:

```text
[-L N(H)] = [-J_res/L^6] = [-J_res]
           =[-2J_prim] in K^*/K^{*2}.                    (23)
```

Thus the old `L` exponent is the even value `1-7=-6`, and `H` is absent from
the norm square class.  This reflection left open whether a selected
degree-27 `x`-coordinate separates the generic points.  THM-3495 closes that
gate with a full-degree squarefree specialization and an off-grid determinant,
so (23) is now the actual eliminant square class, not merely the norm class.

Still open or out of scope after THM-3495:

- the depth-four norm divisor and an all-level newest-factor law;
- the corresponding norm tower for THM-3438's quartic weighted map `G`;
- arbitrary Keller-map, atom, conjugacy, or counterexample classification;
- `JC(2)`, `DC(1)`, and `DC(2)`;
- any LRC, tournament, or physical-current transfer.

This obeys the MISTAKE-236 guardrail: a precise fixed-map divisor theorem is
not an exact classification of Jacobian counterexamples.

## Reproduction

Exact companion:
[jc_level3_global_exact_pole_residue_20260816.py](../04-computation/jc_level3_global_exact_pole_residue_20260816.py).

```bash
python3 04-computation/jc_level3_global_exact_pole_residue_20260816.py
python3 -O 04-computation/jc_level3_global_exact_pole_residue_20260816.py
python3 04-computation/keller_level_three_norm_slice_probe_20260815.py
```

The new companion pins the inherited slice script/output hashes, the source
`H` artifact and coefficient ledger, the 527-term `P` ledger, the exact
hostile points, and semantic SHA-256
`5c12ecc88284a89667bcbbbbcaa505d221c04a4ccf70a1f9dfdc0456439f570c`.
It writes no result artifact because
this fork was restricted to distinctive files under `04-computation` and
`07-reflections`.
