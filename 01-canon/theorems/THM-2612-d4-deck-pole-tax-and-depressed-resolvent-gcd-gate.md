---
id: THM-2612
title: "D4 deck pole tax and depressed-resolvent gcd gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For a hypothetical complex planar Keller map of generic degree four and
  geometric monodromy D4, the quartic root field has a unique nontrivial
  automorphism. It is a regular free deck involution over the finite-etale
  complement of the Jelonek set, but its induced birational involution of
  A2 must have a nonempty pole divisor supported over the Jelonek boundary.
  Otherwise it would extend to a polynomial involution of A2; Smith gives a
  fixed point, while global etaleness and F o tau=F make any such fixed point
  force tau=1. For a depressed primitive quartic, the involution and its
  exact denominator norm are explicit. Separately, at every tame odd-residue
  divisor, positive quartic-to-resolvent order-index tax is possible only
  when that divisor divides both q and p^2-4r. The gcd gate is necessary, not
  sufficient. No boundary component, pole order, D4 exclusion, degree-22
  mixed-support closure, split/even descent, 2-adic order raising, JC(2), or
  DC(2) follows.
source: klein-2026-07-28-d4-jelonek-boundary
depends_on:
  - THM-2465-g1-exclusion-package-for-degree-four-twojet-keller
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
related:
  - MISTAKE-297
  - HYP-9027-twojet-disc-jelonek-odd-exponent-law
script: 04-computation/jacobian_d4_deck_pole_tax_thm2612.py
output: 05-knowledge/results/jacobian_d4_deck_pole_tax_thm2612.out
script_sha256: ad6f6f02a119b9cb049274da99fdd47cf7f2fa824068292d7b45618753faf8cf
output_sha256: 62370d58cc8a47df041609359e2befa72e0fde6c5cf916ca771a28264cd90a88
hash_basis: working-tree bytes (LF)
---

# THM-2612 -- the D4 involution must pay at the omitted boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

MISTAKE-297 repaired the old degree-four Smith argument: a deck action on the
finite-etale open source cannot be applied to the contractible completed
source unless the action first extends across the omitted divisor.  In the
`D4` lane that failure is not merely a warning.  The unique root-field
involution cannot extend, so it must acquire a genuine divisorial pole over
the Jelonek boundary.

This gives a lawful replacement for the discarded argument.  It does **not**
exclude the `D4` lane: a pole is precisely what allows the open involution to
avoid the fixed-point contradiction.

## 1. Exact statement

Let

```text
F=(P,Q): X=A^2_C -> Y=A^2_C                              (1)
```

be a Keller map with `det JF in C*`, generic degree four, and geometric
monodromy `D4` in its transitive action on four sheets.  Put

```text
K=C(Y),                    L=C(X),
A_F=the nonproper-value set of F,
V=Y\A_F,                   U=F^(-1)(V).                 (2)
```

The inherited purity statement in THM-2465 gives a connected finite-etale
cover

```text
F|_U: U -> V.                                               (3)
```

Then:

1. `Aut_K(L)=C2`.  Its nonidentity element `tau` extends uniquely to a deck
   involution of (3).
2. As a birational self-map of `X`, `tau` has a nonempty pole divisor.  Every
   pole prime lies in

   ```text
   X\U=F^(-1)(A_F).                                         (4)
   ```

   Equivalently, for affine coordinates `x,y` on `X`, some prime divisor
   `E subseteq F^(-1)(A_F)` has

   ```text
   min(v_E(tau^*x),v_E(tau^*y))<0.                          (5)
   ```
3. If a primitive element `w` has depressed minimal polynomial

   ```text
   f(T)=T^4+pT^2+qT+r in K[T],                              (6)
   ```

   and `z in K` is the unique rational root of the squared-pair resolvent

   ```text
   S(U)=U^3+2pU^2+(p^2-4r)U-q^2,                            (7)
   ```

   then, in `L`,

   ```text
   D=2w^2+p+z,
   s=(q-2zw)/D,
   s^2=z,
   tau(w)=-w-s,                                             (8)

   Norm_(L/K)(D)=S'(z)^2 != 0.                              (9)
   ```

   Formula (8) includes the `z=q=0` boundary; for `T^4-2`, it is simply
   `tau(w)=-w`.
4. Let `R` be a discrete valuation ring with uniformizer `pi`, residue
   characteristic different from two, and fraction field of characteristic
   zero.  Suppose (6) is monic and `R`-integral, the generic quartic is
   separable, and the local inertia is tame.  Let `i_4,i_3` be the index
   lengths of the quartic polynomial order and the squared-pair-resolvent
   order in their normalizations.  If

   ```text
   i_3-i_4>0,                                               (10)
   ```

   then in fact `i_3-i_4=1` and

   ```text
   pi | q,                    pi | (p^2-4r).                 (11)
   ```

   Thus, on any factorial odd-residue coefficient model, the positive-tax
   support is contained in

   ```text
   V(gcd(q,p^2-4r)).                                        (12)
   ```

The pole conclusion and the gcd conclusion are distinct.  The theorem does
not say that the pole divisor equals (12), or that every divisor in (12)
carries either a pole or positive index tax.

## 2. The field automorphism is canonical

Write

```text
D4=<rho,sigma | rho^4=sigma^2=1, sigma rho sigma=rho^(-1)>   (13)
```

in the vertex action on a square, and let `H` be one vertex stabilizer.
Then

```text
|H|=2,          |N_D4(H)|=4,
N_D4(H)/H=C2.                                             (14)
```

The usual normalizer formula for a separable root field gives

```text
Aut_K(L)=N_D4(H)/H=C2.                                    (15)
```

The nonidentity coset sends a vertex to its opposite vertex.  Moreover the
subgroup interval

```text
H < J < D4                                                (16)
```

has exactly one member, namely `N_D4(H)`.  This recovers the corrected
THM-2465 statement that a `D4` quartic root field has one proper
intermediate field, while avoiding the false maximal-stabilizer claim from
MISTAKE-297.

Since `U` is the normalization of `V` in `L`, the field automorphism extends
uniquely to a deck involution of (3).  It is nontrivial and hence acts freely
on the connected finite-etale cover.

## 3. Why a boundary pole is forced

The deck involution is regular on `U`, so every divisorial pole of its two
coordinate functions is supported in `X\U`.

Suppose, for contradiction, that it has no divisorial pole anywhere on
`X=A^2`.  The affine plane is normal.  A rational function regular at every
codimension-one point of a normal affine variety is globally regular, so
both `tau^*x` and `tau^*y` are polynomials.  Since `tau^2=1` in `L`, this
extends `tau` to a polynomial involution of all of `A^2`.

Now the corrected Smith step is legal.  The complex affine plane is
mod-`2` acyclic, so Smith fixed-point theory gives a point

```text
a in A^2(C),                    tau(a)=a.                      (17)
```

The deck identity on the dense open set extends as a polynomial identity:

```text
F o tau=F                                             on A^2. (18)
```

Because `F` is Keller, it is a local analytic biholomorphism at `a`.  Shrink
to a neighborhood `W` on which `F` is injective, and then use continuity of
`tau` to choose a smaller neighborhood
`W_0 subset W intersection tau^(-1)(W)`.  Equations (17)--(18) give
`tau=1` on `W_0`, hence everywhere by polynomial identity.  This contradicts
(15).  Therefore the pole divisor is nonempty, proving (4)--(5).

This is exactly the missing hypothesis in the old Smith argument.  Smith
does not exclude the open `D4` action; it proves that extending that action
without a pole is impossible.

## 4. The involution's functional form

For a depressed quartic, the three squared pair sums are the roots of (7).
In the `D4` lane the pairing action has one fixed point and one two-cycle, so
`S` has one rational root `z` and one irreducible quadratic factor.

Put `D,s` as in (8).  Exact reduction modulo `f(w)=S(z)=0` gives

```text
(q-2zw)^2-z(2w^2+p+z)^2=0,                                (19)

f(-w-s)=0,                   s(-w-s)=s,                    (20)

tau^2(w)=w.                                                   (21)
```

The denominator is generically genuine and has an invariant norm:

```text
Res_T(f(T),2T^2+p+z)-S'(z)^2
   =-8(p+z)S(z).                                             (22)
```

Since the quartic is separable, so is the squared-pair resolvent, and hence
`S'(z)!=0`; thus (9) also proves `D!=0` in `L`.  The resulting root cannot
equal `w`: otherwise `s=-2w`, so `4w^2=s^2=z in K`, contradicting
`[K(w):K]=4`.  It is therefore the nonidentity element of (15), not the
identity embedding.  Equations (19)--(22) prove (8)--(9).

The norm formula is a useful functional description, but it is not yet a
component-by-component pole formula on `A^2`: the primitive coordinate,
the depressed coefficients, and `z` may themselves carry boundary
valuations.  Any future exact pole-order theorem must retain those sidecars.

## 5. The tame gcd gate

THM-2598 proves, at a tame divisor, that the raw quartic and resolvent
polynomial discriminants agree and that

```text
i_3-i_4=(d_4-d_3)/2,                                       (23)
```

where `d_4,d_3` are the discriminant exponents of the normalized quartic
and pairing covers.  The complete inertia table is

| inertia on four roots | `d_4` | action on pairings | `d_3` | tax |
|---|---:|---|---:|---:|
| identity | 0 | identity | 0 | 0 |
| transposition | 1 | transposition | 1 | 0 |
| double transposition | 2 | identity | 0 | 1 |
| three-cycle | 2 | three-cycle | 2 | 0 |
| four-cycle | 3 | transposition | 1 | 1 |

Thus positive tax means double-transposition or four-cycle inertia.

Pass to a strict henselization and reduce the four integral roots modulo
`pi`.  Roots in one inertia orbit have the same residue.

- For a four-cycle all four roots reduce to one value `a`.  Depression gives
  `4a=0`, hence `a=0`, and the residual quartic is `T^4`.
- For a double transposition the residual multiset is `(a,a,b,b)`.
  Depression gives `2(a+b)=0`, hence `b=-a`, and the residual quartic is

  ```text
  (T-a)^2(T+a)^2=(T^2-a^2)^2.                              (24)
  ```

Both residual polynomials have

```text
q=0,                       p^2-4r=0.                        (25)
```

This proves (11)--(12).

## 6. The gcd gate is not sufficient

Over `C[[t]]`, consider the tame local quartic algebra

```text
f_t(T)=((T-1)^2-t)((T+1)^2-t^2).                            (26)
```

It is depressed, with

```text
p=-t^2-t-2,
q=2t(t-1),
r=(t-1)^2(t+1),
p^2-4r=t(t^3-2t^2+9t+8).                                  (27)
```

Thus `t` divides both entries of the gcd gate.  Its squared-pair resolvent is

```text
(U-4)(U^2-2U(t^2+t)+t^4-2t^3+t^2).                         (28)
```

Both raw polynomial discriminants have valuation three.  The only field
ramification is the first quadratic's single transposition, so both
normalized discriminant exponents are one and

```text
i_4=i_3=1,                    i_3-i_4=0.                     (29)
```

This is a local-algebra hostile to any coefficient-only converse.  It is not
a `D4` Keller witness and proves no sufficiency or nonsufficiency theorem
after adding extra global hypotheses.

## 7. Exact verification

Run

```bash
python3 04-computation/jacobian_d4_deck_pole_tax_thm2612.py
python3 -O 04-computation/jacobian_d4_deck_pole_tax_thm2612.py
```

The companion independently checks:

- the `D4` point-stabilizer normalizer and subgroup interval;
- all three rational identities in (19)--(22);
- the complete five-row inertia and index-tax table;
- the two positive-tax residual quartics;
- the exact local hostile (26)--(29); and
- the `T^4-2` zero-resolvent-root boundary.

All truth-bearing checks use explicit `require` calls and remain active under
optimized Python.  Normal, optimized, and stored transcripts byte-match.

## 8. Scope and remaining frontier

The gain is an exact boundary invoice:

```text
D4 root-field symmetry
   -> unique involution on the finite-etale open
   -> unavoidable pole somewhere over A_F,                  (30)

positive tame quartic/resolvent index tax
   -> divisor lies in V(gcd(q,p^2-4r)).                      (31)
```

The information still missing is equally exact:

- which Jelonek component carries the pole;
- the pole order and whether it meets the positive-tax support;
- an integral affine formula for the depressed coefficients and rational
  resolvent root on the physical source;
- a contradiction with Keller geometry.

The tame proof deliberately excludes residue characteristic two.  It does
not address integral `2`-adic order raising.  It also does not address the
split/even short-edge descent or any degree-twenty-two mixed-support
stratum.  No `D4` witness is constructed or excluded, and `JC(2)` and
`DC(2)` remain open.

**QED.**

The independent hostile audit rederived the normalization/deck-extension
step over the finite-etale Jelonek complement; checked freeness and the
support of every pole prime; and verified the normal-affine
codimension-one intersection argument before applying Smith theory.  It
made the local-injectivity shrink in Section 3 explicit and checked that a
fixed point forces the polynomial involution, rather than merely its germ,
to be the identity.

The same audit independently enumerated the `D_4` point-stabilizer
normalizer and the complete strict subgroup interval; derived the depressed
factorization signs in (8), including `q=z=0`; used (9) to certify the
denominator is nonzero; and checked that the formula is the nonidentity
root-field automorphism.  It rederived every tame inertia row, the
double-transposition and four-cycle reductions, the direction of the index
tax, and the zero-tax gcd hostile.  Finally it replayed ordinary and
optimized companions against the stored transcript and recovered the
declared LF hashes
`ad6f6f02a119b9cb049274da99fdd47cf7f2fa824068292d7b45618753faf8cf`
and
`62370d58cc8a47df041609359e2befa72e0fde6c5cf916ca771a28264cd90a88`.
No pole component, pole order, gcd converse, `D_4` exclusion, mixed-support
closure, wild residue-two claim, or Jacobian conclusion entered the audit.
