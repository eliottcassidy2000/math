---
id: THM-3877
title: "Sign-kernel transfer through a torsion-free quadratic resolvent"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For a finite
  separable field extension with sole reduced field-discriminant divisor of
  multiplicity one, the normal-closure inertia is a transposition.  The
  discriminant quadratic base change absorbs that inertia.  If the normalized
  quadratic resolvent has p-divisible units and no p-torsion in its class
  group, the remaining sign kernel cannot have a cyclic quotient of order p.
  If these conditions hold for every prime, the sign kernel must be perfect.
  Consequently the THM-3874 three-cusp resolvent excludes every degree-three
  and degree-four field with that quintic as sole simple branch divisor.  The
  theorem is a transfer obstruction, not a construction and not JC(2).
source: jc_quartic_c3_construct / post-THM-3874 sign-kernel lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The audit rederived the
  permutation-representation Artin conductor of a non-Galois field, used
  faithfulness of the normal-closure action to turn conductor one into a
  literal order-two transposition, and checked inertia intersection after the
  sign base change.  It independently verified the Krull divisor/Kummer
  argument for every prime, the perfect-kernel conclusion, and the complete
  degree-three/four quantifiers, including the C4 and A5 boundaries.  The
  proof is field- and presentation-independent.  The exact companion verifies the tame permutation-conductor
  test, all degree-three and degree-four sign kernels, their cyclic
  abelianizations, the absence of transpositions in the regular C4 action,
  and the perfect A5 equality boundary.  Normal and optimized runs byte-match
  the frozen transcript.
depends_on:
  - THM-3874-three-cusp-quadratic-k3-affine-class-group
related:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3869-three-cusp-square-residual-cardano-line-ramification
script: 04-computation/jc2_sign_kernel_cyclic_quotient_transfer_thm3877.py
output: 05-knowledge/results/jc2_sign_kernel_cyclic_quotient_transfer_thm3877.out
script_sha256: 2b0f8b43d6d415e382404246d9b00ff604f8ae08928ce99e6ed07f83cb19d4dc
output_sha256: b8dd1bf05cc59218248b9b4a05f2312d1f1b8d59ac1c8c19e5e5a9655e6378ba
semantic_sha256: 5f7d99f84623a015d9c0a537e84ac24148b14081e1dd08e59f9ae97287bbac40
hash_basis: raw LF bytes
---

# THM-3877 -- the sign cover transfers simple branch into an unramified kernel

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero.  The core statement is
valid over any normal affine `k`-domain; no polynomial presentation of the
cover is required.

## 1. The reusable transfer theorem

Let `R` be a normal finitely generated `k`-domain, let `K=Frac(R)`, and let
`L/K` be a finite separable extension of degree `d`.  Write `M/K` for its
normal closure and

```text
G=Gal(M/K) subset S_d                                           (1)
```

for the faithful transitive action on the `d` embeddings of `L`.  Assume that
the field-discriminant divisor of the normalization of `R` in `L` is a
nonempty reduced divisor `D`, with coefficient exactly one on every
irreducible component and no other codimension-one branch.

Put

```text
H=G intersection A_d,        E=M^H,                            (2)
```

and let `A` be the normalization of `R` in the discriminant quadratic field
`E`.  Then:

1. `M/E` is unramified at every height-one valuation of `A`.
2. Fix a prime `p`.  If

   ```text
   A^*/(A^*)^p=0,                 Cl(A)[p]=0,                  (3)
   ```

   then `H` has no quotient isomorphic to `C_p`.
3. If `(3)` holds for every prime, then `H` is perfect.

The statement concerns codimension-one ramification.  It does not assert that
the corresponding normalization is etale at singular codimension-two points.

## 2. Simple field discriminant means one transposition

Fix a height-one prime `q` of `R`, a prime of `M` above it, and its inertia
group `I`.  Residue characteristic zero makes the ramification tame.  For the
permutation representation `(1)`, the field-discriminant exponent is its tame
Artin conductor

```text
v_q(Disc(L/K))=d-#orbits(I on {1,...,d}).                     (4)
```

At a component of `D`, the left side is one.  A nontrivial cyclic permutation
group with `d-1` orbits has exactly one orbit of length two and fixes every
other letter.  Hence

```text
I=<tau>,                     tau a transposition.             (5)
```

Outside `D` the exponent is zero, so `I=1`.  In particular the sign character
of `G` is nontrivial, and `E/K` in `(2)` is genuinely quadratic.

For a Galois subextension, inertia is obtained by intersection with the
corresponding subgroup.  Thus the inertia of `M/E` is

```text
I intersection H=<tau> intersection A_d=1.                   (6)
```

Equations `(5)-(6)` at `D`, together with trivial inertia away from `D`, prove
part 1.

## 3. Units and class torsion detect every cyclic quotient

Suppose that a surjection

```text
chi:H -> C_p                                                   (7)
```

exists.  Its kernel is normal in `H`, so the fixed field
`N=M^(ker chi)` is a connected cyclic degree-`p` extension of `E`.  By part 1
it is unramified at every height-one valuation of `A`.

Because `mu_p subset k`, Kummer theory writes

```text
N=E(gamma^(1/p))                                               (8)
```

for some `gamma in E^*`.  Codimension-one unramifiedness implies

```text
div_A(gamma)=p Z                                               (9)
```

for a Weil divisor `Z`.  Therefore `[Z] in Cl(A)[p]`.  Under `(3)`, write
`Z=div_A(beta)`.  The element

```text
u=gamma/beta^p                                                (10)
```

has zero valuation at every height-one prime.  Since `A` is a normal
noetherian domain, it is the intersection of those valuation rings in `E`;
applying the same statement to `u^(-1)` shows `u in A^*`.  The first condition
in `(3)` gives `u=v^p`.  Hence `gamma=(beta v)^p`, so `(8)` is trivial, a
contradiction.  This proves part 2.

If `(3)` holds for every prime and `H` were not perfect, the nontrivial finite
abelian group `H^ab` would have a quotient `C_p` for some prime `p`.  Part 2
excludes it, proving part 3.

## 4. Degrees three and four

The transfer is already decisive in low degree.

For `d=3`, transitivity and the nontrivial sign character force `G=S_3`, so

```text
H=A_3=C_3.                                                    (11)
```

For `d=4`, transitivity gives `4 | |G|`, while the sign character gives
`[G:H]=2`; hence `H` is a nontrivial even-order subgroup of `A_4`.  The group
`A_4` is solvable, so `H` is solvable.  A nontrivial finite solvable group
cannot be perfect: otherwise its derived series would never reach `1`.
Consequently `H^ab` has a quotient `C_2` or `C_3`.  This avoids any dependence
on a classification of transitive quartic groups.  As an exact control, that
classification yields

```text
S4: H=A4, H/V4=C3;       D4: H=V4, H -> C2;                 (12)
```

while the regular `C4` action contains no transposition and cannot satisfy
the simple-discriminant hypothesis.  The even groups `A4,V4` have square
discriminant and are excluded before `(2)`.

The boundary is sharp for this mechanism.  In degree two, `H=1` and the
quadratic branch cover itself is allowed.  For full symmetric degree five,

```text
G=S5,                         H=A5=[A5,A5],                  (13)
```

so the cyclic-quotient transfer is silent.  It does not say that such a cover
exists.

## 5. The THM-3874 specialization

For the exact three-cusp quintic of THM-3854, THM-3874 proves that the normal
quadratic resolvent `Q=Spec(A)` satisfies

```text
A^*=k^*,                         Cl(A)=Z.                    (14)
```

Since `k` is algebraically closed, `(3)` holds for every prime.  Applying
Sections 1-4 proves that no separable field extension of degree three or four
over `k(x,y)` can have that quintic as its sole simple field-discriminant
divisor.  This is field-level: it includes, but is not limited to, normal
finite-flat cubic and quartic orders.  In degree four it is precisely the
`V4`/resolvent-`C3` obstruction suggested by `S4/V4=S3`, with the sign `C2`
layer separated first.

THM-3874 additionally turns a depressed-cubic identity
`P^3-Q^2=Delta H_0^2` into an unavoidable totally ramified factor of `H_0`.
That exact Cardano payment is not reproved here.

## 6. Scope and reproduction

The theorem does not obstruct covers whose sign kernel is perfect, branch
divisors with higher discriminant multiplicity, extra branch components,
wild ramification, or cyclic layers supported only after changing the
quadratic resolvent.  It constructs neither a Keller map nor a Jacobian
counterexample; `JC(2)` remains open.

Run

```bash
python3 04-computation/jc2_sign_kernel_cyclic_quotient_transfer_thm3877.py
python3 -O 04-computation/jc2_sign_kernel_cyclic_quotient_transfer_thm3877.py
```

and compare both streams byte-for-byte with
`05-knowledge/results/jc2_sign_kernel_cyclic_quotient_transfer_thm3877.out`.
