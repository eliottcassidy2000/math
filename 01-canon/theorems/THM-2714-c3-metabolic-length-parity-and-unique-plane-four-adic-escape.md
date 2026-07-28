---
id: THM-2714
title: "C3 metabolic length parity and unique-plane four-adic escape"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  Let a finite two-primary discriminant module with isometric C3
  action admit a C3-stable metabolizer.  Its standard isotypic sector is a
  finite module over the unramified quadratic DVR O=Z_2[omega], and its
  total O-length is even.  Hence a two-elementary standard sector has even
  THM-2708 Hermitian nullity; if it carries the quartic plane, the nullity
  is at least two.  At nullity one the only possible standard sector is
  O/2^(2r), its metabolizer is 2^r O/2^(2r), and the canonical plane is at
  THM-2695 secondary level one exactly for r=1 and level two for r>=2.
  This is a necessary discriminant-module gate, not geometric realization,
  reflection completion, or a general A4/S4/JC2/DC2 conclusion.
source: a4-resolvent-next-gate-scout-2026-07-28
depends_on:
  - THM-2695-secondary-kummer-bockstein-picard-divisibility-spectrum-and-prime-alignment-boundary
  - THM-2708-c3-hermitian-gain-holonomy-discriminant-gate
  - THM-2711-c3-stable-discriminant-metabolizer-and-mod-four-isotropy-gate
related:
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
  - THM-2700-danielewski-s3-resolvent-standard-plane-exclusion
  - THM-2703-c3-boundary-tree-arm-determinant-standard-plane-gate
script: 04-computation/c3_metabolic_length_parity_thm2714.py
output: 05-knowledge/results/c3_metabolic_length_parity_thm2714.out
script_sha256: 0e824fc5012a1fa0d26a54f8db2d20f85dc404059c9413f23d80e1098ae4b251
output_sha256: 55355e4814b952c2bf8ce374a7df0c6e7df2d946961e3568c50e2be0bcc7a53a
hash_basis: LF-normalized bytes
---

# THM-2714 -- C3 metabolic length parity and unique-plane four-adic escape

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2708 sees how many standard `F2[C3]` planes occur in discriminant
two-torsion.  THM-2711 adds the discriminant pairing and the stable
metabolizer forced by a full-rank unimodular completion.  Combining them
does more than exclude the symplectic `D4` plane: it constrains the entire
two-adic height of the standard sector.

The reason that two and three interact so rigidly is literal arithmetic:

```text
Z_2[C3] = Z_2 direct_product Z_2[omega],
omega^2+omega+1=0,
Z_2[omega]/2 = F4.                                        (1)
```

The second factor is an unramified quadratic discrete valuation ring.  A
metabolizer has square-root order, while every module over this ring has
order a power of four.  Their compatibility forces an even length.

## 1. The integral standard sector

Let `A` be a finite two-primary abelian group equipped with

```text
b:A times A -> Q_2/Z_2                                  (2)
```

a perfect symmetric bilinear pairing.  Let `sigma` act isometrically with
`sigma^3=1`.  Since three is a unit in `Z_2`, the idempotents

```text
e_0=(1+sigma+sigma^2)/3,
e_1=1-e_0                                                 (3)
```

give an exact decomposition

```text
A=A_0 direct_sum A_std.                                  (4)
```

The first summand is the fixed sector.  The second is a finite module over

```text
O=Z_2[T]/(T^2+T+1)=Z_2[omega].                           (5)
```

Because `sigma` is an isometry, `e_0` and `e_1` are self-adjoint.  Hence the
two summands in `(4)` are orthogonal, and the restrictions of `b` to them
are perfect.

Let `H subset A` be a `C3`-stable metabolizer:

```text
H=H^perp.                                                 (6)
```

It decomposes under `(3)` as

```text
H=H_0 direct_sum H_std,
H_std=e_1 H.                                              (7)
```

Orthogonality of `(4)` and `(6)` gives

```text
H_std=H_std^perp inside A_std.                            (8)
```

Thus the standard sector is independently metabolic.  This step is why a
global metabolizer cannot hide a standard deficit in the trivial sector.

## 2. Even `O`-length

The structure theorem over the DVR `O` writes

```text
A_std = direct_sum_(i=1)^s O/(2^(n_i)),
n_i>=1.                                                   (9)
```

Put

```text
N=sum_i n_i=length_O(A_std).                              (10)
```

The residue field of `O` has four elements, so

```text
|A_std|=4^N.                                              (11)
```

Equation `(8)` gives

```text
|H_std|^2=|A_std|,
|H_std|=2^N.                                              (12)
```

But `H_std` is itself a finite `O`-module.  Therefore its order is `4^K`
for some integer `K`.  Comparing with `(12)` proves

```text
N=2K is even,
v_2(|A_std|)=2N is divisible by four.                     (13)
```

This parity is necessary for every `C3`-stable metabolizer, independently
of the detailed discriminant form.  It is not visible in `A[2]` alone when
higher two-power torsion is present.

## 3. The two-elementary doubling theorem

Each summand in `(9)` contributes one copy of the standard plane to
`A[2]`.  Equivalently,

```text
s=dim_F4(A_std[2]).                                       (14)
```

For a boundary lattice in THM-2708, the right side is exactly

```text
s=nullity_F4(B),                                          (15)
```

where `B` is its Hermitian gain-holonomy matrix.

Suppose the standard sector is two-elementary.  Then every `n_i=1`, so

```text
N=s.                                                       (16)
```

Equations `(13)` and `(16)` prove

```text
A_std two-elementary and metabolic
   ==> nullity_F4(B) is even.                              (17)
```

If a quartic character plane is required inside the metabolizer, the
nullity is nonzero.  Therefore

```text
nullity_F4(B)>=2.                                         (18)
```

This strengthens both earlier binary gates in this scope.  Singularity of
`B` only gives nullity at least one; mod-four isotropy only says that one
visible plane has zero restricted pairing.  When the standard discriminant
sector has exponent two, a full stable metabolizer forces a second standard
plane.

In particular, the `D4` root lattice and THM-2708's balanced `3K3` voltage
lift each have a single two-elementary standard sector.  They remain sharp
abstract carrier controls, but neither can be the entire full-rank boundary
lattice in the THM-2711 completion scope.

## 4. Exact classification when only one plane is visible

Now assume

```text
dim_F4(A_std[2])=1.                                       (19)
```

Then `(9)` has exactly one invariant factor:

```text
A_std=O/(2^n).                                            (20)
```

By `(13)`, `n` is even.  Write

```text
n=2r,                    r>=1.                            (21)
```

Every `O`-submodule of the cyclic module `(20)` is an ideal quotient.
The square-root order in `(12)` uniquely forces

```text
H_std=2^r O/(2^(2r)) ~= O/(2^r).                         (22)
```

The only standard plane in `A[2]` is

```text
W=A_std[2]=2^(2r-1) O/(2^(2r)),                           (23)
```

and `(22)` contains it.  Thus a nullity-one completion is possible only by
moving the missing dual plane into higher two-adic depth:

```text
A_std has exponent 2^(2r)>=4,
|A_std|=16^r,
W subset 2 A_std[4].                                      (24)
```

This is the unique-plane **four-adic escape**.  It is the strongest possible
repair of `(17)`: one visible plane is allowed only when its invisible dual
composition factor lies above it in the same cyclic `O`-tower.

## 5. Exact secondary Kummer level of the escape

In the full-rank independent rational-surface application of THM-2711,

```text
Pic(U)_std=H_std ~= O/(2^r).                              (25)
```

The canonical quartic plane is its order-two socle.  THM-2695 identifies
the secondary obstruction with

```text
Pic(U)[2]/2 Pic(U)[4].                                    (26)
```

Equation `(22)` makes the answer exact:

```text
r=1:
  H_std=O/2,
  2 H_std[4]=0,
  W is secondary level 1 (mu_4-nonliftable);

r>=2:
  H_std=O/(2^r),
  W subset 2 H_std[4],
  W is secondary level 2 (mu_4-liftable).                 (27)
```

There is no level-zero branch here because the THM-2711 hypotheses already
place the plane in `Pic(U)[2]`.  Notice also the distinction between ambient
and Picard divisibility.  Equation `(24)` always makes `W` divisible in
`A_std`; equation `(27)` makes it divisible inside the metabolizer only from
`r=2` onward.

## 6. Sharp algebraic controls

### 6.1 One elementary plane fails

The standard part of the `D4` discriminant group is `O/2`.  Its `O`-length
is one, so `(13)` fails.  Equivalently, a metabolizer would have order two,
but every nonzero `O`-submodule has order four.  This rederives the stable
line obstruction in THM-2711 without choosing a basis for its symplectic
form.

### 6.2 Two elementary planes pass

Let

```text
A_std=(O/2) direct_sum (O/2)                              (28)
```

with the orthogonal trace pairing, and put

```text
H_std={(u,u):u in O/2}.                                   (29)
```

Then `H_std=H_std^perp`, it is `C3`-stable, and it contains one diagonal
standard plane.  This is the standard-sector form of THM-2711's doubled
`-2`-triple lattice.  It attains equality in `(18)`.

### 6.3 Every even cyclic height passes algebraically

For each `r>=1`, give

```text
A_r=O/(2^(2r))                                            (30)
```

the perfect pairing

```text
b_r(x,y)=Tr_(O/Z_2)(x conjugate(y))/2^(2r) mod Z_2.       (31)
```

Then

```text
H_r=2^r O/(2^(2r))                                        (32)
```

is `C3`-stable and self-orthogonal.  The cases `r=1,2,3` realize both sides
of `(27)` exactly.  These are finite discriminant-form controls; no claim is
made that every `(30)` is realized by an effective boundary configuration
or a Keller resolvent.

## 7. Exact companion

Run

```text
python 04-computation/c3_metabolic_length_parity_thm2714.py
python -O 04-computation/c3_metabolic_length_parity_thm2714.py
```

and compare both transcripts with

```text
05-knowledge/results/c3_metabolic_length_parity_thm2714.out.
```

The script uses exact integer and `F4` arithmetic with explicit failure
raises.  It verifies:

- the ring law and perfect trace pairings in `(30)`--`(31)` for `r=1,2,3`;
- the exact self-orthogonal subgroups of orders `4,16,64`;
- the `C3` orbit on the three nonzero vectors of every visible plane;
- the secondary nonliftable/liftable split at `r=1` versus `r>=2`;
- the two-elementary diagonal metabolizer `(28)`--`(29)`;
- all `3,002` invariant-factor profiles with at most six factors and
  exponents at most eight, including the parity law and the unique-factor
  exponents `2,4,6,8` in that box.

Normal and optimized transcripts agree after LF normalization and match the
stored output.  The companion contains no Python `assert` nodes.

## 8. Scope and next decisive test

The theorem is a necessary finite-module gate.  It does not prove:

- that a metabolic abstract module is realized by irreducible boundary
  curves on a rational surface;
- that the required plane carries the `C2` reflection completing `C3` to
  the specified `S3` action;
- that a cyclic `O/2^(2r)` tower comes from the graph quartic rather than an
  unrelated cover;
- anything for non-full-rank boundaries where unit squareclasses return;
- any physical LRC homology, endpoint, owner, carry, or phase statement;
- a general `A4`, `S4`, planar Jacobian, or Dixmier conclusion.

For an `S3` action the reflection acts semilinearly by
`omega <-> omega^2`.  Even length survives this extra symmetry, but existence
of a reflection-stable metabolizer is an additional condition.

The cheapest next test is now finite and sharper than another discriminant
comparison:

```text
1. compute the two-primary Smith module of a surviving boundary lattice;
2. split its C3-standard O-primary invariant factors;
3. reject odd total O-length;
4. at nullity one, identify r and compare the forced THM-2695 level (27);
5. only then test semilinear reflection and effective realization.        (33)
```

This is a precise realization of the user's `2`/`3` clue: the order-three
symmetry packages binary classes into `F4`, and unimodularity forces those
packages to occur with even total two-adic depth.

QED (candidate; independent hostile audit pending).
