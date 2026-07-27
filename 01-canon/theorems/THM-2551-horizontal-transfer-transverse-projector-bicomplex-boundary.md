---
id: THM-2551
title: "Integral tensor intersection: horizontal C91 transfer preserves the transverse all-unit kernel"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On the abstract product carrier of THM-2544's relation residues and
  THM-2548's C91 mapping torus, the seven-step transfer commutes with the
  unrestricted/all-unit projectors and preserves their kernels exactly.
  All ranks, saturated images, sharp positive controls, and the Hall hostile
  are integral.  SCOPE: canon has not constructed a lawful physical current
  identifying the relation target, deck root, and semantic arrival roots;
  no row is removed and LRC(14) remains OPEN.
source: codex-2026-07-27-holotopy
depends_on:
  - THM-2544-all-unit-target-projector-kernel-and-lawful-image-intersection-obstruction
  - THM-2548-seven-step-c91-transfer-and-full-norm-separation
  - THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2543-augmentation-norm-relative-phase-local-system-dichotomy
script: 04-computation/lrc14_transfer_projector_tensor_intersection_thm2551.py
output: 05-knowledge/results/lrc14_transfer_projector_tensor_intersection_thm2551.out
script_sha256: 9eea54d9ed8f374589b263198aa89fc586809a50fa7f4938c968a4e18662c646
output_sha256: 10ffad586382d3a8a0fb7a58607d83d360ec5236414586b789b985290a364e61
hash_basis: LF-normalized bytes
---

# THM-2551 -- horizontal transfer preserves the transverse projector wall

This theorem asks the cheapest possible compatibility question between two
recent exact carriers:

```text
THM-2544: relation-residue current --J--> all-unit target current;
THM-2548: C91 chart table         --D--> seven-step transferred table.
```

On their honest product they do not help one another.  They form a commuting
square, and the all-unit kernel is an invariant horizontal subsystem.  This
is the precise algebraic content of the remaining "mixed 2-cell" debt.

The word *holotopy* below is only mnemonic: `D` is a transfer, not a
square-zero differential.  Every assertion is an ordinary statement about
split free abelian groups.

## 1. Integral splittings on the two axes

Use THM-2544's notation

```text
Omega=K_91,                 n=|Omega|=91^8,
R=Z^Omega,
U,J:R -> Z^G,               |G|=169.                 (1)
```

Assume the septimal support has size at least two.  For every `q in G`,
THM-2544 supplies a unit point `u_q`, a nonunit point `v_q` in the same
target fibre, and

```text
h_q=delta_(u_q)-delta_(v_q),
(U,J)delta_(v_q)=(delta_q,0),
(U,J)h_q=(0,delta_q).                                  (2)
```

Put

```text
E_U=<delta_(v_q):q in G>,       E_J=<h_q:q in G>,
K=ker(U,J).                                               (3)
```

Subtracting the two displayed sections from any vector proves the integral
direct sum

```text
R=K direct-sum E_U direct-sum E_J,                        (4)

rank K=n-338,             rank E_U=rank E_J=169,
(U,J)|(E_U direct-sum E_J)=(id,0) direct-sum (0,id).       (5)
```

In particular `R=ker J direct-sum E_J` and
`ker U=K direct-sum E_J`.  No field or dimension argument is being used.

On the horizontal axis let

```text
T=Z^(C_7 x C_13),
D=D_a=1+S_a+...+S_a^6,
A_clk=ker D,                 L_eq=im D.                    (6)
```

THM-2548 gives the split exact sequence

```text
0 -> A_clk -> T --D--> L_eq -> 0,
rank A_clk=6,              rank L_eq=85,                  (7)
```

with primitive image and Smith form `1^85 direct-sum 0^6`.

## 2. The commuting square and exact intersection

On

```text
W=R tensor T                                                   (8)
```

define

```text
Dbar=I_R tensor D,
Jbar=J tensor I_T,              Ubar=U tensor I_T.          (9)
```

Write `D_G=I_(Z^G) tensor D`.  Then

```text
Jbar Dbar = D_G Jbar,            Ubar Dbar = D_G Ubar.      (10)
```

The split `R=ker J direct-sum E_J` makes the intersection literal:

```text
im Dbar intersection ker Jbar
  =ker J tensor L_eq
  =Dbar(ker Jbar).                                         (11)
```

Thus horizontal transfer preserves the entire transverse all-unit kernel;
it neither leaks out of it nor creates an all-unit current from a kernel
input.  Its rank is

```text
85(n-169)=399714648472864920.                              (12)
```

The composite has

```text
ker(Jbar Dbar)
 =ker J tensor T  direct-sum  E_J tensor A_clk,

rank ker(Jbar Dbar)=91n-169*85
                    =427929800129774046,                   (13)

im(Jbar Dbar)=Z^G tensor L_eq,
rank im(Jbar Dbar)=169*85=14365.                           (14)
```

The image in `Z^G tensor T` is saturated.  Equivalently all nonzero Smith
invariants are `1`, and the target cokernel is free of rank

```text
169*6=1014.                                                (15)
```

For the designated nonzero-target projector
`J_*=pr_(G minus {0}) J`, replace `169` by `168`:

```text
rank(im Dbar intersection ker J_*bar)
  =85(n-168)=399714648472865005,
rank(J_*bar Dbar)=14280,
rank ker(J_*bar Dbar)=427929800129774131,
rank coker(J_*bar Dbar)=1008.                              (16)
```

## 3. The unrestricted and all-unit outputs remain independent

Apply the full split (4).  The joint composite satisfies

```text
ker(((U,J) tensor I)Dbar)
 =K tensor T
   direct-sum E_U tensor A_clk
   direct-sum E_J tensor A_clk,                            (17)

im(((U,J) tensor I)Dbar)
 =(Z^G direct-sum Z^G) tensor L_eq.                        (18)
```

Hence

```text
rank image=338*85=28730,
rank kernel=91n-338*85=427929800129759681,
target cokernel free of rank 338*6=2028.                   (19)
```

More sharply,

```text
Jbar(ker Ubar intersection im Dbar)=Z^G tensor L_eq.       (20)
```

So even after transfer, `U` and `J` remain independently prescribable on
every transferred profile.  The partial transfer does not couple them.

## 4. Exact survival criterion and sharp controls

Equation (10) and THM-2548 give the iff criterion

```text
Jbar Dbar x=0
  iff Jbar x belongs to Z^G tensor A_clk.                  (21)
```

After tensoring with any characteristic-zero splitting field,

```text
T=T_const direct-sum T_root direct-sum A_clk,
dim=(1,84,6),                                                (22)
```

where `T_root` is the sum of all root-charged C91 characters.  Thus every
already all-unit root-charged signal survives in either order; an all-unit
kernel input remains invisible in either order.  The criterion is sharp:

- `p=sum_q delta_(u_q)` has `Jp=1_G`; tensoring it with any nonzero
  `a_clk in A_clk` gives a `J`-visible vector killed by `Dbar`;
- `c_host=sum_q delta_(v_q)` is nonnegative, has `Uc_host=1_G` and
  `Jc_host=0`; with `e=delta_(0,0)`,

```text
x=Dbar(c_host tensor e)=c_host tensor De                   (23)
```

  is nonnegative, has full unrestricted target support and every one of the
  `84` root-charged C91 modes, but still has `Jbar x=0`;
- for arbitrary real `0<=b(q)<=a(q)`, THM-2544's positive current
  `c_(a,b) tensor De` realizes the transferred outputs `(a tensor De,
  b tensor De)` exactly.

The full 91-step norm is strictly worse: it retains only the constant C91
character, exactly as in THM-2548.

## 5. Simultaneous projector success still does not imply arrival

The remaining semantic roots are a third factor.  Let

```text
A=delta_(0,0)+delta_(1,1),
B=delta_(0,1)+delta_(1,0)                                  (24)
```

in `H=Z^(F_13 x F_13)`.  They have identical head and later marginals, but
diagonal masses `2` and `0`.  Let `epsilon:H->Z` be total augmentation.
Tensor either with `p tensor De`, where `Jp=Up=1_G`.  After applying
`epsilon` to the semantic factor, the two nonnegative product currents have
identical unrestricted and all-unit target outputs and identical complete
root-charged transfer profiles; before augmentation their head and later
one-point marginals also agree.  Their THM-2545 Hall diagonal differs from
full to zero.  The un-marginalized `H`-valued outputs are, of course,
different.

Therefore even simultaneous `U`, `J`, and `D` success cannot recover the
semantic coupling from product observables.  This is the three-factor form
of THM-2545's aligned/swap hostile.

## 6. Holotopy interpretation and exact boundary

The horizontal transfer is a partial contraction only along the C91 axis.
The submodule

```text
ker J tensor T                                               (25)
```

is transverse and `Dbar`-invariant, with horizontal image exactly
`ker J tensor L_eq`.  A rescue must therefore impose a lawful non-product
coupling or restriction not visible in this ambient tensor diagram.  A
root- or semantic-dependent projector with a nonzero mixed commutator is one
operator realization of such a 2-cell, but is not asserted to be the only
possible realization.

This conclusion is type-sensitive.  THM-2544's relation target is
`q in F_13^2`; THM-2548's deck root is `r in F_13`; THM-2545's head and
later semantic roots are `h,b in F_13`.  No canonical theorem identifies
them.  Ordinary residue translation may also create zero coordinates and
does not preserve the all-unit mask.  The tensor product above is an exact
ambient boundary, not a constructed physical LRC current.

No row is excluded.  LRC(14) remains open.

## 7. Exact companion

The dependency-free companion checks all rank arithmetic, the integral
split blocks, all `91` C91 transfer characters over `F_547`, the positive
kernel and clock-augmentation controls, and the aligned/swap Hall hostile.
Run

```bash
python3 04-computation/lrc14_transfer_projector_tensor_intersection_thm2551.py
python3 -O 04-computation/lrc14_transfer_projector_tensor_intersection_thm2551.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_transfer_projector_tensor_intersection_thm2551.out
```

after LF normalization.  Every executable check raises explicitly under
optimized Python.

## 8. Independent hostile audit

An independent audit reconstructed the integral decompositions from the
literal `v_q,u_q,h_q` sections and a chosen lift of `L_eq`; no flatness or
field-dimension inference is used.  It rederived every intersection, kernel,
image, rank, and free-cokernel formula, and independently computed the
`91 x 91` transfer Smith form `1^85 direct-sum 0^6`.  It also checked the
`84/6` Fourier split, nonnegative kernel hostile, clock-augmentation sharp
boundary, and the semantic aligned/swap observation map.  The audit required
the explicit semantic augmentation in Section 5 and broadened Section 6 from
"mixed commutator" to any lawful non-product coupling.  Normal, optimized,
and stored transcripts agree exactly after LF normalization; the hashes in
frontmatter are independently reproduced.
