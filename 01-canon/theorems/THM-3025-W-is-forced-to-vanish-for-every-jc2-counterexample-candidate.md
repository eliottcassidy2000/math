---
id: THM-3025
title: "JC(2): the subleading correction W vanishes outside the one-root leading-form locus"
status: >
  PROVED + VERIFIED-EXACT / REPAIRED UNDER MISTAKE-422 / AWAITING
  INDEPENDENT HOSTILE AUDIT.  In the setup of THM-3016, write
  J0=Jac(P,Q) in k*, j=Jac(H,Q_(m-1)), and
  W=P_(n-1)-kappa H^(a-b)Q_(m-1).  THM-3016 proves
  j Jac(W,H)=0; its multiplier is j, not J0.  If j!=0, the common-form and
  coprime-degree argument forces W=0 whenever H has at least two distinct
  roots.  If j=0, the same argument first forces Q_(m-1)=0 and then, through
  the subleading Keller equation, P_(n-1)=0.  Thus K>=2 implies W=0 in both
  branches.  The earlier conclusion survives, but its one-line division did
  not; see MISTAKE-422.
source: death-star-2026-07-31-coinC2
depends_on:
  - THM-3016
related:
  - HYP-9070
  - THM-3003
external:
  - "Abhyankar--Moh: embeddings of the line in the plane (one place at infinity)."
  - "Jung--van der Kulk: structure of the plane automorphism group."
script: 04-computation/jc2_W_forced_thm3025.py
output: 05-knowledge/results/jc2_W_forced_thm3025.out
---

# THM-3025 -- `W=0` outside the one-root locus

## 1. Typed setup and the corrected rigidity identity

Let `k` have characteristic zero and let `(P,Q)` be a planar Jacobian pair.
Use distinct notation for the two Jacobians:

```text
J0 := Jac(P,Q) in k*,                         (the full Keller constant)
j  := Jac(H,Q_(m-1)).                         (a subleading form)
```

As in THM-3016, suppose

```text
P_n=c H^a,  Q_m=c' H^b,  deg H=g=gcd(n,m),
n=ga,       m=gb,         gcd(a,b)=1,
kappa=ca/(c'b),            a>=b,
W=P_(n-1)-kappa H^(a-b)Q_(m-1).
```

The degree-`n+m-3` Keller equation is

```text
c a H^(a-1) j = -c' b H^(b-1) Jac(P_(n-1),H).            (L1)
```

THM-3016's Pluecker calculation proves

```text
j Jac(W,H)=0.                                             (R)
```

The multiplier in `(R)` is **not** `J0`.  The previous version divided by it
as though it were the nonzero Keller constant; MISTAKE-422 records that
error.  The proof must split on `j`.

## 2. Common-form degree lemma

If nonzero homogeneous binary forms `F,G` satisfy `Jac(F,G)=0`, then they are
powers of a common homogeneous form `R`.  In particular `deg R` divides both
`deg F` and `deg G`.

Apply this with `G=H`, `deg H=g`.  Since `g|n` and `g|m`,

```text
gcd(g,n-1)=gcd(g,m-1)=1.                                 (C)
```

Consequently, if either a nonzero form of degree `n-1` or a nonzero form of
degree `m-1` has zero Jacobian with `H`, the common form has degree one and
`H` is a power of one linear form.  Equivalently, `H` has one distinct root
direction (`K=1`).

## 3. The two branches

### Branch I: `j!=0`

Equation `(R)` gives `Jac(W,H)=0`.  Here `deg W=n-1`; if `W!=0`, the lemma
and `(C)` force `K=1`.  Therefore

```text
K>=2 and j!=0  =>  W=0.                                  (B1)
```

### Branch II: `j=0`

If `Q_(m-1)!=0`, then `Jac(H,Q_(m-1))=0`; the lemma and `(C)` again force
`K=1`.  Hence `K>=2` forces `Q_(m-1)=0`.  Equation `(L1)` then gives

```text
Jac(P_(n-1),H)=0.
```

If `P_(n-1)!=0`, the same lemma, now in degree `n-1`, forces `K=1`.
Therefore on the `K>=2` locus both subleading forms vanish, and so does `W`:

```text
K>=2 and j=0  =>  Q_(m-1)=P_(n-1)=W=0.                  (B2)
```

Combining `(B1)` and `(B2)` proves

```text
K>=2  =>  W=0.                                           (T)
```

The mirrored proof for `b>=a` exchanges `P` and `Q`.

## 4. Sharp boundary and tower consequence

The restriction is sharp at the level of the homogeneous equation:

```text
H=(x+y)^2, deg W=3: W=lambda(x+y)^3,
H=(x+y)^3, deg W=5: W=lambda(x+y)^5.
```

Thus `K=1` genuinely permits nonzero `W`; excluding that locus for a planar
counterexample uses the cited one-place-at-infinity theory, not this proof.

On `K>=2`, `(T)` discharges the `W=0` hypothesis in THM-3016's displayed
cross-term identity.  When `a>b`, it gives

```text
Jac(P_(n-1),Q_(m-1))
 =kappa(a-b)H^(a-b-1)Q_(m-1)j.
```

When `a=b` the factor `a-b` makes the cross term zero and the expression is
handled separately; no negative power of `H` is intended.  Any further
Euclidean-tower claim retains the exact scope and branch qualifications of
THM-3016.

## 5. Verification and scope

The companion computes exact nullities for five `K>=2` hostile forms and two
`K=1` controls, checks both coprimality identities on a finite hostile box,
and verifies the subleading product identity on symbolic samples.  It is a
referee for the algebra, not a substitute for the common-form lemma.

This theorem constrains the leading homogeneous layers; it does not decide
JC(2).  It makes no electrical, arithmetic-row, or higher-dimensional
reduction.  The conclusion of the pre-repair file is retained, while its
false division is withdrawn.

Referee: `jc2_W_forced_thm3025.py` -- `ALL CHECKS PASSED`.
