---
id: THM-3615
title: "LRC pointed root-difference lift flag"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Retaining root difference before marginalization uniquely lifts the
  THM-3602 centered flag to K2sharp<R4sharp<Q6sharp with dimensions 2<4<6.
  Difference marginalization is an isomorphism on both the raw and centered
  six-carriers, but it collapses one line of a hidden rank-two raw-to-centered
  augmentation.  That line exactly accounts for both anomalous marginal
  intersection dimensions.  This is a static theorem over one pinned finite
  field after branch digit has been summed; no chronology, current,
  characteristic-zero transfer, row exclusion, or LRC(14) result is claimed.
source: kps-s188 + agent Godel / THM-3602 continuation, 2026-08-21
audit: >
  PASS -- independent hostile reconstruction reproduced the full parent census,
  direct root-difference inversion, both marginal isomorphisms, the complete
  sharp flag and intersection ledger, the rank-two augmentation and its lost
  marginal line, and all affine, reversal, Fourier, coordinate-swap, and
  duplicate-state negative controls.  The companion independently reconstructs
  the THM-3602 flag and the clean-room four-way pointed tensor, verifies both
  marginal isomorphisms, every sharp rank/intersection, the rank-two
  augmentation and its one-dimensional marginal kernel, slice and Fourier
  ledgers, hostile affine/reversal/transpose controls, and the failure of the
  tempting duplicate-state physical identification.  Normal and optimized
  runs exit zero and are identical to the stored transcript after LF
  normalization.
depends_on:
  - THM-3602-lrc-centered-a4-flag-inside-pointed-six-carrier
  - THM-3593-lrc-common-a4-anova-graph-flag
related:
  - THM-3585-lrc-common-a4-channel-plane-and-centering-complement
script: 04-computation/lrc_pointed_difference_lift_flag_thm3615.py
output: 05-knowledge/results/lrc_pointed_difference_lift_flag_thm3615.out
script_sha256: 25a37a623f1d955f2e0b6723e0ab103bac32c9ccfabc1d0d9146123be0f6b7c7
output_sha256: 493cf285d16881fbb7018724b47f44dbcbc183b984aa7602c30043628d0a5f40
semantic_sha256: 7994b4ff0428249ac299e6a90517de7e36d0329603f7473d648e199677b1e9c9
affine_atlas_sha256: b90a1366ecfb55134816025974c931f2555a2b6cf538e57c674e6b46bdb593c2
hash_basis: LF-normalized bytes for files; canonical JSON for semantic ledgers
---

# THM-3615 -- LRC pointed root-difference lift flag

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This restores exactly one coordinate discarded in THM-3602 and
identifies what its marginalization destroyed.  The result is stronger than a
rank coincidence but remains static and branch-summed.

Work over the pinned field and two thirteen-element coordinate sets

```text
F=F_p,       p=755373809845391722745761,
D=F13(root difference),                  T=F13(relation).          (1)
```

The shared primitive thirteenth root is the one pinned in THM-3602.  Let
`mathcal T(x,r,s,t)` be the independently reconstructed clean-room four-way
tensor, where `x` is one of the six pointed states, `r` is the branch digit,
`s in D` is root difference, and `t in T` is relation.  The six states are

```text
(0,0), (1,0), (1,6), (3,6), (3,12), (2,12).             (2)
```

## 1. Sharp raw and centered carriers

Sum the branch digit but retain difference and relation:

```text
g_x(s,t)=sum_(r in F13) mathcal T(x,r,s,t),
P6sharp=span_F{g_x:x in (2)} subset Fun(D x T,F).         (3)
```

Let difference marginalization and relation centering be

```text
mu(f)(t)=sum_(s in D)f(s,t),
C_T(f)(s,t)=f(s,t)-1/13 sum_(u in T)f(s,u),
Q6sharp=C_T(P6sharp).                                    (4)
```

The two marginal maps are exact isomorphisms onto the THM-3602 carriers:

```text
mu:P6sharp  -> P6,                 ranks (6,6), kernel 0,
mu:Q6sharp  -> Q6,                 ranks (6,6), kernel 0.          (5)
```

Their canonical graph digests are respectively

```text
09c232a61bd3bd308cc4b66e5c22b8ff4c177cade263b6fea12b519dd355b2c2
0291a20c81074a28474f986e526522397d0e252c2bdb0ad8600462a4aeb59d26. (6)
```

Because THM-3602 proves `K2<R4<Q6`, use the second isomorphism in `(5)` to
define the unique typed lifts

```text
K2sharp=(mu|Q6sharp)^(-1)(K2),
R4sharp=(mu|Q6sharp)^(-1)(R4).                           (7)
```

Then

```text
K2sharp subset R4sharp subset Q6sharp,
dimensions          2 < 4 < 6.                          (8)
```

The canonical sharp-space digests are

```text
K2sharp  7b73626f894ea6730fde769fbe57abbc3dca13819c3ab38e3c62aa6fbe2a95af
R4sharp  f0ca0f513a106e883962ea97f1df6ed56ddf584472a5186d02275c92f3570936
P6sharp  6014f99bee469da11163238341e815f6901d1460ec2f3d1ac62515b65470fda5
Q6sharp  167fb13d3fea42a269881e4be9bd7d2e7405bc5a19c50a7d01c0f9c5ec0c0240. (9)
```

## 2. The hidden rank-two augmentation

The raw and centered sharp carriers are not the same subspace.  Their exact
intersections are

```text
R4sharp intersect P6sharp = K2sharp,                 dimension 2,
P6sharp intersect Q6sharp,                          dimension 4,
R4sharp intersect Q6sharp = R4sharp,                dimension 4. (10)
```

Relation centering restricts to an isomorphism

```text
C_T:P6sharp -> Q6sharp.                                  (11)
```

Write its inverse uniquely as

```text
j(q)=q+LambdaSharp(q),
LambdaSharp(q) in Fun(D,F) tensor <1_T>.                 (12)
```

Thus the correction is constant in relation but may depend on root
difference.  It has the exact ledger

```text
A2sharp=im LambdaSharp,                    dimension 2,
ker LambdaSharp=P6sharp intersect Q6sharp, dimension 4,
rank(LambdaSharp|R4sharp)=2,
ker(LambdaSharp|R4sharp)=K2sharp.                         (13)
```

THM-3593 has its own rank-two map `U:R4->E2` with kernel `K2`.  Equations
`(5),(7),(13)` therefore induce the well-defined quotient isomorphism

```text
E2 -> A2sharp,
U(r) |-> LambdaSharp((mu|Q6sharp)^(-1)r).                (14)
```

This is an exact linear identification of quotient planes.  It is not a
physical identification of the state and difference coordinates.

That qualification is forced by an exact hostile control.  Since the six
pointed generators are independent, write a vector of `R4sharp` in the
ordered coefficient coordinates of `(2)` and set

```text
Delta(c)=(c_1-c_2,c_3-c_4).                            (15)
```

It is tempting to identify these two duplicate-state imbalances with the two
difference profiles in `A2sharp`, which would force the point correction
pattern

```text
(0,+H_1,-H_1,+H_2,-H_2,0).                            (16)
```

The exact reconstruction rejects this identification.  The pattern residual
has rank two, and

```text
rank(Delta|R4sharp)=2,
rank(Delta|K2sharp)=2,
dim(ker Delta intersect K2sharp)=0,
rank(Delta,LambdaSharp)=4.                              (17)
```

Thus `ker Delta` is not `K2sharp`, and `LambdaSharp` does not factor through
`Delta`.  Equation `(14)` identifies two abstract quotient planes only; it
does not supply the missing physical coefficient/address intertwiner.

Now marginalize the augmentation itself.  Exactly one line disappears:

```text
rank(mu|A2sharp)=1,                 dim ker(mu|A2sharp)=1,
rank(mu LambdaSharp|R4sharp)=1,     dim ker=3.            (18)
```

This one lost line explains both marginal anomalies from THM-3602:

```text
upstairs: dim(R4sharp intersect P6sharp)=2,
          dim(P6sharp intersect Q6sharp)=4;
downstairs: dim(R4 intersect P6)=3,
            dim(P6 intersect Q6)=5.                      (19)
```

The kernel line in `(18)` has digest

```text
bbf00a18f8ec98a9a07661f5df8c9413d3d493faac29ed5d459710332a9c994e, (20)
```

and `A2sharp` itself has digest

```text
33697ba84b3d4a3d6f64ea4577bcef00f5a3be511a3de0cbe9ff3b1bdf6914d3. (21)
```

## 3. Difference support, Fourier support, and hostile covariance

Every sharp space vanishes identically on the same-root slice:

```text
f(0,t)=0.                                                (22)
```

Their physical support is all 156 cells with `s!=0`.  Difference-slice ranks,
ordered by `s=0,...,12`, are

```text
K2sharp: (0,1,1,1,1,1,1,1,1,1,1,1,1),
R4sharp, P6sharp, Q6sharp:
         (0,3,3,3,3,3,2,2,3,3,3,3,3).                  (23)
```

At each fixed relation value, restriction to the difference fibre has full
space rank `2,4,6,6`, respectively.  For the two-dimensional Fourier
transform

```text
fhat(a,b)=sum_(s,t)f(s,t) zeta^(as+bt),                  (24)
```

the exact support counts are

```text
space      Fourier cells     relation frequencies
P6sharp        169           all b
Q6sharp        156           b!=0
R4sharp        156           b!=0
K2sharp        130           b notin {0,+1,-1}.          (25)
```

In particular, Fourier support alone does not distinguish `R4sharp` from
`Q6sharp`.  The hostile covariance atlas gives:

- point reversal has permutation `(5,4,3,2,1,0)`, graph rank six, and
  symmetric/antisymmetric ranks `(3,3)`;
- no map `(s,t)->(-s,at+b)`, even with an overall sign, realizes that reversal;
- identity is the only stabilizer among relation-affine maps and separately
  among difference-affine maps; every other one of the 155 maps has zero
  intersection with each of `Q6sharp,R4sharp,K2sharp`;
- the relation- and difference-affine histograms coincide entrywise, but this
  is not swap symmetry: every sharp space has zero intersection with its
  `(s,t)->(t,s)` transpose.

The affine flag-atlas digest is

```text
b90a1366ecfb55134816025974c931f2555a2b6cf538e57c674e6b46bdb593c2. (26)
```

## 4. Exact boundary and next lawful bridge

The companion reconstructs all spaces from disjoint parents and verifies the
displayed maps, bases, intersections, slices, transforms, and affine atlases.
Reproduce with

```bash
python -B 04-computation/lrc_pointed_difference_lift_flag_thm3615.py
python -B -O 04-computation/lrc_pointed_difference_lift_flag_thm3615.py
```

This theorem restores root difference but not branch resolution: the branch
digit `r` was summed in `(3)` before `P6sharp` was formed.  The inverse of `mu`
therefore supplies a unique static difference profile for each vector in
`R4`, but no branch-resolved realization compatible with the source/channel
coefficients of THM-3593.

The next lawful object is a typed lift into

```text
Fun(point x branch_digit x root_difference x relation,F) (27)
```

whose branch marginal is `R4sharp` and which intertwines owner/address and
temporal-entry data.  Without that domain-level intertwiner there is no legal
chronology, physical current, integer or characteristic-zero transfer,
forbidden LRC row, or LRC(14) consequence.  The present theorem claims none.

**QED.**
