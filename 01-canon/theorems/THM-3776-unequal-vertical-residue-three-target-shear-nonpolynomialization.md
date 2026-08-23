---
id: THM-3776
title: "Unequal vertical residues obstruct three natural rational target shears"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE / PENDING INDEPENDENT HOSTILE AUDIT.
  Two distinct nonzero simple-pole coefficients over one common smooth
  vertical parameter cannot be made simultaneously regular by a rational
  target-shear word q-p or q-p-q.  The variable-exchanged p-q and p-q-p
  statement holds when p, rather than q, is the common uniformizer.  The
  proof candidate includes all finite, polar, slope-zero, and higher-contact
  branches and arbitrary regular tails.  Starting in the opposite orientation
  can regularize the two completed DVRs locally; an exact three-shear witness
  keeps the scope sharp.  No claim may be used before independent audit.
source: root + jc_sparse_direct_search / 2026-08-23
audit: PENDING INDEPENDENT HOSTILE AUDIT.
depends_on: []
related:
  - THM-3709-cohn-alternating-two-by-two-elementary-decoration-nonentry
  - THM-3719-cohn-aligned-identity-shortening-laurent-nonentry
  - THM-3758-quadratic-radial-carrier-rational-exact-split-fibre-nonentry
  - THM-3771-arbitrary-profile-cubic-radial-carrier-log-canonical-nonentry
  - THM-3773-quadratic-rational-keller-cover-degree-and-target-word-nonpolynomialization
script: 04-computation/jc2_unequal_vertical_residue_target_shear_thm3776.py
output: 05-knowledge/results/jc2_unequal_vertical_residue_target_shear_thm3776.out
script_sha256: 2e077eeeda339b1813f410058a79b9516bf1d8eb1ff8d8d18b1b87baf1c0fdf1
output_sha256: ad268b187407cd58151a4cfbbe745effb598d89a9f55c0c0a02579de01a70553
semantic_sha256: 7d797e83e0a6d9fc742fcaac77720b6aaa0687f397c0ae5f21ad90949c80f29b
hash_basis: raw LF bytes
---

# THM-3776 -- unequal vertical residues survive three natural target shears

**RESERVED / PROVISIONAL PROOF CANDIDATE / PENDING INDEPENDENT HOSTILE
AUDIT.**  This is a local theorem about rational target operations.  It is
not an all-word classification.

Let `k` be a field of characteristic zero and let `K/k` be a field equipped
with two discrete valuations `nu_1,nu_2`.  Suppose the completed valuation
rings contain one common element `q` as a normalized uniformizer and

```text
Ohat_i = kappa_i[[q]],
p=a_i/q+b_i+O(q),                  a_i in k*,
a_1!=a_2.                                                (1)
```

The residue fields `kappa_i` may be larger than `k`, and the regular tails
`b_i+O(q)` are arbitrary.  Only the two leading coefficients are required to
belong to the common coefficient field.

For `F in k(s)` define the rational target shears

```text
S_q(F):(p,q) |-> (p,q+F(p)),
S_p(F):(p,q) |-> (p+F(q),q).                            (2)
```

All evaluations below are assumed defined in `K`.  Then:

1. no word `S_p(G) S_q(F)` can put both output coordinates in both
   `Ohat_1` and `Ohat_2`;
2. no word `S_q(H) S_p(G) S_q(F)` can put both output coordinates in both
   completed valuation rings.

Equivalently, unequal nonzero vertical principal coefficients obstruct every
natural alternating rational target word of length at most three that starts
by shearing the common regular parameter `q` with the polar coordinate `p`.
After exchanging `p` and `q`, the exact dual statements hold for `p-q` and
`p-q-p` when

```text
p is the common uniformizer,          q=a_i/p+O(1).      (3)
```

The hypotheses are realized by the two split components in THM-3758: the
rational primitive has opposite nonzero principal coefficients.  They also
apply to every pair of distinct nonzero addresses in THM-3771.  Thus neither
reciprocal nor Mobius target corrections can repair those two divisors in
two or three naturally oriented steps.

## 1. Two shears: reciprocal terms only translate the invoice

Write `t=q`.  Set

```text
q_1=q+F(p),                    p_1=p+G(q_1).            (4)
```

If `q_1` is regular at both valuations, then `F` has no pole at infinity.
Its Laurent expansion there has the form

```text
F(s)=f_0+f_1/s+O(s^-2).                                (5)
```

Substitution of `(1)` gives, independently of every regular tail,

```text
q_1-f_0=c_i q+O(q^2),
c_i=(a_i+f_1)/a_i.                                     (6)
```

If `c_i=0` for either valuation, `q_1-f_0` has order at least two.  A
rational `G` regular at `f_0` cannot cancel the simple pole of `p`; a pole of
order `m>=1` in `G` produces order at most `-2m`, which also cannot cancel a
simple pole.  This is the slope-zero or higher-contact exit.

Otherwise `G` can cancel the simple pole only if it has a simple pole at
`f_0`:

```text
G(y)=g_-1/(y-f_0)+O(1).                                (7)
```

The coefficient of `q^-1` in `p_1` is

```text
a_i+g_-1/c_i.                                          (8)
```

Regularity therefore forces

```text
g_-1=-a_i c_i=-(a_i+f_1).                              (9)
```

One common `g_-1` cannot satisfy `(9)` for `a_1!=a_2`.  This proves the
two-shear statement.

The requested special operations are visible in one line.  For

```text
F(p)=lambda/p
```

condition `(9)` reads `g_-1=-(a_i+lambda)`; choosing
`lambda=-a_i` merely enters the fatal higher-contact branch.  For a Mobius
function

```text
F(p)=(alpha p+beta)/(gamma p+delta),                   (10)
```

`gamma=0` either makes `F` constant or leaves `q_1` polar.  If `gamma!=0`,

```text
f_0=alpha/gamma,
f_1=(beta gamma-alpha delta)/gamma^2,                  (11)
```

so the Mobius step again translates every `a_i` by the same number and
cannot equalize two of them.

## 2. Three shears: the polar branch repeats the same difference

Now allow

```text
q_2=q_1+H(p_1)                                         (12)
```

and suppose the final pair `(p_1,q_2)` were regular at both valuations.  If
`q_1` is regular, Section 1 already says that `p_1` cannot be regular.  It
remains to treat `q_1` polar.

Let its pole order be `d`.  Then `F` has pole order `d` at infinity.  To
cancel the order-one pole of `p`, the leading polar order `e` of `G` at
infinity must obey

```text
d e=1.                                                  (13)
```

Both orders are positive integers, so `d=e=1`.  Write the required Laurent
jets as

```text
F(s)=A s+B+C/s+O(s^-2),                 A!=0,
G(y)=-y/A+D+E/y+O(y^-2).                              (14)
```

The leading coefficient `-1/A` is forced by cancellation.  Substitute the
full expansion

```text
p=a_i/q+b_i+d_i q+O(q^2).                              (15)
```

All `b_i,d_i` terms cancel, leaving

```text
p_1=ell+c_i q+O(q^2),
ell=D-B/A,
c_i=(E-C-a_i)/(A a_i).                                 (16)
```

The finite value `ell` is the same on both divisors.  If `c_i=0`, a pole of
`H` at `ell` has order at least two in `q` and cannot cancel the simple pole
of `q_1`; a regular `H` cannot cancel it either.  Otherwise `H` must have a
simple pole

```text
H(x)=h_-1/(x-ell)+O(1).                                (17)
```

Since `q_1=Aa_i/q+O(1)`, cancellation in `q_2` forces

```text
h_-1=-A a_i c_i=a_i+C-E.                               (18)
```

Again the two requirements differ by `a_1-a_2`.  This contradiction proves
the three-shear statement.  It also shows exactly why arbitrary regular
tails cannot help: they disappear before the first surviving coefficient in
`(6)` and `(16)`.

The dual assertion `(3)` follows by exchanging the coordinate names in the
same completed-DVR calculation.  This is a statement about regularity, so no
Jacobian-sign convention is hidden in that exchange; each shear in `(2)` is
itself determinant one.

## 3. The opposite starting orientation really can escape locally

It would be false to apply the theorem to a `p-q-p` word while retaining the
original hypothesis `(1)`.  A first shear of the polar coordinate by the
uniformizer can delete one residue and separate the two divisors before the
next reciprocal correction.

The smallest exact witness takes two embeddings into `k((t))`:

```text
q=t,
D_1: p=1/t,                       D_2: p=3/t.           (19)
```

Apply

```text
p_1=p-1/q,
q_1=q+p_1/(p_1+1),
p_2=p_1-1/(q_1-1).                                      (20)
```

Then

```text
D_1: (p_1,q_1,p_2)=(0,t,-1/(t-1)),
D_2: (p_1,q_1,p_2)=(2/t, t+2/(t+2), 1/(t+1)).          (21)
```

Both final functions `(p_2,q_1)` are regular in both completed DVRs.  This
does not produce a polynomial pair: the denominators in `(20)` introduce
new global divisors that the two local charts do not see.  It does prove that
the orientation and the word-length scope above are sharp, and it identifies
the next live state as a global divisor ledger rather than another local
residue calculation.

## 4. Exact controls and relation to the Cohn canon

The companion named in the metadata verifies the Laurent jets `(5)--(9)`
and `(14)--(18)` with arbitrary constant and linear regular tails; reciprocal
and all nondegenerate Mobius branches; the slope-zero order jump; finite and
polar order exits; several unequal residue pairs; the exact positive local
witness `(19)--(21)`; and the variable-exchanged dual.  Normal and optimized
executions must byte-match the frozen transcript before promotion.

THM-3709 and THM-3719 study alternating elementary decorations of source
Jacobian matrices over a polynomial ring.  The present operations instead
are rational symplectic changes in the two-dimensional target field, and the
invariant is a pair of completed-DVR principal coefficients.  The mechanisms
are analogous Euclidean words, but neither theorem implies the other.

The next counterexample-oriented question is now precise: can the locally
successful opposite word `(20)` pay every newly introduced denominator
divisor on an actual THM-3758 or THM-3771 surface?  **QED, conditional on
independent hostile audit.**
