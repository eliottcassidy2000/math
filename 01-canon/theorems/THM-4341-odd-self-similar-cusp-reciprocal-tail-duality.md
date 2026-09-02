---
id: THM-4341
title: "Odd self-similar cusp reciprocal-tail duality"
status: >
  PROVED FORMAL-LOCAL + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For every
  odd cusp order m=2g+1 and split order 1<=r<m in characteristic zero, the
  normalized exceptional tail has genus g-floor(r/2), while the persistent
  A_(r-1) defect has delta floor(r/2). The involution r<->m-r reciprocates
  the Newton slope and swaps those two quantities. Complementary valuation
  excesses in common sigma-units have constant product, and the sharp
  integral conductor-buffer threshold is B=g. Oriented split types and their
  reciprocal-pair quotient have exact natural-number indexings in doubled
  and ordinary triangular blocks. No isomorphism between reciprocal germs,
  global Keller extinction, even-order cusp statement, or positive-
  characteristic extension is claimed.
source: root + hostile-referee agent / planar-Jacobian next-sharp wildcard, 2026-09-02
depends_on: []
related:
  - THM-4289-a23-blowdown-observer-kahler-dualizing-quotient
  - THM-4340-u-zero-repeated-cusp-planar-jacobian-extinction
mistake_firewall:
  - MISTAKE-540
primary_script: 04-computation/odd_self_similar_cusp_reciprocal_duality_thm4341.py
primary_output: 05-knowledge/results/odd_self_similar_cusp_reciprocal_duality_thm4341.out
primary_script_sha256: 3bb06e3c9de90c96fd2b110f02d25fd723b5ea1112d6a769bad023aafd6eb25b
primary_output_sha256: d5cc1aba9c3b205c21140e94cd9df87f49b5e59489cf67c0fc39646da0f20afc
independent_audit_script: 04-computation/odd_self_similar_cusp_reciprocal_duality_independent_audit_thm4341.py
independent_audit_output: 05-knowledge/results/odd_self_similar_cusp_reciprocal_duality_independent_audit_thm4341.out
independent_audit_script_sha256: 9e1d9b62f1839e97e0ba4405177dc33fb057b1b236fcfe6c79f4cc61e984c2bc
independent_audit_output_sha256: d07e9b50dd8f93292103bfde90b4e06f56dc2f1f4c3566527fc59a8eeda1bc40
hostile_referee_script: 04-computation/odd_self_similar_cusp_reciprocal_duality_hostile_referee_thm4341.py
hostile_referee_output: 05-knowledge/results/odd_self_similar_cusp_reciprocal_duality_hostile_referee_thm4341.out
hostile_referee_script_sha256: 2975f895eceefaee3f5636a0f6ec0e6eb87a900573509deecefd463ccd5d0515
hostile_referee_output_sha256: c9c95a416f7f8d4e5f0b4e4b076201385332a4e76c704a12aab7d6a4eb693fda
hash_basis: raw LF bytes
audit: >
  PASS AFTER SCOPE REPAIR. The dependency-free primary path checks 62,750 oriented types through
  g=250, honest integer base-change weights, genus/delta reciprocity,
  differential products, both index bijections, and scope hostiles. The
  import-independent SymPy path constructs 6,480 tail polynomials over Q(c),
  checks squarefreeness and odd degree, independently inverts both indexes,
  and exhibits failures when m is even, c is zero, or the characteristic
  divides a branch exponent. A third hostile path refutes geometric
  isomorphism, naive plane-closure smoothness, full-order reciprocity,
  even-order extension, and globally injective bare odd-square labels.
  The theorem retains only the repaired numerical and formal-local claims.
  Normal and optimized streams byte-match all three frozen outputs.
---

# THM-4341 -- Odd self-similar cusp reciprocal-tail duality

**PROVED FORMAL-LOCAL + VERIFIED-EXACT + INDEPENDENTLY AUDITED. THIS IS A
LOCAL DEGENERATION AND INDEXING THEOREM, NOT A NEW JACOBIAN-CONJECTURE
STRATUM.**

## 1. Formal-local statement

Let `k` be an algebraically closed field of characteristic zero, let

```text
m=2g+1>=3,                  1<=r<m,                      (1)
```

and let `u(rho)` be a unit of `k[[rho]]`. Over `k[[t]]`, consider the
completed germ

```text
C_(m,r): x^2=z^m-(tz)^r u(tz)
        =z^r[z^(m-r)-t^r u(tz)].                         (2)
```

Then, after a finite base change, simultaneous normalization of the
horizontal `A_(r-1)` factor, and an integer weighted blowup, the reduced
nontrivial exceptional curve is

```text
E_(m,r): Y_0^2=Z^epsilon[Z^(m-r)-c],
epsilon=r mod 2,                         c=u(0)!=0.       (3)
```

It is connected and smooth, its smooth projective model has one place above
infinity, and

```text
g(E_(m,r))=g-floor(r/2)=floor((m-r)/2).                  (4)
```

The horizontal singularity removed before `(3)` has

```text
delta(A_(r-1))=floor(r/2),                               (5)
```

so tail genus plus persistent delta is exactly the original cusp delta `g`.

On the parameter set `(1)`, the numerical classification involution

```text
iota(r)=m-r                                                (6)
```

has the exact effects

```text
lambda_r=r/(m-r)       -> 1/lambda_r,
g(E_(m,r))             -> delta(A_(r-1)),
delta(A_(r-1))         -> g(E_(m,r)).                    (7)
```

This is a duality of the classified invariants and scales. It does **not**
assert that the two completed germs or their tail curves are isomorphic.

Finally, if `t=sigma^d` and a relative residue has the form

```text
omega=sigma^k z^B dz/x,                                 (8)
```

its order excess beyond `k`, measured in the common pre-clearing
`sigma`-valuation, is

```text
chi=B+1-m/2,                    Delta_r=chi d lambda_r.   (9)
```

For fixed `(m,d,B)`, reciprocal contributions satisfy

```text
Delta_r Delta_(m-r)=(chi d)^2.                           (10)
```

For integral `B`, the coefficient in `(9)` is positive for every split
exactly when `B>=g`; at the boundary it is `1/2`. Thus `B=g` is the sharp
integral buffer threshold. If also `k>=0`, every resulting total order is
positive.

## 2. Honest base change and normalized tail

No fractional root-coordinate cover is required. Make the dominating base
change and substitution

```text
t=tau^[2(m-r)],       z=tau^[2r]Z,       x=tau^[mr]Y.   (11)
```

All three terms in `(2)` have `tau`-weight `2mr`. Dividing by that common
power and setting `tau=0` gives

```text
Y^2=Z^r[Z^(m-r)-c].                                     (12)
```

Write `r=2q+epsilon` and normalize the persistent square factor by
`Y_0=Y/Z^q`. This is `(3)`. The normalized total chart is

```text
Y_0^2=Z^epsilon[Z^(m-r)-u(tau^(2m)Z)],                  (12a)
```

whose special exceptional curve occurs reduced. Its branch polynomial has
degree

```text
d_branch=m-r+epsilon.                                   (13)
```

Because `m` is odd, `(13)` is odd. The factor `Z^(m-r)-c` has only simple
nonzero roots in characteristic zero, and when `epsilon=1`, zero is also a
simple root. Hence `(3)` is squarefree. Its smooth projective model has one
place above infinity and genus

```text
(d_branch-1)/2=g-floor(r/2).                            (14)
```

The finite `Z=0` chart consists of smooth points after normalization. The
local curve has one boundary place above `Z=infinity`; identifying its global
attachment requires the complement chart and is deliberately not claimed
here. Any further toric or cyclic-quotient resolution inserts rational
chains only. Equations `(5)` and `(14)` now give the delta partition.

For the complement `m-r`, the elementary odd-sum identity

```text
floor(r/2)+floor((m-r)/2)=g                              (15)
```

proves the swap in `(7)`. The balanced-ray slope in `(11)` is
`r/(m-r)`, so `(6)` reciprocates it.

The oddness and unit hypotheses are load-bearing. If `m=2g` is even,
`d_branch` is even, the smooth projective model has two infinity places,
tail genus plus persistent delta is `g-1`, and `r=g` is fixed by
reciprocity. If `c=0`, the displayed branch polynomial is not squarefree.
In positive characteristic dividing `m-r`, its nonzero roots can be
inseparable; for example `Z^3-1=(Z-1)^3` in characteristic three.

## 3. Reciprocal valuation-excess law and sharp buffer

On the rational scale suggested by `(11)`, measured in `sigma`-orders,

```text
ord(z)=d*r/(m-r),             ord(x)=(m/2)ord(z).        (16)
```

Pulling `(8)` to the relative tail chart contributes one copy of
`ord(z)` from `dz` and subtracts `ord(x)`. This gives

```text
ord_E(omega)=k+(B+1-m/2)d*r/(m-r),                       (17)
```

which proves `(9)`. This is the excess beyond `k` in a common, possibly
fractional, pre-clearing `sigma`-valuation. Replacing `r` by `m-r` replaces
the last fraction by its reciprocal, proving `(10)`. Full orders including
`k`, or orders measured after separately minimal base changes, need not have
constant product.

Since `m=2g+1`,

```text
B+1-m/2=B-g+1/2.                                        (18)
```

For integral `B`, `(18)` is positive exactly at `B>=g`. At `B=g-1` it is
`-1/2`, so with `k=0` every tail contribution is negative. This proves
sharpness of the coefficient threshold. A positive `k` may rescue a
particular below-threshold row; no stronger converse is claimed.

The threshold is the formal operation-level bridge to THM-4340: `B=g` is
also the normalization-conductor exponent in the `z` coordinate, and its
remaining half-step is exactly what pays every self-similar tail scale.

## 4. Exact natural-number indexings

Let

```text
S={(g,r):g>=1, 1<=r<=2g}.                               (19)
```

The oriented index

```text
n(g,r)=g(g-1)+r                                         (20)
```

is a bijection `S -> N_(>0)`. Indeed, fixed `g` occupies the consecutive
block

```text
2T_(g-1)+1, ..., 2T_g,              T_j=j(j+1)/2.       (21)
```

and adjacent blocks meet without gaps or overlap. The reciprocal involution
reflects this block:

```text
n(g,r)+n(g,m-r)=2g^2+1.                                 (22)
```

There is a second exact index after quotienting by reciprocity. Put

```text
j=2r-m,                 h=(|j|+1)/2,
N(g,{r,m-r})=T_(g-1)+h.                                 (23)
```

Here `j` is a nonzero odd integer, `1<=h<=g`, and `N` is a bijection
`S/iota -> N_(>0)`. The sign of `j` is the one-bit sidecar that recovers the
oriented split:

```text
r=g+h       if j>0,                 r=g+1-h if j<0.     (24)
```

Moreover `j^2=(2h-1)^2`. With the genus block `g` retained, replacing the
`h`-th positive odd square by the natural number `h` is lossless on the
reciprocal quotient. Globally the triangular block in `N` is required, and
on the oriented set the sign bit is required as well. This is the precise
version of the odd-square compression suggested by the repository's sequence
work.

The reduced rational `lambda_r` is a valid Stern--Brocot address for the
exceptional ray, and `(6)` sends it to its reciprocal. That address alone is
not a complete cusp type: `(m,r)=(3,1)` and `(9,3)` both give `lambda=1/2`.
The total order `m`, orientation, coefficient `c`, and residue data remain
mandatory sidecars.

## 5. Why the tournament carrier stops at a preorder

For a multi-term degeneration, compare terms by their valuations. Away from
equalities this gives a transitive orientation by numerical order. At
equality, however, the initial polynomial `(12)` appears and carries the
tail's genus. Arbitrarily orienting the tied pair deletes precisely that
polynomial and can turn a positive-genus tail into an invisible event.

Thus the lawful finite carrier is a valuation preorder together with the
initial polynomial on each tie class. An ordinary tournament is neither
required nor lossless here. This observation is an exact information-loss
statement, not a transfer of tournament theorems to cusp geometry.

## 6. Reproduction and boundaries

Run both certificates in normal and optimized modes and byte-match the
frozen outputs:

```bash
python3 -B 04-computation/odd_self_similar_cusp_reciprocal_duality_thm4341.py
python3 -B -O 04-computation/odd_self_similar_cusp_reciprocal_duality_thm4341.py
python3 -B 04-computation/odd_self_similar_cusp_reciprocal_duality_independent_audit_thm4341.py
python3 -B -O 04-computation/odd_self_similar_cusp_reciprocal_duality_independent_audit_thm4341.py
python3 -B 04-computation/odd_self_similar_cusp_reciprocal_duality_hostile_referee_thm4341.py
python3 -B -O 04-computation/odd_self_similar_cusp_reciprocal_duality_hostile_referee_thm4341.py
```

The primary audit checks the formulas through `g=250`; the independent
symbolic path constructs every tail through `g=80` over `Q(c)`. The hostile
referee owns the excluded stronger claims and forced the scope repairs in
Sections 1--4. The proof there is uniform in `g,r`; finite replay is hostile
corroboration, not induction.

This theorem supplies a formal-local classifier and exact indexing scheme.
It does not supply a geometric isomorphism under reciprocity, a global owner
atlas, attachment multiplicities, component multiplicities,
proper-flat degree ledger, or seam entry needed for a new planar Jacobian
extinction theorem. Even `m`, `c=0`, and positive characteristic remain
outside its statement.
