---
id: THM-2082
title: "Rank-one code-wheel blindness and the translated-prime-grid terminal test"
status: >
  PROVED. A primitive scalar speed row has, at every prime p, the two-term
  evaluation-code enumerator 1+(p-1)z^w_p, where w_p counts the p-nondivisible
  speeds. Hereditary primitivity is exactly w_p>=2 for every p, so its
  one-deletion CRT wheel is empty. A translated p-grid nevertheless gives a
  rigorous terminal escape whenever its noncarrier union bound leaves one
  residue. An unbounded family in terminal sizes 7--10 freezes the entire
  prime code profile and passes the divisor, hereditary, quarter, height, and
  scalar-fold filters. A p=17 pair with identical code enumerators but 0
  versus 12 safe grid residues proves that Hamming support cannot replace
  projective residue incidence. This is a method boundary and a conditional
  escape tool, not LRC(14).
source: codex-2026-07-22-LRC-code-wheel-blindness
depends_on:
  - THM-2069
  - THM-2077
  - THM-2080
  - THM-2081
related:
  - THM-765
  - THM-2073
  - THM-2086
script: 04-computation/lrc14_rank_one_code_wheel_blindness_codex_20260722.py
output: 05-knowledge/results/lrc14_rank_one_code_wheel_blindness_codex_20260722.out
script_sha256: 7aab4e6ee7c3d43fad2bace64a2318b88be1b37fbffa95ee4628aec55e1abbd0
output_sha256: 0f048f47e44fc287c80a0c3825f5a5abf457f7e5cdc72dc388b364276f22caf8
hash_basis: working-tree bytes (LF)
---

# THM-2082 -- rank-one wheel blindness and translated prime grids

Put

```text
D_q={t in R/Z:||qt||<1/14},
E_h={t in R/Z:||ht||<1/7},
G_Q=intersection_(q in Q)(R/Z minus D_q).              (1)
```

This theorem separates two roles that were conflated in attempts to use
THM-2069 on a terminal scalar core. The code/cogirth wheel is an exact
**deletion-primitivity** carrier. It does not retain the projective residue
placement that decides whether a rational grid meets `G_Q minus E_h`.

## 1. Exact rank-one specialization

Let `Q={q_1,...,q_s}` be a positive integer row with `gcd(Q)=1`. For a prime
`p`, put

```text
w_p=#{i:p does not divide q_i}.                         (2)
```

The THM-2069 evaluation code is

```text
C_p={lambda(q_1,...,q_s):lambda in F_p},
W_(C_p)(z)=1+(p-1)z^(w_p),
cogirth(C_p)=w_p.                                      (3)
```

Consequently a nonzero rank-one projective direction survives deletion of
`k` coordinates exactly when

```text
w_p>k.                                                 (4)
```

In particular,

```text
Q is hereditarily primitive
  iff w_p>=2 for every prime p.                         (5)
```

When (5) holds, every one-deletion index is one and THM-2069's squarefree
one-deletion modulus is

```text
Q_1=1.                                                 (6)
```

### Proof

Multiplication by a nonzero `lambda` changes no zero coordinate, so every
nonzero codeword has precisely the support counted by (2). This proves (3),
and (4) is THM-2069's deletion theorem.

Because the whole row is primitive, `w_p>=1`. If `w_p=1`, deleting its unique
`p`-nondivisible entry leaves a row with gcd divisible by `p`. Conversely, if
some deletion has gcd divisible by `p`, at most the deleted entry is
`p`-nondivisible, so `w_p=1`. This proves (5). Every retained deletion row in
the hereditary case has gcd and hence index one, proving (6). QED.

This also explains why the density statement of THM-2069 has no rank-one
content: the primitive covectors of `Z` are only `+1` and `-1`. In an LRC
scalar core, height is a new arithmetic phase coordinate, not a parameter
sampled by that density theorem.

## 2. Translated-prime-grid escape lemma

For a real number `x`, write

```text
dist(x,14Z)=min_(n in Z)|x-14n|.                       (7)
```

Let `p` be prime, let `z` be a positive integer divisible by `p`, and let
`c` be an integer. Call `q` a **p-carrier** when `p|q`, and put

```text
w=#{q in Q:p does not divide q},
u_p=ceil(p/7),                 v_p=ceil(2p/7).          (8)
```

Assume every carrier in `Q` obeys

```text
dist(cq/z,14Z)>=1.                                     (9)
```

Consider the `p-1` translated nonzero grid points

```text
t_a=a/p+c/(14z),              a=1,...,p-1.             (10)
```

There are two escape criteria.

1. If `p` does not divide `h` and

   ```text
   w u_p+v_p<p-1,                                      (11)
   ```

   then some `t_a` lies in `G_Q minus E_h`.
2. If `p|h`,

   ```text
   dist(ch/z,14Z)>=2                                   (12)
   ```

   and

   ```text
   w u_p<p-1,                                          (13)
   ```

   then again some `t_a` lies in `G_Q minus E_h`.

Thus either criterion contradicts terminal containment `G_Q subset E_h`.

### Proof

If `p|q`, then `qa/p` is integral and

```text
||q t_a||=||cq/(14z)||>=1/14                           (14)
```

by (9). Every carrier is therefore safe for every `a`; equality is allowed
because `D_q` is open.

If `p` does not divide `q`, multiplication by `q` permutes the nonzero
residues modulo `p`. A translated `p`-grid has at most `ceil(p/7)=u_p` points
in any open circle arc of length `1/7`. Hence each noncarrier excludes at
most `u_p` values of `a` through `D_q`.

When `p` does not divide `h`, the same argument puts at most `v_p` grid
points in the guard arc, whose length is `2/7`. The union bound (11) leaves
one of the `p-1` candidates outside all these events. When `p|h`, equation
(12) instead makes every candidate guard-safe, and (13) leaves a candidate
after the `Q` exclusions. QED.

At the apex prime `p=7`, this becomes a particularly small conditional
classifier. If `7` does not divide `h`, at most three noncarriers close by
(11). If `7|h`, at most five close by (13). Taking `c=1` and a selected
`7`-carrier `z in Q`, the only escape from the lemma is explicit ratio
incompatibility:

```text
dist(q/z,14Z)<1 for some 7-carrier q,
or, in the carrier-guard branch, dist(h/z,14Z)<2.       (15)
```

This is useful leverage because it returns an actual speed ratio or an
actual safe grid point. It is not a uniform theorem: arbitrary terminal rows
can have too many noncarriers or can take the ratio branch (15).

## 3. A frozen unbounded family

Let

```text
L=lcm(2,3,...,14)=360360,
B_j=L 32^j,                         j>=1,
Q_(s,j)={1,2,...,s-1,B_j},          7<=s<=10,
h=13.                                                     (16)
```

This family has unbounded terminal height and simultaneously satisfies the
following exact properties.

1. It is divisor-complete through `14` and hereditarily primitive.
2. Its smallest multiple of four is `4`, while `B_j>7*4`, so it satisfies
   THM-2077's first recursive quarter-escape alternative.
3. The THM-2078 relative-height filter holds:

   ```text
   (13-s)h<2(s+1)B_j.                                  (17)
   ```
4. Its entire prime-code profile is independent of `j`:

   ```text
   w_p=(s-1)-floor((s-1)/p)+1_(p does not divide L).    (18)
   ```

   The minimum cogirths for `s=7,8,9,10` are respectively

   ```text
   3,4,4,5.                                            (19)
   ```

   Hence both the one- and two-deletion CRT wheels are empty throughout the
   family.
5. It passes the scalar containment invoice extracted from THM-2080.

To state item 5, for `q in Q` put

```text
g=gcd(q,h),          a=q/g,          b=h/g,
F(r)=r_bar(14-r_bar),       0<=r_bar<14,
C(q,h)=F(b+2a)-F(b-2a).                                  (20)
```

THM-2080's exact fold formula becomes

```text
I(q,h)=measure(D_q intersect E_h)
      =2/49+C(q,h)/(196ab).                            (21)
```

If `G_Q subset E_h`, the union bound on the guard complement necessarily
gives

```text
sum_(q in Q) C(q,h)/(ab)<=20(s-7).                     (22)
```

For (16), the left side is

```text
s=7       -103/65,
s=8       -103/65,
s=9       -231/130,
s=10      -733/390.                                   (23)
```

Thus every row passes (22), including its sharp rank-seven boundary.

Nevertheless every row has the same explicit strict escape:

```text
t=3/31 lies in G_(Q_(s,j)) minus E_13.                 (24)
```

Indeed `32=1 mod 31` and `L=16 mod 31`, so `B_j=16 mod 31` for all `j`.
The least absolute residues of `3q mod 31` for `q=1,...,9` are

```text
3,6,9,12,15,13,10,7,4,                                (25)
```

and the far speed has residue `14`. All are at least `3`, while
`3/31>1/14`. The guard has least residue `8`, and `8/31>1/7`, proving (24).

### Proof of the remaining family claims

The far speed is divisible by every integer through fourteen. Deleting any
entry other than `1` leaves `1`; deleting `1` leaves the coprime pair `2,3`.
This proves divisor completeness and hereditary primitivity. The quarter and
height assertions are immediate from `B_j>=L*32`.

The prime divisors of `B_j` are exactly those of `L`, because the only prime
dividing `32` is already a divisor of `L`. Among `1,...,s-1`, precisely
`floor((s-1)/p)` entries are divisible by `p`; this proves (18) for every
prime, not merely a tested range. Its minimum occurs at `p=2`, giving (19).

For (21), substitute THM-2080's reduced variables into its one-sided fold
identity and clear denominator `196ab`. If `G_Q subset E_h`, then
`E_h^c subset union_q(D_q intersect E_h^c)`. Comparing measures gives

```text
5/7<=s/7-sum_q I(q,h),                                 (26)
```

and (21) turns (26) into (22). Direct reduction modulo fourteen gives (23);
the far term is zero because `B_j/13` is divisible by fourteen. QED.

The conclusion is deliberately limited but decisive: no terminal height
contradiction can follow from only the complete prime code/cogirth profile
together with divisor completeness, hereditary primitivity, quarter escape,
the relative-height inequality, and the scalar fold invoice. Those data are
constant or uniformly admissible while `B_j` tends to infinity.

This does **not** refute THM-2081. Its relative Hunter edges retain the
three-frequency residue placement that the scalar and code filters discard,
and THM-2086 now detects the rank-seven member of this family uniformly by
its lacunary relative-Hunter cone.

## 4. Identical Hamming data, opposite rational-grid behaviour

The support blindness already appears at one prime. Put

```text
Q_cov=(1,360360,19,37,89,107,109,901093),
Q_gap=(1,360360,103,137,239,307,409,901171).           (27)
```

Both rows are pairwise coprime, divisor-complete through fourteen,
hereditarily primitive, and satisfy the odd quarter escape because their
last entry exceeds `(5/2)360360`. Modulo `17` their residues are

```text
Q_cov: (1,11,2,3,4,5,7,8),
Q_gap: (1,11,1,1,1,1,1,1).                            (28)
```

Every coordinate is nonzero, so both rank-one codes have exactly

```text
W(z)=1+16z^8,                  cogirth=8.               (29)
```

But on the nonzero `17`-grid,

```text
#{a in F_17^x:a/17 in G_(Q_cov)}=0,
#{a in F_17^x:a/17 in G_(Q_gap)}=12.                   (30)
```

For `Q_cov`, (28) contains one representative of every antipodal nonzero
class. For each nonzero `a`, some product `aq` is `+1` or `-1` modulo `17`,
and `1/17<1/14`. For `Q_gap`, only the projective classes of `1` and `11`
occur. Exactly `a=1,3,14,16` are dangerous, leaving

```text
a=2,4,5,6,7,8,9,10,11,12,13,15.                       (31)
```

This proves (30). The code enumerator knows the support size eight; it does
not know how those eight labelled coordinates occupy projective residue
classes.

## 5. What transfers to the active problems

The theorem gives a precise division of labour.

* **LRC(14).** Use the scalar code wheel only to sieve deletion primitivity.
  The next live coordinate is residue incidence localized to the guard
  complement: translated grids for a conditional exit, and THM-2081's
  projective three-frequency edges for the all-height rank-seven gate. The
  family (16) rules out trying to manufacture height from the older scalar
  sidecars alone.
* **Gaussian Moment / DvdK.** The same carrier boundary appears there.
  Existence or cardinality of balanced channels is sufficient in the unique
  channel stratum, but multiple coincident channels require coefficient phase
  and radial ownership. Hamming support is analogous to channel support, not
  to the signed channel sum. This theorem supplies no arbitrary-coefficient
  noncancellation result.
* **Extremal codes.** THM-2069's `[72,36,16]` gate likewise needs the labelled
  incidence of minimum cocircuit supports, not only `A_16`. Equations
  (27)--(30) are a small rank-one model of that exact loss.
* **Finite atlases.** The translated-grid lemma is constructive: failure
  produces the explicit ratio obstruction (15). It can therefore branch an
  atlas into a certified rational phase or a finite carrier-ratio packet,
  unlike a density-one statement.

These are mechanism transfers, not implications between the conjectures.
In particular, neither NC2/GMC(2) nor the code wheel currently proves
LRC(14).

## 6. Assumption challenge and Tournament Analysis

Before choosing the carrier, the natural candidate vertices were runners,
gaps, fixed circle sections, section boundaries, wall-crossing events,
residues, code supports, matroid cocircuits, and proof obligations. The
rank-one code support preserves exactly the predicate

```text
"which primes divide all but at most k specialized speeds?"               (32)
```

It destroys positivity, phase height, projective residue placement, safe-set
components, endpoint owners, and guard-relative pair intersections. For the
translated lemma, the faithful vertices are the nonzero `p`-grid residues,
labelled by which danger or guard event excludes them. This preserves the
rational escape predicate but forgets the intervals between grid points.

Ordering residues by clearance gives a transitive tournament: zero directed
cycles, singleton strongly connected components, and one Hamiltonian path
after ties are broken by residue label. The two rows in (27) can share every
Hamming fingerprint while differing as in (30), so that tournament is only a
priority list. The faithful object is the labelled residue-event incidence
hypergraph, or downstream the weighted restricted-event graph of THM-2081.

## Exact referee

The companion checks the family for `s=7,...,10` and `j=1,...,5`, including
`3,360` prime-profile identities through `p<=997`, all exact fold values, and
`190` phase inequalities. It audits both translated-grid branches on forty
prime/control rows and verifies the two `p=17` incidence sets exactly. Normal
and `python -O` transcripts byte-match and end in `RESULT=PASS`. QED.
