---
id: THM-2168
title: "Three-target second-depth majorization and the evaluation-shear boundary"
status: >
  PROVED. A three-band cover of the asymmetric source rectangle has the
  following odd-prime strengthening of THM-2148: either one target is the
  signed narrow source character, or at least two targets inherit the full
  determinant valuation of the source lattice. The mechanism is a complete
  classification of two translated finite-group bands after the unique
  descended target becomes inactive. In the all-independent (3,5,0) LRC
  profile, the largest and second-largest divided-blocker 13-adic depths
  dominate the largest and second-largest aligned depths. Rational guard-line
  incidence has the exact residual census (0,0), (0,1), (4,3), or (5,3);
  the latter two are precisely the still-open scalar 4+3 and 5+3 tails.
  Exact daughter pairs are the only singleton depth exceptions and give
  independent height-13 two-term relations. Bounded scalar relation rank does
  not close the remaining profile: evaluation-kernel shears preserve every
  W_H and the mod-13 profile while moving determinant and polar data without
  bound. This is a strict reduction and loss ledger, not an elimination of
  (3,5,0) or a proof of LRC(14).
source: codex-2026-07-24-LRC-Fano-rank-height-synthesis
depends_on:
  - THM-2080
  - THM-2138
  - THM-2139
  - THM-2148
  - THM-2164
related:
  - THM-2073
  - THM-2140
  - THM-2169
script: 04-computation/lrc14_three_target_second_depth_referee_thm2168.py
output: 05-knowledge/results/lrc14_three_target_second_depth_referee_thm2168.out
script_sha256: 1ae4e345188c8b40897d02c446d2edfb4c0aed44541535314064bb8d9db00077
output_sha256: 582b9f2a99bac4eaca784b41e5be3af2aaa69d54a4b21842ad208962aa73a7bf
hash_basis: working-tree bytes (LF)
---

# THM-2168 -- three-target second-depth majorization

Let `Gamma` be a rank-two character lattice,

```text
K=Hom(Gamma,R/Z),             D_w={X:||w.X||<=1/14}.  (1)
```

Take rationally independent nonzero characters `a,b` and three nonzero
targets `w_1,w_2,w_3`.  Suppose

```text
{X:||a.X||<1/7 and ||b.X||<1/14}
                         subset D_(w_1) union D_(w_2) union D_(w_3). (2)
```

For every odd prime `p`, put

```text
A=nu_p(det(a,b)).                                     (3)
```

Then at least one of the following holds:

```text
some w_j=+b or -b;                                    (4)

at least two distinct j satisfy
nu_p(det(a,w_j))>=A.                                  (5)
```

There is a more informative integral alternative behind (5).  Let

```text
L=ker(a) intersection ker(b),                         (6)
T={j:w_j|L=0}.
```

THM-2148 gives `T!=empty`.  If (4) fails, then either

```text
|T|>=2,                                               (7)
```

or, after relabelling `T={1}`, the other two restrictions are the same
nonzero order-two character and

```text
2w_2, 2w_3, w_2-w_3 in Za+Zb.                        (8)
```

Thus the factor two in (8) is the only integral loss, and it is invisible
to every odd-prime valuation.

## 1. Two translated bands that cover a finite group

We first prove the finite statement used after a descended target becomes
inactive.

> **Two-band lemma.** Let `F` be a finite abelian group, let `chi,psi` be
> nontrivial characters of `F`, and let `theta,phi in R/Z`.  If
>
> ```text
> F={x:||chi(x)+theta||<=1/14}
>      union {x:||psi(x)+phi||<=1/14},                 (9)
> ```
>
> then `chi=psi` is an order-two character and the two displayed sets are
> its complementary kernel cosets.

If a character has image order `m`, a translated closed arc of length
`1/7` contains at most

```text
floor(m/7)+1                                          (10)
```

points of its `m`-grid.  For `m>=3`,

```text
floor(m/7)+1<=m/7+1<m/2,                              (11)
```

while for `m=2` the maximum is exactly one point.  Hence each set in (9)
has size at most `|F|/2`, with equality only for an order-two character.
Coverage forces equality for both sets and forces the two half-size sets
to be disjoint.

An order-two character band, when nonempty, is one coset of its index-two
kernel.  Cosets of two distinct index-two subgroups intersect: the two
independent characters give a surjection to `C_2 x C_2`, on which every
pair of prescribed values occurs.  Disjointness therefore forces the two
kernels, and hence the two nonzero characters, to agree.  The two cosets
are complementary.  This proves the lemma, including closed endpoints.

## 2. Proof of the local depth theorem

The finite group `L` in (6) is the common kernel of the surjection

```text
(a,b):K -> (R/Z)^2.                                   (12)
```

By THM-2148, at least one target kills `L`; equivalently `T` is nonempty.
Exact annihilator duality says

```text
L^perp=Za+Zb.                                         (13)
```

If `|T|>=2`, taking determinants with `a` in (13) shows that
`det(a,b)` divides the determinants of two targets.  This proves (5).

Suppose now that `T={1}`.  Write

```text
w_1=ma+nb,                    m,n in Z.                (14)
```

If the whole source rectangle in (2) lay in `D_(w_1)`, the asymmetric
one-target lemma of THM-2139 would give

```text
w_1=+b or -b.                                         (15)
```

Assume (4) fails.  There is then a source phase `(x,y)` at which the first
band is inactive.  On the full `L`-fibre over that phase, the bands of
`w_2,w_3` must cover `L`.  Their restrictions are nontrivial because
`T={1}`.  The two-band lemma gives

```text
w_2|L=w_3|L!=0,              2(w_2|L)=0.              (16)
```

Consequently `2w_2`, `2w_3`, and `w_2-w_3` kill `L`, and (13) proves (8).
Taking determinants with `a` gives

```text
det(a,b) divides 2det(a,w_2),
det(a,b) divides 2det(a,w_3).                         (17)
```

Because `p` is odd, (17) implies (5).  This completes the local proof.
Notice that no assertion about which complementary order-two cosets are
chosen is needed; translated-fibre coverage itself supplies the two-target
cover at the one phase where `w_1` is inactive.

## 3. The rational guard-line census

Use the notation of THM-2148's LRC application:

```text
c_*j=13u_j,                    j=1,2,3,
c_1,...,c_5 aligned with the guard g modulo 13.       (18)
```

The five source containments are

```text
{||g.Y||<1/7,||c_i.Y||<1/14}
                  subset union_(j=1)^3 D_(u_j).       (19)
```

Let

```text
k=#{i:c_i in Qg},             r=#{j:u_j in Qg}.       (20)
```

We first prove the exact surviving census

```text
(k,r) in {(0,0),(0,1),(4,3),(5,3)}.                  (21)
```

If `k>0`, THM-2148's dependent clause gives `r>0`.  Let `alpha` generate
the saturated guard line.  The guard and the `k+r` line terminals have
distinct positive scalar coefficients; the guard coefficient is odd.  If

```text
2<=k+r<=6,                                            (22)
```

THM-2080 supplies a positive-measure set of scalar phases on which the guard
and all these line terminals are strictly safe.  On each corresponding
`alpha`-fibre, the other `8-k-r` terminal characters restrict nontrivially.
Their danger bands have fibre measure `1/7`, and

```text
8-k-r<=6.                                             (23)
```

Their union therefore misses positive fibre measure.  Fubini contradicts
the cover.  It follows that `k+r>=7`, so the caps `k<=5,r<=3` give

```text
k>=4,                    r>=2.                        (24)
```

The apparent case `(k,r)=(5,2)` is also impossible.  Its scalar line consists
of a thirteen-unit guard, five thirteen-unit aligned coefficients, and two
positive-valuation blocker coefficients.  THM-2138's all-depth `5+2`
annulus theorem supplies a positive-measure scalar escape.  The sole
remaining off-line blocker occupies only `1/7` of each circle fibre, so
Fubini again gives an escape from all eight terminals.  Hence `k>0` leaves
only `(4,3)` and `(5,3)`.

If `k=0` but `r>=2`, the same THM-2080/Fubini argument applies with the
`r` line blockers and `8-r<=6` off-line bands.  Therefore `r<=1`, proving
(21).  The last two cases in (21) are exactly the scalar `4+3` tail with one
off-line aligned label and the fully scalar `5+3` tail.  Neither is closed
by THM-2138.

## 4. The `(3,5,0)` second-depth invoice

First assume `k=0`, so all five `g,c_i` pairs are rationally independent,
and put

```text
delta_i=det(g,c_i),             epsilon_j=det(g,u_j). (25)
```

Order the valuations, using `nu_13(0)=infinity`:

```text
a_1<=...<=a_5,     {a_i}={nu_13(delta_i)},
beta_1<=beta_2<=beta_3,
                   {beta_j}={nu_13(epsilon_j)}.        (26)
```

For every aligned label, the local theorem gives the dichotomy

```text
u_j=+/-c_i for some j;                                (27)

or at least two blockers have valuation at least
nu_13(delta_i).                                       (28)
```

The LRC cocharacter makes every `u_j` and `c_i` positive, so the minus sign
in (27) is impossible.  We call

```text
u_j=c_i                                                (29)
```

an **exact daughter pair**.  Distinctness of the terminal values makes
daughter pairs a matching: one divided blocker cannot be the daughter of
two aligned labels, and one aligned label cannot have two equal divided
blockers.

Apply (27)--(28) to the two aligned labels of largest 13-adic depth.  If
either is not a daughter, it alone puts two blockers at least as deep.  If
both are daughters, matching puts their depths on two distinct blockers.
Together with THM-2148's original divisibility edge for the deepest label,
this proves

```text
beta_3>=a_5,                  beta_2>=a_4.             (30)
```

For the original blockers `c_*j=13u_j`, equation (30) becomes

```text
largest two blocker depths >= a_5+1, a_4+1.           (31)
```

In the line cases `(4,3)` and `(5,3)`, at least four aligned determinants
and all three divided-blocker determinants vanish.  Thus

```text
a_4=a_5=beta_2=beta_3=infinity,                       (32)
```

and (30) remains uniformly true.

Equation (31) is complementary to THM-2140's real polar certificates.
THM-2140 gives upper height gates in the one-blocker profiles; the present
three-blocker argument gives two depth obligations.  There is deliberately
no determinant-magnitude conclusion: a target on the guard line has
determinant zero and therefore infinite 13-adic depth, while contributing
no Archimedean lower bound.  No valid upper polar bound for the selected
targets follows from a three-band union, so (31) alone does not force
collinearity.

## 5. Exact daughters and bounded relation rank

Each daughter pair is stronger than a scalar coincidence:

```text
c_*j-13c_i=0 in Gamma.                                (33)
```

For `d` daughter pairs, the `d` relations (33) have disjoint paired supports
and are linearly independent, all with coefficient height `13`.  Whenever
the displayed terminal labels are coordinates of the original speed row,
or receive the same common dyadic lift as in THM-2073, (33) gives the same
support-two relation in the thirteen-speed lattice.

Thus:

1. two or three daughters already give two or three independent height-13
   relations;
2. with exactly one daughter in a distinct zero-safe thirteen-speed row,
   THM-2164's `dim_Q W_105>=2` supplies a second actual relation of height
   at most `105` independent of (33).

This is the precise interface with THM-2144/2164.  It does not say that the
second scalar relation vanishes as a character relation in `Gamma`.

## 6. Why scalar rank does not control the Fano carrier

The missing map is exposed by an exact shear symmetry.  Let

```text
h:Gamma->Z
```

be the LRC evaluation cocharacter and choose `0!=k in ker(h)`.  Keep `g`
fixed and, for arbitrary integers `N_i,M_j`, set

```text
c_i'=c_i+13N_i k,
u_j'=u_j+M_j k,                c_*j'=13u_j'.          (34)
```

Then

```text
h(c_i')=h(c_i),          h(c_*j')=h(c_*j),            (35)
c_i'=c_i mod 13,         c_*j'=0 mod 13.              (36)
```

Consequently every scalar relation lattice and every bounded span `W_H`,
for every `H`, is unchanged.  Positivity, distinctness of evaluated speeds,
and the blocker/aligned mod-13 profile are also unchanged.  The same is true
of the entire thirteen-card deletion-relation deck in THM-2169: deleting a
coordinate and then evaluating still sees exactly the same twelve scalar
values.

On the other hand, with `d=det(g,k)`,

```text
det(g,c_i')=det(g,c_i)+13N_i d,
det(g,u_j')=det(g,u_j)+M_j d.                         (37)
```

Here `d!=0`: otherwise `g` would be rationally parallel to `k` and hence
`h(g)=0`, contrary to the positive guard value.  Thus the determinant
magnitudes can be moved without bound while all scalar relation data stay
fixed.  Integral-span incidence and real polar bodies can change as well.

This gives the exact connection ledger:

```text
source:       bounded scalar relations among evaluated speeds;
map:          evaluate every character by h;
preserved:    W_H, positivity/distinctness, special mod-13 types;
              all bounded scalar relation data on every deletion;
destroyed:    transverse character coordinate, determinant ideals,
              exact Fano incidence, polar facets, and band coverage;
needed sidecar:
              a lifted relation sum in Gamma, or its determinant with g;
decisive test:
              the independent shears (34).                           (38)
```

The shears are not claimed to preserve the cover predicate.  Their role is
sharper: they prove that THM-2164's scalar conclusion, by itself, contains no
bound on the determinant/polar coordinates needed to finish (30)--(32).

## 7. Exact hostile controls

The local proof and the five-label majorization received an independent
hostile audit.  The audit specifically found and repaired the zero-determinant
branch: determinant zero is infinite valuation, not positive Archimedean
size, and THM-2080 is required to classify it.

The companion script independently enumerates all translated closed
radius-`1/14` character bands on every invariant-factor group

```text
C_m x C_n,                    m|n, n<=18.             (39)
```

Across `57` nontrivial groups, `3,414` nontrivial characters, and `76,876`
translated band states, it finds `84` covering half-fibre pairs.  Every one
uses the same order-two character and complementary cosets, exactly as the
two-band lemma requires.  Translation probes include every arc endpoint and
one point in every complementary interval, so the subset sweep is exact,
not a phase mesh.  The script also prints a concrete shear family with fixed
evaluated row `(1,2,39)` and determinants drifting from `(13,1)` to
`(13013,2001)`.  Normal and optimized Python transcripts are identical.

Reproduce with

```text
python3 04-computation/lrc14_three_target_second_depth_referee_thm2168.py
```

The finite sweep is a hostile control, not a dependency of the all-group
proof.

## 8. Scope

The theorem narrows the all-independent `(3,5,0)` branch to a two-deep-blocker
invoice, isolates exact daughter locks as its only singleton escape, and
reduces every nontrivial guard-line branch to scalar `4+3` or `5+3`.  It also
names the missing transverse lift for using the new bounded relation rank.
It does not:

- exclude the scalar `4+3` or `5+3` branches;
- give an upper Archimedean bound matching (31);
- lift THM-2164's second scalar relation to a zero character sum; or
- eliminate `(3,5,0)`.

The faithful carrier remains THM-2148's labelled bipartite source/blocker
incidence graph, now with edge type (integral, order-two, or daughter),
13-adic depth, and determinant magnitude.  No intrinsic binary orientation
is present, so a tournament would discard the load-bearing sidecars. QED.
