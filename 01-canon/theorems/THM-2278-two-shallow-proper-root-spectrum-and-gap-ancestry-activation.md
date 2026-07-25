---
id: THM-2278
title: "Two-shallow proper-root spectrum and gap-ancestry activation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. On every strict
  first-depth-one scalar profile, the two THM-2273 shallow-owner flows have a
  common-time labelled root word whose every nonzero C_13 character has
  pointwise squared norm at least (5/21)^5 on the common image.
  Consequently every nonzero global mod-13 residue carries labelled Fourier
  energy at least 712123625/2124823693768776, and the same character is
  active on every deepest-successor safe gap hit by the common image. The
  labels retain the later middle-owner root address, not the earlier shallow
  root digits or one exact integer Fourier atom; no scalar profile is
  excluded.
source: codex-2026-07-25-two-shallow-proper-root-spectrum
depends_on:
  - THM-2269-marked-expiration-root-spectrum-and-branch-state-no-go
  - THM-2273-shallow-owner-flow-and-deep-successor-gap-spread
related:
  - THM-590-apex-7-cyclotomic-gap
  - THM-2276-shallow-owner-residue-aligned-crossing
script:
  - 04-computation/lrc14_two_shallow_proper_root_spectrum_thm2278.py
  - 04-computation/lrc14_two_shallow_proper_root_spectrum_referee_thm2278.py
output:
  - 05-knowledge/results/lrc14_two_shallow_proper_root_spectrum_thm2278.out
  - 05-knowledge/results/lrc14_two_shallow_proper_root_spectrum_referee_thm2278.out
script_sha256:
  - e0b87c60290c9e78c3a7a52af1ee45ce8f6607ed1e5e10382f6ca60f9b4403bd
  - 0fc7e795b2852476477e6eb9713816317d2f82d2044be39ba3dc724284da558a
output_sha256:
  - d07b1c021be7985279c08efca6df9f5958fb7619a85dbe404136182e418ea534
  - 22cefb6e3374f55ae2fc4d7eb268a16264556a95fc4b2f41cf41b7fc1f5ae415
hash_basis: working-tree bytes (LF)
---

# THM-2278 -- two labelled shallow flows fire on every root character

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The marked spectrum in THM-2269 occurs at one owner's expiration. THM-2273
instead transports two shallow exclusive-owner flows to one common time and
keeps them outside a named deepest-successor comb. The missing composition is
to retain both shallow labels at that common time.

The resulting object is not one unlabelled root mask. It is a two-component
root word

```text
(outside-middle-owner ancestry, inside-middle-owner ancestry).          (1)
```

Each active component is a nonempty proper subset of the thirteen roots.
Cyclotomic irreducibility then forces every nonzero character to fire. An
algebraic-norm argument makes the assertion quantitative without enumerating
the `8,190` proper masks.

## 1. The two shallow flows at their common time

Use THM-2273's strict scalar notation

```text
T(x)=13x mod 1,

c_1=13u_1,       c_2=13^b u_2,       c_3=13^c u_3,

2<=b<c,          5<=c<=19,           13 does not divide u_1u_2u_3.       (2)
```

Its disjoint shallow exclusive-owner pieces `E_1,E_2` satisfy

```text
E_1 subset D_(c_1) minus (D_(c_2) union D_(c_3)),

E_2 subset D_(c_2) minus (D_(c_1) union D_(c_3)).                       (3)
```

At time `b`, put

```text
G_1=T^b(E_1),                G_2=T^b(E_2),

Y_j=T(G_j),                  Y=Y_1 union Y_2.                            (4)
```

Exact transport of the middle blocker gives

```text
G_1 subset D_(u_2)^c,        G_2 subset D_(u_2).                         (5)
```

In particular `G_1,G_2` are disjoint even though their terminal images
`Y_1,Y_2` may overlap. Equation (4) is exactly THM-2273's common-time image:

```text
Y=T^(b+1)(E_1 union E_2).                                                (6)
```

That theorem supplies the uniform mass floor

```text
measure(Y)>=Y_0:=5696989/76962600>1/14.                                 (7)
```

It also puts the whole image outside one named deepest successor. With

```text
s=c_3/13^(b+1)=13^(c-b-1)u_3,
```

one has

```text
Y subset D_s^c.                                                          (8)
```

## 2. A rational floor for every character of every proper 13-mask

We first prove the finite cyclotomic lemma used below.

> **Proper-root lemma.** Let `S` be a nonempty proper subset of
> `F_13`, let
>
> ```text
> zeta=exp(2*pi*i/13),
> A_k(S)=sum_(r in S)zeta^(kr).
> ```
>
> Then, for every `k=1,...,12`,
>
> ```text
> |A_k(S)|^2 >= delta_13:=(5/21)^5=3125/4084101.                         (9)
> ```

### Proof

The complement of `S` has the negative of the same nonzero character sum,
because

```text
sum_(r in F_13)zeta^(kr)=0.
```

We may therefore replace `S` by its complement and assume

```text
n=|S| in {1,...,6}.                                                       (10)
```

Put

```text
alpha_l=sum_(r in S)zeta^(lr),
x_l=|alpha_l|^2,                         1<=l<=6.                         (11)
```

No `alpha_l` vanishes. Indeed, its `0/1` coefficient polynomial has degree
at most twelve. If it vanished at a primitive thirteenth root, irreducibility
of

```text
Phi_13(X)=1+X+...+X^12
```

would force the mask to be empty or full.

The cyclotomic norm pairs conjugate factors and gives

```text
product_(l=1)^6 x_l
 =Norm_(Q(zeta)/Q)(sum_(r in S)zeta^r)
 in Z_(>=1).                                                            (12)
```

Finite Parseval gives

```text
2 sum_(l=1)^6 x_l=13n-n^2,

sum_(l=1)^6 x_l=n(13-n)/2<=21.                                          (13)
```

Fix one `x_j`. By AM--GM on the other five factors,

```text
1
 <=product_l x_l
 <=x_j((sum_(l!=j)x_l)/5)^5
 <=x_j(21/5)^5.                                                         (14)
```

Thus `x_j>=(5/21)^5`. Every nonzero `k` is one of the six conjugate pairs,
so (9) follows. QED.

The constant in (9) is a transparent rational floor, not the exact minimum
over proper masks. The proof uses only the prime cyclotomic norm and
Parseval. At a composite root degree, a proper union of cosets can have a
vanishing nonzero character.

## 3. The labelled common-time root word

For `0<=y<1`, label the thirteen roots of `y` by

```text
x_r(y)=(y+r)/13,                         r in F_13.                      (15)
```

Define the two masks and their occupancies:

```text
m_(j,y)(r)=1_(G_j)(x_r(y)),

n_j(y)=sum_r m_(j,y)(r).                                                 (16)
```

A unit danger comb meets every endpoint-free thirteen-root fibre in one or
two sheets. This follows geometrically because its arc length `1/7` lies
strictly between `1/13` and `2/13`; multiplication by the unit `u_2`
permutes the root labels. Combining this fact with (5) gives

```text
y in Y_1:       1<=n_1(y)<=12,

y in Y_2:       1<=n_2(y)<=2.                                           (17)
```

Thus every active mask in (16) is nonempty and proper. For `1<=k<=12`,
put

```text
M_(j,k)(y)=sum_(r=0)^12 m_(j,y)(r)zeta^(-kr),

M_k(y)=(M_(1,k)(y),M_(2,k)(y)).                                         (18)
```

The proper-root lemma applies separately to every active label. THM-2269's
two-root calculation gives the stronger `4/169` floor on the second label.
Pointwise,

```text
|M_(1,k)(y)|^2>=delta_13 1_(Y_1)(y),

|M_(2,k)(y)|^2>=(4/169)1_(Y_2)(y).                                      (19)
```

Since `4/169>delta_13`, equations (7) and (19) imply, for every nonzero
`k`,

```text
|M_k(y)|^2
 =|M_(1,k)(y)|^2+|M_(2,k)(y)|^2
 >=delta_13 1_Y(y),                                                      (20)

integral |M_k(y)|^2dy
 >=delta_13 measure(Y)
 >=712123625/12572921264904.                                            (21)
```

This is the common-time composition: all twelve root characters fire
simultaneously on the full THM-2273 image, while the vector coordinate
remembers which shallow flow supplied the root.

## 4. Global Fourier residue classes, with labels retained

Let

```text
f_j=1_(G_j),

f_hat_j(q)=integral f_j(x)exp(-2*pi*i*qx)dx.                            (22)
```

The root-character calculation from THM-2269 applies componentwise:

```text
integral |M_(j,k)(y)|^2dy
 =169 sum_(q congruent k mod 13)|f_hat_j(q)|^2.                          (23)
```

Summing (23) over the two labels and using (21) gives, separately for every
`k=1,...,12`,

```text
sum_(q congruent k mod 13)
  (|f_hat_1(q)|^2+|f_hat_2(q)|^2)

 >=delta_13 Y_0/169

 =712123625/2124823693768776.                                           (24)
```

Equation (24) is a labelled direct-sum statement. It does **not** assert the
same floor for `1_(G_1 union G_2)=f_1+f_2`, because its two Fourier
coefficients can cancel. The distinction is load-bearing.

## 5. Every forced deepest gap carries every root character

Up to endpoints, the complement in (8) is the disjoint union of the `s`
safe gaps

```text
I_r=[
 (r+1/14)/s,
 (r+13/14)/s
],                              r in Z/sZ.                              (25)
```

Let

```text
R(Y)={r:measure(Y intersection I_r)>0}.                                 (26)
```

Integrating the pointwise inequality (20) on one gap gives, for every
`r in R(Y)` and every `k=1,...,12`,

```text
integral_(I_r)|M_k(y)|^2dy
 >=delta_13 measure(Y intersection I_r)
 >0.                                                                    (27)
```

Thus the *same* set of gap vertices is active for every nonzero root
character. THM-2273's mass/capacity calculation supplies

```text
|R(Y)|
 >=ceil(
   39878923 * 13^(c-b-1) * u_3
   / 461775600
 ).                                                                     (28)
```

In particular, every character is active on at least

```text
c=b+1:       one gap;
c=b+2:       two gaps;
c=b+3:       fifteen gaps;
c=b+4:       190 gaps.                                                   (29)
```

All `135` nonadjacent strict rows force at least two character-active gaps;
the fifteen `b=2` rows force at least fifteen. On an adjacent row, THM-2273
gives the sharper alternative

```text
u_3<=11,

or every nonzero character is active on at least two gaps.               (30)
```

This decorates THM-2273's gap-ancestry graph by all twelve nontrivial
characters. It is naturally a labelled bipartite carrier—shallow source
components versus deepest safe gaps—not a tournament.

## 6. Exact-frequency landing: what is preserved and what is still lost

The vector word (18) preserves:

```text
the E_1/E_2 shallow ancestry label;
the final root sheet from time b to time b+1;
inside-versus-outside membership for the middle unit owner u_2;
the deepest-successor gap index r.                                      (31)
```

It does not preserve:

```text
the earlier c_1/u_1 root digit from time one;
the intervening b-1 root digits or their Hasse transition address;
one particular integer q inside a mod-13 residue class.                 (32)
```

This is the precise interface with the THM-2276 proof candidate. Its
valuation-one crossing is naturally divided at the shallow `c_1/u_1`
address. Equation (24), by contrast, is a sum of squared coefficients over
an entire residue class at the later middle-owner address. Even if the
divided crossing has residue `k`, (24) proves only

```text
some labelled q congruent k mod 13 is energetic,
```

not

```text
f_hat_1(K/13)!=0
or
f_hat_2(K/13)!=0.                                                       (33)
```

The gap activation in (27) does not change this logical direction: a
positive local `L^2` norm is not a prescribed global Fourier coefficient.
Consequently this theorem alone does not shrink THM-2276's hard multiplier
bank.

The cheapest decisive next test is now finite once a coefficient row is
fixed. Refine the exact gap-ancestry graph by the missing `b-1` root digits,
attach to each edge the phase of the proposed divided crossing, and ask
whether the exact coefficient in (33) can vanish while all twelve
inequalities (27) hold. A surviving zero gives the required stopping
witness; infeasibility would be the first lawful exact-frequency landing
lemma.

## 7. Hostile boundaries

Three losses cannot be suppressed.

1. **Do not merge the labels.** Fibrewise it is possible for the first mask
   to occupy all roots outside `D_(u_2)` and the second to occupy all roots
   inside it. Their union is the full thirteen-root mask, whose every
   nonzero character is zero. The vector norm in (20) remains positive.
2. **Do not drop primeness.** Proper masks at composite root degree can be
   character-zero unions of cosets. Irreducibility of `Phi_13` is the exact
   reason (12) is positive.
3. **Do not turn gap count into gap mass.** Equation (27) is positive on
   every hit gap, but THM-2273 gives no uniform positive mass for each
   individual gap. Only their total image mass and number are controlled.

The theorem excludes no scalar profile and proves no successor cover or
blocker return. Its advance is a quantitative, common-time,
ancestry-labelled spectrum on the same deepest safe gaps forced by
THM-2273.

## 8. Exact reproduction

Run

```bash
python3 04-computation/lrc14_two_shallow_proper_root_spectrum_thm2278.py
python3 -O 04-computation/lrc14_two_shallow_proper_root_spectrum_thm2278.py

python3 04-computation/lrc14_two_shallow_proper_root_spectrum_referee_thm2278.py
python3 -O 04-computation/lrc14_two_shallow_proper_root_spectrum_referee_thm2278.py
```

The primary companion checks the norm/Parseval rational certificate,
THM-2273's full `150`-profile mass bank, the exact Fourier floors, and all
gap-count boundaries. The referee independently evaluates the cyclotomic
resultant of all `4,095` complement-reduced masks and constructs every
endpoint-free unit-comb root cell for all twelve nonzero residues modulo
thirteen. Both companions use explicit raising checks, so normal and
optimized modes execute the same audit and produce byte-identical
transcripts. Independent audit separately checked the time indices,
proper-root argument, Fourier factor `169`, gap localization, constants,
and all scope warnings. QED.
