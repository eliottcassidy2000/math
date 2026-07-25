---
id: THM-2269
title: "Marked expiration root spectrum and branch-state no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENT EXACT REFEREE. For any exclusive-owner
  flow at its expiration, each occupied terminal fibre has one or two marked
  predecessor sheets. Every nonzero C_13 Fourier mode of that rooted mask is
  pointwise nonzero and has squared magnitude at least 4/169. Consequently
  every nonzero mod-13 residue class carries global Fourier energy at least
  4*image_mass/28561, while the total nonconstant fibre energy is at least
  twelve times the expiration-image mass. THM-2263 makes the every-residue
  floor at least 15041431/501746795550 on every strict profile and
  5229541/501746795550 on every repeated-first profile. A four-bit proposed
  successor state is not exact: marked ancestry and current selected-blocker
  service are distinct bits. A rational branch satisfying every forward
  pointwise cover clause still has permanent zero marked-core future, so
  sibling-sheet ancestry/gluing remains necessary. A marked residue class
  containing THM-2145's integer cross-cut frequency need not contain that
  frequency in the marked support; no profile is excluded.
source: codex-2026-07-25-marked-root-spectrum
depends_on:
  - THM-2255-valuation-separated-pair-cap-and-exclusive-owner-mass
  - THM-2263-thirteen-adic-gap-pair-spectrum-and-profile-sharp-owner-floor
related:
  - THM-2261-expiration-image-surjectivity-and-one-core-carrier-no-go
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2144-anisotropic-selberg-kraft-relation-box
  - THM-2145-two-block-spectral-crossing-and-6-plus-7-carry
  - THM-2166-hybrid-core-smoothing-low-carry-crossing
  - THM-2267-static-owner-coverage-is-flag-and-transition-holonomy-is-a-cut-kernel
  - THM-2268-two-shell-private-owner-trident-and-raw-carry-cocycle-no-go
  - THM-2270-simultaneous-balanced-cut-relation-and-six-uniform-orientation
  - THM-2275-mixed-scalar-relation-and-guard-blocker-crossing
script:
  - 04-computation/lrc14_marked_expiration_root_spectrum_thm2269.py
  - 04-computation/lrc14_marked_expiration_root_spectrum_referee_thm2269.py
output:
  - 05-knowledge/results/lrc14_marked_expiration_root_spectrum_thm2269.out
  - 05-knowledge/results/lrc14_marked_expiration_root_spectrum_referee_thm2269.out
script_sha256:
  - 605f96b1429aa77b1562f2f66d923b5512469a9c746100d98f10519670424ed2
  - d43a2ce8358faf4734ba3729942ccd2d86af2f63fda8ac6ab30e96ed3346c776
output_sha256:
  - 4f8e2b564d052e4ec451e9177d2e5a5b8c4ee0fb89af581278fcd2510a2c53a6
  - bdf1d4c488869b36f76757cabca2aab833383c345bedca3aedf5b9c58d19bea3
hash_basis: working-tree bytes (LF)
---

# THM-2269 -- the marked predecessor has a nonconstant root spectrum

Put

```text
T(x)=13x mod 1,
D_s={x in R/Z:||sx||<1/14}.                          (1)
```

Let `E` be a measurable exclusive-owner stratum with labelled blocker

```text
c=13^lambda u,             lambda>=1,       13 does not divide u. (2)
```

Define the last pre-expiration set and its image by

```text
F=T^lambda(E),
B=T(F).                                               (3)
```

Exact owner transport gives

```text
F subset D_u.                                        (4)
```

THM-2255 and THM-2263 give a labelled `E` of this form in every remaining
scalar profile, with a quantitative lower bound on `measure(B)`.

## 1. The rooted marked-predecessor word

Represent `y in R/Z` by `0<=y<1` and label its thirteen roots

```text
x_r(y)=(y+r)/13,                      r in F_13.      (5)
```

The marked root mask and occupancy are

```text
m_y(r)=1_F(x_r(y)),
n(y)=sum_(r in F_13)m_y(r).                          (6)
```

Then

```text
B={y:n(y)>0}.                                        (7)
```

Because `F subset D_u` and a unit danger comb occupies at most two roots of
every thirteen-root fibre,

```text
n(y) in {1,2}                  for y in B.           (8)
```

This is the first finite object that survives THM-2261's image
surjectivity. The unmarked image bit only records `n>0`; the rooted word
records which one or two owner-predecessor sheets actually reached `y`.

## 2. Every nonzero root character fires

Let

```text
zeta=exp(2*pi*i/13)
```

and use the unnormalized finite Fourier transform

```text
M_k(y)=sum_(r=0)^12 m_y(r) zeta^(-kr),
0<=k<13.                                             (9)
```

If `n(y)=1`, every `M_k(y)` is a root of unity and is nonzero. If
`n(y)=2`, with marked sheets `a!=b`, then

```text
M_k(y)=zeta^(-ka)(1+zeta^(-k(b-a))).                 (10)
```

For `k!=0`, equation (10) cannot vanish. Otherwise a thirteenth root of
unity would equal `-1`; squaring would give a nontrivial thirteenth root
whose square is one, impossible because thirteen is odd and prime.
Therefore

```text
y in B,  1<=k<=12       implies       M_k(y)!=0.     (11)
```

There is also a uniform quantitative form of (11). For a two-sheet mask put
`d=k(b-a) mod 13`, so `d` is one of `1,...,12`. Then

```text
|M_k(y)|^2=|1+zeta^d|^2
          =4 cos^2(pi*d/13)
          >=4 sin^2(pi/26)
          >=4/169.                                  (11a)
```

The first minimum occurs at `d=6,7`. The last inequality is the elementary
chord bound

```text
sin(t)>=2t/pi,                0<=t<=pi/2,            (11b)
```

at `t=pi/26`. A one-sheet mask has `|M_k|^2=1`, so the same lower bound holds
there. Thus, separately for **every** nonzero root frequency,

```text
integral |M_k(y)|^2 dy >=(4/169)measure(B),
                                      1<=k<=12.      (11c)
```

Finite Parseval gives the exact pointwise energy

```text
sum_(k=0)^12 |M_k(y)|^2=13n(y).                     (12)
```

Since `M_0=n`, the nonconstant energy is

```text
sum_(k=1)^12 |M_k(y)|^2
 =13n(y)-n(y)^2
 =12                         if n(y)=1,
 =22                         if n(y)=2.              (13)
```

In particular,

```text
sum_(k=1)^12 integral |M_k(y)|^2 dy
 >=12 measure(B).                                   (14)
```

By pigeonhole, one fixed nonzero root frequency satisfies

```text
integral |M_k(y)|^2 dy>=measure(B).                  (15)
```

This is an exact positive root-sheet channel. The elementary two-variable
LP behind (14) is sharp: if `a` and `b` are the measures of fibres with
occupancy one and two, respectively, then

```text
measure(B)=a+b,
total nonconstant energy=12a+22b
                            =12(a+b)+10b.            (16)
```

With only an image-mass floor, the minimum occurs at `b=0`.

## 3. Translation to global Fourier residue classes

Let `f=1_F` and use circle Fourier coefficients

```text
f_hat(q)=integral f(x)exp(-2*pi*i*qx)dx.             (17)
```

For a trigonometric polynomial, substituting (5) into (9) and summing the
root characters gives

```text
M_k(y)
 =13 exp(2*pi*i*k*y/13)
    sum_(h in Z) f_hat(k+13h)exp(2*pi*i*h*y).        (18)
```

The identity extends to every `f in L^2(R/Z)` by density. Parseval on the
`y` circle now gives

```text
integral |M_k(y)|^2dy
 =169 sum_(q congruent k mod 13)|f_hat(q)|^2.        (19)
```

Equations (14)--(15) become

```text
sum_(13 does not divide q)|f_hat(q)|^2
 >=(12/169)measure(B),                               (20)

for some k in {1,...,12},

sum_(q congruent k mod 13)|f_hat(q)|^2
 >=measure(B)/169.                                   (21)
```

More strongly, (11c) and (19) give, for every `k=1,...,12`,

```text
sum_(q congruent k mod 13)|f_hat(q)|^2
 >=(4/28561)measure(B).                              (21a)
```

Thus the marked pre-expiration set cannot be spectrally supported only on
multiples of thirteen, or even omit any nonzero residue class. This is
stronger than saying its terminal image is large: every nonzero
ramification residue class carries positive quadratic mass.

## 4. Profile-sharp quantitative floors

THM-2263 improves THM-2255's pair ledger and expiration bounds. In every one
of the `150` strict rows it supplies a labelled flow with

```text
measure(B)>=L_strict:=15041431/70270200.             (22)
```

Hence (20)--(21) give

```text
sum_(13 does not divide q)|f_hat(q)|^2
 >=15041431/989638650,                               (23)

some k!=0:
sum_(q congruent k mod 13)|f_hat(q)|^2
 >=15041431/11875663800.                             (24)
```

In fact (21a) gives the following floor in **each** of the twelve nonzero
residue classes:

```text
for every k!=0:
sum_(q congruent k mod 13)|f_hat(q)|^2
 >=15041431/501746795550.                            (24a)
```

For every one of the `15` repeated-first rows, the corresponding image floor
is

```text
L_repeat:=5229541/70270200,                          (25)
```

and therefore

```text
sum_(13 does not divide q)|f_hat(q)|^2
 >=5229541/989638650,                                (26)

some k!=0:
sum_(q congruent k mod 13)|f_hat(q)|^2
 >=5229541/11875663800.                              (27)
```

The corresponding every-residue floor is

```text
for every k!=0:
sum_(q congruent k mod 13)|f_hat(q)|^2
 >=5229541/501746795550.                             (27a)
```

The strict image floor also exceeds the one-core two-window mass from
THM-2261 by

```text
15041431/70270200-6055/28561
 =24332839/11875663800>0.                            (28)
```

Equation (28) is now a visibly nontrivial gap, but THM-2261 proves that the
one-core set is not a target. Equations (23)--(24), rather than that invalid
support inclusion, are the surviving quantitative information.

## 5. The proposed four-bit successor state misses a bit

At a terminal phase `y`, a tempting rooted state assigns each sheet only

```text
(marked ancestry, two other owner bits, guard/unit literal).  (29)
```

This is not an exact full-cover state. The marked bit `m_y(r)` says that the
sheet lies in `F=T^lambda(E)` and therefore descended from the selected
owner. The scalar cover at the **current** phase instead uses the three
original blocker bits

```text
s_y(r)=1_(D_c)(x_r(y)),
o_(1,y)(r), o_(2,y)(r).                              (30)
```

After division through the owner depth, `m=1` does not imply `s=1`.
Consequently a direct labelled-coordinate cover state must retain the five
named predicates

```text
(m,s,o_1,o_2,ell),                                   (31)

ell=1_(C_H minus union_i D_(q_i)),

ell<=s or o_1 or o_2.                                (32)
```

Identifying `m` with `s` erases exactly the expiration event. Moreover,
iterating (31) exactly cannot retain `ell` as one bit: THM-2201 shows that
the literal alone is not autonomous under the labelled updates. THM-2201
proves exact Hasse-jet incidence updates once the transition address is
supplied; it does not by itself supply the missing full transition theorem.
A conservative faithful candidate therefore retains the ancestry bit, all
three current owner bits, all six labelled guard/unit sheet words, their
rooted Hasse jets, the transition address, and sibling gluing. Any compression
must prove that it preserves the same incidence updates and cover data.

## 6. Exact hostile branch, including every forward cover clause

The distinction in Section 5 is realized uniformly on every strict profile.
For

```text
2<=b<c
```

take

```text
H=1,
x=79/338,

(q_1,...,q_5)=(2,339,677,1015,1353),

(c_1,c_2,c_3)=(13,13^b,13^c).                       (33)
```

At time zero,

```text
x in C_1,
x notin union_i D_(q_i),
x in D_(c_1),
x notin D_(c_2) union D_(c_3).                      (34)
```

Indeed `2x=79/169` is safe, while each of the other unit speeds is
`1+338i` and sends `x` to `x+79i`. The blocker calculations are those of
THM-2261.

The forward orbit is

```text
x_0=79/338,
x_1=T(x_0)=1/26,
x_t=1/2                         for every t>=2.      (35)
```

Define the literal and current blocker bits along this branch by

```text
ell_t
 =1_(C_H)(x_t) product_i(1-1_(D_(q_i))(x_t)),

O_(j,t)=1_(D_(c_j))(x_t).                            (36)
```

They satisfy

```text
t=0:        ell=1,       (O_1,O_2,O_3)=(1,0,0);
t=1:        ell=0,       (O_1,O_2,O_3)=(0,0,0);
t>=2:       ell=0,       (O_1,O_2,O_3)=(0,0,0).     (37)
```

At `t=1` the guard is off because `||1/26||<1/7`. At `t>=2`, the guard is
on but `q_1=2` is dangerous at `1/2`, so the unit literal is off. Therefore
the exact pointwise cover clause

```text
ell_t<=O_(1,t) or O_(2,t) or O_(3,t)                 (38)
```

holds at every time on the forward branch.

Nevertheless the selected normalized-core word after expiration is

```text
1_(D_1)(x_t)=0                    for every t>=2.    (39)
```

At the immediate marked predecessor `x_1`, the divided owner bit
`1_(D_1)(x_1)` is one, while the current selected service bit
`1_(D_13)(x_1)` is zero. This is the promised `(m,s)=(1,0)` witness.

Thus any LP whose atoms are single forward branch words and whose only
global-cover constraint is (38) admits this Dirac branch as a feasible point
while retaining the marked-to-zero obstruction (39). The witness does **not**
satisfy or assert the scalar cover on all sibling roots or all circle phases.
It proves the exact stopping boundary: a full-cover successor argument must
glue sibling sheet words and their ancestry, not merely enforce the cover
down the marked orbit.

## 7. Relation and cut sidecars: the missing map

THM-2145 proves that every labelled `6+7` cut of a strict LRC(14)
counterexample has a shared nonzero integer Fourier frequency `K`, with
bounded coefficient representations on both sides; THM-2166 supplies its
low-carry form. THM-2144/2164 supply complementary bounded rank, and
THM-2270 supplies one relation crossing all balanced cuts simultaneously.
If `13` does not divide `K`, equation (21a) now guarantees marked energy
somewhere in the residue class of `K`; there is no remaining
choice-of-residue mismatch. But that still does not land at the particular
integer frequency `K`:

```text
THM-2269:
  every nonzero residue k mod 13 supports quadratic Fourier mass of 1_F;

THM-2145:
  some nonzero integer K is a common frequency of two safe-block
  approximants.                                      (40)
```

Neither theorem proves

```text
f_hat(K)!=0,
or that the cut coordinates map to the selected scalar owner ancestry.   (41)
```

If `13` divides `K`, one must first descend to
`K/13^nu_13(K)` at the correct root depth. In either case the missing
composition theorem must retain this valuation address and show that a
bounded cross-cut frequency lands in the marked spectral support rather than
merely in an energetic residue class.

THM-2275 makes this valuation issue unavoidable on the quotient-faithful
scalar cut. Its crossing frequency has the form

```text
K_scalar=a_H H+sum_i a_i q_i
        =-sum_j b_j c_j !=0.                         (41a)
```

Every scalar blocker `c_j` is divisible by `13`, so

```text
13 divides K_scalar.                                 (41b)
```

Thus that particularly well-labelled crossing certificate always begins in
the ramified zero residue. Its exact interface with (21a) is the normalized
integer

```text
K_scalar/13^nu_13(K_scalar) mod 13,                  (41c)
```

but THM-2275 does not transport this descended frequency to `1_F`. This
identifies the missing map more sharply: preserve the scalar
guard/unit-versus-blocker labels while descending the frequency through the
same selected-owner root address.

THM-2267 gives the complementary transition target: once exact transported
pieces and eligibility sets are glued, a binary owner cut can carry positive
min-cut energy. THM-2268 globally forces all three private labels and at
least three owner switches, but does not put those torsion witnesses in the
Haar-positive ancestry flow `F`. Equations (23)--(27) provide positive
Haar-scale marked energy; THM-2267/2268 provide owner-cut and switch
structure. The missing bridge is a common ancestry/cut gluing which makes
them act on the same pieces.

## 8. Scope and reproduction

This theorem proves an every-residue marked root-sheet Fourier channel and
the named-coordinate/forward-branch no-go. It does not prove a nonzero
owner-switch tax, align a THM-2145/2271 relation with the marked spectrum, or
exclude any scalar profile.

Run

```bash
python3 04-computation/lrc14_marked_expiration_root_spectrum_thm2269.py
python3 -O 04-computation/lrc14_marked_expiration_root_spectrum_thm2269.py

python3 04-computation/lrc14_marked_expiration_root_spectrum_referee_thm2269.py
python3 -O 04-computation/lrc14_marked_expiration_root_spectrum_referee_thm2269.py
```

The primary companion checks every one- and two-sheet mask, the sharp
two-variable energy LP, the angle certificate behind the uniform
every-character bound, all fractions, and the hostile branch on all `150`
strict profiles. The independent referee rederives finite Parseval by
character orthogonality, independently checks the angle certificate and
global residue-class transform on exact rational trigonometric polynomials,
and reconstructs the witness from congruences. Both normal/optimized
transcript pairs are byte-identical. QED.
