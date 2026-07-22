---
id: THM-2053
title: "Rank-two relation planes have explicit geodesic and transverse-grid LRC terminals"
status: >
  PROVED from the standing cited LRC for at most twelve distinct nonzero
  speeds. The torus-geodesic estimate and the adjacent-normalized-column
  construction are elementary. Every rational two-plane in Q^13 containing a
  primitive positive row has torus margin at least 1/13. The resulting exact
  sufficient gate is max_i|a z_i-b u_i|<=(a^2+b^2)/91; failure of this gate is
  only uncertified, not unsafe. An intrinsic transverse
  coordinate also gives M(v)>=1/13-R/(2N), where N divides an exposed pair
  sum, and the exact residue deck D_N(m) gives a finite rational phase-height
  conductor independent of all longitudinal coefficients. This makes every
  THM-2052 plane finite in parameter space, but does not enumerate all residual
  templates or prove LRC(14).
source: codex-2026-07-21-DC2-LRC14-termination
depends_on:
  - THM-2052
  - LRCUpTo13
  - "T. Sungkawichai and T. Trakulthongchai, Eleven, twelve, and thirteen lonely runners, arXiv:2604.23906v1 (preprint)"
related:
  - THM-1002-pair-sum-denominator-bound-and-the-bounded-gap-case
  - THM-2055
  - HYP-4342
  - HYP-4346
  - HYP-3456
  - HYP-8846
script: 04-computation/lrc14_transverse_residue_deck_referee_codex_20260721.py
output: 05-knowledge/results/lrc14_transverse_residue_deck_referee_codex_20260721.out
script_sha256: bbe1d4f81fb6a1269488cfbd99c9391800ecbeb26ba64c2ac2345c135f442b2
output_sha256: 26339f905a40615cfdb51548078cf7269c054ab25688e8e7c189573354561bd7
---

# THM-2053 -- large directions in every relation plane

## 1. Parameter planes and their torus

Let `U` be a two-dimensional rational subspace of `Q^n`, and choose a
`Z`-basis `u,z` of the saturated lattice `U intersect Z^n`. Write

```text
c_i=(u_i,z_i) in Z^2,
F_U(x,y)=min_i ||u_i x+z_i y||,
M_T(U)=max_((x,y) in (R/Z)^2) F_U(x,y),                 (1)
L(U)=max_i sqrt(u_i^2+z_i^2).
```

For a primitive nonzero parameter direction `d=(a,b) in Z^2`, put

```text
v(d)_i=a u_i+b z_i,
M(d)=max_(t in R/Z) min_i ||v(d)_i t||.                 (2)
```

Thus `M(d)` is the ordinary one-dimensional lonely-runner maximin of the row
selected by `d`, while `M_T(U)` permits the two parameters to move
independently.

## 2. General primitive-geodesic rate

For every primitive `d=(a,b)`, define

```text
E_U(d)=max_i |a z_i-b u_i|/(2(a^2+b^2)).               (3)
```

Then

```text
M_T(U)-E_U(d) <= M(d) <= M_T(U),                       (4)
E_U(d)<=L(U)/(2 sqrt(a^2+b^2)).                        (5)
```

This is the exact all-direction version of HYP-4342's `(1,N)` rate. The first
bound retains the oriented area between the target direction and each column;
the second is its round Lipschitz envelope.

### Proof

Let

```text
C_d={(at,bt) mod Z^2:t in R/Z}.
```

The right inequality in (4) is immediate because `C_d` is a subset of the
two-torus and `F_U(at,bt)` is exactly the objective in (2).

For the left inequality, put `n=(b,-a)`. Since `gcd(a,b)=1`, the subcircle
`C_d` is the kernel of the primitive character

```text
(x,y) |-> bx-ay mod Z.                                  (6)
```

Given any torus point `X`, choose a lift to `R^2` and an integer `m` with

```text
|n dot X-m|<=1/2.
```

Moving `X` in the normal direction by

```text
-(n dot X-m)n/||n||_2^2
```

lands in the kernel (6). Hence every torus point is within Euclidean distance

```text
1/(2||n||_2)=1/(2 sqrt(a^2+b^2)).                      (7)
```

of `C_d`.

More precisely, the normal displacement is `-s n/||n||^2` with `|s|<=1/2`.
Its change in the `i`th character is at most

```text
|c_i dot n|/(2||n||^2)
 = |b u_i-a z_i|/(2(a^2+b^2)).                         (8)
```

Apply (8) to a maximizer of `F_U` and take the maximum over `i`; this proves
the left inequality in (4). Cauchy--Schwarz gives (5). QED.

The proof uses only elementary torus geometry. It does not invoke LRC and it
does not assume that the coordinates `v(d)_i` are positive or distinct.

## 3. A checkable generic floor

Call `U` **full-support-repeat exposed** if some primitive direction
`d_0=(a_0,b_0)` has

```text
p_i=a_0u_i+b_0z_i !=0 for every i,
|p_i|=|p_j| for some i!=j.                              (9)
```

For `n=13`, the standing cited LRC through thirteen total runners implies

```text
M_T(U)>=1/13.                                          (10)
```

Indeed, the thirteen nonzero integers in (9) have at most twelve distinct
absolute values. Apply the cited theorem to that set of at most twelve speeds.
It supplies a time `t_0` at which every distinct speed has distance at least
`1/(k+1)>=1/13`. Repetitions and sign changes do not change circle distance,
so

```text
F_U(a_0t_0,b_0t_0)>=1/13,
```

which proves (10).

Combining (4) and (10) gives the exact anisotropic sufficient gate delivered
by this estimate:

```text
max_i |a z_i-b u_i| <= (a^2+b^2)/91
                         ==> M(d)>=1/14.               (11)
```

Indeed, (11) says `E_U(d)<=1/182`, and

```text
1/13-1/14=1/182
```

so equality is allowed: LRC needs a weak `1/14` witness. By (5), the convenient
round corollary is

```text
sqrt(a^2+b^2)>=91 L(U)  ==>  M(d)>=1/14.              (12)
```

There is an exact geometric form of the parameters not certified by (11).
Put `Jc_i=(z_i,-u_i)`. Failure of (11) is equivalent to membership in one of
the `2n` open disks

```text
||d-(91 sigma/2)Jc_i||^2 < (91^2/4)||c_i||^2,
                  i=1,...,n,  sigma in {+1,-1}.        (13)
```

Indeed, after choosing `sigma` to match the sign of `Jc_i dot d`, the
inequality is just `||d||^2<91 sigma Jc_i dot d` completed to a square. Every
circle in (13) has the origin on its boundary. Thus the exact uncertified
carrier for this gate is a union of tangent disks in the integer parameter
lattice, not the whole round envelope from (12). Membership in this carrier
does not imply that the corresponding row is unsafe.

## 4. Adjacent normalized columns expose the repeat direction

The word "generic" can be removed constructively, without an external
direction theorem. Assume `U` contains a primitive row `v` with every
coordinate positive. Let

```text
L=sum_i Z c_i
```

be the rank-two lattice generated by the coefficient columns. Evaluation at
the actual row is a homomorphism

```text
phi:L->Z,       phi(c_i)=v_i.                          (14)
```

It is primitive because the `c_i` generate `L` and `gcd_i v_i=1`.
In the coordinates of Section 1 one in fact has `L=Z^2`: the row lattice
`U intersect Z^n` is saturated, so Smith normal form says that the gcd of the
`2 x 2` column minors is one, equivalently the columns generate the full
coefficient lattice. We retain `L` to make the transverse construction
intrinsic.

Normalize the columns into the affine line `phi=1`:

```text
p_i=c_i/v_i in L tensor R.                            (15)
```

There are at least two distinct `p_i`, since the columns span a rank-two
space. Order the distinct points on this affine line and choose representatives
`p_i,p_j` of two adjacent values. Put

```text
w=c_i+c_j.
```

Then

```text
w/(v_i+v_j)=(v_i p_i+v_j p_j)/(v_i+v_j).              (16)
```

lies strictly between the two adjacent normalized columns. Hence `w` is
parallel to no `c_k`: positive proportionality would make the normalized
points equal, and negative proportionality is impossible because both
`phi(w)` and `phi(c_k)` are positive. In particular `c_i,c_j` are independent.

A primitive covector in `Hom(L,Z)` perpendicular to `w` is therefore nonzero
on every `c_k`, while its values on `c_i,c_j` are opposite. In a lattice basis
of `L` this is an integer parameter direction, so it is exactly the
full-support repeat projection (9). Thus (10), and consequently the exact and
round terminals (11)--(12), hold on every genuine positive two-plane. Duplicate
normalized columns cause no problem: adjacency is taken among distinct values.

The only external input is the standing cited LRC result for twelve distinct
nonzero speeds. In particular, no planar direction-count theorem is needed.

## 5. Exact transverse-grid terminal

The adjacent pair gives a sharper, intrinsic terminal. Write

```text
w=g w_0,
```

where `g>=1` and `w_0` is primitive in `L`, and extend `w_0` to a lattice
basis `(w_0,eta)`. Write

```text
c_k=a_k w_0+m_k eta.
```

Since no column is parallel to `w`, every `m_k` is nonzero. From
`c_i+c_j=g w_0` one has

```text
m_i=-m_j=:b!=0.                                       (17)
```

Moreover `gcd_k m_k=1`: the images of the generating columns `c_k` generate
`L/Zw_0`, which is infinite cyclic with coordinate `m_k`.

On the transverse circle

```text
K={X in Hom(L,R/Z):X(w_0)=0},
```

put `theta=X(eta)`. The restricted objective is

```text
F_K(theta)=min_k ||m_k theta||.
```

All speeds `m_k` are nonzero and the pair in (17) repeats one absolute value,
so there are at most twelve distinct absolute speeds. Settled LRC for at most
twelve speeds gives

```text
max_theta F_K(theta)>=1/13.                            (18)
```

Put

```text
N=phi(w_0)>0,       M=phi(eta),       R=max_k |m_k|.
```

Primitivity of `phi` gives `gcd(N,M)=1`. The actual orbit `t phi` meets `K`
exactly at `t=ell/N`, and its `eta`-coordinates are

```text
ell M/N       (ell mod N).
```

They form the complete `N`-grid because `M` is a unit modulo `N`. The function
`F_K` is `R`-Lipschitz. Writing `M(v)=max_t min_k||v_kt||`, nearest-grid
rounding in (18) proves

```text
M(v)>=1/13-R/(2N).                                    (19)
```

Consequently

```text
N>=91R  ==>  M(v)>=1/14.                              (20)
```

If `N>91R`, the witness is strict; reducing `ell/N` then gives
`q_exit(v)<=N`. At equality only the weak LRC witness is asserted.

The complementary region is finite for each coefficient template. Choose the
sign of `eta` so `b>0`, and write

```text
c_i=A w_0+b eta,       c_j=(g-A)w_0-b eta.
```

Positivity of the two exposed speeds forces

```text
-AN/b<M<(g-A)N/b.
```

Thus `N<91R` leaves finitely many positive integers `N`, and for each `N`
only finitely many coprime integers `M`. It also carries exact modular data:

```text
v_k=a_kN+m_kM == m_kM (mod N),
v_i+v_j=gN.                                           (21)
```

Hence `N` is a divisor of an exposed pair sum, exactly the ruler type from
[THM-1002, the pair-sum denominator bound](THM-1002-pair-sum-denominator-bound-and-the-bounded-gap-case.md),
and the complete denominator-`N` row is determined by the transverse template
`(m_k)` up to multiplication by the unit `M`.

## 6. Exact residue deck and phase-height conductor

The actual information at the denominator `N` is sharper than the Lipschitz
bound (19). For an integer `x`, let `|x|_N` be the least absolute residue
modulo `N`, and define the transverse residue deck

```text
D_N(m)=max_(ell mod N) min_k |ell m_k|_N/N.           (22a)
```

At the actual time `t=ell/N`, equation (21) gives

```text
||v_k t||=|ell m_k M|_N/N.
```

Since `M` is a unit modulo `N`, multiplication by `M` permutes the complete
residue grid. Therefore

```text
max_(ell mod N) min_k ||v_k ell/N||=D_N(m),
M(v)>=D_N(m).                                         (22b)
```

This quantity is independent of every longitudinal coefficient `a_k` and of
the unit `M`. A hypothetical counterexample must satisfy

```text
D_N(m)<1/14.                                          (22c)
```

It also has exact divisor monotonicity. If `d|N`, the denominator-`d` grid is
a subset of the denominator-`N` grid, so

```text
D_N(m)>=D_d(m).                                       (22d)
```

Thus the bad moduli form a down-set under divisibility. In particular, if
`d|N`, `d<=14`, and every `m_k` is nonzero modulo `d`, then the single grid
point `1/d` gives

```text
D_N(m)>=D_d(m)>=1/d>=1/14,                            (22e)
```

and the entire `N`-fiber is safe.

There is an exact finite conductor for every template. Put

```text
G(m)={theta in T:min_k||m_k theta||>=1/14},
lambda(m)=maximum length of a connected component of G(m).
```

The set is nonempty by (18). Every closed circle arc of length at least `1/N`
meets the complete `N`-grid. Hence

```text
D_N(m)<1/14  ==>  N<1/lambda(m).                      (22f)
```

All boundary points of `G(m)` belong to the finite rational wall set

```text
theta=(14r+1)/(14|m_k|) or (14r-1)/(14|m_k|) mod 1.   (22g)
```

Consequently `lambda(m)` is exactly computable by rational interval
arithmetic, as is the sharp phase-height conductor

```text
H(m)=1+max({N>=1:D_N(m)<1/14} union {0}).             (22h)
```

Every `N>=H(m)` is safe. No monotonicity in the usual order of integers is
claimed: for the AP template below, `D_13=1/13` but `D_15=1/15`. Divisibility
monotonicity (22d) is the valid statement.

The conductor need not be found by scanning every `N`. Split each circular
component of `G(m)` into a closed interval `[alpha,beta]` with rational
endpoints. The `N`-grid meets it exactly when

```text
ceil(N alpha)<=floor(N beta).                         (22h')
```

After putting the endpoints over a common denominator `Q` and writing
`N=Qh+s`, both sides of (22h') are linear in `h` plus a correction depending
only on `s mod Q`. Thus the bad moduli are an exact finite union of initial
segments in residue classes modulo `Q`. This imports HYP-3456's AP84
floor-count operation into the transverse atlas while preserving the sampled
loneliness predicate; endpoint owners must be retained separately if later
Euler transport needs them.

The component-length argument also improves the uniform constant whenever the
template has fewer distinct absolute speeds. Let

```text
s=|{|m_k|:1<=k<=13}|<=12.
```

Settled LRC for these `s` speeds gives a transverse maximizer of height at
least `1/(s+1)`. Since the objective is `R`-Lipschitz, its `1/14` superlevel
component has length at least

```text
lambda(m)>=(13-s)/(7(s+1)R),
D_N(m)<1/14  ==>  N<7(s+1)R/(13-s).                  (22i)
```

For `s=12` this recovers `N<91R`; already for `s<=11` it improves the terminal
to `N<42R`.

### Exact arithmetic-progression deck

More generally, for every `k,N>=1`,

```text
D_N({1,...,k})=floor(N/(k+1))/N.                      (22j)
```

Write `N=(k+1)q+r`, `0<=r<=k`. For any `ell`, the `k+1` residues
`0,ell,...,k ell` have a circular gap at most `floor(N/(k+1))=q`; the
difference of its endpoints is `j ell` for some `1<=|j|<=k`, proving the
upper bound. Taking `ell=q` attains it: for `1<=j<=k`, both `jq>=q` and
`N-jq>=q+r`.

At target height `1/(k+2)`, this deck is good exactly when `q>=r`. Thus the
sharp eventual conductor is

```text
N>=k(k+1),                                            (22k)
```

and `N=k(k+1)-1` is bad. For the thirteen-entry transverse template
`m=(1,-1,2,...,12)`, this gives

```text
D_N(m)=floor(N/13)/N.
```

Its largest bad modulus at height `1/14` is `155=13*11+12`, and every
`N>=156` is safe. The exact conductor `156` is far smaller than the generic
Lipschitz cutoff `91R=1092`.

For `k=11`, formula (22j) also proves the AP grid rungs previously found by
`04-computation/lrc_gridmax_window_search_opus_S403.py`: values such as
`2/27`, `3/40`, `3/41`, and the subsequent `g/N` windows are forced by
`floor(N/12)/N`, not merely discovered by search. This is a closed-form
explanation of that scout, not a new global LRC implication.

## 7. Consequence for the THM-2052 atlas

THM-2052 places every hypothetical primitive counterexample in the kernel of
an eleven-dimensional bounded relation code. Each such kernel has dimension
at most two, and only finitely many kernels occur.

For a one-dimensional kernel, its saturated integer lattice has only one
primitive generator up to sign, so that atlas cell is already a single exact
row. For a two-dimensional kernel, Section 4 supplies the repeat direction;
choose a saturated basis. A primitive positive speed row has a primitive parameter
direction: if `(a,b)=g(a',b')` with `g>1`, every speed is divisible by `g`.
The residual not certified by (11) consists only of primitive parameters
violating that sufficient gate. By (12), it is contained in the finite disk

```text
a^2+b^2 < (91L(U))^2.                                 (23)
```

Section 6 first discards every entire fiber with
`D_N(m)>=1/14` and leaves only the finite bad-modulus down-set inside
`N<91R`, together with the exact pair-sum and residue data (21). Inside those
bad fibers, THM-2055 deletes nonvertex determinant owners from the determinant
test and intersects each normal-fan sector with one tangent disk from (13).
The corresponding non-hull columns remain in the residue deck and in the LRC
row. Every remaining row can then be decided exactly by the linked pair-sum
denominator theorem.
Therefore the entire THM-2052 atlas has an explicit finite parameter-space
terminal **without harvesting a twelfth bounded relation**.

THM-763 already gives a global but enormous finite speed ceiling
`sum_i v_i<=91^12`; this theorem does not claim first decidability. Its new
content is a plane-intrinsic uniform safe tail and the structured tangent-disk
residual (13), obtained without enumerating every speed row below that ceiling.

What remains is computational/structural rather than infinitary: construct or
compress the finite atlas, choose controlled saturated bases, and discharge
the bad residue fibers intersected with THM-2055's owner-labelled tangent-disk
sectors. This theorem does not claim those lattice points have been enumerated
and does not prove LRC(14).

### Calibration on the exact one-tail plane

Take `u=(1,2,...,13)` and `z=e_12`, so `v(1,b)` replaces speed `12` by
`w=12+b`. For `b!=0`,

```text
max_i |z_i-bu_i|=13|b|.
```

Consequently (11) becomes

```text
b^2-1183|b|+1>=0,
```

which holds for every integer `|b|>=1183`. The gate therefore reduces this
infinite line to a finite interval. HYP-2896 then does much better on that
interval and beyond: its divisibility walls `12|w` and `14|w` give `q=12`,
`q=14`, or an explicit affine binding phase. This is the intended division of
labor: geodesic geometry removes the projective tail, while a resonance fan
discharges the finite arithmetic core.

## 8. Coefficient-height control on the transverse template

Specialize to a two-plane `U=W_(Q,3)(v)^perp` from THM-2052, with
`Q=91^6`. For each `k` outside the exposed pair `{i,j}`, the guaranteed
height-`Q` relation on the triple `{i,j,k}` is a column identity

```text
A_k c_i+B_k c_j+C_k c_k=0,
|A_k|,|B_k|,|C_k|<=Q,       C_k!=0.                   (24)
```

The last inequality follows from independence of `c_i,c_j`. Put
`L_ij=Zc_i+Zc_j` and `D=[L:L_ij]`. The eleven cosets `[c_k]` generate
`L/L_ij`, while (24) says that `[c_k]` is killed by `C_k`. Therefore

```text
D<=product_(k outside {i,j}) |C_k|<=Q^11.             (25)
```

In the basis `(w_0,eta)`, the determinant of `c_i,c_j` is `-gb`, so

```text
D=g|b|.                                               (26)
```

Reducing (24) modulo `Zw_0` gives

```text
(A_k-B_k)b+C_km_k=0,
```

and hence `|m_k|<=2Q|b|`. Thus

```text
R<=2Q|b|,       gR<=2QD<=2Q^12.                      (27)
```

For a still-unresolved row, `N<91R`. Combining this with (21), (24), and
(27) yields

```text
v_i+v_j=gN<182Q^12,
max_k v_k<182Q^13.                                   (28)
```

The numerical box (28) is much weaker than THM-763's existing `Q^2` total-
height bound. Its value is structural: the finite coordinate is the exposed
pair-sum divisor `N`, and the entire denominator-`N` residue row is the fixed
transverse template `(m_k)` up to a unit.

## 9. Quantifier repair and assumption challenge

HYP-4346 says that if infinitely many independent directions share finitely
many scale-only templates, then the plane collapses. Applied naively, it only
produces *some* escaping direction and has the wrong quantifier for a specified
target row. Equations (4), (11), and (12) are the repair: the repeat direction
plus the settled lower-dimensional theorem raises the two-torus to `1/13`, and
then every sufficiently long **specified** target direction inherits the
`1/14` bound with an explicit loss.

This also challenges the earlier assumption that the only terminal above
rank eleven is a twelfth bounded relation. Rank twelve is still a valid
maximal-minor terminal, but the torus-geodesic region is already a terminal at
rank eleven. The next decisive problem is not to enumerate the whole round
envelope (23). First use the exact bad decks `D_N(m)<1/14`; then use THM-2055's
support norm, rational normal fan, and tangent sectors (13) inside those fibers.
The surviving carrier must retain the bounded triple-code shape, pair-sum
divisor `N`, modular template `(m_k mod N)`, and active hull owner.

An oriented-matroid quotient alone is too coarse. Take columns

```text
c_r=(1,r)  (0<=r<=11),       c_K=(1,K),       K>11.
```

They have the same rank-two chirotope and order type for every `K>11`, and the
column lattice is saturated because `det(c_0,c_1)=1`. The covector
`phi(x,y)=x+y` gives the positive distinct row `1,2,...,12,K+1`, and its
normalized columns retain `c_0,c_1` as an adjacent pair. For that pair,
`w=(2,1)` and in the unimodular basis
`((2,1),(1,1))` the transverse coordinates are

```text
m_r=2r-1,       m_K=2K-1.
```

Thus `R` is unbounded while the oriented matroid is fixed. Any valid atlas
compression must retain Plucker/determinant magnitudes or an equivalent height
sidecar, not only signs and incidences. This is consistent with the repo's
MISTAKE-185 warning against the old HNF shortcut.

## 10. Exact referee

The stored referee is an arithmetic regression, not an independent proof of
the general saturation, transverse-grid, or conductor arguments above. It
checks the AP deck law for every `1<=k<=12` and `1<=N<=300`, divisor
monotonicity on four hostile templates, every available small-divisor exit,
the nonmonotone pair `D_13,D_15`, the largest AP bad modulus `155` within that
range, and the fixed-chirotope/unbounded-`R` family. The proof of (22j)--(22k),
not the finite scan alone, establishes the all-modulus conductor `156`. The
referee passes both

```text
python 04-computation/lrc14_transverse_residue_deck_referee_codex_20260721.py
python -O 04-computation/lrc14_transverse_residue_deck_referee_codex_20260721.py
```

and byte-matches the stored output. The frozen SHA-256 hashes are recorded in
the frontmatter.
