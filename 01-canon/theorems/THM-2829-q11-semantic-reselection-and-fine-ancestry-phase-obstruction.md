---
id: THM-2829
title: "Q11 semantic reselection, Bockstein carry selector, and fine-ancestry obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the fixed
  THM-2806 atom,
  exactly 63 of the 567 semantic cells repair the q=11 six-factor loss
  while retaining the weighted carrier, address/root cylinder, delayed
  semantic values, and both right endpoint scalars.  The first failure is
  the outer THM-2584 U_Q ancestry gate: its nonzero q-support is {7},
  whereas the nonzero E3 support is {3,11}.  Their nonzero strict fibre
  product is empty, but retaining a relative path phase gives an exact positive
  six-point translation-arrow 2-pullback.  The residue arrow lifts through
  0 -> Z/13^5 -> Z/13^6 -> F_13 -> 0 with an ancestry borrow.
  THM-2820's conditional normal label beta(q)=2q uniquely selects q=11;
  its canonical +1 borrow leaves 65098 common sheets, and the resulting
  carry finite difference is the only nonzero response mod 13 among the
  six arrows, with 91-unit magnitude 514.  Canon supplies no physical action
  of that label on U_Q.  The unique
  collapsed rational circulant repair is signed, and no nonzero
  nonnegative circulant repair exists.  No row exclusion or LRC(14)
  conclusion is claimed.
source: root/lrc-q11-semantic-reselection-2026-07-28
depends_on:
  - THM-2584-b-word-depth-five-absolute-deep-root-tensor
  - THM-2782-semantic-arm-right-wing-local-unit-and-endpoint-deck-boundary
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-2813-affine-lift-transvection-and-projective-horn-decoder
  - THM-2820-boolean-idempotent-rigidity-and-norm-top-cotangent-jet-no-go
related:
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2549-future-pullback-target-neutrality-and-cemetery-hall-boundary
  - THM-2555-natural-extension-sheet-charge-and-future-digit-boundary
  - THM-2554-translation-quotient-root-displacement-and-endpoint-swap-parity
  - THM-2697-filtered-affine-handoff-germ-category-and-base-signature-holotopy-boundary
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
  - THM-2792-cyclic-unit-intertwiner-and-positive-naturality-boundary
script: 04-computation/lrc14_q11_semantic_reselection_ancestry_thm2829.py
output: 05-knowledge/results/lrc14_q11_semantic_reselection_ancestry_thm2829.out
script_sha256: 6c29f820bea97a00a5fae3a808b0fb741afed5468e12b289ffae2e9f95c40f98
output_sha256: 8bf43f8570d1fded3fd02257246f465e3adbd3c513ce4302477c13b84aa19712
secondary_script: 04-computation/lrc14_q11_ancestry_e3_path_circulant_thm2829.py
secondary_output: 05-knowledge/results/lrc14_q11_ancestry_e3_path_circulant_thm2829.out
secondary_script_sha256: a7cf624d891b4b3f7147d77c701e60779f16b94d121bf7690bb9eef4e5b1434e
secondary_output_sha256: 56e666b7f5873264a3ed96873d09789a7027ed43134458f64d069f2a01bf433f
hash_basis: LF-normalized bytes
---

# THM-2829 -- q11 is the unique normal path and carry response

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2820 finds a fixed-cell obstruction: the target-only
successor--transvection orbit has one residue retaining all fixed factors
and another retaining both right endpoint origins.  Its own statement
warns that this is not bank-global.  Reselecting the semantic cell repairs
that loss exactly at the endpoint residue `q=11`.

The repair reaches farther than the fixed factors.  It retains the weighted
carrier, the full address/root cylinder, delayed carry data, and both
endpoint scalars.  The first exact failure is instead the outer
depth-five `U_Q` ancestry gate inherited from THM-2584.

That failure has a useful holotopy refinement.  The strict shared-residue
fibre product over the discrete residue set is empty away from the base
point, but the translation-arrow coupling is positive.  Its lift to the
full ancestry/future-digit state carries a nonsplit Bockstein borrow.
Moreover, the normal transvection label already carried by the source
address selects the unique nontrivial path at `q=11`; the borrow response
itself isolates this arrow modulo thirteen.  What remains missing is not
positivity, a path label, or a coefficient-level selector.  It is a
physical natural transformation feeding that label into the `U_Q`
sampling phase without changing the semantic outer word.

## 1. Exact q11 semantic reselection

Use the fixed THM-2806 atom and scales

```text
p=13,
D=p^5=371293,
M=p^6=4826809,

T=297836897838480,
U=T/p=22910530602960,

I=[142004992589460,142005019034340),
J_0=I+431933040,
J_q=J_0+qU.                                             (1)
```

The half-open convention in `(1)` is load-bearing only at finitely many
strict-comb seams.  Every mass and coefficient statement below is an
almost-everywhere statement; no positive endpoint margin is inferred at a
right boundary.

Let the semantic-cell coordinates be `(s,t,c)` in the full
`9*9*7=567` bank used by THM-2806.  Direct exact containment gives

```text
{cells containing I,J_0,J_11 in all six native factors}

  ={0,1,2,3,8,9,10,11,12}
      x {5,6,7,8,9,10,11}
      x {1}.                                             (2)
```

Thus there are exactly

```text
9*7=63                                                   (3)
```

repairing cells.  The lexicographically first witness is

```text
(s,t,c)=(0,5,1).                                        (4)
```

Equations `(2)--(4)` quantify over all semantic reselections of this one
fixed THM-2806 atom.  They do not quantify over all THM-2806 atoms or over
the `165` LRC profiles.

On `(4)`, the clock-one source at `I`, target at `J_0`, and target harmonic
eleven at `J_11` are each one literal weighted interval of weight

```text
w=27581135604.                                          (5)
```

Pulling the harmonic-eleven target back by

```text
431933040+11U                                           (6)
```

returns exactly the source weighted atom at `I`.  This is a
target-carrier harmonic/cross-role translate.  It is not covariance of one
common target under the THM-2813 lift.

## 2. Address, delayed semantic, and endpoint data all survive

Before reselection, the same physical centre has adjacent source and target
address labels

```text
n_s=6716=8 mod13,              n_t=6715=7 mod13.        (7)
```

For the THM-2813 relative affine lift

```text
A_tau(n)=n+tau D(n-7 mod13)             mod M,           (8)
```

put `tau=9`.  Then

```text
A_9(n_s)=3348353,
A_9(n_t)=6715,

n'_t=A_9(n_s)-1=3348352.                                (9)
```

The translated common centre at `J_11` is therefore described by the
adjacent pair

```text
(n'_s,n'_t)=(3348353,3348352)
             =(8,7) mod13.                              (10)
```

Both full address/root cylinders lie in `J_11`.  But `(10)` is not
`(A_9(n_s),A_9(n_t))`; it is the successor--transvection commutator
boundary.  The physical displacement is

```text
7*(A_9(n_s)-n_s)/M
 =7*9/13
 =4+11/13.                                              (11)
```

Thus `q=11` and `tau=9` satisfy

```text
q=7tau,                    tau=7^(-1)q=2q.              (12)
```

The selected interval also retains the exact delayed semantic rows

```text
source carry 12: (0,103478815440),
target carry  6: (0,103478815440).                      (13)
```

At each of the right endpoint origins

```text
(0,0), (12,0),                                           (14)
```

both `q=0` and `q=11` have one piece of mass

```text
26444880                                                (15)
```

and the same two finite-field values

```text
(231164267889491750,630230755085920022).                (16)
```

The endpoint phase exponent is divisible by the endpoint order, so the
equality in `(16)` is exact rather than numerical.  Hence semantic-cell,
carrier-weight, address/root, delayed-semantic, and right-endpoint losses
all occur strictly after the positive repair `(2)`.

## 3. The first obstruction is the outer U_Q ancestry gate

Let `Q` be THM-2584's fixed Boolean word set

```text
Q=build_set(PAT_QB,ZELL),             ZELL=0.            (17)
```

The depth-five inverse-transfer label is

```text
a in Z/DZ.                                              (18)
```

For a normalized source point `x` in `I`, the q-reselected target gate has
the form

```text
G_q(a)
 =1_Q((x+7/M+q/13+a)/D).                               (19)
```

The exact companion enumerates all `D=371293` labels for every
`q in F_13`.  Its target counts are

```text
q=0: 66099,
q=7: 65652,
all other q: 0.                                        (20)
```

Intersecting with the source ancestry sheet gives

```text
q=0: 66099,
q=7: 65612,
all other q: 0.                                        (21)
```

The first common label in both nonempty residues is `a=59162`.  An
independent whole-cylinder audit, rather than a midpoint-only test,
reproduces `(20)--(21)`.

The macro `E3` factor on the same allocation orbit has support

```text
supp(E)={0,3,11}.                                      (22)
```

Therefore

```text
supp(E)\{0}={3,11},
supp(A)\{0}={7},                                       (23)
```

where

```text
A=66099 delta_0+65652 delta_7.                         (24)
```

The nonzero strict shared-q fibre product is empty.  In particular, the
positive physical q11 packet has no target `U_Q` label at all.  Changing
the semantic cell, clock, delayed carry, or endpoint origin cannot repair
this earlier gate.

Nor does a static replacement by either of the other canonical outer
words repair it while preserving the word.  On the same q11 cell, an exact
whole-cylinder scan gives

```text
word       source labels       q11 target labels       common labels
QA              0                    10787                    0
QB          66099                        0                    0
QAB         10786                        0                    0.
```

The nearest positive static cospan is directed:

```text
QB(source) -> QA(target):
  449 whole-cylinder labels,
  first 59306, last 311961,
  every label =156 mod169.
```

It preserves the q11 local cell, delayed rows, and endpoint scalars, but
changes the semantic word `{b}` to `{a}`.  Taking `QA union QB` recovers
the same 449 labels only by forgetting which word survived.  Thus the
remaining sidecar cannot be a different fixed outer mask; it must include
a lawful word-transition arrow or preserve the original word.
The 449-label cospan has the same formal directed shape as THM-2461's
atom-to-word matrix; that theorem explains why it cannot be replaced by
same-time Boolean idempotence, but does not identify this cospan with its
prescribed-clock current or supply this particular edge.

The congruence `156=0+13*12` is a useful hostile against digit
misidentification, not an address decomposition.  The label `a` is a
depth-five inverse-branch coordinate.  In THM-2584's exact chart

```text
X_(v,a)=(y+v+13a)/(13D),
```

the predecessor carry is `floor(13{DX})=v` and the deep digit is
`2v+epsilon`; both are independent of `a`.  Concretely, the canonical
label `a=59162` has `a mod169=12` and source carry twelve, while the first
cross-word label `a=59306` has `a mod169=156` but simultaneously sees
source carry twelve and target carry six.

There is a second exact alias hostile.  `PAT_QA/PAT_QB` indices `7,8` are
the fixed speed slots

```text
W_7=13^3,                  W_8=2*13^5,
```

not address residues.  At `q=0` the actual source/target address edge is
still `8 -> 7`, but all `66099` common sheets retain `QB -> QB`, not
`QB -> QA`.  Hence neither matching pair of numerals supplies a torsor
map.  The first missing law is an address-edge-to-semantic-word natural
transformation.

## 4. The translation-arrow 2-pullback is positive

The empty nonzero strict fibre product in `(23)` is not the end of the
structure.
Define the translation action groupoid

```text
T=F_13 acting on F_13,

objects: q in F_13,
arrow (q,h): q -> q+h.                                 (25)
```

Because the translation action is free and transitive, there is exactly one
arrow between every ordered pair of objects.  Regard `supp(E)` and
`supp(A)` as discrete groupoids over `T`.  Their groupoid 2-pullback, or
equivalently the corresponding comma/path object, has objects

```text
(q,h) with q in supp(E), q+h in supp(A).
```

Give such an arrow the product weight

```text
F_A(q,h)=E(q) A(q+h).                                  (26)
```

This explicit groupoid is what "path" means below.  Without it, `(26)`
would be only a weighted correspondence and must not be called the
homotopy pullback over the discrete set `F_13`.

The weighted arrow support, retaining `h in F_13`, is exactly

```text
(q,h;weight)=

(0,0;66099),   (0,7;65652),
(3,10;66099),  (3,4;65652),
(11,2;66099),  (11,9;65652).                           (27)
```

The two affine sections are

```text
h=-q,                         q+h=0,
h=7-q,                        q+h=7.                   (28)
```

The marginals are

```text
sum_h F_A(q,h)=(66099+65652)E(q),

sum_(q+h=r) F_A(q,h)=3A(r).                            (29)
```

Replacing the target coefficient `65652` by the fixed-representative common
coefficient `65612` gives the equally positive residue-marginal coupling

```text
F_C(q,h)=E(q) C(q+h),

C=66099 delta_0+65612 delta_7.                         (30)
```

In particular

```text
F_C(11,9)=65612>0.                                     (31)
```

Equation `(31)` uses the same depth-five representative `a` on both sides.
It is an exact coefficient correspondence, not yet the natural-extension
lift of the arrow.  Wrapped arrows change that representative.

Thus the strict pullback over the discrete residue set has only the base
point `q=0`.  Its nonzero part is empty because it permits only identity
arrows `h=0`, not because the translation groupoid has no positive ancestry
arrow.

### 4.1. The natural lift has a Bockstein borrow

Package the depth-five ancestry label and the future residue as

```text
n=13a+q in Z/13^6 Z.
```

Adding the arrow phase `h`, represented in `{0,...,12}`, gives

```text
n+h
 =13(a+kappa(q,h))+(q+h mod13),

kappa(q,h)=floor((q+h)/13).                            (31a)
```

The sign in `(31a)` agrees with the `+h` numerator convention below.  The
six arrows in `(27)`, in their displayed order, therefore have canonical
lifted common-sheet weights

```text
(66099,65612,65579,65612,65579,65098),                 (31b)
```

whereas suppressing the borrow gives

```text
(66099,65612,66099,65612,66099,65612).                 (31c)
```

The lifted-minus-same-section response is

```text
(0,0,-520,0,-520,-514).                                (31d)
```

Modulo thirteen this becomes

```text
(0,0,0,0,0,6).                                        (31e)
```

Thus the carry finite difference is supported only on the
beta-selected arrow `(q,h)=(11,9)` modulo thirteen.  Its absolute response
is

```text
514,                    gcd(514,91)=1.                 (31f)
```

In particular, the natural q11/h9 lift is still positive, with `65098`
common sheets.  The reverse orientation `a -> a-1` would instead give
`65619`; this proves that orientation and borrow data are load-bearing.
The exact companion checks `(31b)--(31f)` at the left endpoint, midpoint,
and right endpoint minus one of the whole interval.

This is the nonsplit extension

```text
0 -> Z/13^5 Z -> Z/13^6 Z -> F_13 -> 0.
```

Thirteen lifted generator steps return to the same future residue but move
`a -> a+1`.  Hence the residue translation groupoid has no sheet-preserving
`F_13` splitting: the carry is its Bockstein wall.  The unit response in
`(31f)` is the first exact `U_Q` coefficient diagnostic of that wall.
It remains an auxiliary relative-factor count, not a physical frequency,
canonical current, or proved Bockstein intertwiner.

## 5. The normal decoder selects q11 uniquely

Equation `(12)` supplies the conditional THM-2820 normal label

```text
beta(q)=tau=2q.                                        (32)
```

On the nontrivial ancestry section in `(28)`, compatibility asks

```text
beta(q)=7-q.                                           (33)
```

Therefore

```text
2q=7-q,
3q=7,
q=11,
beta(11)=9.                                            (34)
```

Among the nonzero E3 residues `{3,11}`, q11 is uniquely selected:

```text
q=3  -> q+beta(q)=9,       not an ancestry residue;
q=11 -> q+beta(q)=7,       the nontrivial ancestry residue. (35)
```

This is structural rather than a scan accident.  More generally, if a
successor commutator has

```text
q=u tau
```

and a slope-one fine connection uses `h=tau`, then landing on ancestry
residue `r` requires

```text
(u+1)tau=r,

tau=(u+1)^(-1)r,
q=u(u+1)^(-1)r,                                        (36)
```

whenever `u+1` is a unit.  Here `(p,u,r)=(13,7,7)`, so `(36)` gives
`(tau,q)=(9,11)`.

Equations `(25)--(36)` construct a positive weighted 2-pullback in the
explicit translation action groupoid and a canonical candidate path label
on the retained source address.  On the natural lift, that label necessarily
uses the `+1` ancestry borrow and has the positive weight `65098` from
`(31b)`.  Its mod-thirteen finite difference is the unique response in
`(31e)`.  These facts still do not construct a physical action of that
label on `U_Q`.

## 6. A9 does not already supply the missing phase

The distinction is exact.  The physical q11 movement in `(11)` has already
used `A_9`.  It moves the residue-eight source and fixes the residue-seven
target:

```text
A_9(6716)=3348353,               A_9(6715)=6715.        (37)
```

The reselected target address `3348352` is obtained by restoring
adjacency after moving the source.  It is not the image of the old target
under `A_9`.  Thus the q11 table is already the full
successor--`A_9` commutator boundary, not a table awaiting another
automatic `A_9` action.

The API in `(17)--(19)` has no address or transvection argument.  Adding

```text
G_(q,h)(a)
 =1_Q((x+7/M+(q+h)/13+a)/D)                            (38)
```

is a new target-relative action.  The unreduced numerator in `(38)` is
load-bearing: when `q+h>=13`, it includes the `+1` borrow in `(31a)`,
not the same-section table `(30)`.  At `q=11` it has the conditional values

```text
h=2: q+h=13, target/common counts 66099/65579;
h=9: q+h=20, target/common counts 65652/65098.          (39)
```

But neither value is produced by `(37)`.

The two natural attempts to hide `(38)` in an existing gauge both fail.

First, the THM-2806 representative change obeys

```text
build_set(PAT_QB,ell+sW)(y)
 =build_set(PAT_QB,ell)(y+s/13).                        (40)
```

Inside the depth-five inverse transfer, `(40)` only permutes

```text
a -> a+s13^4                                           (41)
```

and leaves `q` fixed.  Hence the q-support in `(20)` is
representative-gauge invariant.

Second, the THM-2813 correction is divisible by `D`, so

```text
A_tau(a)=a                         mod D.               (42)
```

Placing the q11 physical shift directly on the ancestor variable similarly
reindexes

```text
a -> a+11*13^4=a+314171                              (43)
```

rather than creating `h=9`.

Finally, THM-2584's deep-root pushforward makes the sheet `a` disappear
because its `C/D=2` factor is integral.  There is no retained output slot
on which the normal label could act after that marginalization.

The formal one-column carry wall suggested by THM-2788 has the same useful
role asymmetry as `(37)`: a depth-six analogue would move residue eight
and fix residue seven.  THM-2788 proves no factor/current intertwiner from
that wall to `(38)`, and its actual theorem is at the modular address
level.  It is shape inspiration, not a supplied map.

In the language of THM-2697, the distinction is a base-signature
distinction.  Equations `(40)--(43)` are vertical fibre decorations: they
leave the ancestry base residue `q` unchanged.  Equation `(38)` changes the
induced base map `q -> q+h`.  A successful repair must therefore add a new
normal-phase arrow or 2-cell, not another coordinate on the old empty
base.

THM-2611 identifies the exact conditional consumer.  If the address-normal
labels and the fine `U_Q` phases were supplied as principal `C_13` fibres
and one equivariant fixed-action identification between them were physical,
then a single basepoint would determine the whole phase map; `(34)` would
select its q11 value.  The positive arrow in `(31)` supplies neither that
principal ancestry action nor the equivariant identification.  Thus the
remaining debt is existence and provenance of one bibundle arrow, not a
further holonomy calculation.

THM-2549 identifies the scale of that arrow exactly.  Its old target action
on a depth-five inverse-transfer sheet changes

```text
a -> a-theta*13^4,
```

which is the coarse reindexing in `(41)` and leaves the q-spectrum fixed.
By contrast, its genuinely future target action

```text
z -> z-theta/13
```

pulls back at depth five to a shift by `-theta/13^6`.  With
`theta=-h`, this is precisely the extra base movement in `(38)`.
THM-2549's carry-corrected unit-role formula `b_w=head root` is therefore
the closest positive algebraic sidecar, but that theorem explicitly does
not identify the old-sheet root with a genuinely target-active future root.
THM-2829 isolates the same missing identification on one exact q11 packet.

THM-2555 makes the digit horizon and the borrow precise.  Here
`a in Z/13^5 Z` is the ancestry sheet, while `q` is the first depth-six
future digit.  The old action changes the top ancestry digit; `q -> q+h`
is the new future action, and `(31a)` is its compulsory sheet borrow.
Consequently the positive count `65098` is evidence for the exact
natural-extension response, but not evidence that the LRC packet carries
the required future action.  One must still identify that future digit
with the retained normal label on the same semantic current.

Likewise, in THM-2554's translation quotient, the arrow phase here is the
displacement

```text
d=(q+h)-q=h.
```

The table `(26)` is a weighted displacement-arrow table.  THM-2554 proves
that a nonzero displacement is not semantically forced without a common
physical root torsor and a genuinely target-active field.  Equation `(34)`
selects `d=9` conditionally; it does not supply those hypotheses.

## 7. Collapsing the path forces signs

The path fibre in Section 4 cannot be collapsed by a nonnegative circulant
map.  Work in

```text
Q[C_13]=Q[z]/(z^13-1)
```

and put

```text
a=66099,               b=65652,

A(z)=a+bz^7,
E(z)=1+z^3+z^11.                                       (44)
```

Since seven is a unit modulo thirteen,

```text
B(z)=sum_(j=0)^12 (-b)^j a^(12-j) z^(7j)
```

satisfies

```text
A(z)B(z)=a^13+b^13=:Delta,                             (45)
```

where

```text
Delta=
880705707102199130841131933042872776034368258206185041400058691.
                                                               (46)
```

Thus `A` is a unit and the unique rational circulant coefficient map taking
`A` to `E` is

```text
K(z)=E(z)B(z)/Delta.                                   (47)
```

In residue order `0,...,12`, the numerator coefficients in `(47)` have
the exact sign pattern

```text
(+,+,+,+,+,+,+,-,-,-,-,-,-)                          (48)
```

and

```text
augmentation(K)=1/43917.                               (49)
```

The full integer numerator is recorded in the exact secondary output.

There is also a support-only proof that signs are unavoidable.  Any
positive kernel coefficient at shift `s` contributes the support pair

```text
{s,s+7}.                                               (50)
```

No pair in `(50)` lies inside

```text
{0,3,11}.                                              (51)
```

Therefore no nonzero nonnegative circulant kernel can map the ancestry
profile into the E3 support.  The positive path coupling `(26)` escapes
this no-go precisely because it retains `h`; after `h` is forgotten, the
unique linear circulant repair must cancel.

The carry selector `(31d)` does not evade this positivity boundary.  It is
the signed finite difference between two positive count tables.  Its
91-unit value at `(11,9)` is therefore a sharp detector of the missing
lift, not itself a nonnegative transport kernel.

This is the same naturality boundary isolated abstractly by THM-2792:
invertibility over the group algebra does not imply a positive physical
intertwiner.

## 8. Connection contract and exact scope

The proof-complete candidate mechanism is:

```text
source:
  one fixed THM-2806 atom, its 567 semantic-cell bank, the THM-2813
  source normal address, and the THM-2584 outer U_Q ancestry gate;

target:
  the q=11 full-factor/endpoint carrier and the six-point
  translation-arrow ancestry coupling;

physical map already present:
  source A_9 address -> q=11 successor commutator carrier,
  with semantic cell reselected to (0,5,1);

preserved:
  all six native factors, weighted carrier mass, address/root cylinder,
  delayed carry values, two right endpoint origins and scalars, and the
  source normal label tau=9;

destroyed by strict q matching:
  every nonzero common U_Q ancestry label;

strongest positive survivor:
  the same-section residue coupling F_C(11,9)=65612, the canonical
  borrow-aware lift of weight 65098, and the unique normal-path equation
  beta(11)=9; its carry response is 6 mod13 with 91-unit magnitude 514;

missing sidecar:
  a lawful target-relative action beta -> h on the pre-marginal U_Q
  sheet, including its ancestry borrow and compatible with the physical
  common cospan and the original QB semantic word;

cheapest decisive test:
  construct one physical map realizing the +1-borrow G_(11,9) in (38)
  on the selected q11 atom without replacing its target, endpoint origin,
  semantic word, or ancestry label by a separately chosen object.         (52)
```

The theorem does not identify `J_11` with `A_9` of the original common
target, construct a physical common beta, retain `a` after THM-2584's
pushforward, turn the auxiliary unit carry response `(31f)` into a physical
current, provide the QB-to-QA word transition, or turn the signed map
`(47)` into a positive current.
It gives no canonical owner-word completion, row exclusion, ledger
decrement, or LRC(14) conclusion.

## 9. Exact companions

Run

```text
python 04-computation/lrc14_q11_semantic_reselection_ancestry_thm2829.py
python -O 04-computation/lrc14_q11_semantic_reselection_ancestry_thm2829.py

python 04-computation/lrc14_q11_ancestry_e3_path_circulant_thm2829.py
python -O 04-computation/lrc14_q11_ancestry_e3_path_circulant_thm2829.py
```

Normal and optimized modes byte-match, respectively,

```text
05-knowledge/results/lrc14_q11_semantic_reselection_ancestry_thm2829.out
05-knowledge/results/lrc14_q11_ancestry_e3_path_circulant_thm2829.out.
```

The primary companion pins the THM-2806/2782 physical constructors,
enumerates all 567 semantic cells, checks the selected weighted carrier,
address/root and delayed-semantic data, replays both endpoint fields, and
enumerates every one of the `13^5` outer ancestry labels at all thirteen
residues.  It also checks the borrow-aware six-arrow weights and their
whole-cylinder endpoint stability, the unique mod-thirteen carry response,
and the three canonical outer words plus the 449-label QB-to-QA cospan.

The secondary companion proves the exact group-ring inverse and sign
pattern, exhausts all translated support pairs, constructs both six-point
residue path couplings and their marginals, verifies the unique normal
section, checks the nonsplit thirteen-step deck shift and carry selector,
and checks the address-role, coarse-gauge, and depth-reduction boundaries.

An independent immutable audit replayed both companions in normal and
optimized modes against their stored outputs, rederived the fixed
`63/567` scope, all ancestry and outer-word counts, the natural borrow,
the unique modular carry response, the normal selector, the address/gauge
boundaries, and the signed-circulant no-go.  It also independently checked
whole-cylinder rather than midpoint-only support.  The audit's one typing
repair—separating ancestry `a`, delayed carry/root, factor slots, and
address residues—is installed above.

**QED.**
