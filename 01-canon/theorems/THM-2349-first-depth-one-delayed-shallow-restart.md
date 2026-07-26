---
id: THM-2349
title: "First-depth-one delayed shallow restart and the repeated two-owner refinement"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. On every one of the
  165 first-depth-one scalar rows, the universal depth-one blocker's
  exclusive-owner set has positive measure and a finite coefficient-
  dependent delay produces a positive literal terminal-word subset.
  On the 15 repeated-first rows the second depth-one owner has the same
  property. More generally, every positive rational step subset of a
  shallow carrier disjoint from the deepest danger comb has, in every
  nonzero shallow root character, a marked c_3-edge whose multiplier is
  coprime to 91. Thus the THM-2327 mixed word/deepest-comb/bare triangle,
  THM-2334 169-twist current, and THM-2343/2344 inverse-correlation
  boundary hold on all 165 rows at a finite row-dependent delayed clock.
  Previously covered positive strict strata retain their stronger
  prescribed clock. No grouped current is proved nonzero off the zero
  target, no scalar row is excluded, and LRC(14) remains open.
source: codex-2026-07-25-repeated-first-twin-fork
depends_on:
  - THM-2138-all-depth-unit-annulus-extremal-law
  - THM-2234-first-depth-one-private-two-owner-mass-and-two-step-expansion
  - THM-2291-repeated-owner-bv-mixing-and-delayed-blocker-handoff
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2323-primitive-fixed-colour-cross-correlation-and-same-gauge-word-alignment
  - THM-2325-prescribed-target-gain-full-lattice-91-unit-needle-bank
  - THM-2326-vertexwise-septimally-primitive-c3-degree
  - THM-2331-two-sided-septimal-address-embedding-in-marked-current
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2343-deep-comb-affine-target-catalyst
  - THM-2344-correlation-inverse-rigidity-and-aligned-tooth-twist-hostile
related:
  - THM-2327-two-colour-marked-unit-c3-triangle
script: 04-computation/lrc14_repeated_first_twin_fork_thm2349.py
output: 05-knowledge/results/lrc14_repeated_first_twin_fork_thm2349.out
script_sha256: 86e17c6f256261a424b6001b2e41f3c3387e80d85843eea10132dbdfbf359d2c
output_sha256: 4985eac48395a05b45b0c57da36162b84fdef22e769129ac0e0de1f23a0d7062
hash_basis: working-tree bytes (LF)
---

# THM-2349 -- every first-depth-one row restarts at a delayed shallow word

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The prescribed-clock route into THM-2327 covered positive strict
shallow-owner strata, but it did not cover every alternative strict
resonance arm and it did not cover the `15` repeated-first rows. The scalar
cover already forces a stronger owner fact which is independent of those
branches:

```text
the universal depth-one exclusive owner has positive measure;
on a repeated-first row, both depth-one owners do.                 (1)
```

Once (1) is retained, ordinary BV mixing supplies a later positive return
for the universal owner, and for either owner on a repeated row. The exact
blocker-word partition can be applied at that later time. THM-2327's
existence proof needs only a positive literal word subset, finite step
complexity, shallow/deep valuation separation, and source/deep
disjointness. Its prescribed clock was used for a uniform height invoice,
not for nonvanishing.

The consequence is a frontier unification:

```text
all 165 rows:
  marked 91-unit triangle at a finite row-dependent delayed clock,
  hence the same 169-twist/inverse-correlation obstruction;

previously covered positive strict strata:
  retain the stronger prescribed controlled clock;

15 repeated-first rows:
  either shallow owner may start the delayed construction.          (2)
```

This does not decrement the scalar ledger. It removes missing marked
incidence as a separate obstruction on repeated and alternative resonance
arms.

## 1. An abstract shallow-carrier triangle

Put

```text
D_n={x in R/Z:||n x||<1/14}.
```

Let `c,d` be positive integers satisfying

```text
nu_13(c)=1<nu_13(d).                                  (3)
```

Let `e,f` be nonzero rational-valued step functions with rational
breakpoints such that

```text
0<=f<=e,

support(e) subset D_c,
support(e) intersection D_d=empty                   (4)
```

up to null endpoints. Let `J_e,J_f` be the numbers of nonzero jumps of
`e,f`. For every prescribed

```text
kappa in F_13^*
```

there are integers `X,Y,m` such that

```text
Y=X+m d,
gcd(m,91)=1,

nu_13(X)=nu_13(Y)=1,
X/13=Y/13=kappa                 mod 13,              (5)

f_hat(X)e_hat(X)e_hat(Y)!=0,                         (6)

f_hat(X)
 (1_(D_d))_hat(m d)
 conjugate(e_hat(Y))!=0.                            (7)
```

One may take

```text
|m|<=13 J_e J_f+7J_e-13.                            (8)
```

The bound is only a convenient finite invoice; no sharpness is asserted.

### Proof: the edge avoiding thirteen

Write

```text
g=gcd(c,d),       a=c/g,       D=d/g.                (9)
```

Equation (3) gives

```text
gcd(a,D)=1,
13 does not divide a,
13 divides D.                                       (10)
```

Apply the Perron operator for multiplication by `g`:

```text
F=P_g f,                   E=P_g e.                 (11)
```

Then `0<=F<=E`, both functions are nonzero rational step functions, and

```text
support(F),support(E) subset D_a.                   (12)
```

Perron transport cannot increase the number of jumps.

Use THM-2323's arithmetic-comb fixed-colour cross-correlation theorem with

```text
N=13D.                                              (13)
```

Choose a unit class `K_0 modulo D` satisfying

```text
(g/13)K_0=kappa                  mod 13,             (14)
```

and put

```text
K_1=K_0+D.                                          (15)
```

Such a class exists by CRT. The prime divisors of `13D` are exactly those
of `D`, so both `K_0,K_1` are primitive modulo `N`. They have the same
nonzero residue modulo thirteen.

THM-2323 gives

```text
0<=h_0,h_1<=J_e J_f-1,

F_hat(K_i+N h_i)E_hat(K_i+N h_i)!=0.                (16)
```

Pulling (16) back through (11), put

```text
A=g(K_0+N h_0),
B=g(K_1+N h_1).                                     (17)
```

Both `A,B` are common `f/e` atoms of shallow grade one and root character
`kappa`. Moreover

```text
B-A=t d,

t=1+13(h_1-h_0),
13 does not divide t,
|t|<=13J_eJ_f-12.                                   (18)
```

Thus `AB` is a marked edge carrying the first CRT colour.

### Proof: the edge avoiding seven and the triangle

Apply THM-2326's general modulated-disjointness lemma to the bare step
function `e`, the marked atom `A`, and the deep comb `D_d`. Equations
(3)--(4) verify its hypotheses. It supplies

```text
C=A+s d,

e_hat(C)!=0,
1<=s<=7J_e-1,
7 does not divide s.                                (19)
```

Because `d/13=0 mod 13`, the shift `s d` preserves the grade-one root
character of `A`; thus `C` has that same root character.

The elementary two-colour triangle says that one of

```text
t, s, s-t
```

is a unit modulo `91`. If it is `t`, use the edge `AB`; if it is `s`,
use `AC`; otherwise use `BC`. At least one endpoint of the selected edge
is `A` or `B`, where both `f` and `e` are nonzero. Orient that marked
endpoint as `X` and the other as `Y`. Equations (5)--(6) follow from
(17)--(19), and (8) follows from the triangle inequality.

Finally,

```text
(1_(D_d))_hat(m d)=sin(pi m/7)/(pi m)!=0            (20)
```

because `7` does not divide `m`. This proves (7) and the abstract lemma.
QED.

## 2. Every row has a positive depth-one owner

Return to any of the `165` first-depth-one scalar rows. Order the three
blockers by valuation:

```text
nu_13(c_1)=1,
1<=nu_13(c_2)<nu_13(c_3).                            (21)
```

The repeated-first rows are exactly those with `nu_13(c_2)=1`. Write

```text
A_0=C_H minus union_(i=1)^5 D_(q_i),

E_j
 =A_0 intersection D_(c_j)
       minus union_(r!=j)D_(c_r),

R_j=A_0 minus D_(c_j).                               (22)
```

We first claim, on every row,

```text
measure(E_1)>0.                                      (23)
```

Suppose, for example, that `E_1` were null. At almost every point of
`A_0` serviced by `c_1`, one of `c_2,c_3` would then also be active.
At every other point of `A_0`, the same conclusion is immediate.
Consequently

```text
C_H subset
 union_(i=1)^5 D_(q_i) union D_(c_2) union D_(c_3)
                                                               a.e. (24)
```

This is a five-unit/two-positive-valuation scalar cover. THM-2138 proves
that the entire such tail is empty. Thus (24) is impossible.

On a repeated-first row, `c_2` also has exact depth one. Swapping `1,2`
in the same deletion argument gives

```text
measure(E_2)>0.                                      (23a)
```

This explains why the tempting repeated-row twin-fork branch is empty. If
both shallow exclusive pieces vanished, then outside `D_(c_3)` the two
shallow danger bits would indeed be forced to agree. But either nullity
already contradicts THM-2138 before that quotient is taken.

On every row THM-2234 applies with `j=1`; on a repeated row it also applies
with `j=2`, because each displayed source has exact depth one:

```text
measure(R_j)>=eta,

eta=2593/90090>0.                                    (25)
```

The owner positivity in (23)--(23a) is qualitative; no useful uniform
lower bound is claimed here.

## 3. BV mixing creates a delayed literal word

Take `j=1` on an arbitrary row, or either `j=1,2` on a repeated-first
row, and put

```text
e_j=measure(E_j)>0.
```

Let `P` be the normalized Perron operator of

```text
T(x)=13x.
```

The boundary of `E_j` lies in the nine scalar comb banks. If

```text
S=H+sum_i q_i+sum_r c_r,
```

then

```text
Var(1_(E_j))<=2S.                                    (26)
```

The BV contraction used in THM-2291 gives

```text
P^k 1_(E_j)>=e_j-2S/13^k                            (27)
```

almost everywhere. Choose any finite

```text
k>=2,
13^k>=4S/e_j.                                       (28)
```

Then (25)--(28) and transfer duality give

```text
measure(E_j intersection T^(-k)R_j)
 =integral_(R_j) P^k1_(E_j)
 >=e_j eta/2
 >0.                                                (29)
```

For the other two blocker labels `a,b`, use THM-2305's exact word
partition

```text
R_j
 =Q_(j,{a}) disjoint_union Q_(j,{b})
                    disjoint_union Q_(j,{a,b})      (30)
```

up to null sets. This identity is independent of the chosen time. Hence
some nonempty word `sigma` satisfies

```text
E_(j,sigma,k)
 :=E_j intersection T^(-k)Q_(j,sigma),

measure(E_(j,sigma,k))>0.                           (31)
```

It is a literal rational step subset of `E_j`, with exact presentation

```text
1_(E_(j,sigma,k))(t)
 =1_(E_j)(t) 1_(Q_(j,sigma))(13^k t).               (32)
```

Thus the transported word dilation is

```text
R=13^k,                    13 divides R.             (33)
```

The jump count in (32) is finite but depends on `k`, and therefore on the
actual coefficient row through `e_j`.

## 4. All-row marked unit triangles

Apply the abstract lemma with

```text
c=c_j,
d=c_3,
e=1_(E_j),
f=1_(E_(j,sigma,k)).                                (34)
```

Exclusivity gives

```text
E_j subset D_(c_j),
E_j intersection D_(c_3)=empty,                    (35)
```

and (21) gives the valuation separation (3). Therefore, for every
nonzero shallow root character `kappa`, there are `X,Y,m` satisfying
(5)--(7). In LRC notation,

```text
(1_(E_(j,sigma,k)))_hat(X)
(1_(D_(c_3)))_hat(m c_3)
conjugate((1_(E_j))_hat(Y))
 !=0,                                               (36)

Y=X+m c_3,
gcd(m,91)=1.                                        (37)
```

This is exactly THM-2327's marked word/deepest-comb/bare triangle, now on
an arbitrary first-depth-one row and at the delayed clock `k`.

On a repeated-first row the theorem gives this construction for **either**
shallow label, not only for the mass-maximizing owner selected by
THM-2291. On every row, what is lost relative to the prescribed-return
strata is their usable uniform multiplier bound.

## 5. Propagation to the current frontier

Equation (32) has the same atomic structure used downstream of THM-2327:
nine present owner factors, nine word factors transported by `R=13^k`,
and nine bare endpoint factors. The only change is the exponent `k`.

The proofs propagate as follows.

1. **THM-2331 termwise address embedding.** Its finite-field lift is
   independent of the word clock, and its atomic application may again
   set every transported-word harmonic to zero. Under the same
   THM-2325 septimal-support hypothesis, every prescribed all-`91`-unit
   target address therefore occurs termwise in the delayed marked
   current.

2. **THM-2334 relation current and target DFT.** Its fixed-address
   summability theorem allows arbitrary integer `R`. In the target
   quotient, (33) makes the transported word target-neutral. The delayed
   current consequently has the same exact `169` coordinate twists and
   Parseval variance criterion.

3. **THM-2343 deep translation.** The comb leg still has the nonzero pure
   deepest-axis charge

   ```text
   p=(0,m) in F_13^2.
   ```

   Hence zero-only full landing is still exactly the inverse-character
   line.

4. **THM-2344 correlation inverse.** All present, word, and bare factors
   remain centered interval indicators or complements, while the word
   shift `R ell_i/13` is integral by (33). The Hermitian law, endpoint
   cross-correlation, shifted-convolution-inverse boundary, off-centre
   covariance test, and aligned-axis hostiles all remain valid.

Thus neither a repeated row nor an alternative strict resonance arm has a
separate missing-incidence obstruction. Every row reaches the same grouped
cancellation problem:

```text
break the shifted endpoint convolution inverse
by two-coordinate or cross-axis structure.                         (38)
```

The conclusion is not a profile exclusion. No off-zero target coefficient,
all-unit grouped aggregate, terminal-component phase, bounded visible
address, or LRC(14) closure is proved. The exact scalar ledger remains
`165`.

## 6. Exact companion

The companion freezes:

- the complete `165=150+15` valuation-profile ledger;
- the complete `c=5,...,19` repeated-profile mass-cap ledger;
- the exact depth-one target floor `2593/90090`;
- the Boolean twin-fork implication which THM-2138 pre-empts;
- all `6,552` CRT two-colour pairs modulo `91`;
- `96` normalized primitive root-pair controls with extra prime factors;
- `65,520` deep-shift grade/root checks.

Reproduce with

```bash
python3 04-computation/lrc14_repeated_first_twin_fork_thm2349.py
python3 -O 04-computation/lrc14_repeated_first_twin_fork_thm2349.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_repeated_first_twin_fork_thm2349.out
```

byte-for-byte after LF normalization.

## 7. Independent audit and scope repair

The independent audit rederived the THM-2138 deletion argument, the
`K_0,K_1` primitive lift with arbitrary extra prime factors, the exact
grade/root calculation, the marked two-colour triangle, and both normal
and optimized companion transcripts. It also caught a material
overstatement in the candidate: existing THM-2327 did not cover every
alternative strict resonance branch. The repaired proof uses the
universal depth-one owner directly on every row, so all `165` reach the
delayed construction without relying on that branch classification.

For completeness, the THM-2326 vertex

```text
C=A+s c_3
```

has the same grade and root as `A` because
`nu_13(c_3)>nu_13(A)=1`; this is the grade-preservation step used by all
three edges of the triangle.

The referee also checked the THM-2234 floor, BV transfer normalization,
Perron support and jump transport, multiplier bound, and every downstream
hypothesis used from THM-2331/2334/2343/2344. It found no conflict with the
fixed-time guardrail in THM-2291 or the corrected `150+15` split in
MISTAKE-251.
