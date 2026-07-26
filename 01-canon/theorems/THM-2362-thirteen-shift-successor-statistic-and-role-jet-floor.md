---
id: THM-2362
title: "Thirteen-shift successor statistic and role-jet floor"
status: >
  PROVED + VERIFIED-EXACT + CORRECTED AFTER INDEPENDENT HOSTILE AUDIT.
  centered danger d=1_(||x||<1/14), target-coordinate shifts obey
  sum_s d(y+s/13)=2-d(13y), not 2-d(y). Thus every nonnegative weight
  supported in d has a nontrivial shifted-danger Fourier mode with real
  part at least 11rho/156 and nonzero-mode energy at least
  121rho^2/2028; a weight supported in the complement has corresponding
  floors rho/156 and rho^2/2028. The exact extra coordinate is the
  successor overlap rho_+=int w(y)d(13y). Inverse-root probes
  d((y+s)/13) instead have an unordered count controlled by d(y), while
  their nonzero-mode sum also needs a chosen sheet anchor. On a pure
  THM-2305 word the named danger factor is redundant by the scalar cover,
  so its actual first-jet role has the danger floor. Fork dangers and
  complement roles are not redundant. This proves a positive role-mass
  jet, not survival of the complex marked current, a nonzero relation
  target, a scalar-row exclusion, or LRC(14).
source: codex-2026-07-25-thirteen-shift-successor-statistic
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
related:
  - THM-2232-same-core-signed-eigen-markov-dual-exclusion
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2356-finite-field-chirp-gram-tomography-and-bockstein-pairing
  - THM-2364-anchored-corner-forces-mixed-deep-blocker-colour
script: 04-computation/thirteen_shift_successor_role_jet_thm2362.py
output: 05-knowledge/results/thirteen_shift_successor_role_jet_thm2362.out
script_sha256: a7309b8ec09556f5dac1e0a6eb1afdf45fb01a76395adb06d69c8396c6dbc81c
output_sha256: 5856d80d6ddf2cf5b1383380ff54f3e64ec5d3d83d3047a953abeba61953b8d6
hash_basis: working-tree bytes (LF)
---

# THM-2362 -- a target shift sees the successor bit

**PROVED + VERIFIED-EXACT + CORRECTED AFTER INDEPENDENT HOSTILE AUDIT.**

The thirteen-root count has two similar-looking but inequivalent forms.
The distinction is a scale:

```text
target-coordinate factor shift:
  d(y+s/13) sees d(13y);

inverse-root probe:
  d((y+s)/13) sees d(y).                            (1)
```

Conflating them loses the successor/expiration bit. Keeping that bit gives
uniform nonzero Fourier-mode floors and a correctly typed pure-word jet.

## 1. Exact pointwise shift counts

Put

```text
d(y)=1_(||y||<1/14),

g(y)=1-d(y),                         y in T,

s in F_13.                                           (2)
```

Away from the finite set of interval endpoints,

```text
sum_(s in F_13)d(y+s/13)=2-d(13y),

sum_(s in F_13)g(y+s/13)=11+d(13y).                 (3)
```

Indeed, apply the exact thirteen-root count

```text
sum_s d((Y+s)/13)=2-d(Y)
```

at `Y=13y`. The complement identity follows by subtracting from thirteen.
The `26` open cells in one period are the common refinement cut by either
side of (3), so endpoint conventions affect only a null set.

Let `w>=0` be integrable and put

```text
rho=int_T w(y)dy,

rho_+=int_T w(y)d(13y)dy.                           (4)
```

The second quantity is the exact successor overlap. It is not determined
by the current-scale support bit `d(y)`.

## 2. A danger role has a uniform nontrivial mode

Assume

```text
support(w) subset D={d=1},

rho>0.                                              (5)
```

Define the real shifted-role masses and normalized finite transform

```text
M_s=int w(y)d(y+s/13)dy,

Mtilde(k)=1/13 sum_s M_s zeta^(-ks),

zeta=exp(2*pi*i/13).                                (6)
```

Then `M_0=rho`, while (3)--(4) give

```text
Mtilde(0)=(2rho-rho_+)/13.                          (7)
```

Finite Fourier inversion at `s=0` therefore gives

```text
sum_(k!=0)Mtilde(k)
 =M_0-Mtilde(0)
 =(11rho+rho_+)/13.                                (8)
```

The right side is positive real. Among the twelve nonzero characters,
some `k` satisfies

```text
Re Mtilde(k)
 >=(11rho+rho_+)/156
 >=11rho/156.                                      (9)
```

Cauchy--Schwarz also gives

```text
sum_(k!=0)|Mtilde(k)|^2
 >=(11rho+rho_+)^2/2028
 >=121rho^2/2028.                                  (10)
```

Every nonzero character is primitive because thirteen is prime. Equality
in the last bounds requires `rho_+=0` and all twelve nonzero modes to be
equal as complex numbers. Since their sum is positive real, their common
value is positive real. These are the sharp algebraic boundaries.

## 3. A complement role has its own floor

Assume instead

```text
support(w) subset D^c,

rho>0,                                              (11)
```

and put

```text
N_s=int w(y)g(y+s/13)dy.
```

Now `N_0=rho`, and (3) gives

```text
Ntilde(0)=(11rho+rho_+)/13,

sum_(k!=0)Ntilde(k)=(2rho-rho_+)/13.                (12)
```

Since `0<=rho_+<=rho`,

```text
some k!=0 has
  Re Ntilde(k)>=(2rho-rho_+)/156>=rho/156,          (13)

sum_(k!=0)|Ntilde(k)|^2
 >=(2rho-rho_+)^2/2028
 >=rho^2/2028.                                     (14)
```

Here the uniform boundary occurs at `rho_+=rho`; simultaneous equality in
the Cauchy floor also requires all twelve nonzero modes to be the same
positive real number. Thus a current-scale danger or complement role
always moves under its genuine `s/13` factor translation, but its exact
response remembers the successor bit.

## 4. Why the inverse-root constants do not transfer

For the different probe

```text
d_s^root(y)=d((y+s)/13),
```

the almost-everywhere pointwise laws are

```text
sum_s d_s^root(y)=2-d(y),

sum_s (1-d_s^root(y))=11+d(y).                     (15)
```

The strict-open endpoints are a finite null exception. Formula (15)
determines the zero-character average, but not the sum of the nonzero
characters. Choose a measurable real lift `ytilde` of the circle
coordinate and put

```text
rho_D=int w(y)d(y)dy,

R_s=int w(y)d((ytilde+s)/13)dy,

A=R_0.
```

Then

```text
Rtilde(0)
 =(2rho-rho_D)/13,

sum_(k!=0)Rtilde(k)
 =A-(2rho-rho_D)/13.                               (15a)
```

For the complementary root profile

```text
C_s=int w(y)(1-d((ytilde+s)/13))dy,
```

the exact formula is

```text
sum_(k!=0)Ctilde(k)
 =(2rho-rho_D)/13-A.                               (15b)
```

The current support mass `rho_D` does not determine the chosen
root-sheet anchor `A`. Under the centered lift, `s=0` is a danger anchor
on `D` but a zero complement anchor on `D^c`; `s=6` is a universal
complement-safe anchor, at the cost of rephasing the Fourier inversion.

With the standard `[0,1)` lift, two exact pointwise profiles expose the
missing coordinate:

```text
y=99/100:
  d(y)=1,
  sum_s d((y+s)/13)=1,
  d(y/13)=0,
  sum_(k!=0)Rtilde(k)=-1/13
  for unit point mass;                              (15c)

y=1/2:
  g(y)=1,
  sum_s g((y+s)/13)=11,
  g(y/13)=0,
  sum_(k!=0)Ctilde(k)=-11/13
  for unit point mass.                              (15d)
```

Each profile is constant on a sufficiently small open interval, so
positive integrable weights give the same scaled hostile. Reindexing the
root sheets changes which entry is called `s=0`; only the unordered count
and zero-character average in (15) are intrinsic. An inverse-root
nonzero-mode theorem therefore needs a separately specified measurable
section and an anchor hypothesis.

The separation already appears at

```text
y=1/100.
```

Here

```text
d(y)=1,                    d(13y)=0,

sum_s d(y+s/13)=2,

sum_s d((y+s)/13)=1.                              (16)
```

Thus no endpoint convention or averaging argument can identify the two
operations.

## 5. Pure-word typing and the fork boundary

Use THM-2305's scalar-cover notation. For a selected owner `j` and the
other blockers `a,b`, the pure word is

```text
Q_(j,{a})
 =A_0 intersection D_(c_a)
       minus (D_(c_j) union D_(c_b)).               (17)
```

The scalar cover implies, almost everywhere,

```text
A_0 intersection D_(c_j)^c intersection D_(c_b)^c
 subset D_(c_a).                                   (18)
```

Hence the displayed `D_(c_a)` factor in (17) is logically redundant:
removing it does not enlarge the word. Push the remaining positive word
mass forward by its `c_a` coordinate. The resulting weight is supported
in `D`, so (9)--(10) apply to translation of the **actual named danger
factor**. In THM-2337's factor-coloured typing, this is a genuine
nontrivial `beta_a` response of the pure-word mass. The pure `{b}` case
is symmetric.

The statement is sharp at forks. After removing `D_(c_a)` from

```text
Q_(j,{a,b})
 =A_0 intersection D_(c_a) intersection D_(c_b)
       minus D_(c_j),                               (19)
```

the remaining `D_(c_b)` already satisfies the cover, so it does not force
`D_(c_a)`. Both values of the `a`-danger bit occur in the exact Boolean
truth table. Likewise, removing a target complement from a base owner
word is not redundant. Those roles receive only an auxiliary
duplicate-factor probe unless extra geometry is proved.

The successor term in (4) has a literal dynamical meaning after the word
is transported: it measures overlap with the next thirteen-adic danger
scale. Expiration information may set or bound it. Discarding it is
exactly the forbidden current-bit/successor-bit collapse.

## 6. Scope and noncomposition

The theorem proves

```text
positive pure-word role mass
  -> a nonzero factor-coloured first-jet harmonic.                 (20)
```

It does not prove that the full complex word/deep/bare marked current has
the same harmonic. Endpoint Fourier coefficients and terminal phase can
cancel after the positive mass is inserted. Moreover THM-2337's exact
gauge can change `beta_a` while fixing the relation target, and
THM-2356's vertical tensor hostile absorbs arbitrary jet-only data into
`B(z)`.

Thus (20) does not force a nonzero target `q`, an all-`91`-unit aggregate,
a scalar-row exclusion, or LRC(14). Its gain is the exact successor
coordinate and a lawful pure-role marginal which a future
target--jet coupling may consume.

## 7. Exact companion

The dependency-free companion uses `Fraction` arithmetic to:

- exhaust the `26` open cells in both pointwise count identities;
- check the rational scale hostile (16);
- sweep fourteen exact values of `rho_+/rho` through every identity and
  lower bound in (7)--(14);
- verify the inverse-root zero-character counts and the two
  lift-dependent negative nonzero-mode hostiles (15c)--(15d); and
- exhaust the seven nonempty blocker truth assignments, proving pure
  danger redundancy and fork/complement nonredundancy.

Reproduce with

```bash
python3 04-computation/thirteen_shift_successor_role_jet_thm2362.py
python3 -O 04-computation/thirteen_shift_successor_role_jet_thm2362.py
```

Both transcripts must match

```text
05-knowledge/results/thirteen_shift_successor_role_jet_thm2362.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

The first independent audit correctly checked the target-shift theorem
but inferred inverse-root nonzero-mode sums from the unordered count
without checking the chosen `s=0` sheet. A second hostile audit supplied
(15c)--(15d). The repaired theorem retains only the a.e. inverse-root
counts and explicitly records the missing section/anchor coordinate.
The target-shift identities (3), successor formulas (7)--(14), pure-word
typing, fork boundary, and THM-2364 remain unchanged. Repaired normal and
optimized transcripts, stored output, LF hashes, and documentation routing
are checked below. QED.
