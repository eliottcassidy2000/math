---
id: THM-3718
title: "LRC complete-atom orbit-defect saturation and semantic boundary"
status: >
  PROVED + VERIFIED COMPOSITION.  For any hypothetical scalar-cover
  realization among the current 165 first-depth-one rows, THM-2436 places
  the counterexample in nu_7(c_3)<=M, and THM-2445 supplies a rational
  nonnegative lawful THM-2365 table T_0 with positive
  circulant-complement energy D_0.
  THM-2452 selects one of 128 complete matched Boolean atoms whose table H
  has D_H at least D_0/16384.  Therefore both THM-3713 detector banks fire;
  one rational three-site defect reaches all twelve nontrivial deep
  residues and, residue by residue, an exact deep frequency coprime to 91
  together with an ordinary frequency.
  The THM-2530/2531 anchored path cone further transfers some detector
  defect, with amplitude loss at most 23, to a literal
  occupied-to-excluded-target selector class.  The atom, selector slope,
  target, address, and sufficiently late word clock are adaptive.  The
  excluded target is not a semantic arrival, so no scalar row is removed
  and LRC(14) remains open.
source: codex-lrc14-20260822 / THM-2452 to THM-3713 composition audit
audit: >
  PASS -- two independent line audits verified the 165-row conditional
  quantifier; literal T_0/T_omega projection, normalization, interval and
  diagonal-zero types; all bare and delayed-word constants; all-twelve
  cyclotomic landing and per-residue extraction; the two-gauge selector
  telescope and factor-23 loss; and the semantic non-arrival boundary.
  Normal and optimized exact companions reproduce the stored transcript.
depends_on:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2436-punctured-ninety-one-stalk-repeated-step-spectrum
  - THM-2445-twenty-four-cell-graft-owner-conditioning
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
  - THM-2530-anchored-deep-gram-cone-and-lossless-skew-target-refinement
  - THM-2531-prime-necklace-guard-boundary-selector
  - THM-3713-lrc-orbitwise-deep-offset-detector-and-successor-quotient-hostile
related:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-3731-lrc14-valuation-owner-equivariance-and-semantic-packet-boundary
script: 04-computation/lrc_complete_atom_orbit_selector_bridge_thm3718.py
output: 05-knowledge/results/lrc_complete_atom_orbit_selector_bridge_thm3718.out
script_sha256: e2b09c343ae9f3859d7367dab51e003c5fb73da83042943d95dceee31bc1ab2f
output_sha256: c7e4491298e527d07034de8b222cb9f388802fe5fab71b2a9e5156c263fec1bf
hash_basis: raw LF bytes
---

# THM-3718 -- the adaptive complete atom already has an orbit defect

**PROVED + VERIFIED COMPOSITION.**  This theorem closes an abstract
nonvanishing question left visible by THM-3713, but only after an adaptive
complete-atom selection.  It does not identify that atom with the canonical
exclusive-owner atom or an owner-supported word/root/arrival packet.

## 1. Exact inherited universe

Let a hypothetical scalar-cover realization lie among the current 165
first-depth-one rows.  THM-2436 proves that every such counterexample lies
in its complete post-closure branch

```text
nu_7(c_3)<=M.                                             (1)
```

THM-2445 chooses a maximal-septimal-depth label `q_*`, uses it as the actual `c_3`
graft in a THM-2309/2365 owner pivot, and constructs

```text
T_0:F_13^3 -> Q_(>=0),

D_0=||(I-P)T_0||_2^2>0.                                 (2)
```

Here `P` is literally THM-2365's normalized circulant projection

```text
(P H)(r,s,t)=13^-2 sum_(a,b)H(r+b,s+a,t+b),              (3)
```

not an analogous projection on another table.  The deepest present factor
is the complement of the probe at `r=t`, so `T_0(t,s,t)=0`.

THM-2452 partitions the moving endpoint into 128 rational complete Boolean
atoms, copies the selected atom to the bare endpoint by pointwise
idempotence, and obtains one matched table `H=T_omega` satisfying

```text
H:F_13^3 -> Q_(>=0),
H(t,s,t)=0,
D_H=||(I-P)H||_2^2>=D_0/128^2=D_0/16384.                (4)
```

The norm and normalization in (4) are exactly those in THM-2365 and
THM-3713.  In particular, the selected table is a lawful rational interval
tensor to which THM-3713's universal fiberization, detector-frame, and
rational-saturation conclusions apply.

## 2. Both complete orbit detectors fire

Put

```text
h_u(s,t)=H(t+u,s,t).
```

THM-3713 gives the sharp frames

```text
E_orb>=4 sin(pi/13)^2 D_H,

E_3:=13^-3 sum_(u,s,t)|M_u(s,t)|^2
   >=16 sin(pi/13)^4 D_H,                              (5)

M_u(s,t)=h_u(s,t)+h_u(s+1,t)-2h_u(s,t+1).
```

Equations (4)--(5) imply that some coordinate edge and, separately, some
three-site defect are nonzero.  Since there are `2*13^3` coordinate edges
and `13^3` three-site defects, there are witnesses `delta` and `Delta` with

```text
|delta|^2>=sin(pi/13)^2 D_0/8192,

|Delta|^2>=sin(pi/13)^4 D_0/1024.                      (6)
```

No row-uniform absolute floor is claimed: `D_0` depends on the inherited
row and chosen graft.

## 3. One defect reaches every deep residue

Fix a centre `(s_0,t_0)` at which the three-site profile

```text
d_u=M_u(s_0,t_0)                                      (7)
```

is nonzero.  Lawful diagonal zero makes `h_0=0`, hence `d_0=0`.  All entries
are rational.  THM-3713's anchored cyclotomic argument therefore gives

```text
sum_u d_u zeta^(a u)!=0              for every a!=0. (8)
```

Its exact Fourier typing is

```text
calH(a,b,q)=B(a,b,q-a),                               (9)

d_hat(a)=13 sum_(b,q)
  (1+zeta^(-b)-2zeta^(-q)) calH(a,b,q)
  zeta^(-b s_0-q t_0).
```

The multiplier vanishes only at the zero target.  Thus for each `a!=0`
there is a nonzero target `(b,q)` with

```text
B(a,b,q-a)!=0.                                       (10)
```

Because `H` is the actual rational interval tensor, THM-2365's absolutely
iterated extraction then supplies, for every `a!=0`, some ordinary
frequency `X_a` and deep frequency `m_a` such that

```text
m_a=a mod 13,
gcd(m_a,91)=1,
A_(X_a,m_a)(b,q)!=0.                                 (11)
```

The target, `X_a`, and `m_a` may all vary with `a`; (11) does not preserve
the old marked triangle or give one common all-coordinate-unit address.

## 4. The defect has a literal selector refinement

THM-2530 applies cellwise to the same selected table.  Let `K_(s,t)` be its
anchored deep-window Gram matrix.  Its diagonal is exactly

```text
h_j(s,t)=K_(s,t)(j,j),                 1<=j<=12.     (12)
```

Write `alpha_j,beta_j` for the singleton and adjacent-pair path masses, with
`beta_0=beta_12=0`.  THM-2531's two slope gauges give the nonnegative
occupied-to-excluded-target selector masses

```text
gamma_j^+=alpha_j+beta_j,
gamma_j^-=alpha_j+beta_(j-1),                       (13)
```

and hence the exact inverse identity

```text
h_j=gamma_j^+ + sum_(k<j)(gamma_k^+-gamma_k^-).     (14)
```

Let `L` be either a coordinate edge difference or the three-site mask in
the `(s,t)` variables.  It commutes with (14).  If `Lh_j` is nonzero, one
of the at most `2j-1<=23` terms on the right obeys

```text
|L gamma_k^epsilon|>=|Lh_j|/23.                     (15)
```

Applying this to the two witnesses in (6) gives a literal selector edge and
a literal selector three-site defect (not necessarily in the same selector
class) with

```text
|delta_sel|^2>=sin(pi/13)^2 D_0/4333568,

|Delta_sel|^2>=sin(pi/13)^4 D_0/541696.             (16)
```

Fix one witness from (15), including its sign `epsilon`, and define the
selector defect profile

```text
e_k^epsilon=L gamma_k^epsilon,       k in F_13.      (16a)
```

It is rational, nonzero, and `e_0^epsilon=0`.  Hence it has all twelve
nontrivial `C_13` Fourier colours by the same anchored cyclotomic argument.
This is a selector-current colour statement; it is not another
relation-current extraction beyond (10)--(11).

The two signs in (13)--(14) are two chosen target-coordinate selector
gauges `tau=+1,-1`; they are not both inherited physical guard slopes.
Forgetting their orientation destroys the inverse.  Moreover, THM-2531's
selected arc ends at the excluded/safe target root.  That root has zero
target charge and is not a semantic arrival event.

## 5. Delayed words

For a fixed positive rational terminal word `Q_term`, there is

```text
k_0=k_0(row,omega,Q_term)
```

such that the same selected complete atom has positive drift at every
allowed clock `R=13^k`, `k>=k_0`.  The qualitative nonvanishing and colour
conclusions of Sections 2--4 then apply.  Under THM-2452's stronger explicit
half-error threshold one has

```text
D_word>=mu(Q_term)^2 D_0/65536.                     (17)
```

Consequently the squared word-bearing edge, three-site, selector-edge, and
selector-three-site floors are respectively

```text
mu(Q_term)^2 sin(pi/13)^2 D_0/32768,
mu(Q_term)^2 sin(pi/13)^4 D_0/4096,
mu(Q_term)^2 sin(pi/13)^2 D_0/17334272,
mu(Q_term)^2 sin(pi/13)^4 D_0/2166784.              (17a)
```

Without the half-error threshold, only eventual positivity and the
qualitative consequences are asserted.  The theorem does not place the
defect at the prescribed first-expiration clock.

## 6. Exact gain and exact remaining boundary

The proved composition is

```text
hypothetical cover realization on the 165-row frontier
  -> positive isolated lawful graft table
  -> one adaptive complete matched atom with D_H>0
  -> nonzero edge and three-site orbit defects
  -> all twelve deep residues and fresh exact (X,m), gcd(m,91)=1
  -> one literal occupied-to-excluded-target selector defect. (18)
```

What (18) deliberately does not retain is:

- a row-independent atom, target, magnitude, address, or late clock;
- owner support on the canonical exclusive-owner atom, or the THM-2305
  semantic word/repair label;
- a common physical root section or the physical guard orientation;
- positive charge at the excluded target, hence semantic arrival;
- the old marked triangle or the all-nine-coordinate unit projector.

Thus the pure abstract nonvanishing half of THM-3713's target is discharged
on an adaptively selected complete atom.  The live theorem is now a
transport/typing statement: couple one such selector defect to the canonical
owner-supported temporal/root packet, or give the excluded target a lawful
positive arrival meaning, without losing its target character and deep
residue.

No scalar row is removed, and LRC(14) remains open.

## 7. Exact companion

The standard-library companion checks the two-gauge selector telescope on
all 169 target cells for both coordinate edges and the three-site operator,
and checks every denominator in (4), (6), (16), (17), and (17a):

```bash
python -B 04-computation/lrc_complete_atom_orbit_selector_bridge_thm3718.py
python -B -O 04-computation/lrc_complete_atom_orbit_selector_bridge_thm3718.py
```

Both runs reproduce

```text
05-knowledge/results/lrc_complete_atom_orbit_selector_bridge_thm3718.out
```

byte for byte.  QED.
