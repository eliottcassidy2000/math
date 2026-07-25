---
id: THM-2343
title: "Deep-comb affine target catalyst and inverse-character boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. On every
  positive strict shallow-owner word stratum covered by THM-2327, the
  deepest-comb phase in THM-2334's 169-twist current is the character of
  a nonzero pure deepest-axis vector p. Removing that known character
  defines a phase-free semantic response K, and full target coefficients
  obey A_H(q)=A_K(q-p). Consequently a constant K lands entirely at the
  nonzero target p. The full current is supported only at zero exactly
  when K is a scalar multiple of the single inverse character chi_(-p);
  its total nonzero-target energy is the squared distance from that
  character line. Translating the pure-other, pure-deep, and fork masks
  partitions every phase-free target except -p into affine sets of sizes
  12,12,144; the phase-free zero target lies in the pure-deep set. The
  theorem does not exclude the inverse character, force the word-matching
  affine energy positive, retain an all-91-unit/visible aggregate, exclude
  a scalar row, or prove LRC(14).
source: codex-2026-07-25-deep-comb-affine-catalyst
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2327-two-colour-marked-unit-c3-triangle
  - THM-2334-relation-residue-current-and-character-twist-pushforward
related:
  - THM-2329-boundary-triple-rerooting-and-transverse-gain-obstruction
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
  - THM-2340-owner-word-anova-target-landing
  - THM-2344-correlation-inverse-rigidity-and-aligned-tooth-twist-hostile
script: 04-computation/lrc14_deep_comb_affine_catalyst_thm2343.py
output: 05-knowledge/results/lrc14_deep_comb_affine_catalyst_thm2343.out
script_sha256: 09d6d41dd1f2988f39af997efd7f322617edae89cdb31806f1a7469ba9a31f4b
output_sha256: 5d8e769557f576637b6485a55b9061594de6fdd2c5619fe2404d2748f02abe26
hash_basis: working-tree bytes (LF)
---

# THM-2343 -- the deepest comb is a target translation

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The independent audit rederived the deep-axis typing, all Fourier signs,
the affine translation law and inverse-character iff, the Parseval distance
identity, the explicit coordinate escape tests, and the `12+12+144`
partition.  It also replayed the exact companion under ordinary and
optimized Python with byte-identical output.

THM-2334 reduces nonzero-target landing to nonconstancy of a `13 x 13`
twist current `H`.  That criterion treats the deepest-comb phase and the
endpoint/word response together.  Separating them makes the hostile
boundary much thinner.

The deep leg contributes a known nontrivial character:

```text
H(ell)=chi_ell(p)K(ell),                            (1)
```

where `p` is a nonzero pure target-axis vector.  Fourier modulation is
translation, so the final target current is the phase-free current shifted
by `p`.  Thus a phase-free response concentrated at zero is useful, not
hostile: the comb moves it to the desired nonzero deepest axis.

The genuinely bad zero-only full current requires the exact inverse
character `K=c chi_(-p)`.  Approximate constancy, arbitrary cancellation,
or a general zero-mode concentration is not enough.

## 1. The deepest leg has a nonzero pure target charge

Use THM-2327 on a positive strict shallow-owner word stratum.  Let `j` be
the selected shallow owner and let `i_3` be the labelled coordinate of the
deepest blocker:

```text
w_(i_3)=c_3,
j in {1,2},
gcd(m,91)=1.                                        (2)
```

The two target blockers in THM-2309's quotient are the labels other than
`j`.  Since `j` is shallow, `c_3` is one of them.  Relabel the ordered
target coordinates as

```text
(o,d)=(other target, deepest target).
```

Modulo thirteen,

```text
K_13
 =L_13 direct_sum span(e_o,e_d),
G=K_13/L_13 isomorphic to F_13^2.                  (3)
```

Because `13|c_3`, the coordinate vector `e_(i_3)` lies in `K_13`, and
(3) identifies its quotient class with `e_d`.  Put

```text
p=pi(m e_(i_3))
  =(0,m mod 13) in G.                               (4)
```

Equation (2) gives

```text
p!=0.                                               (5)
```

Thus the deepest comb has an intrinsic pure-axis target charge.  This is
not a cosmetic orientation: the valuation order names the deepest label.

## 2. Demodulating the full semantic twist current

Let `G^` be the character group of `G`.  THM-2334 writes the boundary
value of the full semantic target twist as

```text
H(ell)
 =d_hat(m) chi_ell(p)
   L^ell(X)conjugate(F^ell(Y)),                     (6)
```

where `L^ell` retains every terminal-word factor, `F^ell` is the bare
endpoint, and the coordinate shifts defining both coefficients have
ordinary `L^1` boundary limits.

Define the phase-free response

```text
K(ell)
 :=d_hat(m)L^ell(X)conjugate(F^ell(Y)).             (7)
```

Then (1) holds.  Both `H` and `K` are functions on `G^`, not on a choice
of character representative: `chi_ell(p)` is well-defined because
`p in G`.

At the trivial character,

```text
K(0)=H(0)
 =d_hat(m)(1_(E_Q))_hat(X)
   conjugate((1_E)_hat(Y))
 !=0                                               (8)
```

by THM-2327.  Hence `K` itself is nonzero.

Use normalized inverse transforms

```text
A_H(q)=1/169 sum_(ell in G^)conjugate(chi_ell(q))H(ell),

A_K(r)=1/169 sum_(ell in G^)conjugate(chi_ell(r))K(ell).   (9)
```

Substituting (1) into (9) gives the exact affine landing law

```text
A_H(q)=A_K(q-p).                                   (10)
```

No limiting interchange remains: these are finite transforms of the
already existing boundary values.

## 3. Constant phase-free response is a positive control

Suppose

```text
K(ell)=c!=0                    for every ell.       (11)
```

Then finite orthogonality gives

```text
A_K(r)=c 1_(r=0),
A_H(q)=c 1_(q=p).                                  (12)
```

By (5), the full current lands at a nonzero pure deepest-axis target.
Thus "`K` is constant" is the opposite of THM-2334's bad "`H` is
constant" boundary.  The deepest comb catalyses the zero phase-free mode
into a nonzero target.

More generally a phase-free pure character

```text
K(ell)=c chi_ell(r_0)
```

lands the full current at the single target

```text
q=p+r_0.                                           (13)
```

The deep leg translates the complete target support, not only its zero
fibre.

## 4. The zero-only boundary is one inverse character

Equation (10) gives the exact equivalence

```text
A_H(q)=0 for every q!=0

 iff A_K(r)=0 for every r!=-p

 iff K(ell)=c chi_ell(-p) for every ell,             (14)
```

where `c=A_H(0)`.  The last step is finite Fourier inversion.

This can also be stated as an exact distance.  Parseval and (10) give

```text
sum_(q!=0)|A_H(q)|^2
 =sum_(r!=-p)|A_K(r)|^2                            (15)

 =1/169 sum_ell |K(ell)|^2
  -|1/169 sum_ell chi_ell(p)K(ell)|^2

 =1/169 sum_ell
   |K(ell)-c chi_ell(-p)|^2.                        (16)
```

The final equality is orthogonal projection onto the unit-norm character
line spanned by `chi_(-p)`.  Therefore

```text
some nonzero full target survives
 iff K is not the inverse deep character line.      (17)
```

This is strictly sharper than asking vaguely whether endpoint twists
vary.  They may be nonconstant and still lie on the bad character line;
they may be constant and then land successfully at `p`.

In the ordered dual coordinates `ell=(s,t)` and with
`p=(0,m)`, the entire bad line is explicitly

```text
K(s,t)=c zeta^(-m t).                               (17a)
```

It is constant in the other-target character `s` and is one fixed
geometric progression in the deep character `t`.  Consequently any one
of

```text
K(s,t)!=K(0,t) for some s,t,

K(0,t+1)!=zeta^(-m)K(0,t) for some t,

K(s,t)^13!=K(0,0)^13 for some s,t                  (17b)
```

is already a certificate of nonzero target survival.  The converses of
the individual tests in (17b) are not asserted; the complete iff remains
(17a).  This turns the next computation from a generic variance search
into a rigid character-matching audit.

## 5. Affine word-support masks

The three nonzero target loci in the ordered `(o,d)` coordinates are

```text
S_o ={(x,0):x!=0},
S_d ={(0,y):y!=0},
S_od={(x,y):x!=0,y!=0}.                            (18)
```

By (10), their phase-free preimages are

```text
S_o-p
 ={(x,-m):x!=0},

S_d-p
 ={(0,y-m):y!=0},

S_od-p
 ={(x,y-m):x!=0,y!=0}.                             (19)
```

They remain disjoint and have sizes

```text
12,12,144.                                         (20)
```

Together they partition

```text
G minus {-p}.                                      (21)
```

The excluded point `-p` is exactly the inverse-character coefficient in
(14).  A second useful asymmetry is

```text
0 belongs to S_d-p,
0 does not belong to S_o-p or S_od-p.               (22)
```

Hence a phase-free zero-target aggregate already lands in the correct
full target-support type when the terminal word is the pure deepest-target
word.  For the pure-other and fork words, a nonzero affine component is
still required.

The exact word-matching energies are

```text
E_o =sum_(r in S_o-p)|A_K(r)|^2,
E_d =sum_(r in S_d-p)|A_K(r)|^2,
E_od=sum_(r in S_od-p)|A_K(r)|^2.                  (23)
```

They satisfy

```text
E_o+E_d+E_od
 =sum_(q!=0)|A_H(q)|^2,                            (24)
```

and the appropriate terminal word lands iff its corresponding energy in
(23) is positive.  THM-2340's row/column/interaction ANOVA is the
equivalent formulation after remodulating by `chi_p`.

## 6. Response and tournament audit

The carrier is an affine action, not a tournament:

```text
source:
  phase-free full semantic response K;

map:
  multiply by the known deepest character chi_p;

effect:
  translate every target coefficient by p;

preserved:
  all coefficient values, total energy, target multiplicity, word
  factors, and the ordered other/deep coordinates;

destroyed by forgetting the deep phase:
  which affine point is the bad inverse coefficient;

cheapest decisive test:
  the one-character distance (16), followed by the three affine energies
  (23).                                             (25)
```

The deepest label provides a lawful direction, but it does not orient a
fork between two equivalent targets.  A fork remains the `144`-point
affine mixed component.

The coin-checksum analogy is also exact at the response level.  The deep
character supplies a deterministic nonzero translation.  As its unit
residue varies abstractly, the twelve values traverse the punctured
deepest axis; THM-2327 supplies one fixed such value for the current.
Since the target group has odd order, this is a `13`-cycle translation
rather than an involutive antipode.  Its useful consequence is support
relocation, not binary orbit bisection.

## 7. Sharp controls and remaining boundary

Two one-line controls attain both extremes:

```text
K=c:
  A_H=c delta_p,              nonzero pure-deep landing;

K=c chi_(-p):
  A_H=c delta_0,              zero-only hostile.    (26)
```

These are exact finite-group responses, not claims that either is
realized by a canonical interval current.  They show that the character
in (14) is the sharp failure boundary.

The live analytic statement has now narrowed to:

```text
rule out K in span(chi_(-p))                       [nonzero target],

or more precisely prove E_sigma>0 in (23)          [word support],

then retain an all-91-unit/bounded visible address
and THM-2303 terminal-component phase.              (27)
```

No scalar profile is excluded.  The exact ledger remains `165`; the
repeated-first and alternative resonance branches remain outside the
THM-2327 input; and LRC(14) remains open.

## 8. Exact companion

The companion checks the complete target-translation orthogonality
kernel, all twelve nonzero deepest phases, the affine `12,12,144` mask
partition, and the special membership (22).  Exact rational controls
verify the inverse-character zero-only hostile, the constant-response
deep-axis catalyst, and the projection-energy identity (15).  Every
load-bearing check raises explicitly under ordinary and optimized Python.

Reproduce with

```bash
python3 04-computation/lrc14_deep_comb_affine_catalyst_thm2343.py
python3 -O 04-computation/lrc14_deep_comb_affine_catalyst_thm2343.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_deep_comb_affine_catalyst_thm2343.out
```

byte-for-byte after LF normalization.
