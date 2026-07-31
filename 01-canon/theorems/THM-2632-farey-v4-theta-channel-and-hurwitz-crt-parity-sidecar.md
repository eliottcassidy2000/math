---
id: THM-2632
title: "Farey V4 theta channels and the full Hurwitz C6 CRT-lift sidecar"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The six mod-two
  residue classes of oriented unimodular Farey
  flanks form the S3 permutahedral C6.  Its inversion-set embedding is the
  isometric Q3 minus the two cyclic tournaments; its three theta cuts are
  the three V4 parity/resolvent channels, and its edges are exactly
  {L,R} times those channels.  THM-2056 defect parity selects one channel
  unless its owner residue is (1,1).  The simultaneous integral 2/3
  quotient is PSL2(Z/6)=S3 times A4, of exponent six and abelianization C6.
  Every normalized seven-edge Hurwitz norm returns its chosen root step in
  this quotient, and the six lifts a+13k exhaust all six C6 shadows.  This
  shadow is integral-lift data, not a PSL2(F13)-conjugacy invariant: the
  original Hurwitz pair returns a mod-two root shadow of order three, while the
  displayed normalizing matrix is singular modulo two and six.  The
  affine S4 detector h^7 h^-1=h^6 recovers exactly the six four-cycle
  cocycles, but no physical LRC connector or Keller exclusion follows.
source: codex-2026-07-27-farey-v4-crt-parity
depends_on:
  - THM-2056-kelvin-polar-farey-defect-certificate
  - THM-2440-sharp-two-comb-centred-window-radius
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2597-six-vertex-bicycle-modular-abelianization-cycle
  - THM-2603-hurwitz-projective-root-owner-atlas-and-nonabelian-seven-edge-trivialization
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2619-psl2f13-seven-edge-norm-minimal-projective-kernel-and-retained-target-obstruction
  - THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary
related:
  - THM-2626-paley-borel-projective-frame-torsor-and-physical-c13-boundary
script: 04-computation/farey_v4_theta_hurwitz_crt_thm2632.py
output: 05-knowledge/results/farey_v4_theta_hurwitz_crt_thm2632.out
script_sha256: f84ec9f806d960170868d144c510bbd5caa2da326d9414169445785f75a1cb16
output_sha256: 050c8f3ac8c2d6827ea0a23703c6118fefa02656fa78930251fb933bfdc1c761
hash_basis: LF-normalized bytes
---

# THM-2632 -- two and three meet in a sixfold CRT lift, not in branch count

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The
dependency-free companion checks every finite statement below with explicit
optimized-mode-safe guards.  The LRC paragraphs are typed
consequences and boundaries of the cited theorems; they remove no row.

There are two genuinely different six-state objects in the modular picture.
The first is the permutahedral hexagon of orders on the three nonzero
directions of `V_4`.  The second is the cyclic CRT lift coordinate
`C_6=C_2 x C_3`.  Forgetting edge colours makes their underlying graphs look
the same.  Keeping the actions shows why they cannot be identified without a
new sidecar.

## 1. The Farey residue graph is an intrinsic three-cut partial cube

Let `U=[u v] in GL_2(Z)` be an oriented unimodular flank and reduce its
columns modulo two.  They are distinct nonzero vectors in

```text
V=F_2^2={0,a,b,c},                         c=a+b.          (1)
```

Define

```text
pi(U)=(u_bar, u_bar+v_bar, v_bar).                          (2)
```

The six possible values of (2) are exactly the six total orders of
`{a,b,c}`.  The two positive Farey child matrices act on the right:

```text
L_bar=[1 1;0 1]:  (x,y,z) -> (x,z,y),
R_bar=[1 0;1 1]:  (x,y,z) -> (y,x,z).                       (3)
```

Consequently the residue graph is

```text
Cay(S_3,{(12),(23)})=C_6.                                  (4)
```

Choose a reference order only to write coordinates.  Send a permutation to
its three pairwise comparison bits.  Right multiplication in (3) flips one
bit, and graph distance equals Hamming distance.  The image is

```text
Q_3 minus {(0,1,0),(1,0,1)}.                               (5)
```

The six retained vertices are the six **transitive tournaments on three
channel vertices**.  The two omitted vertices are the two cyclic
tournaments.  The six residue states are not themselves a tournament.

The coordinate gauge in (5) is not intrinsic, but its three Djoković--Winkler
theta cuts are.  A comparison of two nonzero directions has a unique omitted
third direction, their sum.  By THM-2606 this is exactly one of the three
equivalent data

```text
nonzero V4 direction <-> parity covector <-> quartic 2+2 resolvent channel.
                                                                    (6)
```

Each theta class has two edges: one `L` edge and one `R` edge.  Hence

```text
Edges(C_6)  <->  {L,R} x {three theta channels}.            (7)
```

Left multiplication by `GL_2(F_2)=S_3` is regular on the six states and
commutes with the right child moves.  It preserves the alternating `L/R`
edge colour and has two edge orbits of size three.  This is not the regular
cyclic action of THM-2597.  Indeed `S_3` has no element of order six; a
one-edge rotation of the uncoloured hexagon swaps `L/R` and is not a modular
left action.  The two hexagons coincide only after controlled forgetting.

## 2. THM-2056 parity chooses one theta channel, not a gate sign

For the THM-2056 defect

```text
F_w(d)=||d||^2-91 w.d                                      (8)
```

one has exactly

```text
F_w(d) mod 2 = ((1,1)+w_bar).d_bar.                        (9)
```

If `w_bar!=(1,1)`, the right side is a nonzero covector.  It vanishes on one
nonzero direction and equals one on the other two.  Thus every mod-two Farey
triangle (2) has defect parity pattern

```text
{0,1,1},                                                   (10)
```

and its zero direction names one theta channel in (6).  When
`w_bar=(1,1)`, all three defects are even and no channel is selected.
Equation (9) retains only parity.  It does not determine the sign of (8),
propagate a hull owner, or replace the Gram cross term required by THM-2596.

The Pythagorean ternary branches give a sharp warning.  Their two-dimensional
parameter matrices reduce modulo two to

```text
{J,J,I},                  J=[0 1;1 0],                     (11)
```

so two branches send every Farey residue state to its antipode and one stays
put.  The three Berggren matrices are all identity modulo two.  Thus the
ternary cover is not a `C_3` action.  In THM-2596 it is the complete prefix
code `{rho,lambda rho j,lambda^2}`: symbolic Kraft weights are
`1/2+1/4+1/4=1`, while the actual interval lengths are `1/2,1/6,1/3`.
It covers the same interval but is neither equal-area nor a ternary
permutation.

## 3. The full simultaneous shadow is `C6`, not one parity bit

Put

```text
Q_6=PSL_2(Z/6Z).                                           (12)
```

Chinese remaindering and `-I=I` over `F_2` give

```text
Q_6 = SL_2(F_2) x PSL_2(F_3)
    = S_3 x A_4,                           |Q_6|=72.       (13)
```

Both factors have exponent six.  Moreover

```text
Q_6' = C_3 x V_4,                  Q_6^ab=C_2 x C_3=C_6.  (14)
```

For the standard parabolic

```text
U=[1 1;0 1],                                               (15)
```

its class has order six, intersects (14) trivially, and maps onto the whole
abelianization.  Hence

```text
Q_6=(C_3 x V_4) semidirect <U>.                            (16)
```

Conjugation by `U` inverts the `C_3` factor and cycles the three nonzero
elements of `V_4`.  This is the rigorous co-occurrence of the binary and
ternary grammars: they are the order-two and order-three projections of one
parabolic lift coordinate, not two branch arities.

Every group element `x in Q_6` satisfies `x^7=x`.  Therefore, for every
`g in Q_6` and every integer `a`, the two normalized ordered words obey

```text
(U^a g)^7 g^(-7)       = U^a,
g^7 (g^(-1)U^a)^7      = U^a.                              (17)
```

This is stronger than the successful-norm atlas: closure modulo thirteen is
not needed for (17).  When one of THM-2603/2619's thirty normalized words
does close modulo thirteen, its six integer exponent lifts

```text
a, a+13, ..., a+65                                       (18)
```

all give the same mod-thirteen norm, while `13=1 mod 6` makes their shadows
in (17) all six distinct powers of `U`.  Thus a mod-thirteen root exponent
has a full sixfold CRT lift coordinate.

The `C_2` fixed-section projection recovers THM-2622's pure-linear boundary:
even `a` fixes four `V_4` points and odd `a` fixes two.  On
`V_4 x P^1(F_3)`, the two fixed counts for `a=0,...,5` are

```text
(4,4),(2,1),(4,1),(2,4),(4,1),(2,1).                      (19)
```

These counts detect `gcd(a,6)` but identify `a` with `-a`.  They also assume
zero affine cocycle; THM-2622 shows that an affine lift may instead have no
fixed section.

## 4. The integral normalization is load-bearing

The normalized `U,g_t` pairs are not an invariant of the order-seven
conjugacy class relative to an unspecified integral lift.  For the original
THM-2603 Hurwitz matrices

```text
A=xy=[-1  2; 3 -7],       C=[x,y]=[ 5 -17;-17 58],        (20)
```

the descending seven-edge word is

```text
product_(i=6..0) A^i C A^(-i)=A^7(A^(-1)C)^7.             (21)
```

The same exponent-six argument reduces (21) to **its own root step `C`**.
Modulo two,

```text
A_bar=[1 0;1 1],          C_bar=[1 1;1 0],
ord(C_bar)=3,             ord(U_bar)=2.                    (22)
```

So the original relation does not have the normalized `U` parity shadow.
After THM-2603's unimodular chart `P`, the final matrix used there to
normalize the pair modulo thirteen is

```text
Q=[4 5;0 10],             det(Q)=40=1 mod 13,              (23)
```

but it is singular modulo two and noninvertible modulo six.  It cannot
transport the CRT shadow.  (The complete normalizing matrix is `QP`.)  In the
full quotient `Q_6`, both `C` and `U` have order six: (22) distinguishes only
their mod-two shadows.  Equations (20)--(23) are the sharp boundary:

```text
the seven-edge norm remembers the chosen integral root step;
mod-13 projective conjugacy does not choose that step's CRT lift.           (24)
```

Consequently the six choices in (18) are not a free LRC certificate.  A
physical use needs an integer ancestor already typed by the dynamics.

## 5. The affine seventh-power detector recovers one missing `V4` cocycle

The linear `S_3` shadow forgets the affine term in

```text
AGL_2(F_2)=V_4 semidirect S_3=S_4.                         (25)
```

For an affine permutation `h`, define

```text
delta_7(h)=h^7 h^(-1)=h^6.                                (26)
```

The class census in `S_4` proves

```text
delta_7(h)!=1  iff h is one of the six 4-cycles.           (27)
```

For a four-cycle, `delta_7(h)=h^2` is the nonzero pure translation in the
unique direction fixed by its linear reflection.  That direction is exactly
the omitted diagonal of the fixed THM-2606 channel.  A transposition with the
same linear shadow has `delta_7=1`.  Thus (26) separates the two affine
cocycles that the cubic-resolvent `S_3` shadow conflates.

The boundary is exact.  Identity and the three pure double translations all
have `delta_7=1`; and `h^13=h` for every `h in S_4`.  Equation (26) classifies
an already supplied affine quartic action.  It neither constructs one from
the LRC Hurwitz word nor excludes a Keller monodromy group.

## 6. Ordinary gracefulness spends one of the six edge addresses

The hexagon (4) is a partial cube but is not **ordinary Rosa-graceful**: the
sum of required edge differences `1+...+6` is odd, while modulo two the same
sum counts every cyclic vertex label twice.  Deleting any one edge gives
`P_6`, which is ordinary graceful (for example, successive labels
`0,5,1,4,2,3` give differences `5,4,3,2,1`).

By (7), the six possible graceful repairs are the six `(L/R,theta)`
addresses.  Left `S_3` has two orbits of three, so it fixes no repair.  This
is ordinary gracefulness, not the distinct theta-graceful notion from the
partial-cube literature.  It supplies no invariant owner or root choice.

If an affine `V_4` torsor is given both a distinguished origin and a Farey
order of its three nonzero directions, the data determine a full affine
frame: `4*6=24=|AGL_2(F_2)|`.  THM-2606's Feuerbach label supplies the origin
on its sign torsor; (2) supplies only the direction order on a separately
identified residue torsor.  No theorem identifies those two torsors, so the
numerical product does not manufacture a geometric transfer.

## 7. Two-comb radius is not Farey-quasiconcave

Use THM-2440's closed/almost-everywhere normalized two-comb radius.  The
fractions `1/13` and `1/14` are Farey neighbours, with mediant `2/27`, but
direct tooth endpoints give

```text
r_cl(1,13)=15/182,
r_cl(1,14)=15/196,
r_cl(2,27)=5/126 < min(15/182,15/196).                    (28)
```

For example, the central `2`-tooth ends at `1/28`, and the first positive
`27`-tooth ends at `15/378=5/126`; the next gap stops the component.  Hence
the radius is not quasiconcave under the Farey mediant.  In the literal-open
convention, the seam point makes `r_o(1,13)=1/14`, not `15/182`.

THM-2056's quadratic defect is controlled on an acute Farey cone because it
retains the Gram cross term.  Equation (28) proves that its fan propagation
cannot be transferred to the two-comb radius scalar without a tooth/seam
sidecar.

## 8. Transfer contract and decisive next test

The exact new dictionary is

| source | map | preserved | lost / required sidecar |
|---|---|---|---|
| oriented Farey flank | reduction mod 2 | three theta cuts and `L/R` edge colour | integral Gram, owner, height |
| THM-2056 defect | (9) | one parity channel or the zero boundary | defect sign and cone propagation |
| integral Hurwitz word | reduction mod 6 | chosen root step and its full `C6` lift | physical choice of integral ancestor |
| affine quartic action | `delta_7` | four-cycle translation cocycle | double translations and geometric source |
| two-comb pair | radius (28) | connected central coverage | Farey propagation |

After choosing an orientation of the three theta channels, (7) is a six-set
and `C_6=C_2 x C_3` is another.  A typed physical integer ancestor could map
between them by matching its mod-two coordinate to `L/R` and its mod-three
coordinate to the owner-selected theta channel.  That equality is **not**
present in current canon.  It is the cheapest decisive LRC test suggested by
this theorem.  Choosing a lift after seeing which edge closes is only gauge
selection.

No LRC row is excluded, no positive chronological whole-head transition is
constructed, and no Jacobian or Dixmier conjecture is settled here.

## 9. Exact reproduction

Run

```bash
python 04-computation/farey_v4_theta_hurwitz_crt_thm2632.py
python -O 04-computation/farey_v4_theta_hurwitz_crt_thm2632.py
```

The companion checks the six-state partial-cube/tournament atlas; all
`(L/R,theta)` edges and ordinary graceful deletions; every THM-2056 parity
row; the ternary antipode/stay shadows; all `72*6*2=864` `Q_6` norm laws;
the derived subgroup, semidirect action, fixed profile, thirty normalized
mod-thirteen successes and all six lifts; the original-pair normalization
hostile; all `24` affine `S_4` maps and `delta_7`; and the exact three-radius
Farey hostile.  After LF normalization, normal and optimized output must
byte-match the stored transcript and end in `ALL CHECKS PASSED`.
