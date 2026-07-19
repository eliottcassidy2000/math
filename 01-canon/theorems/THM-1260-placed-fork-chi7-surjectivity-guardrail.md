---
id: THM-1260
title: A PLACED SHARP TWO-WALL FORK HAS SURJECTIVE chi7 COLOUR -- every compact toothpick rung, both binary phase sides, and all 128 Fano colour words survive locally
status: PROVED (explicit all-rung sharp same-provider fork; all four target/provider chi7 pairs; arbitrary seven-label speed-colour word; binary centered phase order; exact second blocker on the correct marked side; compact ratio box; dependency-free 85,760-row referee; sorry-free Lean arithmetic core). This is a local placed-interval guardrail, not a six-comb cover, an embedding in a deletion-minimal cover, or a simultaneous blocker-cycle realization
source: codex-2026-07-19 placed-Fano/chi7 probe
depends_on: [THM-1156, THM-1247, THM-1248, THM-1252, THM-1254, THM-1256]
related: [THM-1179, THM-1234, THM-1262, HYP-7678, HYP-7870]
script: 04-computation/lrc14_placed_fork_chi7_surjectivity_thm1260.py
output: 05-knowledge/results/lrc14_placed_fork_chi7_surjectivity_thm1260.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCPlacedForkChi7Surjectivity.lean
script_sha256: bce7538d48df87b78bcf864442335b6156b94bcfa835415749cb0a6dc0dbc709
output_sha256: 218c23c21150f91ca65b87bbc4e21a835813bbaf2a6688c28a9a61e24dc14482
formalization_sha256: 786c201095e69658658a2e2e52cca9a31b40540da3c4c332d1229633e70ff432
---

# THM-1260 -- placed-fork `chi_7` surjectivity

THM-1156 proves a real `chi_7` law for **zero** tooth seams.  THM-1247 then
shows that the contracted `q=15` Fano stalk and speed `chi_7` are independent
before off-grid germs are placed.  THM-1252--1256 now place a blocker target,
both of its wall seams, its centered phase, and the next marked blocker on one
chronological carrier.  The natural question is whether that extra placement
finally forces a Fano line-colour law.

It does not at one fork.  In fact the sharpest positive same-provider fork is
already colour-surjective.

## 1. Exact statement

Strip the full `7`-adic factor and write

```text
epsilon(v)=Legendre_7(v/7^nu7(v)) in {+1,-1}.          (1)
```

Fix arbitrarily:

```text
1<=r<=335,                 delta in {0,1},             (2)
eta=(eta_c,eta_i,eta_j,eta_h,eta_b,eta_u,eta_v)
       in {+1,-1}^7.                                  (3)
```

There are a carrier `c`, six distinct faster speeds

```text
i,j,h,b,u,v>c,                                        (4)
```

and integer tooth addresses with all of the following properties.

1. Their seven speed colours, in the named order, are exactly `eta`.
2. The half-gap centered spoke of source `i` is `t_i=1/2`; the even target
   `j<i` is strictly dangerous there.
3. The designated target tooth `J` is flanked in chronological order by two
   distinct `h`-teeth:

   ```text
   h -> j -> h.                                       (5)
   ```

   Their address return is exactly `r`.  Both raw overlaps are the sharp
   quantum

   ```text
   1/[14 lcm(j,h)].                                   (6)
   ```

4. At the centered spoke `t_j`, a further named label `b` is exactly
   dangerous.  Its designated blocker tooth lies strictly on the lower side of `J` when
   `delta=0` and strictly on the upper side when `delta=1`.  Thus tooth order
   agrees with the binary phase order of THM-1256.
5. Every speed lies in the compact THM-1233 box

   ```text
   c<d<2345c.                                         (7)
   ```

The construction gives a literal interval-chain model of the local predicates
seen on the common marked tooth carrier: a placed two-wall fork and the next
designated blocker tooth.  It does **not** assert that this model embeds in a
six-comb cover of the whole slow gap, nor that six such fragments close to one
blocker cycle.

## 2. The sharp fork realizes all four endpoint colours

Put

```text
e=1  if r is odd,                 e=2 if r is even.   (8)
```

Choose a large integer `D` and define

```text
j=eD,
h=e[(7r-1)D+1]=(7r-1)j+e.                            (9)
```

Then

```text
gcd(j,h)=e,                    h-(7r-1)j=e.           (10)
```

If `r` is odd, choose `D` even.  In both parity branches, `j` is even and
`h` has the same parity as `r`.  Hence

```text
N=j/2,
M_-=(h-r)/2,                  M_+=(h+r)/2             (11)
```

are integers.  Let `J` be the `N`-th `j`-tooth and let the outer intervals
be the `M_-`-th and `M_+`-th `h`-teeth.  Direct endpoint subtraction gives

```text
R_(h,M_-)-L_(j,N)=e/(14jh),
R_(j,N)-L_(h,M_+)=e/(14jh).                           (12)
```

Since `lcm(j,h)=jh/e`, both values in (12) are exactly (6).  The two outer
teeth are separated, the intersections are proper, and their address
difference is `r`.  Thus (9) realizes every rung in (2) at the smallest
possible positive gcd-sheet detuning.

Modulo seven, the reduced speeds satisfy

```text
j/e ==D,                     h/e ==1-D.               (13)
```

Both `e=1` and `e=2` have positive quadratic character.  The four rows

```text
(epsilon(j),epsilon(h))       D mod 7
(+,+)                            4
(+,-)                            2
(-,+)                            6
(-,-)                            3                  (14)
```

prove surjectivity.  This is the key obstruction to a local Fano law: even
the **sharp** overlap equation leaves all four ordered colour pairs alive.
For an exact zero seam the reduced equation is `D+H==0 mod 7`, which flips
the character as in THM-1156.  For the positive sharp seam it is instead

```text
D+H==1 mod 7,                                         (15)
```

and (14) shows that (15) has no binary-character content.

## 3. Complete seven-colour freedom and binary landing

Take `D` sufficiently large in its row of (14), preserving the required
parity.  Choose odd steps `1<=s_c,s_i<=13` so that

```text
c=j-s_c,                         i=j+s_i              (16)
```

have the prescribed colours `eta_c,eta_i`.  This is always possible because
the seven odd steps `1,3,...,13` traverse every residue modulo seven.  Put

```text
k=(c-1)/2.                                            (17)
```

Then the centre of the `k`-th carrier gap is `1/2`, and both `c` and `i` are
half-integral there.  Thus `t_i=1/2` is an exact centered carrier spoke,
while the even target `j` is integral and hence strictly dangerous.

The target clock `Q_j=c+j` is odd.  Its two nearest centered numerators are

```text
P_j=k+j/2+delta,                  delta=0 or 1.        (18)
```

Consequently

```text
2P_j-Q_j=-1  if delta=0,
2P_j-Q_j=+1  if delta=1,                              (19)
```

so `t_j=P_j/Q_j` lies respectively below or above `1/2`.  This is exactly
THM-1248's binary relative-address digit and THM-1256's phase-side law.

To continue the marked blocker path, choose

```text
b=mQ_j,                         m in {2,3}.            (20)
```

At `t_j`, the `b`-tooth of address `mP_j` is centered exactly on the phase.
Since `epsilon(2)=+1` and `epsilon(3)=-1`, one of the two choices in (20)
gives the prescribed `eta_b`, irrespective of the colour of `Q_j`.  Moreover
that tooth is disjoint from `J` and lies on the side prescribed by `delta`:

```text
1/(14j)+1/(14b)<1/(2Q_j)=|t_j-1/2|.                  (21)
```

Both marked teeth lie strictly inside the same carrier gap.  Equation (21)
is placement, not a claim that the intervening corridor is already covered;
in a hypothetical full cover THM-1256 would fill it with the corresponding
chronological subword.

Finally choose `u,v` successively above `j+28` in any nonzero residue classes
of the prescribed characters.  They are distinct from (16), (20), and each
other.  For the explicit referee choice `j>=100003`, all named speeds fit
below `h`; at the terminal rung,

```text
h<=2344j+2,
c>=j-13,
2345c-h>=j-2345*13-2>0.                              (22)
```

Smaller rungs only improve the margin.  This proves (7) and completes the
construction for every word (3).

## 4. What the guardrail rules out

Label the seven named speeds by the points of any Fano plane.  Equation (3)
realizes every one of its `128` binary point-colourings for each fixed rung
and binary phase side.  Therefore none of the following can be a consequence
of one placed fork plus binary landing:

* a prescribed parity or product of speed `chi_7` colours on a Fano line;
* a monochromatic-line prohibition;
* a Paley orientation or nontransitive tournament fingerprint determined
  only by the seven speed residues;
* a shorter toothpick-rung alphabet on same-colour or cross-colour pairs.

There is no conflict with THM-1156.  The two wall events in (12) are positive
overlaps, not zero seams.  Nor does this refute THM-1179's quantitative
same-colour **global pair-density** floor; the present stalk is local and its
sharp overlap can be arbitrarily small after projective dilation.

The line-shape itself explains the failure.  A distinct-provider target is
the open `V`

```text
h_- -- j -- h_+,                                      (23)
```

whose outer teeth are nonconsecutive and do not overlap by minimality.  A
same-provider target degenerates the two outer vertices, as in (5).  Neither
local object is a simple three-edge Fano line.  The missing third edge can
only come from **owner reuse at another word position**.

## 5. Sharpened next lemma: seam-digit triangle closure

The faithful residue is not a point colour but an oriented handoff digit.
For a chronological overlap from the `(a,m)` tooth to the `(b,n)` tooth, put

```text
W=b(14m+1)-a(14n-1)>0,
g=gcd(a,b),             A=a/g, B=b/g, w=W/g.          (24)
```

Then the exact local law is

```text
w==A+B mod 14.                                        (25)
```

The zero-seam `chi_7` flip is the missing-digit fibre `w=0`; THM-1260 lives
on the sharp positive fibre `w=1`.  Quotienting (24) to `epsilon(a),epsilon(b)`
destroys the coordinate that distinguishes those fibres.

The first genuinely untested Fano carrier is therefore a **triangle of
handoff occurrences**, not three runner labels.  A useful next lemma must
prove a dichotomy of this exact form:

```text
some nonbacktracking fork h_- -> j -> h_+
has an h_-/h_+ handoff occurrence elsewhere in the same word,

or the triangle-free owner-transition support forces enough repeated
star/backtrack lcm weight to close the full THM-1253 invoice.             (26)
```

In the first branch, retain the three digits `w` and their absolute tooth
addresses and test (25) around the resulting owner-reuse circuit.  In the
second branch, use the literal occurrence multiplicities rather than Fano
colour.  A line law that does not retain the seam digits, occurrence
locations, and shared-owner transport is ruled out by this theorem.

THM-1262 resolves the first cycle-length test after this guardrail: a
coherent blocker two-cycle cannot use the adjacent marked-inversion branch
and must open an aligned corridor through a third owner.  That theorem
produces one protected bridge, not yet the closing outer-owner edge in (26).
Thus the sharpened Fano question is incidence among **several exported
third-owner bridges** (or the accumulated triangle-free star debt), rather
than the colour of either endpoint of one bridge.

## 6. Tournament and alternate-carrier audit

On the six runner labels, speed order is the transitive tournament with score
histogram `(0,1,2,3,4,5)`, no directed triangles, six singleton SCCs, and one
Hamiltonian path.  It preserves only the compact ordering and loses (11),
(12), (18), and (24).

On the three local tooth positions, chronology is again a path.  The useful
binary switch is phase side `delta`, but (14) proves it is independent of the
speed-colour word.  Fano points and lines, Paley residues, centered phases,
individual teeth, wall events, seam digits, owner-pair occurrences, and
matroid circuits were all challenged as vertices.  The smallest carrier not
refuted here is

```text
(one irredundant tooth word;
 marked blocker-cycle positions;
 every oriented wall occurrence (a,m;b,n;W,g,w);
 owner-reuse incidence between separated positions).                   (27)
```

It preserves the cover predicate, exact phase landing, and the fibre (25).
The bare Fano/`chi_7` quotient preserves only the true zero-seam statement
and the global density floor; it destroys the positive seam digit and word
location.

## 7. Exact replay and scope

The dependency-free referee checks all

```text
335 * 2 * 128 = 85,760                               (28)
```

rows: all compact rungs, both binary digits, and every seven-label colour
word.  Each row verifies the two sharp endpoint overlaps, gcd/lcm sheet,
tooth order and separation, carrier containment, speed descent, exact next
blocker clock, binary marked-tooth side, and compact ratio margin.  It also
checks the four-row character table (14) and the transitive runner-tournament
fingerprint.  Every proof-critical check uses an explicit `RuntimeError`
guard, and the referee parses its own AST and rejects any remaining Python
`assert` node.  Normal and optimized runs are byte-identical.

The sorry-free Lean core kernel-checks (14), the detuning identity (10), both
wall formulas (12), the binary address split (19), exact next-blocker clock,
and terminal compact margin.  Selection of arbitrary residue-class
representatives and interval/tooth assembly remain the explicit paper and
exact-referee layers; there are no proof placeholders or `native_decide`
calls.

The frozen artifact hashes are

```text
script         bce7538d48df87b78bcf864442335b6156b94bcfa835415749cb0a6dc0dbc709
output         218c23c21150f91ca65b87bbc4e21a835813bbaf2a6688c28a9a61e24dc14482
formalization  786c201095e69658658a2e2e52cca9a31b40540da3c4c332d1229633e70ff432
```

THM-1260 does not prove or disprove a Fano law for an **entire simultaneous
coherent blocker cycle**.  It proves the strongest local no-go presently
available: placement of both wall germs, sharp lcm mass, all finite toothpick
rungs, binary phase order, and one further marked blocker still leaves the
full seven-point speed-colour cube free.  The remaining Fano probe begins at
the owner-reuse fibre product (26), not at another point-colour census.
