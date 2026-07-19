---
id: THM-1240
title: THE CENTERED CARRIER-SPOKE BLOCKER CYCLE — every fast speed has a deep positive sum-beat witness, so every six-cover induces a directed third-support cycle
status: PROVED (all-scale centered pair-sum theorem; deep positive slack and nonzero mixed-clock lift; blocker functional-cycle law; Kakeya-cut clock doublet; dependency-free referee; sorry-free Lean arithmetic core)
source: codex-2026-07-19-S78 continuation with tangent-stalk agent
depends_on: [THM-1217, THM-1219, THM-1241, THM-1236]
related: [THM-1156, THM-1238, THM-1239, HYP-7870]
script: 04-computation/lrc14_centered_carrier_spoke_blocker_cycle_thm1240.py
output: 05-knowledge/results/lrc14_centered_carrier_spoke_blocker_cycle_thm1240.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCCenteredCarrierSpoke.lean
script_sha256: 2248eb59a1548c3bb00cd11a9db3c39d0097cae5218638530b3f8ca9f834296b
output_sha256: 093dea6f05a54893a79095c6da68343af8803866744660b882d9aeae97210394
formalization_sha256: 125aa219863e989898c0e3c4895838a75d473d95587e037ea53928fddcf83a32
---

# THM-1240 — centered carrier-spoke blocker cycle

## 1. The centered carrier-spoke theorem

Let

```text
G=G_k(c)=[(14k+1)/(14c),(14k+13)/(14c)]               (1)
```

be a complete closed safe gap of a positive integer speed `c`.  Its centre is

```text
t0=(k+1/2)/c.                                           (2)
```

For every integer speed `d>c`, put

```text
q=c+d.                                                  (3)
```

Choose an integer `p` nearest `qt0` and set `t=p/q`.  Then

```text
t in int(G),
||ct||=||dt||>=d/[2(c+d)]>1/4.                         (4)
```

Writing

```text
D=q||ct|| in Z_(>0),                                   (5)
```

the active-pair slack obeys the macroscopic lower bound

```text
14D-q>=6d-c>0.                                         (6)
```

Thus every carrier--fast spoke is a deep, strictly positive pair-sum edge.
Unlike the fast-fast edge in THM-1238, it can never fall into the centered
common-zero singleton.

## 2. Proof

Nearest-integer approximation gives

```text
|p-qt0|<=1/2,
|t-t0|<=1/(2q).                                        (7)
```

Since `q=c+d>2c`,

```text
1/(2q)<1/(4c)<3/(7c),                                  (8)
```

and `3/(7c)` is the half-length of `G`.  This proves strict interiority.

The pair-sum identity is

```text
dt=p-ct,                                               (9)
```

so the two circle distances agree.  Moreover `ct0=k+1/2`, and (7) gives

```text
||ct||=1/2-|c(t-t0)|
       >=1/2-c/(2q)
       =d/(2q)>1/4.                                    (10)
```

Multiplying (10) by `q` gives `D>=d/2`; hence

```text
14D-q>=7d-(c+d)=6d-c>0,                                (11)
```

proving (4)--(6).

## 3. The nonzero mixed clock

The reduced spoke period is

```text
Q=(c+d)/gcd(c,d)>=3.                                   (12)
```

Now place `d` among a six-speed packet `d1,...,d6` and set

```text
d0=gcd(q,c,d1,...,d6),             L=q/d0.             (13)
```

Then `Q|L`.  The beat numerator in (7) is nonzero modulo `L`.  Indeed, if
`p=Lr`, then `q=d0L` and `c=d0c'` give

```text
cp/q=c'r in Z,                                         (14)
```

contradicting the `>1/4` distance in (4).  Each carrier spoke therefore
supplies a nonzero master-clock obligation without any primitive-relation or
equidistribution assumption.  In particular it is unaffected by the
nonprimitive beat-localization failure recorded in MISTAKE-182.

## 4. Six-cover implies a blocker cycle

Suppose the strict danger combs of

```text
c<d1<d2<d3<d4<d5<d6                                  (15)
```

cover `G`.  Apply the centered theorem to every `di`.  At its phase `ti`, both
`c` and `di` are safe.  Since the six fast combs cover `G`, at least one label

```text
b(i) in {1,...,6} minus {i}                            (16)
```

is dangerous at `ti`.  Choose one such label for every `i`.  The map

```text
b:{1,...,6}->{1,...,6},              b(i) != i         (17)
```

is a loopless functional digraph.  Iterating from any label visits seven
vertices in a six-element set, so two iterates repeat; taking the first
repeat yields a directed cycle of length between two and six.

This cycle is the first nontransitive relation forced directly by the
remaining slow-gap geometry.  It records third-support incidence at six
different deep phases, not an arbitrary orientation of speed order.

## 5. The Kakeya-cut blocker doublet

Let `j` be the macroscopic path cut from THM-1241:

```text
d_(j+1)-d_j>d6/210.                                    (18)
```

The two adjacent carrier-spoke denominators are

```text
q_-=c+d_j,                  q_+=c+d_(j+1),
q_+-q_-=d_(j+1)-d_j>d6/210.                            (19)
```

Each clock has a nonzero deep active point and a nonempty blocker set among
the other five fast labels.  Hence the cut emits a separated blocker doublet:
either one label blocks both clocks, or two distinct blocker labels occur.
Within the full functional graph, every directed cycle is either contained
on one side of the cut or crosses it in both directions.

This does not yet contradict coverage, but it couples the scale-free Kakeya
cut to phase-bearing mixed clocks and a finite incidence pattern.

## 6. Tournament and carrier audit

On runners, speed order again gives the transitive tournament with score
histogram `(0,1,2,3,4,5)`, no cycles, singleton SCCs, and one Hamiltonian
path.  It destroys the new theorem completely.  On blocker obligations, the
observable is

```text
i -> j  iff j is dangerous at the selected deep spoke phase ti. (20)
```

Choosing one outgoing edge produces the forced functional cycle.  Across all
`5^6=15,625` loopless choices, the exact histogram of the shortest selected
cycle is

```text
length 2: 8265,  length 3: 4360,  length 4: 2160,
length 5: 720,   length 6: 120.                         (21)
```

Completing unspecified blocker pairs to a tournament by speed order is only
a gauge; the forced cycle subgraph is the proof-bearing part.

We challenged runners, fast-fast gaps, carrier boundary roles, pair-sum
beats, master residues, blocker labels, wall events, and proof obligations as
vertices.  The faithful vertices are the six **carrier-spoke obligations**

```text
(qi,pi mod Li,Di,14Di-qi,blocker set),                  (22)
```

with the slow-gap centre and cut index as sidecars.  This carrier preserves
deep pair safety and third-support truth.  It loses coverage between the six
sampled phases, so a closing theorem must combine its cycle with the
positioned mass/address constraints of THM-1237 and THM-1239.

## 7. Verification and scope

The dependency-free exact referee checks `553,945` centered spoke rows,
including negative gap indices and ratios through `d<20c`; every row has
depth `>1/4`, slack at least `6d-c`, reduced period at least three, and a
nonzero residue.  It exhausts all `15,625` loopless blocker maps and checks
`138,606` adjacent spoke-clock separations.  Normal and optimized outputs are
byte-identical.

The Lean module kernel-checks the interior scale, deep-distance and slack
inequalities, nonzero-clock divisibility core, six-label orbit collision,
loopless fixed-point exclusion, and cut-clock identity.  The standard nearest
integer selection and danger-set interpretation remain explicit paper
providers; there are no proof placeholders.

Frozen hashes are

```text
source         2248eb59a1548c3bb00cd11a9db3c39d0097cae5218638530b3f8ca9f834296b
output         093dea6f05a54893a79095c6da68343af8803866744660b882d9aeae97210394
formalization  125aa219863e989898c0e3c4895838a75d473d95587e037ea53928fddcf83a32
```

THM-1240 proves the blocker cycle but not its impossibility.  It does not by
itself prove six-comb noncoverage or LRC(14).
