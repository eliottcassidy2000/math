---
id: THM-1243
title: THE RESONANT TOOTHPICK RAY HAS A UNIFORM GLOBAL REROUTE — the local seven-crack blocker star is punctured by an explicit parity phase of depth greater than 3/28
status: PROVED (all m>=27 explicit thirteen-speed certificate; exact parity residue ledgers and 3/28 depth; alternate-cell separation; dependency-free exact referee; sorry-free Lean arithmetic core)
source: codex-2026-07-19-S78 continuation
depends_on: [THM-1239]
related: [THM-1236, THM-1240, THM-1242, MISTAKE-185]
script: 04-computation/lrc14_resonant_toothpick_global_reroute_thm1243.py
output: 05-knowledge/results/lrc14_resonant_toothpick_global_reroute_thm1243.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCResonantToothpickGlobalReroute.lean
script_sha256: d513823df4322095815f7bec360900785a9a24e3f757fc12fabc4d48b14242ba
output_sha256: 4c371aed42014f5331c197dfb2b05619ba0570fb1a5cbc2b40e32d474e9d39e8
formalization_sha256: dad6863f46cdaccac90303bf4591f5ad3dccc25bfe8aa421e5462631c5937521
---

# THM-1243 — resonant toothpick ray global reroute

## 1. The uniform global certificate

For `m>=27`, put

```text
a=7m+1,
V_m={1,2,3,4} union {a,a+1,...,a+7} union {14a}.       (1)
```

This is a set of thirteen distinct positive speeds.  Define

```text
q=14m+9=2a+7,
p=3m+4+(m mod 2),
s=p/2=floor((3m+5)/2),
t=p/q.                                                  (2)
```

Then the exact lonely depth is

```text
min_(v in V_m) ||vt||=s/q>3/28>1/14.                  (3)
```

Thus every packet in (1) satisfies LRC(14) with a uniform margin greater
than `1/28` above the target.  In particular this dispatches the complete
infinite `m>=42` ray used in THM-1239: although the resonant speed `14a`
strictly covers all seven curvature cracks in one selected `a`-gap, the
whole packet has the explicit global witness (2).

## 2. The master-clock word

The numerator `p` is even in both parity branches.  Since `q=2a+7`, for
`r=0,...,7` one has

```text
(a+r)p
 =((q-7)/2+r)2s
 ==(2r-7)s                         (mod q).             (4)
```

The eight consecutive speeds therefore become the symmetric odd-multiplier
word

```text
-7s,-5s,-3s,-s,s,3s,5s,7s.                            (5)
```

The four small speeds give `2s,4s,6s,8s`, while the resonant blocker gives

```text
14ap== -98s                         (mod q).             (6)
```

Equations (4)--(6) retain the exact kernel and its master period.  They do
not use a low-height resonance truncation, so the finite-interval inverse
failure in MISTAKE-185 is irrelevant.

## 3. Exact parity ledgers

Suppose first that `m=2h`, so `h>=14`.  Then

```text
(q,p,s)=(28h+9,6h+4,3h+2).                            (7)
```

In the speed order of (1), the thirteen least residue numerators are

```text
6h+4, 12h+8, 10h-3, 4h-7,
7h-5, 13h-1, 9h+6, 3h+2,
3h+2, 9h+6, 13h-1, 7h-5,
14h-97.                                                (8)
```

For example, the last entry follows from

```text
98s=11q-(14h-97),                                     (9)
```

and the displayed remainder is between zero and `q/2`.  Direct comparison
shows that every entry in (8) is at least `3h+2`, with equality at `a+3`
and `a+4`.

If `m=2h+1`, then `h>=13` and

```text
(q,p,s)=(28h+23,6h+8,3h+4).                           (10)
```

The corresponding ledger is

```text
6h+8, 12h+16, 10h-1, 4h-9,
7h-5, 13h+3, 9h+12, 3h+4,
3h+4, 9h+12, 13h+3, 7h-5,
14h-139.                                               (11)
```

Here

```text
98s=11q-(14h-139).                                    (12)
```

Every entry in (11) is at least `3h+4`.  At the boundary `m=27`, the small
speed `4` and blocker `14a` also tie; afterward only `a+3,a+4` minimize.
This proves the exact equality in (3).

The depth estimate is not asymptotic.  Its cleared numerator has the fixed
parity values

```text
28s-3q=29                 if m is even,
28s-3q=43                 if m is odd.                 (13)
```

Hence `s/q>3/28`, proving all of (3).

## 4. The reroute really changes address cells

THM-1239's erased carrier gap is

```text
G_m(a)=[(14m+1)/(14a),(14m+13)/(14a)].                (14)
```

For `m>=27`,

```text
right endpoint of G_m(a)<1/6,
t=p/q>1/5.                                             (15)
```

The first inequality clears to `64<14m`; the second clears to
`0<m+11+5(m mod 2)`.  Thus (2) is not a point secretly surviving inside the
locally erased crack complex.  It moves by a macroscopic amount to another
carrier address cell.  This is the precise global reroute demanded by the
guardrail in THM-1239.

## 5. Structural, Kakeya, and tournament audit

The toothpick self-similarity has two different clocks.  On the selected
Kakeya needle, multiplication by `14a` reproduces the seven shrinking cracks
at odd integer addresses.  Globally, the conjugate clock `q=2a+7` turns the
same consecutive packet into the antipodal word (5), and parity selects a
section whose central residues stay deeper than `3/28`.  Local replication
therefore does not imply global trapping: the self-similar stalk has an
explicit transverse escape section.

The quartet `{1,2,3,4}` and the `j=4` continuum/Fano language are not needed
to find the escape.  Their four residues are present in (8)--(11) and are
safe directly.  Likewise THM-1156's `chi_7` seam bipartition does not govern
the local blocker star, because that star consists of strict containments,
not exact two-owner seams.  The useful remnant of the Fano probe is the
lesson that an address atlas must retain metric residue data.

For Tournament Analysis take the pairwise observable to be the difference
between least-residue depths and break ties by speed.  The resulting gauge is
transitive, with score histogram

```text
(0,1,2,...,12),
```

zero directed triangles, singleton SCCs, and one Hamiltonian path.  It loses
the involution `r<->7-r`, the master-clock identity, parity, and the terminal
multiplier `98`; it has no proof power.  We challenged runners, selected
gaps, address cells, wall events, residues, multiplier modes, blocker stars,
and proof obligations as vertices.  The faithful carrier is

```text
(q,s; 2s,4s,6s,8s; +/-s,+/-3s,+/-5s,+/-7s; -98s),    (16)
```

together with parity and the selected-gap address.  It preserves the global
lonely predicate and records exactly what the local crack incidence erased.

## 6. Verification and scope

The dependency-free exact referee checks every `m=27,...,100000`: `99,974`
thirteen-speed residue ledgers, including all `99,959` rows of THM-1239's
one-blocker range.  It checks (3), (8), (11), (13), and the alternate-cell
separation.  Normal and optimized outputs are byte-identical.

The Lean module kernel-checks the master-clock congruence (4), terminal
congruence (6), every unique entry of both parity ledgers, the exact depth
excesses, and the two address-separation inequalities.  Identifying those
positive representatives with the circle norms is the explicit paper layer;
there are no proof placeholders or `native_decide` calls.

Frozen hashes are

```text
source         d513823df4322095815f7bec360900785a9a24e3f757fc12fabc4d48b14242ba
output         4c371aed42014f5331c197dfb2b05619ba0570fb1a5cbc2b40e32d474e9d39e8
formalization  dad6863f46cdaccac90303bf4591f5ad3dccc25bfe8aa421e5462631c5937521
```

THM-1243 closes one exact resonant address ray, not arbitrary blocker words,
six-comb slow-gap coverage, or LRC(14).  The next global target is to prove
that every high-degree blocker address word admits a comparable transverse
clock, or else force incompatible clocks around THM-1240's blocker cycle.

