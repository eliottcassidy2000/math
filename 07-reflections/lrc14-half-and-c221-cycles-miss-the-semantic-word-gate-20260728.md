# The half and `C_221` cycles miss the semantic word gate

**Status: VERIFIED-EXACT stopping certificate on the canonical typed row.**
The result below compares the proved THM-2698 half-odometer cycle and the
finite-exact transverse `C_221` design target with the literal THM-2305 words
and THM-2461 prescribed-clock coupling.  It proves that neither displayed
cycle is a semantic word cospan.  It is not a row exclusion and does not prove
LRC(14).

## 1. Inheritance and the exact test

On the canonical row

```text
(H,q1,...,q5,c1,c2,c3)
 =(1,14,27,40,53,66,13,13^3,2*13^5),
```

THM-2701 writes the literal singleton terminal words for source `c1` as

```text
Q_a = A0 intersect D_(c2) intersect D_(c1)^c intersect D_(c3)^c,
Q_b = A0 intersect D_(c3) intersect D_(c1)^c intersect D_(c2)^c.   (1)
```

The unshifted `c1` complement in (1) is semantic, not cosmetic.  THM-2461
identifies the actual source-to-word object as the prescribed-clock span

```text
E_j <- E_j intersect B^(-k_j) Q_(j,sigma) -> Q_(j,sigma),
k_j=lambda_j+1.                                             (2)
```

Here the three clocks are `2,4,6`.  A local rail/root/clock edge does not
establish either endpoint of (2).

## 2. Both delayed cycles fail at the same literal bit

The delayed points of the two carriers are

```text
THM-2698:       11/24 -> 11/24,
transverse C221: 4/17 <-> 13/17.                           (3)
```

At every point in (3), exact arithmetic gives:

```text
guard safe; all five ordinary speeds safe;
c2 dangerous; c3 safe; c1 dangerous.                       (4)
```

The `c1` distances are `1/24` in the half carrier and `1/17` in the
`C_221` carrier, both strictly below `1/14`.  Thus each point is the
source-clock-debt version of target `a`, but none lies in literal `Q_a` or
`Q_b`.  This is an interior factor failure, not an endpoint convention.

The physical event centres fail even earlier as THM-2305 sources.  None of
the four centres belongs to any exclusive `E_j`.  Ordinary six-step dilation
does send them to

```text
11/24, 11/24, 4/17, 13/17,
```

and each of these is the deepest-owner fork

```text
Q_(3,{1,2}).                                               (5)
```

But the corresponding starting centre is not in `E_3`.  Hence (5) is a
terminal label without the branch-specific source leg of (2); it does not
define a cospan.

## 3. The half phase is invisible at every prescribed clock

Let

```text
B0(y)={13y},                 B1(y)={13y+1/2}.              (6)
```

Iteration gives

```text
B1^k(y)={13^k y+(13^k-1)/24}.                              (7)
```

Since `13^2=1 mod24`, the affine constant in (7) is integral for every even
`k`.  In particular,

```text
(13^2-1)/24=7,
(13^4-1)/24=1190,
(13^6-1)/24=201117,

B1^k=B0^k                    for k in {2,4,6}.              (8)
```

Equation (8) is the decisive semantic boundary for THM-2698.  The half-turn
changes odd intermediate states and thereby repairs the raw delayed packet,
but it changes none of the three THM-2305 endpoint maps.  Consequently every
atom-to-word matrix `W_(k_j)(tau,sigma)` built only from the same source and
terminal sets is identical under `B0` and `B1`.  The phase cannot create a
new prescribed-clock span.

## 4. Even the one-step literal `B1` language has no recurrence

The half phase also fails the weaker attempt to use consecutive literal
terminal words as its states.  Two sparse implications suffice.

First, if `s in Q_a`, then `c2 s` is dangerous.  Since `B1^2=B0^2`,

```text
c1 B1^2(s)=13*13^2 s=c2 s mod1.                            (9)
```

Thus state `j+2` is `c1`-dangerous and cannot lie in either literal word.

Second, suppose `s,B1(s) in Q_b` and put `z=B1^5(s)`.  Because

```text
z={13^5s+1/2},
```

the two target-`b` teeth become exactly

```text
D_2(z),                     D_26(z).                       (10)
```

The half shift vanishes after multiplication by the even speeds.  The
half-open contraction lemma used in THM-2701 applies unchanged: `D_2(z)`
puts `z=c+e`, `c in {0,1/2}`, `-1/28<=e<1/28`; `D_26(z)` forces
`-1/(28*13)<=e<1/(28*13)`; hence

```text
|14e|<1/26<1/14,
```

so speed `14` is dangerous at state five.

Every six-letter word is therefore empty.  The exact certificate partition
is

```text
first a at state 0 -> c1 failure at 2:    32 words;
first a at state 1 -> c1 failure at 3:    16 words;
first a at state 2 -> c1 failure at 4:     8 words;
first a at state 3 -> c1 failure at 5:     4 words;
first four states b -> speed-14 at 5:      4 words.        (11)
```

Thus the literal `B1` language is nilpotent by six states and has no
recurrent strongly connected component.  This is not a vacuous one-state
failure: the strict rational point

```text
y=12894291/80000000
```

realizes the four-state word `bbba`, with exact minimum margins

```text
10260037/560000000,
 1272987/560000000,
 6115129/280000000,
 3083243/560000000.                                      (12)
```

The companion does not decide whether a five-state `B1` word exists, so the
claim is the proved six-state zero, not an asserted sharp index.

## 5. THM-2644 sees only its two sharp hostile boundaries

The natural `C_221` phase labels of the transverse cycle are

```text
m0=114,                   m1=107=-114 mod221.             (13)
```

If these labels are made into one nonnegative group-algebra weight
`mu=delta_114+delta_107`, then

```text
M=2, E=2, delta=M^2-E=2, R=(mu*mu)(0)=2, mu(0)^2=0.       (14)
```

Thus `R=delta` exactly: this is THM-2644's inverse-pair equality case, not
its strict fixed-branch gate.  A single phase `delta_114` is pure but does
not return, since `2*114=7 mod221`.  Moreover the actual stalk update is
`nu -> 13nu+m`; multiplication by `13` is not invertible modulo `221`.
So (14) is only a projected backtracking audit, not a physical convolution
transition on a common odd torsor.

For THM-2698, both edges project to the unique nontrivial `C_2` label.  On
that quotient `A=delta_1` is pure and `A^2` returns, but

```text
1=-1 mod2,                    A*=A.                        (15)
```

This is precisely THM-2644's sharp even-group hostile: same-oriented return
and reverse/backtracking are indistinguishable, and the returning branch is
not the identity.  The full affine edges are chronological, but use two
state-dependent translations and still lack the lawful common semantic
middle fibre required by THM-2644.

## 6. Consequence and next target

The comparison leaves one clean diagram:

```text
local half/C221 packet cycle
  -> positive rail/present/root/unit/clock support
  -> delayed target-a tooth with an exact c1 debt
  -/-> literal THM-2305 endpoint
  -/-> prescribed-clock atom-to-word cospan
  -/-> strict odd-torsor return gate.                       (16)
```

The half phase is therefore not the missing semantic operation.  The
remaining constructive exits are the same typed ones exposed by THM-2701:

1. build a common-ancestry mapping cone which absorbs the guard/source-clock
   debt and maps it to a target-active current; or
2. re-root at an actual exclusive owner and rebuild the rail/present/root/unit
   packet at its prescribed clock.

Another local periodic-point search in the unchanged literal endpoint
category cannot produce recurrence, by (11).

## Reproduction

```bash
python3 04-computation/lrc14_phase_cycle_semantic_gate_probe_20260728.py
python3 -O 04-computation/lrc14_phase_cycle_semantic_gate_probe_20260728.py
```

Both runs byte-match

```text
05-knowledge/results/lrc14_phase_cycle_semantic_gate_probe_20260728.out
```

SHA-256:

```text
script  a34618cf8d2d7266db44c750bb56ef4c40922888f92ec0f31a4619049b09872a
output  e6e5c6db41aad1e20c977eda3499875b34ede56c56a35b281eb79e2c068714d6
```

The companion uses `Fraction` arithmetic and explicit optimized-mode guards.
Normal and optimized transcripts are identical.
