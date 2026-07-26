---
id: THM-2417
title: "Middle-depth-two relaxed orbit dual obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Under
  THM-2391's primitive final-lane scalar-cover hypotheses, the
  septimal top depth M=2 is impossible. On a positive-measure family
  of generic 343-orbits both C_3 and c_3 are safe. After guard gauge,
  the top word is one residue row Q_t and the other six rows contain
  210 nonguard obligations. Every lower q word lies in an exact
  490-mask family. Each physical blocker lies in a 2401-mask
  relaxation obtained by retaining the lawful quotient-blocker slope,
  applying d -> 13^(-1)d modulo 343, and freeing its physical phase.
  Exact nonnegative integer duals exclude all seven top offsets and
  both contiguous/parity quotient partitions. Thus the whole M=2
  alternative is empty, including positive-clean packets; combined
  with THM-2415, every survivor of the last septimal lane has M=1.
  This does not decrement a thirteen-adic profile row, produce a
  canonical target address, close the last lane, or prove LRC(14).
source: codex-2026-07-26-middle-depth-two-dual
depends_on:
  - THM-2391-blocker-caged-septimal-single-layer-address-reduction
related:
  - THM-2392-clean-toothpick-or-bounded-cross-ancestor-cage
  - THM-2405-two-level-septimal-sheet-independence-and-middle-depth-two-cage-elimination
  - THM-2414-thirteen-skew-septimal-word-transport-and-local-stopping-atlas
  - THM-2415-last-lane-septimal-depth-two-cap
  - THM-2420-compositional-thirteen-root-final-septimal-lane-exclusion
  - MISTAKE-264
script: 04-computation/lrc14_middle_depth_two_relaxed_orbit_dual_thm2417.py
output: 05-knowledge/results/lrc14_middle_depth_two_relaxed_orbit_dual_thm2417.out
script_sha256: b84dcdeca3ee6934fa7c4632862096b2bd38e94bf35bdc66f241b01e9a6dc367
output_sha256: 8174eafd729c20aded84ac10db770e97eed86bbe651100fe3258285bae93a2e0
hash_basis: working-tree bytes (LF)
---

# THM-2417 -- the middle-depth-two relaxed orbit dual obstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2391 reduces every survivor of the last septimal lane to one common
primitive lower layer. At depth `M=2`, its quotient blockers have only
two possible top-row set-types. The missing operation is to keep the
full `343`-orbit after transporting those quotient slopes to the
physical blockers.

Doing so gives a finite relaxation, not a finite sampling:

```text
lawful physical M=2 orbit
  -> exact lower-q top containment
  -> exact quotient-blocker support type
  -> lawful blocker step transport d |-> 13^(-1)d
  -> forget every remaining phase correlation
  -> one of fourteen relaxed six-mask cover problems
  -> exact weighted dual contradiction.                         (1)
```

Because the relaxed family contains every physical tuple, its
infeasibility removes the full `M=2` alternative.

## 1. Last-lane hypotheses and a simultaneous-safe orbit

Retain THM-2391's primitive scalar cover

```text
C_H subset
  union_(i=1)^5 D_(q_i)
  union D_(c_1) union D_(c_2) union D_(c_3),                     (2)

D_v={x:||vx||<1/14},                  E_H={x:||Hx||<1/7},
```

where `c_j=13C_j`. In the final septimal lane,

```text
nu_7(q_*)=M,

nu_7(H)=nu_7(q_i)=nu_7(C_1)=nu_7(C_2)
       =nu_7(c_1)=nu_7(c_2)=0,                 q_i!=q_*,

nu_7(C_3)=nu_7(c_3)>M.                             (3)
```

Assume for contradiction that

```text
M=2.                                                   (4)
```

Put `N=7^3=343` and consider additive orbits

```text
O_z={z+n/N:n in Z/NZ}.                                (5)
```

Both `C_3` and `c_3` are divisible by `343`, so their danger states are
constant on every orbit. Each danger comb has measure `1/7`; hence the
set on which both are safe has measure at least

```text
1-1/7-1/7=5/7>0.                                      (6)
```

Choose an orbit in this simultaneous-safe set after deleting:

- the scalar cover's null exceptional set and its `343` translates;
- THM-2391's null exceptional set and its `343` translates; and
- every strict comb endpoint used on the orbit.

The remaining base set still has positive measure. Thus all following
word identities are literal pointwise identities on the chosen orbit.
This is a physical-root orbit; no base frequency is divided by seven,
so the base/root error in MISTAKE-264 does not occur.

## 2. Exact orbit words

Write

```text
Z_343=Z/343Z,

P(a,d)={a+jd:0<=j<49} subset Z_343,       7 does not divide d.
                                                               (7)
```

Index (5) by the `H`-gauge. Since `H` is a seven-unit, this is a
permutation of `Z_343`. A further cyclic translation makes the guard
word

```text
G={0,1,...,97}.                                        (8)
```

The strict interval defining `E_H` has length `2/7`, so these are
exactly its `98` orbit points.

Because `q_*=49u` with `7` not dividing `u`, its danger word is one
complete residue row

```text
Q_t={t+7j:0<=j<49}                                    (9)
```

for some `t in {0,...,6}`. Every offset can occur, so all seven must be
retained. Since `98=14*7`,

```text
G intersection Q_t={t+7j:0<=j<14}.                   (10)
```

Every lower seven-unit danger word has exactly `49` points. If its
speed is `v`, its step in the `H`-gauge is

```text
d_v=H v^(-1) mod 343,                                 (11)
```

so its word is exactly one `P(a,d_v)`. The base phase supplies `a`.

Define the six-row obligation set

```text
U_t=Z_343 minus (G union Q_t).                         (12)
```

It has

```text
|U_t|=343-98-49+14=210.                               (13)
```

## 3. The relaxed lower-q family

THM-2391 proves on every generic `C_3`-safe top row that every lower
`q_i` word is contained in the guard:

```text
Q_t intersection D_(q_i) subset G intersection Q_t,
                                                    q_i!=q_*.    (14)
```

Accordingly define

```text
A_t={
  P(a,d):
  a in Z_343, 7 does not divide d,
  P(a,d) intersection Q_t subset G
}.                                                       (15)
```

Exact enumeration gives

```text
|A_t|=490                         for every t.          (16)
```

The possible steps modulo `49` are precisely

```text
d mod 49 in {+1,-1,+2,-2},                            (17)
```

recovering THM-2391's four slopes. The family (15) allows the four
labelled lower masks to choose their phases independently and even to
repeat a mask. Both choices enlarge the physical universe.

## 4. Quotient supports and physical blocker relaxation

The two quotient blockers partition (10). THM-2391 proves that their
supports have one of the following two set-types:

```text
contiguous:
  S_(t,0)={t+7j:0<=j<7},
  S_(t,1)={t+7j:7<=j<14};

parity:
  S_(t,0)={t+7j:0<=j<14, j even},
  S_(t,1)={t+7j:0<=j<14, j odd}.                     (18)
```

For either support put

```text
Delta_t(S)={
  d in Z_343^x:
  P(a,d) intersection Q_t=S for some a
}.                                                       (19)
```

These are steps of the **divided** blockers `C_i`. If

```text
d=H C_i^(-1) mod 343,
```

then the corresponding physical blocker `c_i=13C_i` has step

```text
H c_i^(-1)=13^(-1)d mod 343.                          (20)
```

Retain (20), but forget the actual phase transport. Define the physical
relaxation

```text
B_t(S)={
  P(b,13^(-1)d):
  b in Z_343, d in Delta_t(S)
}.                                                       (21)
```

For every `t`, type, and side, exact enumeration gives

```text
|Delta_t(S)|=14,                  |B_t(S)|=2401.       (22)
```

The quotient steps reduce modulo `49` to `+/-1` in the contiguous case
and `+/-2` in the parity case. Translating a support changes its start,
not its step, so the two side families in one type have the same mask
set after the phase in (21) is freed.

This is the load-bearing inclusion:

> The actual physical blocker word has exactly the step (20) and some
> actual start `b`; therefore it lies in (21). Allowing all `b`
> independently discards the thirteen-skew affine digit and all
> common-base phase coupling. It can add masks but cannot remove a
> physical mask.

Thus (15) and (21) form a genuine relaxation of every lawful physical
orbit, not a same-fibre identification of `C_i` with `c_i`.

## 5. Why six masks must cover all 210 obligations

Modulo endpoints,

```text
C_H=E_H^c.
```

Therefore (2) says that the guard together with the eight ordinary
danger masks covers the entire orbit. On `U_t`:

- the guard is absent by (12);
- `q_*` is absent by (9); and
- `c_3` is absent by the simultaneous-safe choice.

Consequently every physical survivor would induce masks

```text
A_1,A_2,A_3,A_4 in A_t,

B_0 in B_t(S_(t,0)),       B_1 in B_t(S_(t,1))       (23)
```

such that

```text
U_t subset A_1 union A_2 union A_3 union A_4
              union B_0 union B_1.                   (24)
```

Notice that `C_3`-safety was used to obtain (14) and (18), while
`c_3`-safety was used in (24). Equation (6) is what licenses both
conditions on one positive-measure orbit family.

## 6. Exact weighted-dual obstruction

For a nonnegative integer weight `w` on `U_t`, put

```text
W=sum_(x in U_t)w(x),

alpha=max_(A in A_t) sum_(x in A intersection U_t)w(x),

beta_i=max_(B in B_t(S_(t,i)))
                 sum_(x in B intersection U_t)w(x).   (25)
```

If (24) held, nonnegativity would give the necessary inequality

```text
W<=4alpha+beta_0+beta_1.                              (26)
```

The exact companion stores integer weights and reconstructs every mask
in (15) and (21) before taking the maxima. It obtains:

| type | `t` | `W` | `alpha` | `beta_0` | `beta_1` | RHS of (26) | margin |
|---|---:|---:|---:|---:|---:|---:|---:|
| contiguous | 0 | 204 | 35 | 30 | 30 | 200 | 4 |
| contiguous | 1 | 204 | 35 | 30 | 30 | 200 | 4 |
| contiguous | 2 | 735 | 127 | 113 | 113 | 734 | 1 |
| contiguous | 3 | 513 | 89 | 78 | 78 | 512 | 1 |
| contiguous | 4 | 735 | 127 | 113 | 113 | 734 | 1 |
| contiguous | 5 | 204 | 35 | 30 | 30 | 200 | 4 |
| contiguous | 6 | 204 | 35 | 30 | 30 | 200 | 4 |
| parity | 0 | 1998 | 339 | 320 | 320 | 1996 | 2 |
| parity | 1 | 9978 | 1702 | 1584 | 1584 | 9976 | 2 |
| parity | 2 | 1765 | 301 | 280 | 280 | 1764 | 1 |
| parity | 3 | 1911 | 326 | 303 | 303 | 1910 | 1 |
| parity | 4 | 1765 | 301 | 280 | 280 | 1764 | 1 |
| parity | 5 | 9978 | 1702 | 1584 | 1584 | 9976 | 2 |
| parity | 6 | 1998 | 339 | 320 | 320 | 1996 | 2 |

Every row strictly reverses (26). Hence none of the fourteen relaxed
cover problems has a solution. By Sections 1--5, no physical `M=2`
orbit exists. This contradicts (4).

We have proved the candidate conclusion

```text
M=2 alternative under THM-2391: empty.                (27)
```

This is the whole alternative, not only THM-2405's no-clean branch and
not only one of THM-2392's orbit outcomes.

## 7. Positive and hostile controls

The exact verifier retains the following boundaries.

1. **Nonvacuous mask geometry.** There are `50,421` distinct unit-step
   `49`-progression masks. Every `A_t` has `490` masks. Both quotient
   support types have `14` oriented steps per side, and every relaxed
   physical blocker family has `2401` masks. Thus the obstruction is
   not caused by an empty support convention.

2. **Full hostile phase freedom.** The physical blocker phase in (21)
   is arbitrary and independent for the two blockers. The verifier
   therefore tests a strict phase relaxation rather than silently
   identifying divided and physical words.

3. **Four-mask boundary.** Replacing `4alpha` by `5alpha` makes every
   stored certificate inequality noncontradictory. This does not assert
   that a five-mask cover exists; it records that the displayed duals
   use the exact four-lower-`q` inventory.

4. **THM-2414 live `W=8` atlas.** Its top speed is `q_*=7`, so
   `M=1`, and its high speed has septimal depth two. It is outside the
   `M=2`, `N=343` universe. Moreover it is a local stopping packet, not
   a scalar cover. THM-2417 therefore does not contradict it.

5. **THM-2414 excluded `W=7` hostile.** Its top speed is `49` and hence
   does have `M=2`, but the verifier reconstructs the two strict
   full-bin failures: at address `s=6` the top and guard are dangerous
   while both divided low blockers and `C_3` are safe; at `s=13` the
   top and off-layer `q_5` are dangerous while the guard and all
   divided blockers are safe. Thus it does not induce a tuple admitted
   to the relaxation because it violates the proved premises (14) and
   (18).

6. **THM-2415 phase-grid hostile.** THM-2415's `M=2` example only shows
   that its particular thirteenth-digit orbit-hitting mechanism stops
   at depth two. It is not a six-mask scalar cover and does not satisfy
   (24). The present obstruction uses the full lower-mask inventory,
   so the two results have compatible sharp boundaries.

## 8. Consequence and remaining frontier

**Current-frontier postscript.** Equation (28) was the sharp consequence
of THM-2417 at promotion time. THM-2420 subsequently excludes this
remaining `M=1` alternative and empties the final septimal lane. The
paragraphs below record THM-2417's standalone boundary, not the current
open frontier.

THM-2415 proves uniformly, including positive-clean packets, that every
survivor has

```text
M in {1,2}.
```

Combining it with (27) gives the structural sharpening

```text
every survivor of the last septimal lane has M=1.      (28)
```

Thus the four-slope/Fano cross of THM-2391 is now an excluded
intermediate object. The live residue is its adjacent-layer binary
address: the two quotient blockers choose the two guard addresses, and
the four lower `q` labels choose among them with repetition. The
thirteen-skew transport and canonical owner/target alignment remain
open there.

No thirteen-adic profile is removed merely by (28); the `165`-row
ledger is unchanged. The theorem does not identify a canonical
expiration target, prove the required all-coordinate address, close
the `M=1` lane, or prove LRC(14).

## 9. Exact companion

Run

```text
python3 04-computation/lrc14_middle_depth_two_relaxed_orbit_dual_thm2417.py
python3 -O 04-computation/lrc14_middle_depth_two_relaxed_orbit_dual_thm2417.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_middle_depth_two_relaxed_orbit_dual_thm2417.out
```

after LF normalization. The dependency-free verifier:

- reconstructs all `50,421` progression masks;
- reconstructs every lower-`q` and transported-blocker family;
- checks all family sizes and the lawful `+/-1,+/-2` slope classes;
- recomputes all fourteen exact dual maxima;
- checks reflection-related offsets independently against the rebuilt
  families; and
- replays the THM-2414 `W=8`/`W=7` boundary controls.

Every executable check uses an optimization-safe `require` function.
No Python `assert` is used.

## 10. Independent hostile audit

An independent set/frozenset implementation, sharing no bitset engine with
the companion, reconstructed:

- the simultaneous `C_3/c_3`-safe `343`-orbit and all translated null-set
  removals;
- the single `H`-gauge words `G,Q_t,U_t`;
- every lower-`q` containment and both quotient support types;
- the same-gauge physical transport

  ```text
  H C_i^(-1) -> H(13C_i)^(-1)=13^(-1)H C_i^(-1) mod 343;
  ```

- all `14` quotient steps and `2401` transported blocker masks per side;
  and
- all fourteen dual maxima in the table of Section 6.

Every independently reconstructed margin is in `{1,2,4}` and is strictly
positive. Normal, optimized, and stored transcripts byte-match, and the
frontmatter hashes agree. The audit also replayed the THM-2414 `W=8/W=7`
and THM-2415 phase-grid hostiles and confirmed that none satisfies the
six-mask scalar-cover premises. No flaw was found.
