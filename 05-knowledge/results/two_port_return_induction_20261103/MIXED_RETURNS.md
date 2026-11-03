# Mixed returns: positive growth defeats exact ternary reset cycles

**Status: PROVED four-macro return obstruction and a separate seven-macro
extension; FINITE-EXACT endpoint and actual-source inverse controls. A closed
constructor and the unrestricted Hamilton-count spectrum remain OPEN.** Input commit:
`aa7514b988295fe0caac7fc06359a28b330dccf3`. No external priority claim is made.

The two newly inherited singular return channels do change the 3-adic bulk
level. Mixing them with b and q2 nevertheless cannot close a finite bank of
fixed seed interfaces under exact ternary resets. The obstruction is an
integer normalized height, not a finite residue scan. It applies to arbitrary
return lengths and arbitrary actual seeds in the stated model.

There is also a complete short reverse algorithm for the *physical integer
bulk candidates* of this return language. Its four projective image intervals
are disjoint, so each reverse step has at most one possible last macro.
Actual source realization remains a separate obligation.

## 1. Inheritance, maps and the live concept board

Closest positive mechanisms:

- [RETURN_BANKS, Theorem R](../adaptive_state_induction_20261103/RETURN_BANKS.md)
  proves acyclicity when all return quotients are 3-units. Its exact edge model
  is retained below. Its actual q5-to-b/q2 channel switch remains a possible
  escape from a *fixed* correction line, not a closed return bank.
- [Adaptive INCOMING I2](../adaptive_state_induction_20261103/INCOMING.md)
  gives `Q_s78^2=3U`, `det U=49`, and the mod-49 row-lattice gate.
- [Return-channel INCOMING I2](../return_channels_20261103/INCOMING.md)
  gives `D=L2 b L2`, `L2=r24 a a a`, `Q_D=3R`, `det R=-100`, and the
  mod-100 gate. D costs 19 vertices. The old unit-quotient hypothesis fails
  for these channels; its conclusion must not simply be reused.
- [Ternary-repair THEORY R3](../ternary_repair_clock_20261103/THEORY.md)
  supplies four compact actual repairs and 30 actual endpoints through four
  levels, with nonzero carries. These are useful test objects, not a closed
  exact reset language.

Canonical hostile: the actual H31/95/97 seeds have bulk on the same ray but
fail both individual lattice gates. Corrected near miss: an integral
quotient inverse or a residue-clock repair is not an actual source inverse.
Least-used sidecar: the ordered positive cone of the saturated bulk, together
with common 3-content rather than just its residue.

| Live concept | Initial hypothesis / cheap decisive test | Outcome |
|---|---|---|
| Mixed singular returns, anchor | Combine the mod-49 and mod-100 channels | A cone separation excludes all lengths into the original seeds |
| Positive bulk geometry | Retain the ratio lost by valuation | Four disjoint last-macro intervals give a reverse decoder |
| Content versus growth, niche | Normalize away every factor of 3 | A positive integer height proves all finite exact banks acyclic |
| Repaired endpoint carries | Move to actual off-ray representatives | Compact endpoints still reject; some deeper endpoints pass bulk tests only |
| Actual source inverse | Test each surviving full state in a declared bank | Twelve candidate states all absent in the checked 66-letter bank |
| Changed channel/native moves, wildcard | Leave the exact one-line reset model | Still a live escape; no identification with the two-port alphabet is assumed |

The perspective generator was `actual marked sources × saturated bulk ×
content/ratio/height × mixed returns/carries × correction-line quotient ×
unbounded word length`. Each round below changes what the next one tests.

Conventions are inherited from
[inverse65 INVERSE I1](../inverse65_zeckendorf_20261103/INVERSE.md): states
are row vectors `v=(h,x,y)`, words act on the right in **processing order**,
and literal source-component order is reversed. The source word and vertex
port are retained separately. The quotient map is

```text
tau=(2,-1,0),
B=[[0,1],[0,2],[1,-1]],
beta=vB=(Y,X)=(y,h+2x-y),       P B=B Q.
```

This map preserves exact returns modulo the error line and destroys the
charge along tau, adjacency and source membership. B is saturated: its first
and third rows have determinant -1. Actual nonempty states obey
`h>=1`, `0<=x<=y<=h+x`, and `h<=x+y`; the physical proof uses
`(H-E,E,C)=(h+x-y,y-x,x+y)` and the injection `C>=H` in the same source.
Consequently

```text
Y>0,  X=(h+x-y)+x>=0,  X<=3Y.                       (M1)
```

The formal empty unit `(1,0,0)` has bulk `(0,1)` and is not an actual seed.
We do not impose the old bank-specific upper bound `x,y<=16h` in (M1).

## 2. Round one: a cone test beats the individual lattice gates

The four allowed macros and their exact quotient data are:

| Macro | Q | det Q | tau gain | Y-growth lower bound c |
|---|---|---:|---:|---:|
| b | `[[5,8],[2,3]]` | -1 | -1 | 5 |
| q2 | `[[11,18],[4,7]]` | 5 | 1 | 11 |
| s78 | `[[77,152],[28,55]]` | -21 | 1 | 77 |
| D | `[[715050,1009641],[268050,378483]]` | -900 | -4 | 715050 |

For a nonnegative row beta, the first coordinate of beta Q is at least cY.
All four matrices are strictly positive and preserve the cone in (M1).

The actual processing seeds are

```text
U =r5 a a,             v=(31,14,25), beta=(25,34);
U2=a a q1 a a,         v=(95,41,75), beta=(75,102);
U4=q0 a a a a a,       v=(97,40,75), beta=(75,102).
```

Let `L(Y,X)=25X-34Y`. The coefficients of `L(beta Q)` are respectively
`(30,7)`, `(76,39)`, `(1182,423)`, and `(929325,348375)`. They are all
strictly positive. Every nonzero nonnegative bulk sent through a final macro
therefore has L>0; all three target seed bulks have L=0.

**Proposition M2.** For v equal to any of these three states, there is no
actual state w, integer m and nonempty word r over `{b,q2,s78,D}` satisfying

```text
w P_r = 3v + m tau.                                  (M2)
```

Multiply (M2) by B and use the final macro to obtain the contradiction.
This excludes all word lengths before integrality or source membership is
tested. The separator even excludes a formal-empty parent. Empty returns
are not excluded: `v(U2)=3v(U)+tau` and `v(U4)=3v(U)+2tau` are inherited
actual positive controls.

First failed implication: *a singular return can absorb a factor of 3, so
mixing singular returns should permit an exact reset*. It can absorb that
factor while sending the bulk into the wrong projective region. The surviving
statement is the inherited level-change law, with its altered normalized bulk.

## 3. Round two: normalized height forbids every finite exact bank cycle

For a nonzero integer row beta define

```text
ell(beta)=min(v_3(beta_i): beta_i!=0),
K(beta)=Y/3^ell(beta).                                (M3)
```

On the physical cone, K is a positive integer. It is unchanged by multiplying
beta by 3. For any nonsingular integer matrix Q,

```text
0 <= ell(beta Q)-ell(beta) <= v_3(det Q).              (M4)
```

The lower bound follows from integrality. For the upper bound, use
`(beta Q) adj(Q)=(det Q) beta`: the left side has common valuation at least
`ell(beta Q)` and the right side has valuation `ell(beta)+v_3(det Q)`.

For the four macros, `(c,v_3(det Q))` is
`(5,0),(11,0),(77,1),(715050,2)`. Thus

```text
K(beta Q) >= [c/3^v_3(det Q)] K(beta) > K(beta),       (M5)
```

with certified lower factors `5,11,77/3,79450`. Every nonempty mixed return
strictly increases K in its *forward matrix* direction.

**Theorem M3 (exact fixed-seed return banks).** Let every interface i have an
actual nonempty fixed seed state v_i. Permit edges i->j with arbitrary fixed
words r in the four-macro alphabet and arbitrary integer m, provided

```text
v_j P_r=3v_i+m tau.                                  (M6)
```

If r is nonempty, `K(beta_j)<K(beta_i)`. If r is empty, K is unchanged and
`ell(beta_j)=ell(beta_i)+1`. Consequently every finite such bank is acyclic,
even if its reset words have arbitrarily large lengths. Along any path,
the number of nonempty returns is at most `K(beta_start)-1`, even without
assuming the bank is finite. Empty resets alone can use an unbounded family
of seed states; the theorem does not exclude that different model.

**Proof.** Quotient (M6) gives `beta_j Q_r=3 beta_i`. Apply (M5) and the
scaling invariance of K. For an empty word use common valuation directly.
A cycle with a nonempty return strictly lowers K; a cycle of empty returns
strictly raises ell. Neither can close. The count bound is integer descent.
Continuation-language compatibility is still needed to use any edge for
construction; omitting it only makes this obstruction stronger.

An independent check uses eigenvalues. A k-edge cycle gives
`beta Q_word=3^k beta`, multiplying the edge return matrices in reversed
cyclic order. If the total word is nonempty, positive first-coordinate growth
gives `3^k>=product(c)>3^v_3(det Q_word)`. But an integral 2-by-2 matrix with
integer eigenvalue `3^k` has other eigenvalue `trace(Q_word)-3^k`, an integer;
its nonzero determinant is divisible by `3^k`, a contradiction. If the total
word is empty, `beta=3^k beta` is impossible. The physical Y>0 justifies the
growth division; even for a nonzero nonnegative vector with Y=0, strict
positivity would rule out such an eigenrelation immediately.

The general mechanism works at any prime p for nonsingular integer matrices
preserving a cone with positive integer coordinate Y, whenever each generator
has `Y(beta Q)>=cY(beta)` and `c>p^v_p(det Q)`. It is a sufficient obstruction,
not a characterization. Equality or failure of that inequality leaves the
argument silent. The full d-dimensional adjugate proof of (M4) works as well.

This is stronger than the inherited 3-unit theorem on this particular alphabet.
It is **not** a theorem about arbitrary words in the old 65 or enlarged atom
banks: D is one complete macro, and cutting it apart allows intermediate
error lines. Nor does it cover carries, source-dependent encodings, changing
tau, two-port births/discards or unrestricted native transformations.

**Corollary M3a (the other three named incoming loops).** The same acyclicity
and nonempty-edge bound hold for the enlarged seven-macro alphabet
`{b,q2,s78,D,L6,L45,L13}`. The additional actual processing words are the
three loops in [return-channel THEORY T4](../return_channels_20261103/THEORY.md):

| Macro / actual processing word | Q | det Q | tau gain | c / 3^v_3(det Q) |
|---|---|---:|---:|---:|
| L6 = r24 q4 | `[[720,1101],[270,413]]` | 90 | 6 | 80 |
| L45 = r19 r32 r48 q4 | `[[1020432,1560769],[366879,561148]]` | 5985 | 45 | 340144/3 |
| L13 = r19 r19 a a r18 | `[[211965,378780],[76220,136205]]` | 81225 | 13 | 70655/3 |

Their vertex costs are 9,19,17. Multiplying their inherited calibrated atom
matrices verifies both `tau P=gain*tau` and `PB=BQ` entry by entry. Their
Q matrices are strictly positive, preserve (M1), and have common determinant
valuation 2; the displayed lower height factors are all greater than one.
Thus exactly the proof of (M5)-(M6) applies to arbitrary mixed words in the
seven macros. No finite word search supplies this extension. Neither the
projective decoder nor the endpoint inverse results below use the three
added loops: disjoint image intervals are a separate four-macro property.

## 4. Round three: repaired endpoints and a deterministic bulk inverse

Write rho=X/Y. The images of the physical cone `0<=rho<=3` are disjoint
closed intervals:

| Final macro | Image interval for rho |
|---|---|
| D | `[71503/50640,336547/238350]` |
| b | `[17/11,8/5]` |
| q2 | `[18/11,39/23]` |
| s78 | `[317/161,152/77]` |

These follow by evaluating the fractional-linear action at 0 and 3; its
derivative has the sign of det Q and is never zero. Exact rational comparison
shows the displayed order and strict gaps.

The four compact repairs have bulks `(239,460)`, `(239,432)`, `(281,432)`,
`(260,467)`. Both H293 clock twins have bulks `(218,271)` and `(204,383)`.
All six ratios lie in the gaps. Thus **none admits a nonempty exact return
(M2) from any actual source**, for any integer m, using the four macros.
This does not negate their actual carried replacements: a relation
`u=3v+j tau+7 kappa` has bulk `beta(u)=3beta(v)+7kappa B`, and the extra
bulk term is precisely what (M6) forbids.

For arbitrary target v the following is a complete reverse algorithm for
possible *physical integer bulks* of nonempty returns:

1. Start from `gamma=3beta(v)`. Its ratio determines at most one possible
   final macro. If it is outside all four intervals, stop.
2. Multiply by that Q inverse. Reject a nonintegral result or one outside
   `Y>0, 0<=X<=3Y`; otherwise retain this suffix as a possible bulk return.
3. Continue with the predecessor. Its Y is at most the old Y divided by 5,
   so termination is finite, with logarithmically many steps in Y.

Every actual return supplies each of these predecessors, so none is pruned.
The candidates do not certify an actual parent. For a bulk `(Y,X)`, all
physical full-state candidates can be written

```text
w=(Y+g,(X-g)/2,Y),
max(X-2Y,-X,1-Y) <= g <= floor(X/3),
g = X (mod 2),             Y+X odd.                 (M7)
```

The interval is finite. Thus one can restore the missing charge, test (M2)
exactly, and then invoke an independently justified actual-source membership
procedure for a declared parent language. For any fixed integer m, invert
the full 3-by-3 matrix directly instead. The script records these candidates
without calling them witnesses.

The cheap test of the *published* 30 tree endpoints has only four targets
passing any nonempty physical bulk inverse:

| Target H / exact state | Possible suffixes, in processing order | Parent bulks |
|---|---|---|
| 881 / `(881,692,885)` | b | `(315,540)` |
| 2597 / `(2597,2078,2641)` | b | `(903,1704)` |
| 2615 / `(2615,1992,2578)` | b; q2 b | `(924,1557)`; `(48,99)` |
| 2669 / `(2669,971,1801)` | b | `(651,1074)` |

The other 26 stop outside every last-macro cone. For offsets `m=+/-1,+/-2`,
the first three b suffixes produce twelve physical integer full states; all
twelve are absent in the **selected 65 plus s78** actual-source language.
The q2 b candidate at H2615 and the b candidate at H2669 have no physical
integer full-state inverse at those four offsets. These are finite candidate
tests with an all-word membership algorithm, not global absence claims.

The script rechecks the inherited inverse's C16 invariance, physical cone,
strict norm `5P(6,9,5)^t>=9(6,9,5)^t`, H monotonicity and nonsingularity for
all 66 letters before using those guards. Extending an old 65-bank negative
verdict without those checks would be invalid. No negative result here
classifies the entire 31733-letter repair bank or arbitrary marked tournaments.

Positive controls matter. For any macro word r, put `v=v(U r)` and `w=v(U2)`.
Then `wP_r=3v+m tau`, where m is the product of the tau gains in r. The
receipt tests r=b,q2,s78,b q2,s78 b,D and recovers their actual parent bulk.
For example `v(U b q2)=(3325,2797,3331)` and the return from U2 has m=-1.
These are inherited transported seed identities, not new closure or coverage.
They demonstrate that the decoder and obstruction leave genuine nonempty
acyclic edges alive.

## 5. What remains live, and reproduction

The failed inference now has two precise forms: paying the valuation cost
does not pay normalized positive height, and obtaining an integer bulk inverse
does not realize its charge as an actual source. The strongest survivor is a
constructive supply of individual transported ternary identities and finite
carried repairs, with explicit source words and clocks.

The next proof interface must explain how a legal operation changes something
excluded by (M6): a bulk carry, correction line, unbounded seed parameter, or
source-dependent reencoding. A carry must be propagated as the full vector
`kappa P_w`, not merely its residue or first coordinate. The common two-port
observer supplies a different exact continuation language; a map from a
golden source replacement to its named-port commands must still be proved.
No equality of numerical dimensions or a native endpoint conversion supplies
that map automatically.

Reproduce from this directory:

```sh
python3 mixed_returns.py > mixed_returns.json
python3 -O mixed_returns.py > /tmp/mixed_returns_optimized.json
cmp mixed_returns.json /tmp/mixed_returns_optimized.json
```

The script imports no producer. It reads the hash-pinned 960-chart atlas and
the existing repair receipt, multiplies exact integer matrices, uses exact
rational inverses and interval comparisons, and tests 39 named target entries
(with repeated published representatives retained). It performs no H-support
census. Independent adjacency recursion checks all three original seed states
through their H/singleton/double-singleton queries, and all 7! permutations
check the H31 seed. Positive inverse controls retain actual words. Normal and
optimized Python outputs agree; guards do not depend on `assert`. Root's
independent proof audit checked the cycle direction, Y=0 boundary, empty-word
case and integer-eigenvalue argument. A separate spectrum-agent audit checked
the adjugate valuation bound, normalized-height direction, unbounded-bank
bound and empty-reset exception, with no corrections. Input hashes and all
small exact results are in [mixed_returns.json](mixed_returns.json).
