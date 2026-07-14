---
id: THM-772
title: Folded unit-grid divisor transfer in the two- and three-sheet equality packets
status: PROVED (elementary unit-grid and safe-interval argument; exact shell certificate independently reproduced)
source: codex-2026-07-14-S7
depends_on:
  - THM-769
related:
  - THM-523
  - THM-593
  - THM-765
  - HYP-6820
verification: 04-computation/lrc13_folded_unit_grid_divisor_transfer_codex_S7.py
  (+ 05-knowledge/results/lrc13_folded_unit_grid_divisor_transfer_codex_S7.out)
---

# THM-772 — Folded unit-grid divisor transfer

## Statement

Let `A` be a primitive tight twelve-speed set.

### A. Two-sheet equality packet

Suppose a binding scale gives THM-769's two-exception packet

```text
A=2U union {x,y},       |U|=10,       x,y odd.
```

Then:

1. `U` is primitive;
2. `U` contains a multiple of every `m=2,...,12`;
3. `U` has no multiple of 13;
4. if `B=max(U)`, then

   ```text
   x,y<=11B,             min(x,y)<=11B-36;
   ```

5. with `mu=M(U)` and `rho=(mu-1/13)/B`,

   ```text
   1/(xy)+2rho <= 2/(13x)+2/(13y).                        (A*)
   ```

If `A` is in the deep branch in the global sense that it contains a
13-multiple, at least one of `x,y` is such a multiple.

Thus the quotient core is itself a primitive divisor-complete object, not an
arbitrary loose ten-set.

### B. Three-sheet equality packet

Suppose instead that the three-exception equality edge is

```text
A=3U union {x,y,z},     |U|=9,        3 does not divide xyz.
```

Then:

1. `U` is primitive;
2. `U` contains a multiple of every `m=2,...,11`;
3. `U` has no multiple of 13;
4. if `U` has no multiple of 12, then, after permutation,

   ```text
   {x,y,z} mod 36 = {epsilon_1*2, epsilon_2*10, epsilon_3*14},
   epsilon_i in {+1,-1};
   ```

5. if `B=max(U)`, then `x,y,z<=10B`.

The eight independent sign choices in Part B(4) are all locally feasible at
the unit fractions of denominator 12.  This theorem does not assert that any
extends over the full loose set `G_U`.

## 1. Unit fractions turn a missing divisor into a colour CSP

For `s in {2,3}`, write the equality packet as

```text
A=sU union F,       |F|=s,
G_U={tau:phi_U(tau)>1/13}.
```

By THM-769, at every `tau in G_U` each exception owns exactly one of the `s`
lifts and the `s` colours are all represented.  In particular every exception
is individually eligible:

```text
||w tau||<=s/13.                                           (1)
```

Fix `2<=m<=12` and suppose `U` contains no multiple of `m`.  For every unit
`a mod m` and every `u in U`, multiplication by `a` leaves `u` nonzero modulo
`m`, so

```text
||u a/m||>=1/m>1/13.
```

Hence every unit fraction `a/m` belongs to `G_U`.  If `N_w(a)` is the unique
nearest integer to `aw/m`, (1) becomes

```text
least_absolute_residue(aw mod m) <= floor(sm/13),          (2)
j_w(a)=-w^(-1)N_w(a) mod s,                               (3)
{j_w(a):w in F}=Z/sZ.                                     (4)
```

The inverse in (3) exists because every exception is off-sheet.  The phase
data have exact period `sm` in `w`: adding `sm` changes
`w(a+jm)/(sm)` by the integer `a+jm` on every lift.  Thus (2)--(4) are a finite
residue problem modulo `sm`, with closed eligibility at the endpoint.

## 2. The exact local shell table

The complete table is as follows.  Residues are modulo `sm`; `+/-R` means the
union of the signed classes, and `P` is the number of unordered colour-complete
packets allowing repeated residue classes.

| `m` | eligible `s=2` residues | `P_2` | eligible `s=3` residues | `P_3` |
|---:|---|---:|---|---:|
| 2 | empty | 0 | `+/-2` | 0 |
| 3 | `{3}` | 0 | empty | 0 |
| 4 | empty | 0 | `+/-4` | 0 |
| 5 | `{5}` | 0 | `+/-5` | 0 |
| 6 | empty | 0 | `+/-{1,5,7}` | 8 |
| 7 | `{7}` | 0 | `+/-7` | 0 |
| 8 | empty | 0 | `+/-8` | 0 |
| 9 | `{9}` | 0 | empty | 0 |
| 10 | empty | 0 | `+/-10` | 0 |
| 11 | `{11}` | 0 | `+/-11` | 0 |
| 12 | empty | 0 | `+/-{2,10,14}` | 8 |

This table can be checked directly from (2)--(4).  For `s=2`, all even moduli
have no eligible odd residue.  At an odd modulus the only all-unit eligible
class is `m mod 2m`; both exceptions then have the same colour word, so no
two-colour packet exists.  Consequently absence of any divisor `m=2,...,12`
contradicts persistent two-sheet ownership.  This proves A(2).

For `s=3`, all local packets vanish except at 6 and 12.  If `U` had no
6-multiple, it would also have no 12-multiple and both unit grids would apply.
The denominator-6 eligible residues modulo 18 are all odd,
`+/-{1,5,7}`, whereas the denominator-12 eligible residues modulo 36 are all
even, `+/-{2,10,14}`.  No exception is eligible on both grids, so a joint
packet is impossible.  Thus `U` has a 6-multiple and, by the zero rows of the
table, a multiple of every other `m=2,...,11`.

At modulus 12, the eight colour-complete packets are exactly one independent
sign choice from each of the three pairs `+/-2`, `+/-10`, and `+/-14`.  This
proves B(2) and B(4).  The verification artifact independently enumerates nearest
integers and colours and was reproduced by a second implementation evaluating
the lifted times directly.

## 3. The quotient cores are primitive

Let `d=gcd(U)`.  The loose set `G_U` is nonempty by the settled
lower-dimensional LRC bound, and it is invariant under `tau -> tau+1/d`.
Fix `tau in G_U`.  For an exception `w`, its phases over this `d`-shift orbit
form a `D`-grid, where

```text
D=d/gcd(d,w).
```

By individual eligibility, the whole grid lies in the single closed arc
`||z||<=s/13`, of length `2s/13<1/2`.  A circular `D`-grid with `D>=2` needs an
arc of length at least `1/2`.  Therefore `D=1` and `d|w` for every exception.
It follows that `d` divides every speed of `A=sU union F`.  Primitivity of `A`
forces `d=1`, proving A(1) and B(1).

At the chosen binding point, THM-769 also makes `p/13` a local maximum of
`phi_U` of height exactly `1/13`.  A 13-multiple in `U` would have clearance
zero there.  This proves A(3) and B(3).  If `A` globally belongs to the deep
branch, some speed of `A` is a 13-multiple; since none lies in `sU`, it must be
an exception.  Existence of an `s=2` maximizer alone does not rule out
additional shallow maximizers, so the explicit deep hypothesis is needed for
this last conclusion.

## 4. Safe-component width bounds every exception

Put `B=max(U)` and `mu=M(U)`.  The `B`-Lipschitz property of `phi_U` gives,
around a maximizer, an open connected arc `I subset G_U` of length at least

```text
2(mu-1/13)/B.                                             (5)
```

For a core with `12-s` speeds, settled LRC gives

```text
mu>=1/(13-s).
```

Every exception is eligible throughout `I`.  Its eligibility teeth are
separated closed intervals of length `2s/(13w)`, so the connected arc `I`
lies in one tooth.  Comparing lengths in (5) gives

```text
2(1/(13-s)-1/13)/B <= 2s/(13w),
w <= (13-s)B.                                             (6)
```

For `s=2` this is `w<=11B`; for `s=3` it is `w<=10B`, proving A(4)'s first
bound and B(5).  The strict sign in `G_U` and closed sign in the eligibility tooth are
why the non-strict conclusion in (6) is valid.

### Two further taxes in the two-sheet packet

Lift the safe arc of radius `rho=(mu-1/13)/B` around a core maximizer to the
real line.  The containing `x`- and `y`-teeth have centres `k/x` and `l/y`,
where `k,l` have opposite parity.  Since `x,y` are odd, `ky-lx` is odd and
nonzero.  The intersection of the two teeth contains an interval of radius
`rho`, so

```text
1/(xy) <= |k/x-l/y|
        <= 2/(13x)+2/(13y)-2rho.
```

This is (A*).  In particular, using `rho>=2/(143B)`,

```text
13+52xy/(143B) <= 2(x+y).                                (7)
```

For the sharper one-exception bound, choose a global pair-event maximizer
`p/q` of `U` **in lowest terms**.  The exact rational-maximizer theorem gives
`q<=2B`.  If `q` divided both odd exceptions, then `q` would be odd and the
nearest integers `px/q` and `py/q` would have the same parity, contrary to
opposite sheet ownership.  Hence some `w in {x,y}` is not divisible by `q`.
Because `gcd(p,q)=1`,

```text
||wp/q||>=1/q.
```

Containment of the radius-`rho` safe arc in the `w`-tooth yields

```text
1/q+w rho <= 2/13.
```

Using `q<=2B` and `rho>=2/(143B)` gives

```text
w<=11B-143/4,
```

and integrality proves `min(x,y)<=11B-36`.  Reduction of `p/q` is essential;
the phase gap `1/q` need not hold for an unreduced pair-sum presentation.

## Tournament / information-preservation audit

The exact carrier above has vertices `(m,a,j)`: a divisor modulus, a unit
fraction, and a sheet-colour obligation.  The verification also forms a
tournament on the 22 coarser modulus obligations.  Its pairwise observable is
local obstruction strength; switching from raw eligible-class counts to
eligible density flips 43 edges.  Both gauges are transitive (score histogram
`0,...,21`, no directed cycles, singleton SCCs, one Hamiltonian path after the
fixed tie rule).

This explicitly challenges the assumption that runners should be vertices.
The modulus tournament is telemetry only: it destroys the unit-fraction
columns and, in particular, cannot see the joint denominator-6/12 splice that
proves the 6-divisor.  The obligation hypergraph preserves the LRC predicate.

## Frontier effect

The last two-exception deep branch is now narrower than “ten even plus two
odd.”  Its quotient is a **primitive ten-speed set containing multiples of all
eleven moduli 2 through 12**, and both odd tighteners lie below eleven times
its maximum.  The three-sheet equality edge has the analogous primitive
divisor transfer and one explicit modulus-36 residual shell.  Eliminating
these persistent covers uniformly remains open.
