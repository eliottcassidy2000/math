# Thirty-six, signed gluings, and exact half-step clocks

Status: **PROVED**, elementary identities and permutation constructions; **FINITE-EXACT**, the stated independent controls. The displayed Collatz cycles are explicit certificates, not an exhaustive all-height classification. No convergence or prime-distribution claim is made. Date: 2026-09-25.

## 1. Recovered claim and inherited boundary

The exact earlier occurrence is [level11_short_20260922.md](level11_short_20260922.md), sections “What the two gluings actually do” and “Radix and diagonal statements”: the six selected representatives `±{1,5,17}` give `6²=36` ordered pairs, with 18 contributed by each triple of row labels. It explicitly distinguishes those representatives from the 20 vertices of the six exhibited odd-only cycles. It also distinguishes a full relation matrix from a tournament's 15 arcs. The [current synthesis](collatz_procgen_20260922_synthesis.md), section 3, summarizes that count as “arithmetic true, no map.”

The construction below adds an actual map involving 36. It uses a **half-step in the ordinary dynamical clock**, where the cycle through 17 has period 18. It does not retroactively identify the earlier 36 ordered pairs with these new 36 states.

Closest proved mechanism: [signed cycles](arithmetic_braids2_20260917_signed_cycles.md), sections 1 and 3, supplies signed conjugacy and the ordered valuation clock. Canonical hostile: the three minus cycles persist under sign reflection, so sign gluing alone does not merge their basins. Corrected near miss: minima `1,5,17` are not cycle lengths or vertices of a three-cycle. Least-used sidecar here: the cycle's clock and the relative phase between two labelled copies. [THM-4472 — four-vertex reading of 3n plus/minus 1](../../01-canon/theorems/THM-4472-four-vertex-reading-of-3n-plus-minus-1.md) further separates sign reflection from time reversal.

Live board: sign conjugacy; odd/shortcut/ordinary clocks; cycle square roots; copy-phase torsors; arithmetic inverse fibres. The last operation belongs to the parallel lane and is explicitly separated in section 6.

## 2. Three clocks and two different gluings

For `sigma=±1`, define

\[
C_\sigma(n)=\begin{cases}n/2&n\text{ even},\\3n+\sigma&n\text{ odd},\end{cases}
\quad
T_\sigma(n)=\begin{cases}n/2&n\text{ even},\\(3n+\sigma)/2&n\text{ odd},\end{cases}
\]

and, on nonzero odd integers,

\[
U_\sigma(n)=\frac{3n+\sigma}{2^{v_2(3n+\sigma)}}.
\]

All three satisfy `P_(-sigma)(-n)=-P_sigma(n)`. This is conjugacy preserving forward time, with both the parameter and the state negated. It does not equate the two positive systems.

For any of these clocks, write `P_sigma` for its restriction to positive states. The signed gluing

\[
G_\sigma(\epsilon m)=\epsilon P_\sigma(m),\qquad m>0,\quad\epsilon\in\{\pm1\},
\tag{1}
\]

is two reflected copies of **one** positive system. In the ordinary clock its odd rule is `3n+sigma*sgn(n)`.

- `G_+` uses `3n+1` on positives and `3n-1` on negatives: two positive-plus copies.
- `G_-` uses `3n-1` on positives and `3n+1` on negatives: two positive-minus copies, with the six exhibited cycles under discussion.

The actual signed map `C_+` has a different assignment: its negative half reflects the positive-minus system. These three choices should not be conflated.

If an odd-only cycle has `L` odd states and valuation word `(k_1,...,k_L)`, put `K=sum k_i`. Its shortcut period is `K`; its ordinary period is `K+L`. Each ordinary odd step contributes one multiplication and then `k_i` halvings; the shortcut combines that multiplication with the first halving.

| Minimum | Odd states | Valuation word | Odd period `L` | Shortcut `K` | Ordinary `L+K` |
|---|---|---|---:|---:|---:|
| 1 | `(1)` | `(1)` | 1 | 1 | 2 |
| 5 | `(5,7)` | `(1,2)` | 2 | 3 | 5 |
| 17 | `(17,25,37,55,41,61,91)` | `(1,1,1,2,1,1,4)` | 7 | 11 | 18 |

Thus the displayed positive-minus locus has respectively 10, 15, or 25 states; its reflected gluing has 20, 30, or 50. None of these whole-locus counts is 36.

## 3. A global half-step map with two visible clock sections

**S1 (PROVED).** For any map `P` from a positive state set into itself, define on two signed copies

\[
R(m)=-m,\qquad R(-m)=P(m)\quad(m>0).
\tag{2}
\]

Then `R²=G`, where `G(epsilon*m)=epsilon*P(m)`. Indeed `R²(m)=P(m)` and `R²(-m)=-P(m)`. Equivalently, on `(m,phase)` the map is `(m,0)->(m,1)->(P(m),0)`.

This is an exact global construction for each of the positive Collatz maps, without assuming convergence. Every `P`-cycle of period `m` becomes an `R`-cycle of period `2m`. Under `R²`, its positive and negative sections are two separate copies of the original cycle. A trajectory has a periodic tail under `R` exactly when its projected trajectory has one under `P`; the extra clock therefore preserves the unresolved convergence obligation.

For ordinary minus dynamics the 18-cycle is

```text
17,50,25,74,37,110,55,164,82,41,122,61,182,91,272,136,68,34.
```

Equation (2) interleaves it with its negative copy:

```text
17,-17,50,-50,25,-25,...,34,-34,17.
```

This is a genuine **36-cycle of the specified half-step map**. The other two exhibited cycles have half-step periods 4 and 10, so the entire reflected locus still has 50 states. A sign flip and a clock phase are part of (2); its arrows are not ordinary Collatz arithmetic arrows.

## 4. Why some clocks require duplicate cycles

**S2 (PROVED).** A permutation `P` of a finite set has a permutation square root `R` on the same set if and only if, for each even `m`, the number of `m`-cycles of `P` is even.

To prove necessity, square a cycle of `R`: an odd-length cycle remains one cycle of that length, while a cycle of even length `2m` splits into two `m`-cycles. Thus every even-length `P`-cycle must be paired with another of equal length. Conversely a lone odd cycle has a square root, and any two equal-length cycles can be interleaved, as follows.

For an odd cycle `A_i`, with indices modulo `m` and `P(A_i)=A_(i+1)`, the unique square root that stays in this cycle is

\[
R(A_i)=A_{i+(m+1)/2}.
\tag{3}
\]

For two labelled cycles `A_i,B_i` of the same length, every interleaving square root is, for a unique `a mod m`,

\[
R(A_i)=B_{i+a},\qquad R(B_j)=A_{j+1-a}.
\tag{4}
\]

There are exactly `m` choices. Completeness follows also from `RP=PR`, which forces a fixed shift along each source cycle; squaring forces the two shifts to sum to one. The relative phases form a torsor for `Z/mZ`: a group acts freely and transitively, while an origin depends on the selected markings.

If adding disjoint periodic cycles is allowed while retaining the original permutation unchanged, the smallest square-root completion adds exactly one extra `m`-cycle for every even `m` occurring an odd number of times. No smaller addition can repair that multiplicity.

| Original displayed locus | Required added cycles | Minimal completed size | Number of roots on that labelled completion |
|---|---|---:|---:|
| Odd-only: lengths `1,2,7` | one 2-cycle | `10+2=12` | 2 |
| Shortcut: lengths `1,3,11` | none | 15 | 1 |
| Ordinary: lengths `2,5,18` | one 2-cycle and one 18-cycle | `25+20=45` | 36 |

These are minimal completions of the **specified finite periodic locus**, not classifications of square roots of the entire Collatz functional graph. Transient branches pose additional extension conditions.

## 5. Thirty-six is also an exact phase-choice count

**S3 (PROVED).** The minimal ordinary completion has cycle type

\[
2^2\,5^1\,18^2
\]

on 45 labelled states. Formula (4) gives 2 choices for the pair of 2-cycles and 18 choices for the pair of 18-cycles. The single 5-cycle has the unique root (3), namely a shift by 3. Therefore this completed permutation has exactly

\[
2\cdot18=36
\]

square roots. Every root has cycle type `(4,5,36)`. Its parameter set is the relative-phase torsor `Z/2Z × Z/18Z`.

The full reflected 50-state ordinary locus has type `2² 5² 18²`. Its two 5-cycles can remain separate, with one choice, or interleave, with five choices. Consequently it has exactly `36*(1+5)=216` roots: 36 with root cycle type `(4,5,5,36)` and 180 with type `(4,10,36)`.

There is a useful gauge boundary. For a paired block let reflection exchange aligned labels, `J(A_i)=B_i`. A root (4) commutes with `J` precisely when `2a=1 mod m`. This is impossible for even `m`. Thus the required roots on paired even cycles necessarily choose a phase that breaks this reflection symmetry. Retaining both copies does not make the half-step choice canonical under every symmetry.

These statements give two specific new meanings of 36—a cycle period and a count of relative-phase roots—with their maps and state spaces. They do not identify either with the earlier 36 representative pairs, a 36-state arithmetic residue automaton, or prime statistics.

## 6. Distinction from the arithmetic inverse-fibre half-step

The parallel arithmetic construction uses `(n,sigma)->(2n+sigma,-sigma)`. Its square is the sibling map `n->4n+sigma`, and its cube is `(n,sigma)->(8n+3sigma,-sigma)`. Directly,

\[
3(2n+\sigma)-\sigma=2(3n+\sigma).
\]

After dividing out powers of two, the original and transformed pairs share an odd target. That operation moves through inverse fibres and changes the additive parameter. It is not (2), which is a square root in a chosen forward time clock. Their “half-step” meanings are different and both need their side coordinates. No equality between their cycle periods or their relative-phase choices is inferred.

## 7. Verification and scope

Run

```text
python 04-computation/experiments/thirtysix_signed_20260925.py
python -O 04-computation/experiments/thirtysix_signed_20260925.py
```

The [script](../../04-computation/experiments/thirtysix_signed_20260925.py) uses explicit checks. It directly follows each of the three supplied starts at each clock; verifies signed conjugacy, the two gluings, and (2) for all integers from -2001 through 2001 where defined; checks the valuation clocks; independently enumerates the square-map fibres of every permutation on `n=0..8` (46,234 permutations in total); constructs all 36 and all 216 advertised roots and checks their squares and cycle types; checks every paired reflection phase at lengths 2, 5, and 18; and checks the smaller odd-only and shortcut completions. Its [frozen output](thirtysix_signed_20260925.out) retains every cycle witness. Normal and optimized runs give identical output. An independent agent proof audit passed the square-root criterion, phase formulas, minimal-completion hypotheses, and the counts 36 and 216.

The unrestricted theorems here concern signed conjugacy, clock suspension, and finite permutations. The cycle catalog is only the explicitly exhibited locus. Neither extra copies nor the half-step lift cause a trajectory to reach a root that its original trajectory did not reach.

The independent audit also checked the global identity `R²=G` for arbitrary positive self-maps (including noninvertible ones), and the paired-block reflection obstruction `2a=1 mod m`: **PASS**.
