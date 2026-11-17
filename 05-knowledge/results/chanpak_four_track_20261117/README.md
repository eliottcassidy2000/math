# Four ordered prime tracks in the noncrossing sevenfold construction

The [six-reversal exterior](../chanpak_increasing_square_20261117/README.md) has read `αP₂+βH` at every step. The seven positive beta coefficients preserve any strict ordering within a fixed two-cover fiber, not just a two-state ordering. Here is a literal four-track witness with two independent internal edits.

Use the ordinal sum A→B→C→singleton as its variable block. In outgoing-mask notation:

| Block | Rows | Optional flip | Hamilton count(s) | Two-cover count |
|---|---|---|---|---:|
| A | (134,124,152,161,169,133,61,66) | (1,4) | 485/491 | 6968 |
| B | (22,68,104,35,238,67,137,47) | (0,2) | 359/383 | 5528 |
| C | (22,84,120,115,96,67,1) | none | 77 | 1108 |

Both H and P₂ multiply under ordinal sum: each labeled path slot simply concatenates its restrictions to the factors in order, permitting an empty restriction. Thus all four settings have P₂=85,358,334,464, while their distinct H values in increasing order are 13,406,855; 13,572,713; 14,303,135; and 14,480,081. The two horizontal changes commute on disjoint arcs; all exterior changes are disjoint from both. The exterior has eight vertices outside the module, giving 32 vertices total. All resulting tournaments are strong by the same contraction argument as for the two-track example.

| Setting (A,B) | Prime start | Sevenfold end |
|---|---:|---:|
| 00 | 5,804,433,777,827 | 40,631,036,444,789 |
| 10 | 5,804,434,607,117 | 40,631,042,249,819 |
| 01 | 5,804,438,259,227 | 40,631,067,814,589 |
| 11 | 5,804,439,143,957 | 40,631,074,007,699 |

Since each state has the same P₂ and each exterior beta is positive, all four counts remain in this strict order at every intermediate stage. Each of the four vertical tracks is increasing by the established all-source inequality. This confirms independence of the two literal internal edits, not arbitrary control over their absolute prime outputs.

## A target-side inverse test, and its limit

For this ordinal grammar an assigned target array p_ij has the form K+5c a_i b_j, with K=68P₂ and c=H(C). Its mixed difference is D=p00+p11−p10−p01=5c(a1−a0)(b1−b0). If the two nonzero moves increase H then D>0. The matrix remaining after subtracting K must have rank one, so K is forced:

    K=(p00*p11−p10*p01)/D.

Here D=55440 and K=5804366743552. A candidate must have K a positive multiple of 68 and all residuals positive; the ratios of adjacent residuals then force the ratios of the inner Hamilton counts, before any physical-source test. A nondegenerate additive rectangle (D=0) cannot come from this particular product grammar. Moreover, for fixed-order fixed-P₂ edits the inner changes are multiples of 6, so necessarily D is divisible by 180. The actual quartet is a constructive example, not a selector for arbitrary quartets.

## Minimum order of a one-arc fixed-P₂ repair

A separate literal enumeration checks every labeled tournament through order 7 and every one-arc neighbor, including weak tournaments: none preserves P₂ while changing H. The order-eight one-arc repair therefore has the least possible order. Two edits do suffice at order seven, for rows (26,28,73,112,100,71,3) and (42,28,73,112,101,70,3), with H=153 and 159 at common P₂=1812. This example preserves the coordinate at the endpoints only.

Replays: `python replay.py` checks all marked block profiles, the offset and the four primes by trial division. Compile `min_order.cpp` with C++20 for the separate minimal-order enumeration.
