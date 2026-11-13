# Audit of the H11/H13 fibres and padding theorem

**Status: VERIFIED.** Object: [arithmetic/MODULES.md](../arithmetic/MODULES.md),
Theorems C and D. Auditor: repair lane, 2026-11-13.

The full-fibre reduction is sound. Busch's primary Theorem1 gives the residue-
dependent strong-core lower bound, each branch at least 5^((n−1)/3); therefore
the cited weaker inequality H³≥5^(n−1) is valid. Since13³<5⁵, a strong H11 or
H13 component has at most five vertices. The SCC product forces exactly one
such component and singleton remaining components. All possible neutral
padding positions are included because H and P₂ multiply over ordinal sums.
The primary source was independently read at
[Busch, Theorem1, pp.3–4](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/).

Independent labelled enumeration of all 1,024 five-vertex tournaments found
120 with H11, all P₂=100, and 120 with H13, all P₂=104. Every matching H was
also independently counted by all 5! permutations. This agrees with the
author's complete small-head certificate. The H13 target obstruction follows:
both45P₂±13 retain a factor13 and exceed the only seven-prime possibility91.
The H11 direct prime target9011 is correctly distinguished from a fixed-digit
gate from source2989; its canonical digit parent is3003.

Theorem D is also sound. Ordinal padding multiplies P_j by j^t, and conversion
between nonempty and possibly empty ordered bins has integer coefficients.
Subtracting the digit only changes the base-one coefficient, which is an
integer multiple of H(M0) minus d. If that coefficient is nonzero, the real
exponential sum has at most J−1 distinct real zeros. Divide by its lowest
nonzero mode and apply Rolle inductively; this proves the stated finite-hit
bound, including sparse mode sets. An infinite set of nonnegative integer
solutions forces all coefficients to vanish and therefore H(M0)|d.

This controls fixed finite frames and one fixed nonunit base module under
ordinal padding. It does not exclude a one-time Bezout match, changing modules,
multiple independently variable slots, or unbounded new frames. The note
states these boundaries correctly. No repair to either theorem is requested.
