# Independent audit of Session236 offset transport

Verdict: **PASS WITH PROVENANCE CORRECTIONS**. No false identity, sign error,
parity error, or overbroad mathematical quantifier was found in the submitted
`offsets/PROOF.md`, `offsets/replay.py`, or `offsets/results.json`.
The claimed newness of two mechanisms needs narrowing as described below.
Chan–Pak completeness remains OPEN. This audit did not inspect the other
Session236 lanes or establish any assertion about unrestricted port orbits.

## Required provenance corrections

1. **The singleton-right-block offset transport is inherited near-join algebra.**
   In `01-canon/theorems/THM-4865-the-sole-exit-tower-is-tribonacci-h-k-equals-h-k-1-plus-h-k-2-plus-h-k-3-for-every-rooted-base-the-lineage-c2-law-at-any-vertex-and-the-three-root-response-numerations-attain-every-odd-number-in-69-702189.md`,
   part (d) explicitly credits HYP-9453's law `h+2^k B(T,t)`.
   HYP-9453's full filename is
   `05-knowledge/hypotheses/HYP-9453-the-near-join-sweep-toward-chan-pak-completeness-density-of-b-sets-two-strong-realisers-for-every-attained-odd-count-at-least-15.md`.
   Although the coverage conjecture there is OPEN, the near-join identity is
   separately explicitly PROVED. Its B is exactly this packet's R=E+O: partition
   into a rooted starting path and a Hamilton path of its complement. Its
   tournament `T -> M_k -> z`, reversing the edge between z and the old root,
   is exactly the B=singleton instance here. Cite this before the examples and
   call the general right-block formula a corollary/extension of the inherited
   near-join and seam calculus, not a newly discovered offset mechanism.

2. **Strict cyclic drift positivity and its run proof have a direct precedent.**
   `05-knowledge/results/signed_responses_slope_classification_clock_barrier_20261105.md`,
   proof of statement 3, lines 166–170 in the audited version, gives
   `Delta_cyc = sum_i epsilon_i (U_i-epsilon_(i-1)) > 0` for a mixed cyclic word,
   explicitly including runs crossing the word boundary, and credits Session198.
   It is the source note of
   `01-canon/theorems/THM-4882-signed-one-apex-responses-slope-and-parity-classification-session198-clock-merged-positive-count-barrier.md`.
   The exact dictionary is `U_i=W_i+1` and `epsilon_i=1-b_i`.
   On a zero digit b_i, its contribution is
   `U_i-epsilon_(i-1)=W_i+b_(i-1)`; on a one digit it is zero.
   Thus its Delta_cyc is exactly the present D, not merely an analogous drift.
   The old source writes `U'=2U-epsilon`; this becomes `W'=2W+b` under the same
   dictionary. Its cyclic formal orbit lies in (0,1), the present orbit in (-1,0).
   The all-length claim follows from the identical cyclic run argument; no
   period bound occurs in that argument. Replace “New all-length obstruction”
   and “strengthens a finite search ... to a proof” with an explicit inherited
   positivity/rederivation statement. The constant-offset and arithmetic-
   progression exclusions can be offered as the present corollaries and the
   coordinate dictionary as useful synthesis. I did not independently find the
   original Session198 files in this checkout, but the available November 5
   source states the relevant formula and run proof directly.

Neither correction changes any equation or computed value.

## Independent all-order verification

### Padding and bridge without a seam identity black box

For one appended sink s, every rooted path either stays in the old tournament
or is an old rooted path followed by s. These contribute equal old complement
counts, because a remaining sink is last in every Hamilton path. The second
choice reverses order parity. Therefore E'=O'=E+O, H'=H. Iterate to obtain (1).
In particular `R(T_m,u)=2^m R(T,u)`.

For m>=1 the old u cannot end a Hamilton path of T_m: there is a sink after the
old tournament and that sink has no outgoing arc. A Hamilton path of the
ordinal sum T_m=>B therefore never uses u->w, so reversing that arc removes
none of the `H(T)H(B)` old paths.

Every new Hamilton path uses w->u exactly once. It has the unique four-part
form

    [path in T_m outside P] [path Q in B ending at w]
    [path P in T_m starting at u] [path in B outside Q].

The first and last parts may be empty. Conversely each such choice is a legal
path, since all other crossing edges point from T_m to B. The number of choices
is `R(T_m,u) * (q_e+q_o)`. Hence directly

    H(J_m)=H(T)H(B)+2^m R(T,u)(q_e+q_o).

This also explains precisely why the m=0 extrapolation is invalid: u may then
end an old Hamilton path, which is destroyed by the reversal. The generic
un-padded formula instead subtracts
`(# paths of T ending at u)*(# paths of B starting at w)` from
`H(T)H(B)+R(T,u)(q_e+q_o)`.

### Offset and coefficient barriers

The constant term C is the actual ordinal-sum count H(T=>B), so its positivity
and exclusions 7,21 are exactly the inherited all-order spectrum theorem
`01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`.
They are a restriction on the family parameters, not an exclusion of targets
with a different scale or decoration. The proof correctly states this boundary.

A=1 forces both positive integer response factors to equal 1. An out-neighbor
of u supplies a second positive rooted-path summand; hence R=1 forces u a
sink, and then R=H(T-u). H=1 iff transitive gives the claimed rigidity. The
converse argument on B is identical. Thus A=1 implies C=1, with no omitted
empty-tournament exception because a port is required on each side.

### Periodic words

The bit recurrence, affine block coefficients, and exact solution (6) are
consistent. For a block of length k, the formal periodic starting point is
`-d/(2^k-1)`; before each subsequent bit it is the negative fraction for a
cyclic rotation. All such values lie strictly between -1 and 0 for mixed words.
For a maximal zero run of length l preceded by 1, let the starting point be -r.
The run contributes `1-(2^l-1)r`; because the endpoint -2^l r remains in (-1,0),
this is greater than 2^(-l). Summing over cyclic runs proves D>0. A run crossing
the chosen cut is allowed: the formal orbit and previous-bit parity are cyclic,
so reordering the finite sum changes nothing.

For the explicit fixed-clock obstruction, set h_m=H_m. Then

    h_(m+2)-(q+1)h_(m+1)+q h_m = -(q-1)D != 0.

An exponential plus constant with base q has left side zero. On the subsequence
m=a+jt, q becomes q^t and the linear coefficient becomes tD>0, so the same
argument applies. If some other fixed positive dyadic base is proposed, the
positive exponential leading term forces that base to be q (or q^t on the
subsequence) by the ratio limit. This does not exclude arbitrary nonperiodic
words, varying blocks, or a nonlinear redefinition of the clock.

The first-block parity guard is essential and correctly stated. Leading zero
bits do not change the represented H and do not invalidate either recurrence.
The all-zero and all-one cases have D=0 as claimed; an odd seed before a zero
block first needs the explicitly permitted initial-block discard.

## Reproducible independent computation

Run from the repository root:

    python3 05-knowledge/results/session236_chanpak_offsets_orbits_20261109/offsets_audit/independent_check.py

No repository code or submitted verifier is imported. Hamilton paths are counted
by a memoized first-vertex deletion recursion (not the submitted endpoint DP).
Root responses enumerate partial permutations. One-apex counts are independently
compared with increasing even-size chain subsets whose first/last bits are 1/0,
which count the odd cycles through the apex. Period drift is calculated from
cyclic integer rotations; affine block coefficients are separately recovered
from increments at two admissible starting integers, not the submitted
coefficient accumulator.

`results.json` records:

- 116 padding checks on all labelled rooted tournaments of orders 1–3, m=1–4;
- 1,682 bridge checks on all ordered pairs of those roots, m=1,2;
- 511 independently counted apex tournaments, all words of lengths 0–8;
- all 131,070 binary words of lengths 1–16: positive drift iff mixed;
- 2,046 explicit cyclic-cut invariance checks for all lengths 1–10;
- hostile singleton m=0 and equal-H/unequal-response controls.

The supplied replay was also executed and its full JSON was byte-identical to
the submitted `offsets/results.json`; the reproduced output is retained as
`supplied_replay.json`. This confirms the packet's stated 87 padding, 1,682
small bridges plus 50 random bridges, 2,046 periods, 254 short apex counts and
both concrete rays. Finite tests are checks, not the proof of all-order claims.

## Correction verification

The revised `offsets/PROOF.md` applies both required provenance corrections above
in the status, inheritance section, section headings, detailed coordinate
dictionary and closing scope paragraph. The full correction diff was inspected.
The revised proof therefore receives **PASS** within this audit's mathematical
and provenance scope. The original verdict is retained above to preserve the
correction lineage. No mathematical formula changed.
