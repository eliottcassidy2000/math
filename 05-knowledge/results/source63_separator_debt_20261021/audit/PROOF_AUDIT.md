# Independent audit: source63 separator debt

Date: 2026-10-21. Verdict: **PASS within the explicitly stated scopes.**

Audited producer note: `SOURCE63_SEPARATOR_DEBT.md`, SHA256
`8077d5e3552f93bca4dda19f5314bdadc26556396ea1df8658c4e9ca9466a22a`.
The verdict applies to that frozen version and its manifested input packet.
No mathematical correction or producer-file change was required. This is an
independent mathematical/code audit, not a formal proof assistant certificate
or a claim of historical novelty.

## 1. Independence, provenance, and computational universe

The audit checker does not import or execute the producer's verifier. It:

1. Parses the word and large exact fractions displayed in the note.
2. Derives affine carries and signed emissions from unrolled common-denominator
   sums, rather than importing the producer's state routines.
3. Reconstructs the inherited sufficient depth-18 test on exponents 1 through
   1024 and obtains exactly the same 63 residual exponents.
4. Discovers every retained source trajectory directly from the actual integer
   `125^j`, without using the supplied valuation word to drive the trajectory.
   The newly discovered complete word must equal the decompressed inherited
   word byte for byte and have its retained SHA256. The actual child, suffix,
   source/child registers, first half-descent clocks, and first ROOT clocks
   are separately checked. No terminating representative is substituted.
5. Checks exact exponent guards for four terminal valuations per template;
   performs denominator-preserving division modulo 15, 35, and 105; and replays
   the displayed complement's 149 odd shortcut branches modularly.

All 14 source-inventory hashes were also checked against actual Git object
contents at revision `c4a7875652a5ea516e93a31c3c6410d2204773b2`. The frozen
`guarded_children_105_20261018g.json` is byte-identical to the original at that
revision, with SHA256
`c08611140beb4db0545e67de0d9ff8f066d6bfbd909b29ca6e9b716432800822`.
This authenticates correspondence to the pinned repository content, not the
real-world authorship of a Git commit. No remote or repository file was changed
by this audit.

Normal and optimized independent runs each perform **1,204,846 active checks**
and produce identical receipts, SHA256
`dd99d61117672ac15548759e2c7d64f7a1db9757ac738694ec10d9c9f3c94636`.
Checks use explicit exceptions, so optimization does not disable them.
All branching arithmetic is integral or rational. No floating-point decision
or materialization of `2^H`, `N_1`, or an enormous lifted `125^J` is used.

The retained finite source totals are independently recovered as:

- 575,652 odd steps;
- 1,152,705 shortcut steps;
- 1,728,357 ordinary unaccelerated steps.

Twenty-seven of the 63 supplied seeds have different first literal and first
odd-checkpoint half-descent shortcut clocks. The audit retains both, rather
than rounding the earlier literal event to its later odd endpoint.

## 2. Signed coordinate, exact guard, and typing

For one letter `a`, the shifted coordinate `z=n-1` transforms as

    z' = (3/2^a) z + (4/2^a - 1).

Thus a valuation-one letter emits positive charge 1, valuation two emits zero,
and a letter at least three discards debt `1-4/2^a`. Propagating these emissions
through their actual ordered suffixes gives the stated nonnegative `P,D` and
the exact signed difference `P-D`. At an actual endpoint `h`, the positive
envelope lies strictly below `n-1` exactly when `D<n-h`. This is not an
absolute-value identity, a native signed-row energy identity, or a payment
of Gram cross terms.

The source guard modulo `2Q`, not just modulo `Q`, is essential. For a positive
word, the full affine congruence with odd endpoint is equivalent to the exact
valuation word. One can inductively remove its first step: divisibility gives
an integral first quotient, and the next odd terminal condition fixes its odd
bit; repeating yields every prescribed exact valuation. Conversely an exact
word has an odd endpoint, giving `A n+B=Q mod2Q`. Since `A` is odd, the displayed
inverse modulo `2Q` is lawful. Positivity follows at every native step from a
positive source.

The audit preserves the distinction between literal source 63 and the count
of 63 residual exponent templates. The mod15/mod35 observer has only 105
compatible joint states. Its state `(1,1)` is not the actual integer ROOT.
Every ROOT claim audited here ends at the integer 1, at its first arrival.

## 3. S1: separator normalization and its exact limitation

Fix `theta=(2/3)^r` and a single common `c`. A separator of `q>=r` ones has
multiplier `eta=(3/2)^q`, debt zero, and `eta theta>=1`. Consequently

    D_usv <= c(lambda_usv-theta)
             + c lambda_v(1-eta theta)
           <= c(lambda_usv-theta).

Both factor multipliers are at least one, so the multiplier condition is also
preserved. This proves unrestricted finite composition with changing
generators and separator lengths, provided the fixed common cone conditions
hold. It does not license changing `theta,c` arbitrarily between joins.

For the separator of exactly `r` ones, appending that separator to a word
multiplies both its lambda and debt by `eta=1/theta`. Its inequality becomes
`D_us<=c(lambda_us-1)`, and the map intertwines the separator product with
ordinary concatenation. This normalization is algebraic. The note correctly
does not insert that separator into source63's actual orbit.

For a nonzero-debt generator with `lambda>=1`, the least admitted constant is
`D/(lambda-theta)`. A finite choice of `r` with constant at most one, or
strictly below one, exists precisely when `D<lambda`: theta tends to zero but
never reaches it. At `D=lambda`, every finite constant is greater than one.
Multiplication by any growth-only separator preserves the sign of
`lambda-D`. The hostile `(3,1,1,1,1)` and the failure of ordinary `u u` in the
specific source63 cone are valid boundary controls, not general no-go results
for every possible controller.

## 4. S2 and the all-depth source quantifier

For a legal cone prefix with lambda at least one, a next actual descending odd
step has `mu lambda<h/n<1`, because the full unshifted affine carry is positive.
Thus its valuation is at least two. For `c>0`, the strict multiplier inequality
gives the note's strict debt estimate. Valuation two yields debt less than `c`;
a valuation at least three yields debt less than `1+c`, hence at most the
strict bound 2 when `c<=1`. Since positive odd endpoints differ by at least
two, payment follows, with strict margin greater than `1-c` if `c<1`.
For `c=0`, the debt is just the final one-step debt, less than one. The separate
zero case correctly avoids multiplying a strict inequality by zero.

The common source63 constant is exactly

    14060731968178560/14425381885981321 < 1.

The ordinary cone constant is greater than one but no greater than the
inherited terminal-packet allowance `20/13`. Therefore this result is not
misrepresented as the first available proof paying source63's own exit.

Every prefix of each chosen generator expands, as do one-letter separators.
For any specified finite composed prefix, choosing a sufficiently large final
valuation makes the full affine multiplier less than one. Its exact cylinder
has infinitely many positive representatives. All sufficiently large such
representatives satisfy the descending affine inequality, and all earlier odd
checkpoints remain above the source. This proves

    for every finite chosen macro word, infinitely many legal paid sources,

not a single finite source executing all depths, and not a policy by which
every original source can choose its next valuation. The equality controls up
to 64 copies verify the displayed 2,176 odd steps, 257 debt resets, 3,449
shortcut steps, and 1,038-digit source. These controls support the formulas;
the all-depth result rests on the composition proof.

## 5. Source63, fixed-word lifts, and the new N_t family

Independent native discovery from 63 yields precisely the displayed 33-letter
prefix followed by valuation 4 and child 61. Its first odd and literal
descent clocks are respectively `(34,54,88)` and `(54,88)` in the stated clock
types. Its exact first-exit debt and positive margin are

    D = 215402584767891/140737488355328,
    2-D = 66072391942765/140737488355328.

Every proper odd-prefix multiplier and the coefficient after three divisions
of the final run exceeds one. The full word contracts and already descends
at the least positive cylinder source 63. This proves the claimed first
arrival clocks throughout `63+2^55 t`, `t>=0`. The affine and positive-envelope
barriers match the displayed fractions, with their different units retained.

The actual suffix from 61 is `(3,1,1,5,4)`. The original 63 first reaches ROOT
in 39 odd, 68 shortcut, and 107 ordinary steps. A positive lift in the same
39-letter cylinder reaches `1+2*3^39 t`, not ROOT. The audit explicitly checks
such a lift, preventing propagation of the named ROOT certificate to the
whole first-exit cylinder.

For `H=2*3^33`, elementary lifting gives the exact order of 2 modulo `3^34`.
The base identity `61*2^54-B=63*3^34` then proves every displayed `N_t` is an
integer; its odd numerator divided by an odd denominator proves oddness.
It increases from 63, so it stays positive. Its affine endpoint is the odd
integer 61, giving the complete exact terminal guard, not just divisibility.
The exact word is therefore `u,(4+Ht)` and the same actual suffix to ROOT is
available at the same actual child.

After four divisions in the final run, the value is `61*2^(Ht)`. The displayed
strict gap has base value `2*3^34>0` and increases thereafter. Three divisions
still have expanding coefficient. Hence the first literal clock remains
54/88, whereas the first odd shortcut clock is `54+Ht` and the ROOT clocks
are `68+Ht` shortcut and `107+Ht` ordinary. This proves the separated-clock
claim symbolically. Modular controls do not pretend to materialize `N_1`.
For `t>0`, these are new sources, not any replacements for supplied 125-powers.

## 6. S4: all63 guards, sharp common envelope, and ROOT boundary

For every template the audit independently checks the anchor divisibility,
the minimum-final-valuation half inequality, all prefix coefficient tests,
and `0<j0<2^K`, with `K=A0+c0-2`. The exponent isometry

    v2(125^j-125^j0)=2+v2(j-j0)

preserves the native prefix and permits actual terminal valuation `c>=c0` on
the original exponent class. Every positive exponent in that class is at
least `j0`. The affine half-descent inequality therefore persists as the
source increases. No affine endpoint or terminating-source representative
replaces the actual `125^j`.

For `c>=3`, the debt is `1+(3D_v-4)/2^c`; at `c=2` it is `3D_v/4`. The c=2
boundary is correctly separate. These formulas give
`sup D(c)<=max(1,D(c0))`. The positive envelope barrier is
`1+3P_v/(2^c-3lambda_v)`, with positive denominator at `c0` and thereafter,
and is nonincreasing in c. Thus exact finite coefficient bounds at the 63
minimal valuations prove the all-terminal-valuation statements. There is
no hidden valuation cutoff.

The displayed common barrier and maximal minimum-valuation debt are exactly
reproduced. Both maxima occur at template 638, whose actual anchor valuation
is 8. The common barrier is therefore the least common barrier for these
real affine envelopes, not merely an upper estimate. It need not be the
least threshold for a guard-aware discrete argument, a ROOT height bound, or
a signed-energy bound. All actual family sources exceed 9 and all these
debts are less than 2, so the inherited half descent pays the scalar charge.

For the native odd-modulus observer, `ord_105(2)=12`, so retaining c modulo12
is sufficient. The full c is still needed for exact division and clocks.
All prime-register computations retain the denominator before division and
the mod5 compatibility condition. All 63 supplied original child suffixes
are authenticated and rediscovered at their original sources. None of this
grounds the infinitely many other children of the exponent classes.

## 7. Infinite complement and stopping boundary

The largest sufficient literal half clock is 149. Even when the actual final
valuation is larger, the c0-th division has already halved the source, so this
bound is on the literal clock, not the possibly later odd clock.

The displayed residue modulo `2^149` satisfies
`3*125^j+5=0 mod2^151`. It forces 149 successive odd shortcut branches after
the common two initial steps. The explicit affine formula is above half of
the source at its starting point and increases on those odd branches.
Neither initial step halves the source. Hence no source in this positive
exponent progression halves by shortcut151. The audit parses this displayed
residue, checks its congruence, replays the branches, and tests disjointness
from every new guard and the baseline atlas.

The progression is infinite and unhandled by these particular finite atlases.
It is not a divergent trajectory or an exponent satisfying the all-ones
constraint at every depth. The arbitrary-L argument obstructs a fixed
uniform half-descent horizon in the 125-power family, not source-dependent
horizons or future stronger coverage mechanisms.

Accordingly the audited results leave universal source entry, all lifted-child
grounding, arbitrary exponent acceptance, and native signed ROOT-energy
cancellation **OPEN**, exactly as the producer says. The scalar debt proof
does not discharge any of those separate obligations.
