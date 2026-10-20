# Independent audit of the one-swap Diophantine closure

20 October 2026. **Verdict: PASS. No mathematical gap found in the stated
primitive one-swap exclusion.** This verdict combines the previously audited
critical content gate with the source application, analytic reduction,
continued-fraction completeness, and finite residual checks examined here.
It is not an assertion about unrestricted Collatz cycles or convergence.

## Audited version

The actual replacement packet is `session13-collatz-diophantine/`, specifically
`certify_one_swap_closure.py`, `closure.json`, `closure.optimized.json`,
`ONE_SWAP_EXCLUSION.md`, and `PRIMARY_SOURCES.md`.

The older `certify_log_cutoff.py` and `cutoff.json` still described rate
`1/2500`, cutoff `10^14`, and threshold `30000` when first inspected. They are
not the audited replacement and are not dependencies of this verdict.

Checked verifier SHA-256:

    2b33270c2db9053feac46bc9a88668558ae814f997caeac84267233bb37531fb

Checked normal and optimized closure receipt SHA-256:

    2d0c6c4e30d07f5d9727d599c40e7b65c0e4fcc8d8557a152a92fc1c20f4faea

Both receipts are byte-identical. Re-executing the verifier independently
reproduced the same receipt.

## 1. Imported bound and its application

I independently opened the primary paper and checked its hypotheses,
definition of the multiplicative form, and real-field version of Theorem 9.4,
printed page 16:
[Bugeaud, Mignotte, Siksek](https://arxiv.org/pdf/math/0403046).

The application is to the nonzero multiplicative form
`Lambda=2^q*3^(-p)-1`, not directly to `q log 2-p log 3`. Take the real field
`Q`, degree `1`, two algebraic numbers `2,3`, integer exponents `q,-p`,
and `B=q`. The required heights are bounded by `A1=7/10,A2=11/10`.
All hypotheses hold. The source's real-field constant specializes to

    C=1.4*30^5*2^(9/2)*(7/10)*(11/10).

Since `2^(9/2)<23`, certified by `512<529`, one has
`C<602494200<10^9`. The bound therefore gives
`log Lambda>-10^9(1+log q)`. No recalled irrationality-measure constant,
general Collatz verification threshold, or minimum-cycle theorem is used.

The already audited norm bound gives the necessary upper estimate directly:
if the word represented a positive integer cycle, `D=g` and

    0<Lambda=D/3^p<(7/9)^(p/2)<exp(-p/9).

Thus there is no missing factor caused by switching between the
multiplicative form and its logarithm.

## 2. Global bound and Legendre threshold

Combining the two estimates and `q<2p` yields

    p<9*10^9*(1+log(2p)).

The function `p/(1+log(2p))` is increasing for `p>=1`: its derivative has
the sign of `log(2p)`. At `p=10^12`, the certificate proves
`log(2p)<30`, and `9*10^9*31<10^12`. The monotonicity therefore excludes
every `p>=10^12`, rather than just the tested boundary.

For the additive gap `lambda=log(1+Lambda)`, one has
`0<lambda<exp(-p/9)`. The inequality

    2p*exp(-p/9)<log 2

is verified at `p=100`. Its left side is decreasing for every `p>9`, so it
holds throughout the required range. Consequently

    0<q/p-log(3)/log(2)<1/(2p^2)

for `p>=100`. Coprimality makes `q/p` reduced. Legendre's criterion forces
an ordinary continued-fraction convergent; there is no omitted
semiconvergent case.

## 3. Interval validity and complete convergent enumeration

The production verifier's positive atanh series has the claimed tail bound:
in the omitted terms, replacing each denominator by its first omitted
denominator gives the stated geometric upper bound. Thus the log intervals
are rigorous. Dividing positive interval endpoints in opposite order gives a
rigorous interval for `alpha=log(3)/log(2)`.

Each continued-fraction digit is accepted only when both endpoints have the
same integer part. Subtracting that digit and reciprocating in reverse order
then encloses the next complete quotient. The intervals remain positive.
The recurrence starts with the usual convergent initial values. It includes
every convergent, not a selected numerical list.

There are exactly 18 convergents with `100<=p<10^12`. Nine have a rigorously
negative gap, and nine have positive gap lower bounds exceeding
`2^(-min(floor(p/9),256))`. This last quantity is at least `exp(-p/9)`,
because `min(floor(p/9),256)<=p/9` and `log 2<1`.

The terminal convergent is

    8573543875303 / 5409303924479.

Its denominator exceeds the global cutoff. All subsequent denominators
increase, so the required range is exhausted. Gap comparisons use exact
rationals, not the displayed decimal approximations of older drafts.

## 4. Independent implementation and residual words

The separate checker `verify_diophantine_independent.py` imports no code from
the production verifier. It uses a different logarithm enclosure:

    -log(1-t)=sum_(k>=1) t^k/k,

with 240 terms and geometric remainder. It obtains `log 2` at `t=1/2`,
and `log 3` by adding the bounds for `log 2` and `log(3/2)`, the latter
at `t=1/3`. It independently reproduces the complete continued-fraction
prefix, all 18 candidates and their exclusions, and the terminal convergent.

For `2<=p<100`, it computes the critical length by the integer bit-length of
`3^p`, rather than by a logarithm interval. It constructs every actual cyclic
swap and calculates carry right-to-left, independently of the production
verifier's left-to-right affine composition. It reproduces:

- all 98 small weight entries and their coprimality decisions;
- exactly 59 admissible critical shapes;
- exactly 3548 unequal cyclic-adjacent swap occurrences;
- each shape's exact maximum content and strict inequality `g<D`.

The production certificate also compares its two carry implementations.
Cut-crossing swaps are present because each loop includes the edge joining
the last and first positions. Duplicate descriptions are retained, not lost
through a potentially incomplete necklace deduplication.

The theorem's `p=1` case is independently elementary: the carry is a power
of `2`, while `D` is odd, giving `g=1<D`. For `p>=2`, the content gate first
forces the critical clock; the large and small ranges then exhaust all
possibilities. The `p=2` critical shape is noncoprime. The `D>1` condition
and the ordinary trivial cycle are handled consistently.

## 5. Reproduction and conclusion

From the workspace root:

    python3 session13-collatz-resultants/refinement/audit/verify_diophantine_independent.py
    python3 -O session13-collatz-resultants/refinement/audit/verify_diophantine_independent.py

After relocating the files, pass `--closure-report /path/to/closure.json` to
select the production receipt explicitly.

The outputs saved as `diophantine.independent.json` and
`diophantine.independent.optimized.json` are byte-identical and report PASS.
They are corroborating certificates, while the all-weight steps are the
symbolic norm gate, the cited transcendence bound, monotonicity, and Legendre's
criterion.

The claimed exclusion follows for the specified coprime, single-swap
mechanical family. No unresolved mathematical condition remains in this
closure argument beyond its explicitly cited and independently source-checked
classical theorem. The frozen earlier analytic audit files were not modified
while performing this closure audit.
