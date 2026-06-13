# HYP-2222: LRC n=14 Uses a Pillai-Fixed Gcd-Shell Clock

**Status:** OPEN LRC proof-use hypothesis with a proved clock-classification
lemma and exact finite scout.

## Claim

The aliquot fixed-point lesson transfers to LRC through the shell modulus

```text
C = 2n - 1.
```

For odd `C`, define the antipodal gcd-shell mass

```text
A(C) = sum_{1 <= a <= (C-1)/2} gcd(a, C).
```

Equivalently, if `P(C)=sum_{r=1}^C gcd(r,C)` is Pillai's arithmetical
function, then

```text
A(C) = (P(C)-C)/2.
```

The LRC clock is **Pillai-fixed** when

```text
A(C) = C,
```

or equivalently `P(C)=3C`.  This is the gcd-shell analogue of a perfect-number
fixed equation such as `sigma(N)=2N`.

## Proved Lemma

The only odd Pillai-fixed clocks are

```text
C = 15 and C = 27.
```

So the only LRC total denominators with `2n-1` Pillai-fixed are

```text
n = 8 and n = 14.
```

Proof sketch: `P(C)/C` is multiplicative, with local factor

```text
1 + a*(1 - 1/p)
```

at `p^a || C`.  The fixed condition is `P(C)/C=3`.

If `3` does not divide `C`, then one odd prime power cannot solve the equation,
and the two smallest possible odd-prime factors already give

```text
(9/5)*(13/7) > 3.
```

If `3^a || C`, then:

- `a=1` forces the remaining factor to be `9/5`, hence a single `5`;
- `a=2` leaves a ratio below every possible odd-prime factor;
- `a=3` gives `C=27`;
- `a>3` overshoots.

Thus the odd solutions are exactly `15=3*5` and `27=3^3`.

## LRC Reading

For `n=14`, `C=27` has gcd-shell strata

```text
gcd 1: 9 shells
gcd 3: 3 shells
gcd 9: 1 shell
```

and the weighted carrier mass is

```text
9*1 + 3*3 + 1*9 = 27 = C.
```

The AP row is shell-fixed.  The known `Vstar` floor row

```text
{1,2,3,4,5,6,7,8,9,10,11,13,24}
```

is not shell-fixed: it moves shell `12` to shell `3`.  But it is fixed by the
coarser gcd/divisor carrier, because both shells lie in the gcd-`3` stratum.
So `Vstar` is a genuine divisor-carrier fixed row, not a raw shell
transversal.

## Evidence

S646 computes the Pillai-fixed clocks through `C<=999` and finds only
`[15,27]`, matching the proved lemma.

It also runs an AP single-swap scan through `n=14` using the THM-369 pair-sum
pinch oracle.  Every tight non-AP single-swap in the scan preserves weighted
gcd-shell mass:

```text
n=5:   2 -> 7    shell 2 -> 2, gcd 1 -> 1
n=6:   2 -> 9    shell 2 -> 2, gcd 1 -> 1
n=8:   6 -> 12   shell 6 -> 3, gcd 3 -> 3
n=14: 12 -> 24   shell 12 -> 3, gcd 3 -> 3
```

There are no below-floor rows in the scan.

## Proposed LRC Move

Use the gcd-shell carrier as a first obstruction before expensive lift/carry
work:

```text
mass-changing row -> strictly loose
mass-fixed row -> must pass pair-pinch + carry/owner no-leak checks
```

For `n=14`, this would focus the residual proof around the AP/Vstar
divisor-carrier fixed pocket.  The next lemma to prove is:

> Once the `C=27` gcd-shell carrier is fixed, the pair-pinch and carry-owner
> labels force either AP, Vstar, or a strict witness above `1/14`.

That is not yet proved.  The contribution here is a new fixed-clock invariant
and a finite necessary condition that turns the perfect-number fixed-point
language into an LRC proof target.

## Tournament Analysis

S646 uses proof-carrier lenses as vertices, not runners.  Considered alternate
vertex sets include runners, residues, antipodal shells, gcd strata, pair
pinches, boundary events, carry fibers, and proof obligations.  The chosen
quotient preserves tight/floor evidence and destroys raw speed order.

The majority tournament is transitive with one Hamiltonian path.  Its ranking is

```text
pillai_fixed_clock
> gcd_shell_carrier
> single_swap_zero_defect
> mod3_doubling_law
> pair_pinch_oracle
> carry_lift_conservativity
> scalar_gap_M
> raw_tuple_scan
```

Fingerprints: `score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}`,
`directed_3_cycles=0`, singleton SCCs, and `hamiltonian_paths=1`.

**See also:** HYP-2221, HYP-2220, HYP-2217, HYP-2177, HYP-2167,
THM-401, THM-369, `04-computation/lrc_pillai_fixed_clock_s646.py`,
`05-knowledge/results/lrc_pillai_fixed_clock_s646.out`,
`07-reflections/lrc-pillai-fixed-clock-carrier-s646.md`.
