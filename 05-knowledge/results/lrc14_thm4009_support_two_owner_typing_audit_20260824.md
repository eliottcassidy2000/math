# THM-4009 support-two owner typing in the THM-3910 branch

**Status: PROVED arithmetic classification + FINITE-EXACT hostile audit.**
This is an information-loss repair, not an LRC(14) closure.  It refines the
suggestion in the THM-4009 Fourier/LP audit to compare the `47` short
support-two ratios with THM-3910's `17` surviving external pair types.

## 1. Three owner classes

For two positive integer speeds `x,y`, the primitive relation supported on
those coordinates is, up to global sign,

```text
(y/g)e_x-(x/g)e_y,                 g=gcd(x,y).          (1)
```

Thus THM-4009's bound `||a||_2^2<=195` applies to the **reduced speed
ratio at the labelled support**, not to an unlabelled pair of numbers.  In
THM-3910's row

```text
(u_1,...,u_11,tp,tq),
```

a support-two relation has exactly one of three owners:

```text
body/body,              body/extra,              extra/extra. (2)
```

For extra/extra ownership the common scale `t` cancels.  Of the `17`
THM-3910 types, `15` have reduced pair norm square at most `195`; the only
failures are the scale-one types

```text
(8,21): 8^2+21^2=505,          (9,11): 9^2+11^2=202.   (3)
```

This excludes only **extra/extra ownership** in those two cells.

## 2. Uniform body/extra cutoff

Let `U=max u_i` and pair a body speed `u<=U` with the extra speed `tp`.
The larger primitive coefficient is

```text
tp/gcd(tp,u).
```

Every integer coefficient in a vector of squared norm at most `195` has
absolute value at most `13`.  Therefore body/extra ownership requires

```text
tp <= 13 gcd(tp,u) <= 13u <= 13U.                      (4)
```

Consequently `tp>13U` rules out every body/`tp` support-two relation without
inspecting the body.  The analogous statement holds for `tq`.

Combining `(3)--(4)` gives a correctly typed conditional survivor statement:
if the body has no short body/body relation, the external pair has reduced
norm square above `195`, and both `tp,tq>13U`, then THM-4009's forced short
relation has support at least three.  None of these side conditions holds
automatically in THM-3910.

## 3. Canonical hostile to an owner-free intersection

For the AP11 body `(1,...,11)`, the fixed relation

```text
2*(speed 1)-1*(speed 2)=0,             norm square=5  (5)
```

survives for every external pair and every scale.  Hence all `17` types have
a THM-4009-short relation even when their extra/extra ratio is incompatible.
The exact checker replays this hostile on all `17*133=2261` rows with
`11<=t<=143`, the complete AP11 range in which a body/extra relation could
still meet `(4)` for `p=1`.

Therefore the proposed `47 x 17` calculation is valid only as an
**owner-labelled discovery product**.  It is not a mathematical intersection
of two type lists.  A decisive atlas must retain the two support labels, their
body/extra roles, component incidence, and the endpoint-owner word already
identified by the Fourier audit.

## 4. Exact audit

The companion performs the `47`-ratio and `15/2` type censuses, the AP11
hostile replay, and `68,007,600` body/extra cutoff checks.  Normal and
optimized transcripts agree and match the frozen output.

```text
python -B 04-computation/lrc14_thm4009_support_two_owner_typing_20260824.py
python -B -O 04-computation/lrc14_thm4009_support_two_owner_typing_20260824.py
```

Artifacts:

```text
script_sha256 = d00b3b6206818e8c8a7d8bbba5b153473783a98e341e611dd5d3e90e381e060f
output_sha256 = e1b3cf909d1f10bc7adca26ae5044dc492a29b33a01551badf9aa390f361ffa1
```

**OPEN:** use the owner-labelled decomposition to cross THM-4009 with the
actual THM-3910 body fibres.  The short relation alone still loses common
loneliness time, phase, arrival, and weak endpoint ownership.  LRC(14)
remains open.
