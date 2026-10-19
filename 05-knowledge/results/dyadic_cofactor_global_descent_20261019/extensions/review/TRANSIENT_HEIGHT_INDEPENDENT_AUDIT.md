# Independent audit of transient-height scaling

Date: 19 October 2026.

## Verdict and scope

**Pass.** The stated universal excursion envelope and the limit superior

    limsup (M(n)-3n)/n^(2/3) = 3

over positive even starts tending to infinity follow from the previously
audited cofactor-map descent theorem. No repair is required.

The separately proposed sharpening also passes:

    M(d²(d-1)) = 3d³-2

for every power of two `d>=8`.

The proof concerns the auxiliary cofactor map, with no unrestricted aliquot
or Collatz height conclusion. Native prime-domain validity is not asserted
along the resonant family.

Audited source:

`TRANSIENT_HEIGHT_SCALING.md` before its documented status-only publication update

SHA-256:

`a1d3cfb1c9093da6a2c2519201ef72e0856ca16996f303c705e4df119acb77da`

The final reviewed source includes the exact-maximum sharpening. The earlier draft hash was
`d76080861ff03be0f330607a9713daac4e529b035f540221dd4d66e5ff5f368b`.

## Local excursion review

An exact useful identity is

    3n-2-F(n) = (p-1)(2c+1) >= 0.

Thus the first image always has the claimed basic ceiling. Immediate descent
and the non-power-of-two `q` case are settled by this identity and their
already-verified smaller endpoint.

For a power start `n>=32`, the three consecutive images are
`3n-2`, `3n/2+2`, and `3n/4+4`; each is at most `3n-2`.

When `q=d>=4` and `c>d`, write `e=c/d>=2`. The first image satisfies

    n_1 = n+2c-d <= (1+2/(d-1))n <=5n/3.

The second-image difference is correctly given by

    n_2-n = e(d-1)+d.

It is at most `n=e d(d-1)` because the difference between the right and left
sides is `e(d-1)²-d>0` for `e>=2`, `d>=4`. Hence `n_2<=2n`.
Both ceilings are strictly below `3n-2` in these cases, since `n>=24`.
The remaining middle-branch image descends. In the high branch, either the
third image descends or it is `n+d<2n`, followed by descent. Consequently no
nonresonant intermediate value was omitted.

For the first resonance, the block reaches `d²`, then `3d²-2`, then the two
smaller power-start values. Its maximum is `3d²-2`. For the second resonance,
the initial values are

    d³-d², d³+d²-d, d³, 3d³-2,
    3d³/2+2, 3d³/4+4.

For `d>=8`, the maximum is `3d³-2`. The preceding first image is below that
ceiling as well. For `d=4`, the actual exceptional/six-step chains supply
the stated maxima 46 and 190. The other finite exceptions have the claimed
maxima 22 or 46.

## Propagation review

The verified descent algorithm partitions the actual orbit into finitely many
blocks whose starting checkpoints strictly decrease, ending at a fixed point
or finite exception. Therefore every block starts at `m<=n`, where `n` is
the original start. Its intermediate points need not decrease, but their
local ceilings have all just been bounded.

At a first resonance `m=d(d-1)`, the inequality `d<=sqrt(n)+1` follows from
`(d-1)²<=m<=n`. Hence

    3d²-2 =3m+3d-2 <=3n+3sqrt(n)+1.

At a second resonance `m=d²(d-1)`, the inequality `(d-1)³<=m<=n` gives
`d<=n^(1/3)+1`. Therefore

    3d³-2 <=3n+3n^(2/3)+6n^(1/3)+1.

The latter expression dominates the first-resonance expression for `n>=1`,
already because `n^(2/3)>=sqrt(n)`. It also dominates all nonresonant
ceilings. Including the terminal exceptional tails gives exactly

    M(n)<=max(46,3n+3n^(2/3)+6n^(1/3)+1).

After subtracting `3n` and dividing by `n^(2/3)`, the right-side asymptotic
upper bound tends to 3. The fixed finite term 46 causes no asymptotic issue.

## Lower family and exact maximum

For dyadic `d>=8`, the start `n_d=d²(d-1)` has exact cofactor `d²` and odd
part `d-1`. Its next input has exact cofactor `d`, and the next is the power
`d³`. Thus the stated consecutive prefix is genuine and reaches `3d³-2`.
The resulting ratio is bounded below by

    (3-2/d²)/(1-1/d)^(2/3),

which tends to 3. This proves the matching limit superior without requiring
any assertion about later maxima.

The proposed exact-maximum strengthening is also valid. After five steps,
let `m=3d³/4+4`. For `d>=8`, `m<d³`. The nonconstant universal-envelope term
at `m` is larger than 46, so

    M(m) <=3m+3m^(2/3)+6m^(1/3)+1
          <(9/4)d³+3d²+6d+13
          <3d³-2.

For the final strict inequality the gap is

    (3/4)d³-3d²-6d-15
     =(3/4)d²(d-4)-6d-15
     >=3d²-6d-15>0.

It is already 129 at `d=8`. Since the verified prefix has maximum
`3d³-2`, no subsequent point can improve it. This proves the exact formula.
There is no circularity: the universal envelope was proved using block
ceilings, not exact full-orbit maxima for this family.

## Independent executable controls

Files:

- `audit_transient_height.py`
- `transient_height_audit_results.json`
- `transient_height_audit.out`
- `transient_height_audit.optimized.out`

All arithmetic and comparisons in the checker are exact. The checker uses
explicit failures rather than Python assertions, and normal/optimized runs
produce identical PASS output.

Controls passed:

- Every local excursion from the 2,000,000 even starts through 4,000,000.
- 162,000 additional structured cases with dyadic exponents through 180.
- Independent complete-orbit maximum calculations for every even start
  through 1,000,000, visiting 500,011 distinct nodes.
- Exact resonant maxima for all `d=2^t`, `3<=t<=256`, including 768-bit starts.
- 300 reproducibly seeded additional complete orbits with cofactor
  exponents and odd-part bit lengths through 1,024.

For numerical envelope checks the program compares the exact maximum with
the smaller integer quantity obtained by replacing `n^(1/3)` with its floor
in the increasing expression `3n+3x²+6x+1`. Passing this stronger finite check
therefore verifies the manuscript's envelope without floating-point error.

These controls supplement the proof; they do not supply its universal or
asymptotic quantifiers. No novelty or historical priority is established.
