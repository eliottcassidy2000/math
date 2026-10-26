# Independent audit: arbitrary-modulus Collatz clocks

2026-10-26. **PASS; no mathematical correction requested.** This review
covers C1–C5, the positive-source and refinement boundaries, and the explicit
Berggren vector/determinant example in [the producer note](../arithmetic/README.md).
The reviewed candidate SHA256 is
`682f926e35f7c6252e9af2fc5c73d07dcb7bee166876e71e40cd45ff967a8ad4`.
It does not certify Collatz convergence, a halting reduction, twin infinitude,
or the spectral portions of THM-4808. Reviewed files and any status/link-only
promotion are bound in [MANIFEST.json](MANIFEST.json).

## Proof review

**C1, including even moduli, passes.** For a proposed residue arrow, write
`n=r+Mk`. The congruence makes `D=(3^e r+e-2s)/M` integral. The desired
output congruence is equivalent to `D+3^e k` being even, which selects one
parity of `k`. When `M` is odd this also forces `n` to have parity `e`;
when `M` is even the original congruence already forces `r` to have parity
`e`. Thus no parity condition is missing in either case. Translating the
source by `2M` gives arbitrarily large positive witnesses.

There are two labelled outgoing arrows at every residue, including `M=1`.
For odd `M` the two labels can share an endpoint, but only at the stated
residue. For even `M` the label is fixed by the residue and the two endpoints
differ by `M/2`. Therefore the obstruction to a same-level deterministic
quotient holds for every `M>1`; the singleton quotient is correctly exempted.

**C2's full finite-path iff passes.** After a fixed prefix of length `j`
with `o` odd steps, changing a source within its class by `M*2^j*k` changes
the actual current value by `3^o*M*k`. The one-step argument therefore
selects exactly one parity of `k` for the next prescribed labelled arrow,
because `3^o` is odd. This proves the unique source class modulo
`M*2^L`, with no coprimality assumption on `M` and `3`. Every source in the
class realizes the whole path, and every realization lies in the class.
Consequently all finite relation paths lift to positive integers, not just
their individual arrows. The representative is allowed to depend on the
whole finite path; no infinite fixed-source conclusion is inferred.

The affine carry recursion and concatenation formula C3 have the correct
order: the suffix multiplies the prefix carry by its power of three, while
the suffix carry is multiplied by the prefix power of two. Their evaluation
requires the parity cylinder. The duration `L+o` is the literal ordinary
Collatz duration, since each odd shortcut consists of one odd step and the
following even step. It is not an unspecified machine runtime.

**C4's iff and its stated notion of covering pass.** If `N/M` is odd,
`gcd(N,M*2^L)=M`. CRT therefore gives one source class modulo `N*2^L`
for every prescribed fine starting residue and coarse path. C2 identifies
this class with exactly one fine path. If `N/M` is even, the two fine
outgoing arrows project to the same coarse arrow: one coarse arrow has
two lifts and the other has none. This holds whether the coarse modulus
is odd or even. Labels are retained, so coincident coarse endpoints do not
invalidate the counterexample. The theorem concerns outgoing path lifting;
it does not assert the incoming-degree property of a two-sided graph cover.

The shrinking-precision deterministic arrow is also correct: adding
`2^L*M` preserves input parity and changes the output by a multiple of
`2^(L-1)*M`. There is no illicit division by two at a fixed even modulus.

**C5 and the sieve composition pass with the printed bounds.** Counting
the integers `k` with `0<a+kQ<=H` gives exactly the difference of floors.
The `a=0` case correctly counts the first positive source as `Q`, and
`Q=H` still permits at most one source. The inherited height gate
`H+2<q^2`, with the actual-prime equality exceptions retained, makes the
two-entry sieve complete on `[3,H]`. Filtering the exact source class by
this integer predicate is legitimate. It is not a predicate on residues
alone, as the `3` versus `9` example demonstrates.

**The infinite-path hostile passes.** The all-zero path has source classes
`0 mod M*2^L`; the all-one path has classes `-1 mod M*2^L`. Each finite
prefix has the stated least positive realization. An integer realizing
the entire all-zero word must be zero; one realizing the all-one word
must be minus one, since `3^L(n+1)/2^L` stays integral. These are actual
nonpositive fixed points. This is the precise failure of exchanging the
finite-prefix existence quantifier with one positive source for all time.

**The Berggren example is independently proved and passes.** The printed
root-first word successively sends `(2,1)` through `(5,2)`, `(4M+1,2)`,
`(8M,4M+1)` and `(16M+2,4M+1)`. There is exactly one determinant-minus
branch. The final vector is primitive, of opposite coordinate parity and
in the positive Euclid chamber for every `M>=1`, while reducing to the root
vector modulo `M`. This proves loss of the determinant bit from the vector
observer. It does not lose that bit from the full transport matrix at odd
`M>1`, where `1` and `-1` are distinct. The note explicitly keeps these
observers separate and does not promote the pre-audit THM-4808 wholesale.

## Classical ancestry and source boundary

I read Lagarias's author-linked CECM pages on
[finite parity vectors](https://www.cecm.sfu.ca/organics/papers/lagarias/paper/html/node4.html),
[the 2-adic extension](https://www.cecm.sfu.ca/organics/papers/lagarias/paper/html/node10.html),
and [generalized decidability results](https://www.cecm.sfu.ca/organics/papers/lagarias/paper/html/node12.html).
The visible prose supports the producer's classical attribution and its
distinction between 2-adic, positive-integer and generalized-system claims.
Some mathematical displays are absent from the extracted HTML; none is
needed to validate the self-contained formulas audited above. The note
makes no classical-Collatz Turing-completeness claim.

The referenced growing-sieve definition and W2–W3 were checked directly.
They use `T(x)=(x(x+2),2(x+1),(x+1)^2+1)`, so the signed value `x=-1`
indeed yields `(-1,0,1)`. This is an inherited boundary example under a
different map, not a dynamical conjugacy or a common arithmetic theorem.

## Independent computation

[check.py](check.py) imports no repository or producer code. Its finite
universe is every modulus `1..12`, every length `0..4`, and all **2,418**
literal source classes. It enumerates congruence-defined paths separately
from integer iteration, checks positive and negative translates of each
class, compares interval floors against literal positive trajectories, and
replays native Collatz steps. It additionally tests every divisor pair
`M|N<=24` and all starting residues, including paths through length three:
there are **49 covering** and **35 noncovering** pairs.

Controls include the height-five obstruction, all-zero/all-one least
positive sources, labelled parallel arrows modulo seven, shrinking dyadic
precision, ordered affine composition, all Berggren family parameters
`M<=24`, and three height-complete sieve/path products. The checker uses
explicit runtime checks that remain active under optimization.

Normal and optimized runs are byte-identical: **89,397 checks, PASS**.
The producer's separate replay reproduces its **489,131**-check receipt.

```bash
python3 -I -B check.py > normal.json
python3 -I -B -O check.py > optimized.json
cmp normal.json optimized.json
```

The finite runs corroborate the formulas and their failure boundaries;
the all-parameter verdict rests on the mathematical arguments above.
