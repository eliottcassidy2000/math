# Independent audit: four-window positive amplitude repair at (-27,1,15)

**Status: independent analytic/source audit PASS, with independent exact
certificates.** The producer's fixed-row representation theorem is accepted
without a mathematical repair. Its coefficients are common positive real
algebraic numbers on the four selected positive square-root branches.
Neither an all-eight-root identity nor a rational odd-atom representation
is claimed or possible in the stated nonnegative class. General Laurent
noncancellation and an all-even certificate at this row remain outside the
result.

## Inheritance and source normalization

I read the proof and certificate first, rebuilt the source and determinant
signs independently, and only then read the producer program. The mathematical
inheritance is the original carried A2B3 row and the coupled-window theorem
in `05-knowledge/results/creative_20260906_laurent_bridge.md`, relative to
**THM-4440**,
`01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md`.
The earlier h=3 integer-phase separator and h=4 single-window obstruction
are scope controls, not evidence that the mixed-window representation fails.
The selected positive square roots and the allowed coefficient field must
remain explicit.

Literal charge enumeration for all masses from 1 through 28 yields returns
only at 14 and 28. This also has an elementary proof: all three charges
are one modulo 14 and `(1,12,1)` is a zero-charge count triple. The complete
first five and doubled ten fibers are exactly

```
(1+j,12-3j,1+2j), j=0,...,4;
(2+j,24-3j,2+2j), j=-1,...,8.
```

Their multinomial coefficients reproduce the original first polynomial and
the complete carried target. In particular the phase is `gamma^2`, not a
new relative parameter on one grouped response. At `gamma=i sqrt(s)`, the
actual constant terms are `gamma P_raw(-s)` and `gamma^2 Q_raw(-s)`.
The normalizations `P_raw(-s)=2002 p(s)` and `T=s Q_raw(-s)` are correct;
the lower carry becomes a negative nonzero constant in `T`.

The referee reconstructs the full beta coefficients by a north/east path
dynamic program, then builds `H` symbolically. It verifies the exact
selected coefficient `[u^9]H=-s P_raw(-s)`, differentiates the literal
polynomial nine ways, and performs ordinary squaring and coefficient
extraction. Every full polynomial `J_r`, including its degree, every
remainder modulo `p`, and the full target agree with the certificate.
No coefficient or carry is omitted. Exact Sturm counting also confirms
the five roots of the displayed alternating beta polynomial are positive,
so the square-root pullback is real-rooted in the carrier variable.

## Four-dimensional certificate and branch scope

The four atoms have `(amplitude degree, derivative order)`
`(0,0),(1,7),(1,8),(3,2)`. Their positive normalizers are respectively
`1,32920473600,131681894400,5184`, exactly `(9)_r^2`.
At the four positive roots `r_i` of the actual first polynomial, let
`z_i=sqrt(r_i)>0`. Dividing each row by `J_0(r_i)<0` gives the stated
four-by-four interpolation matrix with first column exactly one.

Each listed positive rational bracket changes the sign of `p`, and the
four intervals are disjoint. Four distinct roots exhaust a quartic, hence
each root is simple and there are no missing roots. The referee refines
the frozen brackets by **64 additional exact bisections**, then computes
fresh square-root enclosures with denominator `2^210`, checking them by
squaring. This is an independent interval certificate, not floating-point
root approximation.

For the determinant calculation, the producer expands 24 signed
permutations with exact rational intervals. The referee instead rounds
every operation outward to denominator `2^240` and performs **interval
Gaussian elimination**, selecting only pivots whose intervals exclude zero.
Every true matrix operation is enclosed despite discarded dependencies.
The product of the pivot intervals, with row-swap signs, encloses the
true determinant. Additional bracket refinement compensates for this
different path's looser dependency bounds; no producer input is changed.

The independent enclosures verify all 36 original window signs, `T(r_i)<0`,
and the following exact rational bounds, shown as integer numerators over
`10^6`:

| Quantity | Lower | Upper |
|---|---:|---:|
| `det A` | -371065 | -371064 |
| replacement numerator 0 | -365451 | -365450 |
| replacement numerator 1 | -156491 | -156490 |
| replacement numerator 2 | -1845480 | -1845479 |
| replacement numerator 3 | -1916 | -1915 |
| `beta_0` | 984869 | 984870 |
| `beta_1` | 421733 | 421734 |
| `beta_2` | 4973466 | 4973467 |
| `beta_3` | 5162 | 5163 |

In particular `A` is invertible and all four Cramer ratios are strictly
positive. Their definitions produce a **single common coefficient vector**
at the four roots, rather than a separately chosen coefficient at each root.
Cramer's rule therefore proves the exact displayed representation.
Every original ray is negative at each root, and the amplitude powers and
normalizers are positive; the represented actual target is consequently
negative. The known target sign is not being rediscovered as a general
theorem. The new content is this mixed-window positive representation
despite the earlier single-window coefficient obstruction.

## Rational coefficients and all eight roots

The primitive numerator `11p(s)` reduces modulo three to a nonzero scalar
multiple of `s^4+2s+2`. The referee checks all three potential linear roots
and all nine monic quadratic divisors over that field. None divides it,
so the quartic is irreducible over the rationals. The norm of a root is
`1/11`, a nonsquare rational number. A square in the quartic field would
have square norm, so that root is not a square in its field. Adjoining
its square root doubles the degree to eight, proving `p(z^2)` irreducible.
SymPy's exact rational factorization independently confirms the degree-eight
irreducibility; the analytic argument does not depend on that extra check.

At a positive root `z_i`, every nonzero nonnegative odd-power atom is
strictly negative. Subtracting an alleged identity at `z_i` and `-z_i`
isolates twice the odd part, so such an identity on all eight roots cannot
contain any nonzero odd atom. This applies equally to negative integer
amplitude powers, since the roots are nonzero and parity is unchanged.

If all coefficients of an identity at one positive root were rational,
clearing a sufficiently large **even** power of `z` would produce a rational
polynomial divisible by the irreducible `p(z^2)`. It would then hold at all
eight roots, contradicting the odd-part sign argument. Thus the four
coefficients in the certified representation cannot all be rational, and
the stated rational positive odd-atom exclusion is valid. This does not
exclude an all-even representation, nor does it make any assertion about
uniform parameter dependence or minimal atom count.

The incoming quadratic-anchor result is complementary: it concerns a
different model class and its boundary. Choosing algebraic coefficients
from these four ordered roots supplies no finite-phase statement uniform
over that model. This audit accepts the fixed-row certificate independently
of that contextual comparison.

## Frozen evidence and reproduction

The referee imports SymPy and the standard library but no producer code,
repository theorem engine, numerical root solver, or optimization package.
Its finite universe is one genuine support, all 15 first/doubled channels,
all nine original windows, four positive roots, and five determinants.
Source and certificate are read beside the referee; the producer output
also supports its eventual `05-knowledge/results` filing location.

```
python -B 04-computation/continuing2_20260906_laurent_sparse_amplitude_audit.py
python -B -O 04-computation/continuing2_20260906_laurent_sparse_amplitude_audit.py
```

Both final runs pass **644 always-active gates**, with byte-identical LF
output. The final ten gates make the full target/window polynomial degrees
explicit; the earlier preliminary run had 634 gates.

```
referee source 6f196a8d7ea2e1086d1ed4514e63b9ee0a56f2f6ce4731cb987d05489e7c7814
referee output ded2bc088c181204eca38460e392db4068688a0371ee2e573055707c96b0ef06

producer source 2afc7df14b92e70a700a9a06a17b34ee2331eabde89902dcf408f597dccb9e28
producer output 7b95d9102c18309ce245d9017f3456b0f75447117f1ec6ac363393c63f183777
certificate     df8471af42220b5b26f7e4ff78dd9611d7999c0578e36603a6a7a33ba3c0afa2
```

No producer artifact, repository file, or Git state was changed.
