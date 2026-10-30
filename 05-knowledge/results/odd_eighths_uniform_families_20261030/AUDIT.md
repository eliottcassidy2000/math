# Validity and reproducibility

**PROVED scoped / FINITE-EXACT; same-agent proof review.** No independent
reviewer or Lean build is claimed. Source hashes certify inheritance bytes,
not mathematical correctness or priority.

## Proof audit

L1 requires actual complete windows throughout the stated radius chamber.
Closed endpoints are retained, and all tails use one common arrival time.
For C_q the hypothesis is inherited from Session114's Farey proof. L1 is
not a formula for radii at or below1/(q+2).

L2 assumes2m>s(q); the missing tight branch is not filled by the false
substitution m/(qm+s). The iff concerns preservation of a known maximum.
A forbidden attachment may still satisfy LRC at a smaller radius. Formula(5)
tests actual added integer speeds at both actual maximizing points. Repeated
speeds have their ordinary set interpretation.

L3 resolves ties by MAXIMUM speed, as required by physical radius relaxation.
It retains both directed sides, then uses reflection. The shrinking-mass
family changes the actual speeds while preserving residues and the maximum;
this is a hostile to a mass floor based only on the latter data.

E1 distinguishes parity-prefix cylinders from complete valuation cylinders.
The source address modulo2^A omits the final endpoint-oddness bit. Conditioning
always uses the same source. Cross-cost orthogonality is proved only where
the high conductor forces it. Finite-cost factorization precedes passage to
the infinite sums, justified by uniform convergence for each fixed m.

E2's exact bandwidth is a statement about this ordered Gram family. Its
nonzero outer edge explicitly decodes the conditioned source residue. It is
not a lower bound on every nonlinear representation and does not reconstruct
an unbounded integer from a fixed finite residue.

E3 applies the finite-support energy identity to S_M with M>=h. The coefficient
period retains c q^-M, and at h=3 its square-root sign also retains M parity.
The result concerns conditioned Haar expectation. Neither an almost-everywhere
limit nor these finite conditional identities are promoted to fixed-seed H1.

E4 tests product preservation on the FULL finite cyclic character group,
not merely a list of labels. Its forced unit comes from a generator. The
two-letter obstruction has the exact valuation, and all higher word terms
give the sufficiency. Odd offsets are allowed, and the plus/minus reflection
does not preserve the positive integer domain.

## Exact checker universes

`check.py` is standard-library-only. Checks use explicit exceptions and
remain active under Python -O. It does not read its expected receipt.

-60 conditional Gram universes: q=3,5,7,9; h=1..4; EVERY odd c modulo2^h;
  m,n<=6; halving cost cutoff T=8. Direct character pairs use forward shortcut
  replay on residues, while the predicted factorization uses suffix vectors
  and binomial prefix counts. Arithmetic is exact in Q[zeta_(2^h)], represented
  by rational coefficients modulo X^(2^(h-1))+1. No floating-point equality
  or FFT approximation is used.
-180 infinite-formula remainder checks compare the full Gram sum with the
  boundary-polynomial formula at M=h,h+1,h+3 in those60 universes. The all-M
  identity is proved analytically; these are algebraic controls, not independent
  enumeration of infinitely many characters.164 outer-band entries check
  the explicit nonzero source decoder. The h=1 diagonal and adjacent entries
  replay the inherited identities through m=12.
-252 full character-permutation universes: q,p in{3,5,7,9,11,13}, D=1..7.
  Every character is constructed.180 of these (D<=5) additionally test every
  ordered additive pair independently.30 unequal-multiplier pairs check the
  explicit first-depth antipodal defect.405 odd-offset universes use
  q,p in{3,5,7}, b in{-1,1,3}, c in{-1,1,5}, D=1..5.
-320 joint-tail first-arrival decisions: q=3..12, eight listed finite banks
  including singleton and common-divisor banks, and four exact radii per
  chamber. The actual full-speed survivor is compared with the separate
  tail first-arrival inequality.
-1432 arbitrary-speed attachment decisions cover q=3..12, m=2..4 satisfying
  the strict guard, and every v=1..2D. All surviving maximizing points are
  checked, not merely existence.84 multi-tail banks check the physical mass
  derivative by exact interval lengths at a small rational epsilon.12 further
  controls enumerate every opposite-slope rational clock for the GLOBAL
  maximum, independently of the top-window theorem.
-All64 labelled four-vertex tournaments replay the inherited3/8,1/8,1/8,3/8
  class masses; strong cases are checked by literal reachability.
  The original twelve-vertex81/27 control is also replayed from its pinned
  checker:67 common cyclic triples and regular binary counts1620/1728.
  No new tournament minimality or A000568 term is asserted.

Eight additional bridge controls at m=2,3,6,10 and both M parities recover
incoming115's physical harmonic-LRC mass from our conditional Gram energy.
The maximizing points, their tail phases and their nearby interval lengths
are computed literally, retaining multiplicity. This is an independent
target check of PHASE_BRIDGE.md, not a replay of the inherited count formula.

Normal and optimized runs produce identical JSON: **4778 local gates plus
seven inherited tournament gates**. Reproduce from the repository root:

```text
python -X utf8 -B 05-knowledge/results/odd_eighths_uniform_families_20261030/check.py
python -X utf8 -O -B 05-knowledge/results/odd_eighths_uniform_families_20261030/check.py
python -X utf8 -B 05-knowledge/results/odd_eighths_uniform_families_20261030/seal.py
```

Two runtime inputs are imported: the prior exact interval/clock enumerator
and the prior81/27 checker. The manifest pins them and fails if their live
bytes change. The former is used as an independent geometric computation,
not to generate expected answers from the new formulas.

## Scope and provenance

No claim that the source-address partition, h=1 energy, single-tail spectrum,
or elementary cyclotomic character was invented here. The new packet supplies
the all-depth conditional factorization, exact periodic remainder, arbitrary-bank
extension with its exact mass law, and the guarded phase-response decoder.
Incoming115 published the same multiplier-compatibility classification during
this session and receives precedence; E4 is a replay with an odd-offset
extension. Its harmonic-LRC mass formula is also inherited in the bridge.
These are same-agent proofs with independent finite
algorithms, not externally refereed major-conjecture resolutions.
