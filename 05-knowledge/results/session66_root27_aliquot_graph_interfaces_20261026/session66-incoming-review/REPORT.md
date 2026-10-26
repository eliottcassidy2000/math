# Incoming KP-S4: correction overlay and independent finite certificates

**Reviewed commit:** `44abcaf22b78125a11e039054f0d0ba9e4978f08`.
**Reviewed note:** `05-knowledge/results/sheets_root27_rank_untouchables_20261025.md`.
**Date:** 2026-10-26. **Status:** independent, scoped audit; incoming historical files are preserved unchanged. This overlay, not the unqualified historical claims, is the endorsement boundary.

## Verdict and retained contribution

The incoming note usefully connects shifted modular tables, source-specific carry cylinders, a Haar extinction probability, and distinct-Goldbach witnesses. Its **root-27 target, a second odd untouchable, and existence/value of the untouchable density remain OPEN**. The carry reformulation does not provide a decreasing global rank.

- **PROVED after qualifications:** the carry/LCP statement, full-match valuation formula, fixed-size frontier extremality, and Haar probabilities below.
- **CITED plus a finite exact bridge:** every odd target `9 <= m <= 4*10^18+1` is touchable. This implies the incoming weaker bound, but it uses more than ordinary Goldbach alone.
- **FINITE-EXACT independently reproduced here:** `150232` untouchables through `10^6`, all six incoming power-of-ten counts, the first forty values, the subtotal `318` through `3000`, and all nine residue tables. The count is obtained without floating arithmetic or any producer import.
- **HISTORICAL / NOT independently replayed here:** `443731` non-semiprime-source odd targets, `56267` semiprime/prime-square-only odd targets, the `2509482` composite-pair count, and the incoming long orbit/halting censuses. No counterexample to these numerical values was found; they are not silently promoted by our narrower replay.
- **REFUTED inference:** a finite positive share of semiprime-only *nonexceptional* targets does not rule out coverage of all distinct-Goldbach exceptions by other sources. “Goldbach-complete in practice” is not an established equivalence or complexity reduction.

Theorem A has a separate independent supplement at `session66-audit-modular/incoming/APPENDIX.md`. It proves the proposed prime-order compatibility criterion, repairs necessity using the constant term rather than cancelling `3`, states `m>=3`, and diagnoses the incoming `T1b` unit-translation test. That supplement owns the detailed modular audit; the present report does not duplicate its code.

## 1. Provenance and audit method

`SOURCE_PINS.json` pins the exact commit, original repository path, Git blob ID, SHA-256, and preserved local snapshot for the note and four B/C computational files. `check.py` verifies both hashes offline before testing. No original incoming text, experiment, or output was edited, no repository commit was made, and no external message was sent by this auditor.

Inheritance: the closest mechanism is the carry-cylinder equivalence, not a new root-27 proof. The hostile is the fully matching positive source `z=21` with word `1000`. The corrected near miss is mistaking a canonical frontier presentation for a universal monotone rank. The least-used sidecars are the residual state valuation, fixed tree size, and the distinction between a Goldbach exception and a semiprime-only value. Concept board: **carry precision / legal tree frontier / extinction measure / source-complete inverse census / target-support certificate**.

The local processor status check could not run because this sparse checkout has no `.machine-id`; no identity was invented or registered. This has no effect on the pinned mathematics or offline checks. Git publication remains the parent session's responsibility.

## 2. Theorem B: exact repair

Line numbers in this section refer to the sealed incoming note.

### 2.1 Matching words require the residual valuation (lines 105–110)

For a positive integer source `z`, let `kappa(z)` be its carry word under
`T(z)=(3z+(z mod 2))/2`. For a word `w` of length `n`, set
`Q_z(w)=3^n z+sum_i w_i 2^i 3^(n-1-i)` and `h=LCP(w,kappa(z))`.
The exact piecewise formula is

```
if h < n:   v2(Q_z(w)) = h,       rho_z(w) = n-h;
if h = n:   Q_z(w) = 2^n T^n(z),  rho_z(w) = -v2(T^n(z)).
```

The first branch follows because the first differing bit has an odd coefficient after its power of two is removed; later bits contribute only higher powers. The second branch follows by iterating the affine update. Thus `max(rho,0)=n-h` is determined by the truncated LCP, but the untruncated rank is not. One must retain word length even in the mismatch branch, and residual state valuation in the matching branch.

**Minimal complete-tree hostile:** `w=1000` matches both `z=5` and `z=21`, with `n=h=4`, but

```
Q_5(1000)=432=2^4*27,   rho_5(1000)=0;
Q_21(1000)=1728=2^6*27, rho_21(1000)=-2.
```

The incoming numerical checker already correctly separates matching and mismatching branches (rank script lines 138–144); it does not support the stronger prose “only through the first disagreement.” For general 2-adic sources, allow `v2(0)=+infinity`; the finite script's `10^9` sentinel is not a mathematical valuation at zero. All arithmetic tests here use positive integer sources, so `Q>0`.

### 2.2 Completed frontier versus executed prefix (lines 127–130)

A legal proper prefix `u` leaves `q=1+2#1(u)-#0(u)>0` pending slots. The word
`F=u 0^q` is **always a complete full ternary-tree word**: the prefix never exhausts its slots, and its appended zeros exhaust them precisely at the end. The prefix `u` *alone* is complete exactly when its credit first reaches zero.

For example, `u=1`, `q=3`, and `F=1000` is a complete tree with positive current credit. Hence the sentence claiming the completed frontier is a complete tree iff `q=0` is false. Replace it by the distinction above. The certificate equivalence survives:

```
A full tree w has v2(Q_z(w)) >= |w|
iff w is a full-tree prefix of kappa(z)
iff the source credit walk first reaches zero at the end of that prefix.
```

Formal processing of that final leaf and native TM execution are different operations: the separate root-27 audit retains the blocked-final-leaf guard and does not count it as an executed native macro.

### 2.3 Frontier formula also has a full-match boundary (lines 113–114)

For `F=u0^q`, put `t=|u|`, `n=t+q`, and
`r=LCP(0^q,kappa(z)[t:])`. Then

```
if r < q: rho_z(F)=q-r;
if r = q: rho_z(F)=-v2(T^n(z)).
```

At `z=21,u=1,q=3`, the incoming expression gives `q-r=0`, whereas the actual rank is `-2`. On source `7`, the first branch is valid at every noncertificate frontier; asserting it unconditionally for every frontier would silently need control of the open certificate boundary. The **clipped** rank formula `max(rho_z(F),0)=q-r` is unconditional.

### 2.4 A precise survivor: unique optimality at fixed size

The sentence “No sequence of shape changes can lower the rank below the frontier's” lacks a size/cut quantifier. Here is a valid theorem preserving its useful content.

**Fixed-size frontier theorem.** Suppose the source has not halted before executing its `j`th internal-node symbol. For `j=0`, take the empty prefix. Otherwise let `u` end immediately after that `j`th `1`, with all earlier symbols on the legal carry ray. Let `F_j=u0^q`, so `|F_j|=3j+1`. Among all complete ternary trees with exactly `j` internal nodes, `F_j` is the unique minimizer of `rho_z`.

**Proof.** The prefix `u` contains all `j` ones of `F_j`; afterward it contains only zeros. If `F_j` mismatches the carry word, its first mismatch is therefore `0` versus `1`. Any different tree with `j` ones must already differ from `u`: a tree agreeing with all of `u` has used every one and is forced to equal `F_j`. At that earlier difference `F_j` still agrees with the source, so the other tree has strictly shorter LCP and strictly larger rank. If `F_j` fully matches, every other same-size tree has a mismatch and positive rank, whereas `rho_z(F_j)<=0`. This includes the first halt-certificate frontier; no continuation past a halted tree is asserted. QED.

This is **not** an arbitrary-size monotonicity theorem. At source `7`, the canonical size-one tree `0` has rank `1`, and the next canonical tree `1000` has rank `3`. The same scalar can increase under legal refinement. The fixed-size extremum is a useful shape normalization, not a well-founded global descent or proof of root-27 nonhalting. Carry reformulation also does not prohibit new proof methods using structure of the fixed carry word; “the only lever left is the source” is too strong as a methodological no-go.

### 2.5 Haar probability survives, with the least-fixed-point clause

Carry cylinders of length `j` biject with residues modulo `2^j`, so Haar-random 2-adic source bits become independent fair carry bits. A pending slot is replaced by either zero or three slots, each with probability `1/2`. Its extinction probability is the **least** fixed point of
`p=(1+p^3)/2`, obtained by iterating from `p_0=0`; the equation alone also admits `p=1`. The least fixed point is `(sqrt(5)-1)/2`.

An odd source forces the first `1`, leaving three slots, giving `p^3=sqrt(5)-2`. The cylinder with prefix `1110` leaves six slots, giving `p^6=9-4sqrt(5)`. These are exact Haar statements, not a natural-density theorem for positive integer sources, and not a probability assertion about the fixed source `7`. The sibling root-27 audit independently checked these qualifications and the tiny hostiles.

## 3. Theorem C: retained reduction and exact cited boundary

For odd `m>=9`, a distinct-prime partition `m-1=p+q` gives the genuine source `pq`, since `s(pq)=1+p+q`. Therefore an odd untouchable `m>=9` is necessarily a distinct-Goldbach exception after subtracting one. Conversely, such an exception is untouchable **only if every other source family also misses it**. This is the correct reduction; neither direction claims that Goldbach is equivalent to uniqueness of the odd untouchable `5`.

**Primary provenance.** Booker's *Finite connected components of the aliquot graph*, [arXiv v2, PDF page 2, footnote 2](https://arxiv.org/pdf/1610.07471v2), explicitly distinguishes ordinary Goldbach from distinct-prime Hypothesis G and records the bound `p<=9781` over even `6<=e<=4*10^18`. The original computation page defines the minimal first prime, records the `9781` record, and dates completion of the range to April 4, 2012: [Oliveira e Silva's verification page](https://sweet.ua.pt/tos/goldbach.html), “Introduction,” “News,” and “Computational results.” These are cited large computations, not computations rerun here.

For `e>19562`, `p<=9781<e/2`, so that partition is distinct. Our integer checker separately exhibits a distinct pair for **every even `8<=e<=19562`**, not just an ordinary Goldbach pair. Thus all odd `9<=m<=4*10^18+1` are touchable; `1,3,7` have sources `2,4,8`. Any second odd untouchable is consequently at least `4*10^18+3`. This strengthens, rather than retracts, the incoming weaker bound `>4*10^18`.

## 4. Census precision: flaw, exact repair, and unreplayed data

### 4.1 The raw FFT is not an exact certificate

The incoming census script lines 86–98 convert the prime indicator to `float64`, apply `np.fft.rfft` and `irfft`, round with `np.rint`, and subtract the prime diagonal before testing positivity. It supplies no proved forward-error bound, interval enclosure, exact convolution, or complete support witness certificate. Integer *storage after rounding* does not alone justify “every count an exact integer” as a proof statement. Passing OEIS comparisons and a small brute-force check does not repair an unbounded rounding claim.

The mathematical case split itself is sound for `N=10^6`: writing a composite source `n=pk` with `p=spf(n)` gives `k<=N-1`; the script covers small `p`, prime `k`, and composite `k` exhaustively. Its non-semiprime marking path is integer-based. The identified gap is the unqualified computational certification, not a demonstrated incorrect count.

### 4.2 Independent exact support and census

`check.py` uses an exact prime sieve and direct distinct-prime witnesses for all **499997 even sums from 8 through 1000000**, with largest least first prime **523**. Odd sums have distinct-prime support exactly when `e>=5` and `e-2` is prime. This gives the entire Boolean support relevant to the incoming convolution, without needing its pair multiplicities.

The independently implemented even-target census uses the classical parity/source bounds from [Pomerance–Yang, Section 5, PDF page 8](https://math.dartmouth.edu/~carlp/uupaper3.pdf): even sources lie below `2N`; odd sources must be squares. Prime square roots supply `p+1`; composite roots satisfy `a^3<=N^2`. We enumerate even sources through `2000000` and the **3771 odd composite roots through 10000**, by integer divisor sums and factorizations. Odd targets are supplied by the exact pair witnesses, with the small cases handled directly. A separate complete divisor sieve checks all targets through `200` against the independent source bound `(200-1)^2`.

This recovers the incoming total **150232**, split **150231 even and one odd (`5`)**, every listed residue count, all six power-of-ten totals, the first forty values, and the `3000` subtotal. The finite ratio is exactly `150232/1000000`; no limiting density is inferred. The source bounds and parity algorithm are prior art, not a new algorithmic priority claim.

`fft_support_diagnostic.py` additionally replays just the incoming floating support layer and compares it entry-by-entry against independently exact support: **zero mismatches** in both normal and `-O` runs. Its observed maximum rounding residual, about `1.09e-11`, is diagnostic only. Agreement with the exact support is the certificate for this replay; small numerical residual is not the proof. No claim about all convolution multiplicities or a general FFT error bound is made.

### 4.3 Remaining output boundary

The counts `443731`, `56267` and `2509482` remain sealed historical incoming outputs, not newly replayed independent results. The first two classify sources, which our target-support census deliberately does not enumerate completely. In particular, verifying total touchability does not verify the assertion that **every** preimage belongs to a named source family. The finite `11.25%` classification must retain this provenance and must not be used as an asymptotic or exception-set obstruction.

## 5. Rejected untouchable conclusions and replacement wording

Incoming lines 207–209 claim that a finite share of semiprime-only values rules out uniform non-semiprime coverage of distinct-Goldbach exceptions. This changes the quantified set.

Let `E` be the even numbers at least eight with **no** distinct-prime partition, and let `B` be the odd targets at least nine reached by a non-semiprime source. The unresolved alternative route is `E+1 subseteq B`. The reported semiprime-only targets are already touchable by distinct semiprimes, and consequently lie **outside `E+1`**. For example, `9-1=3+5`, `11-1=3+7`, and `17-1=3+13`. Their lack of other sources says nothing against `E+1 subseteq B`. Even positive *asymptotic* density of nonexceptions without other sources would not refute coverage of `E+1`.

**Replacement:** “A second odd untouchable must be a distinct-Goldbach exception with no remaining source. Proving all distinct-Goldbach exceptions are hit by other source families is another sufficient route to uniqueness of `5`; the computed dependence share does not rule out that route. Neither route is completed here.”

Also replace the density-coding explanation at lines 210–212 by: “No tractable measure-preserving coding for the aliquot target event is constructed here. The contrast with the carry-cylinder model does not prove that such a coding cannot exist.” For the shifted Chowla statement, use **positive lower density**, not an existing exact density. No asymptotic density theorem follows from the shift or this census.

## 6. Reproduction and final status

From this directory:

```sh
python3 -B check.py > normal.json
python3 -B -O check.py > optimized.json
cmp normal.json optimized.json
# Optional: NumPy diagnostic, not needed for the integer proof.
python3 -B fft_support_diagnostic.py > fft-normal.json
python3 -B -O fft_support_diagnostic.py > fft-optimized.json
cmp fft-normal.json fft-optimized.json
```

The exact checker has no third-party dependency, no producer imports, no `assert`-only gate, and no network dependency. It checks **32704** binary words over the declared source/length universe, including **288** full matches, and **168** fixed-size frontiers (odd sources `1..63`, at most five internal nodes). Paired final receipts and hashes accompany this report. The report's fixed-size theorem and Haar proof are all-parameter arguments; these finite tests are adversarial controls, not their proofs.

**Scope closed:** exact incoming pin preservation, identified proof/logic repairs, primary distinctness citation, and independent bounded numerical endorsement. **Still open:** root-27 survival/decreasing global rank, any second odd untouchable, all-height source-family coverage, and the actual untouchable density. No historical producer mutation or Git publication was performed by this audit.
