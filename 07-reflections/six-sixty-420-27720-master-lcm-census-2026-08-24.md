# The hidden `6, 60, 420, 27720` master-LCM staircase — 2026-08-24

## Truth surface

This note classifies a repository-wide search for the exact natural numbers

```text
6, 60, 420, 27720.
```

The central identity is **PROVED**:

\[
L_m=\operatorname{lcm}(1,\ldots,m),\qquad
(6,60,420,27720)=(L_3,L_5,L_7,L_{11}).                 \tag{1}
\]

Thus the indices are the first four odd primes.  This is a sampled
initial-segment LCM staircase, not a four-term recurrence.  The full staircase
between the endpoints is

```text
n:    3   4   5   6    7    8     9    10    11
L_n:  6  12  60  60  420  840  2520  2520  27720.       (2)
```

The sample skips the genuine jumps at `12`, `840`, and `2520`; its successive
ratios are `10,7,66`.  Any account which treats it as multiplication by a
fixed or simply prime-valued factor is false.

The executable certificate is
`04-computation/six_sixty_420_27720_master_lcm_census_20260824.py`, with frozen
output in
`05-knowledge/results/six_sixty_420_27720_master_lcm_census_20260824.out`.
It verifies all identities below and independently compares the two formulas
for every prime through `10000`.

```bash
python 04-computation/six_sixty_420_27720_master_lcm_census_20260824.py
python -O 04-computation/six_sixty_420_27720_master_lcm_census_20260824.py
python 04-computation/six_sixty_420_27720_master_lcm_census_20260824.py --scan-corpus
```

## Inheritance pass

- **Closest proved mechanism:** THM-1415,
  `01-canon/theorems/THM-1415-switching-is-the-canonical-star-quotient.md`,
  already identified subset-LCM attraction as the reason small arithmetic
  constructions pile onto `60` and `420`.
- **Canonical hostile example:** THM-4042,
  `01-canon/theorems/THM-4042-prime-sector-ap-cover-exact-clock-and-holonomic-law.md`,
  has the requested four clocks at sectors `5,7,11,13`, but the apparent law
  fails immediately at sector `17`.
- **Corrected near miss:** a raw token match is overwhelmingly contaminated by
  repeated rational denominators, enumeration addresses, theorem numbers, and
  coefficients.  Co-occurrence alone is not a connection.
- **Least-used relevant sidecar:** THM-1105,
  `01-canon/theorems/THM-1105-arithmetic-position-law.md`, uses the same `L_Q`
  in the dual role of annihilating every rational denominator `q<=Q`.

## I. Four full parametric carriers already in the repository

### A. Prime-sector AP owner clocks

THM-4042 is **PROVED + FINITE-EXACT**.  Its exact prime-sector clock is

\[
\Pi_P=\prod_{\ell\le P-2}
 \ell^{\min(\lfloor\log_\ell(P-2)\rfloor,v_\ell(P-1)+1)}.             \tag{3}
\]

Writing

\[
M_P=L_{P-2},\qquad R_P=\operatorname{rad}(M_P),
\]

the census exposes the equivalent master-LCM law, now recorded in THM-4042:

\[
\boxed{\Pi_P=\gcd(M_P,(P-1)R_P)
             =R_P\gcd(M_P/R_P,P-1)}                    \tag{4}
\]

for every prime `P>=5`.  Hence

\[
\Pi_P=M_P\quad\Longleftrightarrow\quad M_P/R_P\mid P-1.             \tag{5}
\]

At the four rows in question,

```text
P                 5       7       11       13
Pi_P              6      60      420    27720
previous prime    3       5        7       11
L_previous        6      60      420    27720.          (6)
```

This explains exactly why the sequence appeared in the current AP work.  It
also gives a decisive hostile control.  At `P=17`,

```text
M_17=L_15=L_13=360360,   R_17=30030,
Pi_17=30030*gcd(12,16)=120120=L_13/3.                   (7)
```

The second factor `3` in `M_17` is destroyed because `3` does not divide
`P-1=16`.  Only its radical copy survives.  Therefore the previous-prime law
and the divisibility tower both fail at the next row.  The exact audit through
all primes `P<=10000` found `Pi_P=M_P` only at `P=5,7,13`.  This is
**FINITE-EXACT**, not an all-prime classification of the equality set.

### B. The LRC plateau constant

THM-1420,
`01-canon/theorems/THM-1420-the-sixty-plateau-in-lrc.md`, gives the structural
LRC constant

\[
L_{\lfloor n/2\rfloor}.                                \tag{8}
\]

The same family contains all four requested values without having been read
as one sequence:

```text
clock       exact n-band
6           6..7
60          10..13
420         14..15
27720       22..25.                                     (9)
```

The stored artifact
`05-knowledge/results/sixty_in_lrc_kps_S128c111.out` also prints the covering
modulus ladder `Q=3 -> 6`, `Q=5,6 -> 60`, `Q=7 -> 420`, and
`Q=11,12 -> 27720`.  THM-1420 is a **VERIFIED-EXACT deflation**, not a proof
of LRC(14); its value here is the parametric replacement of an isolated `60`.

### C. The rational-denominator blocker

THM-1105 contains a **PROVED exact refutation** of a bounded-witness-denominator
proposal.  A speed `L_Q` is zero modulo every `q<=Q`, so it blocks every
rational witness of those denominators.  Instantiating `Q=3,5,7,11` gives
exactly `(6,60,420,27720)`.

This is dual to the AP owner clock but not identical to it:

```text
source             THM-4042 owner phase
target             THM-1105 blocker speed
map                 Z -> Z/qZ for every q below the cutoff
preserved           common quotient/divisibility address
destroyed           owner weights, pole shifts, gap geometry, sign
needed sidecar      AP owner word on one side; danger interval on the other
cheapest test       compare prime-adic clock exponent with v_q(L_Q).       (10)
```

HYP-4040,
`05-knowledge/hypotheses/HYP-4040-dyadic-magnitude-ladder-discrepancy.md`,
pushes the same theorem into a magnitude lower bound: `S_X` contains `L_X`,
every `q<=X` is blocked, the least lonely denominator exceeds `X`, and the
largest speed is `e^{psi(X)}`.  Its lower-bound component is **PROVED**.

THM-1098,
`01-canon/theorems/THM-1098-primitive-finite-atlas-height-obstruction.md`,
THM-566,
`01-canon/theorems/THM-566-lrc14-covering-sets-have-no-uniform-bounded-denominator-witness.md`,
and THM-1100,
`01-canon/theorems/THM-1100-extended-rational-point-sieve.md`, are unscaled,
scaled, and finite-atlas members of this same divisor-loading lineage.  They
are genuine appearances of the cutoff LCM but should not be counted as
independent mechanisms.

### D. The Möbius/Ramanujan `v`-grid moment

THM-879 must be cited by its collision-safe slug,
`01-canon/theorems/THM-879-vgrid-moment-and-moebius-sinc.md`.  It is
**PROVED** and gives the exact restricted second moment

\[
S_v(L)=\sum_{d,e\le L}
 {M(\lfloor L/d\rfloor)M(\lfloor L/e\rfloor)
  \over\operatorname{lcm}(d',e')},\qquad
d'={d\over\gcd(d,v)}.                                  \tag{10a}
\]

If `gcd(v,L_L)=1`, then `d'=d` for every `d<=L`, so `S_v(L)=S(L)` exactly.
If `v` is resonant, divisor depths collapse and the moment is amplified.  The
generic equality clause alone factors through `rad(L_L)`, but the full
function `v -> S_v(L)` does not: it has exact minimal period `L_L`.  Indeed,
`S_0(L)=1`, while for each prime `p|L_L`, the largest-`p`-power block in
`(10a)` gives

\[
S_{L_L/p}(L)=\frac1p.                                  \tag{10b}
\]

Thus no `L_L/p` is a period, and every proper divisor of `L_L` divides one of
those nonperiods.  Therefore the exact minimal periods at cutoffs
`L=3,5,7,11` are

```text
L_L = 6, 60, 420, 27720.                               (10c)
```

This is a fourth full carrier of the cutoff sequence in a genuinely different
Möbius/Ramanujan mechanism.  The target numbers are periods of the whole
moment function, not values of the moment.  Its connection contract is:

```text
source             divisor mode d
target             v-grid mode d/gcd(d,v)
preserved           every mode and the exact second moment when gcd(v,L_L)=1
destroyed           divisor depth under resonance
needed sidecar      gcd(d,v), not only the grid period
cheapest test       compare S_v(L) with S(L) before asymptotics.            (10d)
```

### E. The valuation-depth sidecar at a high cutoff

THM-3848,
`01-canon/theorems/THM-3848-rational-base-prefix-atom-tree-and-lonely-runner-separation.md`,
provides a sharp continuation of the same filtration.  Put

```text
P=product_(prime ell<=121) ell,   L=L_121,
U=84P,                            V=84L.                (10e)
```

Then `rad(U)=rad(V)=P`, but `U=55 mod 121` while `V=0 mod 121`.
The row containing `U` is strictly lonely at `10/121`; every denominator at
most `121` divides `V`, so none can witness its row.  This isolates the exact
coordinate which the quartet carries beyond its primorial shadow:
prime-power valuation depth.

That is the same sidecar exposed algebraically by `(4)`:

\[
\Pi_P=R_P\gcd(M_P/R_P,P-1).                            \tag{10f}
\]

Here `R_P` is the radical/primorial address and `M_P/R_P` is the valuation
depth.  At sector `17`, loss of one `3` from that second coordinate causes the
clock collapse; in THM-3848, retaining the full second coordinate separates a
lonely phase from a completely blocked denominator bank.  This is a genuine
cross-frontier identification of the *missing coordinate*, not merely a
shared number.

## II. The previously unnoticed `3960 -> 27720` phase square

The most surprising cross-session match is not another literal copy of the
whole quartet.  It is the same intermediate clock `3960` arising in two
independent LRC constructions and receiving the same extra sector factor `7`.

### A. One-far exact periodicity

THM-563,
`01-canon/theorems/THM-563-signed-deltaw-periodicity-bound.md`, proves

\[
P(B)=7\operatorname{lcm}(B)                            \tag{11}
\]

for the signed one-far deviation.  The old general-period computation
`04-computation/lrc_periodmax_general_macmini_0621s6.py` and its stored output
contain both

```text
B={0,1,2,3,4,5,6}:       lcm(B\{0})=60,    P(B)=420;
B={0,3,5,6,8,9,11}:      lcm(B\{0})=3960,  P(B)=27720.  (12)
```

The second row is also echoed in
`05-knowledge/results/lrc14_periodmax_type_bridge_codex_20260621.out`.
These are genuine exact phase clocks, not recycled display denominators.

### B. The translated hostile family

THM-1226,
`01-canon/theorems/THM-1226-gcd-period-projective-charge-obstruction.md`, uses

```text
R={1,2,3,5,7,11,13}.
```

Its positive difference set is exactly

```text
{1,2,3,4,5,6,8,9,10,11,12},
```

whose LCM is `3960`.  The hostile translation must also satisfy
`A == 0 (mod 14)`, so the least joint modulus is

\[
\operatorname{lcm}(3960,14)=27720.                    \tag{13}
\]

THM-3451,
`01-canon/theorems/THM-3451-prime-quarter-half-twist-cap-seven-classification.md`,
contains the smaller parallel join

\[
\operatorname{lcm}(60,14)=420.                        \tag{14}
\]

Together these make an exact commuting numerical square:

```text
60    -- join with 14 -->   420
 | x66                         | x66
3960  -- join with 14 -->  27720.                      (15)
```

Both base clocks have gcd `2` with `14`, so the join contributes exactly a
new factor `7`.  What is preserved is the common cyclic address.  What is
destroyed is crucial: THM-563's endpoint-sawtooth phase and THM-1226's
pair-gcd/projective-charge obstruction are different observables.  A future
transfer would need the signed endpoint ledger on one side and the primitive
pair labels on the other.  The cheapest decisive test is to evaluate the
THM-1226 difference packet in the THM-563 sawtooth functional; equality of
periods alone is already exhausted as evidence.

## III. Other genuine LCM clocks and meshes

### A. Consecutive CRT observer

THM-4000,
`01-canon/theorems/THM-4000-centered-base-split-cubic-observer-and-tripotent-crt-atlas.md`,
uses

\[
\operatorname{lcm}(6,7,8,9,10,11)=27720.              \tag{16}
\]

The residue observer lives on `C_27720`.  Exact integer samples instead retain
the product `6*7*8*9*10*11=332640=12*27720`; quotienting to the LCM loses CRT
lift/Smith carry data.  This is a genuine clock but only the last requested
value.  Its separate cubic fact `gcd_z(z^3-z)=6` is not the same mechanism.

### B. Root-of-unity observer lattice

THM-4044,
`01-canon/theorems/THM-4044-sixty-clock-hasse-alias-and-planar-jc-boundary-firewall.md`,
is **PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED**.  Its proof now records
the uniform corollary

\[
\ker O_{k,M}=((P^M-1)^k),\qquad
\ker(O_{k,M}|_{P^2K[P]})=P^2(P^M-1)^kK[P].             \tag{17}
\]

For `M|N`, root sets include forward and kernels include backward.  Since

```text
6 | 60 | 420 | 27720,
```

the requested values form a genuine nested observer lattice, with first pure
residual alias degree `Mk+2`.  The map preserves roots of unity and jet depth;
it destroys AP owners, LRC danger arcs, and the Keller equations.  This is a
useful observer analogy, not a transfer of a counterexample.

### C. Endpoint and fibre grids

- THM-2258,
  `01-canon/theorems/THM-2258-depth-three-uniform-five-charge-spectrum-overlap-exclusion.md`,
  has a genuine common endpoint mesh
  `14*lcm(33,3,44,55,165,198)=27720`.  The endpoint `7608/27720` reduces to
  `317/1155`; the large denominator is the shared grid, not the reduced point.
- THM-3106,
  `01-canon/theorems/THM-3106-projected-k3-z232-exact-screen-and-complete-cell-cardinality-descent.md`,
  has ruler `E=(1,9,10,11,12,14)`, `lcm(E)=13860`, full grid
  `14*lcm(E)=194040`, and seven-fibre step `194040/7=27720`.  THM-2928,
  `01-canon/theorems/THM-2928-critical-seven-comb-grid-tensorization-and-drift-tariff.md`,
  uses the same `194040/7` capacity.  These are one grid family, not
  independent discoveries.
- THM-823,
  `01-canon/theorems/THM-823-hamming-five-common-sheet-deck-boundary.md`,
  hides the same endpoint immediately below its replay modulus:
  `360360=lcm(2,...,12)*13=27720*13`.  This is an exact clearing-clock
  continuation, not another full carrier.
- THM-720,
  `01-canon/theorems/THM-720-looseness-dichotomy-quantified.md`, and
  MISTAKE-096 in `01-canon/MISTAKES.md` use `27720` as a smooth finite-cover
  milestone.  The mistake entry is the important survivor: no fixed
  denominator cover can become a uniform LRC certificate.
- THM-2135,
  `01-canon/theorems/THM-2135-root-profile-invoices-and-first-deep-scalar-tail-closure.md`,
  uses the fact that every integer at most `7` divides `420`; THM-960,
  `01-canon/theorems/THM-960-scale-six-hamming-six-affine-nerve-obstruction.md`,
  has an exact hereditary order-`6` clock.  THM-2367,
  `01-canon/theorems/THM-2367-septimal-root-averaging-graft-and-cover-alignment.md`,
  and HYP-2631,
  `05-knowledge/hypotheses/HYP-2631-lrc14-exact-period-ap-drop-repair.md`, use
  exact `420` endpoint grids.  These are genuine partial members, not full
  sequence carriers.
- HYP-3877,
  `05-knowledge/hypotheses/HYP-3877-regime-C-arithmetic-band-route.md`, uses
  `lcm(2,...,11)=27720` as the exact finite magnitude threshold above its
  bounded arithmetic bank; MISTAKE-095 in `01-canon/MISTAKES.md` corrects the
  scope of the extrapolation.  Again the value is a cutoff, not a uniform
  certificate.

The exact resonance ledger
`05-knowledge/results/lrc_resonance_ledger_kps_S73.out` supplies an especially
useful hostile control: two families identical modulo
`420=lcm(2,...,7)` have different danger measure and different per-denominator
attribution.  Thus even a genuine `C_420` address forgets the statistic that a
loneliness proof needs.

### D. Historical CRT staircase

The historical exact artifact
`05-knowledge/results/two_three_tower_deep.out` prints the accumulator

```text
60 -> 420 -> 840 -> 2520 -> 27720 -> 360360.            (18)
```

This comes from repeatedly imposing `I(b) == I(-1) (mod b+1)`.  It is a true
CRT/initial-LCM staircase and independently confirms that `60,420,27720` are
not isolated constants.  Its historical HYP-920/921 labels are not promoted
here into current canon.

## IV. Why these numbers are statistically attractive

THM-1415's exact subset census on nonempty subsets of `{1,...,12}` of size at
most `6` has `2509` samples:

| target LCM | count | frequency rank |
|---:|---:|---:|
| 6 | 10 | 53 |
| 60 | 123 | 2 |
| 420 | 97 | 5 |
| 27720 | 15 | 44 |

Over all `4095` nonempty subsets, the counts are respectively
`10,132,132,192`.  This quantifies an important prior: LCM constructions on
small integers are naturally pulled toward these plateau values.  It does not
make distinct observables equivalent.  In particular, THM-1415 already
firewalls three unrelated `60`s: the Pisano period, `ord_1001(2)`, and
`|A_5|` have different causes; `A_5` even has exponent `30` and no element of
order `60`.

## V. Partial overlaps and hostile boundaries

### Sun's structured hostile scales

THM-4036,
`01-canon/theorems/THM-4036-sun-2468-energy-and-support-exponent.md`,
chooses `t=24P_k`, with `P_k` the product of the first `k` odd primes.  Its
first scales are

```text
72, 360, 2520, 27720, 360360.                           (19)
```

The `27720` overlap is exact but only finite: this is a primorial-product
construction, while `(1)` is an initial-segment LCM with prime-power
envelopes.  The shared value occurs at `24*(3*5*7*11)` and should be labelled
**PARTIAL OVERLAP**, not a common recurrence.

### Fibonacci and triangular work

The Fibonacci material has exact Pisano prefix matches `pi(4)=6` and
`pi(10)=60`.  The repository has no corresponding recorded Pisano periods
`420` and `27720`; a `420` in the tiling graph is an edge count.  Thus the
Fibonacci connection stops after two terms.  The triangular/Pell lane has an
isolated `420` in HYP-2454,
`05-knowledge/hypotheses/HYP-2454-triangular-power-balance-towers.md`, but no
full quartet.  These are useful hostile boundaries: a two-term clock match is
weak evidence.

There is one further exact but finite Fibonacci/Farey overlap in
`07-reflections/fibonacci-farey-q15-four-state-and-mod6-24-state-bridge-codex-20260814.md`:
the canonical six-speed trimode edge has

\[
1+\frac12+\frac13+\frac14+\frac15+\frac17=\frac{1019}{420}.          \tag{19a}
\]

This denominator is exactly `L_7`; it is a finite root mass, not a Pisano
period or a Fibonacci tree-index harmonic mass.

The tropical Morgan--Voyce product `3*8*21*55=27720` is also isolated; its
neighbours are `1512` and `498960`.  It is an exact factorization, not this
staircase.

The open reference note
`05-knowledge/reference/apery-style-irrationality-framework.md` uses
`B_n=L_n^m b_n`, `L_n=lcm(1,...,n)`, as the denominator-clearing balance.
Formally it therefore carries the quartet at `n=3,5,7,11`, but the targeted
Apéry recurrences are an **OPEN framework**, not proved results.

## VI. Raw search and contamination audit

The semantics-free scan covered every regular file below `04-computation` and
`05-knowledge/results`, excluding the census script and its frozen output:

| token | exact tokens | matching lines | matching files |
|---:|---:|---:|---:|
| 6 | 723417 | 500250 | 24805 |
| 60 | 24825 | 21432 | 5341 |
| 420 | 5695 | 3035 | 870 |
| 27720 | 2028 | 1901 | 101 |

The universe has `38138` files and `638951186` bytes; `38` files contain all
four exact tokens.  Yet `1798/2028 = 88.658777%` of all `27720` tokens are the
literal denominator string `/27720`.  Four large LRC ledgers account for most
of that repetition.  This is why raw co-occurrence was used only to generate
candidates and every promoted connection was then reconstructed from its
formula.

Rejected false positives include:

- `bits=27720`, `T#27720`, and JSON mask `27720`: enumeration addresses;
- AMM coefficient `+/-27720 = binom(56,3)`: an unrelated coefficient;
- `ip_3cycle_fast.out`: self-labelled erroneous output;
- PSL(2,7) `420`: an edge count, not a cyclic action of order `420`;
- the corrected Paley `T_11` count `c_7=3960`, hence `7c_7=27720`, and
  THM-2319,
  `01-canon/theorems/THM-2319-crt-unit-bispectrum-needle-and-mixed-polarization-no-go.md`,
  with its unrelated CRT count `3960`: neither has a map preserving the LCM
  phase observable;
- incoming THM-4049,
  `01-canon/theorems/THM-4049-lrc14-d2-two-phase-residue-firewall.md`, uses
  `60` both as one lift-time numerator over denominator `112` and as a finite
  histogram count.  Neither use is an LCM clock;
- the seductive prefix `L_3+1=7`, `L_5+1=61`, `L_7+1=421` is prime, but the
  fourth term breaks it: `L_11+1=27721=19*1459`.  THM-1226 uses `q=27721`
  only because `q=A+1`, never because it is prime;
- theorem IDs, equation labels, commit fragments, and target-specific ruler
  bounds.

## VII. What the pattern now says, and what it does not

The useful common object is not the four integers alone.  It is the cutoff
operator

\[
X\longmapsto L_X=\operatorname{lcm}(1,\ldots,X),       \tag{20}
\]

viewed through different quotients:

1. an AP owner clock retains selected prime powers of `L_(P-2)`;
2. an LRC spectrum or covering system retains the whole cutoff LCM;
3. a blocker speed maps to zero in every quotient `Z/qZ` below the cutoff;
4. a one-far endpoint ledger multiplies the base LCM by the sector count;
5. a CRT or root-of-unity observer joins local clocks but forgets lift data.

That is the genuine synthesis.  The preserved predicate is simultaneous
periodicity/divisibility in a finite family of cyclic quotients.  The destroyed
information depends on the application: owner labels, pole positions, signs,
gcd packets, exact lifts, or boundary jets.  Those destroyed coordinates are
precisely why the pattern does not by itself advance LRC(14), construct a
planar-Jacobian counterexample, or classify Sun holes.

The strongest new exact result of the search is the master identity `(4)` and
its equality criterion `(5)`.  The strongest newly recovered full carrier is
THM-879's divisor-mode quotient.  The strongest cross-frontier coordinate is
`M_P/R_P`, whose valuation depth simultaneously explains the THM-4042 collapse
and the THM-3848 radical hostile; the `3960 -> 27720` square `(15)` is the
strongest phase-clock coincidence.  The sharp stopping certificate is `P=17`:
the next prime-sector row destroys both the previous-prime rule and the nested
clock chain.  Future work should therefore search for maps preserving a named
sidecar, not for further bare occurrences of the same natural numbers.
