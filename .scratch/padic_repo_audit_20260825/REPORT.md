# Version-pinned audit: `octonion/p-adic-zeta-irrationality`

**Audit date:** 2026-08-25
**Scope:** primary-source repository, manuscript, and numerical certificate;
typed comparison with THM-4056/4057/4059/4061/4068/4071/4076 and existing
irrationality work.
**Canon edits:** none.

## Outcome first

The repository contains a serious, self-contained **research-draft claim** of
22 individual Kubota--Leopoldt irrationalities.  It is not merely a numerical
announcement: the 2,317-line manuscript supplies a proposed proof and names
its imported interfaces.  However, it is a one-day-old, untagged, unreviewed
repository, and its own README asks for independent specialist review.  In
the math repo's status language the correct current classification is:

- **AUTHOR-CLAIMED / UNREFEREED:** the hybrid arithmetic-holonomy theorem and
  the resulting 22 irrationality statements;
- **FINITE-EXACT / independently replayed here:** the 22 final real margin
  inequalities, conditional on the manuscript's formulas;
- **CITED:** Calegari's Eisenstein family, Buzzard continuation, Katz Hasse
  invariants, GPRJ's de Rham torsor, and the Bost--Charles/CDT slope, source,
  zero-estimate, and prime-window interfaces;
- **not claimed:** transcendence, algebraic independence, an irrationality
  measure for the new cells, or a rational/irrational/transcendental
  classification by tournaments.

The new draft is therefore real incoming signal, but it should not overwrite
the older `OPEN` singleton ledger as `PROVED` until the geometric and adelic
interfaces have been independently audited.  The accurate update is:
"previously open in the audited published literature; now claimed by the
unrefereed Long draft at the pinned commit below."

## Immutable source pin and artifact inventory

- Repository: <https://github.com/octonion/p-adic-zeta-irrationality>
- Audited commit:
  [`b46a1770901551961710e155d775aae7c5ea39e7`](https://github.com/octonion/p-adic-zeta-irrationality/commit/b46a1770901551961710e155d775aae7c5ea39e7)
- Default branch: `main`; no tags or releases.
- Git history at audit time: three commits: MIT license, initial source import,
  then compiled PDF addition.  The last two commits carry the same generic
  "initial import" subject; there is no public repair history for the source.
- GitHub metadata at audit time: created 2026-08-25 00:40:36Z, last pushed
  2026-08-25 01:05:55Z, zero stars/forks/issues, not archived.
- Contents: README, 2,317-line TeX manuscript, compiled PDF, one 373-line
  standard-library Python certificate, frozen output, and SHA-256 manifest.
- No CI, test suite, formalization, dependency lock, issue/PR audit, or code
  checking the symbolic geometric arguments is present.

The manuscript labels itself a "Re-audited product-digit repaired research
draft" dated 24 August 2026, but that repair lineage is not represented in
the public Git history.

## Exact claimed theorem and mechanism

The draft claims irrationality of

```text
zeta_2(s), s=3,5,...,29                         (14 cells)
zeta_3(s), s=3,5,...,11                         (5 cells)
zeta_5(3), zeta_5(5), zeta_7(3)                (3 cells)
```

The proof architecture is a high-dimensional Hermite--Padé/Apéry analogue,
not an explicit scalar recurrence.

1. Under the contrary hypothesis `eta=zeta_p(s)/2 in Q`, form

   ```text
   K(q)=eta+sum_(n>=1) sigma^p_(-s)(n) q^n.
   ```

   Calegari's cited result supplies the overconvergent Eisenstein form.

2. Define raw Eichler rows

   ```text
   R_e=D^(s-e)J,  1<=e<=s.
   ```

   The manuscript proves

   ```text
   n^e [q^n]R_e=sum_(d|n,p not|d)(n/d)^s in Z,
   L_N^e R_e mod q^(N+1) integral,
   R_e(q)-ell^(-e)R_e(q^ell) in Z_ell[[q]]  (ell != p).
   ```

3. A canonical BGG lift places the rows in a rank `s+1` algebraic
   Eisenstein--Eichler connection.  Functional monodromy rank gives an
   injective degree-`D` global source and `r=(s+1)D+O(1)` distinct Bost pivot
   orders.

4. At small auxiliary primes, the Hasse invariant makes the normalized
   unipotent coordinate rational modulo `ell` with degree
   `(p+1)ell/12+O(1)`.  Rank--nullity on the **genuine positive global source**
   produces a large kernel on which the unique depth-`s` multiplier is
   divisible by `ell`, saving exactly one `ell` in the determinant lattice.

5. At large auxiliary primes, the calculation stays on the de Rham frame
   torsor.  The raw Dwork relation, a residue-class Taylor no-backflow lemma,
   digitization of the full products, and descending pole-grade induction
   recover the CDT multiplier window without the invalid global scalarization
   explicitly rejected by the draft.

6. Bost's slopes inequality compares the modified source degree with pivot
   heights.  Buzzard continuation gives the finite-place radii; an internal
   modular Jensen formula gives the archimedean collision energy.  Rationality
   would force

   ```text
   M_(p,s)(xi,Y)=s Lambda_p(Y)-(s+1)tau_(p,s)(xi)-C_p(Y)<=0.
   ```

7. The certificate encloses the same quantity strictly above zero in 22
   cases.  The smallest lower endpoint is the `(p,s)=(5,5)` cell,
   `0.131799356827...`.

There is no Apéry recurrence in `n`, no Padé determinant/Casoratian printed as
a finite sequence, and no transcendence criterion.  The recurrence-like
object is the exact Frobenius relation in step 2; the Padé object is the
global degree-`D` source paired with a holonomic horizontal leaf.

## Reproduction audit

At the pinned commit:

```text
normal output lines:    35
optimized output lines: 35
frozen output lines:    35
normal diff:             0
optimized diff:          0
```

All 22 normal-mode positivity assertions passed.  Raw-LF normalization gives
the manifest hashes exactly:

```text
1408c9092d0c34917253c8f520db56853112c58640d15c98d58b76a61a0478f5  script
9e1d8f74eb198c7e7b0c7420a1e8021d7bad4d2bc8014cc939b6b353111c0da0  output
```

On the current Windows checkout Git converted the files to CRLF, so hashes of
working-tree bytes differ; converting CRLF back to LF reproduces the manifest.
Thus the manifest is sound but its README command is not portable to a default
`core.autocrlf` checkout without stating the raw-LF basis.  Also, positivity
and domain checks are Python `assert` statements, so normal mode is the actual
certificate mode; `-O` reproduces the transcript but suppresses those gates.

The program rigorously checks the **last numerical implication only**.  It
does not check BGG gluing, source descent, Hasse degree, torsor Cartier
strictness, Bost/CDT applicability, continuation radii, or modular-Jensen zero
multiplicity.

## A new exact audit: which margins genuinely use the Hasse saving?

Keeping every witness and every other formula fixed, delete only the Hasse
rank--nullity saving.  Since

```text
tau_base-tau_hyb = 2(s-1)J_p(xi)/(s+1)^2,
```

the margin loses exactly

```text
Delta M = 2(s-1)J_p(xi)/(s+1).
```

An exact `Fraction` replay plus the certified real interval shows that only
three of the displayed witnesses change sign:

| cell | hybrid lower margin | formal no-Hasse lower margin |
|---|---:|---:|
| `(2,29)` | `1.5222476232` | `-0.1456124267` |
| `(5,5)` | `0.1317993568` | `-0.9157543697` |
| `(7,3)` | `0.4364130578` | `-0.2503715416` |

The other 19 displayed witnesses remain numerically positive after this one
formal deletion.  This does **not** say those 19 were previously proved: the
torsor repair and all other proof gates remain load-bearing.  It isolates the
three cells for which the new one-power Hasse determinant saving is
numerically indispensable in the manuscript's own packet.

A separate `FINITE-EXACT + FINITE-INTERVAL` companion now proves the complete
global optimization of the displayed formula and certifies the next adjacent
cells still negative:

```text
(2,31): global margin <= -1.943953720741976
(3,13): global margin <= -1.655957196706988
(5,7):  global margin <= -3.841753001089363
(7,5):  global margin <= -3.609764970278347
```

The proof audits every floor and positive-part boundary of the prime-window
function: its derivative joins continuously and is nondecreasing, while the
negative Hasse integral is convex.  Hence the paper's rational `xi_witness`
is the unique global minimizer of `tau`.  In the `Y` variable each collision
layer enters with both value and derivative zero; the derivative of the
margin is continuous and strictly decreasing, so a rational derivative-sign
bracket and one concave tangent bound the global maximum.  The verifier and
frozen transcript are
`next_case_margin_audit.py` and `next_case_margin_audit.out` beside this
report.  Normal and optimized transcripts both match the frozen output.
The companion forces LF on stdout, so this is a raw-byte comparison even on
Windows (not merely a line-normalized comparison):

```text
6a4063e4a041fe30f419d34a594ea8e56b32fcf7315bd4137ca4f0eae71c4a70  next_case_margin_audit.py
cbe529a95384e41943adca5c45e546dc8484431e3e78349b18b151ece35c3f54  next_case_margin_audit.out
```

This is an obstruction for the pinned hybrid-margin formula, not a general
no-go against a different continuation template or irrationality method.

## Delicate proof gates needing specialist audit

No error is asserted here, but these are the first interfaces whose failure
would invalidate the headline theorem while leaving the numerical transcript
unchanged.

1. **BGG and genuine multiplier:** verify the exact negative-weight BGG
   normalization, the claimed single unipotent action, and that `H^{-r}`
   introduces only fixed, uniformly spread-out poles.
2. **Global source descent:** verify descent from a fixed neat cover to the
   coarse `X_0(p)`, uniform saturation/base change in `D`, and the claim that
   all frame/model changes cost only `o(D^2)`.
3. **Hasse degree and kernel:** audit the orbifold-to-coarse degree
   `(p+1)ell/12+O(1)`, integrality of the actual multiplier, and the passage
   from a mod-`ell` kernel to the adelic source lattice.
4. **Fixed-bundle Bost/CDT transfer:** check that the cited scalar source,
   pivot-height, and prime-window formulae tolerate the vector-bundle source
   with exactly the claimed leading coefficients.
5. **Torsor Cartier pole grades:** independently reconstruct the product-digit
   expansion, Taylor associated grades, and no-backflow induction over the
   special-fibre coefficient algebra.
6. **Continuation:** check Buzzard's theorem in negative integral weight and
   the exact `f_p` radii, especially the internal `p=5,7` Newton-polygon
   identifications.
7. **Modular Jensen:** verify the primitive-bottom-row enumeration, orbifold
   multiplicities, and equality rather than a one-sided Jensen bound.

The `(5,5)` margin is only `0.1318`, so constant/sign/normalization errors in
these interfaces cannot be dismissed as numerically harmless.

## Reusable objects and exact connection ledger

### 1. THM-4056 -- the connection is literal

| field | entry |
|---|---|
| source | external raw row coefficient `[q^n]R_e` |
| target | THM-4056 `L_N` valuation-depth clock |
| map | `n -> (v_ell(n))_ell`, row label `e ->` clearing depth `e` |
| preserved | `n^e[q^n]R_e in Z`; `L_N^e` clears every index `n<=N`; new clock coordinates occur at prime powers |
| lost | numerator cancellation, the distinguished-prime deletion `p not|d`, exact minimal denominator, and analytic decay |
| required sidecar | numerator valuations and the exact Dwork residual `R_e(q)-ell^(-e)R_e(q^ell)` |
| cheapest decisive test | tabulate minimal coefficient denominators and Dwork residual valuations for all `e<=s,n<=N`; compare to the worst `e v_ell(n)` clock |

This is the strongest reusable bridge.  It also supplies a new operation for
the repo's Apéry searches: seek a mod-`ell` algebraization of the multiplier
of the deepest recurrence row and use genuine-source rank--nullity to save one
`ell` on a large subspace.  A scalar gcd observed after the fact is not enough.

### 2. THM-4057 -- two lawful rational-edge carriers

**Modular-Jensen carrier.**  The manuscript indexes nonidentity Hauptmodul
collisions by primitive bottom rows `(c,d)` of `Gamma_0(p)`, with
`gcd(c,d)=1`, `c>0`, and `p|c`.  Taking `d mod c` maps these to the primitive
rational arcs `d/c` in THM-4057.

| preserved | coprimality, denominator shell `c`, and the `phi(c)` number of primitive owners |
| lost | the full matrix/coset, translate of `d`, and analytic height `Im(gamma tau)` |
| sidecar | the modular transform or collision height; denominator orientation `p|c` |
| hostile | reciprocal reversal generally sends the selected denominator to the numerator and destroys `p|c` and `cY<1` |
| test | enumerate bottom rows through a finite `c` cutoff and compare shell counts with the natural-edge fibre dictionary |

These rows form a coprimality graph/subset, not a tournament.  Completing all
missing pairs by a gauge would add information not used by Jensen.

**Pivot carrier.**  Let `V_D={n_1<...<n_r}` be the distinct Bost pivot
orders.  Orient every pair from smaller to larger.  This is a canonical
transitive tournament and

```text
sum_i n_i >= 0+1+...+(r-1)=binom(r,2)=# tournament edges.
```

This is exactly the combinatorial lower bound used at manuscript lines
1850--1855.  It preserves rank and filtration order and loses every local
height, coefficient, and arithmetic denominator.  The required sidecar is
the weighted adelic evaluation height.  Consequently the tournament is a
valid proof skeleton but cannot classify the value as rational, algebraic
irrational, or transcendental by itself.

### 3. THM-4059 -- collision shells are its unsigned quotient

At a fixed denominator `c`, THM-4059 retains every unit `d`, its modular
inverse, and Stern sign

```text
epsilon_c(d)=(-1)^D(d/c);  sum epsilon_c(d)=S(c).
```

For odd `c` this is `(-1)^(d+d^(-1))`; for even `c` it is
`-(-1)^((d d^(-1)-1)/c)`.  Keeping this parity split matters because the
`p=2` Jensen shells are even.

The modular Jensen formula unfolds that shell to the scalar multiplicity
`phi(c)`.  Thus the source-to-target map is

```text
labelled inverse-depth unit packet -> cardinality of the unit packet.
```

It preserves exact denominator and owner count and destroys inverse,
Stern depth, and sign.  A Stern-signed Jensen diagnostic is well-defined,
but it cannot replace the positive collision energy: signed cancellation is
not an upper or lower bound for the unsigned zero mass.  The missing sidecar
is total variation/owner-resolved collision height.

Cheapest test: for small active shells, retain each residue `d mod c` before
unfolding, compute both signed and unsigned Jensen contributions, and exhibit
the first cancellation witness.  This tests structure, not a proof saving.

### 4. THM-4061 -- a possible signed-Jensen Fourier laboratory

Both sources use the same primitive unit shell.  THM-4061 retains the modular
inverse and turns the Stern sign into a complete Kloosterman/Fourier sum;
Jensen integrates away the owner and retains only `phi(c)`.  A lawful
connection is therefore:

| source | owner-resolved modular collision shell |
| target | THM-4061 inverse-parity Fourier packet |
| map | decorate each primitive bottom row by `epsilon_c(d)` before summing |
| preserved | denominator, coprimality, and owner label |
| destroyed | positivity after taking the signed sum |
| sidecar | absolute collision measure and Fourier transform of the Jensen kernel |
| test | derive the exact finite Fourier formula for one shell and compare `c=3,5,9,15` hostiles |

This may reveal cancellation patterns but cannot reduce Bost energy without a
separate domination theorem.

### 5. THM-4068 -- local factorization versus global owner data

The common mechanism is not a shared theorem but a precise meta-pattern.
THM-4068 factors the complete CRT kernel while canonical representative
parity refuses to factor.  The draft constructs different Hasse kernels
`K_(D,ell)` at different primes; their local determinant indices add adelically,
but no single global scalar coordinate or common kernel is asserted.

| source | family of prime-local Hasse kernels / local Kloosterman kernels |
| target | adelic source lattice / global Stern packet |
| preserved | local integrality or complete additive character |
| destroyed | common basis, canonical representative, and owner synchronization |
| sidecar | adelic norm/model in the draft; dyadic CRT carry in THM-4068 |
| test | compute actual small-`D` Hasse matrices at two auxiliary primes and compare the individual kernel dimensions with their intersection |

The hostile is immediate: summing local codimensions does not prove a common
rational subspace of that rank.  The draft avoids this by modifying the
adelic lattice prime by prime; any attempted tournament selector must do the
same or retain a common-owner sidecar.

### 6. THM-4071 -- two different p-adic stationary mechanisms

THM-4071 stratifies a prime-power Kloosterman sum by
`v_p(h),v_p(k)` and uses the first surviving top digit plus a two-critical-
class p-adic Morse coordinate.  The draft stratifies a coefficient by its
negative `ell`-adic pole grade, digitizes the full truncated product, and uses
Taylor no-backflow.

| map | valuation filtration -> first nonzero associated-graded digit |
| preserved | prime-power depth, residual unit/nonunit status, and exact vanishing before the first live grade |
| lost | additive phase in THM-4071 versus coefficient multiplier/row ownership in the draft |
| sidecar | THM-4071 affine carry `kappa`; draft full product digits `B_(e,nu)` and fibre coefficients |
| hostile | valuation alone does not contract a lift fibre (THM-4071), and digitizing `Q_e` without `R_e` loses positive valuations (draft) |
| test | generic finite-ring replay of no-backflow alongside the THM-4071 affine lift recursion; perturb the omitted carry/product digits to force failure |

The reusable lesson is exact: **digitize the full product and retain the carry**.
Neither stationary-phase estimate transfers numerically to the other object.

### 7. THM-4076 -- exponent counts are only a quotient of the valuation clock

THM-4076's apex congruence retains `Omega(q)=sum_p v_p(q)`.  The draft and
THM-4056 need the full weighted valuation vector

```text
(v_ell(L_N) log ell)_ell.
```

The map from the vector to `Omega` preserves the total number of prime-power
layers and destroys their primes, logarithmic costs, and cutoff locations.
The minimal hostile is any pair with equal `Omega` but different valuation
vectors, e.g. `9` and `15`.  Therefore the mod-four apex degree can diagnose
the number of valuation events but cannot price an irrationality lattice.

The more abstract resemblance--one inversion fixed-point bit versus one
Hasse `ell`-power saving--is only a meta-pattern (a fixed-locus/codimension
sidecar recovers one discrete unit), not a mathematical map.

## Tournament status and arithmetic classification firewall

Three tournament-like objects should not be conflated.

1. **Rational-edge orientation:** numerator to denominator on primitive
   integer pairs.  It organizes rational approximants and modular bottom rows.
2. **Stern depth tournament:** a finite parity/sign observable on rational
   packets.  It controls shell imbalance, not approximation quality.
3. **Pivot-order tournament:** the transitive total order of nonzero Padé
   pivots.  It records the rank/vanishing filtration term in Bost's inequality.

None classifies rational/irrational/transcendental status.  For any fixed
`p`-adic precision `N`, a prescribed digit prefix is shared by a rational, by
an algebraic irrational in `Q_p` (after a suitable Hensel lift), and by
uncountably many transcendental `p`-adic numbers.  Any finite tournament built
only from those digits, congruences, rational edges, or valuation comparisons
is identical on all three controls.  The separating data in the draft are
global and asymptotic: functional independence, analytic continuation width,
adelic denominator height, and nonzero pivot evaluations.

Even an infinite pattern needs a theorem linking that pattern to height and
nonvanishing.  An orientation count or balanced apex is neither necessary nor
sufficient for irrationality.

## Relation to earlier repo irrationality work

- The Apéry framework's four obligations--denominator clearing, decay,
  nonvanishing, and coefficient growth--reappear here as hybrid adelic lattice,
  analytic width/energy, pivot injectivity/zero estimate, and source degree.
  The draft is a vector-bundle Hermite--Padé realization of the same proof
  architecture.
- The existing `zeta_2(5)` energy scouts already warned that derivative and
  energy must come from the same template and used modular companion preimages
  to audit Jensen mass.  The new draft gives an exact general modular Jensen
  formula; that proposition should be checked against those rigorous
  `Gamma_0(2)` controls.
- The 2026-08-21 singleton ledger listed `zeta_5(3), zeta_7(3), zeta_2(7),
  zeta_3(7), zeta_2(9)` as open in then-audited literature.  All are now among
  the draft's claims.  This is a status-changing claim, not yet a status-changing
  proved theorem.
- The old `zeta_13(3)` frontier is untouched: the draft treats only
  `p in {2,3,5,7}`.  Its Hasse kernel suggests a denominator-thinning
  experiment at 13, but no continuation/energy gate is supplied there.
- The manuscript explicitly repairs the same type of scalarization error that
  repo guardrails emphasize: raw derivative depths do not become globally
  defined scalar coordinates after a nonzero-weight frame change.

## Highest-value next computations

1. **Independent symbolic local audit (highest priority).**  For the smallest
   nonclassical cells `(5,5)` and `(7,3)`, construct the BGG rows and genuine
   multiplier in a computer algebra system, reduce at several auxiliary
   primes, and verify the Hasse degree and kernel rank on the actual global
   source.  Positive and hostile controls: `(2,3)` and one deliberately
   doubled-unipotent scalarization.
2. **Independent pole-grade engine.**  Reimplement Proposition 4.11 from the
   displayed algebra only, using random flat coefficient rings including
   zero-divisor special fibres; verify no-backflow and exhibit failure if the
   full product digit or carry is deleted.
3. **Jensen reconciliation.**  Specialize the paper's formula to the prior
   `Gamma_0(2)` lune/circle controls and compare exact active cosets,
   normalization, and collision multiplicities.
4. **Beyond the now-closed adjacent cells.**  The pinned formula's immediate
   `(2,31),(3,13),(5,7),(7,5)` extension is rigorously negative.  Search must
   now change continuation energy, Hasse codimension, or large-prime cost,
   rather than tune `xi,Y` inside the same packet.
5. **LCM minimal-denominator atlas.**  Measure actual cancellations in
   `[q^n]R_e` relative to `n^e`, grouped by prime-power clock events.  Feed the
   resulting full valuation vectors, not gcd scalars, into the repo's Apéry
   recurrence searches.
6. **Signed Jensen laboratory (diagnostic only).**  Retain bottom-row owners,
   attach THM-4059 Stern signs, and compute Fourier-resolved signed shells.
   Require an unsigned total-variation sidecar before any energy claim.
7. **Pivot tournament audit.**  Extract finite truncated pivot orders and local
   heights.  Verify that the transitive edge count supplies only the
   `binom(r,2)` term while every arithmetic distinction lives in edge/vertex
   weights omitted by the bare tournament.

## Primary links

- Pinned repository README:
  <https://github.com/octonion/p-adic-zeta-irrationality/blob/b46a1770901551961710e155d775aae7c5ea39e7/README.md>
- Pinned manuscript source:
  <https://github.com/octonion/p-adic-zeta-irrationality/blob/b46a1770901551961710e155d775aae7c5ea39e7/hybrid_arithmetic_holonomy.tex>
- Pinned certificate:
  <https://github.com/octonion/p-adic-zeta-irrationality/blob/b46a1770901551961710e155d775aae7c5ea39e7/hybrid_rational_interval_certificate.py>
- Calegari, *Irrationality of certain p-adic periods for small p*:
  <https://arxiv.org/abs/math/0408214>
- Calegari--Dimitrov--Tang, *The linear independence of `1`, `zeta(2)`, and
  `L(2,chi_-3)`*: <https://arxiv.org/abs/2408.15403>
- Calegari--Dimitrov--Tang, *Arithmetic holonomy bounds and effective
  Diophantine approximation*: <https://arxiv.org/abs/2510.04156>
- Graham--Pilloni--Rodrigues Jacinto, published DOI:
  <https://doi.org/10.1112/S0010437X25102479>
