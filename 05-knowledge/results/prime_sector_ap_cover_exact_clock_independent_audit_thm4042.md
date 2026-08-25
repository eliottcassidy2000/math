# Independent hostile audit: THM-4042 prime-sector exact clock

Audited files:

- `01-canon/theorems/THM-4042-prime-sector-ap-cover-exact-clock-and-holonomic-law.md`
- `04-computation/prime_sector_ap_cover_exact_clock_thm4042.py`
- `04-computation/prime_sector_ap_cover_exact_clock_independent_audit_thm4042.py`

## Verdict

**PASS.** The owner classification, eventual local max--min law, closed
winners, cyclic word, exact phase period, sharp common desingularizer,
recurrence, and D-finite/nonalgebraic OGF conclusions are correct. I found no
counterexample to a theorem quantifier.

THM-4042 imports THM-4033's separately proved sharp odd-prime geometric onset
`(P^2+3)/4`. This report checks the imported interface and the downstream
selector/clock theorem; it does not independently re-audit THM-4033's
three-gap/guard proof.

Three precision points found during audit are repaired in the canonical file:

1. local displacements and both `rho` values are explicitly typed in the
   `theta=P x` coordinate, explaining the later `1/P` normalization;
2. least cyclic period is formally a least positive translation period;
3. the top-word formula states the convention `rad(1)=1`.

An older exploratory helper with an invalid pre-onset fallback was not
promoted and is outside this audit. The official target and independent
scripts both reject unavailable terminal tracks instead of repairing them.

## 1. Persistent-owner completeness and eventual onset

Write `theta=P x`.  At a reduced rational `x=a/q`, the orbit in the
`theta`-circle is the uniform grid

\[
\{P k/q:0\le k<q\}.
\]

- If `q<P`, its `q` occupied sectors cannot cover `P` sectors, so the point is
  a persistent noncover owner.
- If `q>P`, the mesh `P/q<1`; every open unit sector contains a grid point,
  giving a finite covering prefix stable under perturbation.
- If `q=P`, primality makes every nonzero numerator a unit.  Times
  `0,...,P-1` cover for a small positive perturbation and times `1,...,P`
  cover for a small negative perturbation.  Their union gives one finite
  prefix stable on a two-sided neighborhood.
- For irrational `x`, density supplies an interior hit in every sector and
  hence a stable finite prefix.

Thus the persistent owner set is exactly

\[
B_P=\{a/q:1\le q<P,\ 0\le a<q,\ (a,q)=1\}.
\]

After choosing disjoint owner neighborhoods, the complement is compact and
has a common stable covering prefix.  Shrinking the owner neighborhoods
leaves a compact annulus, which is covered by another common prefix.  Within
the final neighborhoods, the no-skip condition `q|delta|<1`, together with
retention of each base-sector representative, makes the terminal time

\[
E_s(n)=n-((n-s)\bmod q)
\]

decisive on every track.  This proves the exact local formula for every
sufficiently large `m`.  There are finitely many owners, phases, missing
sectors, and competing affine fractions; strict inequalities between their
distinct leading coefficients therefore stabilize uniformly after another
finite enlargement.  The claimed existence of a single `M_P` is justified.

What is **not** proved is a closed form for the least such `M_P`.  The exact
finite data currently have the following scope:

| `P` | tested window | first equality continuing to window end | status |
|---:|---:|---:|---|
| 3 | `3..16` | 3 | finite suffix |
| 5 | `5..24` | 7 | finite suffix |
| 7 | `7..32` | 13 | finite suffix; the all-tail claim is separately available from the seven-sector theorem |
| 11 | `11..50` | 31 | finite suffix |
| 13 | `13..45` | 43 | finite suffix |

The independently recomputed `P=13` bit string is thirty zeros followed by
three ones, so the stated `43` is the least equality in that window and
agrees through `45`.  It is not by itself a persistent-tail certificate.

## 2. Local max--min law and closed winners

For `q<P` and prime `P`, multiplication by `P a` is invertible modulo `q`.
Consequently the nonzero fractional tracks never lie on an integer wall, and
the occupied sector labels are distinct.  In the `theta` coordinate, moving
a track positively to a missing sector costs

\[
d^+_{s,r}-f_s,
\]

while moving negatively costs

\[
d^-_{s,r}-1+f_s.
\]

Monotonicity and no-skip propagation show that dividing by the last available
time `E_s(n)`, then taking `min_s` and `max_r`, is both necessary and
sufficient.  The possible strict negative endpoint changes no measure.

Sorting the base points as `y_k=P k/q`, the last positive-side missing sector
in a gap has cost

\[
P/q-1-\{y_{k+1}\}.
\]

Its unique maximum is before `y_q=P`, giving coefficient `(P-q)/q` and
track `a s=-1 mod q`.  The negative-side first missing sector has cost

\[
\{y_k\}+P/q-2.
\]

Its unique maximum has `P k=-1 mod q`; the following track satisfies
`a s=1-P^{-1} mod q`.  This yields exactly

\[
s^+=-a^{-1},\qquad
s^-=(1-P^{-1})a^{-1},
\]

with coefficients `(P-q)/q` and `(P-q-1)/q`.

For `q=P-1`, the negative coefficient is zero and the radius really is
identically zero: the sole missing sector is filled immediately on the
negative side while time zero retains sector zero.  The `q=1` calculation
separately gives `P-1` and `P-2`.

The independent checker exhaustively rebuilt the nested leading selectors
for all reduced owners at primes through `19`: 546 side/owner checks passed,
including uniqueness whenever the winning coefficient is positive.

## 3. Word convolution and normalization

On phase `r=n mod q`, a winning track `s` contributes at pole shift

\[
c=(r-s)\bmod q.
\]

As `a` runs over the units modulo `q`, the plus tracks run through `-U_q`, so
`c-r` runs through `U_q`.  For the minus side,

\[
c-r=(P^{-1}-1)a^{-1}.
\]

After multiplying each theta-radius by `1/P`, this is exactly the two-term
word in equation (3).  A denominator-`q` block has shifts only
`0,...,q-1`; hence pole `c` receives exactly the blocks `q>=c+1`, proving the
convolution in (4) and rational phase law (5).

An independent construction from the raw max--min selectors, without
importing either audited file, matched the closed word for all 35 denominator
blocks at `P<=13`.

## 4. Exact phase minimality

For complete precision, define

\[
\delta_{P,q}=\min\{d\ge1:
w_{P,q}(c+d)=w_{P,q}(c)\ \text{for every }c\in\mathbb Z/q\mathbb Z\}.
\]

This minimum divides `q`.  The coefficient-vector period is then
`Pi_P=lcm_q delta_(P,q)`.

The descending-pole proof is valid.  At coordinate `c=P-2`, only the
`q=P-1` block occurs, so every period of the aggregate vector must be a
period of that block.  Subtract it and descend: at coordinate `c=q-1`, the
remaining top block is exactly denominator `q`.  Thus every `delta_(P,q)`
divides any aggregate period.  The lcm is both a period and the least one.
Uniqueness of partial fractions transfers minimality to the rational phase
functions.  Multiplication by the nonzero common polynomial transfers it to
the phase polynomials.

For `q=P-1`, the word is the scaled unit indicator.  Unit status modulo
`P-1` depends exactly on residue modulo `rad(P-1)`.  Conversely, any
translation period missing a prime divisor can be defeated by CRT with a
unit translated to a multiple of that prime.  Therefore the top period is
exactly `rad(P-1)`.

There is also a closed formula for every lower word. Write
`t=P^(-1)-1=(1-P)/P mod q`. At a prime power `ell^e || q`, with
`h=min(v_ell(P-1),e)`, the local image `tU_q` is the exact-valuation set
`v_ell(c)=h` if `h<e`, and the zero class if `h=e`. Its least additive period
is therefore `ell^min(e,h+1)`. CRT multiplies these periods.

The fibre multiplicity is constant on `tU_q`. If `t` is a unit, both terms
of the word are supported on `U_q`, giving period `rad(q)` and the same
formula. If it is not a unit, `U_q` and `tU_q` are disjoint and their two
positive word values cannot coincide: with `s=P-q>=2` and fibre size `K>=2`,
coincidence would say `s=(s-1)K`, forcing `(s,K)=(2,2)`, but `q=P-2` makes
`t` a unit. Hence the level sets cannot swap, and for `2<=q<=P-2`

\[
\delta_{P,q}=\prod_{\ell^e\parallel q}
\ell^{\min(e,v_\ell(P-1)+1)}.
\]

Maximizing prime exponents over `q<=P-2` gives the closed product formula
for `Pi_P` in THM-4042. The independent checker compared this formula with
the directly computed least cyclic period for all `6,027` pairs `(P,q)` with
prime `P<=251`.

The independent profile reproduced

```text
Pi_2=1, Pi_3=2, Pi_5=6, Pi_7=60,
Pi_11=420, Pi_13=27720, Pi_17=120120.
```

The target probe's `minimal_period` tests only divisors of the list length;
this is complete for a cyclic word because every translation symmetry has a
gcd-with-length symmetry and the least positive period divides the length.

## 5. Leading constant and sharp common desingularizer

Summing a word is phase invariant.  The plus and minus masses for denominator
`q` are respectively

\[
\frac{\phi(q)(P-q)}{Pq},\qquad
\frac{\phi(q)(P-q-1)}{Pq},
\]

including the stated `q=1` specialization.  This gives `kappa_P` exactly,
and positivity makes the degree of `Q_(P-1)D_P` exactly `P-2` on every phase.

Every pole `c=0,...,P-2` is active on some phase already in the translate
family of the top word.  All residues at a fixed pole are nonnegative, so an
active pole cannot cancel.  A polynomial clearing every phase must therefore
vanish at every `c` and be divisible by

\[
Q_{P-1}(n)=\prod_{c=0}^{P-2}(n-c).
\]

Since this polynomial works, it is the unique monic common desingularizer of
least degree.  This proves both sharpness and exact phase preservation after
clearing.

## 6. Recurrence and OGF

On each `Pi_P` phase, `Y_P=Q_(P-1)D_P` is a polynomial of degree `P-2`.
The `(P-1)`st difference with step `Pi_P` therefore gives equation (10) for
all sufficiently large `n`.  This is a polynomial-coefficient recurrence
(unused intermediate shifts simply have zero coefficient), so the whole
sequence is P-recursive after adjoining its finite prefix.  Its OGF is
D-finite.

From `D_P(n+1)=kappa_P/n+O(n^-2)` with `kappa_P>0`, the OGF differs from
`kappa_P log(1/(1-z))` by a function bounded as `z` tends to one from below.
An algebraic branch has a Newton--Puiseux expansion and cannot have unbounded
logarithmic growth.  The nonalgebraicity conclusion is valid.  For cleaner
indexing, THM-4042 states the OGF explicitly as an eventual one; a finite
prefix changes neither conclusion.

The independent checker also verified the displayed recurrence exactly from
the closed controller for the sample primes through `17`.

## 7. Prime/composite boundary

Primality is used through `gcd(P,q)=1` for every proper owner denominator.
It supplies both the nonvanishing fractional tracks needed to retain occupied
sectors and the inverse `P^{-1} mod q` in the closed winner.  For composite
sector count, denominators sharing a factor with the sector count introduce
integer-wall tracks.  Negative drift can then delete a formerly occupied
sector, information absent from the prime max--min formula.

The stated hostile is exact:

\[
D_4(20)=299/2907,
\qquad D^{\rm naive}_4(20)=2069/23256,
\qquad D_4(20)-D^{\rm naive}_4(20)=1/72.
\]

This refutes the naive prime formula for composite `4`; it does not say that
no corrected composite theorem exists.  A boundary-orientation and
vanishing-track sidecar is the right missing datum.

## 8. Code scope and reproduction

The target companion uses exact `Fraction` arithmetic throughout. Its
lexicographic `(C,c)` order is the correct order at infinity because

\[
C/(n-c)=C/n+Cc/n^2+O(n^{-3}).
\]

Zero candidates are identically zero and are safely normalized to shift
zero. The target's global-row calculation is limited to `P<=13`, while the
profile through `61` is computed from the proved closed word; the theorem
describes that distinction correctly.

Independent replay:

```powershell
python -B 04-computation/prime_sector_ap_cover_exact_clock_independent_audit_thm4042.py
python -B -O 04-computation/prime_sector_ap_cover_exact_clock_independent_audit_thm4042.py
```

Both outputs are byte-identical and end with:

```text
winner_formula_checks=546
direct_word_block_checks=35
p_adic_period_checks=6027
period_profile={2: 1, 3: 2, 5: 6, 7: 60, 11: 420, 13: 27720, 17: 120120}
finite_prime_windows=P3,P5,P7
composite_P4_m20_gap=1/72
RESULT=PASS
```
