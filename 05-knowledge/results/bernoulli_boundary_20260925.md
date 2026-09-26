# The exceptional B1: boundary jumps, Collatz clocks, and power-sum problems

**Status:** CITED classical Bernoulli identities and existing number-theory
results; PROVED elementary bridges and scoped obstructions, independently
audited; FINITE-EXACT controls with explicit universes. Collatz, the
Erdos--Moser equation conjecture, and Agoh--Giuga remain OPEN. No novelty
claim for the classical identities or their immediate consequences.

The user's seed is that the initial Bernoulli numbers form a microcosm,
while the later odd zeros suggest a symmetric macrocosm. A precise useful
interpretation is **endpoint defect versus symmetric bulk**. The index of
a Bernoulli number is not itself a time or scale parameter; any transfer
must specify the summation, normalization, and retained endpoint data.

## Inheritance and concept board

Anchor: relate Bernoulli boundary terms to actual Collatz halving counts.
Niche: the normalized power-sum equation, where the same endpoint correction
survives in a rigorous asymptotic. Wildcard: the literal exceptional-even
case in Agoh--Giuga.

Closest proved mechanism: the exact affine-block defect and countdown in
[the preceding Kuratowski reframe](kuratowski_reframe_20260925.md), especially
its sections 3--5. Canonical hostiles are the minus cycle `5->7->5`, the
`5n+1` cycle through 13, and arbitrary initial rises. Corrected near miss:
ordinary mean drift does not control every integer orbit. Least-used
sidecar: the value of the periodic Bernoulli function at its discontinuity.

| Live concept | What it retains | Decisive boundary |
|---|---|---|
| Bernoulli reflection | endpoint orientation | odd polynomials need not vanish away from endpoints |
| Periodic B1 jumps | integer divisibility | midpoint convention loses the exact atom |
| Dyadic jump tower | exact halving count and block countdown | changing centers can replenish it |
| Cylinder count | bulk density and endpoint phase | an error below one can contain the only survivor |
| Power sums | endpoint plus all even corrections | B1 alone predicts the wrong scale |
| Agoh--Giuga | parity and prime-factor restrictions | the odd composite case remains |

Prior-work checks recovered [THM-3224, complete LRC orbit Bernoulli
gcd-carry](../../01-canon/theorems/THM-3224-complete-lrc-orbit-bernoulli-gcd-carry-and-owner-hodge-splitting.md),
which already separates full-orbit Bernoulli data from cell ownership and
does not prove LRC. The corrected clipping issue near THM-739 in
[MISTAKES](../../01-canon/MISTAKES.md) similarly forbids exporting a full-circle
formula to an unmarked window. The existing
[Giuga circuit audit](collatz_mod6_20260921_grand_circuit_typing.md) already
derives the prime-factor criterion used below; this is not a new Collatz link.
The repository's HYP-2442 named Erdos--Moser concerns transitive subtournaments,
not the power-sum equation here; the shared name does not identify the problems.

The newly fetched [squares/doubles foundry](collatz_sqdbl_20260925_squares_doubles_foundry.md)
also requires two-place arithmetic and tests both signs and other multipliers.
The Bernoulli recoding below respects those controls; it supplies no automatic
solution of its new no-slow-divergence hypothesis.
Meta-pattern used: **Separate unbounded local support from a height-bounded
modular cover**. No new universal card is promoted.

## 1. Why the first odd coefficient is exceptional

Use the convention

    F(t)=t/(exp(t)-1)=sum_(n>=0) B_n t^n/n!,   B_1=-1/2.

The alternate convention has `B_1=+1/2`; state the convention when switching
between sums including left or right endpoints. The algebra gives

    F(-t)=F(t)+t,
    F(t)+t/2=(t/2)coth(t/2),                              (1)

so subtracting the single odd defect leaves an even function. Another exact
centering is

    exp(t/2) F(t)=t/(2 sinh(t/2)),                         (2)

whose odd coefficients, including the first, all vanish.
The standard generating function is recorded in
[DLMF 24.2](https://dlmf.nist.gov/24.2).

An endpoint derivation is equally informative. Bernoulli polynomials obey

    B_n(x+1)-B_n(x)=n x^(n-1),
    B_n(1-x)=(-1)^n B_n(x).

At `x=0`, the endpoint difference is 1 for `n=1`, and zero for `n>1`.
Reflection then forces every odd `B_n(0)` after `B_1` to vanish. This is an
exact symmetry with one endpoint defect, not an approximate large-index law.
See [DLMF 24.4](https://dlmf.nist.gov/24.4).

**Hostile to discarding the odd sector:** `B_3=0`, but
`B_3(1/4)=3/64`. The odd polynomials are not identically zero. Changing the
phase or origin can restore their values; the even Bernoulli numbers also
do not become a negligible tail merely because their indices are large.

Optional formal perspective: `log F(t)=-t/2+an even series`, with even
formal cumulants `-B_(2r)/(2r)`. These are not the Bernoulli numbers themselves.
Nor is F an ordinary real probability moment-generating function: its
formal second cumulant is `-1/12`. The exact interpretation is a translated
formal even series (or the reciprocal MGF of a uniform variable).

## 2. The B1 jump is a divisibility detector

Let

    psi(x)={x}-1/2,

where fractional part is in `[0,1)`; in particular `psi(integer)=-1/2`.
This is the periodic extension of `B_1(x)=x-1/2` in
[DLMF 24.2(iii)](https://dlmf.nist.gov/24.2#iii).

**PROVED.** For integer `z` and integer `m>=2`,

    J_m(z)=1/m+psi((z-1)/m)-psi(z/m)=1_(m divides z).       (3)

*Proof.* Write `z=qm+r` with `0<=r<m`. For `r>0` the difference of
fractional parts is `-1/m`; for `r=0` it is `1-1/m`. This covers negative
z as well, because Euclidean residues have the same range.

Thus for every nonzero integer z,

    v_2(z)=sum_(s>=1) J_(2^s)(z).                         (4)

The combined summands eventually vanish. More precisely, the first L
summands total `min(v_2(z),L)`. At zero they total L, so the infinite sum
diverges: the exceptional arithmetic case is retained, not regularized away.
The same proof gives `v_p(z)=sum_(s>=1) J_(p^s)(z)` for every prime p.

**Endpoint warning.** If one sets `psi(integer)=0` as in the symmetrized
Fourier/Dedekind convention, (3) is false. It becomes half the sum of the
indicators `m|z` and `m|(z-1)`. For example `z=m` gives one half instead
of one. An almost-everywhere equality of sawtooth conventions is inadequate
when the target is exact integer divisibility.

This is a genuine connection contract: source = divisibility of an integer;
target = an endpoint jump of periodic B1; map = (3); preserved predicate =
exact divisibility, hence valuations after summing scales. Forgetting the
endpoint value loses the atom. No statistical hypothesis is involved.

## 3. Applying the jump tower to Collatz blocks

For the odd Collatz map `U(n)=(3n+1)/2^v2(3n+1)`, (4) gives its exact
halving count. More substantially, take an actual exponent block
`w=(k_1,...,k_r)`, all `k_i>=1`. Write

    F_w(n)=(A n+C)/2^S,   A=3^r,   S=sum k_i,
    E_w(n)=(A-2^S)n+C,

with the signed ordered carry and source congruence retained. The preceding
note proves the exact repetition counter

    R_w(n)=floor((v_2(E_w(n))-1)/S),                       (5)

for positive odd n and nonzero E. At E=0 the block repeats forever.
Substituting (4) makes (5) a count of B1 jumps across dyadic scales.

The stronger scale-by-scale relation is exact. At every legal copy,
`E_w(F_w(n))=A E_w(n)/2^S`, and A is odd, so

    J_(2^s)(E_w(F_w(n)))=J_(2^(s+S))(E_w(n)).              (6)

The block shifts the entire divisibility profile down by S. This recovers
the countdown, not a new global descent theorem. At a change to block v,

    D_w E_v(y)=D_v E_w(y)+Delta,
    D_w=A_w-2^S_w,  Delta=D_w C_v-D_v C_w.                (7)

Equal dyadic valuations on the right can cancel and create many new levels.
The existing reset family `n_H=(2^(H+3)-13)/9`, `H=6j+5`, gives exactly
one `(1,2)` block to `2^H-1`, followed by `H-1` rises. Its first member is
`27->41->31`. The Bernoulli representation retains this obstruction exactly.

For a full actual orbit `x_0,...,x_r`, the exact drift identity is

    log_2(x_r/x_0)
      = r log_2(3) - sum_(j<r) sum_(s>=1) J_(2^s)(3x_j+1)
        + sum_(j<r) log_2(1+1/(3x_j)).                    (8)

The spatial mean of the valuation over positive odd residue classes is 2:
the density with valuation at least s is `2^(1-s)`. This follows since
`3n+1=0 mod 2^s` is one odd residue class. It does not establish the same
average along each deterministic orbit. Formula (8) exposes exactly the
temporal estimate still needed for descent; replacing its double sum by
`2r` without proof is the failed implication.

## 4. Boundary corrections can contain the whole survivor

Write positive odd sources as `n=2j-1`. A length-r exponent word of total
halving count S selects one class `j=a mod B`, where `B=2^S` and
`1<=a<=B`. Its exact count in `1<=j<=N` is

    C_(a,B)(N)=N/B+psi(-a/B)-psi((N-a)/B).                 (9)

*Proof.* The count is `floor((N-a)/B)-floor(-a/B)`; expand each floor.
The endpoint correction has magnitude strictly less than one.

That seemingly small correction can be all the remaining mass. For the
minus cycle `5->7->5`, t copies of `(1,2)` have `B=8^t,a=3`. At N=3,
the count is one for every t, while the bulk `3/8^t` tends to zero and
the correction tends to one. For the `5n+1` cycle `13->33->83->13`, use
`B=128^t,a=7,N=7`. The legitimate plus root 1 gives the same phenomenon
with `B=4^t,a=N=1`. The formula cannot distinguish allowed from forbidden
survivors unless the actual orbit and target are retained.

**PROVED bounded-correction obstruction.** No function

    V(n)=a log n+h(n),  a>0,  sup_(n odd positive)|h(n)|<infinity,

decreases at every positive odd Collatz step n>1. For `n_H=2^H-1`, the
first H-1 steps rise to `y_H=2*3^(H-1)-1`, and `y_H/n_H` tends to infinity.
Hence `V(y_H)-V(n_H)>=a log(y_H/n_H)-2||h||_infinity>0`, contradicting
the sum of the proposed stepwise inequalities.

This rules out finite sums of fixed periodic Bernoulli corrections and
uniformly absolutely summable collections. The unbounded sum (4) evades
this particular obstruction; the previous note's finite-center theorem
still rules out its simplest finite linear combinations with log height.
Adaptive centers, nonlinear combinations, and variable stopping times
are not excluded by either result.

## 5. A constructive next target: retain complete boundary packets

Bernoulli multiplication gives the exact cancellation

    sum_(j=0)^(m-1) psi((x+j)/m)=psi(x).                  (10)

Proof: both sides have the same unit-slope evolution and the same jumps;
or sum the fractional parts and check on one period. Thus all children
of a residue class can be combined without discarding their parent boundary.
This is a useful operation for a proof search, not permission to cancel
an incomplete set of children.

For a plus affine block, a first-descent certificate must retain the guard

    A<B and n>C/(B-A),   as well as the actual source cylinder. (11)

Proposed OPEN research target: construct parameterized packets of
first-descent cylinders whose residual boundary terms under (10) reduce
to an already certified smaller packet, with every cutoff in (11) retained.
The first meaningful success would cover a new infinite residual family
of the existing sibling grammar. Universal packet closure would need an
independent well-founded rank or coverage proof; it does not follow from
the identity. Test the reset family and the nontrivial minus/5x+1 cycles
before interpreting any finite cancellation as progress toward convergence.

## 6. Erdos--Moser: the boundary term survives the macro limit

The power-sum conjecture asks whether

    1^k+2^k+...+(m-1)^k=m^k

has any positive integer solution other than `(m,k)=(3,1)`.
For positive integer k and m>=2, Faulhaber's identity gives exactly

    sum_(j=1)^(m-1) j^k/m^k
      = m/(k+1)-1/2
        + sum_(r=1)^floor(k/2) [B_(2r)/(2r)!]
          * k^(falling 2r-1)/m^(2r-1).                   (12)

At a solution the right side equals one. This is a concrete realization
of the user's intuition: the single B1 term remains `-1/2` after
macroscopic normalization. But the even hierarchy remains essential too.
Keeping only the first two terms predicts `k/m -> 2/3`, which is wrong.

For any hypothetical unbounded family of integer solutions, the actual
limit is `k/m -> log 2`. Here is an elementary justification. Rewrite the
left side as `sum_(h=1)^(m-1)(1-h/m)^k`. Since `(1-h/m)^k<exp(-kh/m)`,
a solution has `k/m<log 2`. The ratio cannot tend to zero: arbitrarily
many fixed h terms would then tend to one. Along any subsequence with
positive ratio limit lambda, a geometric dominating series applies, and
the sum tends to `1/(exp(lambda)-1)`. Equating this to one gives lambda=log2.
All subsequential limits therefore agree.

Equivalently, for this small lambda the complete even Bernoulli hierarchy
resums to

    1/lambda-1/2+sum_(r>=1) B_(2r) lambda^(2r-1)/(2r)!
      =1/(exp(lambda)-1).

The exact power-sum identity was checked independently below; the
asymptotic is also classical. Existing results go substantially further:
Gallot--Moree--Zudilin prove

    k=log2 * (m-3/2-(25/12-3log2)/m+O(m^-2))

and show that a nontrivial integer solution forces `2k/(2m-3)` to be a
continued-fraction convergent of log2. See Theorem 1 and Corollary 1 in
[their paper](https://arxiv.org/pdf/0907.1356). These are CITED prior results,
not new work here. The remaining research target is an arithmetic
incompatibility for all admissible convergents, retaining the endpoint
shift `2m-3`; deriving more asymptotic terms alone is not such a proof.

## 7. Agoh--Giuga: a literal role for the exceptional first coefficient

The conjectural criterion `n B_(n-1)=-1 mod n` characterizes primes.
Here the Bernoulli convention makes the even part immediate:

    n=2:       2 B_1=-1;
    even n>2:  B_(n-1)=0, so the congruence fails.

For odd n>1 the remaining criterion is equivalent to n being squarefree
and, for every prime divisor p,

    p-1 divides n-1,      n/p=1 mod p.                    (13)

These are the simultaneous Carmichael and Giuga restrictions for composites;
see [Borwein--Maitland--Skerritt, Theorem 3](https://www.carmamaths.org/resources/jon/giuga2013.pdf)
and the inherited repository audit cited above. Either condition alone
is insufficient: 30 is Giuga but fails the even test; 561 satisfies the
Carmichael condition but `561/11=51=7 mod11` violates the Giuga condition.
The Bernoulli endpoint insight settles this parity split, not the remaining
odd composite exclusion. No map transporting that exclusion to Collatz
has been supplied.
For an odd squarefree candidate, the individual conditions in (13) are
exactly the jump conditions `J_(p-1)(n-1)=J_p(n/p-1)=1`. Thus the same
endpoint detector has a literal role here too; rewriting the gates does
not establish that no composite can satisfy all of them.

## 8. What was gained, and reproducibility

The useful transfer is now explicit: the exceptional endpoint jump of B1
encodes exact divisibility; its dyadic copies encode Collatz clocks and
the earlier rational-center countdowns. Full residue packets admit exact
Bernoulli cancellation, while an omitted endpoint atom can contain an
entire exceptional orbit. This locates a concrete proof obligation rather
than claiming that symmetry itself proves termination.

Run:

    python3 04-computation/experiments/bernoulli_boundary_20260925.py
    python3 04-computation/experiments/bernoulli_boundary_20260925_audit.py

The [main output](bernoulli_boundary_20260925.out) records exact controls:
degrees 1..40 for endpoint/reflection identities; 129,087 divisibility
checks (`-1024<=z<=1024, 2<=m<=64`); 4097 valuation inputs at 13 scales;
51,216 cylinder counts (`1<=B<=32,1<=a<=B,0<=N<=96`); 552 complete-child
identities; three periodic controls at 20 depths; and 1560 normalized
power sums (`1<=k<=40,2<=m<=40`). No floating point is used in these gates.

The [independent audit](bernoulli_boundary_20260925_audit.out), written
without the main script, checks 19,899 divisor indicators, 2010 truncated
valuations including zero, 1064 power sums, and endpoint-convention hostiles.
Separate agents audited the proof scopes and primary-source connections.
[Manifest](bernoulli_boundary_20260925_manifest.json) records raw LF hashes.
Finite controls corroborate the elementary proofs; they prove none of the
open universal statements. No new theorem ID or Lean claim is introduced.
