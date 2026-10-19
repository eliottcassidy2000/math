# Local Collatz return ranks, exact prime budgets, and why a finite bank can recharge them

**Status, 2026-10-19:** PROVED elementary structural results, with exact
source controls. This note supplies a quantitative obstruction to a proposed
coverage mechanism. It does not assert that no stronger global rank exists,
and it does not prove convergence of every positive integer.

## 1. Inheritance and the new question

The [fixed-translation packet](fixed_translation_descent_20261019.md) bounds
isolated integer collisions by d/2+2^{N−c sqrt(N)+O_d(1)} and supplies an
actual smaller-source receipt when min(a,b)≤2x−d. The unresolved task is
source coverage, not another finite residue census.

We inherit the actual affine Terras pair chain from
[THM-4569](../../01-canon/theorems/THM-4569-the-terras-clock-recurrence-of-two-collatz-orbits-is-unconditional.md),
the isolated-stabilizer boundary from the
[integer stabilizer packet](integer_stabilizer_20261019.md), and the
positive-source/negative-center distinction from
[THM-4594](../../01-canon/theorems/THM-4594-the-maximal-class-decided-collatz-sieve-and-the-sign-barrier.md).
The [finite certificate bank](collatz_finite_bank_certificate_20261007.md)
already separates an available symbolic word from the word followed by a
fixed source. None of those mechanisms is re-audited here.

The new results are:

- exact dual budgets: forward repetition of one branch block spends a
  2-adic valuation, and inverse integer lifting spends a 3-adic valuation;
- every fixed or ultimately periodic inverse return code has only finitely
  many distinct integer lifts of one fixed seed;
- a fixed bank of just two actual adjacent-pair return words can recharge
  either local valuation rank arbitrarily through a change of center;
- the exact source-realization test for relation-state cycles explains why
  nonzero carry and finite-state recurrence alone do not settle coverage.

The anchor is a source-wise rank. The niche is changing its rational center.
The wildcard is the duality between forward parity and inverse integrality.
The hostile is the positive phase pair (2,1); the corrected near miss is
calling separate decreasing local ranks a common decreasing rank.

Throughout T is the Terras map. All positivity and integrality tests below
refer to actual sources, not to uniform counting or Haar measure.

## 2. An integer normalized-carry automaton

For a relation u=3^k v+e between integer endpoints, write

    E=3^{max(0,−k)}e ∈ Z.

The driver bit is β=v mod 2 and the parity discrepancy is σ=E mod 2.
The exact integer transitions are as follows.

For even E:

    β=0: (k,E) → (k,E/2),
    β=1: (k,E) → (k,[3E+3^{max(0,−k)}−3^{max(0,k)}]/2).

For odd E and β=0:

    k≥0: (k,E) → (k+1,(3E+1)/2),
    k<0: (k,E) → (k+1,(E+3^{−k−1})/2).

For odd E and β=1:

    k≥1: (k,E) → (k−1,(E−3^{k−1})/2),
    k≤0: (k,E) → (k−1,(3E−1)/2).                   (2.1)

These are obtained by clearing the powers of 3 in the inherited pair-chain
table. They introduce no rounding. If k<0, an integer source also needs
v≡−E mod 3^{−k}; if k≥0 there is no extra denominator guard. Always retain
v>0 and u>0. The graph state by itself forgets these source conditions.

For an initial integer translation d and an isolated merge history, E is
nonzero throughout and is a 3-adic unit whenever k≠0, as proved in the
preceding packet. If both odd counts are bounded by K, then
|E|<(d+1)3^K and |k|≤K. This gives a finite graph with an odd-count budget,
but it does not turn unrestricted paths into a finite acyclic graph.

## 3. Exact dual budgets for one parity block

Let w be a parity block of length n≥1 with a≥1 odd bits and affine branch
map

    F_w(x)=(3^a x+C)/2^n,
    p=C/(2^n−3^a).                                 (3.1)

The denominator 2^n−3^a is odd and coprime to 3. Thus p lies in both Z_2
and Z_3. Unique factorization ensures the denominator is nonzero.

### Forward parity budget

For an integer source x and an integer j≥0, its first nj parity bits are
w repeated j times **iff**

    v_2(x−p)≥nj.                                   (3.2)

At x=p, interpret the valuation as infinity. In particular, if x≠p, the
number of complete repetitions is exactly floor(v_2(x−p)/n).

**Proof.** The formal iterate satisfies

    F_w^j(x)−p=(3^a/2^n)^j(x−p).

The parity-vector cylinder of w^j is exactly the residue class on which
this affine iterate is 2-adically integral. This is the usual branch-word
bijection, or follows by reducing the numerator congruence at each prefix.
Since p is 2-adically integral, that integrality is equivalent to (3.2). ∎

There is a direct source-height bound. If x≠p, put D=2^n−3^a. Then

    nj≤v_2(Dx−C)≤log_2|Dx−C|.                      (3.3)

The integer in the absolute value is nonzero. A fixed positive integer
cannot follow a negative-center periodic block indefinitely, even though
arbitrarily long finite prefixes have positive integer representatives.

### Inverse integer budget

Let G_w=F_w^{-1}. For an integer seed y and j≥0,

    G_w^j(y)=p+(2^n/3^a)^j(y−p) ∈ Z
        iff v_3(y−p)≥aj.                          (3.4)

**Proof.** Every iterate of G_w has denominator dividing 3^{aj}, so it is
an integer iff it is 3-adically integral. If v_3(y−p)<aj, the second term
in (3.4) has negative valuation and cannot cancel the 3-adically integral
p. If v_3(y−p)≥aj, both terms are 3-adically integral. ∎

Every integer inverse obtained this way follows the specified branch word;
positivity is a separate guard. The local rank v_3(y−p) drops by exactly a
at each inverse step. The forward rank v_2(x−p) drops by exactly n at each
legal forward block.

### Consequence for repeated pair-return constructions

A nonempty oriented return word for a pair with fixed nonzero integer
translation d has the same positive odd count a on both paths. Indeed a=0
would divide the difference by 2^n instead of returning it to d.
Consequently (3.4) applies to every such word.

Starting from one fixed integer seed y≠p, repeating one inverse return
word supplies at most floor(v_3(y−p)/a) integer lifts. If y=p is an integer,
the lift is constant. The same conclusion holds for any ultimately periodic
sequence of inverse return words: compose one period, whose total odd
count is positive, and apply (3.4). A fixed point of that composition gives
only the finitely many sources within one period, not infinitely many
distinct lifted integers.

This is an exact obstruction to the simplest infinite-family construction.
It does **not** rule out an aperiodic sequence from a finite bank.
The all-even inverse map x↦2^n x is an important boundary case: it does
create infinitely many integers, but it scales a nonzero translation and
therefore is not a return word for fixed d.

## 4. A fixed two-word bank can recharge the rank arbitrarily

Use these actual oriented returns of adjacent pairs:

    F(x)=(3x+1)/4,       x≡1 mod 4,
    H(x)=(81x+101)/512,  x≡43 mod 512.              (4.1)

For every source in the indicated cylinder, the upper source x+1 follows
its corresponding word and ends at F(x)+1 or H(x)+1. The lower words have
lengths 2 and 9, with 1 and 4 odd steps. Their centers are

    p_F=1,     p_H=101/431,
    v_3(p_H−p_F)=v_2(p_H−p_F)=1.                   (4.2)

The checker reconstructs these coefficients from actual trajectories; no
formal word is assumed lawful merely because its coefficients look right.

### Arbitrarily large backward recharge

Fix any M≥4. By the Chinese remainder theorem choose a positive integer z
such that

    z≡1 mod 4,   v_3(z−p_H)=M.

Such z exist: prescribe z modulo 3^{M+1} to agree with p_H through level M
but differ at the next digit, then impose the independent mod-4 condition.
Set y=F(z). Since M>v_3(p_H−1)=1, the ultrametric identity gives

    v_3(z−1)=1,    v_3(y−1)=2.                      (4.3)

Thus F^{-1} takes the integer seed y to z, spending one unit of its local
rank. Switching to H changes the available rank from 1 to the arbitrarily
large value M, allowing exactly floor(M/4) integer inverse H blocks.
All these inverse sources are positive: z>1>p_H and the inverse multiplier
512/81 is greater than 1.

For example, M=4 can use y=613, z=817, followed by H^{-1}(817)=5163.
The F-rank was only 2 at the starting seed, while the H-rank after switching
is 4. The construction works for every M, not just this example.

### Arbitrarily large forward recharge

Dually, fix M≥9 and choose a positive integer z with

    z≡1 mod 3,   v_2(z−p_H)=M.

Set x=F^{-1}(z)=(4z−1)/3. It is a positive integer, and

    v_2(x−1)=3,   v_2(z−1)=1.                       (4.4)

The actual F block uses two units of the fixed initial rank 3. At z,
switching to H supplies M units, permitting exactly floor(M/9) successive
H blocks. All are actual positive integer paths. For M=9, one example is
x=2105, z=1579.

These constructions refute the inference that a collection of decreasing
word-specific valuation ranks is already a common well-founded rank.
They do not refute a different global potential. In fact F and H both
strictly lower ordinary positive integer height above their small centers;
this bank already has a real-height rank on its admitted domains. What it
lacks is all-source entry into those domains. The established affine-threshold
sieve covers that distinction; local rank recharge does not solve it.

The constructed seeds vary with M. No infinite inverse chain from a single
fixed seed has been claimed, and these are return histories rather than an
assertion that every constructed pair is an isolated atom.

## 5. What graph cycles actually mean for positive integer sources

A periodic driver word w repeated forever has the unique source p in
(3.1). For a≥1, its actual positive integer realization requires

    2^n>3^a  and  (2^n−3^a) divides C.              (5.1)

Both conditions are necessary and sufficient for that single trajectory;
its upper partner must also satisfy the retained affine integrality and
positivity guards. This is a source test, not a finite-state cycle test.

Three controls make the distinction concrete.

1. **Bounded negative-center loop with positive finite prefixes.** From
   (k,E)=(0,1), driver bit 0 reaches (1,2). Driver bit 1 then fixes (1,2)
   forever. Its driver center is −1; the complete initial lower source is
   −2. For every L, the positive initial lower source 2^{L+1}−2 realizes
   the finite driver prefix 0 followed by L ones. Thus arbitrary positive
   finite prefixes do not supply a positive source for the infinite loop.
   Formula (3.3) gives the exact repetition budget for any fixed source.

2. **A genuine positive nonidentity phase cycle.** The pair (2,1) realizes
   (0,1)→(−1,1)→(0,1), with both normalized carries nonzero 3-adic units.
   The endpoints never merge at equal clocks. A phase-aware receipt can
   handle this known cycle, but the nonzero-carry property itself cannot.

3. **Conditional embedding of every positive Collatz cycle.** If any
   positive T-cycle has two distinct phases, choose them as p+d and p with
   d>0. After its full period both paths return, their odd counts agree,
   and the relation state returns to (0,d) without literal coalescence.
   It is therefore a legal positive nonidentity cycle in this automaton.
   This includes the known 1↔2 cycle; excluding all additional positive
   cycles of the relation automaton already includes the positive Collatz
   cycle-exclusion problem. No additional cycle is asserted to exist.

Finite-state recurrence also permits aperiodic paths that switch cycles.
Repeated vertices do not imply a repeated parity word. The quantitative
local budgets in §3 therefore cannot be applied to every recurrent graph
path without controlling the changes of center. Section 4 shows that such
changes can recharge either valuation budget by an arbitrarily large amount.

## 6. Consequence for the remaining coverage problem

For a hypothetical least nonconvergent positive integer x, the preceding
packet implies that an isolated merge with x+d, for x>d/2, must have
min(a,b)>2x−d; otherwise both sources receive a receipt to z<x. This new
necessary condition does not bound the number of later odd events.

The exact additional work needed is now explicit:

- retain the actual source and its center-relative valuations, not merely
  the normalized carry graph;
- control valuation recharge under arbitrary changes of return word, or
  replace it by a genuinely common source-height rank;
- prove that every unresolved source enters the admitted decreasing bank;
- separately certify phase cycles and test any proposed periodic source
  with the integer and positive guards in (5.1).

The local ranks are real well-founded mechanisms for one word or one periodic
code. The recharge examples explain precisely why they do not yet extend
to a source-covering proof. Haar recurrence, finite-state recurrence and
arbitrarily long positive finite prefixes supply none of these missing
quantifiers.

## 7. Focused checks

The [checker](../experiments/collatz_return_rank_obstruction_20261019.py) performs **3,278
explicit checks**, with identical [normal](collatz_return_rank_obstruction_20261019.json)
and [optimized](collatz_return_rank_obstruction_20261019.optimized.json)
outputs. It checks the actual return coefficients and normalized transition
table, both valuation equivalences, the exact exhaustion of each recharge
budget, positive source guards, and the two explicit graph-cycle controls.

As finite positive controls, F^{-1}(7)=9, H^{-1}(7)=43 and F^{-1}(43)=57
transport the isolated adjacent atom at 7 to actual first merges at depths
12,19,21 respectively, always with odd-count difference −1. These illustrate
legal atom transport and the fixed-seed rank limit, not infinite coverage.
There is no new residue census and no dependency on uncomputed Collatz tails.

```sh
python3 05-knowledge/experiments/collatz_return_rank_obstruction_20261019.py > /tmp/return-rank.json
python3 -O 05-knowledge/experiments/collatz_return_rank_obstruction_20261019.py > /tmp/return-rank-O.json
cmp /tmp/return-rank.json /tmp/return-rank-O.json
cmp /tmp/return-rank.json 05-knowledge/results/collatz_return_rank_obstruction_20261019.json
```
