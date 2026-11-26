# Half-rank averages, the height measure, and signed physical lifts

2026-11-26. **PROVED** elementary height law and moment reductions; **FINITE-EXACT**
physical phase tests. No claim of a proof of the all-curves minimalist
conjecture or of the full marked Moore cover.

## 1. Two rank questions, and an exact measure change

The [October finite-growth paper](../bsd_selmer_growth_20261013a.md) gives a
two-level certificate for Selmer corank at most one and a Monsky matrix
interpretation in residue tournaments. The [twist-pair reading](../five_manuscripts_chowla_collatz_residue_bsd_20261013.md)
retains both local conditions and a uniform passage to a growing family.
The [rank-return session](../overnight_rank_return_20261029/ELLIPTIC.md) makes the
complementary point: individual large-rank certificates do not determine a tail.

These are different ensembles:

- **One base fixed; signed squarefree twists grow.** [Smith](https://arxiv.org/abs/2503.17619)
  proves the half–half distribution of 2-infinity Selmer corank, so the
  rank-at-least-two locus already has density zero in that family by the
  Kummer rank bound. The later [analytic-density and mean-rank manuscripts](https://github.com/openai/math/blob/main/overview.tex)
  claim the half–half analytic density and mean analytic rank 1/2 for
  each fixed base. We do not use those preprint claims to assert a
  theorem about a different ensemble.
- **All curves over Q grow by naive height.** The density-zero and
  half-rank mean are still the general minimalist prediction here.
  [Bhargava–Shankar's bounded mean](https://arxiv.org/abs/1312.7859)
  and their [first 2-Selmer moment](https://arxiv.org/abs/1006.1002)
  are unconditional inputs.

Here is an exact reason the first ensemble does not transfer by averaging.
Represent each curve by its unique short integral pair (A,B) with no prime
simultaneously dividing A to order at least 4 and B to order at least 6,
and H=max(|A|^3,|B|^2). Discard A=0 or B=0; these contribute O(X^(1/2))
models, negligible against X^(5/6). The choice of positive constants
multiplying the two terms in H does not change the limiting law below.

Define a **twist-minimal base** to have B0>0 and no p with p^2|A0 and p^3|B0.
Every remaining curve is uniquely

    (A,B) = (d^2 A0, d^3 B0),    d signed squarefree,
    H(A,B) = |d|^6 H(A0,B0).

In fact min(floor(v_p(A)/2),floor(v_p(B)/3)) is 0 or 1 by minimality;
its nonzero primes give |d|, and sign(B) gives the sign. Conversely any
such base and d give a minimal pair. The restriction B0>0 removes the
sign ambiguity, not the negative twists.

Let M(Y) count the bases up to height Y. Möbius inversion in the box
|A0|<=Y^(1/3), 1<=B0<=Y^(1/2) gives

    M(Y) ~ (2/zeta(5)) Y^(5/6).

For every fixed divisor k the limiting proportion of multiples of
(k^2,k^3) is k^(-5); floors divided by Y^(5/6) are bounded by k^(-5).
The zero/discriminant-zero exclusions are negligible. The same proof
with (k^4,k^6) gives (4/zeta(10)) X^(5/6) for all curves of both B
signs. Equivalently, sum 2 M(X/n^6) over positive squarefree n and use
sum_mu^2(n)n^(-5)=zeta(5)/zeta(10). The latter representation has an
absolute normalized cutoff error O(D^(-4)) when n>D.

Thus the limiting height law, for squarefree n, is

    Pr(|d|=n) = (zeta(10)/zeta(5)) n^(-5).
    Pr(|d|=1) = zeta(10)/zeta(5) = 0.965346496... .

The Euler product makes the divisibility bits at any finite set of
primes *independent in this limiting law*, with Pr(p divides d)=1/(p^5+1).
The sign of d is an independent fair bit by the map B to -B. Neither
assertion implies that the root number or rank is independent of d.

This gives a concrete test of the quantifier. A fictitious exceptional
property "|d|=1" occurs only twice in each fixed-base twist family and
yet has limiting density 96.5% by height. Even a property with global
density one can occur only finitely often on every row: take, for
example, |d|<=H(A0,B0); for each fixed n all but finitely many bases
then have the property, and the uniform tail O(D^(-4)) applies.
Neither example is a claim about real elliptic-curve ranks.

It follows more precisely that vanishing density of rank>=2 in the
height ensemble is equivalent to vanishing density **on every fixed
n slice** (or each of the two signs separately) as the base height grows.
Forward follows from each slice having positive limiting mass; backward
follows by truncation in n and the uniform tail. This reverses the
usual fixed-base quantifier. In particular the twist-minimal slice at
n=1 is the first decisive test, not a negligible subfamily.

## 2. A sufficient arithmetic moment to aim for

Write r for Mordell–Weil rank and use the October growth identity at a
fixed prime p. With torsion removed from full finite Selmer sizes as
there, let a_m=u_(m+1)-u_m, and let s be the p-infinity Selmer corank.
Then a_m decreases to s and is always at least s>=r. If b is the parity
bit of s, a nonnegative exact excess is

    Delta_m = p^(a_m) - p^b .
    Pr(s>=2) <= E[Delta_m]/(p^2-1).

If s is even and >=2 the contribution is at least p^2-1; if odd it is
at least p^3-p. This is the average version of the two-adjacent-level
certificate, without assuming Sha is finite. Thus a useful sufficient
target in the height ensemble is the parity split E[b]->1/2 plus

    lim_(m->infty) limsup_(X->infty) E_X[p^(a_m)] <= (p+1)/2.

This is **not** implied by pointwise eventual stabilization of a_m.
If a low-corank converse is available on the resulting density-one set,
it identifies algebraic rank with b there. The first Selmer moment supplies
the missing arithmetic **mean-tail** step without any high-rank BSD input:
rational points inject modulo 2 into Sel_2, so E_X[2^r] has bounded limsup.
For the unweighted height of Section 1 this boundedness follows by
enclosing its boxes in a fixed enlargement of the height boxes in the
cited theorem and using the X^(5/6) counts. No exact value for the
unweighted moment is needed.
If delta=Pr_X(r>=2) and E_X[2^r]<=C, conditional Jensen gives

    E_X[r ; r>=2] <= delta log_2(C/delta).

For delta=0 the right side is interpreted as zero. Thus once the
rank-zero/rank-one density split is known by height, the mean follows;
rare high ranks cannot carry positive mass in the mean under this bound.
The hard input left is precisely the height-uniform passage to low
corank (or some other proof of the rank density). The standard-height
2-Selmer moment of 3 is not the target persistent moment of 3/2.

There is also a simple second route to density zero **if** the conjectural
identity E_X[#Sel_p]->p+1 can be proved along an unbounded sequence of
primes p. For each fixed p, p^r-1 >= (p^2-1) 1_(r>=2); hence the density
limsup is at most 1/(p+1). Send p to infinity only *after* taking the height
limsup. This route does not supply the parity split by itself.

## 3. On the Moore board, the lift becomes an exact signed graph

The arithmetic points above transfer a method, not a theorem about graphs:
keep a reduced family, and measure which configurations survive a lift
rather than counting reduced configurations as solutions. This extends
the [physical slot test](HADWIGER-NELSON-PHYSICAL-SLOTS.md) and exploits the
[Thue–Morse switches](THUE-MORSE-TRADE-CUBES.md).

Keep the marker at (2,7), and group an exterior branch's admissible
contacts at y and y+28, if its center maps respect the half-turn.
Keep the choices at exceptional branches 27 and 49 unmerged. Rebuilding
the exact board gives **1,064 options** for the same 1,741 physical
contacts. Of these, 677 allow both phases and 387 have one legal phase.
The center capacities on residues modulo 28 are 73 twos, ten ones, one
zero. Each exact residue profile has 52 chosen options.

For a center slot of the representative write s=r+28 eta, and use a bit
b_i for the selected contact's phase. A double-capacity residue shared by
branches i,j demands

    b_i xor b_j = 1 xor eta_i xor eta_j.

A singleton or an option allowing one phase imposes a unary equation.
Encode those using an added fixed-zero vertex. The signed multigraph
on 53 vertices then has 83+f edges, where f is the number of chosen
one-phase options. If it has c components including the fixed-zero
vertex's component, the exact cycle rank is **30+f+c**. The profile
lifts iff every cycle has zero XOR (including parallel edges).
In that case there are exactly 2^(c-1) physical lifts. The proof is
propagation along a spanning forest, with one free root bit in each
unpinned component. If edge labels were independent fair bits, the
probability of consistency would be 2^(-(30+f+c)); this is a control,
**not** an assumption about the actual chart. Counting residue profiles
and counting physical covers therefore have different weights, just as
fixed-base twist density and height density do.

### A short certificate on actual contacts

Among the 677 two-phase options, the independent census finds **607**
pairs that share two or more center residues and already contradict
on parallel XOR checks. Indeed, equality of the two phase differences
would make two distinct physical contacts share two slots, violating
linearity. This forces opposite differences at the shared residues.

There are also **9,449** inconsistent triangles with every pair sharing
just one center residue and the three overlaps using three distinct
centers. No pair alone obstructs these triangles. For example take the
following canonical residues and slots:

| Exterior option (branch,residue) | Three center representatives |
|---|---|
| (3,2) | (2,48,39) |
| (4,2) | (2,35,32) |
| (54,18) | (18,20,4) |

The three equations are b_3 XOR b_4=1, b_4 XOR b_54=0,
and b_3 XOR b_54=0. The three residues are each at capacity two,
and all three options have both physical phases in their domains.
The obstruction is exactly the cycle sum. Each such configuration
gives a valid two- or three-option exclusion in the residue solver.

An exact-at-each-node search using these phase equations rejects
partial selections immediately on contradiction. Bounded runs here
did not produce a full residue-and-phase cover. This is no
infeasibility result for the full board or for the Moore graph.

## Replay

```sh
python height_twist_count.py
python nested_core_check.py --board > board.txt
python antipodal_options.py board.txt > signed.txt
python signed_small_cycles.py board.txt
c++ -std=c++17 -O3 antipodal_solver.cpp -o /tmp/antipodal
/tmp/antipodal phase 1000000 4 < signed.txt > profile.txt
# If a run succeeds, independently validate the lifted physical word:
# python antipodal_options.py board.txt profile.txt
```

The height census compares the formula with direct minimal-pair
enumeration at four low heights and compares the two exact counting
formulas up to 10^18. It counts twist labels, not ranks. Arithmetic
status, the full marked cover, and the existence of the degree-57
Moore graph are separate questions.
