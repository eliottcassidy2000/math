# Forbidden sets, complete coefficient rays, and costed recursion

**PROVED scoped reductions and two complete infinite ray families;
FINITE-EXACT coefficient band; Chan–Pak OPEN.** No finite search below
establishes that the only globally forbidden values are7 and21.

## F1 The negative formulation and its exact stopping test

Let S_n be the raw Hamilton-path counts of tournaments of order at most n,
and F_n the positive odd integers outside S_n. An unselected source or sink
preserves H, so S_n increases and F_n decreases. The conjecture is precisely

```
intersection_(n>=1) F_n = {7,21}.
```

Each F_n is infinite: a fixed finite order has bounded H. The convergence
must be pointwise, or on each fixed finite interval. It cannot mean that a
finite F_n eventually equals{7,21}.

The inherited logarithmic-core theorem gives a rigorous negative stopping
rule. Remove all singleton strong components from a realization of N>1.
The remaining order is at most floor(3 log_3 N); equality is attained by
an ordered chain of directed triangles. Therefore, for integers N>1,

```
N <= 3^(n/3)  implies  [N in F_n iff N is globally forbidden].       (F1)
```

N=1 is separately realized by one vertex. This is a corollary of the
existing [finite-core proof](../spectrum_basins_20261020e/spectrum.md),
not a new logarithmic bound. Busch's exact strong minimum, reproduced
there from his primary paper, gives sharper divisor-specific cutoffs.
Thus global absence of any particular integer has a finite decision
procedure. The unresolved part is a uniform description of all decisions.

For example7 forces a strong factor7 of order at most4, and21 forces a
strong factor7 or21, the latter of order at most6. Complete small inventories
exclude these. This is the inherited all-order exclusion, not extrapolation
from the new order-ten census.

## F2 Separate value factorization from physical endpoint realization

The existing [THM-4094, two-prime-lane completeness](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md)
reduces the conjecture to strong realization of all odd primes except7,
and all7p for odd primes p other than3, with the already realized carry
atoms49,63,343. Its proof uses ordered-SCC multiplication. A least allowed
counterexample beyond the carry atoms must consequently be prime or7
times a prime. This removes many composite cases from a negative search;
it does not construct either remaining lane.

The exact minimum-order recursion is also inherited from
[Session224, first-occurrence fibers](../session224_first_occurrence_20261107/proof/PROOF.md):
if s(d) is the least strong order at raw count d, then

```
b(N)=min_(d|N,d>1) [s(d)+b(N/d)],    b(1)=0 internally.
```

The singleton realizes the external value1. An absent strong factor has
cost infinity. Products are ordered component words when counting types;
the scalar minimum discards those words, their multiplicities and marked
endpoint states. Keep those data when continuing a construction.

## F3 A complete finite coefficient band, with every lower head paid

Write each N>1 uniquely as N=1+c*2^s, with c odd and s>=1. A unit seed
of slope U=c*2^j supplies exactly the entire tail s>=j+1 by the inherited
physical tail/apex construction. Its finite unpaid heads are s=1,...,j.
The pair(c,j), not c alone, is the useful coverage state.

**FINITE-EXACT.** Among all478,684 ordinary strong marked unit seeds of
order10 there are37,471 distinct odd parts of U. Every odd c<=53,557
occurs;53,559 is the first absent odd part in this declared bank. This is
not a nonexistence statement for larger seeds or other tournament classes.

For each c choose an actual seed with minimum j. The26,779 prefix witnesses
are stored and independently recounted in [coefficient_band.py](coefficient_band.py).
Their minimum-valuation histogram is

| j |0|1|2|3|4|5|6|7|8|9|10|11|
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| coefficients |24540|1297|523|231|98|49|23|10|3|2|2|1|

There are3,954 distinct unpaid(c,s) heads. Their largest value is43,131.
Every one appears in the inherited complete order-eleven spectrum except
7 and21. The exact spectrum file is hash-pinned; its least-order entries
are consulted, not inferred from a contiguous-coverage headline.

Consequently the following is proved using exact finite certificates and
the all-length tail identity:

> For every odd c<=53,557 and every s>=1, the integer1+c*2^s is a
> tournament Hamilton-path count unless it is7 or21. A witness has order
> at most s+11; finite heads already have witnesses of order at most11.

This is infinite coverage along a finite coefficient band. Its number of
values below X is O(53,557 log X), so it does not supply positive density
or a proof for all coefficients. It does not enlarge the previously
verified general bounded range of the conjecture.

## F4 A renewable1039 family, and why99 closes its head

Provenance boundary:1039 was already supplied by the inherited ordinary
order-eight seed rows[22,36,8,19,38,77,159,63], O={4,7}, with H=137,
delta=1,U=1039 and no stable port. The new contribution here is the
renewable stable-port family, not the first realization of coefficient1039.

Take the ten-vertex strong marked birth

```
rows=[4,17,26,19,673,15,63,111,255,495], O={2,3,4}, port5.
```

The complete signed port kernel is F_1, and its full slope is
U(a)=98+990*2^a for every a>=1. Hence

```
c_a=49+495*2^a,       c_1=1039,
c_(a+1)=2*c_a-49,    U(a)=2*c_a.                         (F2)
```

Inflating the same port is a physical seed-preserving recurrence, with
source order a+9. It supplies1+c_a*2^s for every s>=2. The first head is

```
1+2*c_a=99*(1+5*2^(a+1)).                               (F3)
```

This factorization gives an all-length physical repair. The five-vertex
seed with rows[2,4,9,19,7], selection{0,1}, has H=9, delta=1 and U=20;
its tail realizes1+5*2^r for every r>=3. A separate seven-vertex tournament
with rows[42,44,9,16,39,76,31] has H=99. Their ordered join realizes(F3)
for every a>=2, of order a+12.

At a=1 the factorization is99*21 and cannot be used:21 is forbidden.
The explicit nine-vertex rows[320,25,35,389,141,27,62,103,246] instead
have H=2079. This is a precise example of a factor route failing while
the target remains realizable. Therefore every head and every tail in
(F2) is now constructed: all a>=1 and all s>=1, with no exceptions.

The occurrence of99 here is algebraic, not an identification with the
separate99*2^k+23 family:2*49+1=99 and2*495=10*99 are exactly the
divisibilities that make the lower-head factorization work.

## F5 The crossing obstruction also gives a complete infinite ray family

[CROSSINGS.md C3](CROSSINGS.md) supplies unit seeds with

```
d_k=5841*2^k-241,    U=4*d_k,    k>=1,
d_(k+1)=2*d_k+241.
```

The stable line gives every exponent s>=3. A new search for the *actual
missing interfaces*, rather than more unit seeds, closes both heads.

* Rows[6,200,50,293,267,83,157,61,231], selection{3,6,7},
  have H=1889, defect−481 and U=11682. A tail of length k realizes
  −481+11682*2^k=1+2*d_k for every k>=1.
* Rows[94,156,24,16,32,79,158,61], selection{1,6,7},
  have H=119, defect−107 and U=1298. A tail of length k+1 followed
  by ordered joining with the five-vertex H=9 tournament realizes
  9*(−107+1298*2^(k+1))=1+4*d_k for every k>=1.

Thus every1+d_k*2^s is constructed for every k>=1,s>=1. The first
head construction has order k+10, the second k+15, and the unit-tail
construction for s>=3 has order k+s+10. These are valid upper bounds,
not minimum-order assertions.

Both interfaces are checked by independent complete endpoint-table
algorithms. [complete_rays.py](complete_rays.py) emits actual adjacency
rows for either family and any positive parameters;53 entire graphs in
its declared grid are also counted literally. The discovery search in
[head_repair.py](head_repair.py) checks all ordinary tournaments through9
and all strong order10 sources admitting a unit selection, while allowing
every selection for the new nonunit interfaces. It finds28 presentations.

The general useful rule is exact: an interface(delta=D,U=A*2^j) realizes
D+A*2^r when r>=j+1. A head family that factors as an already realized
constant times such an affine family is thereby closed uniformly, with
the finite remainder checked separately. Negative defects are permitted;
the final Hamiltonian count is positive because a literal graph is built.

## F6 The99 prime wheel as a costed search policy

The inherited physical family has H_k=99*2^k+23, k>=1, on k+7 vertices,
and recurrence H_(k+1)=2H_k−23. The first25 primes are precisely those
below100. Their exact exponent wheel, from
[the endpoint compiler](../endpoint_coverage_compiler_20261111/CONNECTIONS.md),
has11 active divisibility classes and14 inactive primes. It has period
428,468,040 and226,472,946 surviving residues. These are states of a
divisibility test, not proof of primality or of realizability.

[natural_recursion.py](natural_recursion.py) adds actual construction costs
to this wheel and compares the first30 terms with the inherited complete
least-order bank through11. Two hostile controls decide what a valid
"least resistance" policy must remember:

* H_1=221=13*17: the direct family uses8 vertices; the least-order
  factor witnesses use5+6=11. Factoring first is more expensive.
* H_16=6,488,087=1301*4987: no prime below100 divides it, but joining
  the known factors uses19 vertices, versus23 in the direct family.
  The small-prime wheel can miss a cheaper product route.

The wheel also has infinitely many composite survivors: the inherited
class k=2+8,140,892,760*t, t>=1, survives all25 primes but is divisible
by419. Each number is already realized by the physical ray regardless
of this arithmetic classification.

A rigorous search policy compares available exact graph constructors,
their order costs and remaining endpoint interfaces. It may use the wheel
to prioritize factor branches. A greedy scalar rule or a modular survivor
cannot replace the constructor, and no such rule is proved to terminate
for every target here.

## F7 What a universal proof still needs

The successful mechanism is a division of obligations: renewable unit
seeds for upper rays, nonunit affine interfaces and products for lower
heads, and certified finite exceptions. This improves on requiring one
operation to do all three jobs. It still needs a uniform source-producing
rule for unbounded new coefficients, or a uniform realization theorem for
the prime and7-prime lanes. The finite stable-clone union has only
polylogarithmically many outputs below X, by Session288's proved bound;
even its new cubic family-count polynomial cannot overcome that limitation.

**Named working hypothesis: finite head-interface closure for a recursive
birth grammar.** A useful grammar may have infinitely many birth cores
but only finitely many normalized defect/slope transformations needed
to repair their lower heads. The two families above are positive examples;
the21 factor at a=1 is a necessary exceptional branch. This is OPEN and
does not assert such a grammar already exists. Its cheapest next test is
to compile the complete head obligations of each new birth recurrence and
look for a fixed interface bank; an unrepaired infinite family would
refute a proposed finite bank, not the spectrum conjecture.

**Completed next test in this session:** [FAMILY-COMPILER.md](FAMILY-COMPILER.md)
closes all786 head obligations of all540 binary birth-clock laws through
order10, after adding two auxiliary source graphs outside the unit-bearing
source bank. A further320 higher-mode odd-slope families need no repair.
This is a complete positive test for that library, not a universal grammar.

Connections to Collatz, automata or prime geometry here are exact only
when they retain the affine carry, integer divisibility state, physical
endpoint selector and order cost. The residue automaton and the graph
constructor are different objects. The new results do not prove a Collatz,
Carmichael, Catalan or Zeckendorf completeness claim.
