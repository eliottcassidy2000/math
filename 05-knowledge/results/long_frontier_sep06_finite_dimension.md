# A finite-length obstruction to the asymptotic coefficient bound

**Status: REFUTED by FINITE-EXACT actual normalized lists; INDEPENDENTLY
AUDITED.** Neither the
all-N inequality

    N(R-K3)>=C

nor finite-N optimality of the three-equal-positive/equal-negative family
follows from the [sharp asymptotic theorem](long_frontier_sep06_dimension.md).
Both claims fail already at total lengths four and five. The asymptotic
formula I_N=K3+C/N+o(1/N) is unchanged.

All lists below have p1=p2=1, E>0, and the actual ordered top-two distance.
The quantities R,K3,C have exactly the definitions of the asymptotic note.
In particular C=2.1722010964723645447... and K3=0.35013450126998956116....

## 1. Inheritance and the first failed implication

The nearest proved mechanism is the length-uniform three-atom expansion in
the asymptotic note. Its positive coefficients control lists already near
the global three-atom minimizer. The canonical hostile is the other
zero-energy boundary: one positive atom tends to one and all others to
zero. The [minimizing-sequence classification](long_frontier_sep06_minimizers.md)
proves there

    R -> L=sqrt(2)-2/3 > K3.

This fixed gap excludes the boundary from global asymptotic minimization.
It does not exclude it from winning at small fixed lengths. In fact

    5(L-K3)<C<6(L-K3).                              (1)

These are exact strict radical comparisons, checked with rational bounds.
Thus the boundary itself predicts counterexamples at N=4 and N=5.

The live board is the three-atom limit, the one-atom boundary, total length,
ordered top-two distance, rational normalization curves, and finite versus
asymptotic optimization. The corrected near miss is extrapolating a local
first-order expansion to every point at every fixed N. The missing sidecar
is the competing boundary value. The source-to-target map retains the
actual normalized list along a rational curve; its one-atom limit has a
controlled limiting quotient but loses the strict positive-energy condition.
R is undefined at the endpoint itself. The finite rational controls below
retain positive energy and an actual defined quotient.

## 2. Two explicit nonzero rational counterexamples

At N=4 take

    r=(256,16,1,-16)/257.

The numerator sum is 257 and its square sum is 257^2. All four coordinates
are nonzero, p4<1, and the leading positive coordinates are 256/257 and
16/257. Direct exact evaluation with rational square-root enclosures gives

    21419/10000 < 4(R-K3) <21420/10000,
    21722/10000 < C <21723/10000.                    (2)

Thus 4(R-K3)<C, with no limiting or rounding ambiguity.

At N=5 take

    r=(675,30,1,-15,-15)/676.

Again the numerator sum is 676 and its square sum is 676^2. All five
coordinates are nonzero; the leading positives are 675/676 and 30/676.
The exact certificate is

    21649/10000 < 5(R-K3) <21650/10000 <21722/10000<C.  (3)

Repeated negative coordinates are permitted by the real-list domain.
These examples also refute the bound on lists of length at most N, without
using zero padding or a list of smaller actual length.

The three-equal-positive family has quotient approximately 1.88561808316
at N=4 and 1.16378989075 at N=5. Each actual rational counterexample has
strictly smaller R, verified with disjoint rational intervals. Hence that
family is not a finite-N optimizer in these dimensions. No claim about the
true finite-N optimizer or attainment of its infimum is needed.

## 3. A rational curve explains the obstruction

For any integer N>=4 and k>=1 put

    T=(N-3)(N-2)/2,
    r_N(k)=(T*k^2, (N-3)k, 1, -k,...,-k)/(T*k^2+1),  (4)

where -k occurs exactly N-3 times. This is an actual length-N list. Its
numerator sum is T*k^2+1, and the square sum is

    T^2*k^4+[(N-3)^2+(N-3)]k^2+1
       =T^2*k^4+2T*k^2+1=(T*k^2+1)^2.

Thus p1=p2=1 identically. Every coordinate is nonzero, there are three
positive coordinates, the displayed first two are largest, and E>0.
For fixed N, k tending to infinity gives the one-atom limit and hence R->L.
Equation (1) then proves that each of N=4,5 has infinitely many rational
counterexamples of the form (4). The explicit rows in Section 2 use k=16
and k=15, respectively.

The verifier checks every preceding integer parameter k in these two
declared curves. They are the first violating integer parameters there.
This is not a minimality assertion among all eligible rational lists.

## 4. What the bounded controls do and do not support

For each N=4,...,20, retain precisely two controls: r_N(100) and the
three-equal-positive/equal-negative family of total length N. Exact interval
comparisons give:

* r_N(100) violates N(R-K3)>=C precisely for N=4,5 in this control set.
* The equal-three control satisfies N(R-K3)>C at all these dimensions.
* r_N(100) beats the equal-three control at N=4,5,6,7. The equal-three
  control beats r_N(100) at N=8,...,20.

Consequently finite-N optimality of the equal-three construction is also
excluded at N=6,7 by actual nonzero length-N controls. These two-family
comparisons are not a minimizer census, and they do not locate a global
transition or assert optimality of either family in the surviving dimensions.

The strongest established survivor is the sharp asymptotic theorem, with
its exact first-order equality conditions. A precise remaining question is

    Does N(R-K3)>=C hold for every eligible list of length at most N
    whenever N>=6?                                  (OPEN)

The bounded controls do not refute it. No evidence here proves it either.
Any proof would need a global treatment of the one-atom and two-atom
boundaries as well as the local three-atom response. The cheap question
asked in this pass has a decisive counterexample, so no wider optimization
search is used to conceal the missing global inequality.

## 5. Reproduction and exact universe

[Source](../../04-computation/long_frontier_sep06_finite_dimension.py),
[output](long_frontier_sep06_finite_dimension.out). The source uses only
standard-library integers and `Fraction`. Square roots are enclosed using
integer square roots at forty decimal places; every enclosure is checked
by squaring its rational endpoints. Subsequent interval operations round
outward. Every claimed strict comparison has disjoint rational bounds.

The explicit universe is the parameter head k=1,...,16 at N=4, the head
k=1,...,15 at N=5, and the two fixed controls at each N=4,...,20. Exact
normalization and actual ordered top-two labels are checked for each
rational row. The equal-three controls use their algebraic normalization
formula from the proved asymptotic note. No floating optimizer supplies a
proof input.

```
python3 -B 04-computation/long_frontier_sep06_finite_dimension.py
python3 -B -O 04-computation/long_frontier_sep06_finite_dimension.py
```

Normal and optimized outputs are byte-identical, with 445 explicit gates.

* Source SHA256: `bc3ba159919527741e8b7439823759f13575916549ecdb666b9667f9fb5054f1`.
* Output SHA256: `490e74d700b33019a726f75e78b3672aff13b1e3ecc4737cfe22f83583d9d024`.
* Semantic digest: `1dc1ef2b4c92cf87562a81c1248f567b290475d797d15a5569e83f27b422cf6c`.

The [independent referee](long_frontier_sep06_finite_dimension_audit.md)
checked the complete proof, normalized rational curve, ordered coordinates,
boundary mechanism, and the scope of every finite comparison. Its separate
ordinary polynomial and squared-coefficient calculation recovered
D/E=66305/65793 and 444131/456301 for the two witnesses and confirmed
both exact radical inequalities. The 445-gate normal, optimized, and frozen
outputs and all hashes were independently replayed: PASS. All three primary
artifacts are frozen.
