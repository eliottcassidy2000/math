# Exceptional islands, golden units, and the two sorts of missing coordinate

**PROVED** elementary statements B1–B4 below; B3 compares a candidate scalar
tournament language with the inherited proved square-sum order language.
It does not establish tournament completeness. Classical unit arithmetic
and elementary tree counting are not claimed as externally novel results.

## B1. Why deleting two ternary values creates six roots

In an infinite rooted b-ary tree delete a finite set E not containing its
root. Let c(E) be the number of components induced on E, and e(E) its
number of internal parent-child edges. The surviving forest has

    1+b|E|−e(E) = 1+(b−1)|E|+c(E)                       (3)

roots: the old root, and each child of a deleted vertex which is not itself
deleted. The induced graph on E is a forest, so e(E)=|E|−c(E), proving the
formula. Disconnected deletion sets have a different boundary charge; for
example two isolated ternary nodes create seven roots, whereas two adjacent
ones create six. Deleting the original root requires a different formula
and is excluded here.

On positive odd integers the children 3h,3h+2,3h+4 form a full ternary tree
rooted at1. The coordinate u=(h−1)/2 changes them to 3u+1,3u+2,3u+3;
the unique parent strictly decreases u. Since21 is a child of7, removing
{7,21} leaves exactly

    {1,23,25,63,65,67}.                                  (4)

This is the exact mechanism behind the annotation. It is an arithmetic
address tree. Although each root has an actual tournament witness, legal
source-preserving realization of every ternary branch and its return is
still required to turn (4) into a spectrum proof.

## B2. A special junction between this ternary chain and binary repunits

For positive integers a,b,

    3^a(2^b−1)=2^(2b)−1
       iff (a,b)=(1,1) or (2,3).                         (5)

Cancel2^b−1 to obtain3^a=2^b+1. For b=1 the solution is a=1; b=2 gives5,
not a power of3. For b>=3 reduction modulo8 forces a even. The two even
factors3^(a/2)−1 and3^(a/2)+1 multiply to2^b and differ by2. They must be2
and4, giving a=2,b=3. This elementary proof does not invoke Catalan's theorem.

Thus7→21→63 is the unique instance with b>1 in which a pure ternary chain
starting at a binary Mersenne number reaches the Mersenne number with doubled
exponent. This does not say it is the only relation between powers of2 and3.
Here63=7*9 records both repeated prime support and new valuation. The next
doubling,4095/63=65, leaves the pure3-power chain. That factor65 is exactly
the q+1 companion of q−1=63 at q=64, as in the projective chart of
[PORT_MEMORY](PORT_MEMORY.md), though the chart's operations are different.

## B3. The two exceptional-island patterns have different predicates

Use step2 adjacency on the **candidate** tournament scalar language

    A = {positive odd h} minus {7,21}.

Its components are {1,3,5}, {9,11,13,15,17,19}, and the odd tail beginning23.
The full equality of A with the actual raw Hamilton-path spectrum is OPEN.
The exclusion of7,21 is already proved at every order.

For the ordinary square-sum graph on[1,n], existing constructions and the
small-order checks give the **proved** path-order set

    B = {1,15,16,17,23} union {n>=25}.

Recover the strongest result through
[the archive recovery](../square_sum_archive_stability_20261102/RECOVERY.md),
not by reopening ordinary Hamiltonicity as a conjecture. By step1 adjacency
B has components {1}, {15,16,17}, {23}, and the tail beginning25. The
terminal-path and fixed-endpoint versions have different holes; they must
not be substituted for this predicate. Square-sum cycle existence begins32.

There is a unique order-preserving rank bijection A→B. For h>=23 it is

    n=(h+35)/2,                                         (6)

because rank_A(h)=(h−3)/2, while rank_B(n)=n−19. The constant35 is the
arithmetic consequence of the finite missing prefixes, not a discovered
geometric invariant. The six roots (4) map to

    {1,29,30,49,50,51}.

On the tail, (6) transports a ternary step to

    h→3h+2j   corresponds to   n→3n−35+j, j=0,1,2.       (7)

This is an exact conjugacy of scalar arithmetic on the tails, preserving
rank and branch names. It supplies no map of tournament paths to square-sum
paths, no port transport and no Gerbicz extension at(7). The appearance of
49 is therefore explained without identifying it with the existing base49
square-sum controller. Testing a path-level realization of(7), with named
endpoints retained, is a stronger question than observing the integer49.

## B4. A setting where two bits really do suffice—and where they do not

Let phi=(1+sqrt5)/2. Every unit of Z[phi] has a unique form

    u=sigma*phi^k,  sigma in {−1,+1}, k in Z.              (8)

Here is an elementary proof to fix the scope. An algebraic integer
a+b phi is a unit iff its norm a²+ab−b² is±1. Make a unit positive and
multiply by a suitable phi-power to place v in[1,phi). If its norm is+1,
b=(v−1/v)/sqrt5 lies in[0,1), forcing b=0 and v=1. If its norm is−1,
b=(v+1/v)/sqrt5 lies strictly between0 and1, impossible. This proves(8);
strict positivity and phi>1 give uniqueness.

Observe the integer m=|k|. At every m>0 there are exactly four units,
distinguished by sign sigma and exponent direction tau=sign(k). At m=0
there are just two units±1; take tau=+1 by convention. Multiplication is
implemented by addition of the signed clocks:

    (sigma,k)*(rho,l)=(sigma*rho,k+l).

Conjugation and parity are coupled:

    conjugate(sigma*phi^k)=sigma*(-1)^k*phi^(−k),
    N(u)=(-1)^k=(-1)^m.                                 (9)

So for m>0 conjugation reverses direction and flips sign exactly on odd
levels. This gives the user's positive/negative, even/odd, and growth/decay
duality in one exact carrier. The two-bit fibre collapses at the boundary0.
The fixed integer coordinate is the exponent clock, not the norm.

If instead one observes only |N(u)|, every unit lies above1. Even the
positive units phi^(2j), j>=0, give an infinite fibre; no fixed number of
bits can reconstruct it. This is the missing-clock version of the problem.
Also, addition exits the unit carrier: phi+phi^(-1)=sqrt5 has norm−5.
Thus closure under multiplication does not imply closure under addition.

All units with k nonzero are real quadratic irrationals, whereas±1 are
rational. This family contains no transcendental numbers. Rational/
irrational and algebraic/transcendental are arithmetic classifications,
not involutions automatically paired like sign. The existing
[THM-4088 arithmetic-type blindness](../../../01-canon/theorems/THM-4088-order-tournament-arithmetic-type-blindness-and-lrc-margin-density.md)
proves that their order-tournament alone cannot tell those types apart;
heights and quantitative approximation data can carry what order discards.

## Exact controls

`verify.py` checks(3) on4096 deletion sets, (5) through a,b<=40, unit
coordinates for all norm±1 pairs |a|,|b|<=100, signed clocks−50..50, and the
rank/branch identities on thousands of inputs. `audit.py` recomputes the
forest components by graph search. These are controls of the displayed
proofs, not evidence for either unresolved universal return problem.
