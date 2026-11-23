# A window endgame for the tournament spectrum

## What the analytic contradiction actually uses

The supplied sketch matches the architecture in
[Slijepčević’s outline, §2, Propositions 2.2–2.5](https://hexagonmath.org/pdf/2611.00812v1).
That document is a human-readable reconstruction, not by itself
the verification. The [independent verification record](https://math.mit.edu/~drew/grh/)
reports kernel and statement-comparison checks of the formal release,
including RH and **Dirichlet** GRH (not all automorphic L-functions).
The logical comparison below uses the window principle, not the
automorphic constructions or any new re-verification.

Fix \(\chi\) and an alleged zero with real part \(\beta>1/2\).
The loss constants \(e_D,e_L\) must be fixed before choosing
\(\varepsilon\). Choose that parameter and a test exponent with
\[
 \tfrac12+(e_D+e_L+7)\varepsilon < \theta < \beta,\qquad
 0<\theta<1 .
\]
A fixed multiplicative constant can be absorbed for large \(x\)
because the exponent inequality is strict. The one-sided upper
bound must hold at **every** point of the same arbitrarily late
windows. The contrapositive of windowed Landau forces a point in
**each** sufficiently late window where \(G(x)>x^\theta\).
One window beyond both cutoffs is enough. Neither an average
over \(x\) nor two unrelated sequences of good and bad points
suffice. The slope average can be over auxiliary lattices at
fixed \(x\), provided the explicit formula transfers it back to
a pointwise bound for the same \(G(x)\). The \(\sqrt{x}\) chain
length is used in obtaining that bound, not in the final
contradiction.

## A missing target needs a persistent signal

An absent tournament count need not itself grow into more absent
counts. In fact 7 is absent while 49 and 343 are witnessed;
21 is absent while 63 is witnessed. Thus counting *wholly
empty target windows* is the wrong analogue. One can instead
propagate an obstruction in a **graded** observable.

Let \(K(N)\) be the maximum number of nonsingleton strongly
connected components in a tournament of height \(N\), and
\(-\infty\) if \(N\) is absent. Such a component contributes a
factor greater than one; the factors multiply over the
ordered SCC decomposition. Let \(b_d\) be the number of
strong isomorphism classes of count \(d>1\).
Each \(b_d\) is finite by [Busch](https://doi.org/10.37236/1141).
We use the known exclusions of 7 and 21 and
[witnesses for 49, 63, and 343](../chanpak_midpoint_ogg_chambers_20261118/README.md#3-ogg-prime-chambers-and-two-residual-leaves).

**Prime lane.** For a fixed odd prime \(p\), if \(p\) occurs,
\(K(p^r)=r\) by taking \(r\) prime-count cores. If \(p\) is missing,
every nontrivial component absorbs at least two copies of \(p\);
hence \(K(p^r)\le\lfloor r/2\rfloor\) for every \(r\).
When both \(p^2,p^3\) occur, this upper bound is sharp for every
\(r\ge2\), using 2s and one 3 for an odd exponent. The example
\(p=7\) realizes exactly this masking: all powers from the
square on occur, but their grades are too small.

**Single-seven lane.** Fix an attained odd prime \(p\ne7\).
Put \(\eta_p=\min\{s\ge1:7p^s\text{ occurs}\}\), or \(+\infty\).
If \(r<\eta_p\), \(7p^r\) is absent. Otherwise
\[
 K(7p^r)=r+1-\eta_p .
\]
Only one component can carry the single 7. Its other factors
consume at least \(\eta_p\) copies of \(p\); the residual
\(p\)-powers supply at most one component per copy. A first
witness of \(7p^{\eta_p}\) is itself a strong core after removing
singletons, since a nontrivial SCC factorization would give a
smaller witness on the same ray. Pad it by \(r-\eta_p\) copies
of a strong \(p\)-core to attain equality. In particular if
\(7p\) is missing the grade is at most \(r-1\) on the *entire*
ray, even if eventual heights all occur. The known 21 gap
sets \(\eta_3=2\), hence \(K(7\cdot3^r)=r-1\) for \(r\ge2\).

Higher seven valuation erases that particular evidence. If
\(m\) is supported on attainable primes other than 7 and
\(e\ge2\), then \(K(7^em)=\Omega(m)+\lfloor e/2\rfloor\).
Indeed a pure-seven component consumes at least two 7s,
whereas every other nontrivial component requires at least one
non-7 prime factor. Use witnesses for 49 and 343 for the
seven part and separate prime cores for \(m\). The 49- and
343-witnesses must have strong nontrivial cores by factorization
and exclusion of 7.

### An exact uniform window condition that would close Chan–Pak

Fix \(\lambda>1\). Any ray \(cp^r\) with fixed \(p,c>0\) and
\(p>1\) meets every sufficiently late interval \([Y,Y^\lambda]\).
For **each fixed** odd prime \(p\ne7\), a lower witness
\(K(p^r)>r/2\) at some ray point in an arbitrarily late window
would force \(p\) to occur. For **each fixed** odd prime
\(p\ne3,7\), with \(p\) now known to occur, a witness
\(K(7p^r)>r-1\) would force \(7p\). It is enough to produce
one such \(r\) per relevant fixed prime; a prime varying with
the window would not be sufficient. This is a finite-gap
threshold like the analytic exponent margin, not a density
claim. In the single-seven lane even \(K/r\to1\) can hide the
missing count, so keep the full integer deficit.

These two uniform achievements would finish the conjecture
by ordinal products. With no factor 7, use the prime cores.
With two or more 7s, pack them into 49s and, if needed, a 343.
With exactly one 7, use a \(7p\) seed for some prime \(p\ne3,7\)
dividing the remainder; if none, the remainder is a power of
3, and 63 works provided its exponent is at least two.
The only exceptions left are 7 and 21. No such uniform
witness supply is proved here.

## The exact finite-difference probe

Let \(a_N(n)\) count unlabelled tournaments on \(n\) vertices
with height \(N\). The
[inherited fixed-count size series](../chanpak_support_birth_grammar_20261122/FIXED-COUNT-SERIES.md)
comes from singleton padding: an ordered strong word of total
order \(m\) and grade \(k\) contributes
\(\binom{n-m+k}{k}\) in its polynomial range. The sharp
individual onset is \(n\ge m-k\); Busch gives \(m-k\le3\log_5 N\).
The signed forward difference \(\Delta f(n)=f(n+1)-f(n)\)
kills grades below \(r\) and reads each grade-\(r\) word
with weight one at such orders.

Consequently, for every \(n\ge\lceil3\log_5 N\rceil\) and all
\(r\ge1\), the two lane probes are
\[
 \Delta^r a_{p^r}(n)=b_p^r,\qquad
 \Delta^r a_{7p^r}(n)=r\, b_{7p}\,b_p^{r-1}\quad(p\ne7),
\]
with the second formula using \(b_p^0=1\). The second maximal
grade has exactly one \(7p\)-core among \(r-1\) \(p\)-cores,
and its position supplies the factor \(r\). The formula is
valid even when \(p\) is absent or a lane has zero witnesses.
In positive cases it reads an exact nonnegative integer
from \(r+1\) consecutive order counts, not the size of the
main height or an aliased residue. A lower estimate strictly
above zero after controlling the signed error would be a
literal analogue of the last analytic step. This is an
alternative recognition task, not an algorithmic saving
over enumerating the small strong core bank in general.

## Why mere coverage estimates cannot substitute for the grade

Take the abstract multiplicative monoid of positive odd
integers with \(7,21,101\) removed. Products of retained
integers cannot give the two omitted primes; the only
proper factorization of 21 uses 7. This monoid has bounded
gaps, density one, and every residue class compatible with
oddness represented at bounded gaps for each fixed modulus.
On any such progression use the Chinese remainder theorem
to impose divisibility by a prime square coprime to the
modulus. It nevertheless has a third hole. This is not a
claim about the actual spectrum; it proves that the
multiplicative closure axiom and these window properties
do not force the conjecture. In a physical example an
[inherited all-modulus seam rail](../chanpak_profinite_automata_20261122/ALL-PRIME-RAIL.md)
meets every modular target test at exponent 61 but never
meets the absolute height.

The parallel with the analytic proof is therefore restricted
but useful: find a statistic on which an actual exception
forces a persistent, quantifiable defect, keep the parameter
of that defect fixed, and obtain a bound that rules it out.
Dirichlet GRH gives no such bound on tournament constructions.

## Companion certificates

[Order-uniform rank read](UNIFORM-RANK.md) gives an alternative
to high finite differences: at one explicitly chosen large order,
the logarithm of the class count approximates its SCC grade to
within \(1/4\), uniformly in the target. It also describes the
single-seven grade for a mixed product of primes.

[Collision energy](COLLISION-ENERGY.md) gives a sharp separate way
for a sufficiently strong aggregate inequality to certify a whole
finite window, with an exact subset-sum evaluator and a physical
five-target example. Neither the energy estimate nor the grade lower
witness has been extended to all targets.

Reproduce the finite controls with python3 window-grade-replay.py.
The six literal seed rows are from the inherited Ogg witness bank.
The grade replay builds from just those chosen seeds; it is not a
census of all strong cores.
