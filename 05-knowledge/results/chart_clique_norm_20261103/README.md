# Chart boundaries, exact clique carriers and arithmetic collision fibres

**PROVED scoped all-size theorems / FINITE-EXACT controls / CITED inherited
and classical inputs. Same-agent audit with independent computational paths.**
No new Ramsey bound, general tournament-spectrum proof, knot classification,
or global fence optimum is claimed.

The requested themes produced three new deductions and an independent
geometric verification, with explicit maps and operation scopes:

1. **Complete integer cubic-fibre classification.** For every positive integer
   q, collisions of F_q(t)=q^3+1-t^3+3qt correspond to Eisenstein norms of q.
   There are2[q square] two-point fibres and
   sum_(d|q)chi_3(d)-[q square] three-point fibres. For1729 there are exactly
   eight three-point fibres, all explicitly listed and independently counted.
2. **A literal quartic/Frobenius map and a sharp precision inverse.** Scaling
   w=t/sqrt(6q) puts the cubic trace into the inherited quartic discriminant
   parametrization. The nine q=5 states all have two finite quartic preimages;
   the three-point C=0 boundary is a different chart. Modulo9 the nine count
   states collapse along the three older Solfeggio digit-rotation orbits;
   modulo27 they are recovered. The inverse extends to every precision when
   q=2 mod3, with exactly one extra output ternary digit.
3. **Exact all-size clique-mean carrier.** For r>=3 colours and forbidden
   sizes t_i>=3, fixed independent attachments have minimum rational dimension

       1+(r-1)binom(n,2)+sum_i sum_(k=3)^min(n,t_i-1)binom(n,k)
                         +#{i:n>=t_i}.

   Actual attached-clique probes recover every retained proper flag. A unit
   coefficient minor proves integral precision and exact closed quotients
   modulo every modulus. For(3,3,5) the first relevant dimensions are10,20,39;
   these do not replace tournament dimensions7,28,126.
4. **Shared-boundary verification.** At12 unit edges the exact pseudoforest
   maximum is3 for capped total area and3(2+sqrt3) for product area. A square
   grid gives total area4. A two-face theta graph independently realizes
   Friedman's known product(89+26sqrt13)/16, strictly above the pseudoforest
   ceiling. Incoming Session145 independently proves this bound and stronger
   cactus/theta theorems; our replay is credited accordingly. The two
   geometric problems have different constraints.

Read [ARITHMETIC.md](ARITHMETIC.md), [CLIQUE_CARRIER.md](CLIQUE_CARRIER.md),
[FENCES_AND_PORTS.md](FENCES_AND_PORTS.md), and [BRIDGES.md](BRIDGES.md).
The latter recovers the knots and PRIMEGAME work, proves the maximal inert
scaling criterion, and gives a fully guarded infinite1729 execution family.
The pair JK is indeed the unique reciprocal pair. The three quartets have
prime-support widths5,6,7, while their valuation ranks are all4: the claimed
increase is meaningful for support width, not for every complexity measure.

Incoming Session145 is incorporated: its pooled-parity tournament carrier
has exact order-six minimum571, improving the previous upper bound661.
Session146 now also proves minima2668 at seven and13028 at eight;
56296<=d(9)<=65908. All-order sharpness remains open. Incoming prime-wait work already
classifies the weighted inverse-chart rule and inertia partitions; we reuse
it rather than claiming the quartic boundary classification anew.

## Reproduction

From the repository root:

```powershell
python -B -X utf8 05-knowledge/results/chart_clique_norm_20261103/verify.py
python -B -O -X utf8 05-knowledge/results/chart_clique_norm_20261103/verify.py
python -B -X utf8 05-knowledge/results/chart_clique_norm_20261103/invert.py --q 5 --residue 148 --digits 2
python -B -X utf8 05-knowledge/results/chart_clique_norm_20261103/draw.py
```

The inverse example returns trace2 modulo9, using the output modulo27.
Verification requires SymPy1.14.0; the inverse and SVG generator use only
the standard library. See [AUDIT.md](AUDIT.md) and the immutable source pins.

[Shared-boundary diagram](shared_boundaries.svg) is a numerical illustration;
the radical coordinate and unit-length proofs are exact in the companion.

The primary fence sources are the author's
[total-area gallery](https://erich-friedman.github.io/packing/fence/) and
[product gallery](https://erich-friedman.github.io/packing/product/).
They give best known constructions, not certificates for every optimum.
The April2026 [Small Ramsey Numbers survey](https://www.cs.rit.edu/~spr/ElJC/ejcram18.pdf)
records upper bounds57 and91. Incoming Session145's primary-source audit
supersedes the latter with the computer-assisted preprint bound86 from
[Dauth's neighbourhood-counting manuscript](https://lachy-dauth.github.io/Ramsey_Neighbourhood_Counting.pdf).
That newer source audit is inherited here; the current intervals it supports
are45..57 and61..86. Neither Ramsey bound is improved by this packet.
