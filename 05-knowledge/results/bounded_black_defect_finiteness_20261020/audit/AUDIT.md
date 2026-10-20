# Independent audit: bounded black defect forces bounded strong order

Date: 2026-10-20. Verdict: **APPROVED** for the all-strong finiteness
theorem and its stated leaf-series, recurrence, asymptotic, and strong-core
degree consequences. This report is separate from the earlier paired-black
leaf-reduction audit. 

## 1. Precise scope

For ordinary, not converse-merged, Hamiltonian-path orbits, put
`h=H/|Aut|`, let `b` count blue orbits, and let `Delta=h-b`. The action on
ordered path words is free and preserves blue status, so the number of
black words is exactly `Delta |Aut|`.

This remains valid when the tournament is not self-converse, in which case
`b=0`. Thus the new theorem really covers all strongly connected tournaments,
not only self-converse ones. The singleton is handled separately; order two
has no strong tournament. Removing strongness invalidates finiteness:
transitive tournaments of arbitrary order have defect zero.

The checked manuscript is
`../BOUNDED_BLACK_DEFECT_FINITENESS.md`,
with its companion `SOURCE_AUDIT.md`.

## 2. Cycle quotient and restored minimum degree

For a cycle successor `R`, the group `G=Aut(T) intersect <R>` is cyclic,
has odd order `m`, and equals `<R^q>` with `q=n/m`. There are exactly `q`
distinct global Hamiltonian-path orbits among its cuts, since any
coordinatewise map between two cut words must be the corresponding power
of `R`.

The number of blue cut-orbits is zero or `gcd(2,q)`, using
`r_i r_a=R^(2(i-a))`. Therefore `q<=Delta+2` for every cycle; in the NS
case the stronger `q<=Delta` holds.

The new degree step is sound. The `q` vertex orbits of `G` all have size
`m`. Each induced tournament is vertex-transitive under `G`, hence regular
with internal outdegree `(m-1)/2`. It follows that

    minimum outdegree >= n/[2(Delta+2)] - 1/2.

The proof does not accidentally assume global regularity or transitivity.
Strongness supplies a Hamiltonian cycle by Camion, which is the only
existence input needed for this step.

## 3. Cycle-count upper bound and exact square-root fibers

Each closing Hamiltonian word determines exactly one directed cycle up to
rotation. Selecting one black cut word from each cycle having a black cut
is injective. Hence at most `Delta |Aut|` cycles have a black cut.
Nonclosing black words only enlarge this upper bound.

If every cut of a cycle is blue, adjacent reversal maps yield
`R^2 in Aut`. This is only a necessary condition. The converse must not
be used, particularly for even order.

The full-cycle qualification in the square-root lemma is essential and is
present throughout. For fixed `g`:

- At odd `n`, there is exactly one full `n`-cycle square root when `g`
  is an `n`-cycle, namely `g^((n+1)/2)`, and none otherwise.
- At `n=2m`, there are exactly `m` such roots when `g` has two `m`-cycles,
  and none otherwise. After fixing `x` in one cycle and choosing its image
  `y` in the other, commutation with `g` forces the entire interlacing:
  `R(g^j x)=g^j y` and `R(g^j y)=g^(j+1)x`.

Thus the manuscript's deliberately weaker bound of `n` roots per `g`
is correct. It gives

    HC(T) <= (Delta+n)|Aut(T)|.

This counts actual cycles, not automorphism-orbits of cycles. At even order
divisible by four, there can be no all-blue cycle because `R^2` would have
even order. No converse implication is needed anywhere.

## 4. Published inputs independently checked

I independently queried the primary publisher records.

- Thomassen's primary indexed publisher abstract gives the bound
  `4^k k!` under strongness and minimum outdegree at least `3k+3`.
  The hypothesis is not restricted to regular tournaments despite the
  article title. Source:
  [Hamilton Circuits in Regular Tournaments](https://www.sciencedirect.com/science/article/pii/S0304020808730056),
  DOI [10.1016/S0304-0208(08)73005-6](https://doi.org/10.1016/S0304-0208(08)73005-6).
  Direct publisher-page retrieval returned HTTP 403; the full proof was not
  retrieved. The manuscript accurately labels this a cited theorem verified
  from the primary abstract, not an independently reproved theorem.
- Dixon's primary three-page PDF was opened successfully. Theorem 1 on
  page 503 gives `|Aut(T)|<=3^((n-1)/2)` for arbitrary tournaments, without
  regularity assumptions. The use of Feit-Thompson is explicitly stated.
  Source: [The Maximum Order of the Group of a Tournament](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/682531E83B0A3AA042CFC5852474E914/S0008439500067461a.pdf/div-class-title-the-maximum-order-of-the-group-of-a-tournament-div.pdf),
  DOI [10.4153/CMB-1967-048-9](https://doi.org/10.4153/CMB-1967-048-9).

The cited objects are distinct Hamiltonian circuits, with the ordinary
unrooted convention. The manuscript explicitly fixes its count as cycles
modulo cyclic starting point. Even if a weaker interpretation introduced
an extra factor `n`, the finiteness proof would survive after including
that factor; no exponential-versus-factorial conclusion hinges on it.

## 5. Finiteness comparison

For fixed `d>=0`, let `s=d+2` and

    k=floor((n-7s)/(6s)).

For all sufficiently large `n`, `k>=1` and

    3k+3 <= n/(2s)-1/2 <= minimum outdegree.

The floor and constant seven are correct. Thomassen and Dixon then give

    4^k k! <= HC(T) <= (n+d)3^((n-1)/2).

Since `k=n/(6s)+O(1)`, the logarithm of the left side divided by `n`
tends to infinity; the right-side logarithm divided by `n` tends to
`log(3)/2`. The elementary estimate `k! >= (k/2)^(k/2)` used in the
manuscript is valid for every integer `k>=2`. This supplies a uniform
order bound for all tournaments with `Delta<=d`, proving finiteness of
isomorphism classes.

No explicit cutoff, sharp threshold, or classification is claimed.
The argument is effective in principle because the numerical inequalities
are explicit, but the packet does not compute such a cutoff.

## 6. Leaf consequences

Taking `d=2` makes the strong SC center series `D(z)` a finite polynomial.
Its coefficients and degree remain unknown. In particular this does not
prove that S4 is the only center.

The previously audited exact leaf formula now becomes rational:

    Leaf(z) = 1/[2(1-z-z^3)]
              + [D(z)-(1+z+z^3)/2]/(1-z^2-z^6).

The denominator divides the product of the two displayed denominators.
The product expansion and all recurrence signs in the manuscript are
correct:

    (1-z-z^3)(1-z^2-z^6)
      =1-z-z^2+z^5-z^6+z^7+z^9.

The recurrence is eventual; no onset is supplied while `deg D` is unknown.
Calling this product a known annihilating denominator is justified, but
claiming it is necessarily the minimal reduced denominator would not be.

The asymptotic proof avoids an unnecessary analytic assumption. For the
fixed finite correction polynomial, coefficient comparison gives
`|l_n-a_n/2| <= C a_floor(n/2)`. Concatenation at the prescribed weight
boundary is injective, so
`a_n >= a_floor(n/2) a_ceil(n/2)`. The weight-three block choices prove
`a_(3r)>=2^r`; monotonicity then forces the denominator factor to infinity.
Consequently `l_n/a_n -> 1/2`.

This is a leaf-vertex asymptotic, not a bridge asymptotic.

For strong-core quotient vertices, the manuscript uses the correct two
degree conventions: SC degree is `Delta/2`, while NS merged degree is
`Delta`. Therefore bounded quotient degree implies bounded defect and
bounded strong order. No topological conclusion about connectivity or
bridges follows.

## 7. Fresh finite controls and existing exhaustive evidence

I inspected the independent cycle verifier and its completed exhaustive
order-six output. It checks all 22,890 nontrivial labelled strong
tournaments through six, including 16,800 NS examples at order six,
62,234 Hamiltonian cycles, and 372,582 cuts. Its group, orbit, path,
degree, and square-root assertions implement the proved statements,
not heuristic estimates. I did not repeat that full larger computation.

I separately wrote `check_bounded_defect.py`, importing no producer code.
It freshly checks:

- Every labelled strong tournament at orders three, four, five:
  `2,24,544` tournaments and `2,24,768` cycles respectively.
- Twelve NS order-six controls, one per reported `(H,Aut,HC)` profile.
  The existing census supplies only their masks; all paths, blue words,
  groups, cycles, and inequalities are recomputed.
- Every full `n`-cycle for `2<=n<=9`, grouped by its square, giving fiber
  sizes `1,1,2,1,3,1,4,1`.

All assertions passed. This checks the wider NS scope without introducing
a large new census. Computation validates the mechanism but is not a
dependency of the all-order proof.

One initial README link pointed to `independent/REVIEW.md` instead of
`independent/INDEPENDENT_REVIEW.md`; this was reported and corrected.

## 8. Reproduction and artifact scope

From this publication directory:

    python3 audit/check_bounded_defect.py

The script accepts `--source path/to/verification_results.json.gz` for the
portable compressed census. It uses that input only to select NS6 masks.
Fresh output is in `verification.json` and `verification.out`.
Audited source hashes are in `AUDITED_SOURCE_SHA256SUMS`.

The bounded-defect theorem upgrades the earlier leaf packet's infinite
residual to a finite but unclassified center set. That earlier scoped
reduction remains correct; its old frontier should not be presented as
the latest state after this theorem is adopted.

Portable packaging note: the verifier default source was changed only to the included `../independent/verification_results.json.gz`, and the portable command was rerun successfully. Audited-source checksum paths are relative to the publication directory.
