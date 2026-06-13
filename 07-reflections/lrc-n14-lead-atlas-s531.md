# n=14 LRC lead atlas after repo and web scan (S531)

The outside literature moved the frontier.  In the web convention, `k` is the
number of moving speeds and `k+1` is the total runner denominator.  Our repo's
`n=14` case is therefore `k=13`: thirteen moving speeds plus the stationary
observer.

As of 2026-06-01, the relevant public chain is:

- Rosenfeld proved the case of eight total runners by combining computer
  verification with finite-checking bounds.
- Trakulthongchai refined the sieve and proved nine and ten total runners.
- Sungkawichai-Trakulthongchai pushed the same line to eleven, twelve, and
  thirteen total runners, i.e. `k in {10,11,12}`.
- Thus the repo's `n=14` is exactly the next case after the current arXiv
  frontier.

Web sources scanned:

- https://arxiv.org/abs/2509.14111
- https://arxiv.org/abs/2511.22427
- https://arxiv.org/abs/2604.23906
- https://arxiv.org/abs/2605.27941
- https://arxiv.org/abs/2511.16636
- https://arxiv.org/abs/2603.24784
- https://www.cambridge.org/core/services/aop-cambridge-core/content/view/A51A991DE89B8C9C2E2FF13FBD4501DA/S2050509425101072a.pdf/linearly_exponential_checking_is_enough_for_the_lonely_runner_conjecture_and_some_of_its_variants.pdf
- https://www.quantamagazine.org/new-strides-made-on-deceptively-simple-lonely-runner-problem-20260306/
- https://www.sciencedirect.com/science/article/pii/S1574013725000747
- https://arxiv.org/abs/2111.13688

Repo sources leaned on:

- HYP-1990, HYP-1991, HYP-1995, HYP-1998, HYP-1999, HYP-2001, HYP-2003,
  HYP-2007, HYP-2008, both HYP-2009 files, HYP-2010, HYP-2011.
- THM-363, THM-364, THM-375, THM-376, THM-377, THM-380, THM-381, THM-384,
  THM-385, THM-387.
- `lrc_n14_wacky_attempts_s514.py`, `lrc_n14_permutohedron_s526.py`,
  `lrc_n14_permutohedron_covering_s525.py`, `lrc_n14_metagraph_proof_s524.py`,
  `lonely_runner_scalar_excision_s371.py`, `lonely_runner_k13_scalar_gauge_s367.py`,
  `lrc_polygon_outside_inside_s529.py`, `source_sink_apex_arc_s530.py`.

Rebase note: after the first draft, three concurrent S531 commits landed on
main: recursive apex modularity (HYP-2010), the resonance-debt conjecture
(HYP-2009-resonance-debt-conjecture), and the n=4 parity inside-debt law
(HYP-2011).  Leads 75-87 fold those signals back into the n=14 queue.

## Web-method leads

1. Port the Rosenfeld/Trakulthongchai finite-checking engine to `k=13` and run
   small primes first.  The target is not the full proof at once; it is to see
   whether `J(13,p)` empties after the same lift/projection operations.

2. Exploit that `k+1=14` is composite.  The 2026 paper says prime `k+1` forced
   expensive `c=k+1` lifts, while composite `k+1` can use smaller factors.  For
   `14`, try the diagram `1 -> 2 -> 14` and `1 -> 7 -> 14`, then intersect.

3. Add backward projection to the repo's scalar-gauge quotient.  The web sieve
   keeps only an upper bound on `J(k,p)` under projection; the repo already has
   quotient-normalized scalar cells.  These are the same kind of state space.

4. Compute `I(13,p,1)` survival profiles for Goldilocks primes.  For `k=9`,
   runtime had a middle-prime sweet spot.  For `k=13`, do not start at tiny
   primes or huge primes blindly; sample prime windows and plot survivor counts.

5. Use the Sungkawichai-Trakulthongchai polynomial lemma on the odd core.  It
   directly assumes odd prime `k+1`, so it does not apply to 14, but the mod-7
   CRT quotient has six pairs plus a singleton.  A mod-7 analogue may kill the
   `(1,2,...,6)` core remaining after parity split.

6. Replace the prime `one-to-k` tight tuple with the repo scalar-ramp orbit.
   For `k+1=14`, the tight obstruction is not one residue tuple; THM-363/364
   say scalar ramps are a gauge orbit.  The web eventual-properness proof
   should be recast as "every non-scalar quotient class is eventually proper."

7. Convert the repo's `56` missed-cell scalar-puncture moat into a finite-sieve
   pruning rule.  If a survivor is a one-coordinate half-turn defect, the
   missed cells are explicit positive-margin witnesses; the sieve can discard
   that whole equivalence family.

8. Try powers and composites, not just primes, in the dividing-prime product
   argument.  Rosenfeld noted powers did not help much for `k=7`, but `n=14`
   has strong 2- and 7-adic structure.  Test `p=2^a`, `7*odd`, and `14*odd`
   as divisibility certificates.

9. Use the 2026 bound `product(u_i) < B_k` as a scoreboard for partial progress.
   Every verified `J(13,p)=empty` prime raises a lower product bound.  The
   immediate deliverable is a live "percent of log B_13 covered" table.

10. Search the public code ecosystem for Rosenfeld/Trakulthongchai code and
    adapt data structures rather than reimplementing from scratch.  The first
    repo task is a compatibility wrapper translating their proper/improper
    tuples into this repo's endpoint and microstaircase ledgers.

11. Use the finite-sieve viewpoint to validate HYP-1991.  If all verified
    proper witnesses for `k=13` appear at denominators `14*s` with small `s`,
    it gives an external route to the bounded small-cofactor ansatz.

12. Make "eventual properness" a tournament observable.  Vertices are survivor
    classes modulo `p`; edge `A -> B` if `A` has lower lift survivor burden than
    `B`.  The tie path is scalar-distance then support size.  This gives the
    required Tournament Analysis for the web-sieve port.

## Covering and wall-extremality leads

13. Prove HYP-2003: wall-only implies AP up to scaling.  This remains the
    cleanest n=14 bite: every non-AP set found so far has an open lonely alcove.

14. Strengthen "wall-only" to "no positive-margin microcell."  The S367/S371
    cell system already sees scalar ramps as the only full blockers and the
    coordinate-6 half-turn as the unique best non-scalar near-blocker.

15. Prove the one-coordinate scalar-puncture moat symbolically.  The key
   statement: every scalar ramp with one defect exposes at least the `56=7*8`
   positive-margin cells uniquely blocked by that coordinate.

16. Extend the moat to radius two and beyond with a submodular repair-deficit
    inequality.  Current data: radius one opens at least 56 cells, radius two
    at least 112, and non-reverting repairs create many new misses.

17. Turn missed-cell stencils into a SAT certificate.  The 56 cells are seven
    odd shifts times eight alpha stencils.  Encode "every nonzero normalized
    vector misses one stencil" as a small exact Boolean proof.

18. Use scalar-distance as branch-and-bound cost.  Any branch far from scalar
    should be easier; branches close to scalar are handled by moat lemmas.

19. Classify all possible wall-only equality patterns.  Equality at `1/14`
    should force a regular outside necklace and six hidden diameter ties.

20. Search for sporadic wall-only systems beyond AP.  Literature has sporadic
    tight systems at lower total runner counts; n=14 needs an explicit bounded
    search for non-AP wall-only rows before leaning too hard on AP uniqueness.

21. Use boundary compactification as a finite list of tie paths.  The AP wall
    has six diameter ties; enumerate all tie-resolved source classes with six
    ties and test whether any is realizable by non-AP speeds.

22. Prove "positive open source if any diameter tie is broken unevenly."  In
    HYP-2009 language: disturbing the hidden diameter channel should force the
    outside clasp to open.

## CRT and parity leads

23. Recenter the proof on `14=2*7`: six mod-7 pairs plus singleton `7`.  LRC is
    simultaneous safety of the seven CRT classes.

24. Bound pair-class correlations exactly.  Each nonzero mod-7 class has two
    runners `{a,a+7}`.  Use Jensen's mixed-threshold two-function formula to
    compute exact unequal-threshold class-safe intersections.

25. Upgrade the S524 independence heuristic to a rigorous second-moment bound.
    If the seven class-safe indicators cannot be anti-correlated enough, their
    intersection is nonempty.

26. Use inclusion-exclusion by CRT class, not by runner.  Runner-level
    inclusion-exclusion has 13 arcs; class-level has seven objects with strong
    internal structure.

27. Prove that no primitive system can make all seven class-unsafe sets tile
    the circle.  This is the CRT version of wall-only implies AP.

28. Split odds and evens.  One side has at most six speeds, so the proven
    `LRC@7` result gives witness intervals at threshold `1/7`; the problem is
    only coupling the odd witness set with the doubled even witness set.

29. Use leave-one-out `LRC@13`.  Since 13 total runners are now claimed by
    the 2026 preprint, every 12-speed subset of the 13 moving speeds has a
    `1/13` witness.  The slack from `1/13` down to `1/14` may leave room for
    the omitted runner.

30. Build a 13-subset witness union.  For each omitted runner `i`, compute the
    known/proved witness region for the other 12 at threshold `1/13`; test
    whether the union must intersect the `i`-safe set at threshold `1/14`.

31. Formulate a mixed-threshold induction step:
    twelve runners at threshold `1/13` plus one runner at `1/14` implies all
    thirteen at `1/14`.  This is tailor-made for the 2026 mixed-threshold paper.

32. Use mod-7 polynomial interpolation on the AP residue class.  The prime
    lemma for odd `k+1` says the AP tuple becomes proper for large `p`; for
    n=14, ask for an analogous statement after quotienting by parity.

33. Test whether all hard n=14 survivors project to the scalar-puncture class
    in mod 7.  If yes, the CRT proof reduces to the 56-cell moat.

34. Prove a no-gate/gate dichotomy.  Repo evidence says no `14`-multiple keeps
    unit walls alive, while a `14`-gate exports debt.  This is exactly a sieve
    branch split.

## Gap-race and apex leads

35. Make THM-387 quantitative for n=14.  At every wrap reset, calculate the
    exact right-win inequality.  A proof needs one favorable LS reset.

36. Replace the open race by the HYP-1995 compactified wall ledger.  Initial
    segments are wall-only through many `n`; n=14 may be best attacked by
    forcing a closed LL wall first, then perturbing.

37. Track the dynamic apex arc as an interval exchange.  HYP-2008 says the
    source-sink/apex arc is the observer loneliness gap.  The next proof target
    is: the apex clearance cannot stay below `1/14` for a full lap.

38. Split apex orientation into transitive and wrap branches.  Transitive
    branch should be easy by semicircle logic; wrap branch contains AP-hard
    behavior and should carry the six diameter ties.

39. Use the apex as the lowest harmonic in the resonance expansion.  If the
    apex shell is positive, only higher inside shells can destroy loneliness.
    Bound the inside debt relative to apex credit.

40. Pair gap-race resets with CRT classes.  A reset by a runner in class `a`
    changes the singleton or pair-safe state in a controlled way; the quotient
    race may be only seven-dimensional.

41. Search for a resetter-residue guarantee.  HYP-1995 predicts resetter
    residue plus left-runway source controls wall hits; for n=14 residues
    modulo 7 should force one favorable resetter.

42. Use Rifford-style time bounds locally.  If each parity half gets lonely
    within a bounded number of slowest-runner rounds, the coupling search
    becomes finite over a bounded time window.

43. Construct a "two auctions" proof.  Repo data says left and right gap bids
    are anti-correlated; formalize a no-perfect-anti-correlation lemma for
    thirteen speeds.

## Permutohedron, tournament, and tiling leads

44. Continue HYP-2001 with labelled Hall/Farkas certificates.  Raw handoff
    acyclicity is false; after endpoint-private leaves are peeled, no labelled
    source-avoiding circulation should remain.

45. Translate the web finite-sieve survivors into permutohedral chamber words.
    Maybe the web "improper" classes are exactly the repo's source-avoiding
    chamber circulations.

46. Count the n=14 open Ferrers source menu.  HYP-1999 reaches n=9; computing
    n=14 target signatures would give a concrete finite target inside the
    huge A000568 ambient space.

47. Count the n=14 boundary-compactified source menu.  The open body is round;
    the boundary seam is where the AP wall lives.  This count may be small
    enough to certify wall-only implies AP by target-class exhaustion.

48. Use the 2*Fibonacci circular menu as an automaton, but impose edge-word
    arithmetic.  HYP-1997 showed the raw metagraph is non-reducing; the missing
    constraint is the speed-determined order of wall crossings.

49. Express source avoidance as a forbidden word in the cyclic permutohedron.
    Then try automata/transfer-matrix methods with labels `mod 2`, `mod 7`,
    endpoint depth, and apex orientation.

50. Use tiling recursion rather than free arcs.  HYP-2000 says arc flips are
    dependent.  The n=14 source target may be a monotone statement in the
    recursive sub-ranking tree.

51. Treat the six diameter ties as six special tiles.  Prove any closed
    source-avoiding walk must toggle one diameter tile without paying its
    endpoint debt.

52. Use H as a filter, not a scalar proof.  H plateaus on hard ladders, but it
    can still identify which chamber words are regular-polygon-wall-like.

53. Search for a feedback vertex set after quotienting by realizable edge-word
    constraints.  The raw FVS fails; the arithmetic subshift might restore it.

54. Build a "source-neighborhood" tournament on observer-score layers.  THM-385
    makes observer indegree the blocker count; edges compare which layer has
    more forced repairs.

## Harmonic, zonotope, and analytic leads

55. Push HYP-2007's resonance-order expansion at n=14.  The main term is
    positive; AP equality means inside debt exactly cancels it.  Bound all
    non-AP higher-resonance debt away from full cancellation.

56. Split the harmonic sum by dihedral chord-length channels.  HYP-2009 gives
    n=14 channels `12,12,12,12,12,12,6`; the last channel is diameter ties.
    Treat that last channel separately.

57. Try Kloosterman/Gauss bounds per CRT shell.  Prime cases have square-root
    cancellation; n=14 fails only because shells share moduli.  Use product
    modulus `2*7` estimates rather than generic composite bounds.

58. Use Bedert's Riesz-product improvement as a finite-k majorant.  The global
    bound is weaker than LRC, but for `k=13` a tailored Riesz product might
    certify the missing positive margin outside AP-like resonance.

59. Use Jensen's arithmetic-progression safe-sum formula to average over
    progressions of length 2, 7, and 14.  This directly matches the `n=14`
    divisor structure.

60. Use mixed thresholds as a slack bookkeeping language.  Assign slightly
    stronger thresholds to dangerous coordinates and weaker thresholds to
    already-safe parity classes; prove the mixed vector lies in `MLPS_13`.

61. Reinterpret endpoint-debt product in zonotope covering radius terms.
    Malikiosis-Santos-Schymura and later zonotope papers suggest finite
    checking via LR zonotopes.  The repo's gap-debt invariant may be a lattice
    width/covariance invariant in that language.

62. Test coloopless/cosimple conditions on n=14 candidate zonotopes.  Recent
    zonotope work isolates classes where finite checking is transparent; see
    whether hard n=14 rows fail exactly one condition.

63. Use a dyadic fundamental-domain covering algorithm for n=14 LR zonotopes.
    The shifted-LRC computational paper used dyadic domains for rational
    polytopes; n=14 has natural dyadic plus 7-adic domains.

64. Translate "wall-only implies AP" into a covering-radius equality case.
    Equality cases in convex geometry may be easier to classify than interval
    covers directly.

## Endpoint and adic leads

65. Prove the n=14 gap-debt product cannot go to zero in both coordinates.
    THM-376 shows the hard ladders preserve `gap*debt = 5/11`.

66. Use the product building for primes 2 and 7.  The Bruhat-Tits frontier
    works cleanly for dyadic n=16; n=14 needs a two-prime product frontier.

67. Classify the six even bridge fibers from THM-375.  Each minimum local
    gate cover is forced odd fan plus one even bridge; global proof may be a
    statement that the six bridge choices cannot close cyclically.

68. Combine local bridge fibers with the 56-cell scalar moat.  Both have a
    six/seven/eight numerology; look for an exact incidence map between even
    bridges and the eight alpha stencils.

69. Use endpoint-pressure matroids.  A counterexample needs a leafless,
    private-pivot-free protection matroid, but sampled n=14 rows peel to empty.

70. Re-run pressure after deleting scalar/AP leaves first.  THM-377 says basic
    pressure lifts are acyclic on selected rows; maybe the cyclic core appears
    only after quotienting out AP wall debt.

71. Reverse-generate endpoint cycles, then solve for speeds.  Prior disproof
    sampling by dense gated speed sets leaks.  Build the desired nonpeeling
    endpoint core first and ask whether any primitive speeds realize it.

72. Use small-denominator sieve completion as a debt ledger.  Every time a
    missing modulus is repaired by an `lcm(n,d)` speed, endpoint debt exports
    to a deeper layer.  Prove export strictly increases a well-founded depth.

73. Define a potential `Phi = endpoint_depth - C*log(gap)` and test monotonicity
    under gate repairs.  Hard ladders keep product fixed; repairs may force
    `Phi` upward until a source opens.

74. Use no-14-multiple branch as a base case.  If no speed is divisible by 14,
    unit walls survive.  The whole proof can split: no-gate solved; gate branch
    must export debt and peel.

## Concurrent S531 rebase leads

75. Prove modular `H`-multiplicativity as a pruning theorem for n=14 survivors.
    If a bad class splits into disjoint apex modules, the Hamiltonian-path
    count factors and the quotient is transitive.  Then each module can be
    certified smaller than the full n=14 wall problem.

76. Treat nested apex chains as the true AP detector.  HYP-2010 says disjoint
    modules multiply but nested flips couple; the regular/AP obstruction is a
    nested spine.  Search finite survivors for "nested-only" certificates.

77. Convert recursive apex carving into a three-gap proof target.  Every runner
    carves the observer apex interval into sub-arcs; prove any non-AP recursive
    split leaves positive apex length before all thirteen cuts finish.

78. Add an apex-module coordinate to the permutohedral quotient.  Faces with
    disjoint modules should become product cells; the hard quotient should live
    on a lower-dimensional nested spine.

79. Use modular decomposition before scalar-gauge normalization.  Scalar
    quotients may hide independent apex modules; first factor disjoint modules,
    then run the scalar moat only on irreducible nested cores.

80. Generalize the n=4 parity law to a support-residue matroid for n=14.  For
    even n, Fourier coefficients vanish on multiples of n/2; at n=14 the
    active residues live modulo 7.  Enumerate which speed residue patterns
    forbid all active high-order resonances.

81. Build the n=14 analogue of the n=4 even-sum core.  HYP-2011 says n=4
    difficulty is exactly the residue class where inside debt can exist.  For
    n=14, find the CRT residue cores where all seven active channels can close.

82. Make the active mod-7 harmonic support into a tournament.  Orient residue
    classes by the sign and magnitude of the first nonzero Fourier coefficient;
    compare its cycles with the D_14 chord-channel debts.

83. Use "support width" as a monotone difficulty parameter.  n=4 has one active
    class, n=6 has two, n=14 has six.  A proof may proceed by showing width-six
    resonance still cannot match the AP's perfectly coherent debt unless all
    classes are occupied in AP order.

84. Use the n=4 even-core as a local model for every diameter gate.  Each n=14
    hidden diameter tie resembles a square inside-debt problem; attach an
    odd/even local law to each of the six diameter channels.

85. Turn the resonance-debt conjecture into an optimization objective for
    finite-check survivors.  Instead of only asking whether a survivor is lonely,
    compute its `abs(debt/credit)` and prove the maximum is reached only by the
    AP/scalar-ramp class.

86. Split n=14 debt into pairwise 83 percent plus higher-order residual.  If the
    pairwise Bernoulli part is already below AP by a margin larger than the
    high-order bound, the survivor is killed without full resonance expansion.

87. Unify the two HYP-2009 threads by dictionary: outside clasp debt equals
    Fourier resonance debt after summing hidden chord channels.  Prove the
    D_14 channel ledger is literally a grouping of the resonance-order terms.

## Independent-pair channel addendum

This addendum became HYP-2015 after rebasing over the concurrent HYP-2012
independent-pair metric, HYP-2013 coupling-gap boundary, and HYP-2014
almost-fixed-frame coupling result.

88. Replace support width by independent-pair rank.  The user pointed out that
    the K4 toy model is controlled by two independent pair-arcs under a fixed
    scaffold.  For n=14, record not just which channels are active, but how many
    vertex-disjoint pair states each channel can carry.

89. Track the six-bit diameter state.  The clasp-deleted n=14 diameter channel
    consists of exactly six mutually independent pairs plus one singleton.  This
    is the purest candidate state vector for the hidden endpoint/diameter debt.

90. Search for K4 windows in n=14 handoffs.  A local window should have fixed
    non-window scaffold plus two independent pair toggles whose four states
    determine a marked local tournament class.

91. Quotient by non-diameter scaffold first.  Channels `d=1..6` each have max
    matching 6 but seven possible maximum matchings; the diameter channel has a
    unique maximum matching.  Classify the scaffold, then read the diameter bits.

92. Add matching-shape entropy to the multi-channel parity law.  The obstruction
    is not only how many residue channels survive, but whether their independent
    pair matchings are unique, branched, or locked by endpoint labels.

93. Turn "wall-only => AP" into a finite independent-pair classifier.  Under a
    fixed non-diameter scaffold, prove that every non-AP diameter state opens
    the clasp or exports endpoint/resonance debt.

94. Use matching matroids for the channel state space.  Independent pair sets are
    matchings; cross-channel compatibility should be a matroid-intersection or
    Hall/Farkas problem over labelled endpoint debt.

## Near-term priority queue

1. Implement `lrc_k13_sieve_probe_s532.py`: compute survivor sizes for
   `I(13,p,1)`, then `c=2` and `c=7` lifts, with projection and scalar-gauge
   equivalence.
2. Implement exact CRT-class correlation tables using mixed-threshold safe
   integrals for the seven mod-7 classes.
3. Prove or refute the one-coordinate scalar-puncture moat symbolically from
   the eight alpha stencils.
4. Count n=14 open and boundary source menus in the fixed-safe-arc Ferrers
   model.
5. Build the dynamic apex interval-exchange trace for the known hard ladders,
   labelled by endpoint depth and CRT class.
6. Convert `wall-only => AP` into a finite equality-pattern classifier and run
   it on bounded non-AP families.
7. Add a product-bound scoreboard for the finite-checking route: `log product
   of certified primes / log B_13`.
8. Add `diameter_pair_state` and channel matching-shape features to the n=14
   hard-row audits.

## Most promising synthesis

The best new route is probably hybrid, not single-currency:

```text
finite-check sieve survivors
    -> collapse by scalar-gauge and CRT parity
    -> identify AP/scalar-puncture wall classes
    -> kill them by the 56-cell moat or by endpoint-debt export
    -> use product-bound scoreboard to finish k=13.
```

The geometric translation is:

```text
outside apex/clasp must open
unless the system is AP-like;
if it is AP-like, the scalar moat exposes explicit cells;
if a gate hides those cells, endpoint debt doubles and peels.
```
