---
source: codex-2026-06-03-S580
status: comparative synthesis and proof-program
tags:
  - lonely-runner
  - finite-sieve
  - certificate-compression
  - antipodal-witnesses
  - summand-graph
  - tournament-analysis
---

# Paper Sieve vs Certificate Compression for total n=11,12,13

The question was how much more efficient the repo methodologies could be than
the exact methods from the `n=11,12,13` paper, and how to push the improvement
mathematically rather than only by code.

The answer is: if the certificate quotient is complete, the improvement is not
a constant-factor implementation speedup.  It changes the object of search.
The paper proves exact emptiness of `J(k,p)` across many primes.  The repo
methods try to prove that the only exact fibers worth generating are those that
survive four much cheaper clocks:

```text
n-clock
2n-1 antipodal shell clock
pair-sum/pinch clock
D/U/N endpoint-owner clock
```

## Paper Scale

The 2026 paper uses `k` moving speeds and total runner denominator `n=k+1`.
Its proof table has:

```text
total n=11, k=10: 61 primes, log product 342.716, about 45 minutes
total n=12, k=11: 76 primes, log product 435.047, about 40 hours estimated
total n=13, k=12: 91 primes, log product 547.381, about 40 days estimated
```

The proof uses two kinds of savings before the repo enters:

- symmetry and DFS to compute `I(k,p,1)`;
- lift/project diagrams, with the odd-prime polynomial method avoiding the
  direct `c=k+1` lift for `k=10,12`.

That polynomial shortcut is enormous.  For fixed input set, the direct
`c=11` lift costs `(11/2)^10 = 2.53e7` times a `c=2` lift, and the direct
`c=13` lift costs `(13/2)^12 = 5.69e9` times a `c=2` lift.

## Certificate Scale

The S580 script gives this model:

```text
n  raw sum p^k  full exact ansatz  cert scan ops  raw->cert  final->cert
11 10^27.63     10^38.05           10^4.86        22.77      33.19 orders
12 10^31.27     10^55.01           10^5.26        26.01      49.75 orders
13 10^35.31     10^48.68           10^5.61        29.70      43.06 orders
```

This does not mean we have a proof in `10^5` operations.  It means the current
proof architecture has a plausible target table that small: all antipodal shell
transversals times pair clocks and D/U/N obligations.  The missing theorem is
that this quotient is complete.

## Which Clocks Matter

The `n`-clock matters first.  HYP-2102 reduces LRC(n) to C': if no speed is a
multiple of `n`, `t=1/n` is already a witness; if a speed is a multiple of `n`,
we need to prove that the set is loose.  THM-398/HYP-2103 sharpen the large
branch to a dominance-dodge theorem: any speed bigger than `(n-1)` times the
max of the others is loose.  The real residual is the all-short small-multiple
case.

More importantly, THM-398 gives an interval criterion, not just a dominance
criterion.  Choose any runner `v`.  If the safe set of the other runners,
`G(S\{v})`, has a connected component longer than `2/(n v)`, that interval
cannot fit inside a single danger arc of `v`, so it contains a full witness for
`S`.  For the multiple branch `v=nw`, the only surviving exact work is therefore
the all-short AP-cover case:

```text
n=11, v=11: all-short cap 0.016529; n-clock safe gap 0.074380
n=12, v=12: all-short cap 0.013889; n-clock safe gap 0.069444
n=13, v=13: all-short cap 0.011834; n-clock safe gap 0.065089
```

For larger multiples the cap divides by `w`.  This is a stronger algorithmic
message than "large multiples are easy": most candidates should die by a
one-runner interval sweep before the shell/pair/D/U/N quotient is even built.

S573/HYP-2104 adds the right measure-theory name: `n | v` is the Vitali handoff.
No multiple means the constructive `t=1/n` side, including the measure-zero
worry core; a multiple `v=nw` creates the periodic Vitali-cover family with
centres `k/(nw)` and radius `1/(n^2 w)`.  Its sample computation shows B' proving
`91.5%` of the mult-of-`n` branch at `n=12` and `96.8%` at `n=14`, with zero
all-short tight cases across the sampled `n=6,8,10,12,14`.  So the residual is
not generic finite checking.  It is a Diophantine arc-alignment problem.

The `C=2n-1` clock matters second.  It is the summand graph clock: additive
shells `{a,C-a}` with multiplicative unit action.  Here the paper's arithmetic
and the repo's arithmetic diverge:

```text
total n=11: C=21, units/nonunits = 6/4
total n=12: C=23, units/nonunits = 11/0
total n=13: C=25, units/nonunits = 10/2
```

So `n=12` is the clearest rotated opportunity.  It is composite for the paper's
`k+1` lift, but prime for the repo's `2n-1` antipodal shell clock.

Pair-sum/pinch clocks matter because HYP-2095 says measure-zero rows should
have an unblocked small reduced-sum pair.  This is where exact time is cheap:
one rational pinch can certify what a full ansatz grid would otherwise search.

D/U/N endpoint-owner clocks matter because they preserve what the quotient is
trying to prove.  They remember which endpoint, divisor, unit shell, or
`n`-clock tick is being protected.  Raw runner tuples remember too much and too
little at the same time.

## Which Cases To Ignore Early

Do not start by scanning generic denominators, arbitrary `lp` grids, raw
runner-state tournaments, or complete prime fibers.  Defer them until one of the
primary clocks returns a labelled residual:

```text
no n-multiple? solved by t=1/n
n-multiple but dominant? solved by THM-398/HYP-2103 dodge
chosen v has long G(S\{v}) component? solved by THM-398 interval criterion
C-shell unit complete? normalize through unit action
small pair unblocked? pinch witness
D/U/N core peels? no counterexample core
only then: exact prime-fiber fallback
```

## Creative Algorithmic Upgrades

1. Reverse the paper generator.  Generate endpoint cores and pair-blocker
   certificates first; solve for speed tuples only after a core survives.
2. Add a dodge-first interval oracle.  For each promising `v`, sweep the
   endpoints of `G(S\{v})`; any component longer than `2/(n v)` certifies
   looseness.  The fallback object is no longer a tuple, but an all-short
   cover assignment of tiny safe components to the AP arcs `{m/v}`.
3. Use clause learning across primes.  A paper seed family that fails for
   several primes should become a prime-independent certificate germ.
4. Store D/U/N ledgers in a zero-suppressed decision diagram.  Most obligations
   are absent; ZDDs match the proof object better than vectors in `Z_p^k`.
5. Prove a small-multiple equidistribution lemma.  If C' is proved, the entire
   product-sieve layer becomes unnecessary for fixed `n`.
6. Use the S573 Vitali form to quantify alignment: the centres of all-short
   safe intervals are rational functions of the other runners, and a cover
   forces them to be `1/(n^2 w)`-close to the `1/(nw)` lattice.  That is an
   Erdos-Turan/three-distance target, not a grid-enumeration target.
7. Translate all-short cover assignments into endpoint-owner congruence
   constraints.  A component of `G(S\{v})` has two boundary owners; forcing it
   inside one `v`-danger arc creates pair-sum/divisor/unit relations.  This is
   the bridge from interval search to D/U/N ledgers.
8. Use HYP-2111's certificate calculus as the fallback order.  The exact paper
   fiber should only appear after `G3_all_short_Cprime_residual` fails to route
   into `G8_endpoint_cover_circuit_positivity`.
9. For total `n=12`, write the prime-`C=23` unit-shell proof as a standalone
   lab.  There are no nonunit shells, so every obstruction must be endpoint or
   pair-blocker in nature.
10. For total `n=13`, isolate the two nonunit `C=25` shell holes and test whether
   all paper `k=12` raw seed families project to AP plus those two holes.
11. Treat addition and multiplication separately.  Addition creates antipodal
   shells and pair-sum pinches; multiplication moves unit representatives.
   Exact lifts mix them, which is why they are expensive.

## Assumption Challenge

I explicitly reject "vertices are runners" for this stage.  Candidate vertex
sets considered:

```text
runners, pair-sum cells, antipodal shells, fixed circle sections,
endpoint rows, D/U/N obligations, paper seed families, certificate germs,
prime-lift clauses, Fourier modes
```

The quotient must preserve the predicate "a lonely witness or proof certificate
exists."  It can destroy exact speed order, actual prime residue labels, and
generic time-grid coordinates.  That destruction is a feature if it keeps the
endpoint owner, shell, and pair-blocker data.

Tournament Analysis vertices for the next script should be certificate germs,
not runners.  Edge `A -> B` means `A` has larger residual log-work after applying
the same clocks.  Ties use the Hamiltonian path

```text
n-clock -> C-shell -> pair-pinch -> D/U/N core -> exact fallback
```

The current fingerprint is transitive; cycles should be sought only inside
endpoint-owner fibres or cross-prime clause conflicts.

## Handoff

The most leveraged next computation is a paper-seed translator:

```text
paper seed family -> C-shell multiset -> pair-blocker profile -> D/U/N owners
```

Run it first for public `k=12` raw logs.  If the recurring `12` seed families
all collapse to AP plus the two `C=25` nonunit holes, the `n=13` proof route
becomes a two-hole endpoint theorem rather than a forty-day prime verifier.
