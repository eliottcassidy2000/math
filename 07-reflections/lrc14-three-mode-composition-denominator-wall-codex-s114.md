# LRC14 three-mode composition and the committed-denominator wall

codex-2026-06-22-S114.

The owner's correction changes the right mental model.  The three tournament
recursion modes are not interchangeable sign strings.

```text
Mobius:     A+B+C-D-E-F+G       all sizes
Eisenstein: A+B-C               even sizes
Legendre:   A+B+D-C-E-F+G       odd sizes
```

The Legendre string is a real Venn diagram.  Its corners are `A,D,B` with
sizes `N-1,N-2,N-1`; its pairwise overlaps are `C,E,F` with sizes
`N-2,N-3,N-3`; its triple overlap is `G` with size `N-4`.  The edge unions are
`A+B-C`, `A+D-E`, `B+D-F`; the center/full union is
`A+B+D-C-E-F+G`.  The earlier `(2,0,-2,1)` coefficient vector is only the
cardinality shadow after the size-`N-2` corner and size-`N-2` overlap cancel.

This is exactly the distinction LRC14 needs.  The finite cap node should live
in the corrected Venn geometry: AP/three-gap rigidity is a statement about
which Venn regions are forced to move together under phase rotation.  The
analytic floor node lives somewhere else: exact-period denominator packets and
equidistribution.

The lcm family makes that split unavoidable.  For

```text
S_X = {1,...,11,13,lcm(2,...,X)}
```

every denominator `D<=X` is dead, because the committed speed is zero mod `D`.
So no reduced fraction `a/D` can be a `1/14` witness.  This is a short proof,
not a heuristic, and it refutes every fixed finite denominator-basis closure.

The first denominator above the wall is not simply `nextprime(X)`.  The S114
script finds many mismatches; for instance `X=60` first opens at `67`, while
`nextprime(X)=61`, and `X=110,120` first open at the composite denominator
`121=11^2`.  That points to the right analytic object: a prime-power/residue
opening after the divisor wall, not a prime-only rule.

Incoming mac-mini S45 adds the companion radical filter.  If some prime
`p<=13` divides no speed, then `t=1/p` is already a witness: every phase is a
nonzero residue mod `p`, so its distance from the observer is at least
`1/p >= 1/13 > 1/14`.  Thus a genuinely hard row must kill the primes
`2,3,5,7,11,13` and also kill `14`.  That is the first branch of the Node-3
proof.  The S114 composite-denominator data adds the second warning: once the
small-prime branch is saturated, exact first witnesses live in prime-power
packets, not in a prime-only list.

The practical proof split I would keep:

1. Node 2 finite cap: prove AP/consecutive extremality using the Legendre Venn
   labels.  Do not collapse `D` and `C` just because they have the same size.
   One is a corner, one is an overlap, and that is the finite realizability
   information.
2. Node 3 / Part A: first apply the radical filter (`p<=13` surviving gives
   `t=1/p`), then prove an effective equidistribution lemma for denominators
   just beyond the committed wall.  This should be stated over exact-period
   unit packets with prime-power labels, then intersected with `GOOD cap G_P`
   and charged the `#arcs/Vmax` loss.
3. Mobius remains the common skeleton: divisor inclusion-exclusion
   (`phi=mu*id`) on the analytic side and Venn inclusion-exclusion on the
   finite side.  The scalar cap/floor estimates should see only the final
   labelled packets, not the raw runners.

Assumption challenged: tournament vertices do not need to be runners, arcs, or
even graphs.  Here the productive vertices are recursion modes, Venn regions,
denominator packets, prime-power openings, and proof obligations.  This
preserves the dependency structure of the LRC proof but deliberately destroys
the circle-time geometry, which must be restored by the analytic Part-A node.

Files:

- `04-computation/lrc_lcm_committed_denominator_wall_codex_s114.py`
- `05-knowledge/results/lrc_lcm_committed_denominator_wall_codex_s114.out`
- `05-knowledge/hypotheses/HYP-2901-lrc14-three-mode-composition-and-lcm-denominator-wall.md`
