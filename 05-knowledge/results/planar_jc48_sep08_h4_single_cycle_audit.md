# Independent audit of the H4 single-cycle obstruction

**Status: independent analytic and source audit PASS.**
The support theorem and the stronger equality corollary below are accepted.
The application to an actual mixed-cusp complement remains conditional on
its separate geometric supplier. No actual braid-path certification is
provided by this audit.

Audited artifacts:

- [Primary proof](planar_jc48_sep08_h4_single_cycle.md).
- [Standalone source](../../04-computation/planar_jc48_sep08_h4_single_cycle.py).
- [Frozen output](planar_jc48_sep08_h4_single_cycle.out).

## 1. Exact hypotheses and the general run lemma

The theorem uses four finite permutations, each one nontrivial cycle of a
common length `m>=2`, with fixed points elsewhere. The exact relations are

    aba=bab, bcb=cbc, cdcdc=dcdcd,
    ac=ca, ad=da, bd=db.

I independently reconstructed the general odd-braid lemma with rightmost
composition. If `sigma` fixes all orbit positions
`x,tau(x),...,tau^r(x)` and the odd relation has length `2r+1`, its two
sides at `x` reduce to `tau^r(x)` and `tau^(r+1)(x)`. They cannot be
equal for `x` in the moved support of `tau`. Repeated labels in a short
cycle are permitted; this is essential for the length-five argument at
cycle length two and is explicit in the final proof.

On each nontrivial `tau`-cycle of length `ell`, at least
`ceil(ell/(r+1))` positions belong to the support of `sigma`. Summing gives
the stated support bound and also the stronger sum of individual ceilings.
This argument does not require either permutation to be a single cycle.

For an ordinary braid and exactly half-overlap, every outside position
must be followed by an inside position. Equality of the two cardinalities
forces the converse transition as well. Thus the moved cycle alternates.
With several cycles, the nonnegative slack in each per-cycle inequality
must separately vanish when their total vanishes. The equality sidecar is
therefore justified at the exact scope used in the proof.

## 2. The support reduction is complete

For commuting single cycles of the same length, one moved support is
invariant under the other permutation. A nonempty intersection contains
an entire moved cycle, so the equal support cardinalities force equality.
Otherwise the supports are disjoint. This is the only dichotomy used;
it is not valid for general multiple-cycle permutations.

Write `A,B,C,D` for the four supports. If `A=C`, the terminal odd braid
makes `C` meet `D`; commuting `a,d` therefore forces `D=A`. The first
ordinary braid then makes `B` meet this same support, and commuting `b,d`
forces `B=D`. This already gives the common support.

If `A` and `C` are disjoint, the two ordinary half-support inequalities
fill the whole of `B`. Consequently `m` is even, `B` meets each of `A,C`
in exactly `m/2` points, and `B` is contained in their union. In particular
`c` alternates `C intersect B` and `C minus B`.

Commuting `a,d` forces `D` disjoint from `A`: the equality alternative
would contradict the positive overlap demanded by the terminal odd braid.
Commuting `b,d` similarly forces `D` disjoint from `B`, because `B`
already meets `A`. Every `C intersect B` position is thus outside `D`.
If one alternating `C minus B` position were outside `D` too, the three
consecutive orbit positions centred there would all be outside `D`,
contradicting the length-five run lemma. This also covers `m=2`, where
the predecessor and successor may be the same label. Hence

    I=C intersect D=C minus B,  |I|=m/2,

and `c` exchanges `I` with `C minus D`.

## 3. The final three-image contradiction

For `x in C minus D`, put `v=c(x) in I`. Since `d(x)=x`, the two
length-five words at `x` are `y=d c d(v)` and `c(y)`. Thus `c(y)=y`.
If `d(v)` were in `I`, alternation would put `c d(v)` in `C minus D`,
where `d` fixes it. Then `y=c d(v)` would lie in `C`, contradicting
that it is fixed by the single moved cycle `c`.

Therefore `d(v)` lies in `D minus C`. Now `c` fixes `d(v)`, so the
same identity gives `y=d^2(v)` in `D minus C`. Since `c` maps all of
`C minus D` onto `I`, this holds for every `v in I`. In particular
`I,d(I),d^2(I)` are disjoint: the last two intersections reduce under
the bijection `d` to `I intersect d(I)`, which is already empty.
All three have cardinality `m/2` and lie in `D`, of cardinality `m`.
This is impossible. The support theorem follows for all ambient degrees
and every `m>=2`, without extrapolating a finite census.

## 4. A stronger equality corollary

The same hypotheses force `a=b=c=d`, not only equal supports. This
corollary was derived independently during the audit and is accepted
analytically; it requires no source change.

After support coincidence, restrict to the common moved set. A permutation
commuting with the single cycle `d` is a power of `d`: its value at one
point determines its value at every point by equivariance along that
cycle. Thus `a` and `b`, which both commute with `d`, commute with each
other. The relation `aba=bab` then cancels to `a=b`. Next `ac=ca`
and `bcb=cbc` give `c=a`. Finally `ad=da` and the terminal odd relation
with `c=a` give `d=a`. Every generator is the identity on the complement,
so this is equality as permutations of the original finite set.

Consequently their image is cyclic. If the action is transitive, the
common nonempty moved support is the whole set, and the degree is `m`.
There is no retained fixed label. This conclusion still uses the single
cycle hypothesis; the centralizer of a multiple-cycle permutation need
not be cyclic.

## 5. Source audit and exact universes

I read the entire source. It imports no inherited implementation, and its
`need` checks raise explicit errors independently of Python optimization.
The composition convention is correct: `mul(a,b)` is `a after b`, and
the alternating odd words are the literal displayed relations.

The finite universes and filters are paid as follows.

1. The ordered pair bank contains every pair in `S_d`, `2<=d<=5`, for
   odd labels three, five, and seven. It imposes no support-type filter;
   identities and multiple cycles remain. Only literal braid solutions
   are checked against the run lemma and support conclusions.
2. The relative single-cycle bank fixes the first cycle and ranges over
   every positive overlap subset and every ordering of the second cycle,
   for `2<=m<=7`. Its union degree is `2m-j`. Disjoint support is omitted
   only because it is already ruled out by the proved run lemma. The
   final primary wording now states this inherited filter explicitly.
   Fixing the first cycle by relabelling and using every ordering on the
   second support is complete; duplicate representatives would be harmless.
3. The H4 tuple bank fixes the first cycle and enumerates all remaining
   cycles on each declared ambient set. Its filters are exactly the six
   displayed relations. It does not impose support coincidence or
   transitivity before checking them. The postcomputed orbit uses forward
   generators, which give the full generated orbit for finite permutations.
4. The alternating equality bank fixes `c` standard and the intersection
   at every second position, then takes every ordering of `d` on that
   intersection and the required fresh labels. Its hypotheses are exactly
   the impossible final support configuration; no case is inferred from
   a mere union cardinality.

The frozen full-tuple census contains exactly one accepted tuple for each
of its 21 declared `(m,degree)` rows, consistent with the new equality
corollary. That observation is not its proof: Section 4 proves equality
without a finite cutoff. The source's explicit target gate remains the
support theorem, and its source/output have not been altered for the
analytic strengthening.

I also independently checked the named five-letter permutations with
SymPy's separate permutation implementation: `(123),(345)` satisfy the
five-term relation and fail the ordinary braid while their supports meet
in one point. Thus the old half-support claim really fails in that scope.
The adjacent-transposition chain `(12),(23),(34),(45)` satisfies all
other required relations and is transitive with fixed letters, but fails
the terminal five-term relation. Both are load-bearing hostile controls,
not examples of the theorem's full hypotheses. Equal generators on one
`m`-cycle are valid positive controls in every length.

## 6. Independent replay and pins

I ran fresh ordinary and optimized source replays and compared both with
the frozen output. All three are byte-identical, with **15,111 always-active
exact gates** and **3,650 output bytes**.

Reproduction from the worktree root:

```bash
python3 04-computation/planar_jc48_sep08_h4_single_cycle.py
python3 -O 04-computation/planar_jc48_sep08_h4_single_cycle.py
```

Frozen SHA256 pins independently checked:

- Source, 7,003 bytes:
  `c42119cd358e291f3adf1d1d6b3c049fb44f65470ad30ae60ff54616b92d54b2`.
- Output, 3,650 bytes:
  `968390414754d58917891eb2f864aa8d62ebbf9e7ff2d6d8182c6f1dfaaa34ae`.
- Semantic digest:
  `633a979b411bed9b53142c98fa5a3c4e50d4a6cde4e07bbdafbad93806d92cdd`.

The independent replay files are temporary
`/tmp/h4_single_cycle_orthogonal_normal.out` and
`/tmp/h4_single_cycle_orthogonal_optimized.out`. No producer source,
output, shared navigation, or Git state was modified by this audit.

## 7. Actual geometry remains a separate gate

The abstract argument is complete. To consume it for a Keller map requires
an actual common positive-meridian generating tuple with the prescribed
local pairs, a single nontrivial cycle in the meridian type, and at least
one genuinely retained fixed sheet. Merely finding algebraically
conjugate pairs in a heuristic braid word does not prove that their
retained subsets have the asserted actual local access.

The exact literal six-word H4 presentation was separately checked in a
temporary free-word verifier, including the pre-relation inverse
`old a=x*y*z*y^-1*x^-1`. That verification is only a presentation
statement for the literal word bank. It does not turn numerical paths
into certified paths or supply local-cluster markings. No general
multi-cycle exclusion, mixed-curve closure, or JC(2) claim is made here.

The final primary-text reread accepts the promoted mathematical statement
that all four permutations are equal, the complete centralizer paragraph,
and the explicit positive-overlap filter in its finite-universe description.
No mathematical correction remains. Source/output pins are unchanged;
root owns status promotion and integration.
