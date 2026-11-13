# Independent audit of short tournament-mutation returns

**Status: VERIFIED independent proof review and FINITE-EXACT case-table
reproduction.** This audit concerns sections1–2 of
[the sibling proof](../ports/PROOF.md), not its modular observation code.
It makes no external priority claim. The incoming Session291 inheritance
must be retained separately from validity of the arguments.

The following all-size arguments were checked without using the sibling
enumerator as a premise:

- A mutation is native to the tournament locus exactly at a global source
  or sink: any nonempty in/out-neighbor pair creates an entry0 or2.
- The interval from the first to the last nontrivial SCC is immobile under
  native mutation. The exterior singleton vertices move in their cyclic
  order across this hull, giving a path with `t+1` vertices. Its ordinary
  isomorphism classes are distinct because the count of singleton SCCs
  before the first nontrivial SCC is intrinsic. The transitive cases give
  cyclic rotations labelled and a single ordinary class, including the
  separate orders1 and2 boundary cases.
- Every ordinary tournament has one unique hull and exterior split.
  Consequently `a_n=1+sum_(h=3)^n(n-h+1)k_h`,
  `k_n=a_n-2a_(n-1)+a_(n-2)`, and
  `c_n=a_n-a_(n-1)+1` have the stated directions and scopes. They classify
  a quotient and do not compute unknown A000568 terms without input.
- After a nonnative first mutation at `k`, every other vertex is incident
  with a bad cross-part entry. A second mutation there only negates that
  incident entry; hence a two-step non-backtracking return is impossible.
- Restriction to a vertex set containing every mutation index commutes
  with mutation. The complete three-vertex `k,l,k` calculation forces
  precisely the two extreme two-singleton forms in the proof. In those
  forms the remainder is unchanged and its Hamiltonian count is retained.
- The complete three-vertex `k,l,m` calculation with distinct names forces
  a transitive orientation on those names. If the first mutation was
  nonnative, an incoming fourth vertex exists; each of its four possible
  relations to the other two names leaves a bad matrix entry. Thus the
  local obstruction proves exclusion at arbitrary order.
- Any path of at most three mutations between tournaments splits into
  native steps and the classified excursions, or cancels an adjacent
  repeated mutation. Each piece preserves H. A count-changing return
  therefore requires at least four mutations. Existence at length4 is
  not asserted.

No mathematical repair was needed. The source/sink and native-component
claims are inherited where the current sibling proof identifies them.

## Four-step upgrade

The subsequent [four-step proof](../ports/FOUR_RETURNS.md) also passes
independent review. Its all-size claim is stronger than a bounded census
for a valid reason: every matrix entry after a word supported on `S`
depends only on the induced input on `S` together with the entry's two
endpoints. Thus a bad intermediate entry would already give a counterexample
on at most `|S|+2` vertices. The finite entry implications used for0102 and
0121 therefore do extend to every order.

All five reduced four-letter word shapes are present. For0101 the inherited
full-matrix pentagon and the source/sink criterion give the exact answer.
For0102 and0121 a forced intermediate tournament splits the word into
already proved H-preserving pieces. For0123 the local return condition
forces an extreme TT4 block, so every step is native. For0120 the three
named vertices form a directed triangle whenever an exterior vertex is
present; each exterior vertex lies uniformly on one side, and the
five-vertex test forces all of them onto the same side. The output
reverses and transfers this extreme C3 while preserving the exterior
tournament, so both Hamiltonian counts equal `3H(R)`. At order3 the
explicit six-row table separately preserves H.

The independent checker now exhausts all1,096 labelled tournaments of
orders3,4,5 for the applicable word shapes. It verifies the exact0120 and
0123 return criteria in both directions, the0102/0121 intermediate
implications, the0101 source/sink implication, the full-matrix pentagon,
and the stated0120 output entries. It imports no sibling implementation.

Therefore **every return of length at most4 preserves H**. Any count-changing
return must have length at least5. Existence at5 or any length remains
open; the bounded longer scout is not used in this proof audit.

## Independent finite check

[audit_mutation.py](audit_mutation.py) uses its own mutation implementation
and regenerates all eight sign assignments on three vertices for each of
the words `k,l,k` and `k,l,m`, then all four incoming-fourth-vertex cases.
It does not import the sibling code or its receipt. Its output is
[mutation_audit.json](mutation_audit.json).

The exact return rows are:

```text
k,l,k: (-1,1,1), (1,-1,-1)
k,l,m: (-1,-1,-1), (1,1,1)
```

For the normalized positive `k,l,m` orientation, the final remaining
column is respectively

```text
(-5,-2,4), (-3,0,2), (-1,-2,2), (1,0,0)
```

for `(e,f)=(-1,-1),(-1,1),(1,-1),(1,1)`. These exactly match the proof.

Reproduction:

```text
python -B audit_mutation.py
python -O -B audit_mutation.py
```

The finite rows and four-step implications corroborate the stated local
calculations; the all-size consequence comes from the restriction, SCC
and explicit block-operation arguments above.
