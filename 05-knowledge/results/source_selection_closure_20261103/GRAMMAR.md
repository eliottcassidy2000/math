# A physical cone, a canonical inverse suffix, and the remaining singular core

**PROVED for the inherited 21-entry/28-map grammar:** the two full-rank
commands are recognizable from the output profile, against every other
last command; their maximal suffix has a unique exact inverse. A finite
inverse procedure decides membership of a profile in this actual-source
grammar and returns a literal program. **OPEN:** a selector that produces
every required prime/7p target. No scalar-coverage interval is enlarged.

## G0. Inheritance and the live board

The closest mechanism is [ROOTED_RESETS R1–R3](../two_port_return_induction_20261103/ROOTED_RESETS.md):
the root-child interface, the positive four-profile grammar, and its exact
continuation divisor. Its actual compiler retains physical roots, controller
slot changes, and the final provisional birth. The canonical hostile is the
[marked H5 weak-predecessor example](../two_port_return_induction_20261103/OBSERVER_INDUCTION.md),
where a realizable numerical inverse is outside the strong seed language.
The corrected near miss here is using only nonnegative coordinates as the
physical cone. The least-used sidecar is the insertion of a root **away
from the first position** of a Hamilton path.

The incoming [eight-macro projective decoder](../spectrum31_return_20261103/PROJECTIVE.md)
provides a useful comparison, rather than a dependency: it decodes a
different family by disjoint ratio intervals. Here exact four-dimensional
cone separation recognizes two letters even in the presence of the other
26 grammar commands. Its integer-source guards are retained.

| Live concept | Cheap structural test | Outcome |
|---|---|---|
| Physical positive cone | Insert the root noninitially | The extra inequality h−s≥q holds whenever the root has positive indegree. |
| Full-rank inverse | Compare the 041/403 images | Their images are disjoint; both preserve the odd profile gcd. |
| Singular interference | Compare each of 26 other last letters with both full-rank maps | 54 exact dual certificates recognize the full-rank last letter, including comparison of the two full-rank maps. |
| Actual source entry | Decode until an entry or singular boundary | A unique full-rank suffix, then finite bounded inverse fibres; accepted outputs have literal source programs. |
| Target selection | Use only a bounded number of singular commands | Too few states to contain all primes or all7p values; the singular core cannot have a fixed budget. |

The anchor is constructive source selection. The niche is physical-cone
recognition. The wildcard is the information loss of a singular command:
it becomes the precise location where an all-target inverse must choose.

## G1. A stronger physical cone and a growth bound

For a rooted tournament use the inherited profile

    h=H(T), s=#HP starting at v, q=H(T−v),
    F=#ordered two-slot covers (P,Q), P ending at v, Q possibly empty.

If v has positive indegree, then

    h−s ≥ q.                                           (G1)

Indeed, take any Hamilton path of T−v and select its last vertex that
points to v. Insert v immediately after it. If it is not the final vertex,
its successor is beaten by v by the choice of the last in-neighbor.
The resulting path does not start at v. Deleting v recovers the old path,
so this is an injection into the h−s paths not starting at v. The positive
indegree guard is necessary: a universal source can have h=s=q>0.

Write

    x=(b,t,s,δ)=(h−s−q, q−s, s, F−h−q+s).

For every strong source, x is in

    C={x≥0 : s>0}.                                    (G2)

Here q≥s follows by deleting the first vertex of a path starting at v;
δ≥0 is the inherited endpoint-cover identity; and a strong tournament
has a Hamilton path starting at each vertex. The older coordinates are

    u=(h−s,q−s,s,δ)=(b+t+s,t,s,δ),
    h=b+t+2s, q=t+s, F=b+2t+2s+δ.

This is an integral unimodular change of coordinates. Conjugating all28
inherited maps gives nonnegative integer matrices M_w on x; the receipt
checks every entry. The cone is a necessary physical condition, not a
general realization theorem.

The first provisional birth in any macro has count

    h_first=2b+3t+3s+δ ≥ (3/2)h.                      (G3)

The output parent contains that birth; subsequent births preserve
strongness and increase H. Thus every macro has H_out≥(3/2)H_in. Equality
in (G3) would require b=t=δ=0 and give even h, so actual odd-H sources
have strict inequality at that boundary. Only the stated weak bound is
needed below.

## G2. Full-rank images recognize their own last letter

Set A=041 and B=403. In the x coordinates their matrices are

    A=[[4,1,2,1],       B=[[2,3,1,4],
       [6,3,2,2],          [4,3,2,6],
       [5,3,3,3],          [5,3,2,6],
       [3,1,1,1]],         [2,2,1,2]].

Their determinants are 4 and −4. For an actual source, gcd(x) is odd
because it divides h. Integral forward maps and dyadic inverses imply

    gcd(xA)=gcd(xB)=gcd(x).                            (G4)

Their counts obey

    H(xA)=9b+13t+14s+6δ ≥7h,
    H(xB)=7b+11t+12s+6δ ≥6h.                          (G5)

The following is stronger than distinguishing A from B.

**Last-letter theorem.** For N in {A,B} and any other one of the28 maps M,

    {xM : x∈C} ∩ {yN : y≥0} = empty.                 (G6)

The proof is an exact finite coefficient certificate. For each of the54
ordered pairs (M,N), the receipt gives an integer column λ≥0 such that

    c=M N⁻¹ λ ≤0 coordinatewise,   c_s<0.             (G7)

If xM=yN with x∈C and y≥0, then yλ=x c<0, contradicting yλ≥0.
All coefficients are verified over rational arithmetic. This is a proof
for all points of C, not a test on sampled source states. Examples are

    M=1, N=A:   λ=(0,4,3,1), c=(0,−1,−1,0);
    M=1, N=B:   λ=(3,0,0,1), c=(0,−3,−3,0).

The full list contains only small nonnegative integer λ, with largest
entry11. Scientific/numerical linear programming was used to discover
these columns; reproduction uses only the supplied integers and exact
rational multiplication, with no optimizer dependency.

**Canonical suffix.** Given a profile known to belong to the whole grammar,
test xA⁻¹ and xB⁻¹ for membership in C. If one succeeds, it is the unique
possible last macro, its actual predecessor is forced, and the test may be
repeated. If neither succeeds, the state is an initial entry or a singular
last-macro output. This recovers a unique maximal full-rank suffix and its
core profile. The core may have several singular histories; the theorem
does not make the whole28-letter grammar free.

The source guard is sharp at the entry boundary. The entries with words041
and403 invert under the corresponding full-rank map to

    x=(0,1,0,0), z=(h,s,q,F)=(1,0,1,2).

This is the actual weak TT2 rooted at its sink. It is nonnegative and
integral but has s=0, so is outside C and outside the strong C3 grammar.
These two entries must not be silently decoded into a strong predecessor.

## G3. An exactly decodable primitive family

Start with any of the21 actual entries and use only A and B. This is an
infinite family with a deterministic inverse: test the two inverse cones,
require integral predecessors in C, and stop only at a listed entry.
Each step decreases H by a factor at least6, so the procedure terminates.
An accepted profile comes with its entry and word; replay through the
inherited controller compiler produces the actual strong rooted graph.

All21 entry counts are at most13, whereas the first A/B output has count
at least18. Together with (G6), this proves uniqueness of entry and word.
Distinct accepted profile vectors distinguish rooted isomorphism classes
within this constructed family. The result is not recognition of an
arbitrary supplied graph from its profile; other graphs can share a profile.

The21 entries have orders3,4,5 with multiplicities1,4,16. Both letters add
three vertices. Therefore the generating function for the distinct rooted
classes in this specific family is

    (z³+4z⁴+16z⁵)/(1−2z³).                            (G8)

Exactly20 entries are primitive. The only exception is entry021 with
x=(3,0,3,0), gcd3. By (G4), the primitive subfamily has generating function

    (z³+4z⁴+15z⁵)/(1−2z³).                            (G9)

This is an exact actual-source family and gcd-preserving control system.
It is not itself an all-target construction.

## G4. A finite inverse procedure for the full grammar

When full-rank peeling stops, all genuine choices are among26 singular
maps. They have rank2 or3. For a candidate output y with H(y)=N, a genuine
parent x necessarily satisfies

    xM=y, x∈C∩Z⁴,
    2b+3t+3s+δ ≤ N.                                  (G10)

The bound follows from the first provisional birth, not from an assumed
realization of a formal inverse. It bounds each coordinate by N. Row
reduction of Mᵀ leaves at most two free integer coordinates, so the fibre
can be enumerated finitely. Each surviving candidate has H(x)≤2N/3.

The following recursive algorithm is sound and complete for **profile
membership in this grammar**:

1. Accept a listed entry, retaining its literal source prefix.
2. If either full-rank inverse lies in C, it is the only possible branch
   by (G6). Reject a nonintegral predecessor; otherwise recurse there.
3. If neither does, enumerate (G10) for each singular map and recurse.
   Accept if a predecessor accepts, retaining that macro.

The decreasing H makes the recursion terminate. Completeness follows
because every genuine grammar parent satisfies (G10). Soundness follows
by reversing an accepted chain of actual entry/macro programs. No numerical
profile is accepted merely for being integral, positive, primitive, or
inside a last-command image.

For a scalar target N the candidate profile search is finite as well. Every
entry and every macro output satisfies δ≤3h; for outputs this is the exact
coefficient inequality M[:,δ]≤3M[:,H] checked for all28 maps. Thus enumerate

    b+t+2s=N, b,t,δ≥0, s≥1, δ≤3N,

and apply the membership algorithm. There are O(N³) such outer candidates.
This is a terminating target-directed decision procedure, not an efficient
binary-input complexity claim or a theorem that it accepts every target.
The experiment tests only the declared small structural controls below.

## G5. Singular choices cannot have a fixed budget

The two full-rank maps alone preserve primitivity but cannot contain all
primes or all7p targets. More generally, fix any integer K≥0 and allow at
most K singular macros in a program. If its final H is at most X, (G5)
bounds the number of full-rank letters by L=floor(log₆(X/3)). The number
of program presentations is at most

    21 Σ_{j=0}^K 26^j Σ_{ell=0}^L binom(ell+j,j) 2^ell
       = O_K(X^(log₆2) (log X)^K).                   (G11)

Singular macros only increase H, so interspersing them cannot invalidate
the full-rank length bound. Counting presentations is an upper bound even
when singular words collide.

For completeness, the required prime-count comparison needs no unproved
distribution hypothesis. The central binomial coefficient is at least
4^n/(2n+1). For each prime p, its exponent in binom(2n,n) is a sum of at
most floor(log_p(2n)) terms, each0 or1. Its prime-power contribution is
therefore at most2n, so

    (2n)^(pi(2n)) ≥ binom(2n,n) ≥4^n/(2n+1).

Hence pi(X)=Omega(X/log X), which exceeds (G11) for each fixed K. The same
comparison at X/7 applies to7p. Consequently every fixed singular budget
misses infinitely many primes and infinitely many7p targets. A successful
all-target construction within this grammar must use an unbounded number
of singular choices. This does not rule out a finite-state controller with
unbounded execution or an invariant family that handles those choices.

The constructive problem is now localized: full-rank suffixes have unique
source-preserving inverses; one needs a productive rule on the singular
core, with actual entry and continuation invariants. A finite inverse
search supplies a decision procedure, while an all-target selector remains
the missing theorem.

## G6. Exact receipts and independent graph controls

`grammar_inverse.py` uses only the Python standard library and records the
hash of its inherited matrix receipt. Its finite universe is:

* all28 exact coordinate changes, growth/defect inequalities, and54 dual
  last-letter separation certificates;
* all21 entries and all A/B words of length at most2:147 literal rooted
  graphs, maximum order11, with independent all-subset endpoint DP profiles;
* all26 singular inverse fibres over the C3 entry after one macro, and28
  one-macro membership decisions, each accepted program replayed literally;
* the two actual weak inverse entry boundaries, explicitly rejected as
  strong grammar sources.

The full binary decoder recovers every supplied entry and word exactly.
These are structural and compiler controls, not a support census. Normal
and optimized runs match the frozen JSON byte for byte:

```sh
python3 -I -B 05-knowledge/results/source_selection_closure_20261103/grammar_inverse.py
python3 -I -B -O 05-knowledge/results/source_selection_closure_20261103/grammar_inverse.py
```

Independent proof review passed G1, G2, G4 and G5. A separate Fraction
Gauss–Jordan implementation inverted A/B and recomputed all54 dual
contractions, including their strict signs. The reviewer checked that the
entry test covers all42 entry/inverse pairs and isolates exactly the two
weak TT2 boundaries. That review did not independently replay the147
literal graph controls. The coordinator separately reviewed the generating
functions, gcd/growth arguments and finite inverse bound.

The connection map is rooted source → noninitial-insertion cone → separated
last-command images → a smaller actual predecessor. It preserves rooted
continuation profiles and the literal compiler, loses graph information
outside that continuation interface, and still requires an accepted entry
program. The decisive tests were coefficient separators and the weak TT2
boundary, rather than a larger scalar coverage interval.
