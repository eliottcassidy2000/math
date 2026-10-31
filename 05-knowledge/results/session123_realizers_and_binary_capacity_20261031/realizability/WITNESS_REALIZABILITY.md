# Realizers, witness sections, and a concrete finite-dimensional compiler

2026-10-31. **PROVED elementary deductions with classical categorical and
algorithmic inputs explicitly attributed.** No novelty or open-conjecture
closure is claimed by this packet.

## 1. The precise categorical distinction

Fix a recursive pairing code for natural numbers. A coded set X is regarded
as a partitioned assembly with its actual element code as the singleton set
of realizers. A map is tracked when a partial recursive function takes every
valid input code to the corresponding valid output code. Its tracker need
not halt on invalid codes. These are particular objects and morphisms in
Kleene realizability, not a declaration that our research repository itself
is a realizability topos.

The categorical input is standard. Faber–van Oosten, *More on Geometric
Morphisms between Realizability Toposes*, §1.2, defines assemblies and
tracking, identifies Ass(A) as a quasitopos, and identifies RT(A) as its
regular exact completion; RT(K1)=Eff. In particular the assembly category
is not itself being called a topos. Primary author copy:
https://webspace.science.uu.nl/~ooste110/realizability/geommorRTav.pdf .
Hyland's original *The Effective Topos*,
§7 (the PER representation following Proposition7.1) gives a full faithful
representation by partial equivalence relations and tracked maps; Lemma10.1
identifies maps N→N with recursive functions. We use only equality PERs on
coded subsets and the elementary tracking definition. Primary author copy:
https://www.dpmms.cam.ac.uk/~jmeh1/Research/Oldpapers/hyland-effectivetopos.pdf .

**Section theorem.** Let R⊆N×N be any relation and let π:R→N send(x,w) to x,
with R realized by paired codes. A tracked section of π exists if and only if
there is a total recursive f:N→N with R(x,f(x)) for every x.

**Proof.** From a section's tracker, compute the second projection of its
output on input x. Conversely f gives the tracker x↦〈x,f(x)〉. These maps
are inverse constructions on sections and selectors. No decidability of R
is used. The same statement holds over a coded domain X, with trackers only
required to halt on valid codes for X.

**Search corollary.** If R is decidable and every fiber is nonempty, then
π has a tracked section: test w=0,1,2,... until the first successful witness.
If only recursive enumerability of R is given, dovetail its enumeration
until a pair with first coordinate x appears; the same conclusion follows.
Neither corollary supplies a useful time bound. It is wrong to claim a
blanket choice obstruction for a uniformly decidable total witness relation.

**Classical-bijection hostile.** Let H(e) be the halting bit of machine e
on input e, defined classically, and let R_H={(e,H(e)):e∈N}. Projection is an
underlying bijection and a tracked morphism, but has no tracked inverse or
section: such a tracker would decide halting. Within the assembly category
it is both mono and epi (injectivity and surjectivity suffice), yet not an
isomorphism. It is **not a regular epi**: any monomorphism that is also a
regular epimorphism is an isomorphism. For completeness, if a monomorphism
m coequalizes a,b then a=b; the identity of its domain therefore coequalizes
a,b, and the coequalizer universal property supplies its inverse. This gives
a direct proof here, without requiring a projectivity theorem.

The standard full regular embedding into the exact completion preserves
monomorphisms and regular epimorphisms, not arbitrary epimorphisms. Its
image in Eff remains mono but is not epi: a mono/epi
in a topos is an isomorphism, while the full embedding would supply the
forbidden tracked inverse. Thus *external surjectivity is not internal
surjectivity*. This is a classical obstruction, not a computational test of
undecidability and not a claim that all epi notions coincide.

## 2. Exact safe-time witnesses as integer points

Let n>=1, positive integer speeds v_1,...,v_n, and rational radius r=p/q,
where p,q are positive integers and 2p<=q. Define

    G(v,r)={t∈[0,1]: ||v_i t||>=p/q for all i}.

**Compiler theorem.** G(v,r) is nonempty if and only if the following
integer system in exactly n variables z_i is feasible:

    0<=z_i<=v_i-1,
    q v_j z_i - q v_i z_j <= (q-p)v_i - p v_j   for all i,j.   (1)

For any feasible integer vector, put

    ell=max_i (q z_i+p)/(q v_i),
    h=min_i (q z_i+q-p)/(q v_i).

Then [ell,h]⊆G(v,r), including the possibility ell=h. Conversely every
safe t supplies z_i=floor(v_i t), satisfies(1), and belongs to its interval.
Thus the compiler preserves the *entire* witness fiber through the floor
vector, not merely nonemptiness.

**Proof.** A fixed z_i allows precisely the closed interval
[(qz_i+p)/(qv_i),(qz_i+q-p)/(qv_i)]. The given z-bounds put these intervals
inside(0,1). The n intervals intersect iff every lower endpoint is at most
every upper endpoint. Clearing positive denominators gives(1). Since
0<p/q<=1/2, an allowed interval lies strictly between z_i and z_i+1 in
its phase coordinate; hence its floor really is z_i. This proves every
direction, including equality and repeated speeds.

**Short certificate.** A nonempty fiber always has a witness

    t=(qz_i+p)/(qv_i)

for some i attaining the maximum lower endpoint. Its reduced denominator is
at most q·max(v). The pair(i,z_i), together with the instance, has polynomial
bit length even when n varies; exact verification against every speed is
polynomial in the binary input length. Accordingly the safe-instance decision
problem is in NP. This is not a claim that finding the witness is polynomial
when n is part of the input.

**Fixed-n selection.** For each fixed n, integer feasibility for(1) is
polynomial-time in the binary encoding of p,q and the speeds, by Lenstra's
fixed-dimension integer-programming theorem. A witness vector can be extracted
from a decision oracle by coordinatewise binary search within0..v_i-1,
using at most a polynomial number of calls and preserving dimension n.
Return its maximum lower endpoint. This is therefore a polynomial-time
safe-time selector on feasible instances, for every fixed runner count.

The only imported algorithmic theorem is Lenstra's fixed-dimension integer
feasibility result; the linear system and extraction above are direct
deductions. The original paper is H.W.Lenstra,Jr., *Integer Programming with
a Fixed Number of Variables*, Mathematics of Operations Research8(4)(1983),
538–548, DOI10.1287/moor.8.4.538. The author-hosted text and its introduction
were inspected: https://pub.math.leidenuniv.nl/~lenstrahw/PUBLICATIONS/1983i/art.pdf .
No new general-purpose polynomial-time integer-programming implementation
is supplied here. The executable reference compiler uses exhaustive search
only for bounded validation; it must not be advertised as Lenstra's algorithm.

## 3. What this actually adds to the current LRC result

Session122 proves that every eight-element C⊆{1,...,13} and any five positive
integer attachments have a safe time at1/14. Membership in this input family
is decidable. The compiler above therefore gives a total computable section
of its coded safe-witness projection. With n=13 fixed, it also gives a
polynomial-bit-time selector by the cited theorem. The output has a rational
certificate of denominator at most14·max(v).

This refines an existence statement into a precise tracked witness statement
with a bit-complexity corollary. It neither enlarges the proved family nor
improves the radius, and does not solve unrestricted LRC. A conductor-specific
selector can be far more practical than invoking fixed-dimensional ILP; the
present result is a robust uniform baseline, not a performance comparison.

The map is: common safe time t -> shared integer floor vector z; conversely
z -> full common interval -> its first endpoint. It preserves every original
speed inequality and closed boundary. It forgets a chosen interior point,
but retains its containing interval. The sidecar is the shared z-vector;
independent pairwise existential witnesses cannot replace it.

**Minimal jointness hostile.** At radius1/3, each two-element subset of
{1,2,3} is safe, but the whole set is not. In particular the pair{1,2} forces
t=1/3 or2/3, where speed3 is zero modulo1. This does not contradict(1): its
pair inequalities use one common vector z, whereas pairwise feasibility
allows inconsistent z-values in different pairs.

## 4. Collatz: a different direction of the witness map

The Session122 content graph maps actual valuation words to abstract paths;
not every abstract path has an actual word lift. Clock/length propagation
is a computable *necessary-condition* map, not a section of that abstraction.
Failure of every restored multiplier condition refutes all possible actual
lifts in the three near-tight layers. A surviving clock pair would not
produce a tracked lift: translation, ordered letters and first returns
remain missing. The section theorem identifies the exact extra obligation
without inventing it. No categorical completion creates those lifts for free.

## 5. Computational evidence and limits

`check_witness_compiler.py` independently compares the integer compiler,
a full-comb intersection algorithm, and literal endpoint verification over
its declared bounded universe. It includes the jointness hostile, point-only
closed boundaries, repeated speeds, scaling invariance and reversed labels.
All arithmetic is exact. The halting obstruction and Lenstra complexity are
analytic/cited facts, not claims supported by finite testing.

Normal and optimized Python runs both returned PASS with identical JSON:
1,125 instances, 73,404 integer assignments, 899 feasible instances,
114 point-only instances, 726 explicit matrix checks, and 12 scaling checks.
Reproduce from the packet directory with:

```sh
python3 -B realizability/check_witness_compiler.py
python3 -B -O realizability/check_witness_compiler.py
```

## 6. Explicit loss ledger and hostile obligations

| Source -> target | Map | Preserved predicate | Lost data / sidecar | Hostile test |
|---|---|---|---|---|
| Coded R -> coded N | witness projection | actual verified witness implies instance | witness erased; a section requires one uniform tracker | halting-bit graph is externally bijective and assembly-epi, not regular epi or Eff-epi |
| Safe time -> feasible floor vector | t -> floor(v_i t) | all original phase inequalities jointly | exact point erased; retain full common interval | pairwise-safe {1,2,3} at radius1/3 has no common safe point |
| Feasible floor vector -> rational safe time | maximum lower endpoint | every original inequality, including equality | other safe points omitted, but interval remains as certificate | point-only (1,3) at radius1/2 must not be dropped |
| Actual valuation word -> content/clock abstraction | forget order with declared clock data | necessary content and multiplier conditions | translation, order, first-return lift absent | an abstract path is not automatically an actual word |

The first and last rows are obligations, not claims that a missing section
exists. The middle rows supply actual computable witness maps with stated
complexity and are not categorical relabelings of existential data.
