# Finite seeds, affine boundaries, and the closure that remains

2026-10-05. **PROVED:** an explicit strict first-hit construction inside every
odd dyadic tail; the resulting whole-tail sufficiency theorem; exact forward
and backward fuel for the native L dependency; the finite-head reduction with
its closure hypothesis. **CITED:** the broader arithmetic-progression
sufficiency theorem and the inherited all-length one-member bound.
**FINITE-EXACT:** only the universes in section 6. **OPEN:** arbitrary-source
Collatz completion and any whole-tail completion proposed without its missing
closure proof. No literature-priority claim.

[Program](../../04-computation/experiments/finite_seed_bootstrap_boundary_20261005.py)
and [saved output](finite_seed_bootstrap_boundary_20261005.out).

## 1. Inheritance and the objects kept separate

The nearest proved mechanisms are:

* [THM-4512, coefficient-descent-classes-one-member](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md):
  an all-length threshold bound and at most one potentially unpaid least
  representative per contracting coarse cylinder. Its exact/coarse oddness-bit
  correction remains essential.
* [Frontier family compiler, sections 2 and 6](frontier_family_compiler_20261004.md):
  affine parameter thresholds and completed members built from supplied rooted
  hubs. The free parameter is not thereby universally rooted.
* [Partition cover, section 2](collatz_partition_cover_20261004.md):
  smaller common-future dependencies support well-founded induction only on a
  covered domain with an appropriate closure condition.
* [Paid guard budget, section 3](paid_guard_budget_20261005.md):
  the native L rule on 219 mod256, including its actual common-future words.

Our board is **point certificate; all-parameter law; finite head; closure;
binary guard; ternary fuel**. The anchor is the proposed finite-seed bootstrap.
The niche is the exact L-only ancestral boundary. The wildcard is the
constructive fact that even one whole dyadic tail is a sufficient test set.
The canonical hostile is 9: it belongs to 1 mod4 and descends to 7, which
does not belong to that domain. The corrected near miss is to promote a
descent theorem on a domain into its standalone home theorem. The useful
sidecar is the destination obligation, together with its induction rank.

Throughout, U(n)=oddpart(3n+1) on positive odd integers. A *root certificate*
records a finite actual route to the first occurrence of 1. A common-future
receipt records two actual routes to the same state and transports an
already supplied certificate. Neither object means an assumption that an
unexamined infinite family reaches 1.

## 2. A finite head is not the missing closure theorem

For a positive valuation word w of length j write

    F_w(n)=(Pn+B)/Q,  P=3^j, Q=2^A, B>0.

Assume Q>P. Its coarse odd source class is n=c+Qk, k>=0, where
c=(-B)P^(-1) modQ is the least positive representative. Put
e=(Pc+B)/Q. Earlier valuations are exact and the last may be larger;
the actual endpoint is oddpart(e+Pk), not necessarily e+Pk. Consequently

    k0=max(0, floor((e-c)/(Q-P))+1)

is the exact threshold for *nominal* descent. Every k>=k0 has actual descent.
The finitely many k<k0 are exceptions to this sufficient test, not asserted
failures of actual descent. THM-4512 gives k0<=1 for every contracting word;
this note uses that inherited result without reproducing its analytic bound.

**Finite-bank induction criterion.** Let D be a specified set of positive
odd integers and E a finite set of supplied certified bases. For each
n in D outside E, suppose a proved guarded rule gives a smaller positive
odd m with a checked common future with n. If each such m belongs to D or
has a supplied root certificate, then every n in D is certified by strong
induction on n. The same argument works with an explicitly well-founded
rank in place of size. For a finite bank of contracting coarse cylinders,
the possible threshold heads form a finite E; the criterion still requires
coverage of D and the destination condition.

The proof is direct: the smallest unrooted member of D cannot lie in E,
and its smaller destination is rooted by the stated alternative. The
common-future receipt then roots that member, a contradiction. A closed
cycle of *obligations* is not a supplied base certificate.

For D={n positive odd: n=1 mod4}, every n>1 satisfies
U(n)<=(3n+1)/4<n, and 1 is a certified base. Nevertheless 9->7 leaves D.
This does not prove D is unrooted; it proves these premises omit precisely
the closure needed by this induction. Section 3 shows why a claim to have
rooted this entire D would already solve the global problem.

The fixed-point boundary is also exact: nominal equality occurs at
n=B/(Q-P). If this is a positive odd source in the *exact* word class,
the word is an actual return. A nonroot example would be a nontrivial
positive cycle. Thus the inherited finite census of first coefficient
crossings with only the root exception must not be extrapolated without
an all-length proof. Extra final divisions on a coarse class may destroy
nominal equality as an actual return.

## 3. Every odd dyadic tail is a sufficient set, constructively

**Theorem.** Fix K>=1, an odd c with 0<c<2^K, and a finite cutoff H>=0.
For every positive odd hub u not divisible by 3, there are arbitrarily
large y>H with y=c mod2^K whose actual U route first reaches u after
exactly K odd steps. All its preterminal states can be required to exceed u.

**Construction.** Take K actual steps from c, allowing U(1)=1 solely while
forming this finite template. Let their valuations be w=(a1,...,aK),
with affine carrier (P,Q,B). Then P=3^K, Q=2^A, A>=K, and B is a
3-unit. Choose h>=1 in the congruence class

    Q·2^h·u = B mod3^K,
    y=(Q·2^h·u-B)/3^K.                         (1)

This is one exponent class modulo L=2·3^(K-1). To see existence
elementarily, v3(4^(3^j)-1)=j+1 follows by induction from cubing
1+3^(j+1)t with 3 not dividing t. Hence 2 has order 2·3^(K-1)
mod3^K and its powers comprise every unit. The program lifts one ternary
digit at a time. The factor Q in (1) is required.

Increasing h by L preserves integrality and makes y arbitrarily large.
For each proper prefix i=0,...,K-1, write its carrier (P_i,Q_i,B_i),
including (1,1,0) at i=0. It is enough to require

    y > T=max(H,u,max_i floor((Q_i·u-B_i)/P_i)).             (2)

This is an explicit finite height test, not an unbounded route search.

**Proof of the source and actual guards.** Equation (1) makes y an odd
integer. The original template at c gives Pc+B=0 modQ, so Py+B=0
modQ implies y=c modQ, and therefore y=c mod2^K. The nominal final
endpoint is 2^h u. Divisibility by the full denominator Q gives all earlier
integrality guards; each nonfinal nominal state is odd, since an even
state would make the next numerator 3x+1 odd and prevent its prescribed
positive division by 2. Thus the actual word is

    (a1,...,a_(K-1),aK+h),

and its endpoint is exactly u. The positive prefix carries and (2) ensure
every proper state exceeds u. In particular, a template which padded 1
does not export a padded root route. The output is a strict first-hit-u
receipt; if u has a supplied first-hit-root certificate, it splices safely.

**Sufficiency consequence.** Every positive integer x has a positive odd
forward state v. If 3 divides v, then U(v) is a 3-unit; otherwise use v.
Call the resulting state u. The construction gives a y in the prescribed
dyadic tail with U^K(y)=u. Hence x and y have a common future. If *every*
member of this tail reaches 1, every positive integer does. The converse
is immediate. This is a common-future statement; it does not assert that
the forward orbit of x visits that tail.

**CITED background.** Monks, Monks, Monks and Monks define sufficient sets
by common futures and give an admissible backtracing theorem for arbitrary
arithmetic progressions in Theorem 4.1 of
[Strongly sufficient sets and the distribution of arithmetic sequences in the 3x+1 graph](https://arxiv.org/pdf/1204.3904).
The proof above is a self-contained dyadic-tail specialization, also closely
related to the inherited frontier compiler. It does not use the paper's
printed full-product affine group claim, whose modulus95 counterexample
and repaired global subgroup are recorded in
[the guarded affine lift audit, section 2](collatz_affine_guarded_lifts_20261004.md).
No new theorem ID is assigned here.

**Finite-seed survivor.** A single actually certified 3-unit u roots every
member of the explicitly constructed exponent ray in (1). Thus a finite
seed library can produce infinitely many genuine completed members in
every dyadic tail. For fixed template and hub the sources grow
geometrically, and only O(log X) of this ray lie below X. It is not the
entire free-parameter cylinder, nor a coverage test for an arbitrary
supplied source. The example in the output is

    194179 --(1,4,15)--> 5,
    194179=3 mod8, 194179>10000.

The literal certificate 5->1 may be supplied separately.

## 4. Exact native L fuel and finite ancestral closure

The inherited dependency is

    L(n)=(9n-3)/16,  n=219 mod256, n>0.

It has a positive odd smaller child and a checked common-future receipt;
it is not an actual U edge. With z=7n+3 its conjugacy is

    7L(n)+3=(9/16)(7n+3).                         (3)

The native source guard is exactly v2(z)>=8. Each L step consumes four
binary digits and adds two ternary digits. Therefore the number of
consecutive native L steps from n is exactly

    max(0, floor((v2(7n+3)-4)/4)).                 (4)

When it is positive, the final shifted binary valuation is between 4 and
7 inclusive. For example 1755->987->555 has shifted valuations12,8,4.
Every step pays relative to its own source, so this is a valid finite
dependency chain; its terminal state remains an explicit root obligation.

Backward, the only possible parent of u is (16u+3)/9. It is a native
positive odd L parent iff

    v2(7u+3)>=4 and v3(7u+3)>=2.                 (5)

Indeed z_parent=16z_u/9. Divisibility by 9 makes the parent integral;
the binary condition is precisely its native guard. Positivity and
parent>u are automatic. Each inverse step reduces v3 by two and
increases v2 by four. Thus the exact ancestral depth is

    d_L(u)=floor(v3(7u+3)/2) if v2(7u+3)>=4,
           0 otherwise.                         (6)

For any finite seed set S, its closure under this one inverse rule is
finite, with cardinality at most sum_(u in S)(d_L(u)+1). Different seed
chains may overlap. This is an exact obstruction to an infinite L-only
bootstrap from finite point seeds, not an obstruction to refuelling,
switching controllers, or introducing another family law.

Literal certificates for S={1,123,555} give the exact L-only closure
{1,123,219,555,987,1755}. A finite verification of these bases legitimately
discharges that finite closure. It does not discharge the entire infinite
219 mod256 class. The state-dependent refuelling construction in the
parallel finite-seed tree work deliberately introduces a new operation
and lies outside this theorem.

## 5. Transfer ledger and the retained obligations

| Source -> target | Map and preserved predicate | Lost data / required sidecar | Decisive test |
|---|---|---|---|
| Contracting coarse cylinder -> smaller endpoint | oddpart(e+Pk), valid actual descent above k0 | Destination membership/root proof; last oddness bit | 9->7 leaves1mod4 |
| Certified hub -> dyadic-tail members | Equation(1), strict first-hit hub | Exponent phase, height guard, supplied hub certificate | Literal forward and backward readers |
| Entire dyadic tail -> all positive sources | Common-future construction through a 3-unit hub | No universal home premise has been supplied | One whole-tail home theorem would imply Collatz |
| Native L source -> smaller dependency | Shifted multiplier9/16 | Native guard and terminal root obligation | Exact binary fuel(4) |
| Finite rooted seeds -> L-only ancestors | Unique inverse and ternary fuel(6) | Other operations are not included | Exact finite closure above |

The positive conclusions are constructive completed subfamilies and exact
finite obligations. A standalone infinite-domain root theorem additionally
needs a closed well-founded cover or another proved all-parameter argument.

## 6. Exact reproduction and finite universe

From the repository root run:

    python 04-computation/experiments/finite_seed_bootstrap_boundary_20261005.py
    python -O 04-computation/experiments/finite_seed_bootstrap_boundary_20261005.py

The program uses integer arithmetic only and explicit checks unaffected by
optimization. Its finite controls are:

* all728 units modulo3^K for K=1,...,6: digit-lift logarithms versus an
  independently generated power table;
* all odd residues at binary depths1,...,5, hubs1,5,7,11,13, cutoffs0 and10000,
  and two successive exponent lifts:620 strict tail receipts, independently
  replayed both forward and backward; largest source335 bits;
* all10000 positive odd n<20000: literal native-L repetition and an
  independently iterated inverse-parent test against(4) and(6), plus actual
  common-future words for every applied L step;
* all719 contracting positive words of total valuation cost at most10:
  four coarse parameters each,2876 controls of the nominal threshold and
  actual oddpart endpoint; twelve words have a possible least-source
  exception to the nominal test. This is not an all-length first-crossing
  census or a cycle-exclusion proof;
* finite seed certificates, the explicit domain-closure hostile, the root
  fixed-point boundary, the coarse-versus-exact witness at1, and invalid
  integer/type and 3-divisible-hub inputs.

The all-height statements in sections2--4 rest on their proofs. The finite
universes test implementations and boundaries; their counts do not replace
those proofs or certify untested free parameters.
