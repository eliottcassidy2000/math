# Synchronized births and arithmetic limits of endpoint repair

**PROVED scoped mechanisms; FINITE-EXACT bank classifications. Chan–Pak OPEN.**
Local labels below are packet labels, not newly reserved global theorem IDs.
The inherited definitions and tail construction are in
[the previous packet](../joint_stability_forbidden_20261112/README.md) and
[the full run-tensor proof](../session288_twin_lengths_20261112/stable/PROOF.md).
Here H is the number of directed Hamilton paths, delta=H-W(O), and U is the
selected-start Hamilton/two-path-cover count. A unit seed has delta=1 and
produces H=1+2^t U for every integer t>=1. A stable-clone birth is strong,
has a stable unselected port, and has no stable unselected twin contraction.

## B1. A finite synchronized grammar giving infinitely many births

The 24 explicit records in `synchrony_kernels.json.gz` have two originally
unstable singleton ports a,b and an independent stable singleton port z.
For their positive transitive lengths the complete signed tensor is

    delta(a,b,z)=1+epsilon*(2^a-2^b), epsilon in {1,-1}.

This is an identity of full skeleton tensors, computed in two independent
ways, not an inference from the length-two and length-three discovery tests.
Hence its unit locus is exactly a=b, with z unrestricted. Twelve records
have each sign. Swapping the coupled ports reverses the sign.

Set a=b=m and z=1. Both growing blocks are unstable: increasing only one
of their lengths changes delta by plus or minus 2^m. Every old nonsingleton
unselected block has, when one member is duplicated, defect minus one equal
to either

    4+4*4^m,  or  12-12*2^m+8*4^m.

The first is positive. The second is 8x^2-12x+12 at x=2^m, and is positive
because its discriminant is negative and its leading coefficient positive.
Thus no nonsingleton block is stable. Positive inflation does not merge
distinct maximal marked transitive blocks; this is the inherited weighted
twin-quotient theorem. The free singleton z remains stable. Strongness is
preserved by substituting nonempty modules in a strong quotient. Therefore
every m>=1 produces a strong stable-clone birth of order 8+2m.

For ordinary marked isomorphism, compress the maximal unselected twin blocks
and retain their lengths and all selected singleton marks. Each of these
24 quotients is rigid. There are 12 pairs sharing a quotient. Equality of
total orders forces the same m in both charts; the exact coordinate equations
force m=1. Thus the 24 charts give 12 old order-ten births, and **24 distinct
births at every even order >=12**, with generating functions

    B(z)=12*z^10+24*z^12/(1-z^2),
    B_new(z)=24*z^12/(1-z^2).

The word grammar is finite: choose one of 24 starts at m=2 and repeat the
operation G that adds a twin simultaneously to each coupled block. A further
operation F grows the free port and preserves the seed; birth status above
is asserted only on the z=1 slice. This is a genuine infinite birth grammar,
with fixed quotients, not a complete grammar of all births or all counts.
Its order-twelve states are disjoint from the 834 births obtained by the
declared one-vertex extension grammar. All-length separation and literal
order10,12,14,16 controls are in `synchrony_grammar.json.gz`.

The compressed defect signature is the earlier two-port synchronization
gate multiplied by the free module's H. In arbitrary modules A,B,C it is

    H(C)*(H(A)H(B)+C2(A)H(B)-H(A)C2(B))

after choosing the positive orientation; C2 is the ordered two-path-cover
count. For transitive modules, H=1 and C2(TT_m)=2^m-2. Thus the bridge from
the earlier Farey-type gate to these new births preserves the signed defect,
but discards the frozen-block geometry and U. Those are indispensable sidecars
for birth status and actual Hamilton counts respectively.

One explicit record has rows
`[32,5,33,7,15,474,159,287,607,255]`, selected vertices {1,3},
coupled ports 6,9 and free port4. Its diagonal slope is

    U(m,m,z)=12*4^m+Q(m)*2^z,
    Q(m)=-50+100*2^m+100*4^m+100*6^m.

Here v2(Q)=1. For z>=2m+2, v2(U)=2m+2. Writing z=2m+1+k,
k>=1, its normalized odd coefficient is

    c_k(m)=A(m)*2^k+3,
    A(m)=-25+50*2^m+50*4^m+50*6^m.

Thus an extremely simple synchronized birth rule creates an **unbounded
lower-head invoice j=2m+2**. Finite growth rules do not imply a bounded number
of unpaid exponents. This is the first obstruction the repair language must
address. The displayed coefficient family is not claimed fully repaired.

## B2. Composing extension with synchronization gives46 odd-order families

Every one-vertex orientation, with either new marking, was tested on all24
synchronized order-twelve births. This gives88 new strong births of order13,
disjoint from all1414 in the plain extension layer. For each of the88, choose
one member of each of the two inherited coupled blocks:352 pair tests.
There are172 presentations preserving the diagonal unit identity. After
weighted-quotient deduplication, **46 distinct all-length birth families**
remain. Their complete two-port defects are

    delta(a,b)=1+lambda*(2^a-2^b), lambda in {2,-2,6,-6}.

For each, an actual stable unselected port p has a full three-variable
tensor whose diagonal in a,b is identically1 at every positive p length.
We do not require p to disappear from the off-diagonal tensor. This distinction
avoids discarding valid synchronized families with mixed free-port terms.

Increasing a or b separately at a=b=m changes delta by a nonzero
`+/-lambda*2^m`, so the growing blocks cannot contract stably. For every
remaining old nonsingleton block, the increment defect minus one has a
strictly positive expansion in F_r(m)=r!S(m,r). The exact certificate is
`u_r=sum_(b>=r) binom(b,r)*e_b`; all u_r>=0 and u_1>0. Thus it never
vanishes at positive m. Strongness and the free singleton port survive.
This proves birth status at every positive length, not only the tested ones.

Every retained quotient is rigid. Exact linear weight comparisons decide all
20 same-quotient chart pairs and show no positive-length intersections.
Literal graph-key controls agree at orders13,15,17. Therefore their generating
function is `46*z^13/(1-z^2)`. Together with B1 the new infinite grammar has

    (24*z^12+46*z^13)/(1-z^2),

in addition to its12 old order-ten starts:24 births at every even order>=12
and46 at every odd order>=13. The extension step changes the quotient before
the synchronized rule resumes. This is an explicit composition of two
previously separated operations. It remains a finite family of quotients,
not a universal birth grammar. `mixed_grammar.py` independently recompiles
each whole free-port tensor by endpoint convolution, separate from the signed
walk producer, and supplies the frozen-block and intersection certificates.

## B3. Five independent stable ports, and a local saturation boundary

The order13 tournament with rows
`[128,129,1027,1031,3119,3407,3359,3708,7839,3199,3,1039,3839]`,
selected vertices{2,3}, and unselected ports{4,9,10,11,12} is a strong birth.
Its full five-variable signed Stirling tensor is exactly `a_(1,1,1,1,1)=1`.
All five ports are jointly stable at all positive independent lengths.
The separate endpoint convolution agrees and retains108 positive U terms.
Nonuniform literal expansions of orders15,16,18 give defect1.

With this exact source and selection held fixed, adding one unselected vertex
has18 unit orientations. Fourteen keep the enlarged six-port set acyclic and
the graph strong; ten retain stability of every old port and the new one.
All ten merely place a twin immediately before or after one of the five
ports, hence are descendants. None gives a six-port birth in this declared
source-preserving move. Its transitive port word has letters
`[127,63,31,15,3]` relative to the retained eight-vertex exterior.

Across all eight five-port births in the retained order13 layer, both markings
and all one-vertex orientations give44 strong order14 births:16 have five
stable ports and28 have four. The16 weak source/sink paddings have six ports;
they must not be reported as strong rank-six births. No minimum-size or
global maximum-rank theorem is inferred from these restricted universes.

## V1. Complete valuation classification of the 478 old even higher-mode axes

Write U(a)=sum e_b*b^a in its unique finite exponential form. In this bank
the only odd bases are 1,3,5. Let g=min v2(e_b) over nonzero odd-mode
coefficients, and f(a)=sum_(b odd) (e_b/2^g)*b^a. Even modes vanish modulo
any prescribed 2-power once a is sufficiently large.

If f(a) is odd, v2(U(a))=g eventually. This gives 352 families. Otherwise,
the exact bank coefficients have odd coefficient of3 and even coefficient
of5. At a=2t+p, p=0 or1, write

    f(2t+p)=C1+C3*3^p*9^t+C5*5^p*25^t.

For every r>=0 and integer t>=0, the difference at t+2^r and t has exact
valuation r+3: v2(9^(2^r)-1)=r+3, whereas the 25 term has at least r+4
because C5 is even. Thus a parity branch either has fixed valuation1 or2
modulo8, or has one uniquely lifting zero at every successive precision.
The latter yields arbitrarily high v2(U), since one can choose an arbitrarily
large representative of each residue class, making all even modes vanish
at the requested precision. This proves the classification:

| Eventual behavior | Families |
|---|---:|
| Constant valuation |352|
| Two bounded parity branches |60|
| Unbounded valuation on a parity branch |66|

Explicit bounded cutoffs are at most5. There are 52,736 literal bounded
controls and compatible lifting controls through precision g+32, in addition
to the symbolic argument. Of the 66 unbounded families, 65 use odd bases1,3;
one also uses5. At birth519, port4, the full coefficients are
`[26,-860,1798,1698,240]`; its normalized odd part is
`13+899*3^a+120*5^a`. The even branch has the unique lifting zero.

These are classifications of a named finite collection of infinite axes.
They neither classify all future seeds nor supply their missing heads.

### V1a. All65 two-mode lifting roots are irrational2-adic integers

For each two-mode unbounded branch, the lifting argument gives a unique
t in Z2 with `9^t=-C1/(C3*3^p)`. If t=u/v is rational in lowest terms,
v is odd. Raising the equality to v gives an equality of rational numbers
`R^v=9^u`. Prime valuations imply v divides2u, hence v=1. The real sign
also forces R positive. Thus t is rational iff it is an integer, equivalently
`-C1/C3=3^(p+2t)` as a rational number. Exact rational factorization in
`lifting_roots.py` rejects that possibility for every one of the65 branches.
Their binary digit expansions are consequently not eventually periodic:
an eventually periodic expansion sums to a rational geometric series, and
conversely rational2-adic integers have eventually periodic binary expansions.

This concerns the odd-mode lifting roots, not a claim of randomness or an
automatic identification with another digit sequence. The single mixed3/5
branch is outside this rationality classification. All66 roots remain
computable by the uniform compatible one-bit lifting procedure.

## V2. Why a finite bank of affine transitive templates cannot close 126 families

An **affine TT chart** is a fixed finite tournament template with some vertices
replaced by transitive modules of lengths alpha_i*m+beta_i, alpha_i positive
integers, with all lengths positive for sufficiently large integer m. Fixed
modules may be arbitrary finite tournaments, absorbed into the fixed exterior.
Different charts may have different templates, offsets and slopes.

**Claim.** If the H-count of such a chart grows at most exponentially in m,
its unique exponential expansion has 2-integral coefficients at every odd
base. Here 2-integral means rational with odd denominator after reduction.

If the growing ports contain a directed triangle, take t vertices in each
of its three growing blocks, with t proportional to m. Cycling through the
three blocks allows independent permutations in each, producing at least
(t!)^3 Hamilton paths on that induced subtournament. Hamilton counts cannot
decrease when another vertex is added: insert that vertex in a Hamilton path
at its first permissible position, recovering the old path by deletion.
Thus a cyclic set of growing ports gives superexponential growth.

Otherwise the growing ports are acyclic. The inherited finite run expansion
has integer coefficients in products of F_r(alpha_i*m+beta_i), where

    F_r(l)=sum_(b=1..r) (-1)^(r-b)*binom(r,b)*b^l.

After substitution, each exponential base is a product of b_i^alpha_i and
each coefficient an integer times a product of b_i^beta_i. For an odd base
all b_i are odd, so even negative offsets yield only odd denominators.
Collecting terms preserves 2-integrality. Distinct positive exponential
bases are linearly independent, by Vandermonde.

Now fix an infinite arithmetic progression on which v2(U(a))=j. Its lowest
unpaid head is `1+U(a)/2^(j-1)`. If j>=g+2, some odd-mode coefficient has
negative 2-adic valuation. Restricting to a progression a=Lm+r changes the
bases to b^L and multiplies coefficients by odd b^r, so the defect persists.
Adding the constant1 cannot remove a negative valuation at base1 either.
Consequently no affine TT chart equals this head throughout that progression.

Nor can a **finite union** of charts, with an arbitrary choice of chart for
each parameter, cover the whole progression. Cyclic charts eventually exceed
the target. For an acyclic chart a nonzero difference of finite sums of
distinct positive real exponentials has finitely many real zeros: divide
by the smallest exponential and apply Rolle's theorem inductively. Unless
the profiles are identical, only finitely many integer matches are possible.
The identical case has just been excluded by 2-integrality.

Each of the 60 parity families has a branch j=g+2. Each of the 66 unbounded
families has a residue class modulo4 with j=g+3. The receipt supplies these
126 exact branches and 8,064 independent literal valuation controls. Thus
**126 of the 478 families cannot be fully repaired by any finite bank of
affine TT charts**, even allowing parameter-dependent chart selection.

This does not exclude finite recursive rules that change the core, nonlinear
clock lengths, nontransitive growing modules, or a richer observer. The
obstruction is to this precise representation class, not to Chan–Pak and
not to finite-rule grammars in general. The other 352 families are outside
this eventual obstruction; that is not a proof they admit such repairs.

## V3. A different obstruction to periodically repeating old seam words

The inherited 65-letter seam grammar acts on a three-coordinate state by
integer 3x3 matrices. The exact determinant list has absolute values1 or
at least5, and none equal6. Hence the multiplicative monoid generated by
those determinants does not contain6: two nonunit factors have absolute
product at least25, while a single one cannot be6.

Suppose a fixed word W is repeated m times between any fixed prefix and
suffix, and its scalar output is a sum with all three nonzero modes1,2,3.
Its minimal scalar recurrence has characteristic polynomial
`(x-1)(x-2)(x-3)`. Cayley–Hamilton bounds it by the characteristic polynomial
of the 3x3 word matrix; equality of degrees forces determinant6, a contradiction.

There is a stronger arithmetic-progression version. Among these determinants,
those supported only on primes2,3 satisfy v3>=3*v2. Products preserve this
cone. Therefore no word determinant can be 6^L for any positive integer L.
On a progression a=Lm+r the three modes become1,2^L,3^L, whose required
determinant is precisely6^L. Splitting a family into finitely many residue
branches cannot cure this obstruction. Profiles with at least four distinct
modes already exceed the three-dimensional recurrence bound.

In the declared 478-family bank, **470 have a lowest-head progression
excluded by these rank/determinant tests**; the other eight have only a
constant and a pure3 mode. This includes a repair scheme assigning one fixed
periodic word to each of finitely many arithmetic-progression branches.
It does not assert the elementary real-zero argument for arbitrary matrix
sequences, which may have negative or complex eigenvalues. This is independent
of V2: even wholly integral three-mode profiles meet it. Arbitrary nonperiodic
words and higher-dimensional recursive constructions are not covered.
`repair_obstructions.py` checks the pinned inherited atlas, rather than
assuming its matrices unimodular. The cone is a genuine bridge between
tournament seam composition and prime valuations: the map is determinant,
the preserved operation is multiplication, and the lost data are the actual
state cone and endpoint realizability. A permitted determinant is necessary,
never sufficient, for a physical repair.

## V4. A25-state parity automaton excludes periodic O/Z/J repairs for all478

The inherited physical root-change grammar acts on `(D,E,U,d)` by

    O: (D,0,2U,d),
    Z: (D-E+U,d,2U-d,d),
    J: (D-E+U,D-E,D+U,D).

This is [THM-4904, root-change word fibers](../../../01-canon/theorems/THM-4904-root-change-word-fibers-carmichael-chamber-parity-law-sqrt-density.md),
not a new constructor. Reduce its integer4x4 matrices modulo2 and retain
two flags: whether J has appeared, and whether O or Z has appeared.
The complete reachable semigroup has25 flagged states. Exhaustive transition
closure, independently reimplemented using bit-packed matrix multiplication
and polynomial determinants, gives

| Word class | States | Characteristic polynomial over F2 |
|---|---:|---|
| Empty |1| (x+1)^4 |
| Nonempty O/Z only |4| x^2*(x+1)^2 |
| J only |4| x*(x+1)^3 |
| Contains J and O or Z |16| x^3*(x+1) |

This is an all-word certificate: every generator transition remains in the
listed finite closure. It is not an extrapolation from a maximum word length.

**Modular eigenvalue multiplicity principle.** If an integer matrix M supports
a scalar sequence with distinct integer exponential bases b_1,...,b_k and
nonzero coefficients, Cayley–Hamilton and exponential independence imply
that every b_i is a root of its characteristic polynomial. The monic product
`product(x-b_i)` divides that integer polynomial. Reducing modulo a prime p,
the number of these bases congruent to r cannot exceed the multiplicity of
`x-r` in the reduced characteristic polynomial. This is a necessary spectral
condition; it does not establish physical state or target realizability.

A mixed O/Z/J period therefore cannot support two distinct odd integer
exponential bases: it has only one `(x+1)` factor modulo2. In particular it
cannot support both1 and3^L for any positive integer L.

For an O/Z-only period of length ell>=1, the exact characteristic polynomial
is `x*(x-1)^2*(x-2^ell)`. Indeed d is fixed, E is overwritten by0 or d,
U has multiplier2 at each step, and D retains its own coefficient1. Thus
this case also cannot supply3^L. A J-only period is a power of J, whose
characteristic polynomial is `x*(x^3-x^2-x-1)`. The cubic has one real root
rho>1 and a nonreal conjugate pair of modulus rho^(-1/2)<1 (its discriminant
is-44 and root product1). No nonzero root has modulus1, so no positive power
of J has eigenvalue1. The empty period is constant and cannot help.

Each of the478 even higher-mode families has an infinite exact-valuation
progression on which its lowest unpaid head retains nonzero bases1 and3^L.
For the constant class use its eventual valuation; for the other126 use the
explicit branches of V2. `audit_ozj.py` verifies both coefficients are nonzero
in every case. Consequently **no fixed repeated O/Z/J word, with any fixed
prefix, suffix and source interface, supplies that complete head progression**.
This excludes all478 in this precise periodic architecture, including the
eight pure3 families not excluded by V3. Assigning fixed periods to finitely
many arithmetic-progression branches does not evade the result. Arbitrary
adaptive words, new operations and general recursive grammars remain open.

A stronger tempting guess is false: J-containing words can have eigenvalue1.
The shortest examples have length7; `OOZOZOJ` has characteristic polynomial
`x^2*(x-1)*(x-82)`. The even second mode is exactly what the parity theorem
permits. All28,501 J-containing words through length9 were checked as a
discovery/hostile test;82 have eigenvalue1, none also has eigenvalue3 or9.
The25-state closure, rather than that finite word search, proves the theorem.

The three arithmetic filters are distinct: V2 constrains coefficients,
V3 constrains determinant valuations, and V4 constrains eigenvalue slots.
An ordinary TT run transfer can have eigenvalues1,2,3 and determinant6,
escaping V3/V4, but its integral readout still meets V2. A formal integer
matrix `[[1,4],[0,9]]` with initial vector(1,1) gives `(9^m+1)/2` in its
first coordinate, escaping V2's odd-coefficient restriction. No physical
tournament realization of that matrix is asserted here. It is a concrete
test interface for a richer repair language, not a proof by formal recurrence.

## G1. Exponential growth rank detects active strong components

Call a growing TT port active when it lies in a nontrivial strong component
of the fixed template. Singleton strong components contribute H=1 at every
transitive length and can be removed from the H calculation.

If active ports p,q lie in one strong component, orient p->q and take a
simple return path q,...,p. The walk `p,q,...,p,q` visits both twice and
each intervening vertex only once. Insert every remaining fixed vertex and
every missing port once into this directed walk. Insertion is always possible
by the same first-permissible-position argument used in V2. This produces a
positive skeleton term having at least two runs at both ports. If the ports
lie in different strong components, construct a two-visit cycle at each and
concatenate the completed component walks. The same positivity conclusion
holds. For one active port a two-visit cycle suffices.

Consequently two active ports of positive affine slopes give H>=c*4^m for
large m, and one active port of slope>=2 does too. Cyclic sets of growing
ports already have superexponential growth by V2. Therefore an affine TT
chart with dominant exponential base3 has exactly one active growing port,
of slope1. All other growing ports are H-invisible singleton components.

This turns four finite profile absences into all-order structural barriers.
The targets below are `1+U(a)` for the named old axes:

| Old birth | Original parameter | Profile in new parameter m>=1 | First H |
|---|---|---|---:|
|3, port0|a=m|3+9*2^m+4*3^m|33|
|3, port0|a=m+1|3+18*2^m+12*3^m|75|
|4, port0|a=m|-3+12*2^m+8*3^m|45|
|4, port0|a=m+1|-3+24*2^m+24*3^m|117|

Any chart defined for every m>=1 has nonnegative offset at its unique active
clock. Absorb that offset into a larger initial twin block. Write its count
as G times the count of the pointed active strong component; G is the fixed
product of all other nontrivial component counts. Integer one-port
exponential coefficients force G to divide the target coefficient content,
by uniqueness of the expansion. The possible divisors here are only1 or3.
The active component's first H is at most117. The inherited Busch bound
gives minimum125 at strong order10 and larger minima thereafter. Hence that
component has order at most9. The complete pointed bank through9 excludes
every required quotient profile. `small_head_barriers.py` checks all factors.
This use of the bound was checked against [Busch, Corollary1, p.4](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/)
on12 November2026; the bound itself is CITED, not a new result here.

Thus none of these four families has an affine TT chart starting at the
specified first parameter. Individual counts33,75,45,117 are all attainable.
Later-start charts, multiple charts with exceptional prefixes, and non-affine
recursive grammars are not excluded. This is a structural obstruction to a
family representation, not a forbidden-number claim.

## R1. Five complete repairs among the478, and one partial repair

The following odd coefficient families now have actual constructors for every
a>=1 and every s>=1, with target `1+c(a)*2^s`:

| Old birth / port | c(a) | Repair |
|---|---|---|
|122 /2| -3+12*2^a+150*3^a | Another stable core has exactly this odd slope. |
|334 /1| -3+20*2^a+390*3^a | Another stable core has exactly this odd slope. |
|741 /3| -5+24*2^a+390*3^a | Another odd-slope core for a>=2; c(1)=1213 uses the certified finite band. |
|14 /1| -7+32*2^a+18*3^a | Separate pointed H-profile pays s=1 for a>=2; c(1)=111 uses the finite band. |
|42 /5| -25+25*2^(a-1)+60*3^a for a>=2; c(1)=45 | Separate pointed H-profile pays s=1 for a>=2; finite-band seed at a=1. |

For the last two, the original unit seed pays every s>=2. For the first
three, a genuine odd-slope unit seed pays every s>=1. The finite exceptional
parameters have actual endpoint certificates and any needed finite head
graphs, not an assumption that a small value must be present. The complete
repair kernels are independently recounted by endpoint convolution;164
whole constructed graphs are checked at positive parameter pairs whose
orders are at most18. `repaired_rays.py --family 14 --parameter 4 --exponent 1`
emits the actual adjacency rows of an example.

For birth4/port0, c(a)=-1+3*2^a+2*3^a. The separate6-vertex profile
`-1+6*2^a+4*3^a` pays s=1. The original seed pays s>=3. The s=2 family
remains unresolved, with the early-start chart obstruction in G1.

Thus **five of the478 are completely repaired;473 remain unclosed here**.
Their valuation behaviors are all classified. Classification and repair
are different deliverables; neither is silently substituted for the other.
