# Target height, transitive-module rigidity, and physical Jordan cancellation

**Session282. PROVED, independently reviewed.** Chan–Pak remains
OPEN. The underlying substitution/Stirling and ordinal-product laws are inherited;
new applications concern exact target height, signed-correction rigidity, and
explicit mixed-tape cancellations. No new all-even-gap existence claim is made.

## 1. Inheritance and connection

The closest positive mechanism is Session281's finite continuation experiment
basis. Its hostile is a scalar equality losing the receiver used by later fusion.
Its congruence-complete but power-sparse language is the corrected near miss.
The least-used sidecar here is the number of visits to a transitive module.

Sources: [THM-3121 substitution kernel](../../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md),
[THM-3235 transitive blowup law](../../../../01-canon/theorems/THM-3235-transitive-tournament-blowup-saddle-and-hypergeometric-decimation.md),
[THM-3166 falling-factorial transform](../../../../01-canon/theorems/THM-3166-falling-factorial-order-join-path-colour-transform.md),
[THM-5123 Jordan memory](../../../../01-canon/theorems/THM-5123-common-gap-flag-and-periodic-apex-jordan-bounds.md),
and [THM-1862 SCC product](../../../../01-canon/theorems/THM-1862-order-join-reduction-principle.md).
No novelty is claimed for those general transforms or elementary Vandermonde
algebra. The map is Hamilton path -> module-run skeleton. It preserves the run
assignment and count but its scalar readout loses root interfaces. The live
board is height / runs / cancellation / growing interfaces / physicality / receiver.

## 2. Sharp module growth floor and exact finite target search

Fix tournament Q; replace vertex i by nonempty transitive TT_(m_i), inheriting
all external arcs. Write Q[m]. For a strong component C of size at least3,
let c(C) count directed Hamilton cycles up to cyclic rotation. Camion's theorem
gives c(C)>=1.

**Theorem 1.**

    H(C[m]) >= H(C)+c(C) sum_(i in C)(2^(m_i)-2).       (1)

Proof. Every Hamilton path of C yields one path with each module in one
increasing run. Cut any Hamilton cycle at i to get skeleton i,...,i. Allocate
an arbitrary nonempty proper subset of module i to its first run, its complement
to the second, both increasing. This gives2^(m_i)-2 paths. The repeated label,
cycle and subset distinguish them from each other and the one-run paths. ∎

Since H(C)>=3 and the other SCC factors are positive integers,

    H(Q[m]) >= 2^(m_i)+1                              (2)

for every cyclic quotient vertex i. Equality holds for C3[TT_m,1,1], whose
count is2^m+1. At m=1 the two-run family is empty and the bound is3; zero-size
modules are excluded. Singleton-SCC transitive modules are ordinal factors
with H=1 and are count-neutral.

Thus to search exactly for H<=X, X>=3, in a fixed template, set neutral
module sizes to1 and restrict every cyclic size to

    1<=m_i<=L=floor(log2(X-1)).                        (3)

Use integer bit length at boundaries. For X=1 only acyclic templates contribute1.
The fixed q-vertex template has at most L^q output values below X; a finite
bank has O((log X)^qmax) support and cannot cover an odd tail. Even all labeled
quotients of orders at most q give at most

    sum_(j=1)^q 2^(j(j-1)/2) L^j                      (4)

values. Duplication only weakens this bound. Coverage via transitive templates
must increase quotient complexity, not merely inflate a finite collection.

## 3. Positive compiler and all-size identity certificates

A legal skeleton is a word of quotient labels, adjacent labels distinct and
following Q's arcs. Let a_r count skeletons with r_i visits to module i.
Allocation of its labeled vertices to r_i ordered nonempty increasing runs
has r_i! S(m_i,r_i) possibilities. The inherited run formula is

    H(Q[m])=sum_r a_r product_i r_i! S(m_i,r_i).       (5)

All terms are nonnegative; for fixed m the sum is finite. Together with(3)
this is a complete absolute-target oracle for a fixed template bank, not a
claim of fast complexity or completeness of the bank.

Now suppose k variable modules induce an acyclic quotient and s other
vertices are fixed. Removing the outside vertices from a path leaves at most
s+1 segments. Within any segment the variable labels follow a topological
order, so each is visited at most once. Therefore r_i<=s+1 for every variable
module. A finite skeleton tensor gives(5) for ALL positive module sizes.

For m>=1,

    r!S(m,r)=sum_(j=1)^r (-1)^(r-j) binom(r,j)j^m.

Along affine schedules m_i=a_i n+b_i, with nonnegative integer slopes/offsets
and all sizes positive, group equal bases B=product_i j_i^(a_i) to obtain

    H_n=sum_B c_B B^n,    c_B integers, B positive integers. (6)

Negative offsets can be handled by shifting to a valid initial n. Zero slopes
are allowed. Coefficients may have signs, although(5) is positive. Equal bases
MUST be grouped. If two families' combined base list has size d, agreement at
any d consecutive valid indices proves agreement at every valid index: the
matrix is a nonzero diagonal times the Vandermonde matrix on those bases.
Subtracting a constant at base1 similarly certifies constant gates. This is
finite proof from a known annihilator, not an unsupported fit.

## 4. Signed compensation rigidity and its boundary

**Theorem 2.** Any finite integer signed combination D_n of the preceding
acyclic-variable, affine-size families is either constant or eventually has
absolute magnitude at least exponential. It cannot be a nonconstant polynomial
or a nonconstant bounded/subexponential correction indefinitely.

Proof. After grouping bases, the largest B>1 with nonzero coefficient dominates:
|D_n|>=|c_B|B^n/2 eventually. If none remains, only base1 survives. The same
argument rules out equality with a nonconstant polynomial at infinitely many
unbounded indices. ∎

Nonzero constant differences are allowed. Changing templates with n, cyclic
variable modules, or operations changing external neighborhoods lie outside
this theorem. The next construction is an actual escape, not a contradiction.

## 5. Mixed tapes isolate a physical linear error

Let T(w) be a transitive tape in word order plus an apex, with bit1 meaning
the apex beats that vertex. For nonempty w of length L, directly

 H(T(w))=[last=0]+[first=1]+#adjacent01
            +sum_(j<i,w_j=1,w_i=0)2^(i-j-1).          (7)

Indeed a Hamilton path is an increasing subset before the apex and its
increasing complement after it. If max(before)=i<j=min(after), only adjacent
01 works; if j<i, the intervening positions are freely assigned. Empty before
or after gives the first two terms. For empty w, H=1 separately.

For w=(0^a1^b)^m, a,b,m>=1, R=2^(a+b), D=(2^a-1)(2^b-1), equation(7) gives

 H=m+D sum_(d=1)^(m-1)(m-d)R^(d-1)
  =D(R^m-1)/(R-1)^2+(1-D/(R-1))m.                    (8)

For a=b=s,

 (2^s+1)^2 H(T((0^s1^s)^m))
       =4^(sm)-1+2(2^s+1)m.                         (9)

The multiplier is itself a physical count: two ordinal factors
C3[TT_s,1,1], each having count2^s+1. Hence

    9 H(T((01)^(2m)))=16^m+12m-1,
   25 H(T((0011)^m))=16^m+10m-1.                     (10)

Their difference is2m for every m>=1. The actual product tournaments can
have equal order4m+9 by adding two neutral sources to the first. For m=1
the counts are27,25; for m=2 they are279,275. At m=0 the identity fails:
the singleton boundary gives9-25=-16.

All-even-gap existence is ALREADY known from Sessions199/221/243/244; no
novel existence claim is made. This gives a different explicit signed-mode
extractor, at WORSE exponential height than the known linear-cost pairs.
The mixed tape retains a Jordan/linear term after its exponential modes cancel.
It escapes Theorem2 because zero and one positions are interleaved, not
homogeneous transitive modules of a fixed quotient.

At s=3 there is also the exact positive identity

    H(T((01)^(3m)))=9H(T((000111)^m)), m>=1.           (11)

It does not identify the root interfaces: the apex deletion on the left has
count1, while a tape-apex root in the factor9 ordinal construction has deletion9.

## 6. Full-interface resource and the necessary two-cover distinction

For the balanced tape of(9), its rooted carrier is

 h=1, Q=4^(sm), U=V=1+(4^(sm)-1)/(2^s+1), d=e=0.      (12)

The deletion is transitive. U is1 plus sum2^(L-1-i) over outgoing tape positions;
V is1 plus sum2^i over incoming positions. Alternating-length versions vanish
because last bit is1 and first bit is0. These geometric sums prove(12).

Write F_s(m) for the two-cover count of the FULL tape, not its root deletion.
The two needed exact laws are

 F_1(k)=(8/25)9^k+((4/15)k+178/225)4^k+8/9,
 F_3(m)=(128/1225)729^m
          +((148/315)m+31282/99225)64^m+128/81.        (13)

Finite proof certificate: the apex occupies either ordered cover slot. Its
slot has a before and after subsequence; the other slot is increasing. Track
(b,a), with b the last before bit or2 if empty, and a whether after has begun.
For tape bit t, the three choices are identity (other slot), (b,a)->(t,a)
(before), or (b,a)->(b,1) (after, allowed iff a=1 or t=1). Start(2,0), accept
b!=1, multiply by2. This six-state matrix for period0^s1^s, s=1,3, is annihilated
by (x-9^s)(x-4^s)^3(x-1)^2, as exact matrix multiplication verifies.
The first six F_1 values are2,8,48,336,2576,21072; the F_3 values are
2,128,60672,40933632,29547824384,21516383326464. Formula(13) obeys the same
recurrence and six values, hence holds for all m>=0.

Put P=T((01)^(3m)), Q=T((000111)^m), R=64^m, m>=1. The equal-order deletion
objects D_L=TT6 -> P and D_R=C3 -> C3 -> Q both have6m+7 vertices and count
H=(R-1)/9+2m. By the inherited path-colour law at2, full two-cover counts
multiply across ordinal sums. TT6 has64, C3 has12; consequently
Q(D_L)=64F_1(3m), Q(D_R)=144F_3(m).

Add a universal source root to each. The two actual6m+8-vertex carriers are

 X_L:(H,h,Q,U,V,d,e)=(H,H,64F_1(3m),64F_1(3m),H,0,H),
 X_R:(H,h,Q,U,V,d,e)=(H,H,144F_3(m),144F_3(m),H,0,H).  (14)

Here the beginning mass is the full deletion two-cover count, ending and
end-sign masses are h, and start-sign mass is the inherited zero signed-two-cover
identity for a nonempty tournament. Thus these are physical same-H/same-receiver
twins with sending contrast

 DeltaQ=DeltaU=-(6656/1225)729^m
       +((576/35)m-19232/3675)64^m+512/3.              (15)

This is negative for all m>=1: (729/64)^m>=11m and512/(3*64^m)<=8/3, giving
an already negative linear upper bound after division by64^m. At m=1 the
common H is9 and Q/U are21504 versus18432 (ratio7/6).

These are not interchangeable in every context: a right fusion can detect U.
No unit digit or all-target selector is claimed. The source-root and ordinal
laws are inherited; the resource is the explicit parametric same-order family.
An early draft confused root-deletion Q=R with the full-tournament F_s; literal
profile testing caught this before publication. Equations(13)–(15) retain the
correct coordinate and replace the erroneous64R/144R claim.

## 7. What remains

The exact target oracle bounds fixed-template searches. The finite identity
certificate can synthesize parametric constant gates. Mixed tapes provide
physical signed polynomial corrections outside the rigidity class, and scalar
twins with a controlled continuing-interface contrast. The missing step is
turning such signed relations into one lawful absolute-target selector.
Subtracting counts of two tournaments is not itself a tournament operation.
No result here fills an arbitrary odd interval or completes Chan–Pak.

`verify.py` compares literal graph counts with the skeleton compiler, cycle
floors, tape formulas, exact matrix annihilators and the singleton hostile.
Finite controls accompany the all-order arguments rather than replacing them.
