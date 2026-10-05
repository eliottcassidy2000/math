# Paid Collatz controllers from cyclic words and quadratic portraits

2026-10-04. **PROVED:** a uniform family of paid controllers with two
unbounded parameters, their exact coverage gain against specified banks,
repairs of infinitely many ternary-tower exponent classes, and the finite
portrait boundary below. **INHERITED PROVED:** the proper graph rank,
affine group, prime-adic repeat fuel and anchor-switch boundary.
**FINITE-EXACT:** the declared audits. **OPEN:** universal positive
coverage and complete negative-basin classification. No literature
priority claim is made.

The main gain is an infinite controller family based on repeated words
`1^h 2`. Its negative rational anchor records exactly how many repetitions
are possible; a checked exit then pays the original integer and therefore
the inherited rank. These families add positive-density coverage to the
named binary16/ternary banks, including infinite families of powers of
three that the concurrent early-reroute grammar cannot handle.

## 1 Inheritance and the active concepts

Closest mechanism: [branch tolls and the proper graph rank](collatz_branch_toll_rank_20261004.md).
Canonical hostiles: 7,27,703; a numerical increase after an unpaid anchor
switch; and the negative cycles. Corrected near miss: a finite critical
orbit is not a theorem of global attraction. Least-used sidecars: the odd
cofactor at a precision refill, and a rational anchor determined by a word.

Anchor: paid coverage with unbounded repeat counts. Niche: finite quadratic
portraits and polynomial attempts to combine anchors. Wildcard: what
Burnside's theorem does and does not transfer to the affine word group.
The board is **original rank / guarded word / anchor / binary fuel /
refill cofactor / ternary address / covered domain**.

The old [row and braid audit, section5](collatz_mod6_20260917_row_braid_typing.md)
already proves the three rational PCF polynomial parameters and their
integer portraits. The [affine blueprint audit, section2](collatz_blueprint_20260921_affine.md)
already proves the unrestricted metabelian group. These are recovered
mechanisms, not new discoveries.

While this session was active, incoming commits `30a363ec64`, `d50814bcf9`,
and `81057f0bdf` supplied [certificate memory](collatz_guarded_pumping_memory_20261004.md),
[guarded finite affine lifts](collatz_affine_guarded_lifts_20261004.md), and
[early reroutes](collatz_early_reroute_20261004.md). The first already proves
the two-anchor ultrametric boundary developed independently below; the
second supplies the exact finite group and its legality restoration.
The new paid word families, their density comparison, and the ternary-tower
repairs are the principal additions of this note.

## 2 Exact run lengths give a controller an unbounded register

Use the signed odd map `U(n)=oddpart(3n+1)`. A positive valuation word w
of length r and total halving cost A has formal action

    F_w(n)=(P*n+B)/Q,   P=3^r, Q=2^A,
    rho_w=B/(Q-P).

Its carry B is odd, and Q-P is odd and coprime to3. Thus rho is an odd
2-adic rational. The exact word cylinder is

    n=rho_w mod2^(A+1).                              (1)

Indeed its endpoint must be odd, so `P*n+B=Q mod2Q`; substituting the
fixed-point equation gives (1). Repeating the word keeps its anchor and
multiplies A by the repeat count. Therefore, for an odd integer n!=rho,
the exact maximum number of complete consecutive copies is

    q=floor((v2(n-rho)-1)/A).                        (2)

After those copies the state is

    x=rho+(P/Q)^q*(n-rho),
    v2(x-rho)=v2(n-rho)-q*A in{1,...,A}.              (3)

The next copy fails within its r letters. This is the exact-odd-word
specialization of the inherited binary repeat fuel. At n=rho, when that
anchor is an integer, the correct object is a signed cycle, not an
indefinitely pumpable first-hit root certificate.

The controller description is finite, but its precision and repeat-count
registers are unbounded. It is outside the scope of the incoming theorem
that a sound regular first-hit word language has bounded odd depth.

## 3 A uniform paid exit theorem

Suppose w has negative rational anchor
`rho=-b/d`, in lowest terms, b,d positive odd, and its first valuation is1.
Its rational periodic orbit has exact valuations by (1); all actual
controller inputs and outputs below remain integers. Then `b/d=1 mod4`.
Write r=len(w), A=sum(w), and suppose `3^r>2^A`.
Fix q>=1. On the source stratum

    d*n+b=2^(A*q+1)*t,   t positive odd,              (4)

q complete copies of w are actual, and

    x=F_w^q(n)=(2*3^(r*q)*t-b)/d.

Here n is an integer; the congruence on t imposed by d is retained.
The next actual valuation a is at least2. Define

    a0(r,A,q)=min{a>=2:2^(A*q+a)>3^(r*q+1)},
    L(r,A,q)=a0(r,A,q)+1.                            (5)

**PROVED paid exit.** If `2^(A*q)>=b` and the actual exit has a>=L, its
endpoint y satisfies `0<y<n` and `R(y)<R(n)` for the inherited root rank.

**Proof.** Set c=b/d and `alpha=3^(r*q+1)/2^(A*q+a)<1/2`. Then

    y=alpha*(n+c)-(3*c-1)/2^a.

The size premise in (4) gives n>=c, and c>1/3 since the anchor's first
signed step is negative. Hence y<n. Positivity follows from the actual
positive orbit. Also n=3 mod4, so the original energy is `3(n-1)^2/4`.
For every positive odd y, its energy is at most `3(y-1)^2/4` (with ROOT1
handled separately). Thus strict integer descent pays the original rank.
No child certificate or convergence assumption entered this proof.

For fixed w,q this is one dyadic source cylinder: the q exact copies
followed by a final valuation at least L. Its modulus is `2^(A*q+L)`,
and its density among odd integers is

    2^(-(A*q+L-1)).                                 (6)

The precise terminal exponent need not be fixed. Every legal higher
valuation only decreases the child further. A sharper compiler can use
a0 and check the least source against its child, discarding a finite
initial segment if needed. The saved audit does this for576 exact-exit
families from the three integer cycles; all tested least sources pay.
That finite observation is not promoted to an all-q claim at a0.

## 4 An infinite rational anchor atlas with an all-height gain

For h>=1, choose

    w_h=(1^h,2),   r=h+1, A=h+2,
    rho_h=-(3^(h+1)-2^(h+1))/(3^(h+1)-2^(h+2)).      (7)

These are growing words. The first anchors are `-5,-19/11,-65/49`.
Their reduced numerators b satisfy `b<3^(h+1)<2^(2A)`. Consequently
the paid exit theorem applies for **every h>=1 and q>=2**.

In binary distance, `v2(rho_h+1)=h+1`; in the real metric rho_h also
approaches -1. These two convergence statements describe the same
explicit rational atlas, not attraction of positive Collatz orbits.
The integer remembers which anchor and how much precision it shares
with that anchor.

Adjoin the genuine -17 cycle word `(1,1,1,2,1,1,4)`, with r=7,A=11,
and use every q>=1. Its height premise holds already at q=1.
Call the union of these safe paid cylinders B_new.

**PROVED disjointness.** Different h are distinguished by the initial
run of ones. At fixed h, the first departure from w_h distinguishes q.
The -17 word has prefix `1112 114`; at the same first-run length h=3,
two copies of w_h instead begin `1112 1112`. These banks are disjoint.

They are also disjoint from the sixteen rows in the inherited binary-debt
bank. Their first reset is2, followed by a1, so its debt index must be e=1.
The two e=1 rows require the suffixes `(1,6)` or `(6)` after that reset.
The new words have `(1,2)` for h=1, `(1,1)` for h>=2, and `(1,1)` for
the -17 word. Thus neither old row can apply. Immediate descent and the
initial reset-at-least-three rule also fail there.

Every source in B_new has n=3 mod4 and n!=11 mod32. Filtering out
n=2 mod3 places it in the inherited critical set. This is a genuine
infinite collection of paid critical cancellations.

The exact odd-relative density is the convergent sum

    delta_B = sum_(h>=1,q>=2) 2^(-((h+2)q+L(h+1,h+2,q)-1))
            +sum_(q>=1) 2^(-(11q+L(7,11,q)-1)).       (8)

It is approximately0.00460266542236. The script retains exact rational
bounds, summing h<=24,q<=24 and the -17 terms through q=64. Since L>=3,
the omitted q tail at a fixed A is bounded by
`2^(-A(Q+1)-2)/(1-2^(-A))`. The omitted h tail is bounded by
`(8/7)*2^(-2(H+1)-6)/(1-1/4)`.

These are natural-density bounds. For large h, the source has many
initial ones: every h>H source lies in `n=-1 mod2^(H+2)`. For large
q at a fixed h, it lies in a single deep anchor cylinder. The number
of applicable h,q below a source-height bound is at most quadratic
in log X, by (4) and the numerator bound. Truncate those outer cylinders
first; their density tends to zero. Finite cylinder counts and the
vanishing tails justify (8), not countable additivity alone.

**Named bank comparisons.** Let D2 be the inherited binary residual
density, and delta_T its old ternary sibling coverage. CRT and the
two banks' height-controlled tails give

    new residual density = (D2-delta_B)*(1-delta_T).

It decreases from approximately0.15956213718787 to0.15655262213373.
The newly covered portion outside that named old bank is approximately
0.00300951505414 of all positive odd integers.

The concurrent early-reroute note improves the **original-source**
ternary bank to density between0.4458180914749391 and0.4458180914749392.
Replacing delta_T by that stronger density gives a second exact product
interval in the saved JSON: its residual is approximately0.13268614492279,
with approximately0.00255071390807 newly paid outside that strengthened
original-source comparison bank. Position-dependent reroutes and the other
repository constructions are not exhausted by either comparison.
These are densities of paid dependencies, not completed home proofs.

## 5 Paying some ternary towers that the early reroute cannot enter

The incoming [early-reroute theorem, section6](collatz_early_reroute_20261004.md)
excludes every source3^a from its entire adaptive initial-ones grammar.
After its other elementary dispatch rules, the surviving exponent class
is a=3 mod8. Our h=1 controllers enter that remaining class.

For h=1 and each q>=2, the safe source residue has modulus
`2^m`, m=3q+L(2,3,q), and is3 mod8. The cyclic subgroup generated by3
modulo2^m consists exactly of the units that are1 or3 mod8, and has
order `2^(m-2)`. The order follows from
`v2(3^(2u)-1)=3+v2(u)`. Thus there is a unique exponent residue p_q
modulo2^(m-2) taking3^p_q into that paid cylinder. Lifting one binary
digit at a time computes it without expanding the source integer.

**PROVED entire exponent families:**

| q | Required final valuation at least L | Exponents a of the source3^a |
|---:|---:|---|
| 2 | 3 | 107+128t |
| 3 | 4 | 11+2048t |
| 4 | 4 | 5387+16384t |

Here t>=0. Their actual source words are `(12)^q` followed by the
indicated reset. For example

    3^11=177147 --(1,2,1,2,1,2,4)-->47293 <177147.

Every parameter in each exponent progression pays the original rank,
by the same source-cylinder proof. The q-family domains are disjoint.
All their exponents are3 mod8: for q>=2 the source is -5 mod128,
whose power-three exponent is11 mod32. The union over every q>=2
covers approximately **0.06696427719932 of the exponent class3 mod8**.
This is a density among exponents, not among positive integers. Its
exact relative density is16 times the h=1 contribution to (8).

This is a reciprocal prime-clock construction. Earlier ternary source
guards locate Mersenne exponents using powers of2 modulo3^d. Here binary
controller guards locate ternary-tower exponents using powers of3 modulo2^m.
The modular address chooses a legal paid family; the actual word and
strict endpoint inequality supply the convergence-relevant content.
Neither family covers its entire exponent class.

## 6 What switching anchors really does to the remaining fuel

The inherited two-anchor ultrametric law extends through an affine word
without approximation. Let F(x)=(3^r*x+B)/2^A, and choose odd-denominator
anchors rho,sigma. Put

    K=v2(n-rho),  delta=F(rho)-sigma,
    h=v2(delta),  Knew=v2(F(n)-sigma).

If delta=0, Knew=K-A. Otherwise, when K-A!=h,

    Knew=min(K-A,h).                                (9)

When K-A=h, the two normalized odd summands cancel and Knew>h,
possibly without a uniform upper bound. Thus every strict fuel increase
passes the precise source shell `K=A+h`. This locates the obstruction;
it does not bound the new fuel.

An actual positive example makes the missing coordinate visible:

    n=32t-5 --(1,2)-->x=36t-5,  t positive odd,
    v2(n+5)=5,
    v2(x+1)=2+v2(9t-1).                             (10)

The final precision can be arbitrarily large, on explicit dyadic
progressions of t. At t=1 this is `27 ->41 ->31`, and the new precision
is5. The new odd cofactor is `oddpart(9t-1)`. That expression is a
generalized Syracuse operation on a **coordinate**, not a conjugacy of
the full original dynamics to9x-1 and not a free child certificate.

This is why the controller retains the cofactor as well as the fuel.
An unbounded repeat may be compressed; an anchor refill must still
record which arithmetic expression supplied its new digits and how
the complete route pays the original source.

## 7 Quadratic portraits suggest a construction and expose its limit

The user's finite graphs are exactly the inherited integer preperiodic
portraits of x^2, x^2-1, x^2-2. Their polynomial encodings are

    P_0(x)=P_1(x)=x(x^2-1),
    P_2(x)=x(x^2-1)(x^2-4).

Writing f_c(x)=x^2-c, direct factorization gives

    P_0(f_0(x))=P_0(x)*x(x^2+1),
    P_1(f_1(x))=P_1(x)*x(x^2-2),
    P_2(f_2(x))=P_2(x)*x(x^2-2)(x^2-3).             (11)

These retain the full finite root sets and their pullbacks. They suggest
trying a polynomial of several Collatz anchors in place of one distance.

**PROVED limitation of that specific repair.** Suppose a nonempty
finite root set at each controller mode is transported into the next
root set by **every internal affine edge**, with nonempty forward words.
On a directed cycle the
composite has positive slope3^r/2^A!=1 and acts injectively on its
finite starting set. It must permute that set. But its only periodic
point, even over the complex numbers, is its unique fixed point.
Therefore that set has just one element. Strong connectivity propagates
the singleton conclusion to every mode.

In particular, polynomial identities
`P_target(F(x))=constant*P_source(x)` on infinite guards reduce each
polynomial in such a cyclic component to a power of one linear factor.
They cannot combine incompatible anchors into a common polynomial
fuel. Quadratic maps evade the finite-set argument because they are
not injective: several roots can merge.

The scope matters. This excludes a common polynomial transport identity,
not arbitrary guarded multi-anchor controllers or the actual signed
cycle portraits. The repair pursued here is a finitely described
**infinite rational anchor atlas**, exact valuation guards, and retained
cofactors. The incoming [quadratic escape atlas](quadratic_escape_rank_atlas_20261004.md)
also shows why finite critical portraits do not imply attraction: the
integer points outside these finite quadratic cores escape.

## 8 Burnside provides a decomposition principle not integer coverage

**CITED.** Burnside's theorem says every finite group of order p^a q^b
is solvable. In particular a nonabelian finite simple group needs at
least three prime divisors in its order. Cyclic prime-order simple
groups remain possible. See [Etingof and coauthors, Theorem4.20 and
its proof](https://math.mit.edu/~etingof/replect.pdf).

For Collatz the relevant recovered group is already explicitly

    Z[1/6] semidirect <2,3>,

with translations forming the commutator subgroup and commuting
multiplier exponents in the quotient. For example D(x)=2x and T(x)=3x+1
have commutator D T D^-1 T^-1 equal to translation by1. This proves
metabelian structure directly; it does not need Burnside's finite-order
hypothesis. Modulo m coprime to6 the exact group has order
`m*|<2,3> modm|`, which may have more than two prime factors. The incoming
[finite affine lift theorem](collatz_affine_guarded_lifts_20261004.md)
also retains the common exponent pair across composite moduli.

The useful structural separation is **multiplier exponents / ordered
translation carry / integer guards**. Solvability controls the first
two algebraic layers. It does not direct legal integer paths toward1.
The incoming legal-lift theorem realizes every finite group element
above a supplied certified hub, but does not give a descending route
from an arbitrary prescribed source. Our paid controllers add exactly
that direction and size inequality on their specified domains.

## 9 Remaining target and reproducible scope

The next target is payment through the refill shells that the safe-exit
rules miss. A concrete testbed is the residual part of n=32t-5 with
refill cofactor oddpart(9t-1), keeping the original rank fixed. A second
is the remaining exponents3 mod8 of the ternary tower. The paid
families above handle infinite subdomains; they do not prove that
every residual continuation must eventually enter one of them.

The [script](../../04-computation/experiments/collatz_paid_portrait_controllers_20261004.py)
and its [JSON](collatz_paid_portrait_controllers_20261004.json) and
[stdout](collatz_paid_portrait_controllers_20261004.out) freeze the
comparison bank. The universes include180 words of length at most4,
letters1..4 and cost at most9 on102 signed inputs; five anchors through
five affine words; arbitrary-refill witnesses through precision40;
576 exact-exit families at q1..64; the h<=24,q<=24 rational atlas;
all50000 positive odd sources below100000; and eleven power-three
exponent classes. Large power sources are expanded only when their
exponents are at most10000. Literal replay, affine coefficients,
symbolic all-height guards, and rational density bounds are kept distinct.

Run `python3 -B 04-computation/experiments/collatz_paid_portrait_controllers_20261004.py`
and repeat with `python3 -O -B`; the artifacts must agree byte for byte.
The proofs carry the infinite quantifiers. These finite checks do not
certify universal coverage or identify every negative basin.
