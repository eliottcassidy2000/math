# Divisor gates and bounded-fiber barriers for the twelve forms

Status: **PROVED, independently audited in the stated scope**. This note supplies
all-order elementary proofs and exact replay controls. The twelve-template occupancy assertion remains **OPEN**.
This is not a larger coverage census and proves no new universal occupancy.

## 0. Inheritance and scope

The input is the twelve matrices in `dependencies/MENUS.json`, copied unchanged
from the independently audited compact-central-arithmetic release. The physical
construction and source-bound dependency is the unchanged reviewed
`dependencies/CORE_PROOF.md`, Sections 2--4.

The source coordinates are always

    tau(u)=(g(u), u mod 2, u, 1),  u>=2,
    g(1)=1, g(2h)=g(h), g(2h-1)=g(h)+h-(h mod 2).

The inherited closest mechanism is the complete progression inverse, not
surjectivity of an ambient bilinear form. The hostile boundary remains the
all-height twelve-form gap 2134117. The corrected near miss is the temptation
to select g and parity independently. The less-used sidecar here is elementary
factorization of the actual source-sum coordinates.

Concept board: actual g fibers; parity sheets; rank-one bilinear matrices;
integer divisors; dyadic heights. The anchor is uniform occupancy, the niche
is exact factorization, and the wildcard is a norm-style change of coordinates.
The latter produces a useful classification but no representation theorem.

## 1. Exact odd cores and endpoint equalities

**Claim 1.** For every u>=1, g(u) is odd and 1<=g(u)<=u. Moreover:

1. g(u)=1 exactly when u is a power of two.
2. For every odd a>=3, (a+3)/2 <= g(a) <= a.
3. Equality in the lower bound occurs exactly at a=2^k-1, k>=2.
4. Equality in the upper bound occurs exactly at a=2^k+1, k>=1.

Proof. Parity and positivity follow directly by induction, since the increment
h-(h mod 2) is even. For even u, g(u)=g(u/2)<=u/2. For odd u=2h-1,
if h is odd, g(u)<=h+h-1=u. If h is even, g(h)<=h/2 and
h+h/2<=2h-1 for h>=2. This proves the upper bound.

Stripping all factors of two shows that g(u)=g(a), where a is its odd core.
Every odd a>=3 adds at least two to a positive value, so g(a)>1. This proves
item 1. Write a=2h-1. If h is even then g(a)=h+g(h)>=h+1=(a+3)/2,
with equality exactly when h is a power of two. If h is odd, h>=3 and
g(h)>=3, hence g(a)>=h+2>(a+3)/2. This proves items 2 and 3.

For the upper equality, even h gives equality only at h=2, so a=3.
Odd h gives equality exactly when g(h)=h. Repeatedly descending a->(a+1)/2
therefore reaches 3, and the reverse recursion a->2a-1 gives precisely
3,5,9,17,... . This proves item 4.

In particular, every fixed g fiber consists of finitely many odd cores,
each multiplied by arbitrary nonnegative powers of two. The actual index
restriction u>=2 excludes only the height-zero occurrence of core a=1.

## 2. The J form collapses two dyadic heights to one

Write u=2^r a and v=2^s b with a,b odd and positive. Then

    Psi_J(u,v)=g(a)g(b)+ab-1,                  if r=s=0,
    Psi_J(u,v)=g(a)g(b)+2^(r+s) ab,            if r+s>=1.

This is the literal formula g(u)g(v)+uv-(u mod 2)(v mod 2).
Conversely, the first line is realized exactly when a,b>=3. Every second-line
expression with k=r+s>=1 is realized, except that a=b=1 requires k>=2.
To construct indices, assign one positive height to each core equal to 1,
then distribute the remaining height arbitrarily. This also proves necessity
of the exceptional restriction.

Thus for fixed odd cores the two source heights do not give two independent
additive degrees of freedom in J. They give one dyadic exponential family.
This is an exact set description, not a claim that the family covers a tail.

## 3. Two exact multiplicative gates

Put x(u)=u+g(u), p=u mod 2 and q=v mod 2. Let F_c denote JJZJ when c=3 and
JJZZJ when c=5. Direct expansion of the frozen matrices gives

    F_c(u,v)=c x(u)x(v)-q x(u)-p x(v)+pq,
    c F_c(u,v)=(c x(u)-p)(c x(v)-q)+(c-1)pq.       (1)

These are identities on actual indices of every height.

Define X_e={x(u): u>=2 even} and X_o={x(u): u>=3 odd}. Because g is odd,
X_e consists of odd integers >=3 and X_o of even integers >=6.
The following three cases are an exact if-and-only-if divisor criterion:

- Both sources even: n=cxy with x,y in X_e.
- Exactly one source even: n=x(cy-1) with x in X_e and y in X_o.
- Both sources odd: cn-(c-1)=(cx-1)(cy-1) with x,y in X_o.

Each admissible divisor pair gives actual u,v by finding their source-sum
fibers; either index order is allowed. This is not the weaker assertion that
factors merely have compatible residues.

There is an explicit complete inverse for every source-sum fiber:

- For odd x, test precisely the even integers
  ceil(2x/3)<=u<=x-1 against u+g(u)=x.
- For even x, test precisely the odd integers
  ceil(x/2)<=u<=floor((2x-3)/3), u>=3, against u+g(u)=x.

For even u the inequalities follow from 1<=g(u)<=u/2. For odd u they follow
from (u+3)/2<=g(u)<=u. Testing the displayed equality is necessary; interval
membership alone is not sufficient. `inverse_factor` in `verify.py` implements
all three cases, with these finite exact fibers, and returns actual indices.
There is no hidden upper-height assumption.

### Prime-target obstruction and hostile control

If n is prime, neither even-source case is possible: both factors in the
mixed case are >=3 and >=6c-1, and the both-even case has the nontrivial factor
c. Consequently every prime F_c value must have

    cn-(c-1)=(cx-1)(cy-1),  x,y in X_o,

with both factors at least 6c-1 and congruent to -1 modulo 2c. In particular,
if both n and cn-(c-1) are prime, F_c has no representation of n at any height.

A deliberately scoped hostile example above the occupancy threshold is

    n=2135743,    3n-2=6407227,    5n-4=10678711.

All three numbers are prime, as exact trial division verifies. Thus neither
JJZJ nor JJZZJ represents this n, at any height. It is NOT a gap of the full
menu, because the actual-source witness

    Psi_J(7,191317)=2135743

is checked directly. The example prevents confusing a two-template divisor
obstruction with a counterexample to the twelve-template hypothesis.

## 4. Exact scope of a polynomial factorization strategy

For a nonempty template, its even/even core has shape

    a AB+c Av+d uB+b uv,

where A=g(u), B=g(v). The determinants ab-cd for the eleven nonempty forms are

    J:1, JJ:0, JJJ:-2, JOZJ:3, JJZJ:0, JJJJ:3,
    JOOZJ:5, JJZZJ:0, JOZJJ:-7, JJZJJ:0, JJJJJ:-1.

Hence only JJ, JJZJ, JJZZJ and JJZJJ have a separated product on the even/even
sheet. For a factorization on all formal parity sheets, allowing a rational
multiple of pq to be subtracted, the ONLY nonempty templates are JJZJ and
JJZZJ.

Proof. Subtracting a multiple of pq changes only the parity/parity entry
of the 3-by-3 homogeneous matrix on (A,p,u) and (B,q,v). A separated product
of linear forms is exactly a rank-one matrix. Its (A,u) by (B,v) minor must
therefore vanish. This leaves the four zero-determinant cores above.
For JJ, the parity row on (B,v) is (-1,0), not proportional to its core row
(1,1). For JJZJJ, the q column on (A,u) is (-1,-2), not proportional to
(5,6). No choice of the pq entry repairs these nonzero minors. The other two
factor exactly by (1). This is a formal polynomial classification, not a
no-go theorem for identities special to the actual g sequence or for coupled
or nonlinear changes of indices.

A norm-style coordinate change on the even/even sheet always gives

    a n=(aA+d u)(aB+c v)+(ab-cd)uv.             (2)

This exposes the listed determinant as the residual coupling. Formula (2)
is a bilinear identity, not generally a norm of one algebraic integer, and
it loses the source-digit restrictions if its new coordinates are allowed
to vary freely. No sum-of-squares or norm-surjectivity theorem may be applied
without restoring those restrictions.

## 5. A bounded-g-fiber proof cannot cover all large targets

Let G>=1 be an integer, N>=1, and let L=floor(log2(2N)). The number of distinct
values <=N represented by any of the twelve templates with

    g(u)<=G and g(v)<=G

is at most

    12 G^2 (L+1)(L+2)/2.                       (3)

Proof. Claim 1 gives at most G possible positive odd cores for each index:
core 1 is allowed, and every core a>=3 obeys a<=2G-3. The released physical
product bound is Psi_W(u,v)>=uv/2 for every template. Thus a value <=N
requires uv<=2N. With u=2^r a, v=2^s b, this forces r+s<=L. There are
(L+1)(L+2)/2 nonnegative integer pairs (r,s) satisfying that inequality.
Multiply this upper bound by G^2 possible core pairs and twelve templates.
Excluded singleton indices and repeated values only decrease the count.

For J alone, the height collapse in Section 2 improves the upper bound to
G^2(L+1): one potential odd/odd value and at most L positive-height values
per core pair.

Since the right side of (3) is O_G((log N)^2), no fixed bound on both
source g-coordinates can represent every sufficiently large odd target.
Quantitatively, if every odd n with N/2<n<=N has such a representation, then

    12 G^2 (L+1)(L+2)/2 >= floor(N/4).

This forces the permitted g bound to grow at least on the scale sqrt(N)/log N
in any such interval-covering strategy. It does not require every target,
or every particular witness, to have coordinates that large. It also says
nothing against a construction which keeps only one g coordinate bounded.

## 6. Reproduction and limits

Run from this packet directory:

    python3 verify.py
    python3 -O verify.py

The verifier proves the factorization and its formal classification by exact
coefficient/minor calculations. Finite controls independently compare g,
endpoint equalities, source-sum fibers, and the complete divisor inverse;
trial division checks the hostile prime triple, and a positive J witness
keeps the scope honest. The source-sum controls run through 1024, the g
controls through 4096, and 211 small target checks compare the divisor inverse
with a direct complete inverse. These checks are not an extension of the
large occupancy census. Bounds and all-order claims rest on the displayed
proofs, not extrapolation from the checks.

No literature novelty is asserted. No target above 2134117 has been proved
missing from the full menu. Uniform twelve-template occupancy remains OPEN.
