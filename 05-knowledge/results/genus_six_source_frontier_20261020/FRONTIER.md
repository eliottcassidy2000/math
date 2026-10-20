# The genus-six frontier after the conductor-twelve bound

20 October 2026. **PROVED classification and realization statements; independently audited PASS in the stated scope.**

This packet proves elementary classification and explicit realization statements.
It does **not** exclude generic codimension six and does **not** improve the
audited conductor bound twelve. The planar Jacobian conjecture remains open.

## 1. Outcome and exact scope

The complete symmetric genus-six list is

    <2,13>, <3,7>, <4,5>, <4,6,9>, <5,7,8,9>, <6,7,8,9,10>.

Here genus means the number of numerical-semigroup gaps, not the genus of
the normalization. In the primitive polynomial-algebra setting the
normalization is the affine line, while the algebra can have positive
codimension.

The complete list of nonprincipal plane characteristic sequences, up to
exchanging their first two entries, is

    (2,13), (3,7), (4,5), (4,10,5), (4,6,9), (6,9,4).

All six sequences have **explicit immersed primitive polynomial
realizations**, including both alternate presentations. The last two
symmetric semigroups in the first list are not planar in the
one-place-at-infinity sense. They are not candidates for the primitive
two-generated polynomial source algebra.

Here “immersed” means that the derivative pair never vanishes, equivalently
gcd(p',q')=1; it does not assert injectivity or a locally closed scheme
embedding. “Alternate” refers to presentation types of a given semigroup,
not an assertion that the particular displayed algebras are isomorphic.

Under the inherited primitive monic Keller-source hypotheses, the separate
audited quadratic-semigroup theorem already excludes <2,13>. Thus the
remaining terminal coordinate-degree pairs at generic codimension six are,
up to interchange,

    (3,7), (4,5), (4,10), (4,6), (6,9).

This is a necessary frontier map, not a claim that any of these pairs can
occur in a noninvertible Keller map. In particular, an immersed
one-variable parametrization is not a two-variable constant-Jacobian map.

### Inheritance and live concepts

The closest proved mechanism is the frozen genus-five quartic/sextic
obstruction, together with its exact conductor-twelve synthesis. Its
full degree-separated module basis is retained here.

The canonical hostile is an alternate embedding with no minimal-degree
semigroup generator among its two coordinate degrees. The new concrete
(4,10,5) and (6,9,4) realizations make this danger explicit at genus six.
The corrected near miss is to infer polynomial realizability from the
characteristic axioms alone. The least-used sidecar is the derivative
gcd: it separates realizability from immersion without asserting either
is sufficient for a Keller source.

Anchor: exhaust the genus-six source types. Niche: construct the alternate
presentations. Wildcard: follow their parameter degenerations.
The concept board is symmetric gaps, gcd-drop sequences, terminal
presentations, full integral module bases, immersion, and the first
nonvanishing odd cancellation coefficient.

## 2. Exhaustive symmetric semigroup classification

A symmetric numerical semigroup with six gaps has Frobenius number eleven.
Exactly one of n and 11-n belongs for 0<=n<=11.

If 2 belongs, every nonnegative even integer belongs. Symmetry excludes all
odd integers below twelve, and all integers at least twelve belong. This
gives <2,13>.

Suppose 2 is absent and 3 belongs. Then 6 and 9 belong. Symmetry makes 5
and 8 absent. The element 4 must be absent: otherwise 3+4=7 would belong
while symmetry with 4 excludes 7. Thus 7 belongs and the semigroup is
<3,7>.

Suppose both 2 and 3 are absent. Then 9 and 8 belong. The remaining
independent choices are 4 versus 7 and 5 versus 6. Their four combinations
are exactly

    <4,5>, <4,6,9>, <5,7,8,9>, <6,7,8,9,10>.

Each choice is closed under all sums below twelve, and all larger sums
are already past the conductor. The full gap and Apéry certificates are:

| Minimal generators | Gaps | Apéry set, in residue order modulo the multiplicity |
|---|---|---|
| 2,13 | 1,3,5,7,9,11 | 0,13 |
| 3,7 | 1,2,4,5,8,11 | 0,7,14 |
| 4,5 | 1,2,3,6,7,11 | 0,5,10,15 |
| 4,6,9 | 1,2,3,5,7,11 | 0,9,6,15 |
| 5,7,8,9 | 1,2,3,4,6,11 | 0,16,7,8,9 |
| 6,7,8,9,10 | 1,2,3,4,5,11 | 0,7,8,9,10,17 |

The checker separately enumerates all 4,096 subsets of {0,...,11} and
all 64 complementary-pair selections. Both recover exactly this list.
This finite computation corroborates the complete case split above.

## 3. Complete plane characteristic universe

The primary structural input is the characteristic-sequence description
for curves with one place at infinity. With
d_i=gcd(delta_0,...,delta_(i-1)) and n_i=d_i/d_(i+1), the conditions are:
final gcd one; n_i>1 for i>=2; n_i delta_i in the preceding monoid; and
delta_i<n_(i-1)delta_(i-1) for i>=2. For conductor twelve and a
nonprincipal initial pair, delta_0,delta_1<=13 and later entries<=11.
The conductor is 1-delta_0+sum_(i>=1)(n_i-1)delta_i, and every primitive
semigroup generator must occur. These are cited inputs, source-pinned in
[PRIMARY_SOURCES.md](PRIMARY_SOURCES.md).

The actual coordinate pole orders start the sequence. Thus the first
two entries, not merely the primitive semigroup generators, control
the terminal coordinate degrees. This identification is separately
source-checked in the same source note.

### Why this is a genuinely finite exhaustive test

For nonprincipal positive a,b<=13, put d=gcd(a,b). Write a=du, b=dv,
where u,v are coprime and neither is one. The smaller possible unequal
pair is 2,3, so 3d<=13 and d<=4. In fact the possible d are 1,2,3,4.

After the first pair, each new entry strictly lowers this gcd. For d=2
or 3 at most one further entry can occur. For d=4 there can be two,
through 4 -> 2 -> 1. Therefore every sequence has two, three, or four
entries. A length-three-only scan would not yet be justified. As a
positive control, (8,12,10,5) is a valid length-four sequence, but its
conductor is twenty, not twelve.

The entire Cartesian-product universe has:

    length two:      108 nonprincipal initial pairs;
    length three:  1,188 candidate sequences;
    length four:  13,068 candidate sequences.

Testing the stated axioms and conductor formula yields exactly twelve
oriented sequences:

    (2,13), (13,2);
    (3,7), (7,3);
    (4,5), (5,4);
    (4,10,5), (10,4,5);
    (4,6,9), (6,4,9);
    (6,9,4), (9,6,4).

A second implementation starts from each of the six target semigroups,
recurses only through strictly dropping gcds, and selects by the exact
generated monoid instead of the conductor formula. It gives the same
list. In particular, no length-four sequence contributes at conductor
twelve.

### Elementary nonplanarity certificates

The semigroup <6,7,8,9,10> has five primitive generators. Every possible
characteristic sequence here has at most four entries, so it is not planar.

The semigroup <5,7,8,9> has four primitive generators. A sequence with
four entries must start with gcd four. In its allowed initial universe
{5,7,8,9,10,12,13}, the only nonprincipal gcd-four pair is (8,12), up to
interchange. But 12 is not primitive. That sequence has at most three
slots left for the four required primitives. Every other initial pair
allows at most three total entries. Thus this semigroup is also not planar.

Here “planar” refers to the characteristic degree semigroup at infinity,
not to the different convention for local plane-branch semigroups.

## 4. Why each explicit realization certifies the whole algebra

Let p be a monic polynomial of degree n over a characteristic-zero field k.
Suppose k(p,q)=k(t). Then q is integral over k[p], and its minimal polynomial
over k(p) has degree n. Its coefficients belong to k[p], which is integrally
closed. Consequently 1,q,...,q^(n-1) is a free k[p]-basis of k[p,q].

More generally, when r=q^2-U(p)q-V(p) with U,V in k[p], the ordered list

    1,q,r,qr,r^2,qr^2,...

up to n entries is obtained by a triangular matrix with diagonal entries
one: its successive entries are monic of degrees 0,1,...,n-1 in q.
This includes U=0, used in the explicit realizations, and the nonzero-U
remainder used in Section 6. The particular examples below therefore
give bases of the **whole** algebra, not merely a conveniently generated
subalgebra.

If these basis elements have degrees e_0,...,e_(n-1) in distinct residue
classes modulo n, nonzero summands f_i(p)b_i have different degrees.
No leading cancellation is possible. The degree set is exactly the union
of e_i+nN, and the codimension is

    sum_i (e_i - (e_i mod n))/n.

Primitivity can be certified without solving for t: the intermediate-field
index divides the degree of every nonconstant polynomial in the algebra.
If several exhibited degrees have gcd one, the index is one.

## 5. Six immersed polynomial realizations

All formulas are over any characteristic-zero field. Interchanging p and q
realizes the other orientation.

### The three coprime initial pairs

For (n,m)=(2,13),(3,7),(4,5), take

    p=t^n, q=t^m+t.

The degree gcd is one, the full basis degrees are 0,m,...,(n-1)m,
and hence the semigroup is <n,m>. The only zero of p' is t=0, where
q'=1, so these are immersed. Each has codimension six.

### The standard quartic/sextic pattern (4,6,9)

Take

    p=t^4,
    q=t^6+t^3+t,
    r=q^2-p^3=2t^9+2t^7+t^6+2t^4+t^2.

The degrees 4,6,9 force the fraction field to be k(t). The full basis
1,q,r,qr has degrees 0,6,9,15, with distinct residues modulo four.
It gives exactly <4,6,9> and codimension six. Again p' vanishes only at
zero, where q'=1.

### The alternate quartic/decic pattern (4,10,5)

Take

    p=t^4+t,
    q=t^10+(5/2)t^7+(15/8)t^4+(5/16)t.

The exact identity is

    r=q^2-p^5-(5/64)p^2=t^5/64+5t^2/256.

Degrees 4 and 5 give primitivity. The full basis 1,q,r,qr has degrees
0,10,5,15, so its exact semigroup is <4,5> and its codimension is six.
The Euclidean derivative remainder is q' mod p'=-5/8, proving immersion.

This is a positive hostile against replacing “semigroup <4,5>” by “actual
terminal coordinate degrees (4,5).” The degree-ten coordinate cannot
simply be discarded from the terminal map.

### The alternate sextic/nonic pattern (6,9,4)

Take

    p=t^6+t^4,
    q=t^9+(3/2)t^7+(3/8)t^5-(5/32)t^3-(3/128)t.

The exact identity is

    r=q^2-p^3+(3/16)p^2-(7/1024)p=t^4/2048+9t^2/16384.

The degree gcd of 6,9,4 is one. The full basis 1,q,r,qr,r^2,qr^2 has
degrees 0,9,4,13,8,17, distinct modulo six. Therefore its exact semigroup
is <4,6,9> and its codimension is six.

For immersion,

    p'=2t^3(3t^2+2),
    q'(0)=-3/128,
    q' mod (3t^2+2)=-27/128.

Thus no zero of p' is a zero of q'. This is a positive hostile against
discarding terminal (6,9) when the smallest semigroup generator is four.

The supplementary
[realization note](realizations/ALTERNATE_DEGREE_REALIZATIONS.md) gives
parameter families, their exact nonprimitive and non-immersed boundaries,
and a degeneration from (6,9,4) to (6,9,2). Those identities are
separately checked; the two concrete certificates above suffice for the
frontier classification.

None of these one-variable examples is asserted to satisfy the planar
Keller condition. With coefficients independent of x, their planar
Jacobian is zero.

## 6. The precise failure of the genus-five quartic/sextic mechanism

Center a monic quartic and, for comparison of the k-algebra only, remove
the scalar terms:

    p0=t^4+a t^2+b t,
    q0=t^6+A t^5+c t^4+B t^3+e t^2+h t.

The degree-eleven coefficient of q0^2-p0^3 is 2A. When it is nonzero,
the full basis degrees 0,6,11,17 give codimension seven. Therefore
codimension at most six still forces A=0.

After A=0, subtract (2c-3a)p0q0 from q0^2-p0^3. The degree-nine
coefficient of the remainder is

    beta=2B-3b.

If beta is nonzero, the full basis degrees are 0,6,9,15, and the
codimension is **exactly six**. Consequently the earlier delta<=5 proof
could set B=3b/2, but a delta-six proof cannot do so.

The preceding immersed example p=t^4, q=t^6+t^3+t has b=0, B=1 and
beta=2. It is the smallest cheap certificate that the old
normal-form implication fails at the next boundary. It is a hostile
to that implication, not to the Keller conjecture.

There is a useful differential sidecar. Restore the actual scalar terms:

    P=t^4+a(x)t^2+b(x)t+d(x),
    Q=t^6+c(x)t^4+B(x)t^3+e(x)t^2+h(x)t+f(x).

For J=P_x Q_t-P_t Q_x the top coefficients are

    [t^7]J=6a'-4c',
    [t^6]J=6b'-4B'=-2 beta'.

Thus a hypothetical Keller pair forces beta to be a constant, not zero.
The unresolved quartic/sextic delta-six branch has a genuine nonzero
constant defect. Reusing the previous first-integral proof with beta
silently set to zero would drop precisely this branch. No scalar
subtraction is being claimed to preserve the actual Jacobian.

## 7. Translation to the inherited Keller-source setting

The inherited setting is R=C[x], K=C(x), P,Q in R[y], one coordinate
monic of positive y-degree, K(P,Q)=K(y), and J(P,Q) in C*. The audited
predecessor proves that a nonautomorphism has generic codimension at least
six and generic normalization-conductor degree at least twelve.

If equality occurs, the generic degree semigroup is symmetric of genus six.
The usual constant-leading target reduction preserves K[P,Q] and the
nonzero constant Jacobian. It terminates at nonprincipal positive degrees:
a zero-degree coordinate r(x) would make r'Q_y a unit, forcing Q linear
in y and codimension zero.

The characteristic input applies after scalar extension, which preserves
filtered dimensions and degrees. The normalization remains K[y], with
one rational point at infinity. Hence Section 3 lists every terminal
degree pattern. This uses neither unique embedding nor automatic
globalization of a low-degree generic semigroup generator.

The independently audited all-quadratic-semigroup theorem excludes
<2,m> for every odd m>=3 under these actual monic Keller hypotheses.
It removes <2,13>, not merely the displayed quadratic realization.
All other five terminal pairs remain obligations of this map. In
particular, no further exclusion is inferred from the one-variable
examples or from their derivative gcds.

## 8. Connection contract and stopping boundary

The map is

    actual primitive monic Keller source
      -> generic degree semigroup
      -> bounded characteristic sequences
      -> possible terminal coordinate degrees.

It preserves generic codimension and conductor. Forgetting the
characteristic sequence loses presentation degrees; the two alternate
examples are exact witnesses. Forgetting coefficient dependence on x
loses the Keller equations; immersion cannot restore it. The necessary
sidecars are the actual degree presentation, polynomial coefficient
dependence, and the full Jacobian including scalar terms.

The cheapest next quartic/sextic obligation is the beta in C* branch,
not another proof of the beta=0 case. The other exact target lanes are
terminal (3,7), (4,5), (4,10), and (6,9), with their different full module
bases retained. This packet stops at the auditable frontier. It asserts
no conductor-fourteen result, no classification of all embeddings up
to equivalence, and no statement about arbitrary nonmonic sources or
special-fiber conductor constancy.

## 9. Reproduction

From this directory:

    python3 verify_frontier.py > results.json
    python3 -O verify_frontier.py > results-optimized.json
    cmp results.json results-optimized.json

Both runs pass 129 explicit checks with byte-identical output. The
checker imports no predecessor or realization-worker code. It includes
two semigroup censuses, two characteristic-sequence enumerations, a
length-four positive control, every concrete realization, full
unimodular basis transitions, degree gcds, immersion, exact remainder
identities, and the genus-five-to-six coefficient boundary.
Its resultant-degree checks are supplemental consistency controls, not
independent proofs of primitivity; the degree-gcd field-tower argument
is the load-bearing certificate.

The supplementary realization checker has 48 additional checks and
normal/optimized parity. The independent audit adds 119 checks using
a separate Kunz-coordinate census and exact Apéry computation. Exact
dependency hashes and primary-source limitations are recorded beside
this report; see [the independent audit](audit/INDEPENDENT_AUDIT.md).
