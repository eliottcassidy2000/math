# Exact parity-class responses and physical two-source fiber repairs

Session 200. **PROVED, subject to independent audit:** for every reduced
positive rational a/b with odd b and 1/2<a/b<1, the response aU-bD on strong
unit-interface sources is exactly the parity class a-b modulo 2. This extends
the previous even-numerator theorem. It supplies infinitely many actual
two-source C2 reassociations on explicit rational fibers. **OPEN:** prescribed
positive Hamilton-count coverage, including the full odd spectrum outside
7 and 21.

## 1. Inheritance and physical domain

The closest mechanism is Session 198's even-numerator response theorem,
frozen in `dependencies/EVEN_NUMERATOR_RESPONSE.md`. Its clock starts O,
erasing the incoming endpoint coordinate. With odd numerator the clock
starts Z, and ignoring the incoming endpoint term gives a wrong decrement.
The canonical hostile outside our hypotheses is slope 3/4, whose positive
responses are only 1 modulo 4. The least-used sidecar is E=d at both ends
of a complete periodic clock. Concept board: source words, endpoint memory,
signed response, parity, two-source fibers, positive-count height.

A nontrivial source is a rooted tournament with exactly one root-start
Hamilton path. Its profile (D,E,U,d) records total Hamilton paths, root-ending
Hamilton paths, ordered spanning two-path covers with the nonempty reserved
path starting at the root, and root-deletion Hamilton paths. H(empty)=1.
The physical chronological operations are

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

O adds a sink. Z adds a vertex beating only the old root and retains that
root. J makes the same birth the new root. The transitive edge has profile
t=(1,0,2,1). O/Z preserve transitive root deletion and d=1. A final Z gives
strongness and E=d=1. On strong sources D<=U<=2(D-E); with E=d=1, D,U are odd.

The inherited physical modular compiler says: for every even M>=2 and
odd D0,U0 modulo M, there are strong sources with transitive root deletion,
E=d=1, D=D0 and U=U0 modulo M, and arbitrarily large D. Its literal V(L) has

    (D,E,U,d)=(g(L)+L-(L mod 2),1,2L-1,1), L>=2,

where g(1)=1, g(2h)=g(h), g(2h-1)=g(h)+h-(h mod 2).
`STRONG_SOURCE_INTERFACE.md` and `JOINT_MODULAR_GADGET.md` give the graph
and integer binary-word compiler. This is residue freedom, not independent
exact D,U freedom.

## 2. A cyclic clock for either numerator parity

Fix coprime positive a,b, with b odd and b/2<a<b. Choose k>=1 with
2^k=1 modulo b. Put B=2^k and A=a(B-1)/b. Then A is integral and
B/2<=A<B-1. Let epsilon_i be bit i of A, i=0,...,k-1. The word W uses Z
for one and O for zero, in that order. Its last letter is Z and it is mixed.
Its first letter is Z exactly when a is odd, since (B-1)/b is odd. Put

    R=sum_i epsilon_i 2^(k-1-i),
    C=sum_(j<i) epsilon_j epsilon_i 2^(i-j-1)
       +sum_(i=0,...,k-1) epsilon_(i-1)epsilon_i,
    K=aR-bC,                 epsilon_(-1)=epsilon_(k-1)=1.

Direct iteration on every input with E=d gives

    D'=D+AU-Cd, E'=d, U'=BU-Rd, d'=d,
    (aU-bD)(s.W)=(aU-bD)(s)-Kd.                 (1)

The cyclic adjacent-one sum includes epsilon_0: precisely the incoming E=d
subtraction at the first Z. It vanishes when a is even, recovering the older
formula. Equation (1) is not asserted for arbitrary E when W starts Z.

**Lemma.** K is a positive even integer.

For positivity use the formal periodic cover orbit u_0=R/(B-1),
u_(i+1)=2u_i-epsilon_i. Its binary expansion repeats the epsilon word,
so 0<u_i<1 and the orbit closes. Its formal D increment is

    Delta=sum_i epsilon_i(u_i-epsilon_(i-1))=AR/(B-1)-C=K/b.

This evaluates a word coefficient, not a physical source. For a cyclic run
of j>=1 ones followed by zero, let v be u at that following zero. The sum
of the j deficits 1-u is (1-v)(1-2^(-j)), while the adjacent-one subtraction
is j-1. The run contributes v+(1-v)2^(-j)>0. There is at least one run,
so K>0. Integrality is clear. Apply W to the actual triangle (3,1,3,1):
initial and final responses have the same parity a-b, so K is even.

Examples are 3/5: W=ZOOZ,K=2; 5/7: W=ZOZ,K=4; 5/9: W=ZZOOOZ,K=2.

## 3. Exact full-parity-class coverage

**Theorem 1.** For coprime positive a,b with b odd, the image aU-bD on
strong sources with transitive root deletion and E=d=1 is the entire
parity class (a-b)+2Z if and only if 1/2<a/b<1. In that interval every
permitted target has infinitely many distinct physical source profiles,
each with an explicit finite recipe.

Fix r=a-b modulo 2 and choose ell>=1 with

    c=2^ell(2a-b)>b.                              (2)

Appending O^ell Z to a unit-interface strong input gives

    (D+2^ell U,1,2^(ell+1)U-1,1),

whose response is S=cU-bD-a. To select odd seed residues D0,U0 modulo K
with S=r modulo K, write D0=1+2x and U0=1+2y. The equation becomes

    cy-bx=(r-c+b+a)/2 modulo K/2.

Its right side is integral. Since b is odd and gcd(a,b)=1, gcd(c,b)=1;
Bezout supplies x,y even as an exact integer equation. The modular compiler
realizes these residues with D arbitrarily large and
(c-b)D-a>=r+K. Because U>=D, the actual S>=r+K. Thus h=(S-r)/K is a positive
integer. Append W exactly h times. Equation (1) applies at the true E=d=1
interface each time and gives the final response exactly r. The physical
word retains root-start count one, unit deletion, and transitive root
deletion, and ends strong.

A compressed recipe consists of the modular seed index L, prefix O^ell Z,
word W, and integer repetition h. As seed D tends to infinity, final D is
at least its seed value, so infinitely many distinct profiles occur.
This does not independently prescribe their positive Hamilton counts.

Parity gives necessity of r=a-b modulo 2. For 0<a<=b/2 the strong cone gives
aU-bD<0; for a>=b, U>=D gives aU-bD>=0. Neither image contains a whole
two-sided parity class, proving the iff within the odd-denominator setting.

The complete O/Z orbit of t is

    (D,E,U,d)=(g(n),n mod 2,n,1), n>=2;

odd n>=3 are strong. The modular seeds belong to this orbit. Hence the
equivalent exact arithmetic statement is

    {an-bg(n): n>=3 odd}=(a-b)+2Z.

Conversely every physical d=1 profile belongs to this O/Z orbit. In a
complete O/Z/J history, d cannot decrease, and J replaces d by D. Thus a
J in a history ending at d=1 must occur at D=1, where the source is
transitive with profile (1,0,2^j,1), j>=1. Its J output is
(2^j+1,1,2^j+1,1), also obtained as t.Z^j. Once D>1, no later J can keep
d=1. Replacing that possible J segment by its O/Z recipe leaves the same
profile and proves the converse. This uses the complete physical grammar,
not an assumption that every formal unit-interface tuple is realizable.

The odd-denominator assumption cannot simply be dropped. For 3/4, write
n=2m-1. The recurrence gives 3n-4g(n)=2m-4g(m)+4(m mod 2)-3. For even m
this is 1 modulo 4. For odd m>=3 the strong cone gives g(m)>=m/2+1, so
the response is at most -3. The excluded m=1 gives -1. Every positive
response is therefore 1 modulo 4, and target 3 is absent.

## 4. Target-preserving physical C2 repairs on rational fibers

For two actual nontrivial sources the complete physical C2 form is

    N(v,w)=(D_v-E_v+U_v)D_w-(D_v-E_v+d_v)E_w
           +(D_v+3U_v-d_v)U_w-(U_v+E_v)d_w.

It constructs a strong tournament with exactly two Hamilton cycles. Put
p(v)=D_v+d_v and rho(v)=D_v-E_v+U_v-d_v. The inherited exact defect is

    N(v,Jw)-N(Jv,w)=rho(v)p(w)-p(v)rho(w).          (3)

**Theorem 2.** Let v be any actual nontrivial source with
2<rho(v)/p(v)<3. Reduce p(v)/(rho(v)-p(v)) to a/b. If b is odd, infinitely
many strong unit-interface physical sources w satisfy

    N(v,Jw)=N(Jv,w).                              (4)

Except for at most one w, the unordered pairs of response profiles
{v,Jw} and {Jv,w} differ. These are genuine target-preserving two-source
presentations, not merely scalar symmetry.

Indeed the interval is exactly 1/2<a/b<1. For w=(D,1,U,1), equality of
rho/p in (3) is equivalent to aU-bD=3a+b. This target has parity a-b,
so Theorem 1 provides infinitely many such physical w. Both grafts have
exactly two Hamilton cycles. If the unordered pairs were equal, strict U
growth under J rules out v=Jv, leaving w=v. Discard that at most one
profile. The bound N(x,y)>=3U_x U_y+1 shows their common positive targets
are unbounded too.

A small literal control is

    v=t.ZOOZ=(15,1,23,1),  w=t.ZOZOZ=(31,1,43,1).

Both have rho/p=9/4 and common target 7649. The verifier directly counts
both 15-vertex grafts and their two Hamilton cycles. The old hostile
v=t,w=t.O still gives 41 versus 45:
no universal transfer is claimed.

An example specifically requiring the new odd-numerator theorem is

    v=t.ZJOZO=(17,0,42,3), rho/p=14/5.

Its unit companions must satisfy 5U-9D=24. The 5/9 clock has W=ZZOOOZ,K=2,
and the construction gives exact compressed companions. No such companion
occurs in the complete U<=5000 unit-source bank; this bounded control is
not used to prove existence or nonexistence at larger heights.

### Exact simultaneous-positive-target condition

The constraint aU-bD=3a+b, with gcd(a,b)=1, implies
a(U-3)=b(D+1). Odd D,U and odd b give exactly

    (D,U)=(2a h-1,2b h+3), h>=1.                 (5)

Thus a unit companion exists at a specified h if and only if

    g(2b h+3)=2a h-1.                            (6)

This is a necessary-and-sufficient arithmetic condition, because the full
physical O/Z orbit supplies the source when (6) holds. No integrality or
physical membership condition has been silently removed.

For the fixed left source v=(V,E,H,d), put

    A_v=V-E+3H-2d, B_v=2V-E+4H-d.

Expansion of N(v,Jw) gives A_v D+B_v U+d-H. Therefore its common repaired
positive target is exactly

    n=2(A_v a+B_v b)h-A_v+3B_v+d-H.              (7)

Both A_v and B_v are positive on physical nontrivial sources. For a given
n, equation (7) fixes h uniquely. Theorem 1 proves that infinitely many h
satisfy (6), but it does not prescribe that h. This is the exact remaining
arithmetic condition for this repair family, rather than a claim that two
free coordinates were constructed.

For v=(17,0,42,3), equations (6)-(7) become

    g(18h+3)=10h-1,       n=4952h+421.

An explicit substantially smaller witness than the general modular compiler
uses seed t.ZOZZOZZOOOZ, whose (D,U) are (1463,2639), then appends
(ZZOOOZ)^2. Its final profile is (6003419,1,10806159,1), h=600342,
and common positive target is 2972894005. Both grafts have 34 vertices.
The graph word, binary g recurrence, and two-source formulas verify this
large graph certificate; no literal adjacency-count computation at its
full order is claimed.

## 5. Arithmetic sidecars and the remaining existence obligation

K is invariant under cyclic rotation and bit-complement of its binary
clock. Rotation reindexes the periodic sum Delta. For complement, u_i
becomes 1-u_i and epsilon_i becomes 1-epsilon_i; the difference in Delta
is sum epsilon_i-sum u_i=0, by summing the periodic cover recurrence.
Thus numerator residues related by multiplication by 2 or negation modulo
b have equal K when using the same k. In the formal calculation for a
clock ending zero, the incoming E must be epsilon_(k-1)d, not d.

Denominator or period alone does not determine K: slopes 9/17 and 10/17
both have minimal period 8 but K=2 and K=30. A common Euler period or
totient collision is not itself an exact-response-preserving operation.
The decrement and endpoint interface are indispensable sidecars. No
Carmichael conjecture is assumed or proved.

The Diophantine connection is constructive: Bezout supplies integral
residues, the physical compiler supplies a genuine integral source, and
(S-r)/K is proved a positive integer before it becomes a repetition count.
Diophantine definability alone would not establish any of these nonempty
fibers. No result about integer rings in number fields is imported.

Theorem 2 repairs a forbidden universal reassociation by restricting it
to an explicitly populated physical fiber. It does not start with an
independently prescribed positive target n and preserve that n while
searching for repair. It constructs an unbounded family of n, each with
two physical presentations. D,U may become enormous while their signed
response remains fixed. Universal positive-count coverage stays open.

## 6. Verification

Run `python3 verify.py` and `python3 -O verify.py` serially; receipts must
agree. The standalone checker uses explicit exceptions. It checks clocks
for all reduced interior slopes with odd b<=31, the E=d transfer on a
three-vector basis, positive run decompositions, cyclic/complement laws,
compressed recipes of both parities, independently materialized profiles
when manageable, literal small sources, the new 5/9 companion, and the
two C2 grafts in the displayed collision. It hashes frozen dependencies.
Finite checks corroborate the all-order proof, not saturation extrapolation.
