# Exact odd coverage for rational signed source responses

Session198. **PROVED:** for each reduced rational slope a/b with a positive
even numerator, b a positive odd denominator, and 1/2<a/b<1, every odd integer
is exactly aU-bD on an actual strong source with s=E=d=1 and transitive root
deletion. The slope interval is sharp. This is exact coverage of a signed
response, not exact coverage of positive Hamilton-path counts or of the C2/C3
bilinear image.

## 1. Physical domain and the theorem

A source is a nontrivial rooted tournament with exactly one Hamilton path
starting at its root. Its actual response is (D,E,U,d), using the reviewed
Hamilton, endpoint, ordered-cover and deletion conventions. We use the
physical retained-root operations

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d).                 (1)

They preserve root-start count one and d. Appending either vertex to a source
with transitive root deletion preserves that property. A final Z makes the
source strong and sets E=d. All sources below have d=1, so a final Z gives
exactly s=E=d=1, not just congruences of these interfaces.

The strong-source cone is

    D<=U<=2(D-E), E>=1.                         (2)

It follows directly from the positive grammar: a strong last J or Z satisfies
f<=a in its positive coordinates. Hamilton-path counts D are odd.

**Theorem.** Let a,b be positive coprime integers, with a even and b odd.
The image of aU-bD on strong sources with transitive root deletion and
s=E=d=1 is all odd integers if and only if b/2<a<b.

The easy parity statement is exact: U and D are integers and D is odd, so
aU-bD is always odd. The proof of surjectivity constructs an actual graph
recipe. It never treats independent integer D,U as automatically physical.

## 2. A periodic physical word with a positive arithmetic decrement

Assume b/2<a<b. Choose k>=1 with 2^k=1 modulo b and put

    B=2^k, A=a(B-1)/b.

This is an even integer satisfying B/2<=A<B-1: A>(B-1)/2 and A is an integer.
Let epsilon_i
be bit i of A, i=0,...,k-1. Apply O for epsilon_i=0 and Z for epsilon_i=1,
in that chronological order, and call the resulting word W. Its first letter
is O and its last is Z. Let

    R=sum_i epsilon_i 2^(k-1-i),
    C=sum_(j<i) epsilon_j epsilon_i 2^(i-j-1)
       +sum_(i=1,...,k-1) epsilon_(i-1)epsilon_i.

Directly iterating (1) gives, for any nontrivial actual input source,

    D'=D+A U-Cd, U'=BU-Rd, E'=d, d'=d.          (3)

Indeed before letter i the cover count is 2^i U minus the reversed prefix
weighted by d. A Z adds that cover count minus the preceding ending count.
The first letter is O, so the input E disappears and contributes no boundary
term. The two sums in C are exactly those subtractions.

Since bA=a(B-1), the signed response satisfies

    (aU-bD)(s . W)=(aU-bD)(s)-Kd,
    K=aR-bC.                                   (4)

**Lemma.** K is a strictly positive even integer.

For positivity, consider the formal periodic cover orbit with d=1 and
u_0=R/(B-1). Its binary expansion repeats epsilon_0,...,epsilon_(k-1), so

    u_(i+1)=2u_i-epsilon_i, 0<u_i<1,
    epsilon_i=floor(2u_i), u_k=u_0.

The word is mixed, so all inequalities are strict. The formal ending count
before letter i is epsilon_(i-1), with cyclic indexing. Thus the formal D
increment over a period is

    Delta=sum_i epsilon_i(u_i-epsilon_(i-1))
         =A R/(B-1)-C=K/b.                     (5)

This formal evaluation is only a calculation of the integer word coefficient;
it is not asserted to be a physical source profile.

Split the cyclic word into maximal runs of ones. For a run of j>=1 ones
followed by a zero, let v be the periodic cover value at that following zero.
At the successive one positions, the deficits 1-u have sum
(1-v)(1-2^(-j)). The adjacent-one subtraction contributes j-1. Consequently
this run contributes

    1-(1-v)(1-2^(-j))=v+(1-v)2^(-j)>0

to Delta. Every run is positive and there is at least one, so K>0.
Integrality follows from (4). For evenness, apply the actual word to the
physical unit-deletion seed t=(1,0,2,1). Both signed responses are odd and
d=1, so their difference K is even. In particular K>=2.

The simplest clocks are

    a/b=2/3: W=OZ, K=2;
    a/b=4/5: W=OOZZ, K=2;
    a/b=6/7: W=OZZ, K=4.

## 3. Converting modular source freedom into an exact integer witness

We use the already reviewed STRONG_SOURCE_INTERFACE theorem: for every even
M>=2 and every pair of odd residues D_0,U_0 modulo M, there are arbitrarily
large actual strong unit-deletion sources with

    D=D_0, U=U_0 (mod M), s=E=d=1.

Its literal sources are converses of one-apex gadgets, have transitive root
deletion, and satisfy D>=L for arbitrarily large gadget indices L. This is
physical modular freedom, not assumed independent exact freedom.

Fix an arbitrary odd target r. Choose odd residues D_0,U_0 modulo K satisfying

    2a U_0-bD_0=r (mod K).                      (6)

Such residues always exist. Write D_0=1+2x and U_0=1+2y. After dividing the
even difference by2, (6) is

    2a y-bx=(r-2a+b)/2 (mod K/2).

Since gcd(2a,b)=1, Bezout gives x,y, even as an exact integer equality before
reduction. The modular source theorem now supplies actual D,U with those
residues and arbitrarily large D. Choose one with

    (2a-b)D>=r+K.

Append one O. Its exact signed response is

    S=2aU-bD >= (2a-b)D >= r+K,
    S=r (mod K).

Therefore h=(S-r)/K is a positive integer. Append W exactly h times. By (4)
its final signed response is precisely r. Its deletion sidecar remains1;
its final letter Z gives E=1 and strongness. Root-start count is still one,
and the root deletion is still transitive. This proves the exact theorem.

This construction is a finite, explicit physical recipe: a compiled actual
modular gadget, one O, then a specified positive number of repetitions of a
specified O/Z word. Large repetitions can be recorded in compressed form;
no efficient graph-order or positive-count bound is asserted.

### An equivalent arithmetic image of the one-apex sequence

Let g(1)=1, g(2n)=g(n), and g(2n-1)=g(n)+n-(n mod2) for n>=2.
The exact O/Z orbit of t has profile

    (D,E,U,d)=(g(n),n mod2,n,1), n>=2.           (7)

Indeed O sends n to2n and Z sends n to2n-1, with precisely these g recurrences.
Reversing the two maps reaches2 from every n>=2. Odd n>=3 end with Z and are
strong. The modular seed V(L) also has this form with n=2L-1. Thus the theorem
has the concrete arithmetic corollary

    {a n-b g(n): n>=3 odd} = all odd integers
    exactly when b/2<a<b                        (8)

under its stated coprime even-numerator/odd-denominator hypotheses. Every
integer in this arithmetic image retains the actual graph recipe above.

## 4. Sharpness and the unsolved positive readout

If 0<a<=b/2, (2) gives

    aU-bD <=(2a-b)D-2aE<0.

If a>=b, U>=D gives aU-bD>0 (equality a=b is incompatible with the stated
coprime even/odd pair). Thus outside the open interval the signed image cannot
contain both signs, and cannot be all odd integers. This proves sharpness.

The parity/denominator hypotheses should not be silently dropped. For example,
the interior slope3/4 does not have all odd integers in the arithmetic image.
Write n=2m-1 odd. The recurrence gives

    3n-4g(n)=2m-4g(m)+4(m mod2)-3.

If m is even this is1 modulo4. If m is odd and m>=3, the strong cone for (7)
gives g(m)>=m/2+1, so the expression is at most-3. The excluded m=1 gives-1.
Therefore every positive value3n-4g(n) is1 modulo4; the target3 is absent.
This is a hostile control on a broader rational-slope claim, not an exception
to (8).

The theorem prescribes one signed response only. It does not prescribe D or U
independently, nor a positive count such as alpha D+beta U+gamma. The repetitions
in Section3 can make both D and U enormous while their signed combination is
small. Consequently no all-odd tournament-spectrum, C2 source-pair coverage,
or positive scalar descent conclusion is drawn.

## 5. Every odd offset has actual C2 parent/child realizations

The signed theorem has a direct physical source-pair consequence. Fix
 t=(1,0,2,1). For a strong unit-interface source q=(D,1,U,1), the reviewed
N and M constructions give actual C2 graphs with counts

    n=N(t,Oq)=3D+12U-2,
    m=M(t,Jq)=13D+18U-2,
    m-2n=7D-6U+2.                              (9)

For every prescribed odd offset delta, apply the theorem at a/b=6/7 and
r=2-delta. It produces an actual q with 6U-7D=2-delta, so (9) gives

    m=2n+delta.

There are infinitely many such parent/child count pairs for each delta:
the initial modular seeds have unbounded D, all physical O/Z moves preserve
or increase D, and n>=3D. In particular both binary offsets delta=-1 and+1
have infinite exact physical C2 families, not just congruence families.
Small literal controls are

    q=t . Z O Z=(9,1,11,1):   n=157, m=313;
    q=t . O Z Z=(11,1,13,1):  n=187, m=375.

These cover offsets with a variable parent, not an independently prescribed
parent. In fact eliminating U in (9) gives

    n=17D+2-2delta, hence n=19-2delta (mod34).

For delta=-1/+1 this construction requires parents21/17 modulo34 respectively.
It therefore does not give both children of every parent or universal scalar
coverage. The N-to-M gate change replaces the source structure; no induced
extension of a fixed parent graph is asserted.

## 6. Exact verification and dependencies

Run python3 verify.py and python3 -O verify.py. The checker constructs the
periodic words for all coprime allowed pairs with b<=31, checks the raw transfer
on a coordinate basis, computes the periodic run contribution with rational
arithmetic, and verifies positivity and evenness of K. It uses the physical
joint modular compiler to construct compressed exact recipes for representative
odd targets, checking both modular interfaces and the final exact response.
It also directly builds small physical clocks and evaluates their complete
source profiles by independent subset path counts, then constructs both C2
parent/child controls in Section5 and checks their Hamilton and cycle counts. The finite tests supplement
the all-order proof and do not establish it by saturation.

Dependencies are the previously reviewed positive physical source grammar,
STRONG_SOURCE_INTERFACE, JOINT_MODULAR_GADGET, and BILINEAR_C2, copied unchanged;
the bilinear graph constructor and subset counters are also unchanged. The only
new use of modular freedom is followed by the exact decreasing word (4), with
an actual unit sidecar throughout.
