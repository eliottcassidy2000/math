# Exact square-plus-one intersection of minimum branching defect

**PROVED corollary of the independently reviewed minimum-defect theorem.**
The primality guard remains explicit. This is a structural connection to
Fermat-prime candidates, not an inference from two collections having the
same cardinality.

## 1. The all-degree intersection

Fix a strong rooted tournament (T,r) with d>=3 in-neighbours and the
minimum possible branching defect b=F_r(2)-H(T)=2^d-2. The classification
in [MINIMAL_DEFECT](../metagraph/MINIMAL_DEFECT.md) says that T is exactly

    r -> TT_k -> I -> r, k>=1,

where I=TT_d, except that I=C3 is also allowed when d=3. The transitive
case has H=(2^d-1)(2^k-1)+2; the triangle case has H=9*2^k-3.

**Exact intersection.** Such a rooted source has H=n^2+1 for an integer
n if and only if I=TT_d, k=1, and d is even. Necessarily d>=4, and then
n=+/-2^(d/2).

Proof: if d>=3 and k>=2, both powers of2 vanish modulo4, so

    (2^d-1)(2^k-1)+2 = 3 (mod4),

whereas a square plus1 is1 or2 modulo4. If k=1, H=2^d+1 and n^2=2^d
holds exactly when d is even. The C3 alternative is divisible by3, whereas
n^2+1 is never divisible by3. These arguments exhaust the classified
sources. This is an iff about sources in one exact stratum, not about all
tournaments representing a near-square number.

## 2. Primality singles out the Fermat candidates

For even d write d=2^j u with u odd. If u>1, the identity

    x^u+1=(x+1)(x^(u-1)-x^(u-2)+...-x+1)

with x=2^(2^j) proves2^d+1 composite. Hence a prime value in section1
must have d=2^j, j>=2. Conversely each prime Fermat candidate

    F_j=2^(2^j)+1, j>=2,

has the displayed minimum-defect source on2^j+2 vertices, square root
2^(2^(j-1)), and root defect F_j-3. Thus the prime members of this
classified source family are exactly the Fermat primes with index>=2.
No assertion is made that there are infinitely many, or that the examples
below exhaust the prime members.

| j | d | H | Square root of H-1 | Source order | Root defect |
|---:|---:|---:|---:|---:|---:|
|2|4|17|4|6|14|
|3|8|257|16|10|254|
|4|16|65537|256|18|65534|

These three primalities are independently certified by the full-factorization
Lucas criterion already used in `quadratic_primes.py`. The smaller Fermat
number5 belongs to d=2, where the separate pentagon classification applies;
3 itself is not of the form n^2+1. Those boundaries explain why simply
identifying a five-element collection would hide a mathematical distinction.

There is also an exact map to the aliquot maximum boundary:

    H=2^d+1 prime  iff  A135244(2^d+2)=(2^d+1)^2.

The map retains the same integer H and a certified prime predicate. It
forgets the tournament's root and future allowed operations, and it does
not make the two metagraphs equivalent. The full source construction,
marked root and the primality certificate are the necessary sidecars.

## 3. The nontransitive residual sector at defect6 is excluded too

The completed two-inport, residual-C3 sector in
[the defect proof](../defect/PROOF.md) is the variable-prefix inflation of
the five ready tails. For a prefix P its exact counts are
E Q_P(2)+B H(P), with(E,B) equal to
(3,-1),(6,1),(5,-1),(13,3),(14,3). For
P=TT_a -> C3 -> TT_b, H(P)=3 and Q_P(2)=12*2^(a+b). Hence its five
count rays are

    36D-3, 72D+3, 60D-3, 156D+9, 168D+9,
    D=2^(a+b), a,b>=0.

They are all divisible by3 and therefore contain NO square-plus-one
integer, whether prime or composite. Together with section1 at d=3,
this removes both newly classified nontransitive/minimum-three-inport
sectors from any defect6 near-square target. The remaining two-inport
transitive-residual sector is a distinct obligation; it is not removed
by these arguments.

This is a useful exact filter for a search over source families. It also
shows why finding more fully classified families need not automatically
supply missing prime or square-plus-one targets: a whole source stratum
can be arithmetically excluded.

## 4. Reproduction and scope

`minimum_squares.py` checks the formulas for3<=d<=128,1<=k<=128, gives
independent literal first-vertex Hamilton counts on18 small collars, and
saves the three explicit Fermat-candidate source graphs and Lucas
certificates. The universal iff follows from the classification and the
modular proof above, not from those finite controls. Parent theorem proof
review: PASS by the arithmetic lane; no new theorem namespace is reserved.
