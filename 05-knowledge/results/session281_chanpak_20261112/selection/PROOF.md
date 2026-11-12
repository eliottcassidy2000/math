# Seven physical letters: full congruence coverage and zero-density integer support

**Session281. PROOF DRAFT under independent review. Chan–Pak remains OPEN.**

## 0. Scope and inheritance

Use the [Session280 active-block law](../../session280_arbitrary_switching_20261112/ancestry/PROOF.md):

    M_w(A,W,d)=(hA+bW-d, aA+tW, A),
    t=2^|w|, a=value_2(w)+1, b=1+sum_i(1-w_i)2^i,
    h=g(value_2(w)),
    g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+(v mod2).

Our alphabet is exactly L={10,100,110,1000,1010,1100,1110}; the initial
physical source is the rooted singleton s=(1,1,1). Positive words append
actual one-apex blocks and reverse the seam between successive roots.
Inverses below are finite-ring algebra, NEVER inverse tournament operations.

**Theorem.** For every m>=1 and every residue r containing an odd integer
(all r for odd m, odd r for even m), an L-program from s has H=r modulo m.
For odd m its full packet orbit is exactly

    Prim(m)={(A,W,d) mod m: gcd(A,W,d,m)=1}.

Yet the number of its distinct actual counts <=X is O(X^(399/400)), hence
has zero density among odd integers. The latter is inherited directly from
[Session280's63-letter scarcity theorem](../../session280_arbitrary_switching_20261112/coverage/PROOF.md),
because L is a subalphabet.

**Priority and new scope.** Congruence completeness without integer coverage
is ALREADY proved for fixed-port arithmetic states in
[Session237](../../session237_chanpak_synthesis_20261109/SYNTHESIS.md), crediting
Sessions195–206 and THM-4888/5085. [Session225](../../session225_joint_modular_law_recovery_20261108/PROOF.md)
already gives an all-odd-modulus limiting law for an exactly-two-J ensemble.
The new result here is the fixed SEVEN-letter language, its full primitive
packet orbit, a simultaneous dyadic/odd synchronization proof, and its
coexistence with a proved power-sparse bound. No claim to discovering the
general local/global distinction is made.

The distant sources are [THM-2200's support-hole trichotomy](../../../../01-canon/theorems/THM-2200-convex-semigroup-and-finite-place-support-hole-trichotomy.md),
[THM-4077's Mahler tangent chart](../../../../01-canon/theorems/THM-4077-mahler-denominator19-2adic-tangent-full-shift-isometry.md),
and the Sun additive-bases local/global separation in CURRENT-FRONTIER.
The explicit connection is physical word -> finite residue orbit: it retains
residue and loses absolute height. No theorem about those other systems is
imported. The live concepts are physical initialization, inverse-word cost,
primitive packet, dyadic attractor, CRT synchronization and absolute target.

## 1. Full primitive-packet orbit at every odd modulus

Put T=M_10, B=M_100, C=M_1000 and E_ij(c)=I+c e_ij. Exact identities in
GL_3(Z[1/2]) are

    T=[[3,3,-1],[3,4,0],[1,0,0]],
    B=[[5,7,-1],[5,8,0],[1,0,0]],
    U=M_1010 C^-1=[[1,-1/4,9/4],[0,1,2],[0,0,1]],
    V=M_1100 C^-1=[[1,-1/8,25/8],[0,1,4],[0,0,1]],
    [U,V]=E_13(-3/4),
    U²V^-1=I-(3/8)e_12+(19/8)e_13,
    UV^-2=I-6e_23-(7/2)e_13,
    T^-1 E_13(c) T=E_31(-c).                         (1)

The commutator convention is UVU^-1V^-1. In particular E_13(3) is present.

### 1a. Moduli coprime to six

For gcd(m,6)=1, powers of E_13(-3/4) give every E_13(c). Remove the central
e_13 terms in the next two matrices of(1). Their remaining coefficients
-3/8 and-6 are units, so their powers give every E_12(c),E_23(c).
Conjugation gives E_31(c); elementary commutators give the other directions.
Thus the generated group G_m contains SL_3(Z/m).

Elementary generation/transitivity here can be proved by CRT to prime-power
rings: a primitive column has a unit pivot; determinant-one row operations
eliminate the other entries and normalize the pivot. The same elimination
works for determinant-one matrices. CRT localizes every elementary coefficient.

### 1b. The prime three

`certificate.json` gives positive L-words whose matrices modulo3 are each of
the six E_ij(1). Multiplying those finite words is the base certificate that
the reduction contains SL_3(F_3). No enumeration of an unbounded group enters.
Eight matrices G e_13 G^-1 span sl_3(F_3), for chronological G-words

    (), (10), (100), (110), (1010),
    (100,10), (100,100), (100,1010).                  (2)

The rank-eight certificate is checked exactly. Since E_13(3) is present,
so is I+3G e_13 G^-1. The nilpotent squares to zero, so its3^(j-1)-st power
is EXACTLY I+3^j G e_13 G^-1. Their products realize every correction
I+3^j X modulo3^(j+1) with tr(X)=0. This is precisely the determinant-one
condition at that level. Induction therefore supplies the full kernel of
SL_3(Z/3^k)->SL_3(F_3).

Base lifts must have determinant ONE, not merely determinant one modulo3.
Use commutators: SL_3(F_3) is perfect, since for distinct i,j,l,
E_ij(a)=[E_il(a),E_lj(1)]. Express a base element in these commutators and
lift their factors to the seven-generator group. The lifted commutators
have determinant exactly one in Z[1/2]. The preceding trace-zero corrections
then lift to any desired SL_3(Z/3^k) element. Hence G_(3^k) contains SL_3.

### 1c. Simultaneous odd CRT

Let m=3^k q, gcd(q,6)=1, k>=1. The group first contains a localized SL_3(q),
identity on3^k. Raise the last two unipotents of(1) to powers divisible by
3^k. Their noncentral coefficients are units on q and zero on3^k. Remove
the central terms using powers of E_13(3) whose exponents are divisible by
3^(k-1); coefficients can be freely matched on q. Conjugation and elementary
commutators then give every localized elementary matrix.

Lift an SL_3(F_3) base element by determinant-one commutators as above, then
erase its q-part by that localized SL_3(q). Higher3-adic corrections use
conjugates of E_13(3q): they are identity on q and retain the rank-eight span
at3. This supplies localized SL_3(3^k), hence SL_3(Z/m).

Now s is primitive, every M_w has determinant2^|w|, and SL_3(Z/m) is
transitive on primitive columns. Every inverse generator in the finite group
is a positive power of itself. Therefore the POSITIVE physical orbit is
exactly Prim(m), proving the odd-modulus assertion. This conversion makes
no assertion about the actual height of its physical witness.

In fact every exact packet in this orbit has gcd1: an odd prime common to
its coordinates would pull back through the power-of-two determinant to a
common divisor of s; factor2 cannot divide the all-odd packet.

## 2. Two exponentially sparse scalar sequences fill every dyadic residue

B has characteristic polynomial f(z)=z³-13z²+6z-8, reducing to z²(z-1)
modulo2. Its unique odd2-adic root lambda lifts one bit at a time, because
f' is odd at the odd root. Direct evaluation gives lambda=3 modulo8. Write

    f(z)=(z-lambda)Q(z),
    Q(z)=z²+(lambda-13)z+8/lambda.

Q(lambda)=f'(lambda) is odd. Over R=Z/2^k, polynomial Bezout splits the
B-module into its lambda-line and its Q-primary part. On the latter B² maps
into2 times that module, so B^(2k) kills it. On the former a first-coordinate
normalized generator is

    v=(1,5/(lambda-8),1/lambda),                      (3)

whose coordinates are all odd. The projector is P=Q(B)/Q(lambda); it has
rank one and P²=P. Thus B^r=P whenever r>=2k and lambda^r=1 modulo2^k.

For the physical seeds s0=s and s1=T s=(5,7,1), the amplitudes
alpha_i=(P s_i)_A obey

    alpha_0=1 mod8, alpha_1=7 mod8.                  (4)

Directly B^4 modulo8 gives this projector and these amplitudes. Beyond the
transient,

    (B^j s_i)_A=alpha_i lambda^j mod2^k.              (5)

For k>=3, lambda=3mod8 has order2^(k-2): the elementary valuation identity
v_2(lambda^(2^j)-1)=j+2 follows inductively from v_2(lambda²-1)=3.
Its powers are the index-two subgroup of odd units whose residues modulo8
are1 or3. The two amplitudes(4) occupy different cosets. Thus(5) for the two
seeds covers every odd residue modulo2^k. The k=1,2 cases follow directly
from lambda=3mod4. Both scalar sequences have only O(log X) terms below X.
Their dyadic density is consequently not an integer-tail theorem.

## 3. A finite-monoid sandwich synchronizes both parts

Separate prime-power claims would not prove simultaneous CRT coverage.
Fix m=2^k q, q odd, k>=1. Choose r>=2k divisible by the order of B modulo q
and by the order of lambda modulo2^k. The actual positive program E=B^r is

    E=P modulo2^k, E=I modulo q.                     (6)

Consider the finite residue space

    Omega={u v: u odd modulo2^k} x Prim(q).

Each sandwich E M_w E restricts to a PERMUTATION of Omega. On the dyadic
line it is multiplication by an odd scalar mu_w; on q it is M_w. The scalar
is odd because both P and every active M_w preserve the all-odd packet
modulo2. Let H be the finite group generated by these permutations. Its
dyadic projection is abelian, and its q-projection is G_q.

Therefore every commutator in H fixes the dyadic coordinate, while the
q-projection of [H,H] is [G_q,G_q]. This contains SL_3(Z/q), because that
group is perfect by the elementary commutator identity already used above.
Thus one can act transitively on Prim(q) while fixing the selected dyadic
state. Inverses in this compressed finite permutation group are positive
powers of SANDWICHES, not of the ambient singular dyadic matrices.

First choose i,j in Section2 to obtain the desired dyadic scalar, with j
past the transient. Append E, which fixes that dyadic state. The q-packet is
primitive. Apply a positive sandwich program representing the required
commutator element to send it to(r,1,1) modulo q without changing the dyadic
coordinate. This gives the desired actual physical H modulo m. The q=1 and
k=0 cases were handled separately. This proves full congruence coverage.

## 4. Consequence and the missing height coordinate

The SAME seven-letter language is now provably congruence-complete and
power-sparse. For every finite collection of fixed moduli, combine them by
lcm: every admissible odd residue is occupied. Such tests cannot certify
any particular odd target's absence in this language, although infinitely
many odd targets are absent. Even residues congruent to7 and21 occur, while
the exact counts7 and21 are excluded for every tournament.

This does not forbid a modulus depending on the target together with an
actual height bound. A potential positive transfer would be: produce a
congruence witness modulo a modulus larger than a prescribed interval, AND
bound its true count inside that interval. Our inverse-to-positive-power
conversion does not pay that height bill. Nor does the finite-ring projector
E denote an actual count-preserving graph operation: its real count grows.

The strongest next lesson is therefore operational rather than philosophical:
a residue compiler can now be used as a benchmark which is guaranteed to
pass every fixed modulus but fail an eventual-tail goal. Any proposed
Chan–Pak selector needs a genuinely height-sensitive ingredient beyond it.

## 5. Reproduction

`verify.py` uses standard-library exact arithmetic and explicit exception
gates. It checks rational identities, six elementary words modulo3, the
rank-eight lift, dyadic roots/projectors/cosets, complete small odd packet
orbits, simultaneous residue witnesses and literal physical controls.
`certificate.json` supplies finite certificates, and `results.json` records
controls. Sections1–3, not finite testing, establish the all-modulus result.
