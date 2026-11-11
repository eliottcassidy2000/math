# A unary dyadic controller, Collatz carry readout, and odd-prime digit lifts

**Status: PROVED, conditional only on the inherited physical F3/F5 transfer law; separate FINITE-EXACT controls below. Chan–Pak remains OPEN.** No canon promotion or external priority claim is made.

## 1. Actual source and the question

Use the initialized marked transitive pair with state

    s0=(r,M,q,E)=(-1,3,0,1),     q=H(B)+r.

The inherited THM-5130 transfer law, as specialized in Session264/control/PROOF.md, supplies actual positive tournament operations F_a, a=3,5:

    (r,M,q,E) -> (-r,(a²+1)M+a q−E,(2a−2)M,2r+aM−q).

F5 adds four vertices and moves the selected pair. Closing a source with a transitive tail of t>=1 vertices and one closer gives the ACTUAL scalar count

    H=r+2^t M.                                                   (1)

We ask for digit control INSIDE this initialized grammar. Session199/profinite-source-escape/PROOF.md already supplies scalar modular density for every modulus in another architecture. Thus neither scalar congruence density nor its Collatz corollary below is claimed as new globally. The new mechanism is a single physical letter whose fixed-sign iterations enumerate every odd binary slope residue, with an explicit bit compiler.

Concept board: positive source words / contracting versus unit 2-adic coordinates / exact offset sign / Collatz affine carry / ordinary integer height. The last coordinate remains the unresolved bridge to integer target coverage.

## 2. Main theorem: one letter controls all binary slope digits

Write F for the four-by-four matrix of F5 in the displayed state order and G=F²:

    F=[-1   0   0  0]       G=[ 1   0   0   0]
      [ 0  26   5 -1]         [-2 711 131 -26]
      [ 0   8   0  0]         [ 0 208  40  -8]
      [ 2   5  -1  0]         [-2 122  25  -5].

**Theorem D1 (unary binary controller).** For every k>=1, every odd residue m modulo 2^k, and either exact sign r*, an actual F5-only word starting from s0 has offset r* and slope M=m modulo 2^k. More precisely put e=0 if r*=−1 and e=1 if r*=+1, and b=ceil(k/6). Among the words

    F5^(2(b+j)+e),       0<=j<2^(k−1),                         (2)

EVERY odd slope residue appears exactly once. The source order is at most

    8b+2^(k+2)−6+4e.

A bit-by-bit algorithm finds j without enumerating these words.

### Proof of the 2-adic decomposition

Direct characteristic-polynomial calculation gives

    char_G(t)=(t−1)²(t²−745t+64).                            (3)

The quadratic has, by the elementary simple-root lifting argument modulo2, one even root mu and one odd root lambda in Z_2. Since mu(mu−745)=−64 and mu−745 is odd, v2(mu)=6. Also lambda=745−mu. The factors t−mu and (t−1)²(t−lambda) are coprime over Z_2, because 1−mu and lambda−mu are units. Bezout's identity therefore splits Z_2^4 as the direct sum E plus K, with integral commuting projector P onto E, where G=mu on K and G is invertible on E. This is just polynomial primary decomposition, not a real or complex diagonalization.

Modulo2 the displayed G is idempotent. Its image is the rank-three subspace q=0, and it acts as the identity there. Reducing the direct-sum decomposition modulo2 identifies E/2E with that image. Consequently

    G|E=I+2A

for an integral endomorphism A of E. In particular G²|E=I+4B, B=A+A². Repeated squaring gives, for j>=1,

    (G^(2^j)−I)/2^(j+1) = B  (mod2) on E.                  (4)

Indeed (I+4B)^(2^(j−1)) has first nonconstant term 2^(j+1)B and all higher terms vanish modulo2^(j+2). The j=0 normalized increment is A. Both A and B commute with G, whose reduction on E is the identity.

### Two small calculations supply the initial odd digit

For x_e=P F^e s0, e=0,1, the following computations modulo8 suffice:

    e    x_e mod8          (G−I)x_e mod4   (G²−I)x_e mod8
    0    (7,3,0,1)         (0,2,0,2)       (0,4,0,4)
    1    (1,5,0,5)         (0,2,0,2)       (0,4,0,4).

They require no unknown 2-adic constants: G^4=P modulo8, since (I+2A)^4=I modulo8 on E and mu^4=0 modulo8 on K. Thus use G^4 F^e s0 to compute the table. The M coordinates of Ax_e and Bx_e are odd. Commutation and G=I modulo2 show that the same holds at every G^n x_e. It follows that for every n>=0 and j>=0,

    v2( M(G^(n+2^j)x_e)−M(G^n x_e) ) = j+1.             (5)

This is an exact valuation statement on the projected 2-adic orbit, not on an arbitrary formal source.

The M coordinate is always odd. Equation (5), first for the lowest bit and then by induction, shows that any block of 2^(k−1) consecutive iterates G^n x_e contains every odd M residue modulo2^k exactly once. To transfer to actual integer tournaments, observe

    G^b F^e s0−G^b x_e = mu^b (I−P)F^e s0 =0 mod2^k

for b=ceil(k/6). Thus the literal words (2) have exactly those residues. Their vertex bound follows from two initial vertices plus four per letter. This proves D1.

The bit compiler fixes the sign e and burn-in b. Start j=0. At bit i=1,...,k−1, compare the current M with the requested m modulo2^(i+1). If they differ, add 2^(i−1) to j. Equation (5) says that this toggles precisely the new bit without changing any lower one. Fast modular matrix powering evaluates each query. The compiler outputs the positive run length, exact state, and, optionally, literal adjacency lists. The exponent is not an inverse operation and does not remove old vertices.

**Efficiency boundary.** This is O(2^k) actual vertices, exponential in the number of binary digits. The states themselves generally grow exponentially in the word length. It is a strong exact digit controller, not a small-height integer selector.

### The projected limit is a binary adding machine

**Corollary D1a (exact time-to-slope coordinate).** For either sign e, the map

    n -> M(G^n x_e),     n in the nonnegative integers,

extends to a bijection from Z_2 to the odd elements of Z_2, and

    v2(M(G^n x_e)−M(G^m x_e))=1+v2(n−m).                    (7)

Here n and m may be any distinct 2-adic integers in the extension. Thus `(M−1)/2` is an isometric coordinate for time, and G acts as the binary adding machine in that coordinate. This is a theorem about the projected limit system. No actual finite tournament with a noninteger 2-adic state is postulated.

*Proof.* For ordinary integers write n−m=2^j u with u odd. Telescope u increments from (5). After division by2^(j+1), each increment is odd; an odd number of such terms remains odd. This proves (7). The resulting isometry extends uniquely by completeness. Surjectivity follows because D1's projected blocks fill every odd residue at every level, and compatible nested time residues have a limit. Addition by one is conjugated by this coordinate map, giving the claimed adding machine. Literal finite sources approach these projected coordinates after the explicit burn-in; the projection itself is never used as an illegal physical operation.

## 3. Collatz's carry becomes a literal readout

Let U(N)=(3N+1)/2^v with v=v2(3N+1), for odd positive N. For an odd source slope M, the closing law yields exact cases:

- r=−1, t>=2: v=1 and U(2^t M−1)=3·2^(t−1)M−1.
- r=−1, t=1: v=1+v2(3M−1).
- r=+1, t=1: v=1.
- r=+1, t=2: v=2+v2(3M+1).
- r=+1, t>=3: v=2 and U(1+2^t M)=1+3·2^(t−2)M.

Thus only the two faces (r,t)=(−1,1),(+1,2) can carry arbitrarily deep halving. This is a precise dyadic valuation statement, not a claim that Collatz is a physical tournament operation.

**Corollary D2 (prescribed halving depth in a unary physical family).** Every integer h>=2 occurs as v2(3H+1) for an actual count H from a negative-unit source produced solely by F5. To construct it prescribe

    M=3^(-1)(1+2^(h−1)) modulo2^h

using D1 at r=−1, then close with t=1. The congruence says 3M−1=2^(h−1) modulo2^h, so the valuation is exactly h−1, not merely at least h−1. The depth h=1 is obtained with any negative-unit source and t=2. Dually every h>=3 occurs on positive sources with t=2 by prescribing

    M=3^(-1)(2^(h−2)−1) modulo2^(h−1).

This does not show U(H) is itself a Hamiltonian count, nor that iteration follows a tournament grammar, nor that Collatz converges. Session199's global scalar congruence result already implies unrestricted existence of these scalar valuations; the new result is the unary source realization and explicit compiler.

### Why 7, 21, and 63 are different objects here

The arithmetic identities

    7=2^3−1,   21=(2^6−1)/3,   63=2^6−1,
    3·7=21,    3·21=63,        3·21+1=64

compare two operations: plain tripling and tripling-with-carry. Ordinal joining with a triangle implements plain tripling of an ALREADY EXISTING Hamiltonian count. It cannot start from the excluded scalar value7. The Collatz carry sends21 directly to1, whereas63 lies on the negative dyadic face and first goes to95. The inherited THM-1370 exclusions of7 and21 and THM-4748 Mersenne family at exponent>=4 are inputs, not rechecked here. Therefore63 being a count does not require7 or21 to be counts.

D1 gives a third, distinct role: all of7,21,63 are legal odd slope residue targets. A slope residue is neither a literal slope equality nor a scalar count. The controls compile these targets with both signs and report actual integer slopes, keeping those distinctions visible.

## 4. The odd-prime controller lifts one digit at a time

This extends Session264's prime-field full-state theorem, rather than rechecking its accepted transfer law. Work in signed coordinates y=r(M,q,E), where F_a is an affine map f_a(y)=C_a y+c and each letter flips the exact sign. Let v=(1,−8,−13).

**Lemma O1 (affine digit lift).** Let p be odd and g(y)=L y+t an integral affine map with L=I modulo p. Then for every s>=0,

    g^(p^s)(y)=y+p^s t (mod p^(s+1))                         (6)

for all y, and its linear part is I modulo p^(s+1).

*Proof.* In homogeneous coordinates write g=I+N. The linear block of N is divisible by p and its last row is zero, so N^j is divisible by p^(j−1), and its linear block by p^j. For 2<=j<=p^s,

    v_p(binomial(p^s,j))+j−1 >= s−v_p(j)+j−1 >= s+1

because p is odd. Thus only the j=1 term survives modulo p^(s+1); its linear block vanishes. The case s=0 is the hypothesis itself. The expression p^s t may be reduced using only t modulo p. ∎

**Theorem O2 (every nonexceptional prime power).** For every odd prime p other than5 and761, every k>=1, either exact offset sign, and every desired (M,q,E) modulo p^k, a positive physical F3/F5 word reaches that target from ANY actual unit source.

*Proof.* Session264 proves that g=f5^(p²−1) has L=I modulo p and t=v/10 modulo p, and that the fixed-sign conjugate translation vectors

    v, C3v, C3²v, C5C3v

span F_p³: two indicated determinants are344 and632, with gcd8. It also establishes full reachability modulo p. Suppose a word reaches the requested state modulo p^(k−1). By O1, g^(p^(k−1)) fixes all lower digits and translates the highest digit by v/10; its conjugates translate by the four spanning vectors. Their powers correct all three highest digits. Every inverse used to form a conjugate is replaced by a positive power in the FINITE permutation group on the full residue state plus exact sign modulo p^k. This is valid because both generator determinants are powers of2 and p is odd. Such a positive inverse word has the correct exact sign parity; it is not an inverse over the integers. Induction proves the result. A shortest full-state BFS word has at most2p^(3k)−1 letters, hence adds at most8p^(3k)−4 vertices. ∎

The primes5 and761 remain exceptions to THIS algebraic lift proof. The earlier finite p=5 orbit does not automatically justify all powers5^k.

## 5. A mixed 3-and-7 tower, rather than an unjustified CRT inference

Prime-power transitivity alone need not give simultaneous controllability for a product. Here the common physical F5 clock supplies a direct bridge.

**Theorem O3.** Full-state control at either exact sign holds for EVERY modulus 3^a7^b, a,b>=0, a+b>0. The case a,b>=1 uses one FINITE-EXACT base orbit at modulus21, followed by an algebraic induction.

*Proof.* The single-prime cases follow from O2. The script records an exact BFS over all2·21³=18522 states and every transition, giving the base a=b=1 from s0. Since the positive generators are permutations of this finite state space, a full orbit is strongly connected; hence the base works from any actual unit source.

For induction, the affine map f5^(q²−1) has linear part I modulo q and translation v/10, for q=3,7. Formula(6) implies that f5^((q²−1)q^d) is the identity modulo q^d.

To lift the 3-adic exponent from a−1 to a, a>=2, while preserving 7^b, use

    f5^L,     L=48·3^(a−2)·7^b.

As a power of g3=f5^8, its exponent is 2·3^(a−1)·7^b, so it is identity-linear modulo3^a and its highest-digit translation is

    3^(a−1) (2·7^b/10) v modulo3^a,

with nonzero coefficient modulo3. As a power of g7=f5^48 its exponent is3^(a−2)7^b, so it is the identity modulo7^b.

To lift the 7-adic exponent from b−1 to b, b>=2, while preserving3^a, use

    f5^L,     L=48·7^(b−1)·3^a.

The highest 7-digit translation is 7^(b−1)(3^a/10)v, nonzero modulo7, and the map is the identity modulo3^a. Conjugate and span exactly as in O2, replacing inverses in the finite group modulo the WHOLE product. Conjugation preserves identity on the factor being protected. All these F5 exponents are even, so the exact sign is unchanged. This independently corrects each new prime digit and proves O3.

Thus7,21,63 are not merely matching numbers: they lie on a concrete physical residue-control tower, with63=3²·7 the first nontrivial 3-digit lift of the mixed base21. This does not turn those residues into small exact counts. The controlling words can be very long.

## 6. Connections and remaining obligations

- **Collatz carry -> source digit readout:** map H=r+2^t M; preserve exact halving valuation; lose the source graph under U. Needed sidecar: a physical slope-tripling/carry operation with controlled size. No such operation is claimed.
- **Contracting p-adic eigendirection -> positive unary word:** the projection P is a proof device only; the actual constructor uses the integer burn-in and positive F5 run. Exact sign, modulus and a run-length certificate are retained.
- **Affine odd-prime translation -> physical high-digit lift:** conjugates act only modulo a specified finite ring; the lost coordinate is ordinary height. A finite inverse word expands the graph.
- **7/21/63 -> mixed control tower:** preserve entire source residue state, not literal count equality. The mod21 base is finite evidence; the subsequent induction is algebraic.

Open: Can the source realizing a prescribed digit string be chosen with near-minimal integer M, rather than enormous M in that residue? Can the Collatz critical faces be joined by an actual small source operation, not a numerical substitution M->3M? Can the isolated odd-prime exceptions be removed without a huge full-group enumeration? None of the present results resolves Chan–Pak or Collatz.
