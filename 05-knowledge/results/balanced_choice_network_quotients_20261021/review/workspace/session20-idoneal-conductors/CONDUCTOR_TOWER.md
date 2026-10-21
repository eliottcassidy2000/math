# The exact Eisenstein class tower, two port quotients, and the first carry

21 October 2026. **PROVED:** the explicit dyadic class groups, extension maps, genus coordinates, and the indicated equivariance obstructions. **FINITE-EXACT:** the retained enumerations. The class-group formula and genus terminology are classical; no priority claim is made. The global 65-idoneal classification is not needed or changed.

## 1. Inheritance and the precise question

The inherited packet `05-knowledge/results/divisor_norm_ball_20261021/DIVISOR_NORM_BALL.md`, sections 5–6, already proves that the idoneal integers of the form 3m² are exactly 3,12,48, corresponding to conductors 2,4,8 in the Eisenstein field. Their class groups have orders 1,2,4 and exponent at most two. Its explicit warning that norm-ball points are not ideal classes remains essential.

The present addition is an actual system of canonical maps, not an identification of three numerical lists. It determines what the conductor tower can and cannot explain about two independent deletion bits, including the two different size-two port strata. The hostile boundary is conductor 16, where the missing coordinate is a carry, not an extra independent involution.

Write

    R = Z[ω],  ω²+ω+1=0,  j=1+2ω=sqrt(-3),
    O_f = Z+fR = Z+fZω,
    w = 1+4ω.

Only proper invertible ideal classes are used. In particular the noninvertible conductor ideal is not a class-group element.

## 2. A concrete residue presentation of every dyadic level

The standard conductor exact sequence, together with Cl(R)=1, gives

    G_f := Cl(O_f)
         ≅ (R/fR)^× / ((Z/fZ)^× · image(R^×)).                    (1)

For integral α with gcd(N(α),f)=1, the residue-to-class map in (1) is

    α mod fR  ↦  [αR ∩ O_f].                                   (2)

The required classical input is Conrad, *The conductor ideal of an order*, Theorem 3.8 and exact sequence (5.7), printed pp. 6–7, 16. The following explicit specialization and calculations are proved here. Source: https://kconrad.math.uconn.edu/blurbs/gradnumthy/conductor.pdf

For the dyadic conductors f=2^a considered below, the condition gcd(N(α),f)=1 is equivalent to odd norm. An element of R/fR is a unit exactly when its norm is odd, because X²+X+1 is irreducible modulo 2. There are 3f²/4 units. The subgroup of odd scalar units has f/2 elements; adjoining ω multiplies its order by 3. Thus

    |G_f|=f/2,  f≥2.                                          (3)

There is a particularly useful elementary equivalence test. If α≡1 mod2R, then its class in (1) is trivial precisely when α is congruent modulo fR to an odd scalar. Indeed its expression as a scalar times an Eisenstein unit must reduce to 1 in F4^×, so that unit can only be ±1. Equivalently, for α=A+Bω with α≡1 mod2R, the class is trivial iff f divides B.

### Theorem 1: explicit group and extension coordinates

For every integer a≥2,

    G_(2^a) = <[w]> × <[j]>
            ≅ C_(2^(a−2)) × C2.                              (4)

At a=2 the first factor is trivial. At a=1 the whole group is trivial. For a≥3 the canonical ideal-extension map is

    e_(2^a,2^(a−1)) : G_(2^a) → G_(2^(a−1)),
    (t,e) ↦ (t mod 2^(a−3), e),                              (5)

where (t,e) denotes [w]^t[j]^e. The map G4→G2 is trivial.

**Proof.** Certainly j²=−3, an odd scalar, so [j]²=1. Its coefficient of ω is2, so [j] is nontrivial modulo 4. Write

    w^(2^r)=A_r+B_rω.

Initially A_0=1 and B_0=4. Squaring gives

    A_(r+1)=A_r²−B_r²,
    B_(r+1)=B_r(2A_r−B_r).

Inductively A_r is odd and v2(B_r)=r+2: because B_r is divisible by 4, the second factor 2A_r−B_r has exact 2-adic valuation 1. The equivalence test therefore gives [w] exact order2^(a−2) in G_(2^a). Every power of w is scalar modulo 4, whereas j is not. Hence the two cyclic subgroups intersect trivially. Their product already has order2^(a−1), the full size in (3).

For g|f, extension [I]↦[IO_g] corresponds under (1) to residue reduction modulo g. One can see this directly on (2): the prime-to-f ideal αR contracts to each order, and extension-contraction away from the conductor gives

    (αR ∩ O_f)O_g = αR ∩ O_g.

This proves (5), not just its cardinality. QED.

As a classical specialization this recovers the dyadic part of the inherited idoneal theorem: exponent at most two holds exactly for f=2,4,8. The inherited odd-factor argument is still required to classify all even conductors, not just this dyadic tower.

## 3. Literal reduced forms and canonical kernels

Use primitive positive-definite forms (A,B,C) with discriminant B²−4AC=−3f² and reduced convention |B|≤A≤C, taking B≥0 on the boundary. Formula (2) gives the following exact labels:

| f | residue coordinate | reduced form |
|---:|---|---|
| 2 | 1 | (1,0,3) |
| 4 | 1 | (1,0,12) |
| 4 | j | (3,0,4) |
| 8 | 1 | E=(1,0,48) |
| 8 | w | A=(4,4,13) |
| 8 | j | B=(3,0,16) |
| 8 | wj | C=(7,2,7) |

The multiplication is A²=B²=C²=E and AB=C. Extension gives

    G8 → G4:  E,A ↦ (1,0,12);  B,C ↦ (3,0,4),
    G4 → G2:  both classes ↦ (1,0,3).

Thus the actual conductor kernel is

    K = ker(G8→G4) = {E,A}.                                   (6)

There are two distinct ways the numbers4,2,1 appear, and they must not be conflated:

1. The fibres of maps from the **common carrier G8** to G2,G4,G8 have sizes4,2,1. They are torsors for the respective kernels G8,K,{E}. Except for the whole group, a fibre is not invariant under every translation by G8.
2. The quotient G8-sets G8/{E}, G8/K, G8/G8 have sizes4,2,1 and point stabilizers {E},K,G8. These are the class sets G8,G4,G2, in the reverse conductor order.

A tournament deletion fibre acted on by both toggles is compared with (2), not automatically with (1). Claiming the same action for (1) changes types.

### Independent exact form computation

For α=a+bω with N=N(α) and gcd(b,N)=1, put

    r ≡ −fa b^(−1) mod N.

Then αR∩O_f has basis [N,fω−r]. Choose even B with B≡f+2r mod2N. Its associated form is

    (N, B, (B²+3f²)/(4N)),

which may then be reduced. This gives all displayed labels without choosing them to fit a desired group table.

There is also an independent extension computation. In a form of conductor f choose an equivalent representative with odd A, swapping A,C if needed. For g=f/2 choose the even B' with B'≡B/2 modA. The extension class is represented by

    (A,B',(B'²+3g²)/(4A)).

The checker compares these two paths for every class at conductors 2 through128 in powers of two.

## 4. Norm and genus characters give the two binary axes

For every dyadic conductor f=2^a with a≥3, the norm gives a well-defined homomorphism

    ν_f : G_f → (Z/8Z)^×,   [α] ↦ N(α) mod8.                  (7)

Odd scalar factors change norms by squares, which are1 modulo 8, and Eisenstein units have norm1. In the coordinates of Theorem 1,

    ν_f([w]^t[j]^e)=5^t 3^e mod8.                            (8)

Within this dyadic family, the map is surjective and its kernel is precisely the square subgroup G_f²=<[w]²>. Consequently it identifies the genus quotient with G8, and

    ker(G_f→G8)=G_f².                                        (9)

Calling the square quotient the genus quotient uses Gauss's principal-genus theorem, recalled in Kani, *Idoneal numbers and some generalizations*, section 2.2, equation (2.2), printed p. 200. Source: https://www.labmath.uqam.ca/~annales/volumes/35-2/PDF/197-227.pdf

At f=8 this is an isomorphism on the entire class group:

| class | norm mod8 | χ_−4 | χ_8 | χ_−8=χ_−4χ_8 |
|---|---:|---:|---:|---:|
| E | 1 | +1 | +1 | +1 |
| A | 5 | +1 | −1 | −1 |
| B | 3 | −1 | −1 | +1 |
| C | 7 | −1 | +1 | −1 |

Here χ_−4(n)=(−1)^((n−1)/2), and χ_8(n)=(−1)^((n²−1)/8), for odd n. These are multiplicative sign observables, not a choice of labels on four objects. The conductor quotient G8→G4 is exactly χ_−4 after identifying the nontrivial G4 class with−1.

The class B has additional arithmetic meaning: it is the contraction of the unique ramified prime jR above3. Thus arithmetic supplies two distinguished independent classes:

    A: nonidentity in the conductor kernel;
    B: ramified-prime class.

The second quotient χ_−8 has kernel{E,B}. It is not the actual extension map to O4. There is only one suborder of each conductor, so there are not two different conductor 4 targets hidden here.

### Scope correction: a non-dyadic hostile

The dyadic hypothesis in this section is necessary. An earlier draft said “every f divisible by 8,” improperly broadening the scope of the square-kernel assertion beyond Theorem 1. Although the norm map modulo 8 still exists for such f, its kernel need not be the square subgroup.

At f=40, the residue quotient has 24 classes, its norm-modulo-8 kernel has size 6, and its square subgroup has size 3. A concrete witness is α=1+8ω, whose norm is 57≡1 mod8. The additional well-defined class character (N(α)/5) is −1, so this class cannot be a square. Its reduced form is (25,20,52), of discriminant −4800. The retained checker independently exhausts this quotient and its squares. Thus (8)–(9), the two-character genus description, and the all-level genus-identification claim in this packet are explicitly restricted to f=2^a, a≥3.

## 5. A positive blue-predicate bridge, with its exact scope

The companion balanced-tournament calculation identifies each Hamilton path with a pair of coordinates, each a finite subset or the exceptional symbol ∞. A path is blue exactly when its two coordinates agree. Central-pair deletion forgets the newly inserted bit in each finite coordinate and fixes ∞. Let τ_L,τ_R toggle those two new bits when present.

Fix a **finite blue parent**, so its two old subsets agree. Its four lifts have new bits (ε,η)∈F2² and are blue precisely when ε=η. Choose the zero-bit lift as origin and define

    θ(ε,η) = B^ε C^η = w^η j^(ε+η).                         (10)

This is a group isomorphism from the two-toggle fibre to G8. It sends

    τ_L ↦ B,  τ_R ↦ C,  τ_Lτ_R ↦ A.

Under actual conductor extension G8→G4, the w factor vanishes and the j exponent survives. Consequently

    a lift is blue  iff  θ(ε,η) ∈ ker(G8→G4)
                    iff  ε XOR η = 0.                        (11)

The two blue lifts go to E,A, and the two black lifts to B,C. The three genus characters now have exact path meaning:

    χ_8(θ(ε,η))    = (−1)^ε,
    χ_−8(θ(ε,η))   = (−1)^η,
    χ_−4(θ(ε,η))   = (−1)^(ε+η).                            (12)

Thus a genuine operation-and-predicate-preserving local map exists. It transports both independent bit toggles and the blue/black refinement through the actual arithmetic extension, not merely a four-element cardinality. The conductor kernel is the **diagonal toggle**, not either individual shore toggle.

### What the local map needs and what it loses

For a finite black parent, its two old subsets are already unequal. Every lift remains black, regardless of the new-bit parity. The uniform formula is therefore

    child blue = (parent blue) AND (new bits agree).

The parent-coincidence predicate is necessary side data; conductor parity alone is not a global blue criterion. On exceptional fibres at least one new bit does not exist, so the four-element pointed-group formula does not apply unchanged.

More precisely each deletion fibre is one V4 orbit, with stabilizer

    {1}          if both parent coordinates are finite;
    <τ_L>        if the left coordinate is ∞;
    <τ_R>        if the right coordinate is ∞;
    <τ_L,τ_R>    if both are ∞.

Under (10), these are {E},<B>,<C>,G8, respectively. Hence the two size-two fibre types require the two different quotients G8/<B> and G8/<C>. They are detected by χ_−8 and χ_8 respectively. Neither is the literal conductor quotient G8/<A> in the blue-preserving marking.

No fixed group isomorphism can identify both size-two types with the one conductor quotient: an isomorphism cannot send two different stabilizer subgroups to the same subgroup. This is a structural obstruction, not a count disagreement.

A complete local orbit dictionary therefore retains the inactive-port type and uses

    G8,  G8/<B>,  G8/<C>,  G8/G8,

for the four stabilizer types. A zero-bit section supplies an origin where bits exist. Without a section a free fibre is a torsor, not a canonically pointed group. The conductor flag alone distinguishes A but leaves B,C interchangeable; choosing which is the left toggle is one remaining orientation choice. Its two possible choices are related by the unique nonidentity automorphism preserving that flag.

### Shore reflection and arithmetic symmetry

Shore exchange swaps ε,η. On G8 it becomes

    σ(A)=A,  σ(B)=C,  σ(C)=B,
    σ(w^t j^e)=w^(t+e)j^e.                                 (13)

It preserves the conductor kernel and commutes with extension to G4, while interchanging the two port characters χ_8,χ_−8. Thus it gives a valid symmetry of the three-level abstract class-group tower and of the blue-predicate bridge.

It is **not** induced by complex conjugation, which is inversion and is trivial on G8. Nor does it preserve the ramified-prime decoration: it swaps B, the class of the prime over 3, with C. A claim of a canonical ring/ideal symmetry therefore needs more than this group-level tower automorphism. At the next conductor the obstruction becomes stronger still: σ does not lift even to a commuting group automorphism, as proved below.

## 6. The conductor 16 carry is genuinely nonsplit

The full next-level table is:

| coordinate | conductor 16 form | image at conductor 8 |
|---|---|---|
| 1 | (1,0,192) | E |
| j | (3,0,64) | B |
| w | (13,8,16) | A |
| wj | (7,4,28) | C |
| w² | (4,4,49) | E |
| w²j | (12,12,19) | B |
| w³ | (13,−8,16) | A |
| w³j | (7,−4,28) | C |

Thus

    1 → <w²> → C4×C2 → C2×C2 → 1                        (14)

is the literal extension sequence. It does not split. Every lift of A is w or w³ and has square w²≠1; a homomorphic section would have to lift the involution A to an involution. Equivalently G16[2] maps only to{E,B}, not onto G8.

Writing t=t0+2t1 with bits t0,t1, multiplication by w sends

    (t0,t1,e) ↦ (1−t0, t1 XOR t0, e).

The new coordinate is a carry bit. Two applications change t1; four return to the start. Independent tournament bit toggles cannot satisfy that law while retaining the conductor extension map.

**Shore-reflection lift obstruction.** The automorphism σ of (13) cannot lift to any automorphism of G16 commuting with extension to G8. The subgroup G16[2] is characteristic, and its image is exactly {E,B}. A commuting lift would force σ to preserve {E,B}, whereas σ(B)=C. Thus the three-level blue-preserving shore symmetry itself stops extending at conductor 16.

More generally every adjacent extension G_(2^(a+1))→G_(2^a) is nonsplit for a≥3, by Theorem 1. The elementary-abelian behavior at conductors 2,4,8 therefore has a structural last level, not an observed finite cutoff. Idoneality fails at16 because classes cease to be determined by their genus signs; each genus now contains two classes.

## 7. Two additional firewalls

**Geometric reflection is not class translation.** On ideals the map I↦−I does nothing, and conjugation induces [I]↦[I]^(−1). At conductors 2,4,8 both therefore act trivially on the class group. The nontrivial translations by A and B in the decorated deletion model are different operations. Moreover z∈O_f gives a principal ideal zO_f, hence the identity class whenever that expression is used. The residue representatives w,j used above need not lie in O_f. One must not relabel the divisor norm-ball involutions z↦−z and z↦conjugate(z) as the two class translations.

**Arithmetic conductor is not electrical conductance.** Nothing in this ideal-extension theorem supplies edge weights, a Laplacian, a boundary current, or a response map. A network representation would need those data and a separately proved intertwining law. Counts, class characters, and the word “conductor” supply none of them.

For the path-reversal metagraph there is a further inherited sidecar: deletion changes path edges into bypass edges. The companion tournament packet supplies an exact bypass-set correction and an explicit pair of lifts whose reversed-and-deleted tournaments differ. The class-group quotient alone does not contain that correction, so the result here does not identify the full metagraph or its electrical response.

## 8. Reproduction, scope, and audit

Run

    python3 -B verify_conductor_tower.py > /tmp/conductor.json
    python3 -B -O verify_conductor_tower.py > /tmp/conductor.optimized.json
    cmp /tmp/conductor.json /tmp/conductor.optimized.json
    cmp /tmp/conductor.json verification.json

There are 28,754 live exact checks. The standard-library checker exhausts finite-ring quotient units, their group tables, all primitive reduced forms, contraction labels, and independent form-extension maps for f=2,4,8,16,32,64,128. It checks the repeated-squaring valuation modulo sufficient powers of two through exponent 64. Ordinary and optimized Python output are byte-identical. The all-a conclusion depends on the written induction, not this finite range.

Hostile controls retain the conductor 16 order-four form and its distinct inverse, the square-kernel class, nonsplit lifts, the two different port stabilizers, and the trivial geometric reflection action at conductor 8. A separate f=40 hostile checks the required dyadic scope of the square-kernel assertion. These checks would fail if genus signs were silently treated as complete class coordinates after conductor 8, if ideal reflection were confused with translation, or if the dyadic genus-coordinate conclusion were broadened to all multiples of 8.

No tournament class census, graph-metagraph isomorphism, electrical identity, native dynamic-source identity, or global idoneal completeness is asserted.
