# Two-sided numerical-semigroup saturation and the central residue bottleneck

**Status: PROVED, scoped; physical transfer conditional on the stated actual banks.**
This is an exact bridge to numerical semigroups and Apéry normal forms, not a
Chan–Pak proof. The arithmetic is elementary classical coin algebra, not claimed
as a new numerical-semigroup theorem. The contribution is a complete boundary-hole
model and a local saturation obstruction on the *physical central-fusion carrier*.

## Inheritance and map

Session282 [central convolution](../../session282_chanpak_20261112/central/PROOF.md)
already proves both the exact sum law and an interior conductor bound. We use,
not rediscover, those. Its inherited41/43 h15 gate, new2747/2749 h169 gate and
missing474055 are the canonical hostile. Session285's
[compatible fibers](../../session285_five_web_20261112/tournaments/PROOF.md)
supplies another actual gate14207/14209 with h3357, but no complete digit bank.

A deliberately distant source is
[THM-3703, exceptional-quartic SAGBI and Apéry module](../../../../01-canon/theorems/THM-3703-russell-cylinder-exceptional-quartic-sagbi-module.md).
There a polynomial restriction algebra admits a finite residue module and a
finite conductor-hole list. Here the analogous precise object is the set of
*normalized central count differences*. Preserved: integral additive membership
and an explicit finite boundary obstruction. Lost: absolute count baseline,
physical digit availability, and the continuing carrier. We retain all three.
No Jacobian or algebraic-normalization theorem transfers to tournaments by analogy.

For actual full-central banks indexed i, put counts b_i+2E_i, common deletion
h_i>1, pairwise coprime h_i, D=product h_i and w_i=D/h_i. Use no fixed nonunit
deleting context. The physical count set is

    H0 + 2S,          S = w_1 E_1 + ... + w_r E_r.              (1)

H0 is the literal all-zero baseline. Formula(1) and its preservation of the
remaining fusion packet are inherited. Banks with only a shared receiving
interface are NOT automatically full-central and cannot be substituted here.

## 1. A local bottleneck stronger than primitive gcd

Put R_i=w_i E_i mod h_i. If R_i is a proper residue set, let ell_i be the
longest cyclic run of consecutive residues lying wholly in R_i. Then EVERY
consecutive integer run in S has length at most ell_i. Consequently every
odd-spaced interval in(1) has at most min_i ell_i counts, the minimum over
proper R_i. In particular it has at most |E_i mod h_i| for each unsaturated i.

Proof: modulo h_i all terms except w_i E_i vanish. Project a consecutive run;
a run of length h_i would already see every residue. Shorter runs must lie
in a cyclic run of R_i. Multiplication by w_i permutes residues, so its set
cardinality equals |E_i mod h_i|. QED.

Thus an independent two-choice central gate with distinct digits mod h_i>2
prevents *any* widening beyond two consecutive odd counts, however large or
numerous the other coprime banks become. If h_i>=3 and its two projected
residues are not cyclic neighbors, it prevents even a consecutive pair.
For h15/h169 and E_i={0,1}, R_15={0,4}, so that fixed fusion cannot contain
any adjacent odd counts. This strengthens, but does not replace, the known
single missed474055. Changing the baseline family, sharing a core, noncentral
surgery, or rerooting leaves this theorem's universe and is not excluded.

The exact residue mask is more informative than coefficient gcd. A gcd1 mesh
is necessary for a long integer interval; it does not imply even one unit gap.

## 2. The completed-bank theorem: two finite boundary memories suffice

Assume now E_i=[0,m_i] with m_i>=h_i-1. Define

    Gamma = <w_1,...,w_r>,       Smax=sum_i w_i m_i.

Then the ENTIRE bounded support is exactly

    S = { z in [0,Smax] : z in Gamma and Smax-z in Gamma }.    (2)

This gives positive constructive membership, not only an obstruction.
For any integer z let r_i be its unique digit in[0,h_i-1] satisfying
w_i r_i=z mod h_i. Put R(z)=sum_i w_i r_i. Then

    z in Gamma iff z>=R(z).                                  (3)

Indeed every nonnegative coefficient a_i is uniquely r_i+h_i k_i with
k_i>=0; the excess is D sum_i k_i. Conversely a nonnegative excess multiple
of D can be assigned to any one generator.

For the bounded box write c=(z-R(z))/D and
q_i=floor((m_i-r_i)/h_i)>=0. Session282's exact carry law is

    0<=c<=sum_i q_i.                                         (4)

The canonical residue digits of Smax-z are
r'_i=m_i-r_i-h_i q_i. Hence R(Smax-z)=Smax-R(z)-D sum q_i,
and(3) for Smax-z gives exactly the upper inequality in(4). This proves(2).
A greedy allocation of c among capacities q_i returns actual digit choices
and hence the actual tournament whenever the hypothesized banks exist.

No enumeration through Smax is required: the fixed residue table and its
reflection classify every target. This is the precise limited Apéry bridge.

## 3. Exact holes, not just a guaranteed middle

For r>=2 put

    F=(r-1)D-sum_i w_i,         C=F+1.

The nonnegative gaps of Gamma end exactly at F. This follows from(3): the
largest canonical R is sum_i w_i(h_i-1)=rD-sum_i w_i=F+D, and R-D is the
largest missing integer in that residue. Thus C is the inherited Session282
conductor, now with its sharpness and exact boundary behavior explicit.

Moreover Gamma is symmetric: for 0<=z<=F exactly one of z,F-z lies in Gamma.
Write z=R+cD. The canonical digits for F-z are h_i-1-r_i, because
F≡-w_i mod h_i; their sum is F+D-R. Consequently its carry is -c-1.
Exactly one of c and -c-1 is nonnegative. There are therefore(C/2) gaps.
(C is even because the h_i are odd and pairwise coprime.)

If Smax>=2F+1, the low gaps and reflected high gaps are disjoint. Formula(2)
then gives EXACTLY C missing integers in[0,Smax], independent of how much
the completed banks are enlarged. They are

    Gaps(Gamma) union {Smax-g:g in Gaps(Gamma)}.               (5)

All integers in[C,Smax-C] occur, and (2) also decides every boundary point.
The largest guaranteed central odd interval has Smax-2C+1 values when
nonempty. This counts outputs of this family, not omitted tournament counts.
For r=1 the coefficient w1=1 gives S=[0,m1] with no holes, separately.

**Arithmetic illustration, not claimed physical bank supply:** take
h=(3,5), m=(8,8). The weights are5,3, Smax64, F7 and gaps{1,2,4,7}.
The bounded sum has exactly the eight missing values
{1,2,4,7,57,60,62,63}; all8..56 occur. A hypothetical actual pair of
central nine-banks with these deletion counts would physically realize that
support after applying its own baseline. Existing h3 fibers do not supply
such a bank; Session282 proves a stronger obstruction there. The example
illustrates the compiler, not a construction claim.

## 4. Why this does not silently create an infinite fixed-interface semigroup

For a nonempty deletion D0 of n vertices, every induced subtournament has
at least one Hamiltonian path, including H(empty)=1. Therefore

    Q(D0)=sum_A H(D0[A])H(D0\A) >= 2^n.

Fixing a full packet's h,Q bounds n<=floor(log2 Q). There are only finitely
many such labelled rooted tournaments up to relabelling, and their H values
are bounded by(n+1)!. Thus an infinite central bank at a single fixed packet
cannot exist. A physical normality/conductor program must carry a sequence
of expanding packets and prove the needed digit supply at each scale.
Neither a formal semigroup closure nor signed subtraction in a SAGBI algebra
provides those operations. The relevant missing theorem is a *scale-compatible
physical saturation theorem*, not another gcd computation.

## 5. A shared-core escape is already concrete, but its scale matters

The companion [automata lane](../automata/PROOF.md) finds the actual shared-core
language(A|BE*C)*. Its normalized gaps have a truncated<8,17> description,
with56 persistent holes after sufficient length. This escapes independent-bank
residue projection because the choices are correlated excursions inside one
receiver fiber. It realizes the finite-boundary-memory idea physically.
Its outer mesh2*13^(k-2) remains, however: conductor saturation of normalized
gaps is not absolute odd-target saturation. The reset lane investigates that
separate physical problem. No full Chan–Pak conclusion is drawn here.

## Research outcomes and stopped routes

A structured receiver scan over repeated binary patterns was tried as a cheap
probe for renewable unbounded intervals, but stopped without a theorem or new
witness; no claim of exhaustive absence is made. The productive route was to
replace blind width enlargement with the exact projected-residue condition,
and then compare independent fusion with correlated excursions. The
normalization analogy contributes a finite exact target decoder, while its
physical supply/scale limitations stay visible.

## 6. Mixed moduli do not yet supply a saturated bank: an all-length theorem

A natural positive attempt is to mix the literal five-vertex deletion cores

    D13=(18,8,3,5,14),    D15=(24,17,3,6,12).

All32 root attachments of D13 have(h,Q)=(13,104), and all32 attachments
of D15 have(h,Q)=(15,120). Thus their common ratioQ/h=8 allows correlated
receiver carries, while the deletion types are coprime. This is more general
than fusing already-central switches independently. The entire64-letter
physical alphabet has40 distinct fusion packets; all concatenations of these
packets are actual rooted tournaments, in any order.

**Theorem.** No full-central fiber in this entire mixed-word language contains
a complete set of normalized H-difference residues modulo its deletion h.
This is an all-word-length statement, not a finite lack-of-witness claim.

For a word with a letters of deletion13 and b of deletion15, D=13^a15^b.
Expanding its upper-triangular fusion matrices shows that the scalar is a sum
of three types of terms:

1. D and (H_i-h_i)D/h_i;
2. (U_i-h_i)(V_j-h_j)8^(j-i-1)D/(h_i h_j), i<j;
3. (d_i-h_i)(e_(i+1)-h_(i+1))D/(h_i h_(i+1)).

There are no terms with three marked letters: a path through the4x4 matrix
can visit at most ONE of the two available interior channels, hence at most
two marked letters; channel3 has
zero diagonal so only adjacent letters contribute to type3. This is just
the exact matrix expansion, not a positivity argument. Every term is divisible by

    G=13^max(a-2,0) * 15^max(b-2,0).                       (6)

For lengthk=a+b>=5 at least one exponent is positive, soG>1. It divides D
and every scalar H in every word of that fixed deletion multiset. SinceG is
odd, all normalized differences(H-H0)/2 are also divisible byG and cannot
exhaust residues modD. A fixed full packet determinesa,b fromD, so combining
words of different length or order cannot evade this argument inside a fiber.

For k<=4, the explicit exact packet DP supplies the remaining finite base:

|k|deletion D : maximum number of normalized residues in any central fiber|
|--:|--|
|1|13:1; 15:2|
|2|169:2; 195:2; 225:3|
|3|2197:4; 2535:4; 2925:3; 3375:4|
|4|28561:8; 32955:8; 38025:6; 43875:4; 50625:5|

Every value is strictly smaller than D. At the four levels there are
40,1584,62416,2448564 distinct full packets, respectively. Deduplication is
safe because all future operations in this specified language are functions
of that exact packet and no ordering/terminal constraint is imposed afterward.
The verifier reconstructs every letter from literal subset-path counts and
performs every packet transition. Thus the finite base plus(6) proves the
all-length theorem. Empty words are a trivial single count and not a completed
nonunit bank.

This failed positive attempt is informative: a fixed finite set of prime
contents, integralQ/h and two-interior-channel fusion make content accumulate
faster than correlations can remove it. To make the semigroup compiler useful,
one must change that carrier, inject genuinely new deletion primes, leave
integral-ratio backgrounds, or use a physical reset. Coprime deletion types
alone are not enough. The theorem does not exclude any of those escapes.

### 6.1 A general finite-alphabet deletion threshold, with a useful escape clause

The content mechanism is not special to13 and15. Fix any finite alphabet of
nontrivial rooted tournaments for which Q_i/h_i is an integer and K_i=0.
Permit unit-deletion letters h_i=1. LetP be the finite set of primes dividing
any h_i, let M_p=max_i v_p(h_i), and put B=product_(p in P) p^(2M_p).
For ANY word of deletion D, the scalar count is divisible by

    G_*(D)=product_(p in P) p^max(v_p(D)-2M_p,0).        (7)

The same three-type matrix expansion proves(7), because each summand removes
at most two endpoint deletion factors; all intervening Q/h ratios are integers.
The divisor depends only on D and the alphabet, so it holds even when a full
fiber combines different words, orders, lengths or factorizations of D.
Consequently, onceD>B, G_*(D)>1 and no full-central fiber at that deletion can
saturate its normalized residues. This is a finite-alphabet theorem; no bounded
word length or exhaustive search is used. IfP is empty, D=1 always and the
conditionD>B is vacuous, as it must be.

For the two-core alphabet, B=38025. The sharper count-sensitive(6) and finite
base yield the stronger everywhere-negative result above. The general theorem
also highlights a *nonexcluded constructive route*: add arbitrarily many
unit-deletion letters while keeping D small, so Q and interface geometry can
grow without triggering(7). It does not assert that those growing fibers
actually become complete. Session282's receiver banks pursue precisely a
related low-deletion, growing-Q direction; the distinction between shared
receiver and full-central equality remains decisive.

## 7. A positive low-deletion escape: six full-central values from unit padding

The unit-letter loophole in§6.1 is useful even before full saturation. It
produces SIX actual fourteen-vertex rooted tournaments with common FULL
noncentral packet

    (h,Q,U,V,d,e)=(15,30720,24834,8259,0,15)

and scalar counts

    6689,6691,6695,6697,6701,6703.                         (8)

This is three adjacent odd pairs, NOT six consecutive odds. Its normalized
central digits are{0,1,3,4,6,7}; it does not fill all15 residues.

For a reproducible construction use deletionC5=(6,12,24,17,3), rootincoming
core cutc, prefixtransitive lengtha with incoming binary maskp, and suffix
length8-a with incoming maskq. Fuse these three at their roots, prefix before
core before suffix. Explicit(a,p,c,q) for the counts in(8), in order, are

    (7,12,7,0), (7,12,11,0),
    (3,4,7,8),  (3,4,11,8),
    (1,0,3,34), (1,0,5,34).

The standalone verifier constructs all graph rows, counts literal induced
subset paths and all seven coordinates, and compares with the factorized
fusion calculation. All six noncentral packets coincide, so every common
unit-deletion context on EITHER side preserves their six offsets forever.
Arbitrary contexts multiply them by the product of deletion counts. This is
stronger continuation compatibility than a shared-receiver bank, but less
scalar density than the existing renewable nine-consecutive receiver bank.
No optimality or historical priority beyond a scoped exact-number search is
asserted.

### Why the apparent arithmetic pump stops after one step

LetX be the rootedC5 with incoming cut7. LetA andB be unit-deletion letters
whose incoming length4 masks are1 and8 respectively. Their packets are

    X=(41,15,120,57,78,9,6),
    A=(1,1,16,8,2,0,0),
    B=(9,1,16,15,9,1,1),

written(H,h,Q,U,V,d,e). The exact matrix difference XB-AX changesH,d,e
by+6 and changes no other coordinate. Discardingd,e falsely predicts that
one can slideX along repeated letters and accumulate arbitrary multiples6.

Choose fixed unit contextsL of incoming mask4,length3 andR of mask0,length1.
For every k>=1 consider L A^j X B^(k-j) R, 0<=j<=k. All have a common
noncentral packet, but their scalar differences relative toj=0 are exactly

    0, -6, -6, ..., -6.                                 (9)

Proof: sandwiching the raw difference by unit contextsC,E changes its scalar
by6(d_C+e_E-1). Heree_(B^rR)=1 always, whiled_L=1 andd_(LA^j)=0 forj>=1.
Thus the first slide has size6 and every later slide has size0. The outer
nonempty unit contexts kill thed,e differences; U,V never differed. This
proves(9) for everyk, not just the displayed finite replay. It explains why
the promising six-set supplies more central digits but not an automatic
infinite numerical-semigroup generator.
