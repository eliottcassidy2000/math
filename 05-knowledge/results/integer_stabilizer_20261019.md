# Integer collision atoms: moving stabilizers, resonance denominators, and a finite-height atlas

**Status, 2026-10-19:** PROVED elementary extensions, self-checked by
FINITE-EXACT computations. This is a staged research note, not an independently
audited canon promotion. Universal positive-integer convergence remains OPEN.

## 1. Inheritance and the new scope

The closest inherited mechanism is the affine Terras pair chain in
[THM-4569](../../01-canon/theorems/THM-4569-the-terras-clock-recurrence-of-two-collatz-orbits-is-unconditional.md)
and its Haar-almost-sure identity absorption in
[THM-4581](../../01-canon/theorems/THM-4581-haar-coalescence-affinely-related-collatz-orbits-merge-almost-surely.md).
THM-4569 explicitly excludes a Haar-null set of isolated value coincidences
when identifying literal merges with identity-state hits. We retain this
qualification; the present note studies those rational exceptions directly.
The [signed pair-mass computation](collatz_signed_resonance_pair_mass_20261013b.md),
§4, already distinguishes isolated coincidences from affine identities by
retaining odd-step counts. This is inherited, not newly discovered here.

The [clock holonomy calculus](collatz_clock_holonomy_20261007.md),
[source-aware finite bank](collatz_finite_bank_certificate_20261007.md),
and [effective prefix-mass criterion](collatz_effective_prefix_mass_20261005.md)
already retain actual source addresses, positive guards and clock shifts.
No new proof is obtained by merely replacing Haar measure with an atomic
measure: its required pointwise or mass-decay estimate remains unproved.
The companion [relation-geometry note](golden_signed_relation_geometry_20261019.md) supplies
BS(1,3) period quotients and their prime-power multiplicities.

The new extensions are the complete rational stabilizer formula, its
presentation-dependent absorption boundary, the exact denominator bridge to
period quotients, and a quantitative finite-height theorem for isolated
positive integer collisions. The anchor is an actual integer receipt; the
niche is the moving stabilizer; the wildcard is a branch-valid rational
resonance atlas. The hostile is the out-of-phase pair (1,2). The corrected
near miss is requiring group identity where point stabilization suffices.

Throughout, T is the **Terras map**, T(x)=x/2 for even x and (3x+1)/2 for odd x.
For rational 2-adic integers, parity means numerator parity in reduced form
with odd denominator. No odd-modulus reduction of all of Z_2 is assumed.

## 2. Rational points have affine stabilizers

Let Γ=Z[1/3] ⋊ Z act by (k,e):x↦3^k x+e. For a rational point in reduced form

    v=a/(3^s d),   s≥0, d≥1, gcd(d,3)=1,

its stabilizer is exactly

    Stab_Γ(v) = {(k,(1−3^k)v): k∈ord_d(3) Z}.       (2.1)

Use ord_1(3)=1. The same formula applies to rational points of Z_2, where d
is necessarily odd. Every irrational 2-adic point has trivial stabilizer.

**Proof.** Fixing v forces e=(1−3^k)v. After canceling the reduced numerator,
e lies in Z[1/3] exactly when 3^k≡1 mod d, with negative k interpreted using
the inverse of 3. This is precisely the stated order condition. If k≠0
fixes an irrational point, v=e/(1−3^k) would be rational. If k=0, fixing any
point forces e=0. ∎

For v∈Z[1/3], the stabilizer is the conjugate of the full dilation subgroup
by translation by v. In particular **every integer has nontrivial isotropy**,
whereas a Haar-generic 2-adic point has none. If g(v)=u, all affine
presentations of the same pair form the right coset g Stab_Γ(v).

### Literal merging need not ever hit the identity state

Suppose a lawful pair-chain trajectory has u_t=3^{k_t}v_t+e_t. If u_t=v_t,
then e_t=(1−3^{k_t})v_t is even in Z_2. The parity discrepancy is zero, so
k_{t+1}=k_t. The endpoints remain equal, hence this holds at every later
time. Therefore

- if the merge occurs with k_t=0, then e_t=0 and the state is identity;
- if the merge occurs with k_t≠0, identity absorption never occurs afterward.

For the actual integer pair (45,15), the initial translation presentation
(k_0,e_0)=(0,30) reaches

    (u_6,v_6)=(20,20),   (k_6,e_6)=(−1,40/3).       (2.2)

It then remains on a nonidentity moving stabilizer forever. The *same pair*,
initially presented as (k_0,e_0)=(1,0), reaches identity at time 6. Literal
coalescence is presentation invariant; identity-state absorption is not on
these exceptional rational sources. This is compatible with the inherited
Haar theorem because the relevant exceptional sets are Haar null.

## 3. The correct finite target is a point stabilizer

For an odd modulus M coprime to 3, let A=3^k mod M and b=e mod M, reducing
only coefficients from Z[1/3]. The finite point-stabilizer equation is

    (A−1)v+b=0 mod M.                               (3.1)

Put g=gcd(A−1,M). It has exactly g solutions if g divides b, and none
otherwise. This is the elementary linear-congruence theorem, obtained by
dividing through by g and inverting the remaining coefficient modulo M/g.

For M=11^h and k≠0, the exponent in g is

    min(h,0)                         if 5 does not divide k,
    min(h,2+v_11(k/5))               if 5 divides k.  (3.2)

The first line means g=1. For k=0, g=M. Thus the doubled 11-adic valuation
of 3^5−1=242 controls not only orbit multiplicities but the sizes of affine
fixed-point fibres. This is an exact new use of the inherited valuation
formula, not a new proof of that formula.

For an actual integer endpoint v and a cleared coefficient e=a/3^r with
r≥max(0,−k), set

    D=(3^{r+k}−3^r)v+a = 3^r(u−v) ∈ Z.              (3.3)

If D≡0 mod M and |D|<M, then D=0, hence the endpoints literally merge. The
height test is essential: D=M is a congruence false positive. The finite
character projector M^−1 Σ_{j=0}^{M−1} exp(2πijD/M) tests only the congruence;
it gains exact integer equality only with the height sidecar.

For (2.2), r=1, a=40, v=20 gives D=0 at every modulus, while 3^−1≠1 mod 11^h
for every h≥1. Requiring coefficient identity would reject this valid join.
All counts here use finite uniform counting measure. They are not a pushforward
of Haar measure on Z_2 through a nonexistent natural odd-modulus reduction.

## 4. Branchwise collision dichotomy and the period denominator

Fix an integer translation d>0 and equal Terras depth N≥1. On a joint
lawful source cylinder x≡r mod 2^N, let a,b be the odd counts for x+d and x,
and C_u,C_v their affine constants. Then

    T^N(x+d)=[3^a(x+d)+C_u]/2^N,
    T^N(x)  =[3^b x+C_v]/2^N.                       (4.1)

Thus exactly one of these alternatives occurs:

1. a=b and 3^a d+C_u−C_v=0: the whole cylinder merges by affine identity.
2. a=b and that constant is nonzero: there is no collision in the cylinder.
3. a≠b: there is at most one collision, at the rational candidate

       x*=(C_v−C_u−3^a d)/(3^a−3^b).               (4.2)

The candidate is an actual rational 2-adic collision **iff** it has odd
reduced denominator and belongs to the specified source cylinder. It is an
actual positive integer collision only if it also lies in Z_{>0}. Solving
(4.2) without these guards does not certify a Collatz trajectory.

The reduced denominator of (4.2) divides

    3^{min(a,b)} (3^{|a−b|}−1).                     (4.3)

For a fixed general BS(1,3) presentation u=3^{k_0}x+e_0, clearing powers of 3
gives the analogous factor 3^s(3^{|a+k_0−b|}−1) when the exponents differ.
Consequently the universal period-quotient modulus 3^j−1 has a second exact
role: it bounds the non-3 denominator of an isolated rational collision
whose final dilation exponent is ±j. In particular a denominator with an
11 factor requires 5|j, and its 11-adic exponent is at most 2+v_11(j/5).
Numerator cancellation can make this bound strict.

This is an arithmetic bridge to the rational exceptional set, rather than
an analogy based only on equal orbit counts. It does not assert that every
allowable denominator occurs in a lawful or positive collision.

## 5. A finite-height theorem for isolated positive collisions

**Lemma 5.1 (sharp normalized offset maximum).** For every length-N Terras
branch word with a odd steps and affine numerator C,

    0≤C≤2^{N−a}(3^a−2^a),
    0≤H:=C/3^a≤2^N/6.                              (5.1)

The second maximum is attained by the word whose only odd bit is its last
bit. It is a statement about all branch words; no probabilistic assumption
is used.

**Proof.** If the odd positions are j_1<...<j_a, indexed from zero, then

    C=Σ_{i=1}^a 2^{j_i} 3^{a−i}.

Each j_i≤N−a+i−1, giving the first bound by a geometric sum. After division
by 3^a it is 2^N(2^−a−3^−a). For integers a≥1, 2^−a−3^−a has maximum 1/6
at a=1: its successive difference is −2^{−a−1}+2·3^{−a−1}<0 for a≥1.
For a=0, C=0. The indicated word attains equality. ∎

**Theorem 5.2 (isolated positive integer head).** If d>0, N≥1 and x>0 obey
T^N(x+d)=T^N(x), with different odd counts a≠b, then

    x < d/2 + 2^N/12.                               (5.2)

If a>b, the sharper bound is x<2^N/12. The theorem also holds for positive
rational candidates satisfying the displayed branch equation, even before
the 2-adic branch-validity guard is applied.

**Proof.** Write H_u=C_u/3^a, H_v=C_v/3^b. If a>b, put q=3^{a−b}≥3. The
collision equation gives

    x=[H_v−q(d+H_u)]/(q−1) < H_v/(q−1)≤2^N/12.

If a<b, put q=3^{b−a}≥3. Then b>0, hence H_v>0, and

    x=[d+H_u−qH_v]/(q−1)
      < (d+H_u)/(q−1)≤d/2+2^N/12.

All inequalities use the positive d and the plus-one branch constants.
They do not transfer unchanged to the minus-one dynamics. ∎

**The coefficient 1/12 is asymptotically sharp, uniformly in d, even for
lawful positive integers with d≤2^N.** For N=6m, m≥1, put M=2^N and

    x=5(M−1)/9,   d=(17M+10)/18,   x+d=3M/2.        (5.4)

These are positive integers, and d<M. The upper source has parity word
00...01, with its only odd step last; the lower source has word 110...0.
Indeed its first two iterates are (5M−2)/6 and 5M/4, followed by N−2
halvings. Both paths first meet at 5 at time N, with odd counts 1 and 2.
Their exact height is

    x=d/2+M/12−5/6.                                (5.5)

Thus no smaller universal coefficient of 2^N can replace 1/12 in (5.2),
even with a fixed additive constant and the side condition d≤2^N. This
sharpness allows d to grow with N; it does not exclude a stronger theorem
for a single fixed translation d.

**Consequences.** At every fixed depth, all isolated positive integer merges
lie in an explicit finite initial interval. Every merge above that interval
is a whole-cylinder affine identity. If d≤2^N then every isolated integer
candidate lies below 7·2^N/12<2^N, so it must be the canonical nonnegative
representative of its own cylinder. There is no search over arbitrarily
large integer lifts.

Let I_N(d) be the identity cylinders and R_N(d) the guarded positive integer
roots in (4.2). The literal integer merger set at time N is exactly

    (I_N(d) ∩ Z_{>0}) ⊔ R_N(d).                     (5.3)

Both families increase with N, and their unions remain disjoint: after a
nonidentity literal merge its nonzero dilation exponent is frozen (§2).
The first part is open in the 2-adic topology and measured by the inherited
prefix weights. The second part is a finite set of atoms at each depth and
has Haar mass zero, however important it is for an individual integer.

## 6. Unequal clocks isolate the rational phase repairs

For the same translation d but clocks A,B≥0, with odd counts a,b and affine
constants C_u,C_v, a collision solves

    (2^B 3^a−2^A 3^b)x
       =2^A C_v−2^B(3^a d+C_u).                    (6.1)

If A≠B, the coefficient cannot vanish, by unique factorization at 2 and 3.
Therefore every unequal-clock collision of this fixed translated pair is
an isolated rational root on its branch cylinder. The same conclusion holds
for initial slope 3^{k_0}. It is **not** true for arbitrary rational initial
slopes: u=2v on the even branch has a whole-cylinder identity at clocks 1,0.

The actual pair 233=x+10, 223=x has clocks A=15, B=7, odd counts a=10,b=5,
constants C_u=167983, C_v=211, and common endpoint 425. Dividing (6.1) by 2^7
leaves

    −3159 x=−704457,  hence x=223,
    3159=3^5(2^8−3^5)=3^5·13.                       (6.2)

Its complete actual parity words are saved in the exact output. This is a
branch-valid positive integer isolated resonance, consistent with the
intermediate equal-clock join of 334 and 335. The determinant 2^8−3^5=13
is the exact clock arithmetic here; replacing it with the eleven-cell
modulus would lose the eight-step phase shift.

This explains why a Haar-null branch atlas is a natural carrier for integer
clock repairs. It is not evidence that all positive integers occur in that
atlas with a smaller certified partner.

## 7. Focused exact controls

The companion [script](../experiments/integer_stabilizer_20261019.py) performs **279,578
explicit checks**, not removable Python asserts. [Normal output](integer_stabilizer_20261019.json)
matches its optimized-mode replay byte for byte. It checks:

- both presentations of (45,15), the literal nonidentity tail, and 3,125
  rational stabilizer-membership cases;
- finite affine fixed-point fibres at moduli 11,121,1331, including negative
  dilation exponents and the zero-coefficient case;
- every bit word through length 9 for both offset bounds and equality;
- all 2^N source cylinders for N=1,...,12 and d∈{1,2,3,15,30,61};
- each lawful rational root by direct fraction iteration, each identity or
  parallel cylinder at two independent positive lifts, and every positive
  isolated integer root by an independent complete bounded-head scan;
- the strict height boundary, the actual unequal-clock 223/233 example, and
  the exact integer sharpness family (5.4) at N=6,12,...,60.

Selected complete rows (roots are the lower source x):

| N | d | Identity cylinders | Isolated positive integer roots |
|---:|---:|---:|---|
| 6 | 30 | 0 | 2,10,12,14,15 |
| 10 | 1 | 351 | 7 |
| 12 | 1 | 1476 | 7,8,9,15 |
| 12 | 30 | 49 | 2,4,5,10,12,14,15,22,23,34 |
| 12 | 61 | 15 | 7,8,9,24,29,35,45 |

The integer root counts are not added to Haar cylinder numerators. They are
separate source certificates. In particular the pair (8,7) merges at depth
10 but its initial translation presentation never subsequently hits identity.

Reproduction from the repository root:

```sh
python3 05-knowledge/experiments/integer_stabilizer_20261019.py > /tmp/integer-stabilizer.json
python3 -O 05-knowledge/experiments/integer_stabilizer_20261019.py > /tmp/integer-stabilizer-O.json
cmp /tmp/integer-stabilizer.json /tmp/integer-stabilizer-O.json
cmp /tmp/integer-stabilizer.json 05-knowledge/results/integer_stabilizer_20261019.json
```

## 8. What remains for integer Collatz

A valid route can now test a proposed source receipt against three explicit
carriers: an identity cylinder, a guarded isolated equal-clock root, or a
guarded unequal-clock root. Modulus banks can accelerate exact displacement
checks, but must retain actual parity, endpoint height, source positivity,
and both clocks. The partner must already be certified or genuinely smaller
for induction to close.

The missing theorem is still **uniform source coverage** by such receipts.
It is not supplied by almost-sure absorption, the finite branch atlas, the
11-adic fibre counts, or the height bound whose interval grows with N. Even
universal equal-clock coverage for one fixed translation is false: (1,2)
remain out of phase. Allowing both clocks repairs that obstruction but does
not prove a source-wise well-founded decrease.

The next meaningful target is a source-labelled, clock-aware cover with a
well-founded rank on its unresolved sources, or a proved decay estimate for
a full-support atomic source measure. Both must account for the isolated
stabilizer receipts explicitly. The existing finite source-energy frontier
remains an independent candidate for supplying that missing descent.
