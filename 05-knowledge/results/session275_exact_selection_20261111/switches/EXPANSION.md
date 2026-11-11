# Expanding the controller: parity resolution and a positive causal product gate

**PROVED from the physical fusion carrier and explicit certified inputs.**
These are additional constructive ingredients, not Chan–Pak coverage.

## 1. A mixed one-way bank removes the fixed modulus4

The pure5/9 controller in `PROOF.md` is trapped in one residue modulo4 at
fixed length. Rather than extrapolate from that obstruction, add a second
physical letter bank:

    N=(26,32,3,6,14,29), P=(52,41,2,5,14,28).
    Common (H,h,Q,k,V,e)=(37,9,84,0,43,7).
    (U,d)_N=(59,5), (U,d)_P=(62,4).

Direct subset counting certifies the input data. Thus

    M(P)-M(N)=3E12-E13.

It obeys the same spacer-uniform square-zero law as the5/9 pair. Write C for
the rooted directed triangle. For m>=1 consider the actual words

    T=(A or B)^m fused (N or P) fused C.

All m+1 binary choices remain EXACTLY independent, not approximately so.
The neutral last switch contributes6*5^m: its pairing with C's last column is
3(3-1)-(1-1)=6, multiplied by the preceding root-deletion factor5^m.

For the other slots define

    f_0=1244,
    f_r=(72966*64^r-396*5^r)/59  (r>=1).

Their weights, in order, are5^i f_(m-1-i), i=0,...,m-1.
Indeed the fixed suffix N fused C has(h,V,e)=(9,211,7), so its last-column
entries relevant to the trade are(202,-2,9). Pair this column with
D M(B)^r to obtain the displayed expression; at r=0 retain the parity term.

**Theorem E1.** At every length m>=1 the gcd of all m+1 response weights is
exactly2. Thus this reusable independent controller has no growing-prime
obstruction and no fixed congruence obstruction beyond odd parity.

Proof. All weights are even. The final two weights are1244*5^(m-1) and
30*5^(m-1), whose gcd is2*5^(m-1). For r>=1, f_r is nonzero modulo5 by the
closed formula; f_0=1244 is too. The first weight therefore removes every
factor5 from the common gcd. QED.

At m=1 the baseline is1795 and weights(1244,30); hence the four counts are
1795,1825,3039,3069. At m=2 baseline75055 and weights(79116,6220,150).
These are counts of literal tournaments, with13 vertices for m=1.
The modulus repair does NOT fill the gap between1795 and1825. In fact the
weights are again superincreasing: f_1>60f_0,
f_r=64f_(r-1)+396*5^(r-1) for r>=2, and1244/30>41.
Every successive descending-weight ratio therefore exceeds12.
So this exact selector remains sparse.
This is a useful way to distinguish two repairs: removing an arithmetic
lattice obstruction versus supplying enough overlapping physical digits.

## 2. Converse turns a one-way trade into a controlled nonlinear gate

Converse reverses fusion order and interchanges U,V and d,e while keeping
H,h,Q,k. The converse of the5/9 trade therefore has difference

    D_right=6E24+2E34+4E14,

whereas D_left=6E12+2E13+4E14. For every fusion spacer S,

    D_left M(S) D_right=(36Q_S+4k_S)E14,
    D_right M(S) D_left=0.

Two left trades or two right trades also give zero. Every product containing
three trade differences vanishes, with arbitrary baseline spacers between
those differences. These claims follow directly from triangular multiplication,
not from treating signed matrix differences as physical inverses.

**Theorem E2 (causal quadratic compiler).** An arbitrary word whose variable
slots use the forward5/9 bank or its converse has a Hamilton-count polynomial
of degree at most2 in the independent choice bits. The ONLY nonzero pair
interactions are a forward trade to the left of a converse trade. If S is the
baseline word between those slots and L,R the exterior baseline words, its
coefficient is

    h_L h_R (36Q_S+4k_S)>0.

Every third and higher discrete derivative is zero. All coefficients and the
absolute baseline are computable in fixed4-by-4 state. Consequently the
quadratic polynomial is an exact target-test/constructor for this physical
family, not a perturbative approximation.

Proof. Expand the full matrix product. A two-trade term is the displayed
central matrix times the prefix/suffix h diagonals. A third trade annihilates
it since every trade has zero first column and zero last row. This gives every
coefficient and eliminates every higher term. A nonempty physical spacer has
Q>0 and k=0; an empty spacer is the singleton identity with Q=k=1. Thus all
allowed interactions are positive. QED.

Literal initial gates, ordered as(bit_left,bit_right)=(00,01,10,11):

    empty spacer:       (75,115,115,195), interaction40;
    root-source edge:   (125,225,205,377), interaction72;
    rooted triangle:    (275,475,475,819), interaction144.

For example the edge-spaced tournament implements exactly

    H=125+80x+100y+72xy.

The lower-order terms are paid costs, not deleted or subtracted tournament
counts. Reversing the left/right control types removes the mixed term; it does
not merely change its sign. A periodic arrangement therefore has a precisely
known causal interaction graph rather than an uncontrolled accumulation of
root debts.

## 3. A physical multiplication of independently selected root responses

Let L(x) and R(y) be arbitrary forward one-way-controller words, with baselines
L(0),R(0). Write u(x)=U_L(x)-U_L(0), v(y)=U_R(y)-U_R(0).
Insert any fixed NONTRIVIAL pointed S and form

    F(x,y)=L(x) fused S fused converse(R(y)).

Because the two input differences are confined to a row and a column, the
mixed difference is exactly

    H(F(x,y))-H(F(x,0))-H(F(0,y))+H(F(0,0))
      = Q_S u(x)v(y).                                  (2)

The parity contribution disappears since k_S=0. Each u,v is a linear form in
its own chosen bits. For pure5/9 words,

    u(x)=6 sum_(i=0)^(m-1) 5^i64^(m-1-i)x_i.

Thus (2) is a literal positive multiplication gate on two independently
encoded root-response integers, with exact known linear and baseline costs.
The interface remains in the same fusion carrier, so gates may be explored
without forgetting the state that made them lawful. This does NOT allow an
arbitrary integer as u or v: the sparse digit languages and their linear costs
must still be solved in an inverse target argument. Nor does nesting the
formula automatically implement arbitrary higher-degree circuits inside a
single fixed root; Theorem E2 states the degree bound for this family.

## 4. What the expansion learned

The randomly sampled historical `berggren-positive-cubes-and-amm-causal-horizons-kps-s181.md`
asked whether a frozen coordinate is the resource that must move. Here that
question first freed scalar H, then freed the direction of information flow.
The result is an exact distinction between additive digits, parity repair,
and positive products with paid costs. This is a mechanism-level analogy,
not a transfer of a Pell or AMM theorem.

The S216 eleven-offset substitution bank remains a different operation with
an unresolved two-cover sidecar. The new product gate does not silently fix
that incompatibility. A concrete next target is an interval-overlap theorem
for the causal quadratic form, or a hybrid where the11 offsets pay its linear
terms while preserving its controllable root responses. Exact formulas now
make that a well-posed physical question; a proof of coverage is still absent.

Reproduce with `python3 -B 05-knowledge/results/session275_exact_selection_20261111/switches/expand.py`.

## 5. Concurrent energy bridge: exactly the same6^m can mean very sparse outputs

At incoming main057f7b27c88c, S93's
`01-canon/theorems/THM-5141-one-apex-counts-three-quarters-exponent-reversal-graph-energy.md`
is explicitly **DRAFT / audit pending**. Its signed-digit reversal-graph energy
6^w is used here as a research prompt, not as a proved dependency; none of its
counting-exponent conclusions is imported.

Our pure affine controller permits an elementary independent statement.
Let S_m be its2^m scalar outputs. Every descending weight exceeds twice the
sum of all smaller weights. Hence no nonzero coefficient vector in
{-2,-1,0,1,2}^m has zero weighted sum: the leading nonzero term dominates the
rest. Thus pair differences are classified exactly by their signed digit
vector eps in{-1,0,1}^m, and a vector with z zero coordinates occurs2^z times.
Consequently, for every m,

    |S_m-S_m|=|S_m+S_m|=3^m,
    E(S_m)=#{a+b=c+d in S_m^4}=sum_eps4^(#zeros eps)=6^m.

This is an exact source-target map from Boolean choices to physical scalar
outputs; it preserves the entire signed-digit difference multiplicity, not
just its total energy. The proof is the strong separation of our response
weights, so it neither assumes nor verifies the draft reversal-graph theorem.

Yet these S_m lie in superincreasing gaps and in one residue modulo4. The
energy identity does not imply an anchored interval, much less Chan–Pak.
Conversely a consecutive lattice interval has many carrying identities and
higher additive energy. This sharpens the next research choice: the one-way
controller eliminates every mixed carry and is wonderfully invertible, but
exact filling needs a controlled mechanism that allows neighboring digit
blocks to overlap. The causal product gate supplies one explicitly priced
interaction to investigate; merely optimizing the number of distinct outputs
or reproducing6^m does not replace the overlap obligation.
