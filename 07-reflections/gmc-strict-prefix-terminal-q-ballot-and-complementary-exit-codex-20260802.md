# Strict-prefix product-Gamma transport: ballot roots and a complementary exit

Date: 2026-08-02

Status: **PROVED STRUCTURAL REDUCTION + EXACT COUNTEREXAMPLES + FINITE-EXACT
SCOUT / OPEN ALL-ORDER CONJECTURE.** This note does not prove arbitrary
strict-prefix Hankel sign regularity, product-Gamma width-three goodness,
GMC(2), NC2, LRC(14), or any Jacobian statement.

## 1. Inheritance and audit

The closest proved mechanism on current origin/main is
THM-3056-product-gamma-reciprocal-hypergeometric-and-hankel-reversal.md.
For a pure reciprocal product

~~~text
K_n=product_(ell=1)^Q (beta_ell)_n^(-1),
~~~

positive row clearing turns reversed columns into nested prefix products.
Their coefficient matrix is the path matrix of ordered nonnegative Jacobi
layers, so Cauchy--Binet gives every generalized-Hankel sign.

The canonical near miss is
THM-3065-reciprocal-beta-gap-contiguous-hankel-wall.md. Nonpositive
cumulative exponent prefixes still prove the arbitrary-gap order-two sign,
but a zero prefix admits an exact wrong-sign order-three determinant. Its two
nearest strict repairs survived the finite hostile banks.

The current origin/main file
THM-3079-laguerre-pf-row-transform-and-strict-integer-mesh-terminal-minus-one.md
is still a RESERVED / UNPROVED EMPTY STUB. A separately audited local
candidate at commit cb1d97c512 supplies a Newton--PF proof for terminal
prefix -1; until that promotion reaches the shared truth surface it is
inheritance signal, not a proved dependency here. The present note attacks
the explicitly open terminal mass Q>=2.

The incoming arbitrary-anchored bank was replayed before making the new
reduction.

- Both THM-3110 exact companions run identically under ordinary and optimized
  Python and byte-match their stored outputs.
- The THM-3112 companion does the same.
- gmc_product_gamma_order5_transport_circuit.py does the same. Its exact
  minimal fifth-order circuit and its raw 2 x 2 coefficient sign conflict
  are sound.

That circuit is not the matrix studied below. It forbids an ordinary
coefficient-TP proof for the signed THM-3110 response bank. Here the proposed
matrix is the coefficient array of one strict-prefix Gamma quotient after
integer-mesh row clearing. Confusing the two would turn a useful no-go into a
false global obstruction.

Two proved results arrived on origin/main while the scout was running and
sharpen the connection. THM-3115 turns each normalized low-degree signed
fibre coefficient vector into the boundary of a nonnegative one-chain on the
partition-refinement Hasse diagram; the zero-total-mass face, not the ambient
coefficient cone, is what makes the transport positive. THM-3117 proves every
rank-four macro current has a signed forest-cycle lift while every literal
same-sign lift is impossible. Neither theorem implies the Gamma statement
below. They independently support its typing: retain the inherited face and
allow an internal signed chain before asking for a positive quotient.

THM-3119 then adds a second precise deletion lesson. Raw unweighted block
deletion fails, while the unique factorial diagonal gauge turns it into
genuine size-biased same-label deletion and preserves the rescaled Young
carrier order. This suggests testing diagonal gauges on coefficient levels
for the restricted inverse Jacobi map in Section 6. No such conjugacy is
claimed here; the current positive object remains the inherited divisibility
face.

## 2. Strict prefixes are a literal strand inventory

Let

~~~text
0<alpha_0<alpha_1<...<alpha_N,
alpha_(j+1)-alpha_j in Z_(>0),

H_n=product_(j=0)^N (alpha_j)_n^(e_j),
E_j=sum_(i=0)^j e_i<0.                                  (1)
~~~

Put

~~~text
S_j=-E_j>0,                 Q=S_N=-sum_j e_j.            (2)
~~~

Summation by parts gives

~~~text
H_n=(alpha_N)_n^(-S_N)
    product_(j=0)^(N-1)
      [(alpha_(j+1))_n/(alpha_j)_n]^(S_j).              (3)
~~~

Peel one unit from every strict height. With T_j=S_j-1>=0,

~~~text
H_n=1/[(alpha_0)_n (alpha_N)_n^(Q-1)]
    product_(j=0)^(N-1)
      [(alpha_(j+1))_n/(alpha_j)_n]^(T_j).              (4)
~~~

For the integer gap d_j=alpha_(j+1)-alpha_j,

~~~text
(alpha_(j+1))_n/(alpha_j)_n
 = (alpha_j+n)_(d_j)/(alpha_j)_(d_j).                  (5)
~~~

Consequently, up to a positive constant,

~~~text
H_n=P(n)/[(alpha_0)_n (alpha_N)_n^(Q-1)],

P(n)=product_(j<N)(alpha_j+n)_(d_j)^(S_j-1).           (6)
~~~

This is not merely a convenient algebraic normal form. Read S_j as the
number of active strands above the shape interval
[alpha_j,alpha_(j+1)). One strand crosses every interval and becomes the
reciprocal-Gamma ray at alpha_0; Q-1 terminal strands are recorded as rays at
alpha_N; the remaining S_j-1 strands crossing edge j are the finite
rising-factorial blocks in P.

A birth/death version is equivalent. Since

~~~text
e_0=-S_0,                 e_j=S_(j-1)-S_j,             (7)
~~~

the signed event inventory -e_j has strictly positive cumulative sums S_j.
Stack or queue matching pairs every death with an earlier birth and leaves
exactly Q unmatched rays. Different matchings change the strand sidecar, not
the product in (6).

Equations (3)--(7) are proved identities for every inventory in (1).

## 3. The row-cleared object and its exact ballot law

Fix a column ceiling C. Positive row scaling turns the column with 0<=q<=C
into the polynomial

~~~text
F_q(u)=P(u+q)
       (u+alpha_0+q)_(C-q)
       (u+alpha_N+q)_(C-q)^(Q-1).                     (8)
~~~

Indeed,

~~~text
[(alpha_0)_(u+C)(alpha_N)_(u+C)^(Q-1)] H_(u+q)
 = positive_constant * F_q(u).                        (9)
~~~

All roots of F_q(-u) are positive. More importantly, comparison of adjacent
reversed columns is exact at the original shapes:

~~~text
mult_(alpha_j+q-1; roots(F_(q-1)))
 -mult_(alpha_j+q-1; roots(F_q))
 =-e_j.                                               (10)
~~~

Every other root multiplicity is unchanged. Therefore the cumulative root
surplus through alpha_j is

~~~text
sum_(i=0)^j (-e_i)=S_j>0,                              (11)
~~~

and the total surplus is Q.

The ordinary one-dimensional ballot lemma now gives a monotone matching:
every root parameter of F_q can be injected into a no-larger root parameter
of F_(q-1), with Q roots left over. This is the exact Hall certificate which
an undifferentiated root-majorization slogan was missing.

The same prefix identity explains the already proved order-two survivor. For

~~~text
R(x)=H_(x+1)/H_x=product_j(x+alpha_j)^(e_j),           (12)
~~~

summation by parts gives

~~~text
d/dx log R(x)
 =E_N/(x+alpha_N)
  +sum_(j<N) E_j
    [1/(x+alpha_j)-1/(x+alpha_(j+1))] <0.             (13)
~~~

Thus R is strictly decreasing and every arbitrary-gap order-two Hankel minor
is negative. Integer mesh is unnecessary for (13); it becomes load-bearing
at order three.

## 4. Root matching proves coefficient TP2, but stops at order three

The ballot matching has one exact coefficient consequence. Suppose

~~~text
p(u)=product_i(u+a_i),
q(u)=product_i(u+b_i) product_j(u+c_j),
0<b_i<=a_i,                 c_j>0.                    (14a)
~~~

Then the two coefficient columns of p and q are totally nonnegative. It is
enough to inspect the two elementary root moves. Write h_k=[u^k]h(u), with
h_(-1)=0.

- Replacing one root a by b<=a sends p=(u+a)h to q=(u+b)h, and for k<l

  ~~~text
  p_k q_l-p_l q_k
   =(a-b)(h_k h_(l-1)-h_l h_(k-1))>=0.              (14b)
  ~~~

- Adding one root b sends p=h to q=(u+b)h, and the same minor is

  ~~~text
  h_k h_(l-1)-h_l h_(k-1)>=0.                         (14c)
  ~~~

The coefficient sequence of a negative-root polynomial is PF-infinity: each
linear factor contributes a nonnegative bidiagonal Toeplitz layer. In
particular it is log-concave, so h_(k-1)/h_k is nondecreasing and the last
expressions are nonnegative. Each elementary move therefore makes q_k/p_k
nondecreasing in k. Ratios multiply under composition, proving (14a).
Zero padding above deg p gives strict positive minors against every newly
created top degree.

Apply this to the ballot matching from F_q to F_(q-1). Nonadjacent columns are
compositions of the same moves. Hence the complete monomial coefficient array
A^M in (21) is **proved TN2 for every strict integer mesh and every ceiling**.
This is stronger than merely recovering the already-known order-two Hankel
sign: it locates the exact first open coefficient order at three.

That stopping point is sharp. Consider the three negative-root polynomials
whose positive root parameters are

~~~text
X_0=(1/2,17/10),
X_1=(1/2,17/10,39/10,47/10),
X_2=(1/10,1,12/5,3,3,41/10).                         (14d)
~~~

There is a monotone injection from X_0 to X_1, and one from X_1 to X_2,
always matching an old root to a no-larger new root. Nevertheless, for

~~~text
p_i(u)=product_(x in X_i)(u+x),
~~~

the coefficient rows 1,u,u^2 have determinant

~~~text
det [[u^k]p_i]_(k,i=0)^2
 =-15565881/125000 <0.                                 (15)
~~~

This is not only a coefficient-coordinate hostile. At the increasing
nonnegative nodes u=(0,1/100,1/50), the evaluation determinant is

~~~text
-7111716618502761/31250000000000000000 <0.            (16)
~~~

Hence root matching, even with strictly lower added roots and exact rational
data, does not imply order-three total positivity. The first failed
implication is

~~~text
adjacent monotone root injection
  does not imply a totally nonnegative coefficient or evaluation flag.    (17)
~~~

The strongest survivor is (10)--(11): in the Gamma problem the root moves
occur at integer translates of one fixed birth/death inventory. That
temporal lattice coherence is absent from (14d).

## 5. A concrete production matrix

Reverse the ceiling and put

~~~text
f_t(u)=F_(C-t)(u),                 0<=t<=C.             (18)
~~~

Then the entire column family is the orbit of one fixed operator:

~~~text
f_(t+1)(u)=R_C(u) f_t(u-1),

R_C(u)=(u+alpha_0+C-1)
       (u+alpha_N+C-1)^(Q-1).                          (19)
~~~

The finite strands appear only in the seed f_0(u)=P(u+C). This is the
requested actual production map, not a tournament or analogy:

~~~text
L_C: f(u) |-> R_C(u) f(u-1).                           (20)
~~~

There are two natural coefficient arrays:

~~~text
A^M_(k,t)=[u^k] f_t(u),

f_t(u)=sum_k A^N_(k,t) binom(u,k).                     (21)
~~~

Total nonnegativity of either array gives the weak checkerboard inequality.
The monomial array already gives strictness as well.

- For A^M, the evaluation matrix [u_i^k] is generalized-Vandermonde totally
  nonnegative on increasing nonnegative nodes.
- If selected columns have degrees d_1<...<d_r, the coefficient minor on
  degree rows 0,d_2,...,d_r is triangular with determinant equal to the
  positive constant coefficient of the first column. The matching generalized
  Vandermonde minor is strictly positive, even when the first node is zero.
- For A^N, [binom(u_i,k)] is the Pascal planar-network matrix on increasing
  nonnegative integer nodes. Its total nonnegativity gives the same weak sign;
  a separate nonzero low-degree term is still needed if one uses only this
  basis for strictness.
- Equation (9) supplies only positive row scalings, and reversing the selected
  columns contributes (-1)^(r choose 2).
- Cauchy--Binet writes every reversed generalized-Hankel minor as a sum of
  products of a Vandermonde/Pascal minor and a coefficient minor.

This identifies a precise LGV target:

~~~text
construct one planar network whose boundary measurement is A^N
or A^M, with the finite-strand seed attached to the Q-ray production orbit.
                                                                  (22)
~~~

The proposed vertices are coefficient levels k and reversed times t. The
pairwise observable is literal path incidence; there are no ties and no
tournament gauge. The preserved target is every coefficient minor. Passing
only to sorted root multisets destroys the time-labelled birth/death paths,
as (14d)--(17) demonstrate.

## 6. Why the obvious positive-layer proof fails

Let J(r) be the bidiagonal coefficient map for multiplication by u+r. Adding
a root is a nonnegative Jacobi layer J(r). Deleting a root requires J(r)^(-1),
which is not nonnegative.

The first nontrivial strict example already sees this. Take

~~~text
(alpha_0,alpha_1)=(1,5),          e=(-3,+1),
S=(3,2),                          Q=2.                (23)
~~~

At the root event q=1 -> 0, equation (10) adds three roots at 1 and deletes
one root at 5. The formal elementary-coefficient transfer is

~~~text
T=J(1)^3 J(5)^(-1),

symbol(T)=(1+z)^3/(1+5z)=1-2z+13z^2-... .             (24)
~~~

The negative linear coefficient -2 proves that neither the local root-event
map nor the fixed operator (20) can simply be declared a totally nonnegative
production matrix on the whole coefficient space.

But (24) overstates the failure. The source polynomial F_1 is already
divisible by u+5. Thus J(5)^(-1) is applied only on the face Im J(5), where
it is literal factor deletion. The negative directions of the full inverse
may never be entered by the Gamma orbit.

This is the complementary-exit formulation:

~~~text
prove positivity of the compound maps
  wedge^r[J(births) J(deaths)^(-1)]
only on the nested divisibility faces supplied by the previous cut.       (25)
~~~

Jacobi's complementary-minor identity is the natural local tool: a minor of
an inverse layer becomes, up to the controlled checkerboard orientation, a
complementary minor of the original positive Jacobi layer. Strict ballot
heights say that every deletion has an earlier planar birth available and
that Q boundary paths exit unmatched. What is not yet proved is that these
local complementary choices concatenate without incompatible crossings.

This makes the local-chain/complementary-exit pattern suggested by
[arXiv:1506.00952](https://arxiv.org/abs/1506.00952) mathematically precise
here: the actual chain is (20), the deletion layer is J(r)^(-1) on Im J(r),
and the complementary exits are the Q unmatched ballot paths. No theorem
from that paper is imported, and the analogy has no force beyond this named
map.

## 7. Exact finite signal

The companion
04-computation/gmc_strict_prefix_terminal_q_ballot_transport_scout.py uses
the exact universe

~~~text
alpha_0=1,
one or two integer gaps with offsets at most 4,
all strict heights S_j in {1,2,3},
terminal Q in {2,3},
residual degree at most 5,
column ceiling C=3.                                    (26)
~~~

There are 106 inventories. It verifies:

- 742 direct values of the strand factorization (4)--(6);
- 318 root events (10)--(11), matching 1,527 old roots;
- 318 instances of the fixed production recurrence (19);
- every one of 158,730 minors of A^M, with zero negative minors;
- every one of 158,730 minors of A^N, with zero negative minors; and
- 23,850 independent generalized-Hankel minors of orders two through four,
  all with the checkerboard sign.

There are 61,752 structural zero minors in each coefficient census. They
come from degree support and are not counterexamples. The direct Hankel
minors in the declared universe are all strict.

The script also reproduces the non-strict-prefix control from THM-3065:

~~~text
(alpha_0,alpha_1,alpha_2)=(1/7,2/3,20),
e=(-1,+1,-1),                 E=(-1,0,-1),

det[H_(i+j)]_(i,j=0)^2
 =4914161/84138683904000 >0,                           (27)
~~~

which has the wrong order-three checkerboard sign. Ordinary and optimized
Python are byte-identical to the stored transcript.

## 8. Sharpened conjecture and proof obligations

### Coefficient-network conjecture

Under (1)--(2) with integer shape gaps, for every ceiling C, both coefficient
arrays in (21) are totally nonnegative.

This is stronger than the desired Hankel statement. Its value is that it
names an exact planar-network object. If it is false, the first negative
coefficient minor will reveal which information must be added to the state.

### Strict-prefix Hankel conjecture

Under the same hypotheses, every order r and every pair of strictly
increasing nonnegative integer tuples obey

~~~text
sign det[H_(u_i+v_j)]
 =(-1)^(r(r-1)/2).                                    (28)
~~~

The proved content in this note is the factorization, order-two case,
root-ballot law, production reduction, and both no-go examples. Equations
(26)--(27) are finite exact evidence only.

The cheapest decisive next steps are:

1. compute the compound action of J(b)^A J(d)^(-B) on Im J(d)^B, and
   express each surviving minor by a complementary Jacobi minor;
2. prove the two-shape family symbolically before allowing several nested
   strand intervals;
3. isolate laminar/ribbon strand matchings, where complementary exits cannot
   cross, and prove that subclass by LGV;
4. then test the first non-laminar height profile at ceilings C=4,5,6,
   emphasizing full order-five minors rather than another order-two scan.

A generic Hadamard-product theorem should not be substituted. The proved
Hankel Hadamard closure in
[Fallat--Johnson--Sokal](https://arxiv.org/abs/1612.02210) concerns the
positive Hankel cone, while these reciprocal kernels are already indefinite.
Likewise, a general dual Jacobi--Trudi monomial-positivity claim is not
available: the current replacement
[arXiv:2511.08969](https://arxiv.org/abs/2511.08969) proves ribbon-like
cases, not the general conjecture. Those boundaries make the laminar
subclass in step 3 a lawful target rather than a citation shortcut.

## 9. Relation to the order-five transport circuit

The two order-five observations now have compatible roles.

- The THM-3110 residual bank is one minimal full-support circuit invisible to
  Schur flags through degree four and positive at degree five. Its raw
  coefficient matrix has both signs already at order two.
- The strict-prefix Gamma quotient has a positive finite coefficient signal,
  but its abstract root-event transfer has a negative local inverse
  coefficient already in (24).

Both say that positivity is created after retaining a constrained sidecar.
For THM-3110 the sidecar is the chamber-rooted rank-four circuit. Here it is
the integer birth/death path together with its current divisibility face.
Neither problem is likely to yield to an ambient totally-positive matrix
which forgets how the current state was reached.

The new reusable research move is therefore:

~~~text
when a signed local transition is an inverse of a positive layer,
do not test it on the ambient cone;
identify the inherited image face and move its minors by complementary exit.
                                                                  (29)
~~~

The trigger is an exact divisibility/cancellation inherited from the previous
cut. The counterindication is a generic root-matching family such as (14d),
where no common image face exists.

## 10. Reproduction

~~~text
python3 04-computation/gmc_strict_prefix_terminal_q_ballot_transport_scout.py
python3 -O 04-computation/gmc_strict_prefix_terminal_q_ballot_transport_scout.py
~~~

Both runs must byte-match
05-knowledge/results/gmc_strict_prefix_terminal_q_ballot_transport_scout.out.

LF-normalized SHA-256 at this checkpoint:

~~~text
script  02ff9341a5118f0e68795510020d7cfecf5f1f86c7d8240182228b4c370763c9
output  eb1afd57bfc1db3a5a3340cb3655c1eb3baac16c3b2e46f6f3e17dd879dd425b
~~~
