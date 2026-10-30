# The budget counterexample and an entry-preserving first-arrival repair

**REFUTED:** the previous universal quadratic-phase budget. **PROVED:** the
first-arrival identity, its sufficient envelope, and an infinite family reduced
to the existing sharp five-layer theorem. **PROVED by finite certificates:**
the truncated budget at E=1,...,8, for every finite positive word.
The all-layer charge target and universal Collatz entry remain **OPEN**.

**Incoming frontier update:** Session113 subsequently supplied the sharp
repaired E9 certificate. Section5 credits that stronger result; section6
uses the combined mechanisms to obtain exact infinite-family budgets.

## 1. Inheritance and the exact failure

Use [the preceding COLLATZ.md](../layer_lifts_cyclic_intersections_20261030/COLLATZ.md)
for P=3^d,Q=2^A,F_w(x)=(Px+B)/Q, g=gcd(|P-Q|,B),
z=(B+Q-P)/Q, and chronological concatenation F_uv=F_v composed with F_u.
Appending a positive valuation letter a changes

    (P,Q,B) -> (3P,2^a Q,3B+Q),
    z -> (3z+2^a-2)/2^a.

After the first non-one letter, z>0. An epoch consists of a>1 and L subsequent
ones. Put h=v5(z) at its start, kappa=v5(2^a-2). The epoch is active when
m=P/Q is a square modulo5. The old deficit was

    D_old = sum_active max(h,kappa) - sum_epochs 5a/8
          - sum_inactive 5L/8 - sum_active 5L/(8*5^h).       (1)

The proved inequality S_infinity<=5A/8+D_old remains valid. The proposed
universal D_old<=0 does not.

Let U=(3,1,3,1,1,1,1,1,1,3,1,3) and w=U 1^231 U. Then

    length(w)=255, A=271,
    J=18, B_budget=89231/5000,
    D_old=769/5000>0,
    S_infinity(w)=154 < 1355/8=5A/8.                       (2)

Its eight reset heights are0,2,0,5,0,2,0,9. All epochs are active, and all
reset kappas are0. At the last reset m=4 modulo5 and L=0. There is no return
at any of those nine levels. Quadratic membership guarantees that some future
ones can return; it does not put a return inside this epoch. No global minimality
of this witness is claimed. It has the first failing precision found by the
certified truncated searches below.

This repairs a named hypothesis, not a proved result. The earlier exhaustive
clock<=18 test is unchanged. See the correction in `01-canon/MISTAKES.md`.

## 2. The exact arrival address

For an active epoch and1<=e<=h, set t_e=2*5^(e-1). The element3/2 has order
t_e modulo5^e and generates exactly the units whose reduction modulo5 is a
square. This follows from its order2 modulo5 and v5((3/2)^2-1)=1, then LTE.
There is a unique j_e modulo t_e satisfying

    m=(3/2)^j_e modulo5^e.

Let r_e be the representative of -j_e in [0,t_e). It is the first offset,
starting with offset0 at the reset itself, at which multiplier1 can occur.
The valuation of z stays h along ones. Hence the exact epoch charge is

    sum_(e=1..h) N_e(L),
    N_e(L) = 0                         if L<r_e,
             1+floor((L-r_e)/t_e)      if L>=r_e.           (3)

This follows directly from v5(g)=min(v5(m-1),h); it counts the same source
prefixes, with their original integer clock. Inactive epochs contribute zero.
The initial run of L0 ones contributes sum_e floor(L0/t_e).

The addresses are compatible under reduction and can be computed without
searching a full residue cycle: start with j1=0 or1, then test the five lifts
j_(e-1)+k*t_(e-1), k=0,...,4. Exactly one lifts at each level. This is a finite
algorithm for a supplied word, not a uniform bound on its heights.

Define the attained startup envelope

    C(h,m,L) = sum_(e<=h, r_e<=L) (1-r_e/t_e).              (4)

Each summand lies in(0,1]. Formula(3) gives

    charge(epoch) <= lambda_h L+C(h,m,L),
    lambda_h=(5/8)(1-5^-h).                               (5)

This envelope is never worse than lambda_h L+h. In particular it is zero
when the epoch has no attained return level and L=0. With

    D_arr = sum_active C(h,m,L) - sum_epochs 5a/8
          - sum_inactive 5L/8 - sum_active 5L/(8*5^h),      (6)

we have the all-word sufficient inequality

    S_infinity <= 5A/8+D_arr <= 5A/8+D_old.                 (7)

For the hostile word, D_arr=-60951/5000, comfortably negative. The exact
identity(3) and inequality(7) are proved; **universal D_arr<=0 is OPEN**.
The experiment checks all32767 positive words of clock<=15. It neither proves
that open inequality nor replaces actual source legality with a random model.
Because(4) depends on the completed epoch length, a streaming potential must
retain arrival state or unpaid startup credit; it cannot debit all h levels
when an epoch begins.

## 3. The whole hostile macrofamily is already controlled

For the above U,

    (P_U,Q_U,B_U)=(531441,1048576,3664115),
    Z_U=B_U+Q_U-P_U=4181250, v5(Z_U)=5.                    (8)

Consider w_L=U 1^L U for every integer L>=0. The preceding marked-controller
theorem says

    g(u c^k v) divides C(c,vu),
    C(c,t)=B_t(P_c-Q_c)-B_c(P_t-Q_t).                     (9)

Take c=(1),u=U, and v any prefix V of U. Then C(1,VU)=Z_(VU)>0.
For prefix lengths0,...,12 its valuations are

    5,0,0,2,2,2,2,2,2,2,0,0,5.                          (10)

This is a finite list of explicit nonzero integers, independently computed
from affine composition. Positivity also follows because VU has a non-one
letter. During the middle ones, take v empty in(9); the bound is5. The initial
U has actual content height at most5 as well. Thus every prefix of every w_L
has actual content height<=5. The inherited sharp five-layer theorem gives

    S_infinity(w_L)=S_5(w_L)<=781(40+L)/1250               (11)

for **all** L. The proof uses the exact entry U and chronology VU, not UV.
The finite L=0,...,1249 checks corroborate it; they are not the proof for
unbounded L.

The final defect illustrates the mechanism particularly clearly:

    z_(U 1^L U)=z_U(1+m_U(3/2)^L).                       (12)

There are arbitrarily high final defect heights by the same cyclic-unit
lifting. Whenever the parenthesis is0 modulo5, m_U=1 modulo5 implies the
actual final multiplier is -m_U=4 modulo5. The high defect therefore has no
instantaneous content charge. The controller divisor provides a stronger
uniform safeguard than a height heuristic.

## 4. Eight levels certify; the ninth separates the models

For E fixed, the truncated budget replaces h,kappa by min(h,E),min(kappa,E).
It depends only on (z modulo q,active), q=5^E, a set of2q states. Initial ones
leave this state equal to(0,true); their spent clock can be ignored safely.
In units8q, a one-letter transition incurs

    delta = -5q/5^h_E      if active,
            -5q           otherwise;

and a reset incurs

    delta = 8q [target active] max(h_E(target),kappa_E(a))-5q*a. (13)

`epoch_potential.cpp` constructs a vector V with V<=0, V(seed)>=seed debt,
and V(t)>=V(s)+delta on every edge. This proves the truncated budget for all
finite words by induction. A separate NumPy implementation reads the exported
integer bytes and checks every inequality, seed and tail. The full vector is
reproducible and deliberately not committed as a binary blob.

| E | states | checked letter cap | independently checked edges |
|---|---:|---:|---:|
| 1 |10|12|120|
| 2 |50|16|800|
| 3 |250|20|5000|
| 4 |1250|24|30000|
| 5 |6250|28|175000|
| 6 |31250|32|1000000|
| 7 |156250|36|5625000|
| 8 |781250|40|31250000|

For a>cap, delta<=8qE-5qa. The audited inequality
5q(cap+1)>=8qE-min V makes every larger edge and seed automatic, since max V<=0.
Thus the alphabet is unbounded; only the precision is fixed.

These also imply S_E<=5A/8 for E<=8. They do not give the sharper coefficient
lambda_E, nor control E->infinity. At E=9 the program finds exactly(2), whose
heights are all<=9, so its positive truncated deficit is the actual deficit.
The small quotient is sufficient for this budget consumer but insufficient
for the stronger future-return observer of the preceding packet; there is
no contradiction with that observer's larger state lower bound.

## 5. Merge the stronger incoming phase certificate

Concurrent Session113's
[one-digit repair](../session113_symmetry_consumers_20261030/uniform_content/REPAIRED_CERTIFICATE.md)
retains m modulo5 and a pending first-arrival flag. Its independently checked
certificate gives, with lambda9=488281/781250,

    S_9(w)<=lambda9 A,
    S_9(w)<=lambda9(A-2) if w contains a non-one letter.    (14)

The receipt checks525,000,000 transitions,126 entries, valid-domain closure
and infinite-letter tails. Both bounds are sharp in their declared domains.
The proof and auditor were read here; that full certificate was not redundantly
rerun. It supersedes our weaker E<=8 content estimates as the current finite
precision frontier. Our old-budget certificates remain useful independent
boundary controls, not the strongest current content theorem.

Its low-phase envelope F_r counts levels<=r exactly and postpones the higher
startup cost until the first possible low-phase return. This gives an online
bounded-mode implementation, unlike an arbitrary completed-epoch formula.
The all-E feasibility of that envelope remains OPEN. Our full addresses in
section2 identify every individual arrival, but are not claimed to yield a
smaller or universally paid online controller.

Session113 also supplies R=(5,1^8,5), W=R 1^233 R, whose word W 1^8 W has
clock546, defect height11 and actual content317. It refutes the old budget
by479371/312500. The next section controls the *whole* W 1^L W family.

## 6. Exact periodic charge with a fixed entry and exit

**Theorem (marked single-ones-loop profile).** Fix words u,v, with at least
one non-one letter in u. For every L>=0 the complete untruncated charge
S(u 1^L v) is an affine-periodic function of L, with explicitly computable
period and gain. This assertion concerns the finite word with its literal
entry and exit, not one integer orbit realizing all choices of L.

Let h=v5(Z_u), and for each nonempty prefix V of v put

    c_V=v5(Z_(Vu)),  t_e=2*5^(e-1).

Every Z_(Vu)>0. Let K=max(h,all c_V); if v is empty only h is needed.
For K=0 use T=1 and gain0; otherwise use T=t_K. Write chi_u=1 if m_u is
square modulo5 and0 otherwise. Then

    S(u 1^(L+T) v)-S(u 1^L v)=chi_u lambda_h T.           (15)

The right side is an integer. Fixed initial-prefix charges do not affect
the period. The claimed period need not be minimal.

The reset assumption is essential. With an all-ones entry and no exit, the
successive charges v5(3^n-2^n) are unbounded; a finite affine period for their
partial sums would make those increments periodic and bounded. The initial
zero-defect stratum therefore requires its own density ledger.

**Proof and exact formula.** The fixed prefix u contributes S(u). At its
endpoint, ones have charge at level e<=h precisely at positive indices i
with m_u(3/2)^i=1 modulo5^e. This is a single residue class if chi_u=1,
and empty otherwise. Next, at the endpoint after the additional prefix V,

    m'=m_V m_u(3/2)^L,
    z'=(Z_u/P_u)m'+z_V.

Conditional on m'=1 modulo5^e, z'=0 modulo5^e is equivalent to

    Z_u/P_u+z_V=0 modulo5^e,

whose valuation is exactly c_V: its numerator over the unit denominator
P_u Q_V is P_u Z_V+Q_V Z_u=Z_(Vu). Thus that suffix prefix contributes

    sum_(e<=c_V) [m_V m_u(3/2)^L=1 modulo5^e].            (16)

The terms are either one explicitly lifted residue class modulo t_e or
empty, according to the multiplier's quadratic coset. This formula preserves
every suffix prefix, including phases with no possible charge. Shifting L by
T leaves every suffix indicator unchanged and adds exactly T/t_e returns
in each active middle level. Sum them to obtain(15). QED.

The formula reduces an unbounded-length question to a declared finite
phase bank. It does not give a uniform bound on K or on the number of
controllers needed for changing-word coverage.

### Two exact optimal macrofamily budgets

For U from section1 and W from section5, the phase caps are:

* U: entry height5; thirteen divisors including the empty suffix have
  heights0 four times,2 seven times,5 twice. K=5.
* W: entry height8;254 divisors including the empty suffix have heights
  1 eighteen times,4 234 times,8 once and9 once. K=9; only the full suffix
  attains9. Every initial prefix also has content height at most9.

Consequently the incoming theorem(14) already certifies W 1^L W for every L,
including cases with arbitrarily high final defect. The exact periodic
profile gives still stronger bounds, without needing that global E9 theorem:

| family | clock A | period T | gain per period | optimal uniform bound |
|---|---|---:|---:|---|
| U 1^L U |40+L|1250|781| S<=781 A/1250-6906/625 |
| W 1^L W |538+L|781250|488280| S<=48828 A/78125-1461937/78125 |

Equality holds precisely at L=962 modulo1250 in the first row and
L=491466 modulo781250 in the second. Thus the constants are optimal for
these families and slopes. The second slope is lambda8, even though a
suffix can reach content height9: the ninth level contributes a periodic
terminal correction rather than a continuing density.

The finite maximum checks are exact certificates for these all-L statements.
`macro_profile.py` builds the whole phase profile using(16) and arithmetic
progression counts; a separately implemented C++ program follows the modular
affine trajectory and every suffix letter. All782500 phase values agree.
Literal unbounded-integer gcd checks cover short positive controls, both
known hostiles and the entry words. The hash, unique maximizers and integer
cap histograms are in `macro_profile.json`; no floating point is used.

This is a useful source-preserving return-controller theorem, rather than
another conjecture justified only by a search cutoff. The remaining universal
obligation is payment when entries and exits themselves change without a
bounded, independently certified phase bank.
