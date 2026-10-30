# Joint first arrival and exact recursive resonant echoes

**PROVED scoped, same-agent proof review:** a joint-arrival criterion for
any finite resonant tail bank above the complete near-interval chamber;
an all-q, all-length recursive family with exact maximum and exactly two
maximizing times. **REFUTED and repaired:** the unrestricted local binding
balance in THM-724. General LRC(14) remains **OPEN**.

Write ||x|| for distance to the nearest integer and
M(V)=max_t min_(v in V)||vt||. Speeds are distinct positive integers; the
stationary observer is implicit.

## 1. Inherited complete windows

For q>=3 let C_q={1,...,q-1,q+1}, and let s be the least integer>=2
coprime to q (equivalently the least missing prime). The incoming
[Session114 general-core theorem](../session114_continuation_symmetry_20261030/lrc/GENERAL_CORE.md)
proves M(C_q)=1/q and describes its entire safe set whenever

    1/(q+2) < rho <= 1/q.

At each primitive centre a/q its left/right windows have directional
owners p,n: p=a^(-1) mod q and n=q-p, with an owner1 replaced by q+1.
A side with owner e has length (1/q-rho)/e. The smallest owner is s,
achieved on exactly two reflected sides. This complete-window theorem is
a proved dependency; merely knowing slopes at one maximum would not suffice.

## 2. Any resonant tail bank: exact joint arrival — PROVED

**Concurrent precedence:** `f128982d82` published this same criterion during
this session, as L1 in
[the uniform-tail-bank packet](../odd_eighths_uniform_families_20261030/LRC.md).
It also states the broader complete-window version for arbitrary cores.
The derivation below is a credited independent replay, not a separate new
claim. Sections3–4 extend it to exact recursive maxima after successive
attachments which deliberately destroy the old maximum.

Let B be a finite nonempty set of positive integers and append qB to C_q.
Define the earliest joint safe time

    tau_B(rho)=min{x>=0 : ||m x||>=rho for every m in B},

with infinity if the set is empty. Periodicity lets us use x in[0,1].
For rho in the chamber above,

    M(C_q union qB)>=rho
      iff q rho+s tau_B(rho)<=1.                    (1)

**Proof.** At a centre a/q, a displacement delta in either direction gives
||q m(a/q +/- delta)||=||m x|| for x=q delta. A core side owned by e
allows exactly e x<=1-q rho. Since all e>=s and a side with owner s
exists, existence in any complete window is equivalent to a joint-tail
point x with s x<=1-q rho. The earliest such point is optimal. QED.

This is a source-preserving existence test: it constructs the actual time
on an actual least-owner side. It is not a mean measure bound.
For B={m}, tau=rho/m; (1) recovers the inherited single-tail formula in
its strict top chamber.

**Synchronization hostile.** For B={1,3}, rho=3/10, the latest of the two
individual first arrivals is3/10, while the joint arrival is13/30.
Their individual clearances cannot be optimized separately and then glued.

## 3. Recursive echoes that destroy the previous maxima — PROVED

Fix q>=3 and an integer m with 2m>s. Set

    a_0=m, b_0=0,
    a_(j+1)=q a_j+s,       b_(j+1)=q b_j+1,
    c_j=a_j-s b_j=q^j m,
    D_j=q a_j+s,
    V_j=C_q union {q a_0,...,q a_j}.

Then for every j>=0,

    M(V_j)=rho_j=c_j/D_j.                            (2)

There are exactly two maximizing times. If a*=s^(-1) mod q in{1,...,q-1},
put x_j=(q b_j+1)/D_j; the times are

    t_j=(a*-x_j)/q,       1-t_j.                     (3)

Every newly appended tail q a_(j+1) has clearance zero at both preceding
maximizers, since a_(j+1)x_j=q b_j+1 is integral. Nevertheless,

    rho_j decreases to rho_infinity
       =m(q-1)/(q(m(q-1)+s)) >=1/(q+2),              (4)

and rho_j>1/(q+2) at every finite j. Thus the entire proof stays inside
the complete chamber; it never extrapolates a germ beyond an owner change.

### Closed forms and the common entry threshold

The recurrences give

    a_j=q^j m+s(q^j-1)/(q-1),  b_j=(q^j-1)/(q-1).

These imply (4), strict finite monotonicity, and

    x_j < x_infinity=1/(m(q-1)+s).

To check the lower bound in (4), it suffices that2m(q-1)>=q s.
For s=2, m>=2 and q>=3 give this. For odd s, integer2m>s gives
2m>=s+1, while s<=q-1, so (s+1)(q-1)>=q s.

Define

    theta_j=c_j/((q+1)a_j+s).

The closed forms show theta_j is nonincreasing, and

    rho_infinity > theta_0=m/((q+1)m+s),

because m(q-1)>s. The strict inequality follows from2m>s and q>=3.

### Global upper bound by first-arrival induction

For every rho>theta_0 and every x>=0 satisfying all tail inequalities
||a_h x||>=rho for0<=h<=j, we claim

    x >= (b_j+rho)/a_j.                              (5)

Base j=0: no safe x occurs before rho/m. If (5) holds at j, then

    a_(j+1) x >= (q a_j+s)(b_j+rho)/a_j
                  > b_(j+1)-rho,

where the strict inequality is exactly rho>theta_j. The next safety
condition forbids the open interval(b_(j+1)-rho,b_(j+1)+rho), hence
forces a_(j+1)x>=b_(j+1)+rho. This proves (5) at j+1, without assuming
that x stays in a selected first period of any individual tail.

Suppose M(V_j)>rho_j. Choose rho strictly between rho_j and that maximum;
then rho<=1/q and rho>rho_infinity>theta_0. The complete-window criterion
(1) supplies a joint safe x with q rho+s x<=1. Inserting (5) gives

    (q a_j+s)rho <= a_j-s b_j=c_j,

contradicting rho>rho_j. This is a global upper bound over all times.

### Attainment and uniqueness

At x=x_j, q rho_j+s x_j=1 and a_j x_j-b_j=rho_j.
For every h<=j,

    a_h x_j-b_h>=rho_j

is equivalent to rho_j<=rho_h, hence holds, strictly for h<j.
Also x_j<x_infinity, and

    a_h x_infinity-b_h=m/(m(q-1)+s)<1/2.

Thus all tail phases are in[rho_j,1/2), with no unrecorded fold.
The least-owner core side has clearance rho_j at that displacement, proving
(2) and (3). If any time attains rho_j, (5) and the core deadline force
x=x_j and owner e=s. There are precisely two such reflected sides, proving
that (3) lists all maximizing times. QED.

### Examples and relevance to fourteen runners

- q=3,m=2: V_0={1,2,4,6}, then append24,78,...; maxima1/4,3/13,9/40,...
- q=4,m=2: append8,44,188,... to{1,2,3,5}; maxima2/11,8/47,32/191,...,
  tending to1/6.
- q=6,m=3: first tails18,138,858; maxima3/23,18/143,108/863,...,
  tending to1/8.
- q=11,m=2: C_11 plus22,264 has13 moving speeds and exact maximum11/133,
  at61/133 and72/133.

Generally, q+j+1 is the number of moving speeds. Setting j=12-q for
3<=q<=12 gives explicit thirteen-speed families with an exact margin above
1/14. Their safety may follow from other results as well; the contribution
here is the exact all-length spectrum, complete maximizer set, and recursive
joint-arrival mechanism, not a claim of first proving those individual sets.

## 4. Arbitrary changing multipliers and a digit expansion — PROVED

The fixed recursion is not necessary. Choose **any** sequence of positive
integers K_0,K_1,..., and instead use

    a_(j+1)=K_j(q a_j+s),  b_(j+1)=K_j(q b_j+1).

Keep a_0=m,b_0=0 and the same V_j. Define

    P_0=1, P_(j+1)=q K_j P_j,
    H_j=sum_(h=0)^j 1/P_h.

Then the exact result, with exactly two finite-stage maximizers, is

    M(V_j)=rho_j=m/(q m+s H_j),
    x_j=H_j/(q m+s H_j),
    t_j=(a*-x_j)/q, 1-t_j.                         (6)

Every new speed still vanishes at both preceding maxima. The whole sequence
K_j is unrestricted; no bounded multiplier or periodicity assumption is
hidden in this statement.

**Coefficient identities.** Let H_(-1)=0. Direct induction gives

    b_h=(P_h/q)H_(h-1),
    a_h=P_h(m+(s/q)H_(h-1)),   a_h-s b_h=m P_h.

These yield (6)'s candidate rho and x. Moreover,

    1<=H_j<q/(q-1),

so every rho_j is strictly above the lower limit in(4), hence remains in
the complete core chamber. It is strictly decreasing.

**Global first-arrival proof.** At step h the sufficient induction threshold
from (5) is now

    theta_h=K_h(a_h-s b_h)/(K_h(q a_h+s)+a_h).

For any later finite stage j>h,

    H_j-H_h < 1/((q-1)K_h P_h).

Since a_h(q-1)>s, this gives

    q m+s H_j < q m+s H_h+a_h/(K_h P_h),

or rho_j>theta_h. Thus every test level rho>rho_j exceeds every induction
threshold encountered before j. The same forbidden-interval induction
proves x>=(b_j+rho)/a_j, and the core deadline again gives rho<=rho_j.
This proves the global upper bound, with all changing multipliers retained.

**No old tail has silently folded.** At the candidate x_j,

    a_h x_j-b_h
      = m P_h (H_j-H_(h-1))/(q m+s H_j).

For h<=j, H_j-H_(h-1)>=1/P_h, giving the lower bound rho_j, strict
for h<j. The same sum is strictly less than q/((q-1)P_h), giving

    a_h x_j-b_h < 1/(q-1) <= 1/2.

Thus every old tail is actually safe and unwrapped. Attainment and uniqueness
follow exactly as in section3. This also supplies a simpler independent
verification of that section's fixed-generator formulas. QED.

### Infinite bank and irrational exact clearances

H_infinity=sum_(h>=0)1/P_h converges, and the infinite speed bank
V_infinity=union_j V_j has the exact common-clearance value

    max_t inf_(v in V_infinity)||vt||
       = m/(q m+s H_infinity).                       (7)

The upper bound follows from every finite-stage maximum. The times t_j
converge to t_infinity=(a*-H_infinity/(q m+s H_infinity))/q. For each fixed
speed, continuity and its finite-stage lower bounds show that t_infinity
and its reflection attain(7). We do not claim these are the only infinite-
bank maximizers; finite-stage uniqueness alone would not imply that.

To obtain an explicit irrational example, set n_0=0 and n_h=h! for h>=1,
and K_h=q^(n_(h+1)-n_h-1). All exponents are nonnegative, P_h=q^(h!)
for h>=1, and H_infinity=1+sum_(h>=1)q^(-h!). Its base-q digits have
infinitely many1s separated by unbounded zero runs, so are not eventually
periodic; H_infinity is irrational. The nonconstant rational fractional
linear expressions in(7) and t_infinity are therefore irrational too.
Finite truncations still have the exact rational certificates(6).

This is a literal digit-block construction of infinite compatible times,
not a transfer from an almost-everywhere statement to an integer source.
Its distinction from Collatz is decisive: the LRC time variable lies in a
compact circle. A compatible infinite valuation address need not encode a
positive ordinary integer. Compactness alone cannot supply that missing
Collatz source.

### Reusing an echo as a new attachment context

The incoming uniform-tail-bank L2–L3 classification also extends to every
finite echo, but its coprimality guard must be recomputed. Write

    D=q a_j+s, c=a_j-s b_j,
    h=(alpha s-1)/q, K=alpha a_j-b_j+h,

where alpha=s^(-1) mod q. Then t_j=K/D. Put

    g=gcd(K,D)=gcd(c,q b_j+1),
    D'=D/g, K'=K/g, c'=c/g.

The identity for g follows from qK=alpha D-(q b_j+1), gcd(q,D)=1, and
D=s(q b_j+1)+q c. In particular all displayed quotients are integers and
K' is a unit modulo D'. For any finite added bank Z,

    M(V_j union Z)=rho_j
       iff c' <= v K' mod D' <= D'-c' for every v in Z.        (8)

This is the complete condition, since every echo has exactly the two
reflected maximizers. There are D'-2c'+1 admissible residue classes modulo
D'. Failure strictly lowers the maximum by compactness, not necessarily
below the LRC threshold.

The actual terminal safe-mass coefficient for a passing bank is

    2(1/P+1/N),
    P=max({s} union {v in Z : v=s mod D'}),
    N=max({q a_j} union {v in Z : v=q a_j mod D'}).

At the chosen maximum the two original binding walls are s and q a_j;
all earlier tails have strict slack. The residue condition lists exactly
which new speeds tie each wall, and their largest actual speeds control the
local expansion. This is the incoming L3 proof applied with the reduced
new denominator. Maximum-preserving additions can again drive this
coefficient to zero without moving either maximizing point.

**Coprimality hostile:** q=3,m=2,j=1, K_h=1 gives D=26,K=16,c=6 but
g=2. The true marking is8/13 and rho=3/13. It admits8 residues modulo13,
or16 residues modulo26. Treating K as a unit modulo26 would incorrectly
count15. The iteration has changed the allowable quotient; the gcd is a
required sidecar, even though it was1 at the initial single-tail stage.

## 5. Corrected near miss

THM-724 / `01-canon/theorems/THM-724-single-killer-covering-min-rigidity.md`
formerly asserted M(A union{v})>=mu v/(v+s) from a local binding slope s.
For A={1,6}, mu=3/7, s=1, v=7, this predicts3/8 but the exact maximum is
4/13. Its proof needed the missing whole-core horizon eta. The valid repair
is v min(eta,mu/(v+s)); globally one can use s=max(A) instead.

The separate multi-tail extension fails at core{1,2}, tails3 and12:
M=3/13<1/4. The fastest tail has already folded. These failures explain why
(1) retains a joint arrival and why the proof of (2) establishes a complete
chamber and an unwrapped interval for every old tail.

## 6. Verification and next targets

Run `python3 lrc_audit.py`, also under `-O`. There are684 exact structural
controls for q=3..40, two m values per q, and j=0..8;21 complete safe-window
checks,18 independent opposite-slope clock maxima,480 randomized but seeded
exact joint-arrival equivalence checks, and the three required hostiles. Another486 literal changing-generator
controls and54 complete-window certificates test section4, with multipliers
in{1,2,7} and all length-three choices at the declared q values.
The finite universe and parameters are explicit in the script; randomness
selects test inputs only. All comparisons are rational and exact.

For each complete-window test, emptiness at rho+1/(8 max(V)^2) certifies no
larger maximum: any positive local maximum is an opposite-slope intersection
or a half-integer cusp, with reduced denominator at most2 max(V). Distinct
such rational values differ by at least1/(4 max(V)^2). Exact singleton
safe intervals at rho also certify the complete maximizing set. A second
implementation enumerates all pair clocks for the smaller declared cases.

**OPEN helper:** characterize short joint-arrival recurrences beyond the
proved arbitrary-multiplier family a_(j+1)=K_j(q a_j+s). The induction
identifies the precise condition: the previous arrival must lie strictly
past the left forbidden endpoint of the next chosen integer address.
**OPEN main extension:** retain these address inequalities across changes
of core window and owner; the current proof does not cover arbitrary cores.
