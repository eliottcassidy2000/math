# Height-adaptive prices close the twelfth and thirteenth layers

**Latest PROVED result:** the independently checked E10 certificate in C5
proves, for every finite positive valuation word with original clock A,

    S13(w) <= (305175781/488281250) A,

with additional slack12207031/48828125 whenever w has a non-one letter.
The coefficient is sharp. C6 reduces the next E14 layer exactly to two
native price inequalities. The E12 proof below remains valid and explains
the smaller certificate and its artificial-lift boundary. Neither result
proves an all-layer bound or universal Collatz convergence.

**PROVED, using the finite certificate independently checked below:** for every
finite word of positive valuation letters, with its original integer clock A,

    S12(w) <= (30517578/48828125) A.                       (1)

The coefficient is sharp. If w contains any letter greater than one, then

    S12(w) <= (30517578/48828125) A - 2441401/19531250.     (2)

All seven twelfth-layer price families in the preceding packet are therefore
closed. **OPEN:** a uniform all-layer certificate and universal Collatz. These
are content inequalities for valuation words, not a convergence proof. This
packet has a separate executable certificate checker and same-agent proof
review; it does not claim an independent human or subagent audit.

## Inheritance and the new information budget

The closest proved mechanism is the paid-return lemma in
[return_kernel_core_transplants_20261030/COLLATZ.md, section 3](../return_kernel_core_transplants_20261030/COLLATZ.md).
The preceding [terminal fee theorem](../terminal_fees_phase_covers_20261030/COLLATZ.md)
already proves E11. Its observer-only refinement is the corrected near miss:
retaining another multiplier digit without consuming it did not improve the
fee. The least-used sidecar is that packet's
[height-adaptive observer](../terminal_fees_phase_covers_20261030/ADAPTIVE.md).
Here it charges the actual truncated content and pays each initial phase in
its original clock, instead of using an epoch-density enclosure.

For F_w(X)=(P X+B)/Q=mX+t, Q=2^A, put z=t-m+1 and

    c(w)=v5(gcd(P-Q,B)),
    S_e(w)=sum over nonempty prefixes u of min(c(u),e),
    lambda_e=(5/8)(1-5^(-e)), T_e=2*5^(e-1).

All denominator powers of two are units in the rings used below. A letter a
updates m to 3m/2^a and z to (3z+2^a-2)/2^a. Content is equivalently the
minimum of the 5-adic heights of m-1 and z.

## C1. A nondeterministic refinement that preserves every actual word

Fix E and q=5^E. For h=min(v5(z),E), retain z modulo q and m modulo
5^ell, where ell=max(1,h). The multiplier is always a unit. On a letter a,
compute z' exactly modulo q and h'. If the required multiplier precision
falls, reduce the actual multiplier formula. If it rises, allow **all** unit
lifts compatible with the old multiplier residue. Charge

    min(h',v5(m'-1)).                                    (3)

This is exactly the content of every actual transition represented by that
edge. Induction shows every actual word has a path with the same letters,
integer clock and charges. Artificial paths may also occur. An upper bound
proved on all these paths is consequently valid on the actual words; an
artificial path alone cannot refute an actual-word theorem.

The previous packet's hostile at q=25, (z,m)=(1,3),(1,8), a=10, explains
why the lift must branch: the retained states agree, but the next charges
differ. Choosing a single lift would lose an actual source.

There are (16E+84)5^E/25 states: at h=0 there are16q/5; at each
1<=h<E there are16q/25; at h=E there are4q/5. This is much smaller than
retaining both full residues, but does not assert a fixed-size all-E model.

## C2. Paid entry phases have a finite exact reduction

Every non-one word begins 1^k a with a>=2. Immediately after a,

    z=1-2^(1-a) mod q.

Let ell=max(1,h(z)). The target multiplier depends on k only modulo T_ell.
At slope lambda_(E+2), among all k in one such residue, the least
nonnegative representative maximizes the entry excess. In particular there
is no uncharged arbitrary phase.

To prove this, write k=b T_ell+s, 0<=s<T_ell. The initial charge is
sum_(j<=E) floor(k/T_j). For j<=ell its difference from the charge at s is
exactly b T_ell/T_j. For j>ell, T_j/T_ell is an integer and

    floor((b T_ell+s)/T_j)=floor(b/(T_j/T_ell))
                            <=b T_ell/T_j.

The charge difference is at most lambda_E b T_ell, while the additional
clock is b T_ell. Since lambda_E<lambda_(E+2), the longer representative
cannot improve the excess. The charge on a itself depends only on the target
state, so the same dominance holds after entry. This proof also covers E=1.

## C3. The finite certificate and both infinite tails

Put D=200q, N=125q-5, so N/D=lambda_(E+2). The vector V certifies

    V(target)>=V(source)+D*(3)-N*a,
    V(entry)>=D*S_E(1^k a)-N*(k+a),
    max V<=0.                                           (4)

Here `(3)` in the first line means the charge in equation (3), not the
integer three. Entries enumerate a=2,...,64 and 0<=k<T_ell. Transitions
enumerate a=1,...,64 and every compatible multiplier lift. Every charge is
at most E. Thus the two further inequalities

    max V+D*E-65N <= min V,
          D*E-65N <= min V                               (5)

cover every omitted transition and entry letter. For an omitted entry the
initial all-one excess is nonpositive, independently of its phase. C2 and
(5) therefore cover every non-one word, without a bound on its length.

The producer uses queue relaxation and stratified state decoding. The
separate checker enumerates z first, uses repeated modular halving rather
than inverse-power tables, and checks each inequality without optimization.
It rejects wrong lengths, positive or unreached values, and violated edges,
entries or tails. No assertion statement is needed for its checks.

At E=9 the independent checker reports:

    states = 17,812,500; edges = 1,766,250,000; entries = 150;
    D = 390,625,000; N = 244,140,620;
    min V = -5,468,749,720; max V = -488,281,240;
    V(z=0,m=1 mod5^9) = -1,220,703,020.

Both tails pass. Hence every non-one word whose actual endpoint is identity
modulo5^9 has

    S9 <= lambda11*A-eta,
    eta = 61035151/19531250 = 3.1249997312.                (6)

Every non-one word, regardless of endpoint, also has slack
beta=12207031/9765625=2 lambda11 at slope lambda11. The E6 control passes
10,320,000 edges; replacing its terminal value by zero is rejected by the
independent checker. See the audit JSON files, corruption control and binary
hash receipts. The142,500,000-byte E9 potential is a reproducible build
artifact, not a committed data file.

### The terminal margin is optimal for this abstraction

The word (2,2,3,1,2,2,1,1,2,1,4) has clock21 and coefficients

    P=177147, B=1986245, Q=2097152, B-P+Q=2*5^9.

Its actual S9 is2, with terminal content1. Immediately before its last
letter, z is a unit, so only m mod5 is retained. The last step may choose
the artificial multiplier lift m=1 mod5^9, keeping z=0. This abstract path
has charge10, and price21 lambda11-10=eta exactly. Every potential valid
on all paths in this state graph must admit it. Thus changing weights
alone cannot improve this terminal bound at this E and slope. The actual
integer word is not a counterexample: its multiplier is not identity mod25.

`critical_lift.cpp` finds it among2,097,151 words of clock at most21;
`semantic_controls.py` independently checks its literal coefficients.
The targeted repair is to retain m mod25 also at low heights, which rejects
this last lift, or retain a compatible clock/length condition using
m=3^(word length)/2^A. Neither proposal is claimed to give an all-depth
certificate. This explicitly identifies what the next refinement must fix.

## C4. From terminal debt to all seven families

A primitive10-excursion b has S10(b)=S9(b)+1: its only level-ten return
is terminal. Its endpoint is identity modulo5^9. By (6), every such non-one
excursion has price

    lambda11*A_b-S10(b) >= eta-1 > 2.                    (7)

Cut an arbitrary word at actual returns modulo5^10. The inherited boundary
potential Phi(alpha,beta)=((-alpha) mod5)/5 consumes at most one unit on
each non-one excursion, and at most1/5 on each all-one excursion. Thus a
complete non-one excursion leaves eta-2. An unfinished suffix has no
tenth-layer return, so its S11 equals S9. If it contains a non-one, the
global certificate pays it with slack beta. Otherwise its all-one charge
is already below budget. Actual entering maps are identity modulo5^10,
so the lower-layer tests are the native ones; the next tangent is retained.
Consequently

    S11(w)<=lambda11*A-gamma for every non-one word,
    gamma=min(eta-2,beta)=21972651/19531250>1.             (8)

For every native primitive11-excursion, including each prescribed next
tangent X, its fee at slope lambda12 is at least

    g11(X)>=gamma+(lambda12-lambda11) A_b>1.              (9)

Empty tangent classes have infinite price. Every one of the seven families
requires at most one unit total from at least one such excursion, so all
seven hold. The paid-return lemma gives (1). Its quantitative slack version
gives at least gamma-1 (indeed a slightly better clock-dependent amount),
which is exactly (2). The final unfinished suffix is paid by (8); it is not
silently dropped.

At T12=97,656,250 consecutive ones, S12=(5^12-1)/4=61,035,156.
The ratio equals30517578/48828125, proving sharpness in (1). This word has
a realizable positive source; no infinite all-one positive orbit is inferred.

## Reproduction and scope

Run from this directory:

    c++ -O3 -std=c++17 adaptive_fees.cpp -o /tmp/adaptive-fees
    /tmp/adaptive-fees 9 /tmp/adaptive-fees9.bin
    c++ -O3 -std=c++17 audit_adaptive.cpp -o /tmp/audit-adaptive
    /tmp/audit-adaptive 9 /tmp/adaptive-fees9.bin

The checkpoint contains exact audit receipts. E12 is a proved finite-depth
advance. Neither success at E9 nor the smaller adaptive state count implies
feasibility at every E. The next useful obstruction is an artificial path
that gains a high multiplier lift after earlier digits were forgotten. A
repair must constrain that lift while still admitting every actual source.

## C5. The E10 certificate and the sharp thirteenth layer

The same construction at E=10, R=1 passes its separate checker:

    states = 95,312,500; edges = 9,625,000,000; entries = 150;
    D = 1,953,125,000; N = 1,220,703,120;
    min V = -29,052,734,040; max V = -2,441,406,240;
    V(z=0,m=1 mod5^10) = -7,324,217,040.

See[adaptive10_audit.json](adaptive10_audit.json). Both all-letter tails
pass. The binary is762,500,000 bytes and its hash is recorded separately;
it is a reproducible build artifact, not committed data. Replace9 by10 in
the reproduction commands above.

Thus the non-one terminal fee at slope lambda12 is

    eta10=183105426/48828125 = 3.74999912448,

and the global non-one slack is beta10=61035156/48828125=2 lambda12.
A primitive11 non-one excursion pays at least eta10-1. Applying exactly
the C4 boundary argument one level higher, with the actual tangent retained,
gives global non-one E12 slack

    gamma12=min(eta10-2,beta10)=61035156/48828125>1.

The one-fee lift now proves the displayed S13 bound and non-one margin
gamma12-1=12207031/48828125. This deliberately drops a further positive
clock contribution. The unfinished suffix is paid by gamma12, exactly as in
C4, and is never omitted. Sharpness follows at T13=488,281,250 consecutive
ones, whose S13 charge is305,175,781.

## C6. Retaining native fees leaves only two fourteenth-layer obligations

Do not replace a native excursion fee by the weaker global non-one slack.
Cut a non-one primitive12 excursion b into complete11 excursions.
There is exactly one level12 return, at its end. At least one constituent
is non-one and pays eta10-1 at slope lambda12; all other constituents have
nonnegative fee. Therefore

    delta12(b) >= eta10-2 + (lambda13-lambda12) A_b.

Next cut a non-one primitive13 excursion into complete12 excursions.
Every non-one constituent has the preceding positive fee, and an all-one
constituent has fee1/5. Subtract its single terminal level13 charge. This
gives the deliberately weaker uniform bound

    delta13(b) >= eta10-3 + (lambda14-lambda13) A_b
               > 36621051/48828125 > 7/10.              (10)

These cuts preserve the lower identity tests because their entering maps
are identity in the appropriate quotient. Native non-one excursions
necessarily contain at least one non-one subexcursion; their cuts have no
unfinished suffix. Every original letter and clock unit is counted once.

Let g13(X) be the minimum native non-one13-excursion fee at tangent
X in F5^2, or infinity if the class is empty; v=(4,4) is the all-one tangent.
The inherited exact kernel criterion now says that E14 is equivalent to
just

    g13(0,0) >= 1,                 g13(-v) >= 4/5.       (11)

Indeed two non-one excursions cost more than1 by(10). One non-one plus
two or more all-one excursions also costs more than1. Pure all-one zero
sums require at least five steps and cost at least1. The only remaining
possible underpaid zero sums are a single non-one, with zero or one all-one
step. Conversely the native witnesses for either strict failure in(11)
concatenate to an actual violating word. Minima exist by the inherited
positive clock bound and finite branching below any fixed price.

**OPEN:** both inequalities in(11), all-layer content control, and the
additional arithmetic argument needed for universal Collatz. These are two
native-price questions, not a two-word check. The stronger necessary
[five-form grammar](NATIVE14.md) now retains one primitive11 core, exact
placements and finite core-clock bounds. It does not exclude those cores.
Source-specific ROOT positivity is not asserted by these content estimates.

## C7. Minimum multiplier precision as a targeted refinement

The programs also allow an optional minimum precision R,1<=R<=E. Set
ell=max(R,h) in C1 and C2. The same compatible-lift proof and paid-entry
reduction apply unchanged. The state count is

    4*5^(R-1)*5^E + 16(E-R)*5^E/25.

The first term counts all states with h<R together with h=E; each
intermediate h=R,...,E-1 contributes16*5^E/25. The C3 inequalities and
all-letter tails are unchanged. Both executables reject configurations
whose state count exceeds their integer index type instead of overflowing.

At E6,R2, the independent check passes25,080,000 transitions on352,500
states, with630 paid entries and terminal fee3.9999616. The default-R1
regression reproduces the earlier E6 binary byte for byte. R2 excludes the
C3 fake endpoint because it retains m mod25 even while z is a unit.
Equivalently, m mod25 retains the congruence7*len(w)-A mod20, since3=2^7
mod25. This is exact clock/length memory, not an independent random digit.
Finite R2 success does not assert a uniform all-depth refinement theorem.

At E9,R2, the independent checker passes3,611,250,000 edges on47,812,500
states, with630 paid entries. Its terminal fee is2832031/781250=3.62499968.
This improves the R1 fee but is again **optimal for this abstraction**.
The new extremizer is simply four ones followed by the earlier eleven-letter
suffix. Its clock is25, actual charge5, and compatible artificial terminal
lift raises its abstract charge to12. Its integer coefficients are

    P=14348907, B=43294475, Q=33554432.

The actual terminal content is2; before its last letter the retained
precision is2. Thus25 lambda11-12 equals the certified terminal fee exactly.
The exhaustive fixed-length probe checks all binomial(24,14)=1,961,256
compositions of25 into15 positive letters and finds exactly this witness
under its stated single-final-lift conditions. It is not a classification
of every abstract optimizer.

This is a paid-entry version of the cyclic phase arithmetic in
[Session116's phase-completion theorem](../session116_dimension_owner_towers_20261030/collatz/PHASE_COMPLETION_THEOREM.md),
not a new all-entry optimization theorem. Prefixing ones keeps z=0 until
the suffix begins, so the suffix's z trajectory is unchanged. It rotates the
multiplier by(3/2)^k. The least alignment phases for R=1,...,9 are

    0,4,24,74,1074,6074,24824,56074,56074.

[precision_phase.py](precision_phase.py) charges each prefix by its exact
floor sum, verifies every terminal lift, and checks the first two witnesses
again by literal integer coefficients. At R=8 and9 this suffix already has
actual terminal content9; the witness is no longer an artificial terminal
lift. These explicit prices are upper bounds on the optimal endpoint fee
for their respective abstractions. Only R1 and R2 have matching lower
certificates here. Improving precision repairs one missing digit, but a paid
phase shift can recreate the obstruction at the next precision.
