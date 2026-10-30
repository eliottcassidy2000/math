# Height-adaptive prices close the twelfth layer

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
