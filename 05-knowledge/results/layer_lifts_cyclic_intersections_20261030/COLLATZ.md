# Five-adic potentials, defect epochs, and marked-controller divisors

**PROVED by an exact finite certificate:** the sharp five-layer inequality
below for every finite positive valuation word. **PROVED elementary:** the
untruncated defect bound and marked-controller divisor theorem. **OPEN:**
the all-layer sharp potential, uniform changing-controller debt, and Collatz
convergence. No implication from finite tests to those open claims is used.

## 1. Inherited state, clock, and the fifth certificate

Use the conventions and renewal identities of
[the preceding packet](../anchor_content_flat_quotients_20261030/COLLATZ.md).
For a nonempty word w of positive valuation letters, let

    F_w(x)=(P_w x+B_w)/Q_w,
    P_w=3^|w|, Q_w=2^A, A=sum(w), L_w=P_w-Q_w,
    g(w)=gcd(|L_w|,B_w).

Chronological concatenation satisfies F_uv=F_v composed with F_u. Appending a
letter a updates (P,Q,B) to (3P,2^a Q,3B+Q). Content g is coprime to6.
The state (m,b)=(P/Q,B/Q) modulo5^E retains **both** affine coefficients.
Its transition and reward are

    T_a(m,b)=(3m/2^a,(3b+1)/2^a),
    R_E(m,b)=sum_(e=1)^E [m=1,b=0 modulo5^e].

Write S_E(w)=sum_i min(v_5(g(w_i)),E), where w_i are the nonempty prefixes.
The new exact result is

    S_5(w) <= (781/1250) A.                         (1)

This holds for every finite positive word, at every length, with no additive
constant. It is optimal, including in the universe of words realized by
positive integer sources. The older E=1,2,3,4 coefficients remain
1/2,3/5,31/50,78/125.

**Certificate and proof.** On the 7,812,500 affine states modulo3125, the
integer vector H has minimum0, maximum H(1,0)=5469, and satisfies

    H(T_a s)-H(s) >= 1250 R_5(T_a s)-781a.           (2)

The producer checks all states and a=1,...,20. The independent NumPy auditor
checks the exported bytes by a different transition implementation: exactly
156,250,000 integer inequalities, without floating point. For a>=21,
781a>=16401>5469+1250*5, so (2) is automatic. Summing (2) from the identity
gives 1250 S_5-781A<=H(end)-H(identity)<=0.

The byte hash, zero slack count, equality control and reproduction commands
are in [the receipt](potential5_audit.json) and [AUDIT.md](AUDIT.md).
The vector is reproducible rather than checked in as a31MB binary. Independent
O0 and O2 producers give identical bytes. This is a finite certificate for
an all-word theorem at one fixed prime-power precision, not a bound on E.

Equality is attained by1250 consecutive letters1: its content at prefix i is
3^i-2^i, whose five-adic valuation is0 at odd i and1+v_5(i/2) at even i.
The total truncated reward is781. The positive source2^1251-1 realizes this
finite valuation word. No infinite positive all-ones orbit is asserted.

## 2. An exact lift profile, and what remains unproved

Use unreduced D_E=8*5^(E-1), N_E=5^E-1. For the actual computed potentials,
rescaled to these units, put

    K_E=H_(E+1)-5 H_E composed with reduction,
    sigma_E(s,a)=H_E(T_a s)-H_E(s)-D_E R_E(T_a s)+N_E a.

The inherited exact lift equation is

    K_E(T_a s)-K_E(s)
      >=5D_E [T_a s=identity modulo5^(E+1)]-4a-5sigma_E(pi(s),a). (3)

All four tested lifts have precisely the following set of K values:

    {4-D_E, 8-D_E, ..., 5D_E-4}.                    (4)

| lift | min K | max K = K(identity) | distinct K values |
|---|---:|---:|---:|
| 1 to2 | -4 | 36 | 11 |
| 2 to3 | -36 | 196 | 59 |
| 3 to4 | -196 | 996 | 299 |
| 4 to5 | -996 | 4996 | 1499 |

These are **FINITE-EXACT** profiles, not an induction proof. A concrete
failure of deleting old slack occurs at E=4, s=(2,1564), a=1. The next state
is(3,784), K changes124 to0, and the inequality with old slack omitted has
deficit120. The old slack is280, so retaining its factor5 repairs the
inequality with margin1280. This is a failure for this chosen lift, not an
impossibility theorem for every other potential.

**OPEN strengthened target.** Construct H_E for every E, with identity a
maximum, and corrections satisfying (3) and the uniform range in (4).
The range alone does not imply (3). Feasible potentials at every E would
prove S_E<=(5/8)(1-5^-E)A and hence the untruncated bound

    sum_i v_5(g(w_i)) <=5A/8.                       (5)

The implication is one-way here. The earlier wording calling these
statements equivalent has been repaired; no converse has been supplied.
Even (5) would control only one part of the fixed-source Collatz obligation.

## 3. Untruncated epochs: the location of the missing payment

Define the integral defect Z=B+Q-P and its normalized form z=Z/Q=b+1-m.
Then

    Z'=3Z+(2^a-2)Q,   z'=(3z+2^a-2)/2^a.           (6)

Initially Z=0. It remains zero along letters1; after the first a>1 it is
strictly positive. During every subsequent run of letters1, h=v_5(Z) is
constant. At every prefix,

    v_5(g)=min(v_5(m-1),v_5(z)),                   (7)

with the usual convention for the zero defect. Thus a large defect height
does not by itself create a return.

Index the non-one letters a_j and set

    h_j=v_5(Z immediately after a_j),
    kappa(a)=v_5(2^a-2),
    e_j=max(h_j-kappa(a_j),0).

For a>1, kappa(a)=0 unless a=1 modulo4; in the latter case
kappa(a)=1+v_5((a-1)/4)<=(a-1)/4. Equation(6) and the ultrametric inequality
show that h_j>kappa(a_j) can occur only at an equal-valuation cancellation
between3z and2^a-2. These e_j record valuation jumps, not actual return charge.

**All-word epoch bound.** If o is the number of letters1, then

    S_infinity(w) <=5o/8+sum_j h_j
      <=5A/8+sum_j e_j-sum_j(5a_j/8-kappa(a_j))
      <=5A/8+sum_j e_j-sum_j(3a_j+2)/8.             (8)

**Proof.** Partition into the initial ones and epochs consisting of a_j>1
followed by L_j ones. At level e<=h_j, the multiplier advances by3/2 on an
epoch. Its order modulo5^e is2*5^(e-1), by LTE applied to9/4. Among the
L_j+1 positions starting just after a_j, at most floor(L_j/(2*5^(e-1)))+1
can have multiplier1. Sum over e to get5L_j/8+h_j. The initial ones have
exactly floor(L_0/(2*5^(e-1))) returns at level e, hence cost at most5L_0/8.
Now h_j<=kappa(a_j)+e_j and the displayed kappa estimate prove(8).

**PROVED all-layer subgrammar:** any word with at most one letter different
from1 satisfies(5). Its first defect has h=kappa, so its jump debt is zero.
In fact it satisfies the sharper S_E<=lambda_E A for every E: truncate the
geometric sum in the proof at E and use kappa<=(a-1)/4<=lambda_E a.
More generally the explicit acceptance test
sum e_j<=sum(5a_j/8-kappa(a_j)) suffices for(5). It is not universal acceptance.

**Sharp reason to retain phase and unused clock.** For
w_r=(2,1^(2r),3), the final normalized defect and multiplier are

    z=3(3^(2r)+2^(2r+2))/2^(2r+4),
    m=(9/32)(9/4)^r.

Since9/4 generates1+5Z_5 at each finite precision, there is a unique
r modulo5^(e-1) with (9/4)^r=-4 modulo5^e. Hence h can be arbitrarily large
while the two non-one letters cost only5. But m=2 modulo5 at that endpoint;
further ones alternate m between2 and3. Before the last letter z is a unit.
Consequently **the entire word, and every extension by ones, has zero actual
five-adic content charge**. The first r values for e=1,...,7 are
0,1,21,46,296,921,13421. At r=46 the unrefined jump acceptance test fails
(debt4>25/8), although actual charge is0. This is a false negative of the
sufficient grammar, not a counterexample to(5).

Another concrete quotient failure is w=(3,1^16), A=19. Counting only
v_5(P-Q), truncated at3, gives12>31*19/50; actual content charge is0.
The old trace-blindness theorem already identifies the discarded carry:
[THM-4604](../../../01-canon/theorems/THM-4604-collatz-words-are-the-reducible-locus-carries-are-the-extension-cocycle-and-traces-only-hear-the-tuning.md).
These witnesses specify a useful next refinement: a defect epoch must retain
its multiplier coset and saved clock slack, not just its height.

### A refined budget that repairs the dark-height obstruction

At an epoch's start call it active if m is a square modulo5, equivalently
A+|w| is even at that prefix. This equivalence follows from the Legendre
symbols of2 and3, both-1 modulo5. Ones preserve this flag. Inactive epochs
have zero return charge, however large h is. For an active epoch the earlier
geometric sum is actually lambda_h L+h, where

    lambda_h=(5/8)(1-5^-h), including lambda_0=0.

Thus define the total **active jump debt** and **available budget** by

    J=sum_(active j) max(h_j-kappa(a_j),0),
    B=sum_j 5a_j/8
        +sum_(inactive j) 5L_j/8
        +sum_(active j) ((5/8)5^(-h_j)L_j-kappa(a_j)). (13)

The same proof gives the exact sufficient inequality

    S_infinity(w) <=5A/8+J-B.                        (14)

The initial all-ones run is treated separately by its5/8 bound. The budget
retains unused time in low-height runs, discounted by their actual height.
It is nonnegative before the comparison with J: each kappa is paid by its
own non-one letter. Every dark-height word w_r above has J=0 and is accepted.

**REFUTED (2026-10-30):** the hypothesis J<=B for every positive word fails.
Let U=(3,1,3,1^6,3,1,3) and w=U 1^231 U. Then A=271,
J=18 and B=89231/5000, so J-B=769/5000>0. Its actual content charge
is154<1355/8. Thus neither(14) nor the underlying content target(5) is
refuted. The final height-nine epoch has multiplier4 modulo5 and no
following ones: none of its nine possible return levels is actually visited.
The exact word and first-arrival ledger are in
[the successor packet](../phase_budget_shared_orbits_20261030/COLLATZ.md).
Concurrent [Session113](../session113_symmetry_consumers_20261030/collatz/THEORY.md#6-exact-boundary-the-epoch-budget-fails-already-at-layer9)
independently audits this boundary, gives a second resonant-doubling hostile,
and supplies a [sharp repaired E9 certificate](../session113_symmetry_consumers_20261030/uniform_content/REPAIRED_CERTIFICATE.md).
This was a sufficient target, not an established equivalence. Earlier exact tests
accept all262143 nonempty words of binary clock at most18. A distinct search
finds no positive debt path in the finite truncated-height graphs E=1,2,3
with letters1,...,12. Those graphs do not certify the untruncated budget.
Their states, precision and receipt are in [epoch_budget.json](epoch_budget.json).
The phase gate, discounted clock and inequality(14) remain proved refinements.
The finite checks remain valid in their declared universes; they did not
control higher levels. Retain the failed word as a regression test.

## 4. A fixed-controller theorem at every prime, with the entry retained

For nonempty w and any word t define the integer

    C(w,t)=B_t L_w-B_w L_t.                          (9)

The empty word has P=Q=1, B=L=0. The scalar C is the numerator of the
translation difference between F_t F_w and F_w F_t, up to sign. The faithful
word representation is inherited from THM-4604; the following divisibility
is the new use of that mechanism.

**Marked-entry divisor theorem.** For words u,v, every k>=0 for which the
word u w^k v is nonempty satisfies

    g(u w^k v) divides C(w,vu).                     (10)

If C is nonzero this is an untruncated, all-prime bound, independent of k.

**Proof.** Work modulo g=g(u w^k v), where every affine map is invertible
because gcd(g,6)=1. Content gives
F_v F_w^k F_u=id. Therefore F_w^k=F_v^-1 F_u^-1=(F_u F_v)^-1.
The map F_u F_v=F_vu commutes with F_w modulo g. Subtract the two affine
composition formulas and clear their unit power-of-two denominators to
obtain C(w,vu)=0 modulo g. No source, entry word, or prefix clock was reset.

For a **primitive** word w and a nonempty proper prefix v, C(w,v) is nonzero.
Indeed C=0 means the affine maps commute. Faithfulness then says wv=vw as
words. The elementary commuting-word lemma (strip the shorter initial word
and induct) makes both words powers of a common word. This contradicts
primitivity and0<|v|<|w|. Thus

    g(w^k v) <= |C(w,v)|     for every k>=0.          (11)

At complete period boundaries the complementary exact formula is

    g(w^k)=g(w) (P_w^k-Q_w^k)/(P_w-Q_w), k>=1.       (12)

It follows by applying the same geometric sum to both L and B. Deep content
outside the boundaries is confined to finitely many explicitly computable
prime powers, while the boundaries retain the actual repeated-controller
anchor. For example w=(1,2) has C(w,(1))=-4, so g(w^k(1))=1 for all k.
The nonprimitive word(1,1) has C(w,(1))=0: primitivity cannot be dropped.

For completeness (12) also gives an explicit all-length charge rate for
every fixed controller and prime ell>=5. If ell divides L_w, boundary charge
at repeat j is v_ell(g(w))+v_ell(j). Otherwise set
t=ord_ell(P_w/Q_w), c=v_ell(P_w^t-Q_w^t); charge is0 unless t divides j,
and then is c+v_ell(j/t). Since v_ell(k!)<=k/(ell-1), summing these terms
and adding at most k sum_(proper v) v_ell(C(w,v)) bounds the total charge
through k periods linearly in k, with a computable controller-dependent
constant. This is not the uniform sharp5/8 constant for changing words.

**Resonance and the changing-controller boundary.** In(10), C=0 exactly
when vu and w commute as words (including empty vu). For example u=(2),
w=(1,2), v=(1) gives u w^k v=(2,1)^(k+1): an entry phase can move the true
period boundary. In the nonresonant case the bound depends on the full
retained entry u. A uniform estimate on that entry-dependent valuation is
still required before chaining arbitrarily changing controllers. Replacing
u by the empty word would assert a different, unjustified bound.

**Exact one-phase classification.** Fix any entry u and primitive controller
w of length d. As v ranges over its prefixes of lengths0,...,d-1, at most
one phase can have unbounded content g(u w^k v) as k grows. Such a phase
exists exactly when

    vu=w^r for some integer r>=0.                    (15)

At every other phase the fixed nonzero C(w,vu) bounds the full integer
content. At the exceptional phase write w=vt. Then for all k>=1,

    u w^k v=(tv)^(k+r),                              (16)

so (12) gives the exact content there and it is unbounded with k.

**Proof.** By faithfulness and the commuting-word lemma, C=0 iff vu is a
power of the primitive word w. Its length fixes |v| modulo d, leaving at
most one candidate in the stated range. If r>=1, cancelling v from
vu=w^r gives u=t w^(r-1), which proves(16); r=0 is just u=v=empty.
The geometric sum in(12) tends to infinity, since both P and Q are positive
and at least2. Every nonexceptional phase is bounded by(10). This proves
both directions, not merely a necessary resonance test.

Thus a nonresonant entry can make **all** phases uniformly bounded for one
fixed controller; an exact suffix match instead moves one true period
boundary. These are statements about the family of finite valuation words,
not existence of an infinite positive-integer orbit following that controller.

The [periodic checker](periodic_controls.py) tests complete integer content,
not just its support, against independent direct coefficient iteration.
The existing two-outward-swap mechanical-word content bound remains stronger
on its own specialized family; it is not replaced by this different
periodic/marked-entry quantifier pattern.

## 5. Exact future-return observers cannot quotient away affine states

The incoming Fourier result
[F1–F2](../scaled_count_fft_carriers_20261030/THEORY.md) distinguishes a
closed quotient for one observable from the extra data required by future
named-coordinate queries. Here is the corresponding precise Collatz claim.

**PROVED.** For q=5^E, the positive letters1 and2 generate the whole affine
group AGL_1(Z/q), of size q phi(q)=4*5^(2E-1). Any deterministic state
representation that answers, for every possible future word, whether the
resulting affine map is identity modulo q must distinguish every one of
these states. Any linear weighted representation of this exact indicator
over a field, with fixed letter operators, needs dimension at least q phi(q).

**Proof.** Let T_1(x)=(3x+1)/2, T_2(x)=(3x+1)/4. Their generated group
contains D=T_1 T_2^-1, where D(x)=2x. The two maps DT_1 and T_1D both have
multiplier3 and translations1 and1/2. Their quotient is translation by1/2,
which generates all translations modulo q. The element2 generates every
unit modulo5^E: its order is4 modulo5 and v_5(2^4-1)=1, so LTE gives order
4*5^(E-1). Hence D and the translations generate the whole affine group.
Because the group is finite, inverse generators are positive powers of
generators; positive words in1,2 already reach every state.

For states s!=t, a positive future word for s^-1 sends s to identity and
does not send t to identity. Thus merging them breaks a permitted future
query. For the linear statement, choose one positive prefix for each group
element and one future word for each inverse. The matrix of acceptance
values is the identity matrix of order q phi(q). A d-dimensional linear
representation factors this matrix through dimension d, so d is at least
that order. This is the elementary distinguishing-suffix/Hankel-rank
argument, not a novel general automata theorem.

This has exactly the same factor25 from level E to E+1 as the affine-state
lift. It explains why a lossless return-recognizer cannot replace the pair
of coefficient digits by one digit with no compensating side information.
For example, retained word length and integer clock can already determine
the multiplier; the theorem does not require duplicating that information.
It does **not** prove that an inequality
certificate must list every state, or rule out symbolic potentials or
one-sided overestimates. The pair itself takes only O(E) bits in modular
arithmetic. Also, these are formal valuation-word futures: they are not
all simultaneously legal from one fixed positive integer. A source-specific
controller can exploit its further restrictions, but must state them.

[observer_controls.py](observer_controls.py) independently checks full
reachability and positive inverse words through E=3, plus all distinguishing
pairs at E=1. The all-E result is the algebraic proof above.
