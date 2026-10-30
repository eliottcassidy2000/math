# Terminal fees with a two-layer epoch clock

**PROVED relative to the independently checked certificates and inherited E10:**
all nine eleventh-layer excursion-price families hold. For every finite positive
valuation word w with its original integer clock A,

    S11(w) <= (12207031/19531250) A,                     (E11)

and this coefficient is sharp. If w has any non-one letter,

    S11(w) <= (12207031/19531250) A-5878847/19531250.     (E11+)

The independent E9 F2 certificate checks2,509,000,000 edges,630 entries,
valid-domain closure and all positive-letter tails. **PROVED reduction:**
the next layer is equivalent to seven native price families in section5.
**OPEN:** those seven E12 families, the all-layer inequality, and universal
Collatz. No independent human/subagent audit is claimed; the producer and
certificate checker are separate algorithms, with a same-agent proof review.

The closest proved mechanism is the
[preceding paid first-return lift](../return_kernel_core_transplants_20261030/COLLATZ.md),
which proves the sharp E10 actual-content bound and reduces E11 to nine
native excursion-price families. The corrected near miss is the
[Session113 resonant envelope obstruction](../session113_symmetry_consumers_20261030/uniform_content/THEORY.md).
Its least-used relevant sidecar is the F_r delayed high-return enclosure in
section5 of that same file. We implement r=2 rather than retune false weights.

## 1. The endpoint certificate would settle all nine families at once

Retain the preceding conventions: F_w(X)=m_w X+t_w=(P_w X+B_w)/Q_w,
Q_w=2^A, A=sum of the positive valuation letters,
c(w)=v5(gcd(P_w-Q_w,B_w)), and S_e is the sum of min(c,e) over nonempty
prefixes. Set z=t-m+1 and lambda_e=(5/8)(1-5^(-e)).

**Terminal-fee lemma.** Suppose every word containing a non-one letter and
ending with

    m=1 mod25, z=0 mod5^9

satisfies S9(w)<=lambda11*A-eta, with eta>=2. Then every non-one primitive
10-excursion has delta10>=1. Hence all nine next price families hold and
the sharp E11 all-word inequality follows from the already proved E10 bound.

**Proof.** A primitive10-excursion has no proper return modulo5^10, so
S10(b)=S9(b)+1 exactly: the only tenth-layer contribution is at its terminal
prefix. Its actual endpoint satisfies the displayed (weaker) residue
conditions. Therefore delta10=lambda11*A-S10>=eta-1>=1. Apply the one-fee
lifting theorem of the preceding packet. Every family with b non-one
excursions then has total price at least b, which is at least1-a/5.
No source, block clock or tangent boundary is changed. QED.

The margin2 is a sufficient target; failure of this particular certificate
does not refute the nine necessary-and-sufficient native price conditions.
The certificate in section5 now supplies eta>2, so this implication has been
consumed, not left conditional.

## 2. A failed refinement, and why a richer observer is not enough

First retain m modulo25 but keep the old F1 epoch charges. At E9 and slope
lambda11 the independently checked terminal potential remains exactly

    V(1,0)/D = -15341763/9765625,
    D=390625000.

This has magnitude approximately1.571, below2, exactly as in the smaller
multiplier-mod5 observer at the same slope. The extra coordinate alone has
not improved the payment. The certificate still upper-bounds actual content;
it has not found an actual counterexample to the desired fee.

This negative control checks2,625,000,000 transition inequalities and630
entry inequalities, all positive-letter tails, finite values, nonpositive
maximum and domain closure. See[multiplier9_audit.json](multiplier9_audit.json).
The unchanged F1 charges can pay high-layer content before the multiplier
has actually returned modulo25. The next refinement changes the consumer of
the retained coordinate, not merely its storage.

## 3. The exact two-layer epoch enclosure

An epoch starts with a letter a>1, followed by L ones. Let m be the multiplier
at its start and h=min(v5(z),E). Ones preserve h and multiply m by3/2.
If m mod5 is not1 or4, the epoch has no return at any positive layer. Otherwise
for e<=min(2,h) let d_e be the unique residue0<=d_e<T_e=2*5^(e-1) such that

    m*(3/2)^d_e=1 mod5^e.

Its exact low-layer counts are n_e=max(0,1+floor((L-d_e)/T_e)). For e>2,
the first return offset reduces to d_2 modulo10, hence is at least d_2.
The number of returns at that level is at most
1+(L-d_2)/T_e when L>=d_2, and zero otherwise. Summing proves

    actual epoch S_E <= F2(m,h,L)
      = sum_(e<=min(2,h)) n_e
        + [h>2,L>=d_2] ((h-2)+(lambda_h-lambda2)(L-d_2)).   (C1)

This is the inherited F_r formula specialized to2, not a new assertion that
all its epochs sum below budget. It is valid for every finite epoch length.

An online implementation needs m mod25, z mod5^E, and one pending flag:

* Count levels1 and2 exactly at every letter.
* On a new epoch of height>2 in the active multiplier subgroup, keep the flag
  pending until the first multiplier return modulo25. At that first return,
  charge h-2 and clear the flag. A return at the epoch's first letter counts.
* On later ones in the same active epoch charge the tail density
  lambda_h-lambda2, in addition to the exact low-layer returns.
* A new letter>1 resets the epoch logic. Inactive or height<=2 epochs need
  no pending flag and no tail density.

These online charges sum exactly to(C1). The pending states have m in
<3/2> modulo25 other than1, and h>=3. There are20*5^E ordinary valid states
and9*5^(E-3) pending valid states. The dense serialization reserves29*5^E
slots; unused pending slots are explicitly excluded.

The initial all-ones prefix is paid separately in its original clock:
S_E(1^k)=sum_(j=1)^E floor(k/T_j)<=lambda_E*k. Its actual multiplier modulo25
runs through all ten elements of<3/2>. All ten entry phases are retained;
seeding only the two phases modulo5 would be unsound after refinement.

## 4. Proof-producing computation, with a separate checker

For E<=9 set q=5^E, D=200q, N=125q-5, so N/D=lambda_(E+2).
The producer maximizes enclosed charge minus(N/D) times original clock.
The independent checker verifies the vector V satisfies

    V(target) >= V(source)+gross_charge*D-N*a,
    V(entry) >= gross_entry_charge*D-N*a,
    max V <= 0.

It uses successive modular halving, a separately enumerated multiplier orbit,
and a sum of individual tail densities instead of the producer's inverse
powers and closed density formula. Every finite transition and every initial
phase is checked; no optimality claim about the producer is required.

Every gross charge is at most E. Consequently checking letters1..64 and

    max V+D*E-65N <= min V,
    D*E-65N <= min V

covers all omitted positive letters and entry letters. This is an infinite-
word-length proof from a finite certificate, not a bounded path search.
After any non-one prefix, induction bounds its enclosed excess by its
actual terminal V/D. Initial all-ones payment at lambda_E is at most its
payment at lambda_(E+2). At endpoint m=1 mod25 the flag is never pending.
Thus the terminal state V(1,0) proves the hypothesis in section1 if its
negative margin reaches2 at E9.

The E6 control has V(1,0)=-10546520, D=3125000 and independently passes
20,072,000 edge checks; its margin is3.3748864. This does not extrapolate toE9.
The old F1 extra-observer control itself had a sufficient margin atE6 but
an insufficient one atE9, so monotone extrapolation is explicitly disallowed.

Reproduce the full target run and its independent audit:

    c++ -O3 -std=c++17 epoch2.cpp -o /tmp/epoch2
    /tmp/epoch2 9 /tmp/epoch2-9.bin
    c++ -O3 -std=c++17 audit_epoch2.cpp -o /tmp/audit-epoch2
    /tmp/audit-epoch2 9 /tmp/epoch2-9.bin

Run those commands in this packet's directory. The binary vectors are
reproducible build artifacts, not committed files. The auditors reject
truncated, malformed, excluded-slot or violated-edge data. Finite semantic
controls in[epoch_semantics.py](epoch_semantics.py) compare literal integer
affine coefficients, the closed F1/F2 epoch formulas, and the online F2
charges on21,844 words and the514-letter resonant hostile.

**Scope boundary:** even an E11 result is a fixed-depth content theorem. The
all-layer bound and the additional source-specific convergence obligations
remain separate; neither a finite certificate nor the prospective all-layer
content inequality alone is called a proof of universal Collatz here.

## 5. Audited E11 closure, positive slack, and seven next families

The complete E9 F2 vector passes the independent checker with

    valid states = 39,203,125,
    excluded pending slots = 17,437,500,
    finite edges = 2,509,000,000, entries = 630,
    min V = -5419921740, max V = -488281240,
    V(1,0) = -898826940, D = 390625000.

Thus

    eta = 44941347/19531250 = 2.3009969664,
    beta = -max(V)/D = 2 lambda11,
    gamma = eta-2 = 5878847/19531250 = 0.3009969664.

The vector's453,125,000-byte SHA256 is
`dda483a908a059526c385c70b6b1af80bee9c1acedfa4c7197f390a465a263a8`.
See[epoch2_9_audit.json](epoch2_9_audit.json),
[binary_receipts.json](binary_receipts.json) and
[corruption_control.json](corruption_control.json). A corrupted E6 terminal
value is rejected by an edge inequality. Assertions are not used for checks.

The terminal lemma gives the uniform bound on all25 native tangent classes

    g10(X) >= eta-1 = 25410097/19531250 > 1.             (C2)

An empty class has infinite price and also satisfies(C2). Therefore all nine
families from the previous packet hold; the paid boundary potential proves
(E11) for arbitrary changes of native excursion controller. Its all-ones
extremizer has clock19,531,250 and charge12,207,031, establishing sharpness.

### The unfinished suffix also has to be paid

For a word cut at actual returns modulo5^10, each complete non-one excursion
has fee at least1+gamma. The boundary potential consumes at most1, leaving
gamma. All-ones excursions pay their exact1/5 fees.

The remaining suffix has no return modulo5^10, hence no tenth- or eleventh-
layer charge: its S11 equals its S9. If it contains a non-one letter, the
audited global bound max V=-2 lambda11 D leaves slack beta=2 lambda11.
If it consists only of ones, its cost is nonnegative by the initial run bound.
Consequently every word with a non-one has slack at least
min(gamma,beta)=gamma, proving(E11+). This suffix argument uses the new
certificate, rather than incorrectly recycling the weaker E10 global margin.
The actual entering map is identity modulo5^10, so all lower-layer counts
on these blocks agree with their native counts; the next tangent is retained
by the boundary potential throughout.

### The exact next target is seven families

For a primitive11-excursion let
delta11=lambda12*A-S11 and let g11(X) be its minimum over non-one words
with next tangent X. Equation(E11+) gives g11(X)>=gamma+Delta*A>gamma,
where Delta=lambda12-lambda11=1/(2*5^11)>0. The all-ones excursion still
has price1/5 and increment v=(4,4).

The exact kernel criterion inherited from the preceding packet reduces E12 to

    sum_(i=1)^b g11(X_i) >= 1-a/5,
    sum_i X_i = -a(4,4),

for exactly these seven pairs:

    b=1: a=0,1,2,3;
    b=2: a=0,1;
    b=3: a=0.                                         (C3)

Indeed4gamma>1, 3gamma+1/5>1, 2gamma+2/5>1 and
gamma+4/5>1. Every omitted mixed sequence is already paid; pure all-ones
zero sums need at least five terms and cost at least1. Conversely each
remaining family is part of the necessary native zero-sum criterion.
These seven families remain **OPEN**. The finite boundary does not bound
the lengths of the native words defining its edge prices.

[lift_receipt.py](lift_receipt.py) separately checks all650 boundary-potential
inequalities, the nine closed families, the new slack constants and sharp
all-ones clock, and the seven-family reduction using exact fractions.

The next efficient target is a terminal-conditioned refinement that raises
the native11-excursion fee, rather than a larger unconditional clock margin.
The observer-only negative control shows why a new digit must change the
charge or exclude an actual extremal continuation. A three-layer epoch
enclosure or a certified correction to the existing potential are candidates;
neither has been proved feasible at the next precision here.
