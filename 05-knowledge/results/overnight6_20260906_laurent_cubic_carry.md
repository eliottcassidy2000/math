# A four-channel lower-carry family by quotient characteristic positivity

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [raw-source/resultant audit](overnight6_20260906_laurent_cubic_carry_audit.md)
and [independent polynomial referee](overnight6_20260906_laurent_cubic_carry_referee.md)
accept the complete proof. No theorem ID is claimed.

The proved family statement is: for every integer `g>=11` with
`gcd(g,21)=1`, and every nonzero complex coefficient triple, a Laurent
polynomial on support

```text
(-21, 2g-21, 3g-21)
```

has first nonzero constant-term moment exactly `g` or `2g`. Both alternatives
occur on every support. At each of the three first-row cancellation phases,
the doubled moment divided by the square of the first anchor monomial is
strictly negative. The complete doubled row has eight channels and a lower
carry. The family has unbounded total width `3g`; this is a fixed endpoint-21
family, not a statement about every support with smaller endpoint at most 21.

## 1. Inheritance, corrected near miss and scope

The immediate source is the explicitly open cubic stopping object in
[the carry-corrected trace note](trinomial_trace_sign_empty_core_certificate_sep06.md),
Section 5. Its preceding
[endpoint-15 family](trinomial_width15_empty_core_returns_sep06.md)
closes `(-15,2g-15,3g-15)` by quadratic trace/norm signs. The present first
row is cubic, so a trace and norm alone do not suffice. The additional
second elementary symmetric coefficient is the missing sign coordinate.

The closest canonical tools are
[THM-4436, complete factorial-row roots](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md),
[THM-4432, two first channels with carries](../../01-canon/theorems/THM-4432-trinomial-two-channel-two-rung-noncancellation-with-carries.md),
and [THM-4440, signed duplication SOS](../../01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md).
All are currently proved. None is required for the algebraic sign proof
below, since the needed cubic real-rootedness is checked directly.

The canonical hostile `(-15,1,9)` has positive canonical second-row values
at both first roots; omitting `tau^-1` caused the corrected false claim.
Our entire family exhibits that same sign reversal. The full second row,
its lower carry, and the first anchor are retained throughout.

Exact family strings, `(21,2,3)` and the cubic-family phrase were searched
in relevant theorem and trinomial/moment notes before derivation. The
trace note explicitly labels this family OPEN. The `g=11` specialization
`(-21,1,12)` is an inherited exact finite control, and is not presented as
a new finite discovery. Every parameter already lies outside the old
smaller-endpoint-eight strip: `min(21,3g-21)>=12`. The interior charge `+1`
at `g=11` is not the positive endpoint. The support gaps are `2g,g`, so the
family is not an arithmetic progression. This is a repository coverage
claim, not an external literature-priority claim.

The live concept board is: complete affine fibers; carry normalization;
ordinary-core real-rootedness; quotient multiplication; characteristic
coefficients; parameter-shift positivity. The least-used operation is
replacing a full weighted trace-form certificate by the characteristic
polynomial of the weight operator, while retaining its real spectrum.

## 2. Complete source fibers and the actual moment identities

Write

```text
f(u)=alpha*u^-21+beta*u^(2g-21)+gamma*u^(3g-21),
alpha*beta*gamma != 0,
tau=alpha*gamma^2/beta^3,
X=alpha^(g-10)*beta^9*gamma.
```

A balanced multiplicity vector `(x,y,z)` of mass `m` satisfies
`g(2y+3z)=21m`. The gcd condition forces `g|m`. At masses `g` and `2g`
the complete nonnegative fibers are, respectively,

```text
(g-10+j, 9-3j, 1+2j),       j=0,1,2,3;
(2g-21+j, 21-3j, 2j),       j=0,...,7.
```

All negative-charge counts are strictly positive when `g>=11`.
The first support return is `g`, and no other return mass lies strictly
between `g` and `2g`. The second anchor is twice the first minus the
primitive step `(1,-3,2)`. This is the lower carry; the upper carry is zero.

Using falling factorials, define

```text
P(t)=72t^3+504(g-7)t^2+84(g-7)(g-8)t+(g-7)(g-8)(g-9),
L=(g)_7/9!,                       K=(2g)_14,
Qbar(t)=sum_(j=0)^7 (2g-14)_(7-j)*t^j / [(21-3j)!(2j)!].
```

The literal moment identities are

```text
CT(f^g)=X L P(tau),
CT(f^(2g))=X^2 K tau^-1 Qbar(tau).                 (1)
```

Here `L,K>0` on the stated domain. Every coefficient is retained, including
the second row's extra lower-carry channel.

For `p=P/72`, its discriminant is

```text
disc(p)=(g-8)(g-7)^2 D(g)/1728,
D(g)=74863g^3-1610102g^2+11532575g-27511000,
D(s+11)=74863s^3+860377s^2+3285600s+4167636.
```

Thus `p` has three distinct real roots for every real `g>=11`.
All four coefficients of `P` are positive, so none of its roots is
nonnegative. Write the three strictly negative roots as `lambda_i`.

## 3. Exact characteristic-polynomial certificate

Work in the three-dimensional real algebra `A=R[t]/(p)` at fixed `g`.
Since `p(0)>0`, let

```text
R=t^-1 Qbar mod p,
C = [[0,0,-p0], [1,0,-p1], [0,1,-p2]],
V=R0 I+R1 C+R2 C^2,
det(zI-V)=z^3+B1(g)z^2+B2(g)z+B3(g).              (2)
```

`C` is multiplication by `t` in the basis `1,t,t^2`, and `V` is
multiplication by `R`. Evaluation at the three distinct real roots of
`p` conjugates `V` to `diag(R(lambda_i))`. In particular all eigenvalues
of `V` are real. The exact coefficient definitions are

```text
B1=-Tr(V),
B2=((Tr V)^2-Tr(V^2))/2,
B3=-det(V).                                      (3)
```

The standalone producer derives these by rational polynomial division and
independently expands the determinant in `(2)`. Appendix A gives each
`B_i(s+11)` as an explicit positive-coefficient polynomial divided by a
positive integer. Their degrees are `6,12,18`, respectively. Thus
`B1,B2,B3>0` for all real `g>=11`. The characteristic polynomial in `(2)`
is strictly positive for every real `z>=0`; none of its three real roots
can be nonnegative. Therefore

```text
R(lambda_i)<0, i=1,2,3.                          (4)
```

This is the exact extra coordinate beyond the quadratic trace/norm proof.
For a real-spectrum operator, positivity of all coefficients of
`det(zI-V)` is equivalent to a strictly negative spectrum. Forgetting the
real spectrum would invalidate the implication; for example
`z^3+z^2+z+2` has all positive coefficients but a nonreal conjugate pair
with positive real parts. A trace and determinant alone also fail in
dimension three: eigenvalues `(-5,1,1)` have negative trace and negative
product, but are not all negative (`B2=-9`).

There is a small independent-identity audit interface. Give `g` and `t`
weight one. The coefficients `p_j` have degree `3-j`, so reduction modulo
`p` preserves this weighted-degree bound. Terms of `Qbar/t` with `j>=1`
have weighted degree at most six. In the sole inverse term,

```text
[t^0]Qbar / p0
 =1152(g-10)(2g-15)(2g-17)(2g-19)/21!.
```

Hence the denominator cancels before reduction and
`deg_g R_j<=6-j`. Newton sums have `deg_g Tr(C^l)<=l`, so `(3)` gives
`deg B_i<=6i`. Nineteen distinct rational parameter values therefore
suffice to certify all three Appendix A identities by exact interpolation;
these are degree-bounded identities, not an extrapolation from controls.

## 4. Detection, sign normalization and boundaries

At every zero of the first scalar moment, `(4)` and `(1)` give

```text
CT(f^(2g))/X^2 = K R(tau) < 0.
```

The quotient is real and strictly negative even if `X^2` itself is complex.
There is no equality or common-root case. Thus the first nonzero moment is
`g` when `P(tau)!=0`, and exactly `2g` otherwise. Positive coefficients
give the first alternative. Each first root is attainable, for example
by taking `beta=gamma=1` and `alpha=lambda_i`; these give the second
alternative. The worst detection order `2g` is strictly less than width
`3g`; this is not a sharp-width family.

Changing the first anchor by an integer multiple of `(1,-3,2)` multiplies
the normalized doubled value by `tau^(-2k)>0` at a negative first root.
The sign is therefore anchor invariant. In contrast the uncorrected
canonical second row satisfies `Qbar(lambda_i)=lambda_i R(lambda_i)>0`.
It has the opposite sign at every first root throughout the family.

The hypotheses are necessary for their stated roles. At `g=12`, the
support `(-21,3,15)` has first support return mass four, not twelve:
`(x,y,z)=(1,2,1)` is balanced. The displayed mass-twelve and twenty-four
identities still hold, but their first-return interpretation fails.
At `g=10`, the middle charge is negative, leaving the one-negative
source domain. No assertion here classifies these omitted cases.

The new family is not subsumed by the real-rooted ordinary-core SOS in
THM-4440. At `g=11`, `P` is proportional to
`3t^3+84t^2+42t+1`; its values at `-28,-27` are `-1175,1054`.
There is a cancellation root `tau` in `(-28,-27)`. Normalize the ordinary
cubic core as `H(s)=1+beta*s^2+s^3`, where the real cube root satisfies
`beta^3=1/tau`. Its discriminant is `-4/tau-27<0`. Thus it has nonreal
roots, and its powers do not satisfy the real-rooted hypothesis used by
THM-4440. This is a source-domain obstruction to importing that theorem,
not a counterexample to it. The inherited SOS hostiles `1+s^2` and `s^2`
also retain the real-rootedness and interior-coefficient boundaries.

The connection contract is complete factorial fibers -> weighted quotient
operator via `(1)--(2)`. It preserves every normalized first-root sign and
noncancellation. Characteristic coefficients forget the root labels, but
preserve the all-negative predicate when the real-spectrum sidecar is
retained. The lower carry and first anchor restore the scalar moment sign.
The cheapest hostile tests are `(-5,1,1)` for dropping the middle
characteristic coefficient and `(-21,1,12)` for dropping the lower carry.

## 5. Exact controls, reproduction and remaining object

The companion imports no repository mathematical producer. Eight named
primitive parameters `11,13,16,17,19,23,25,29` are checked by both the raw
charge equation and repeated Laurent multiplication retaining the gamma
degree. Both complete rows and every earlier/intermediate empty return
agree with `(1)`. Independently specialized quotient matrices give positive
signed leading Hermite minors at all eight parameters. The `g=11` minors
exactly recover the inherited trace-note control:

```text
752350432547692,
21236578848830718128306804/3,
942083106929885721679784497035400/27.
```

The corresponding uncorrected canonical forms are positive definite in
each control. This tests the sign reversal without numerical root samples.
Normal and optimized runs pass 135 explicit gates and have identical output.

```text
python -B 04-computation/overnight6_20260906_laurent_cubic_carry.py
python -B -O 04-computation/overnight6_20260906_laurent_cubic_carry.py
```

The primary script and output have the same basename as this report.
The semantic digest is
`490ea9fbea0af011c42e4e2334e2d07653ce52ab7b2dd84db5a37ca887f43881`.
The finite controls are not the proof of the all-height signs; Appendix A
and the exact algebra in Sections 2--3 provide that certificate.

The general carry-corrected negative spectrum remains open. The direct
next family in this progression would be `(-27,2g-27,3g-27)`,
`g>=14`, `gcd(g,27)=1`: first degree four and doubled degree nine,
again with lower carry one. It needs four positive characteristic
coefficients. No computation or all-height claim for that family is made
here. More generally, the source fibers make clear that fixed first-row
degree, rather than raw endpoint height, controls the size of this
certificate; the coefficient degrees and cancellation poles remain the
sidecar needed before generalizing.

## Appendix A. Exact positive-shift certificates

In the following arrays, coefficients are ordered from highest power down
to the constant term. For each `i`, divide the polynomial in `s=g-11`
having the displayed coefficient array by `denominator` to obtain exactly
`B_i(g)` from `(2)--(3)`. Every integer displayed is strictly positive.

```text
B1 denominator = 1
coefficients = [
  21956511799903/3649353012264960000,
  267537505061917/1824676506132480000,
  258708248980799/173778714869760000,
  1838646489428677/228084563266560000,
  22399824402529547/912338253066240000,
  3031962897099883/76028187755520000,
  17098873466993/633568231296000
]
```

```text
B2 denominator = 1
coefficients = [
  50226422680944373519/31074813952297120348441568870400000000,
  113616879901235260073/1553740697614856017422078443520000000,
  10087745819021998248701/6658888704063668646094621900800000000,
  9842145967675989491399/517913565871618672474026147840000000,
  14986214019851670363073549/93224441856891361045324706611200000000,
  14620201527077878112971/15133837963781065104760504320000000,
  49280724294686604537444353/11653055232111420130665588326400000000,
  83332728736136165382389/6133186964269168489823993856000000,
  2940753348365636369347883/92484565334217620084647526400000000,
  592273353663289609850977/11204860800107134741024604160000000,
  2397230679358712154165973/40461997333720208787033292800000000,
  4841024989944210619/120422611112262526151884800000,
  3759382384123836439/301056527780656315379712000000
]
```

```text
B3 denominator = 1
coefficients = [
  141387288101294140031271809/50010707962953418096784756876988735278121222144000000000000,
  1488634097300170056737776259/8335117993825569682797459479498122546353537024000000000000,
  266013922450355569787202539243/50010707962953418096784756876988735278121222144000000000000,
  1240694938040090065142018934647/12502676990738354524196189219247183819530305536000000000000,
  76728325428250169804293930157/59044519436780895037526277304591186869092352000000000000,
  45310373799488330937284618271031/3572193425925244149770339776927766805580087296000000000000,
  11253486543073733988417312426221/117672254030478630815964133828208788889696993280000000000,
  22249516657339498579721339153641/39070865596057357888113091310147449436032204800000000000,
  6788443235038428139948592734664207/2500535398147670904839237843849436763906061107200000000000,
  6525471179554184485896786161788283/625133849536917726209809460962359190976515276800000000000,
  1317151360170792594798884278371443/40593107112786865338299315646906440972500992000000000000,
  412745403974110658559370746734101/5074138389098358167287414455863305121562624000000000000,
  42516058652693891277404272913538851/260472437307049052587420608734316329573548032000000000000,
  50520580098864153210990548449989937/195354327980286789440565456550737247180161024000000000000,
  5150406354879342457552901380932023/16279527331690565786713788045894770598346752000000000000,
  391183167664601075463650413175509/1356627277640880482226149003824564216528896000000000000,
  45808773814158683288343679607/248466534366461626781345971396440332697600000000000,
  141872898165866127059323763/1922657706407143540569939064377216860160000000000,
  10189963546467592898789/732441031012245158312357738810368327680000000
]
```


**Filing:** root integrated this audited report after `f5f0f7f75`;
portable reproduction paths are shown above. The exact producer and
transcript bytes remain pinned by the sixth manifest.
