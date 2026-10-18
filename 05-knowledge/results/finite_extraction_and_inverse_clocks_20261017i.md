# Finite extraction, winding, and clocks: two attached papers under audit

2026-10-17, clock_duality_oct15. Read both attachments end to end, including
the four appendices of the71-page manuscript; rendered equation-bearing
pages of each to check the extracted notation. **PROVED, scoped/self-audited:**
the elementary counterexample, golden critical-value circle, clock-rounding
obstruction and finite golden filter below. **FINITE-EXACT:** small controls.
Numerical complex checks are separately **VERIFIED**. No whole-paper referee
claim. Collatz, six-colour forcing and the planar Jacobian problem remain open.

Sources supplied by the user:

* `2602.09253v1.pdf`: Marinov and Veselinov, *Solvability of meromorphic
  equations in elementary functions*,9 pages, arXiv2602.09253v1.
  [Primary record](https://arxiv.org/abs/2602.09253).
* `finite-extraction.pdf`: *A finite construction of matrix multiplication
  algorithms with exponent approaching9/4*,71 pages, user-supplied manuscript.
  Its theorem statements are manuscript claims unless a component is checked
  directly here or already established in this repo. No publication status
  or external acceptance is inferred from the PDF.

## 1. An explicit counterexample to the first paper's Theorem1.1

That theorem permits any connected open X in C and claims that an elementary
meromorphic f on X with infinitely many critical values has no elementary
inverse. The allowed domain is decisive.

Remove the downward vertical ray from every integer multiple of pi:

    X=C \ union_(k in Z){k*pi-i*t:t>=0}.

The removed set is closed and locally finite. X is open and connected:
each lower strip connects vertically to the upper half-plane. For z in X,
1-exp(2iz) never belongs to the nonpositive real axis. Indeed such membership
would require Re(z)=k*pi and Im(z)<=0, precisely the removed rays.
Consequently the following is one single-valued elementary holomorphic
function on X, using the principal logarithm:

    f(z)=-iz+Log(1-exp(2iz))+i*pi/2-log2.                     (1)

It satisfies exp(f(z))=sin z, and differentiation gives f'(z)=cot z.
For z_k=pi/2+k*pi, all lying in X,

    f'(z_k)=0,  f''(z_k)=-1,  f(z_k)=-k*pi*i.                (2)

Thus the critical values are infinite and distinct. Nevertheless the local
inverse at z=pi/6, w=-log2 is the elementary germ

    g(w)=-i Log(i exp(w)+sqrt(1-exp(2w))).                    (3)

Choose the square-root and logarithm branches giving g(-log2)=pi/6.
The algebraic identity sin(g(w))=exp(w) gives
f(g(w))-w in2*pi*i*Z. This analytic difference vanishes at the base point,
so it vanishes nearby. The same branch-selection argument works at every
regular source point. The square root is elementary via exp((Log)/2).
Equations(1)--(3) refute Theorem1.1 **as stated on connected open domains**.
They do not claim a counterexample to a repaired theorem imposing a global
covering condition on the domain and all inverse continuations.

The missing coordinate is the winding level. Translating z by2pi preserves
sin z but changes this f by-2pi*i. An elementary formula can contain an
unbounded sheet address; finite syntactic length does not bound its history.
For planar Jacobian work, elementary inverse germs likewise do not provide
a globally single-valued polynomial inverse. No Keller implication is made.

Two separate defects in the paper's route reinforce this distinction:

* For f(z)=(exp z-1)^2, a small loop around value0 exchanges
  log(1+sqrt(w))+2pi*i*k and log(1-sqrt(w))+2pi*i*k for **every** k.
  Each critical point is simple, but the global permutation is an infinite
  product of transpositions. A finite local ramification index does not
  imply finite support, as required by the quoted Wielandt theorem.
* exp z has regular monodromy Z and nested blocks nZ. It factors through
  exp(z/n) and the n-th power map for arbitrary n. Repeated decomposition
  has no complexity decrease supplied by the mere fact that exp z is a
  finite expression. In the covering setting, a finite sequence of primitive
  cyclic quotients has finite total degree, whereas exp has infinite degree.

Theorem5.5 adds an infinite-monodromy hypothesis. We do not silently equate
the monodromy of all continued inverse germs with a permutation of fibres
inside the slit domain X; that distinction is one of the problems being
audited. The explicit refutation above is of Theorem1.1.

## 2. Golden arithmetic makes the critical-value problem visible

Let phi=(1+sqrt5)/2 and h(z)=exp z+exp(phi*z). Its critical points and values
are exactly

    z_k=phi[-log(phi)+(2k+1)pi*i],
    h(z_k)=R exp((2k+1)pi*i*phi),  R=phi^(-phi-2).             (4)

Indeed h'=0 iff exp((phi-1)z)=-1/phi. Substitution gives (4), and
h''(z_k)=(1-phi)exp(z_k) is nonzero. Irrational rotation makes the distinct
critical values dense on the circle of radius R. Hence even an elementary
entire function can have a critical-value set that is not discrete. This
directly contradicts the paper's general discreteness assertion in5.1.

The Fibonacci return clock is exact:

    phi*F_n-F_(n+1)=(-1)^(n+1)phi^(-n),
    |h(z_(k+F_n))-h(z_k)|=2R |sin(pi*phi^(-n))|.              (5)

Compare the repo's Q-matrix clock modulo9, with Q^12=-I and Q^24=I,
in [outer_packet_golden_clocks_20261017h.md](outer_packet_golden_clocks_20261017h.md).
The same Fibonacci recurrence controls a finite residue period and arbitrarily
close irrational phase returns. The map from an integer clock to these two
observations preserves recurrence, but discards the integer lift if either
observation is used alone. This is a precise analogy to our band/clock
sidecar, not a conjugacy to Collatz or a proof of inverse non-elementarity
for h.

## 3. What the finite-extraction manuscript actually retains

The read-through ledger:

| pages/sections | mechanism examined | obligation retained |
|---|---|---|
|1--15,1--3|relative schemes, pooled integer recursion, least profile F|source side, padding, bilinearity, positive gain|
|16--23,4|determinant filtration, fixed-type restriction, Fourier separator|shared first leg, sector labels, full-space maps|
|24--29,5|wide row jumps and positive Doob transform|absorption proof before rational occupations|
|29--34,6|integral circulation rounding, exact ledger, entropy loss|every state balance; tiny/zero occupations allowed|
|34--43,7--8|exactification, paid catalyst reuse, explicit rate|degree, copies, catalyst cost, strictly positive Delta|
|44--51,9|rational interpolation, finite seed, bit/precision bounds|fixed seed cost and preprocessing, not practical size|
|51--56,10--11|finite examples, certified logarithm enclosures, questions|numerical versus rational-enclosed outputs|
|56--62,A|exact binary weights, lcm denominator growth|denominator cost and row-sign alternation|
|62--65,B|composed tensor maps on all coordinates|positive-degree terms also have defined reindexings|
|65--67,C|sparse and quadratic exactification|finite support/degree and prime anisotropy|
|67--71,D/references|notation and attribution|distinct meanings of N,K,j and the source references|

The least profile is already in THM-4568,
[the-nine-fourths-growth-lemma-is-exactly-sharp](../../01-canon/theorems/THM-4568-the-nine-fourths-growth-lemma-is-exactly-sharp.md).
Reproving its4/3 growth would be redundant. The relevant additional tools
are the finite ledger, rounding cost, wider jumps, and explicit map extension.
The manuscript's enormous explicit parameters are not an experimentally
constructed competitive matrix-multiplication seed in this session.

One particularly useful algebraic step: if every surviving tensor term has
filtration indices j'_theta>=j_theta and both sides have the same total
index, equality holds in every position. Fixed-type restriction removes
the off-diagonal terms exactly. In our setting this suggests retaining a
nonnegative defect coordinate before identifying terminal states. It does
not license deletion of an arbitrary signed residual.

The separator retains a tentative label g' and the actual sector g; phases
first impose e-e'+g'-g=0, then a nonnegative exponent(g'-g)^2 isolates
agreement. This is a model for authenticated source/port matching. A
sector-dependent change of the shared first leg is permitted **only after**
the separator has supplied independent legs. That is the same logical
restriction as forbidding source-dependent changes of a common Collatz
controller after its owner label has been forgotten.

## 4. Balance-preserving rounding needs a clock lift

Lemma6.1's circulation rounding is valid: integral incidence constraints
allow each rational edge flow to be rounded within its floor/ceiling while
keeping every vertex balanced. This does not preserve an independent
elapsed-time functional.

**Clock hostile (PROVED).** On a two-edge directed cycle, assign flow1/2 to
each edge and clock1 to each edge. The flow is balanced and its total clock
is the integer1. An integer circulation within the rounding intervals is
either(0,0) or(1,1), with clocks0 or2. The augmented incidence/clock matrix
has determinant2. Thus even an integer total clock can be destroyed.

**Repair (PROVED).** Replace state v by(v,t), and each authenticated edge
of integer duration ell by(v,t)->(w,t+ell). Then the clock is the difference
of the vertex potential t. For a finite time horizon, any rounding that
preserves the full lifted boundary automatically preserves total clock,
by telescoping. To preserve each source, use a separate source label too.
If the prescribed boundary masses are not integers, the manuscript's exact
integer-boundary rounding does not apply without a further payment.

More generally, an edge functional is unchanged under every balanced-flow
perturbation iff it annihilates the incidence kernel, equivalently it is a
vertex-potential difference. The positive durations around an unlifted
cycle fail this test. The finite time lift fixes the type issue; it does not
supply a uniform bound on the necessary horizon.

This informs the new [ROOT boundary module](root_sparse_boundary_module_20261017i.md):
round or aggregate only after the injection rows and their clocks have
been included. The paper proves absorption before solving occupations;
solving a formal occupation equation cannot replace our open entry theorem.

## 5. Kathara renewal, fractional drift, and why extraction needs a cutoff

For the user-specified Kathara recurrence a_0=0,a_1=1 and
a_n=sum_(j<n)a_j, a_(t+1)=2^(t-1) counts compositions of t>=1.
At d parts the count is binomial(t-1,d-1). This bridge is inherited from
the previous clock note, not a newly discovered recurrence.

If a letter A has mass2^(-A), every nonempty word of total clock t has
mass2^(-t), so all depths together have mass1/2 at every t>=1. The clock
generating function, including the empty word, is

    1+sum_(t>=1) z^t/2 = (2-z)/(2-2z).                      (6)

It has a pole at z=1. A root-of-unity average of an unbounded series keeps
positive exponents divisible by the averaging modulus; it cannot isolate
the constant term as finite polynomial exactification does. Retain a
literal degree cutoff, a prefix-free stopping rule, or a proved tail bill.

There is nevertheless a genuine square-root drift in a **random real lift**
of the primitive kernel. Put X=x+1/3>0. With independent letters A of law
2^(-A), and a=1 or2 according to its parity,

    X'=(3/2^A)X+2^(a-A)/3,
    E sqrt(X') <= rho sqrt(X)+2sqrt3/7,
    rho=sqrt3/(2sqrt2-1)<1,
    E X'=X+4/15.                                           (7)

The first bound follows from sqrt(s+t)<=sqrt(s)+sqrt(t) and two geometric
series; the second is exact. Thus fractional moments contract while the
first moment has positive drift. The multiplier moment is
3^theta/(2^(1+theta)-1), strictly below1 for0<theta<1 and equal1 at0,1.
Its Cramer threshold1 was already recorded in
[collatz_artin_corrections_20260927_audit.md](collatz_artin_corrections_20260927_audit.md).
The primitive translation4/15 and the explicit square-root constant here
adapt that old mechanism to the current kernel.

This is a useful control, not a source-wise Collatz inequality: real random
letters have forgotten the native integer guard and the matched-clock
two-path readout. The actual energy has weight2^(d+b) and compares endpoints
with equal band. These are precisely the extra coordinates to carry into
any proposed drift-to-energy transfer.

## 6. Golden norm filters: exact compression of two zero tests

AppendixC.3 uses the anisotropy of x^2+y^2 modulo primes p=3 mod4. The
golden version is Q(x,y)=x^2+xy-y^2, the norm of x+y*phi.
If5 is a nonsquare modulo an odd prime p, then

    Q(x,y)=0 modp iff x=y=0 modp,
    (1/p) sum_(ell=0)^(p-1) exp(2pi*i*ell*Q(x,y)/p)
       =1_{x=y=0 modp}.                                   (8)

For bounded integer defects |x|,|y|<p, this tests exact simultaneous zero.
Proof: for y nonzero, x/y would solve r^2+r-1=0, whose discriminant is5;
for y=0, x=0. Fourier orthogonality gives the second identity.
The primes3,7,13 work. **11 fails:** Q(3,1)=11. The Gaussian and golden
filters therefore select different primes; the appearance of11 in both
projects is not permission to reuse the same anisotropy claim.

The Fibonacci matrix sends(x,y) to(y,x+y), negating Q; scalar doubling
multiplies Q by4. Both preserve the zero test. Over3-power rings,
v3(Q(x,y))=2min(v3(x),v3(y)), so Q=0 mod3^(2r-1) iff both coordinates
vanish mod3^r. This is an exact bridge to the previous golden-clock orbit
decomposition. It compresses a predicate, not all the information in a pair.
For ROOT it could authenticate two finite defects (terminal and clock);
it supplies no estimate for how often the defects vanish.

## 7. Reproducible controls and next questions

`finite_extraction_bridges_20261017i.py` independently solves the manuscript's
row2 occupation system, checks a valid rounded circulation, retains the
determinant2 clock hostile, counts all norm zeros modulo2,3,5,7,11,13,17,19,
and checks24 Fibonacci steps exactly. Seventy-digit complex evaluations of
(2)--(5) have residuals below10^-60; these are numerical supplements to the
proofs, not interval certificates.

```bash
python3 04-computation/experiments/finite_extraction_bridges_20261017i.py --out /tmp/paper-controls.json
```

Live questions: can the sparse ROOT injections be balanced in a finite
source-and-clock lift with a uniform energy invoice? Can a coupled contact
circuit make its colour-shift freedom inconsistent? Can finite quadratic
filters organize the joint band/terminal Gram calculation without moving
its original source? These are concrete next tests, not implications of
shared terminology or of the manuscript's matrix-multiplication headline.
