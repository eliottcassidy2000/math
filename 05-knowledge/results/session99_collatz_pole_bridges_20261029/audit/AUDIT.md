# Independent audit of Session99 pole bridges

Status: PROOF REVIEW PASSED for all four lane notes listed below;
FINITE-EXACT controls passed independently. This does not certify a Collatz
convergence theorem, a canonical bijection, or an operator conjugacy.
Source pin supplied and used: `d9d7e33349653b6b65dbad477841c9ed518c081a`.
The inherited Session98 `limits/ALL_POLES_CORE_EXPANSION.md` was inspected
at this pin. The Session94 valuation-tree file was absent from this sparse
checkout; its needed elementary digit/clock facts are rederived below rather
than treated as an independently inspected source.

Reviewed notes: `structural/GEOMETRIC_SLOT_BRIDGE.md`,
`transfer/CARRY_CLOCK_POLES.md`, `poles/POLE_TRANSPORT.md`, and
`readout/POLE_ORDER_READOUT.md`.

## 1. Highest poles and cancellation

The tournament bank with defect k is indexed by partitions d of k via parts
2d+1. Its support is 2k+length(d), at most 3k, with equality uniquely at
(3,...,3). Consequently exactly one summand has order 3k+1. Its numerator
is nonzero at 4^k. Normal convergence of all other banks near that point,
proved in the inherited source by a Gaussian-in-support bound, is essential:
positivity of Taylor coefficients alone would not prove noncancellation at
all continued poles. That stronger normal-convergence argument is present.

For carry moments with J=3k+1, heights h_j=sum_{i>j}m_i are nonincreasing
and h_J=0. At a positive pole 2^-e, e>0, exactly J-1 resonant slots require
h_1=...=h_{J-1}=e. This forces m_1=r-e,m_J=e and all other marks zero.
For each shape its final vertex is a leaf. Thus the unique maximal height
profile has leading coefficient

    c_k binom(r,e) 3^((J-1)(r-e)) 2^(2k) L(2^-e)>0.

Every other profile has lower order there. This proves sharpness without
an invalid blanket claim that all rational summands have positive continued
values. At z=1, the unique top profile has all r marks at the first vertex;
its coefficient is 3^(r(J-1))c_k2^(2k+1)>0. At k=0 carry is identically 1,
so there are no shifted poles. The note correctly excludes unsupported
exact multiplicity assertions for the shifted nonreal poles.

## 2. Rationality and exact clocks

For fixed k,r there are finitely many shapes and weak mark compositions.
Multinomial expansion followed by independent valuation summation proves
rationality. Coefficients of each fixed degree use finitely many positive
valuation words. No convergence or infinite interchange is needed for this
formal derivation. Carry is the raw numerator C, not C/2^D.

The counting variable is z^D. Haar weighting substitutes z/2. A Haar-weighted
normalized-intercept rth moment instead substitutes z/2^(r+1). Pole locations
must be labeled by this convention before they are compared.

There are TWO controllers here. In the structural note an even valuation
branches, with branch probability 1/3, giving the critical full ternary law.
Write a=1+b+2G. The leaf bit b=0 has probability 2/3, the branch bit b=1
has probability 1/3, independently of P(G=g)=(3/4)4^-g. For k branches and
2k+1 leaves the minimal clock is 2k+(2k+1)=4k+1. Conditioning on a particular
shape, or mixing the shapes at fixed k, preserves the independent G law.
In the transfer note only valuations divisible by 3 branch. This probability
is sum_{j>=1}2^-3j=1/7, giving a different, subcritical controller. Its
smallest clock is 3k+(2k+1)=5k+1 and its counting poles are cube roots of
unity. These laws are compatible, but must never be presented as the same
clock or the same stopping rule.

## 3. Source, measure, and symmetry guards

For any odd multiplier b, a finite valuation word of length J and sum D has
carry C_b=sum_i b^(J-i)2^(a_1+...+a_(i-1)). The exact source cylinder is

    n=(2^D-C_b)b^-J mod 2^(D+1).

A reverse induction through the odd inverse branches proves this congruence
forces every individual valuation, not only an odd endpoint. It is one odd
residue class mod 2^(D+1), hence has normalized odd-Haar mass 2^-D. Thus finite
digits are iid geometric for every odd b, after excluding the countable
backward orbit of points where b*n+1=0 in the 2-adics. This measure statement
is not a frequency or convergence theorem for an individual positive source.

The structural slot/vertex bijection requires a fixed core, ambient label
order, and fixed ordered tree shape. It preserves only total marks and their
negative-binomial law. There is no S_(3k)-equivariant map from the conjugacy
class of 3^k to edge-labelled ordered ternary trees when k>=1: the source
stabilizer contains a nontrivial 3-cycle product, while target stabilizers
are trivial. The proof is scoped correctly and does not rule out arbitrary
non-equivariant encodings. Carry and source residues are lost under total
clock projection, as the note's two-word hostile witness demonstrates.

## 4. Strong hostile multiplier control

The odd 5x+1 map has the same Haar digits, critical-tree mark law, geometric
slot factors, and cheaper-controller counting-clock poles. Nevertheless,

    13 --a=1--> 33 --a=1--> 83 --a=5--> 13

is a positive native cycle disjoint from 1-->3-->1 (digits 1,4). The first
cycle's digits are leaves for both controllers, so each auxiliary return
terminates immediately while repeated native dynamics stays in that cycle.
This is an exact counterexample to any inference from those laws alone to
convergence into the cycle containing 1. It is not a statement that 3x+1
has a nontrivial positive cycle. The carry-profile rational algorithm also
generalizes to b by replacing its powers of 3 with powers of b, without
altering the positive pole orders when b>0 is odd.

## 5. Differential operators

For f=z^m/(1-qz)^(m+1), applying (D+1),...,(D+p) successively gives
(m+p)!/m! times z^m/(1-qz)^(m+p+1). Dilation first contributes
2^(-(p-1)m); the stated cycle insertion constant supplies the remaining
loss exponent. The multiplication by z^p occurs last. The derived tower
PDE and the entire core-series pantograph equation agree with these ratios.
The full marker hierarchy retains gcd response data and is not a closed
scalar finite-order ODE. Its projection-loss example (3,7) versus (5,5)
is valid. No dynamical identification follows from differential dilation.

## 6. Independent finite evidence

Run from the Session99 packet root (`session99_collatz_poles`):

    python3 audit/check_independent.py
    python3 -O audit/check_independent.py

Both produce identical JSON. No production checker is imported. Controls:

* 1,020 source cylinders: b in {3,5,7}, word length 1..4, digits 1..4;
  direct native execution at three distinct lifts of each residue.
* 180 coefficient comparisons: k=0..2, r=0..3, D=0..14; direct enumeration
  of chronological valuation words versus independent convolution of the
  rational mark factors (including coefficients that must vanish).
* Pole-bank maximality and leading-amplitude ratio through k=10.
* Both positive 5x+1 cycles above and their disjointness.

Finite controls support, but do not replace, the all-parameter arguments.

## 7. Readout and amplitude lane: final review

The completed readout note is approved on its stated scope. On positive odd
indices, Phi(4^k)=4^oddpart(3k+1) is a conjugacy by its explicit definition.
It does not intertwine the insertion edge k->k+1, and the note says so.
The leading-amplitude ratio cancels to
2^(6k+4)(3k+1)(3k+2). For odd k its second factor is odd, hence subtracting
6k+4 from the amplitude valuation increment recovers v2(3k+1) exactly.
For even k=2r the increment instead measures 1+v2(3r+1); the parity guard
is indispensable and correctly imposed.

Legendre's identity gives v2(D_k)=2k+s_2(k)-s_2(3k), proving the displayed
digital finite difference. This difference is along adjacent integers,
not successive Syracuse sources, so there is no asserted telescoping
potential. Increasing edge 3->5 rules out monotonicity of pole radius and
amplitude. The all-ell singleton-cycle tower has one nonzero rational term
at each distinct pole 2^((ell-1)k), of order ell*k+1, with Gaussian decay
of its coefficients ensuring normal convergence. Thus the ell=5 hostile
uses the actual same analytic construction and not merely similar numerals.

The readout checker was independently rerun both normally and with -O;
outputs are identical (`readout_rerun*.json`). This rerun supplements, rather
than replaces, the independent source/moment tests above.

## 8. Publication path and provenance review

Packet-root relative references to structural/, transfer/, poles/, readout/,
and audit/ remain valid under the intended destination
`05-knowledge/results/session99_collatz_pole_bridges_20261029`.
The readout's `../poles/POLE_TRANSPORT.md` references resolve there.
Checker commands written with bare filenames assume that lane's directory;
the independent audit commands above assume the packet root.
Session98/Session94 shorthand names are provenance labels, not claimed
packet-local links. Source pin is consistent throughout the four notes.
The readout reference to THM-4507 should include its exact canonical path,
`01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md`;
this source-hygiene correction has been sent to the publishing parent.
No external literature priority claim or imported proof is certified here.


## 9. Bounded extension: normalized-intercept moment boundary

`transfer/INTERCEPT_MOMENT_BOUNDARY.md` passes proof review. Its slope bound
U<=3/2 follows from J=3K+1,D>=5K+1. The moment fixed-point equation is the
least nonnegative solution because complete finite-height trees increase to
the almost surely finite subcritical tree; no unsupported selection of a
cubic root is made. For integer r>=1 the pure child-intercept monomials have
coefficient alpha_r=(1+m_r+m_r^2)/7. Dependence of U_i and V_i within a
single subtree is not discarded: cross monomials are bounded using U_i<=3/2
and strictly lower powers of V_i, and the exact mixed-moment formulas retain
E(UV).

Most importantly, completed-tree height truncation is not independent
replacement of absent subtrees by zero intercepts in an otherwise completed
parent. All three children must complete. The note explicitly retains that
indicator and obtains an UPPER BOUND after dropping sibling constraints.
This is sufficient for a uniform geometric bound if alpha_r<1. Under an
assumed finite y_r, the pure terms and strictly positive leaf contribution
give y_r>=alpha_r*y_r+c_r, disproving finiteness when alpha_r>=1, including
the equality case. Thus there is no circular finite-moment assumption.

The threshold test is valid: B_r decreases for r>=1 and is at most 1/21,
so L_r+B_r*x^3-x is strictly decreasing on [0,2]. Its least nonnegative
root is below 2 for r=1,2,3 and above 2 for r=4; alpha_r<1 is equivalent
to m_r<2. Divergence for all integer r>4 follows by restricting to V>=1,
not by an unjustified monotonicity of moments over all V. No noninteger
threshold is asserted. The exact second/mixed equations were checked both
against the original expansion and the independently derived recurrences
below. The mean 7/12 is correctly marked inherited.

### Independent sibling-aware check

Run `python3 audit/check_intercept_independent.py` from the packet root;
normal and optimized outputs agree exactly in `intercept_independent*.json`.
It proves rational barriers m_2<7/5 and m_3<19/10, giving alpha_2<109/175
and alpha_3<93/100; also m_4>=L_4=85536/32767>2. Strict upper barriers use
the same elementary decreasing-cubic criterion but are sharper than the
production check's barrier 2.

Let q,x,m,y,h,t denote respectively completed-height mass, EU,EU²,EV,EUV,
and EV², all including completion indicators. Independent exact recurrences
used by the checker include

    q' = 6/7+q^3/7,
    x' = L_1+x^3/21,
    m' = L_2+9m^3/511,
    y' = L_1/3+x^3/63+y(x^2+qx+q^2)/7,
    h' = L_2/3+3m^3/511+h(m^2+xm+x^2)/21,
    t' = L_2/9+m^3/511+t(m^2+qm+q^2)/7
         +2h(m^2+xm+x^2)/63+2yh(m+x+q)/7.

Here L_j has the note's meaning. Starting all six quantities at zero, four
height levels satisfy exact monotonicity and bounded-slope controls. Setting
q=x=1 in their limits recovers the stated mixed/second formulas. Omitting
sibling factors q in the finite-height y equation gives a STRICTLY larger
quantity, directly testing the dangerous false-equality shortcut. These
finite controls do not themselves prove the infinite-moment theorem; its
proof is the contraction/positive-forcing argument audited above.

The statement that raw positive-integer carry moments diverge even at fixed
K>=1 is also consistent: the root is a branch, its unrestricted valuation
appears as 2^a in a later carry summand, and Haar weight contains 2^-a.
The positive tail therefore already fails to sum for moment exponent 1,
and does so for higher positive integer exponents. This does not apply to
normalized intercepts, nor to counting formal coefficients.
