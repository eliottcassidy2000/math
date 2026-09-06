# Whole midpoint classes: an exact cut, B-1 reversal orbits, and the open sign

Status: **PROVED ANALYTICALLY / INDEPENDENTLY AUDITED** for the exact cut,
nonvacuity and reversal statements (root and independent observer PASS).
Class root geometry is relative to
the cited finite Hadamard theorem. The declared root-value banks are
**FINITE-EXACT** only. Universal class negativity and general trinomial
two-rung noncancellation remain **OPEN**. No new theorem ID.

## 1. Inheritance and precise target

The source is the proved
[virtual Hadamard/midpoint-defect construction](nc2_hadamard_transport_overnight_hexagon_sep05.md),
with its distinction between the actual doubled moment and a virtual
Hadamard square. The closest recovered operation is
[weighted midpoint deletion](nc2_weighted_midpoint_overnight_hexagon_sep05.md).
That proves an all-root nonpositive doubling sign for each individual
auxiliary factor, but its explicit cubic hostile shows that those
factor predicates do not transport through joint Hadamard doubling.

The canonical older hostile is the mixed-sign individual Fourier-pair
example in the
[contiguous-operation note](nc2_channel_contiguous_overnight_hexagon_sep05.md).
The corrected near miss is truncating the complete auxiliary support.
The least-used sidecar here is a **whole midpoint residue or crossing
edge class**, retaining every prefix and suffix in that class. It is
coarser than an individual Fourier pair and finer than the total defect.

The live board is: full path support; residue sections; undershoot at
a cut; reversal; exact coefficient maps; root-position predicates;
factorial/hypergeometric parameter relations. The source-to-target map
preserves the entire actual defect and its carries. Passing from a
class to its real-root count loses its position relative to the first
row; that missing predicate is not restored by the root theorem alone.

Take integers

```text
0<A<B, gcd(A,B)=1, p=B-A, h>=1, x>=0,
0<=r<B, 0<=z<A, y=B*h+r, m=x+y+z, C=p*y+B*x.
```

No positive-charge or first-return assumption is needed for the cut.
Use the complete Laurent rows

```text
F_s(t)=sum_i binom(m,s+A*i)t^i,       0<=s<A,
G_(U,V)(t)=sum_i binom(U+V-A*i,V+p*i)t^i,
G=G_(y,x),     G_2=G_(2y,2x),
H=F_z^2,       P=F_z star G,
V=H star G^2,
Q_raw=sum_j (2m)! t^j/[(2x+p*j)!(2y-B*j)!(2z+A*j)!].       (1)
```

Every sequence is zero unless its underlying path counts are
nonnegative; all its admissible indices, including negative indices
of G, are retained. `star` means ordinary coefficientwise
multiplication, with no binomial normalization. The first P has
support `0..h` and is the complete factorial row. Only
`-1<=j<=2h+1` can occur in `Q_raw`; unavailable carry endpoints are
zero, not inserted. Put `D_raw=Q_raw-V`.

## 2. The exact residue and edge-crossing decomposition

There are three kinds of labelled classes.

For an alpha midpoint residue `s in {0,...,A-1}\{z}`, define

```text
R_s(j)=[t^j]G_2 * sum_(Y congruent s mod A)
              binom(m,Y)binom(m,2z+A*j-Y).                 (2)
```

The binomial convention makes this a finite sum. It counts pairs of
full paths whose alpha component hits a midpoint with the wrong
residue, while the beta component is unrestricted.

For a crossing edge `e=(e_X,e_Y)` equal to E=(1,0) or N=(0,1), put
`s_e=p` or B, respectively. For every `d=1,...,s_e-1`, define

```text
K_(e,d)(j)=H_j * sum binom(i+k,i)
             *binom(X+Y-1-i-k, X-e_X-i),                  (3)
X=2y-B*j, Y=2x+p*j,
sum over i,k>=0 with p*i+B*k=C-d,
                    i<=X-e_X, k<=Y-e_Y.
```

This counts pairs where alpha hits its selected midpoint, while beta
crosses the cut level C through an edge of type e whose starting
level is `C-d`. It retains the **entire** crossing class, not a single
prefix position. Then the exact Laurent identity is

```text
D_raw = sum_(s!=z) R_s
        +sum_(d=1)^(p-1) K_(E,d)
        +sum_(d=1)^(B-1) K_(N,d).                         (4)
```

Proof: classify a pair of unrestricted paths by whether its alpha
path hits the correct midpoint residue. If not, it belongs to exactly
one class (2). Otherwise consider its beta path. Its level functional
`pX+BY` increases strictly along every step. Because `gcd(p,B)=1`,
every integer vertex on level C is one of
`(y-B*l,x+p*l)`. Thus failure to hit a selected midpoint is exactly
crossing the level on a unique edge. Its type and undershoot give
exactly one class (3). These cases are disjoint and exhaustive.
Splitting at the midpoint or crossing edge gives the displayed
binomial products. This proves (4), including all carry terms.

## 3. Every class factors, is nonzero, and has simple negative nonzero roots

For (2), let

```text
s'=(2z-s) mod A,       ell=(s+s'-2z)/A.
```

The full row factorization is

```text
R_s=(t^ell F_s F_(s')) star G_2.                          (5)
```

For a crossing delay d, there is a unique `u in {1,...,B-1}` with
`p*u=d mod B`. Set `v=(p*u-d)/B`; then `0<=v<p`. Define

```text
U_1=y-u,                 V_1=x+v,
U_2=y-B+u-e_X,           V_2=x+p-v-e_Y.
```

All four coordinates are nonnegative: `y>=B`, `u>=1`, and `v<=p-1`.
The crossing factorization is

```text
K_(e,d)=H star [t*G_(U_1,V_1)*G_(U_2,V_2)].               (6)
```

Indeed a prefix at index l and a suffix at index k have combined
endpoint equal to the full endpoint minus the crossing edge when
`j=l+k+1`. Their levels are respectively `C-d` and
`C-(s_e-d)`. This proves the factor t in (6); it cannot be dropped.

Every factor in (5)--(6) is a complete binomial-path row, hence is
real-rooted after its full Laurent shift. Products and monomial
shifts retain this property. Clear the negative powers on the two
Hadamard inputs using the **same** shift. Finite ordinary Hadamard
composition then gives real roots, simple away from zero. The exact
imported statement is
[Brändén, arXiv:math/0403364v2, Theorem 3.8 and Section 7, printed pages 6 and 18](https://arxiv.org/pdf/math/0403364#page=18):
factorial Schur composition preserves real roots when one input has
roots of one sign; dividing its kth coefficient by k! gives ordinary
Hadamard composition and simple roots away from zero. Both statements
were checked directly in the primary paper. No infinite-matrix
Hadamard closure is assumed.

No class is vacuous: **every class has a strictly positive coefficient
at raw index j=1**. For (3), the prefix `(y-u,x+v)` and suffix
`(y-B+u-e_X,x+p-v-e_Y)` just constructed are actual nonnegative path
endpoints. Also `H_1=2 binom(m,z)binom(m,z+A)>0`.

For (2), `s+s'=2z+A*ell`, where `ell in {-1,0,1}`. Split `2z+A`
between the two residue classes by adding a total of `1-ell`
multiples of A. If ell=1, no addition is needed. If ell=0, add A
to a residue at most z. If ell=-1, add A to both residues; each
original residue is then at most `2z-A`, so the result is at most
2z. These choices are nonnegative and at most `m>=B+z> A+z,2z`.
The unrestricted beta endpoint at j=1 is also nonnegative. Thus the
claimed coefficient is positive.

Consequently, after removing its exact zero factor, each class is a
positive constant or a polynomial with all roots simple and strictly
negative. This is an unbounded **class root-geometry** result. It
does not assert a sign at roots of the different polynomial P.

## 4. Exactly B-1 reversal orbits

Reversing each complete path exchanges its two pieces. Coefficientwise,

```text
R_s=R_((2z-s) mod A),
K_(e,d)=K_(e,s_e-d).                                     (7)
```

There are `(A-1)+(p-1)+(B-1)=2B-3` labelled classes. Exactly one of
A,p,B is even, since A and B are coprime and p=B-A. Thus the reversal
involution has exactly one fixed label: the half-residue if A is even,
the half-E-undershoot if p is even, or the half-N-undershoot if B is
even. Therefore there are exactly **B-1 reversal orbits**. In (4),
the fixed class occurs once and every other orbit occurs twice.

This counts label orbits; additional accidental polynomial equalities
are not ruled out. The number B-1 is independent of h and x. It gives
a finite list of whole-class responses for the signed target, but
does not make the degree-h root comparison finite uniformly in h.

## 5. Two exact boundaries that prohibit losing the extra data

The primitive assumption in (4) is necessary unless another family
of midpoint-coset classes is inserted. At

```text
(A,B,h,r,z,x)=(2,4,1,0,0,0),
D_raw=396t+32t^2,
the classes in (2)--(3) sum to 288t+32t^2.                (8)
```

The missing `108t` comes from the beta midpoint `(2,1)`, on level
`2X+4Y=8` but not in the selected period `(-4,2)`. This is a
wrong midpoint coset, not an edge crossing. The repaired nonprimitive
description must retain such additional midpoint classes. This
algebraic hostile is not a primitive first-return claim.

A second obstruction challenges a plausible mixed-square argument.
Let

```text
H(s)=(s+10)^4(s-1)(s-8),
K(s)=(s+10)^4(s-7)(s-100).
```

The positive-root pairs strictly interlace, `1<7<8<100`, and H,K
share the identical negative binomial factor. Nevertheless

```text
[s^2]H=-21200,       [s^2]K=2000,
[s^4](H*K)=1454400000>0.                                (9)
```

Thus common negative roots, interlacing positive roots, and opposite
middle-coefficient signs do not imply mixed signed duplication.
This is an abstract proof-route hostile, not a binomial cut-class
counterexample. The exact prefix/suffix parameter relation in (6)
must be retained in any prospective polarization argument.

## 6. Exact sign banks and the structural stopping point

The companion producer reconstructs (2)--(3) directly from path-count
sums and independently verifies the full-factor convolutions (5)--(6),
the actual multinomial identity, reversal and j=1 nonvacuity. It
imports no repository mathematical implementation.

The first sign bank is the full eligible head with `A<B<=6`,
`h in {1,2,3}`, every r,z, and `x in {1,2,3}`: **638 rows, 4,392
labelled classes, 2,515 reversal orbits**. Eligibility means

```text
gcd(A,B)=1, A*x-p*z>0,
gcd(A*(B*h+r)+B*z, m)=1.
```

It identifies primitive positive-charge first fibres with support
`(-a,mA-a,mB-a)`, `a=A*(B*h+r)+B*z`, and first mass m. The second
bank consists of 300 raw proposals from seed4440 with B=3..12,
h=1..8, all residue ranges and x drawn from
`{1,2,3,5,7,11,23,101,997}`. The same filter retains **140 indexed
rows, 1,628 labelled classes, 884 reversal orbits**. Across the two
banks there are **776 distinct parameter tuples**, not 778.

At every first root in both banks, **each whole class has a strictly
negative raw Laurent value**. The producer isolates all roots of P
exactly, evaluates each class remainder on rational intervals and
refines until its strict sign is certified. Reversal partners are
compared as complete coefficient arrays; one certificate per orbit
therefore covers every labelled class. The lower carry is restored
with the sign of `rho^(-floor(2z/A))`, not discarded. No floating
root estimate, midpoint evaluation or sign-sampling heuristic is used.

These banks motivate, but do not prove, the open assertion

```text
P(rho)=0 ==> every class in (4) has value <0 at rho.        (10)
```

Together with the proved virtual sign, (10) would imply actual
two-rung noncancellation. The stopping point is exact: every class
has a complete factorization and known real-root geometry, but its
root placement relative to P is not proved uniformly. The hostiles
above exclude replacing that placement by generic compatibility or
mixed interlacing.

A second precise open route is a **faithful joint** weighted pencil
`Q_raw+(w-1)P^2`. If its complete Laurent shifts were real-rooted for
an unbounded positive sequence w, the proved square-pencil lemma
would force actual nonpositivity at first roots. The separate
auxiliary weighted networks do not establish this coupled statement.
No such joint network or unbounded preserver is supplied here.

### Incoming alpha completion: an exact overlap map

The independently audited incoming
[three-response specialization](overnight7_20260906_laurent_midpoint_transport.md)
and [all-A=2 alpha completion](overnight8_20260906_alpha_completion.md)
were read in full. They fit this cut but use a different grouping.
At `A=2,B=3,r=0,z=1,x>=1`, write their full factors as O,E,Beta,C,D,
so `F_z=O`, the missing alpha row is `t^-1 E^2`, and
`G_2=Beta^2+2tCD`. With their responses G1,G2,G3, the exact map is

```text
our R_0 = G1+G3,
our K_(N,1)=our K_(N,2)=G2/2; there are no E crossings.
```

Thus their three-response identity is recovered, with the combined
`G1+G3` treated as one full class. The general cut, its B-1 reversal
count and the simple-negative nonzero root geometry of each class
are the additional all-parameter statements here; the special
three-response decomposition is not counted again as new.

More importantly, the incoming A=2 theorem proves uniformly

```text
P(rho)=0 ==> W(rho)<0,
W = alpha_double star G^2,
alpha_double_j=binom(2m,2z+2j).
```

It completes both alpha residues while keeping beta paths that hit
their cut. Its real-rooted carrier preserves the original zero
coefficient. In the displayed A2B3 specialization `W=V+G1`, while
our alpha class is `R_0=G1+G3`, not G1. Accordingly W<0 does not
give an individual class sign in (4). The same actual row admits
the alternative grouping `Q_raw=W+R_skip`, with
`R_skip=alpha_double star (G_2-G^2)`; only the beta-skip payment
remains open in that already-proved A=2 completion. These are
complementary partitions, not contradictory status claims. The
universal **individual** signs asserted in (10) remain open.

For literature inheritance, the factorial splitting behind
[THM-4436](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md)
also matches the hypergeometric block conditions in
[Martínez-Finkelshtein–Morales–Perales, arXiv:2309.10970, Theorem 4.6, printed page 25](https://arxiv.org/pdf/2309.10970#page=25):
negative numerator parameters below `-h+1` and positive denominator
parameters give same-sign real roots. Its statement and proof were
read directly. This is an additional classical-source connection,
not an external novelty claim and not a doubled-mass identity. Its
single-parameter interlacing does not by itself compare m and 2m.

## 7. Reproduction and independent audit

```bash
python3 -B 04-computation/nc2_midpoint_classes_overnight_hexagon_sep05.py
python3 -B -O 04-computation/nc2_midpoint_classes_overnight_hexagon_sep05.py
```

Both full executions pass **144,664 optimization-live gates**, and their
882-byte outputs match each other and the
[stored output](nc2_midpoint_classes_overnight_hexagon_sep05.out) byte for
byte. The two sign banks certify respectively **8,755** and **7,378**
labelled class/root values. The additional structural corner bank has
144 indexed coprime cases: `B=2..5`, every coprime A, `h=1,2`,
`r in {0,B-1}`, `z in {0,A-1}`, and `x=0,1`, including coincident
corner labels. It independently counts the real roots of every class
after removing its exact zero power; constants are allowed.

The checker compares the full shifted and normalized suffix supports,
every suffix coefficient, and the literal edge-count sum with both
prefix/suffix convolutions. In particular it directly tests the factor
t and `j=l+k+1` displayed in (6), not merely an equivalent unlabelled
product. Both hostiles in Section 5 are literal exact coefficients.

Root and a separate observer independently read and passed the full
analytic proof, its scope and the primary finite-Hadamard input. The
observer additionally replayed the entire 144-case structural corner
bank and both hostiles. The complete 778-index sign banks are this
producer's normal/optimized finite certificates, not an independently
reimplemented universal argument. The incoming overlap map was checked
against both complete incoming notes and does not add a sign claim.

```text
source SHA256 fa3e1acdefb9da47882ccbd0b75c3c53af25e134231a3272892bb02aae1c1748
output SHA256 085ae75f15429eda7da56918c5f37b11d609214d2e4bf5cb03caf8349df067a0
semantic trace 27b31c14547cbc70393506f4fcab3b8e289acb42d1c3cfc1fe4228fbd29878f6
```

This is the recorded structural stopping point. The exact cut and
class geometry are proved; their unbounded signs at first-row roots,
the faithful joint weighted pencil and general trinomial doubling
noncancellation are not.
