---
id: THM-3359
title: "Modular C-finite supports, harmonic density, and periodic multiplicative scars"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Modular residue
  supports of every integer C-finite
  sequence are ultimately periodic. Their harmonic logarithmic coefficient
  is the rational cycle density, while the weak and strict multiplicative
  deletion scars have collision-corrected coefficient delta-delta^2/2.
  Tournament walk sequences and regular ternary-language level counts enter
  through finite integer transfer matrices. The harmonic residue does not
  recover tournament orientation, ternary addresses, Berggren ancestry, or
  an arbitrary subset of the natural numbers.
source: codex-2026-08-14-modular-harmonic-supports
audit: >
  An independent proof reconstructed the finite-state bound, Hurwitz/digamma
  remainder, ordered-to-unordered collision asymptotic, transient boundary,
  tournament and generalized-relation typing, and tame shell-clock hypotheses.
  It required explicit all-n recurrence scope, positive index classes, the
  pole iff h>0 boundary, and separation of the Boolean density valuation from
  the nonlinear scar coefficient. Ordinary, optimized, and stored exact
  transcripts agree after those repairs.
depends_on:
  - THM-2438-poisson-newton-ternary-half-and-harmonic-divisor-incidence
  - THM-3315-tournament-cut-switching-centered-coronal-walk-compiler
  - THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction
  - THM-3356-primitive-affine-determinant-shells-parabolic-orbits-and-prime-clock-resultants
related:
  - THM-2005-support-dirichlet-automatic-tournament-atlas
  - THM-3200-fixed-lrc-channel-cleared-overlap-quasipolynomial-and-mass-recurrence-boundary
script: 04-computation/modular_cfinite_harmonic_supports_thm3359.py
output: 05-knowledge/results/modular_cfinite_harmonic_supports_thm3359.out
script_sha256: afb46cf743af9f4f1d0353c896d9d92f1306d466f367385a8a6d388c8f43a1cf
output_sha256: 66caf6683b005772a3401cc472ccaea11b3af1c1bbb6f0a7d49ead33ee06c9af
hash_basis: LF-normalized bytes
---

# THM-3359 -- modular C-finite supports and harmonic scars

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Modular recurrence-support theorem

Let `r>=1`, and let `(a_n)_(n>=0)` be an integer sequence satisfying, for
every `n>=0`,

```text
a_(n+r)=c_(r-1)a_(n+r-1)+...+c_1 a_(n+1)+c_0 a_n,       (1)
```

with integer coefficients. Fix `m>=2` and `S subset Z/mZ`. The relevant
subset of the positive harmonic-series index carrier is

```text
H=H_(m,S)(a)={n>=1 : a_n mod m belongs to S}.             (2)
```

Then there are integers `mu>=0`, `p>=1` such that

```text
mu+p<=m^r,
a_(n+p)=a_n mod m       for every n>=mu.                  (3)
```

If `gcd(c_0,m)=1`, one may take `mu=0`. Let exactly `h` positions of the
eventual `p`-cycle be accepted by `S`, and put

```text
delta=h/p.                                                (4)
```

The ratio in (4) is unchanged if the displayed period is refined. One has

```text
#(H intersection [1,x])=delta x+O(1).                    (5)
```

More precisely, choose the first persistent positive **index** representative
`n_j` of each accepted eventual index class modulo `p`, and put every earlier
accepted index in a finite set `F`. Then

```text
H=F disjoint_union union_(j=1)^h {n_j+kp:k>=0}.           (6)
```

Consequently the support Dirichlet series has the exact expression

```text
D_H(s)=sum_(n in H)n^(-s)
      =sum_(f in F)f^(-s)
       +p^(-s) sum_(j=1)^h zeta(s,n_j/p),                 (7)
```

initially for `Re(s)>1`. Its meromorphic continuation has residue

```text
Res_(s=1) D_H(s)=delta,                                   (8)
```

and it has a simple pole there exactly when `h>0`; for `h=0`, (7) is a finite
Dirichlet polynomial.

and the harmonic partial sums satisfy

```text
sum_(n<=x,n in H) 1/n
 =delta log x+C_H+O(1/x),                                 (9)

C_H=sum_(f in F)1/f
    -(h/p)log p-(1/p)sum_(j=1)^h psi(n_j/p).              (10)
```

In particular,

```text
sum_(n in H)1/n<infinity
 iff h=0
 iff H is finite.                                        (11)
```

Thus an infinite ultimately periodic subset of the natural numbers always
occupies a logarithmically divergent subseries of the harmonic series.

## 2. Proof by the finite recurrence state

Put

```text
X_n=(a_n,...,a_(n+r-1)) mod m.                            (12)
```

The recurrence gives a deterministic transition on the `m^r` states of
`(Z/mZ)^r`. Among `X_0,...,X_(m^r)` two states coincide. Determinism makes
the subsequent orbit periodic, and if the first repetition is
`X_mu=X_(mu+p)`, then `mu+p<=m^r`. The transition determinant is `+/-c_0`.
When `c_0` is a unit modulo `m`, the transition is a permutation, so the
initial state already lies on a cycle and `mu=0`.

Equation (6) follows by reading the first coordinate of this state orbit.
Counting its arithmetic progressions proves (5). For one progression,

```text
sum_(k>=0)(n_j+kp)^(-s)=p^(-s)zeta(s,n_j/p).              (13)
```

Equations (7)--(10) now follow from the residue of Hurwitz zeta at one and
the elementary digamma asymptotic. No analytic property of the recurrence
roots is required.

## 3. The collision-corrected periodic multiplicative scar

Use THM-2438's full positive multiplicative carrier. For an arbitrary
`H subset N_(>0)`, let

```text
I_H(n)=#{h in H:h divides n},                             (14)
```

and let `D_H^(<=)(n)` and `D_H^(<)(n)` be the weak and strict unordered
factor-fibre losses after deleting `H`. THM-2438 proves pointwise

```text
D_H^(<=)=I_H-K_H,
D_H^(<) =I_H-K_H-Q_H,                                    (15)
```

where `K_H` counts distinct hole--hole pairs and `Q_H` the diagonal squares.

For the ultimately periodic support (2), if `delta>0`, then

```text
1/N sum_(n<=N) I_H(n)=delta log N+O(1),                   (16)

1/N sum_(n<=N) D_H^(<=)(n)
 =1/N sum_(n<=N) D_H^(<)(n)
 =(delta-delta^2/2)log N+O(1).                           (17)
```

The coefficient in (17) is not merely the harmonic density. Indeed, with
`A_H(x)=delta x+O(1)`, the ordered hole-pair count is

```text
P_H(N)=#{(a,b) in H^2:ab<=N}
      =sum_(a<=N,a in H) A_H(N/a)
      =delta^2 N log N+O(N).                              (18)
```

Passing to unordered distinct pairs divides the leading term by two, while
the diagonal is `O(sqrt(N))`. Combining (15), (9), and (18) proves
(16)--(17). If `delta=0`, then `H` is finite and THM-2438 instead gives the
exact finite limits

```text
lim_(N->infinity) 1/N sum_(n<=N) I_H(n)
=lim_(N->infinity) 1/N sum_(n<=N) D_H^(<=)(n)
=lim_(N->infinity) 1/N sum_(n<=N) D_H^(<)(n)
=sum_(h in H)1/h.                                        (19)
```

This sharpens the divergent boundary in THM-2438: dense hole--hole collisions
change the leading scar coefficient from `delta` to
`delta-delta^2/2`.

## 4. Finite linear representations: tournaments and ternary trees

Every sequence with an integer matrix representation

```text
b_n=u^T M^n v                                             (20)
```

is C-finite by Cayley--Hamilton, so Sections 1--3 apply.

For a tournament `T`, take its ordinary adjacency matrix `A(T)`. The total
directed-walk sequence

```text
w_n(T)=1^T A(T)^n 1                                      (21)
```

is therefore C-finite. For a cut-switched tournament, THM-3315 gives the
more informative compiler

```text
det(I-x A(T^d))=P(x/2)-(x/2)N_d(x/2).                    (22)
```

The centered denominator `P` alone does not determine the sequence; the
observer numerator `N_d` is a required sidecar.

For a regular language `L subset {-1,0,1}^*`, let `M_(ij)` count letters
sending DFA state `i` to state `j`. Then

```text
b_n=#(L intersection {-1,0,1}^n)
   =e_start^T M^n 1_accept.                              (23)
```

This is the precise common mechanism between tournament walks and regular
ternary-tree level sets: both have finite integer transfer matrices. It is
not an identification of the underlying relations. A directed relation with
a missing pair, a both-way pair, or a tie is not a tournament; its transfer
sequence may still satisfy the theorem.

There are two different support carriers here.  For a unary regular language
`K subset {U}^*`, the accepted **length set** itself is ultimately periodic:
one reads the orbit of the single DFA transition.  For a ternary regular
language, (23) instead says that the **level count** `b_n` is C-finite, so a
modular predicate on `b_n` is ultimately periodic in `n`.  It says nothing of
the base-three integer set represented by the accepted addresses.

The same-language hostile makes this distinction exact.  Let

```text
L_01=1{0,1}^*                                           (23a)
```

in canonical base-three notation.  Its length-`n` count is `2^(n-1)`, so the
predicate `b_n congruent 1 mod 3` accepts exactly the odd lengths and has
harmonic coefficient `1/2`.  But the represented value set has `2^(n-1)`
points in `[3^(n-1),3^n)`, hence block harmonic mass at most
`(2/3)^(n-1)` and finite total mass.  Its gap from `(3^n-1)/2` to `3^n` is
unbounded, so it is not ultimately periodic.  This is the automatic-value
boundary also visible in THM-2005: finite-state address recognition is not a
unary time transition.

Finite Boolean combinations of modular supports remain ultimately periodic.
On a common cycle their harmonic residues are finitely additive, and in
particular

```text
delta(H xor K)=delta(H)+delta(K)-2delta(H intersection K). (24)
```

This scalar valuation knows normalized cycle occupancy, not the cardinality
`h` unless the period `p` is retained. It does not orient an arc. Moreover,
`delta` is the Boolean valuation; the nonlinear scar coefficient is not. For
example, even and odd indices each have scar coefficient `3/8`, while their
union has coefficient `1/2`, not `3/4`.

There is a useful exact Boolean-algebra formulation. The full power set
`P(N_(>0))` is a Boolean ring under intersection and XOR, and

```text
I_harm={A subset N_(>0):sum_(n in A)1/n<infinity}         (24a)
```

is its proper summable ideal. Modulo finite sets, the Boolean algebra of
ultimately periodic supports is exactly the algebra of clopen subsets of the
profinite integers: a residue class modulo `p` becomes a clopen coset of Haar
mass `1/p`. Under this identification, `delta` is normalized Haar measure and
also the residue at `s=1` of (7). Thus the harmonic residue is a genuine
measure on the periodic Boolean quotient, even though the unrenormalized
harmonic weight (24a) is infinite on every nonempty clopen class.

## 5. The complete order-four tournament test

The `2^6=64` labelled tournaments of order four have four isomorphism types:

| type | labelled count | `det(I-xA)` | total walks `w_0,w_1,...` | positive odd-index support |
|---|---:|---|---|---|
| transitive | 24 | `1` | `4,6,4,1,0,...` | `{3}` |
| source over `C_3` | 8 | `1-x^3` | `4,6,6,6,...` | empty |
| `C_3` over sink | 8 | `1-x^3` | `4,6,6,6,...` | empty |
| strong, scores `2,2,1,1` | 24 | `1-2x^3-x^4` | `4,6,8,11,16,22,30,43,...` | `n=3 mod 4` |

For the strong class,

```text
w_(n+4)=2w_(n+1)+w_n,                                   (25)
```

so modulo two the walk sequence is the pure cycle `0001`. Hence

```text
sum_(n<=N,w_n odd)1/n
 =1/4 log N+(gamma+log 2-pi/2)/4+O(1/N),                 (26)

delta-delta^2/2=1/4-1/32=7/32.                          (27)
```

For the transitive type the support is the singleton `{3}`, so the harmonic
mass and both limiting scars equal `1/3`.

The THM-3315 hostile is visible already in this switching class. Its centered
polynomial is

```text
P(z)=1+4z+12z^2+16z^3+8z^4.                             (28)
```

Equal-cardinality singleton switches can have

```text
N_tr(z)=4+12z+16z^2+8z^3,
N_str(z)=4+12z+32z^2+24z^3.                             (29)
```

They therefore have the same centered spectrum and cut size but finite versus
logarithmically divergent odd-index harmonic support. Conversely the source
and sink `C_3` types have identical total-walk sequences by transposition, so
the harmonic observer cannot recover apex direction.

## 6. Fibonacci recurrence and the three Berggren ancestry rays

The Fibonacci state modulo two is the three-cycle

```text
(F_n,F_(n+1)) : (0,1)->(1,1)->(1,0)->(0,1).              (30)
```

Thus `F_n` is odd exactly when `3` does not divide `n`, and

```text
sum_(n<=N,F_n odd)1/n
 =H_N-(1/3)H_floor(N/3)
 =2/3 log N+(2gamma+log 3)/3+O(1/N),                    (31)
```

with multiplicative-scar coefficient `4/9`.

THM-3339 gives the exact ancestry language

```text
L_F=(BA)^* union A(BA)^* union C(BC)^*.                  (32)
```

Its depth-`d` level count is one at even depth and two at odd depth, with

```text
sum_(d>=0)b_d x^d=(1+2x)/(1-x^2),        b_(d+2)=b_d.    (33)
```

Modulo two, the positive accepted depths are the even depths. Their harmonic
coefficient is `1/2` and their scar coefficient is `3/8`. Each of the three
ancestry index rays has coefficient `1/3` and scar coefficient `5/18`; each
of THM-3339's six channel-order states has coefficient `1/6` and scar
coefficient `11/72`.

These values are different because ancestry index, tree depth, and channel
time are different copies of `N`. The index carrier is mandatory data.

## 7. Determinant-shell Boolean clocks

THM-3356 supplies a further exact instance. Let `u` be primitive, `c!=0`,
choose `v` with `det(v,u)=1`, and fix a parabolic residue
`r mod |c|`. Its orbit is

```text
x_(r,t)=cv+(r+ct)u,              t>=0.                   (34)
```

Put `Q(n)=||cv+nu||^2=An^2+2hn+C`. For a tame modulus

```text
M=product_(j=1)^k p_j^(e_j),
p_j distinct odd primes, p_j=1 mod 4, e_j>=1,
gcd(M,Ac)=1,                                             (35)
```

THM-3356 makes the root set `R_M={n:Q(n)=0 mod M}` an affine
`F_2^k` torsor. Since `c` is a unit modulo `M`, the root times on orbit `r`
are exactly

```text
R_(M,r)=c^(-1)(R_M-r) mod M.                             (36)
```

Thus every parabolic orbit sees every Boolean root state once per `M` turns.
For `B subset F_2^k`, let `A_(r,B)` be the corresponding positive orbit times
`t+1`, and let `a_epsilon in {1,...,M}` represent their residues. Then

```text
D_(r,B)(s)=M^(-s)sum_(epsilon in B)zeta(s,a_epsilon/M),  (37)

Res_(s=1)D_(r,B)(s)=|B|/M.                              (38)
```

Since every coordinate of the spinor `Phi(x_(r,t))` is quadratic in `t`,

```text
(E-1)^3 Phi(x_(r,t))=0.                                 (39)
```

If `b_(r,B)` is the period-`M` clock indicator, then

```text
(E^M-1)b_(r,B)=0,
(E^M-1)^3 (b_(r,B)(t) Phi(x_(r,t)))=0.                  (40)
```

This is a typed source for the quasipolynomial recurrence grammar of
THM-3200, not a physical LRC current.

For the U-spine `u=(1,1)`, `c=1`, `v=(1,0)` and
`M=1105=5*13*17`, the eight roots are

```text
23,231,418,431,673,686,873,1081.                         (41)
```

Antipodal conjugation `n -> -1-n` gives four vertices

```text
{23,1081}, {231,873}, {418,686}, {431,673}.              (42)
```

They form a literal `K_4` quotient of the Boolean cube: a vertex has harmonic
residue `2/1105`, an edge's endpoint union has residue `4/1105`, and every
perfect matching has the same all-vertex support and residue `8/1105`.
Therefore the harmonic residue loses all three XOR edge colours. Ordering the
finite constants would only impose a phase-dependent transitive tournament
gauge; it is not an intrinsic `V_4`-invariant orientation.

The tame hypothesis is sharp. With `u=(1,1)`, `v=(1,0)`, and `c=M=5`,

```text
Q(n)=||(n+5,n)||^2 congruent 2n^2 mod 5.                 (43)
```

The two-root cube collapses to one double root. One of the five parabolic
orbits sees it at every time and four see none. Also, weighting orbit `r` by
geometric hypotenuse gives `1/Q(r+ct)=O(t^(-2))` and a convergent branch, unlike the
positive-density harmonic time support. Both hostiles show why the index
carrier and `gcd(c,M)=1` must be retained.

## 8. Equality and failure boundaries

1. Every `A subset N_(>0)` defines a harmonic subseries
   `sum_(n in A)1/n`; finite-state structure is not automatic.  In fact
   `a_n=n^2` is itself C-finite, with `(E-1)^3a=0`, but its **value support**
   `{n^2:n>=1}` has unbounded gaps and reciprocal mass `zeta(2)`.  The theorem
   applies instead to modular **index supports**: for example
   `{n>=1:n^2=0 mod 5}=5N_(>0)` has `delta=1/5`.  Likewise
   `{ceil(k log(k+1)):k>=1}` has density zero and divergent mass.
2. The theorem concerns the **index support** (2), not the set of distinct
   recurrence values. Collapsing repeated values can change the answer.
3. Exact zero drift in `{-1,0,1}^*` is not regular: intersecting with
   `+^* -^*` gives `{+^k -^k}`. Drift modulo `q` is a lawful finite-state
   quotient, but it forgets absolute drift.
4. Ties are genuine states. With all three ternary letters, the number of
   length-`n` words returning modulo three is `3^(n-1)`, which is odd. With
   only `+/-1` it is `(2^n+2(-1)^n)/3`, even for every `n>=1`. Their modular
   harmonic supports are respectively all positive lengths and empty.
5. A transitive `T_4` comparison is too coarse for Fibonacci ancestry. The
   windows `(1,2,3,5)` and `(2,3,5,8)` induce the same order but have opposite
   Cassini sign and different THM-3339 ray states; `(1,1,2,3)` has a tie and
   gives no tournament without an added gauge.
6. Counts lose addresses. Two disjoint ternary languages can have the same
   `3^(n-1)` level count, and harmonic residue retains only their cycle
   density. It cannot reconstruct a branch word, owner, current, or source.
7. Reserved THM-3357/3358 are not dependencies. A future finite weighted
   Berggren transfer matrix would enter through (20), but harmonic periodicity
   would not prove its proposed parent circuit or composite compiler.

## 9. Exact companion and scope

The companion has no Python `assert` nodes and compares direct transfer-matrix
powers with recurrence-state iteration. It exhausts all labelled order-four
tournaments, their cut observers, modular singleton filters, complete
two-state ternary DFAs, the tie/no-tie hostile, Fibonacci parity, Berggren
depth counts, and the U-spine `K_4` clock. Ordinary and optimized transcripts
are byte-identical to the stored result. Exact counts and hashes are printed
there and pinned above.

This theorem proves an index-set and mean-scar classification. It does not
orient a generalized relation, classify arbitrary subsets of `N_(>0)`, recover
Berggren ancestry from counts, construct the reserved weighted/composite
compilers, prove an LRC physical mass, or settle LRC(14).
