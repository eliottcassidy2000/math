# Coupled twin lengths: a physical equality gate and a Farey interface

**PROVED below; independent review pending.** The tensor method is inherited
from THM-3121, `path-cover-walk-content-substitution-kernel`, and Session260's
`constructors/RESEARCH.md`. The specific gate, synchronized family and
arbitrary-module determinant interface are the present contribution. No
literature-priority or Chan–Pak-completeness claim is made.

## 1. The object and its complete two-port response

Let Q have outgoing bit rows `(6,4,8,19,39,15)` on0..5, selected set
O={0,4}, and unselected ports1,5. Replace these ports by nonempty tournaments
A,B with uniform inherited exterior orientation. The selected vertices stay
fixed. For any nonempty T, write C_r(T) for ordered r-nonempty-path spanning
covers, H(T)=C_1(T). Write (h,c,e)=(C_1,C_2,C_3)(A) and (g,d)=(C_1,C_2)(B).
The seed observables are exactly those in
`five_stress_connection_atlas_20261112/SEEDS.md`: W counts ordered nonempty
two-path covers whose first end and second start lie in O, U also allows an
empty first path and requires the second start selected, and δ=H−W.

**Theorem1 (universal module interface).** For every such A,B,

    H = 15hg+6cg+5hd+2cd,
    W = 14hg+5cg+6hd+2cd,
    δ = hg+cg−hd,
    U = 86hg+72cg+18eg+42hd+34cd+8ed.                 (1)

These are physical tournament substitutions. They do not permit arbitrary
positive tuples (h,c,e,g,d).

**Proof/certificate.** Compress each maximal visit to a module to one run.
There are four fixed exterior vertices. The two ports form an acyclic
quotient, so between consecutive exterior vertices a port cannot be
revisited. Every one- or two-path skeleton therefore has at most six visits
to each port. A finite recursion over the used exterior set, last type, and
whether the separator has been used exhausts every skeleton; repeated ports
increase the run indices, while fixed vertices are used once. Adjacent equal
port runs are forbidden because they are one maximal run. There are no
port-only directed cycles, so this recursion terminates without any assumed
length cutoff. The full nonzero coefficient tables are:

| observable | (1,1) | (2,1) | (3,1) | (1,2) | (2,2) | (3,2) |
|---|---:|---:|---:|---:|---:|---:|
| H |15|6|0|5|2|0|
| W |14|5|0|6|2|0|
| U |86|72|18|42|34|8|

Each skeleton of run-content(r,s) lifts in C_r(A)C_s(B) ways, because runs
are ordered by their occurrence. `kernel.py` enumerates this entire finite
object. `verify.py` checks the complete table and independent literal
vertex-level counts, including nontransitive modules. This proves(1) for
arbitrary module size, not only the literal verification sizes. ∎

Q is strong (e.g. it has the spanning cycle0→1→2→3→4→5→0).
Nonempty substitutions preserve strong connectivity, so every displayed
unit seed below is strong and has U>0 by(1).

## 2. Neither length is free, but their equality is renewable

For A=TT_a and B=TT_b, a,b>=1, C_1=1, C_2=2^a−2,
C_3=3^a−3·2^a+3. Thus

    δ(a,b)=1+2^a−2^b,
    H(a,b)=1+2·2^a+2^b+2·2^(a+b),
    U(a,b)=−2·2^a+2·3^a−2·2^b+10·2^(a+b)+8·3^a2^b.  (2)

**Corollary2.** The complete positive-length unit locus is exactly a=b.
At every point(a,a), changing either length alone destroys δ=1; increasing
both lengths by the same positive amount preserves it. Thus two individually
unstable ports can cooperate. This does not contradict the different open
hypothesis that *individually stable* ports are jointly stable.

At a=b=m,

    H_m=1+3·2^m+2·4^m,
    U_m=8·6^m+10·4^m+2·3^m−4·2^m,
    c_m=U_m/2=4·6^m+5·4^m+3^m−2·2^m.                (3)

The registers x=2^m,y=3^m update by(x,y)↦(2x,3y), and the formulas are
polynomial in these registers. U satisfies the exact recurrence
U_(m+4)=15U_(m+3)−80U_(m+2)+180U_(m+1)−144U_m; H satisfies
H_(m+3)=7H_(m+2)−14H_(m+1)+8H_m. These are synchronized physical growth
rules, not arbitrary integer recurrences detached from tournaments.

c_m is odd for every m>=1. The inherited selected-tail construction
therefore gives literal counts1+c_m2^s for every s>=2. It does **not**
supply the missing s=1 head. Also c_m~4·6^m, so this particular coefficient
family has O(log X) members below X and cannot by itself cover all odd
coefficients. No unrestricted exclusion follows.

## 3. Replacing length equality by physical Farey adjacency

Equation(1) says precisely

    δ=1  iff  g(h+c)−hd=1
         iff (h+c)/h − d/g = 1/(hg).                (4)

Thus the two fractions are Farey neighbors, with the first larger.
The equation forces gcd(h,c)=gcd(g,d)=gcd(h,g)=1. Since h,g are odd and
C_2 is even (pair ordered path covers by exchanging the two paths), the
numerator h+c is odd, while d is even. These necessary arithmetic facts do
not imply physical realizability of the module profiles.

There are actual nontransitive examples:

* A rows `(18,20,9,19,36,79,31)` has(h,c,e)=(65,846,4446).
  B=TT_4 has(g,d)=(1,14). Equation(4) is911−910=1. The resulting
  fifteen-vertex seed has(H,U,δ)=(34289,1085398,1).
* A rows `(30,4,8,18,6)` has(h,c)=(5,54); B rows
  `(114,164,73,19,38,76,154,61)` has(g,d)=(589,6950).
  Here589·59−5·6950=1. Both modules are nontransitive. Their seventeen-vertex
  substituted graph and full U are independently counted in `REPLAY.json`.

A complete inherited profile census through module order8 yields eight
ordered nontransitive profile pairs satisfying(4), among1922 distinct
(H,C_2) profiles. This is a **FINITE-EXACT profile statement** conditional
on the inherited class inventory, not eight isomorphism pairs or an all-order
classification. `module_pairs.json` stores representatives.

There is also a useful all-order search bound. Cutting each Hamilton path
at each of n−1 edges injects into ordered two-path covers; hence
C_2(T)>= (n−1)H(T). Given A, (4) gives

    |B|−1 <= d/g = 1+c/h−1/(hg) < 1+c/h,
    hence |B| < 2+c/h.                              (5)

Thus all possible partner modules B for a *fixed physical A* lie in an
explicit finite range. This bounds partner size, not the size of every
possible A or the work of enumeration. In particular:

* A=TT_k forces C_2(B)=(2^k−1)H(B)−1 and |B|<2^k.
* B=TT_k forces C_2(A)=(2^k−3)H(A)+1. If H(A)>1, then
  |A|<=2^k−2; if H(A)=1, A is TT_k by the equation.

These are finite partner-search bounds for each k, not independent
inflation permissions.

The directed physical partner graph also has a strict potential off the
transitive loops. Writing r(T)=C_2(T)/H(T), (4) says

    r(B)−r(A)=1−1/(H(A)H(B)).

Every edge involving a nontransitive module strictly increases r by less
than1. Equal-ratio edges have H(A)=H(B)=1 and are precisely equal-sized
transitive loops. Thus after deleting those loops the partner graph is
acyclic. The size bounds make both in- and out-neighborhoods of each fixed
physical module finite. This is a locally finite directed grammar for
**this substitution interface**, not a grammar for every unit seed.


## 4. The exact renewal boundary of the Farey lift

Put q_A=c+2h and q_B=d+2g, the ordered two-slot cover counts permitting
empty slots. The inherited ordinal-product law is
Q_2(A→B)=Q_2(A)Q_2(B), while H(A→B)=H(A)H(B).
Simultaneously append(or prepend) TT_t to both modules. Then h,g stay fixed
and q_A,q_B multiply by2^t. At an existing unit pair,

    δ_t=hg+2^t(1−hg)=1+(hg−1)(1−2^t).               (6)

For t>0 it remains a unit seed iff hg=1. A tournament has H=1 iff it is
transitive: a nontransitive tournament has a nontrivial strong component,
whose spanning cycle gives at least three Hamilton paths, and strong
components multiply. So simultaneous neutral ordinal padding renews the
chain-diagonal family, but cannot renew either nontransitive example by
this operation. In fact the same obstruction is complete for **independent** ordinal
padding. Let a,b>=0 be the total added transitive lengths on the two sides,
x=gq_A and K=hg−1. At a unit pair hq_B=x+K, hence

    δ_(a,b)=hg+x2^a−(x+K)2^b.                         (7)

If K>0, then x>=2hg>K. For a<=b the last difference is at most−K2^b,
and equality with−K requires a=b=0. For a>b it is
2^b[x(2^(a−b)−1)−K]>0, so cannot equal−K. Therefore a nontransitive
unit pair permits **only(a,b)=(0,0)**. If K=0, the exact allowed locus is
a=b. This includes arbitrary transitive collars distributed before and
after each module because only their total lengths enter Q_2.

Abstract determinant-one transformations must not be confused with
physical realizable profile transformations.

## 5. A complete direction classifier for any fixed acyclic-port tensor

Let a fixed acyclic unselected-port inflation chart have positive-length
response δ(m)=sum_α c_α prod_i α_i^(m_i), with finitely many positive integer
base vectors α. This finite exponential representation follows from the
inherited Stirling skeleton theorem. Subtract the target1 and combine equal
base vectors. Fix a positive integer basepoint b. For an integer direction
v>=0, v≠0, put B_α(v)=prod_i α_i^(v_i) and w_α=c_α prod_i α_i^(b_i).
Then

    δ(b+tv)−1 = sum_B [sum_{α:B_α(v)=B} w_α] B^t.    (8)

**Theorem3.** The ray is unit for every t>=0 iff every bracket in(8)
vanishes. All such directions form a finite union of rational linear
subspaces intersected with the nonnegative integer lattice. An exact finite
classifier is:

1. Enumerate partitions of the finite base set into blocks of zero total
   weight w_α.
2. For every block require all its B_α(v) equal. Factoring the bases makes
   these equations rational linear: for every prime p,
   sum_i(v_p(α_i)−v_p(β_i))v_i=0.
3. Take the union of those solution subspaces, restricting to v>=0,v≠0.

Proof: distinct positive exponentials are linearly independent, by successive
largest-base limits. Any actual cancelling ray supplies its equal-base
partition. Conversely a zero-weight partition with within-block collisions
suffices, even if additional blocks collide. Thus exclusions or strict
chamber inequalities are unnecessary. ∎

If there are s distinct nonzero grouped bases along a non-unit ray, it hits
δ=1 at at most s−1 real values of t, and therefore at most s−1 admissible
integer points. Proof: a nonzero sum of s exponentials with distinct real
exponents has at most s−1 real zeros, inductively divide by the smallest
exponential, differentiate, and apply Rolle. The domain may be restricted
to t>=0 without weakening this bound. This is a ray bound; it is not a
classification of all nonlinear multi-length solution sets.

For(2) at basepoint(1,1), the subtracted weights are+2 at base(2,1) and−2
at base(1,2). The only nonzero direction subspace is v_1=v_2, exactly
recovering the physical synchronization gate. This is a genuine toric
collision criterion, not a numerical guess from diagonal samples.

## 6. Provenance, experiments and limits

The earlier atlas supplied stable ports, isolated C5 resonances and the
stable-clone contraction lemma; it did not justify treating every pair of
unstable ports as independently bounded. The new gate is a counterexample
to that tempting inference. The finite screen tested12300 unselected-port
pairs in strong marked seed representatives through order8 and found14
nonconstant tensors constant on the diagonal; no minimality claim uses the
screen. The smallest displayed core has order6.

Classical Farey determinants and exponential independence are not new.
The bridge carries actual module profiles, endpoint marks and their full
run tensor. It loses arbitrary integer realizability; the sidecars are the
physical module witnesses and continuing U response. Chan–Pak, a general
stable-port theorem, a complete nonlinear length-locus classification and
unrestricted prime-core classification remain OPEN.
