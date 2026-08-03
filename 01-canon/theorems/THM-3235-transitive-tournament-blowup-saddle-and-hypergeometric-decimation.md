---
id: THM-3235
title: "Transitive tournament blowup saddle, hypergeometric decimation, and observer twins"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Tournament substitution obeys a general Schur-complement determinant law.
  For a uniform transitive m-blowup of a strong b-vertex quotient, the unique
  saddle splits each base coordinate into m equal coordinates, the product
  radius is R_P^m/m^(bm), and the factorial-normalized fixed-depth sequence is
  an exact hypergeometric multiple of the m-section of the base sequence.
  Determinant symmetries partition saddle coordinates into owner classes.
  Two explicit nonisomorphic order-six quotients have the same complete
  path-colour profile, product radius, and non-P-recursive class but different
  uniform-lift Hamilton responses, proving that those scalars are not a
  sufficient substitution state.
source: root/creative-reframes/2026-08-03
audit: >
  The assertion-independent exact companion reconstructs both order-six masks,
  both multivariate determinants from principal minors, the full determinant
  symmetry groups of orders 48 and 6, their transitive owner partitions, both
  determinant factorizations, complete backward/path-cover/path-colour
  profiles, 16 direct uniform determinant substitutions, one arbitrary-factor
  block determinant control, four substitution-associativity controls, a
  generic walk-observer split, and Held--Karp lift responses through order 18.
  Normal, optimized, and stored transcripts agree byte-for-byte; the script
  AST contains no assert node.  An independent hostile audit rederived the
  block/Schur law, Neumann equivalence, saddle and decimation formulas,
  recurrence closure, owner partition, both exact twins, and all scope limits.
depends_on:
  - THM-3226-unbalanced-q4-unequal-saddle-and-transcendence-wall
  - THM-3213-tournament-normalized-cyclic-diagonal-and-fast-moving-jet-transform
  - THM-3166-falling-factorial-order-join-path-colour-transform
  - THM-3121-path-cover-walk-content-substitution-kernel
related:
  - THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence
  - THM-3202-c3-repeated-join-moving-jet-formula-and-cfinite-obstruction
script: 04-computation/tournament_transitive_blowup_radius_observer_twins_thm3235.py
output: 05-knowledge/results/tournament_transitive_blowup_radius_observer_twins_thm3235.out
script_sha256: eac84038d8820be4bcc6004ee7da2974bf955f3f21f3e4bd3db3e07d25ae09fa
output_sha256: 689bd2788613967149f00564f49df3272aaf16d353f9edb7d1712c352ee3cdaa
hash_basis: LF-normalized bytes
---

# THM-3235 -- transitive tournament blowup saddle, hypergeometric decimation, and observer twins

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem connects three coordinates previously used separately: the
quotient determinant controlling the positive saddle, the
factorial-normalized repeated-substitution sequence, and the operation state
consumed by a later substitution. Uniform transitive blowup has an exact law
in all three coordinates. The final order-six pair marks the information
boundary: complete path-colour data and the product radius still do not
determine the quotient as an operator.

## 1. Conventions and the general block determinant law

For a tournament \(R\) on \([q]=\{1,\ldots,q\}\), use the row-to-column
adjacency convention

~~~text
(A_R)_(ij)=1  iff  i->j.
~~~

For \(X=\operatorname{diag}(X_1,\ldots,X_q)\), put

\[
 D_R(X_1,\ldots,X_q)=\det(I-A_RX),\qquad
 W_R(X)={\bf1}^{T}X(I-A_RX)^{-1}{\bf1}.                 \tag{1}
\]

The second expression is the nonempty directed-walk observer of THM-3121. It
is a rational function; the first is a multiaffine polynomial.

Let \(P\) be a tournament on \([b]\). For each \(i\), let \(S_i\) be a
tournament on a disjoint vertex set, let \(X_i\) be its diagonal variable
matrix, and form the substitution \(P[S_1,\ldots,S_b]\): every arc \(i\to j\)
of \(P\) becomes all arcs from \(S_i\) to \(S_j\). Define

\[
 D_i=I-A_{S_i}X_i,\qquad
 w_i={\bf1}^{T}X_iD_i^{-1}{\bf1}.                       \tag{2}
\]

Then the exact **block determinant law** is

\[
 \boxed{
 \det(I-A_{P[S_1,\ldots,S_b]}X)
   =\left(\prod_{i=1}^{b}\det D_i\right)
      D_P(w_1,\ldots,w_b).}                              \tag{3}
\]

Equation (3) is initially an identity in the localization where the \(D_i\)
are invertible. After the displayed denominators are cleared, it is a
polynomial identity and therefore holds identically.

### Proof

Let \(D=\operatorname{diag}(D_1,\ldots,D_b)\). Let \(L\) be the block-column
matrix whose \(i\)-th column is the all-one vector on block \(i\), and let
\(R\) be the block-row matrix whose \(i\)-th row is \({\bf1}^{T}X_i\) on
block \(i\). With the stated adjacency convention,

\[
 I-A_{P[S]}X=D-LA_PR.                                   \tag{4}
\]

Sylvester's determinant identity gives

\[
\begin{aligned}
 \det(D-LA_PR)
 &=\det(D)\det(I-D^{-1}LA_PR)\\
 &=\det(D)\det(I-A_PRD^{-1}L).
\end{aligned}                                           \tag{5}
\]

The matrix \(RD^{-1}L\) is diagonal with entries \(w_i\). Substitution into
(5) is exactly (3). No strong-connectivity or uniformity hypothesis was used.
\(\square\)

## 2. Uniform transitive blocks

Let \(T_m\) be the transitive tournament \(1\to2\to\cdots\to m\) and write

\[
 P^{\langle m\rangle}=P[T_m,\ldots,T_m].                \tag{6}
\]

Because \(A_{T_m}X\) is strictly upper triangular,

\[
 \det(I-A_{T_m}X)=1.                                    \tag{7}
\]

Every nonempty subset of \([m]\) supports exactly one increasing directed
path. Expanding the inverse in (1) therefore gives

\[
 W_{T_m}(X_1,\ldots,X_m)=\prod_{j=1}^{m}(1+X_j)-1.      \tag{8}
\]

If the variables in base block \(i\) are
\(X_{i1},\ldots,X_{im}\), (3) becomes

\[
 \boxed{
 D_{P^{\langle m\rangle}}(X)
  =D_P(Y_1,\ldots,Y_b),\qquad
 Y_i=\prod_{j=1}^{m}(1+X_{ij})-1.}                       \tag{9}
\]

Put \(u(x)=\exp(x)-1\). Since \(1+u(x)=\exp(x)\), the exponential coordinate
linearizes every block:

\[
 \boxed{
 D_{P^{\langle m\rangle}}(u(x_{11}),\ldots,u(x_{bm}))
 =D_P\left(u\left(\sum_jx_{1j}\right),\ldots,
           u\left(\sum_jx_{bj}\right)\right).}           \tag{10}
\]

Thus transitive substitution is multiplicative in the raw variables and
additive in the exponential saddle variables.

## 3. Exact saddle and radius renormalization

Assume now that \(P\) is strong on \(b\ge3\) vertices. Let

\[
 \mathcal N_P=
 \{x\in\mathbb R_{\ge0}^{b}:
   \rho(A_P\operatorname{diag}(u(x_i)))<1\}.             \tag{11}
\]

By THM-3226, there is a unique support point

\[
 c=(c_1,\ldots,c_b)
 =\arg\max_{x\in\overline{\mathcal N_P}}\prod_i x_i,
 \qquad R_P=\prod_i c_i,                                 \tag{12}
\]

and it is the unique strictly minimal, nondegenerate positive diagonal
saddle.

For \(P^{\langle m\rangle}\), write \(S_i=\sum_jx_{ij}\). The transitive
diagonal blocks in the Schur complement have nonnegative inverses for all
positive variables. Consequently the nonsingular-M-matrix criterion,
together with (10), gives

\[
 x\in\mathcal N_{P^{\langle m\rangle}}
 \quad\Longleftrightarrow\quad
 (S_1,\ldots,S_b)\in\mathcal N_P.                        \tag{13}
\]

At fixed block sums, AM--GM gives

\[
 \prod_{j=1}^{m}x_{ij}\le (S_i/m)^m,                    \tag{14}
\]

with equality only when all \(x_{ij}=S_i/m\). Maximizing over the sums and
using the uniqueness in (12) proves

\[
 \boxed{c_{ij}^{\langle m\rangle}=c_i/m}                \tag{15}
\]

and

\[
 \boxed{
 R_{P^{\langle m\rangle}}
 =\prod_{i,j}c_{ij}^{\langle m\rangle}
 =\frac{R_P^m}{m^{bm}}.}                                 \tag{16}
\]

This is an equality of exact analytic values, not a numerical saddle fit.
The quotient \(P^{\langle m\rangle}\) is strong, so THM-3226 also supplies
strict complex minimality and a nondegenerate \((bm-1)\)-dimensional saddle.
For every fixed depth \(d\ge1\),

\[
 \frac{F_{P^{\langle m\rangle}[T_r,\ldots,T_r]}(d)}
        {(r!)^{bm}}
 \sim C_{P,m,d}\,
 r^{d-1-(bm-1)/2}
 \left(\frac{R_P^m}{m^{bm}}\right)^{-r},
 \qquad C_{P,m,d}>0.                                    \tag{17}
\]

Recall the canonical definition: if \(A_Pv=\rho v\) and
\(w^TA_P=\rho w^T\), then \(P\) is **Perron-balanced** when
\(w_iv_i=(w^Tv)/b\) for every \(i\). THM-3213 and the Lagrange equations of
THM-3226 make this equivalent to the unique support point being diagonal.
Equation (15) therefore gives

\[
 P^{\langle m\rangle}\text{ is Perron-balanced}
 \quad\Longleftrightarrow\quad
 P\text{ is Perron-balanced}.                            \tag{18}
\]

Moreover, because the radii are positive,

\[
 R_{P^{\langle m\rangle}}\text{ is algebraic}
 \quad\Longleftrightarrow\quad R_P\text{ is algebraic}. \tag{19}
\]

Indeed one direction follows from (16), and the other follows because \(R_P\)
is a root of \(z^m-m^{bm}R_{P^{\langle m\rangle}}\).

For the unbalanced quotient \(Q4\) of THM-3226, whose saddle is
\((s,t,t,s)\), the blowup has the two repeated scales \(s/m,t/m\) and

\[
 R_{Q4^{\langle m\rangle}}=R_{Q4}^{m}/m^{4m}.            \tag{20}
\]

Thus uniform blowups transport the Q4 arithmetic wall exactly; they do not
resolve the open transcendence of \(R_{Q4}\).

## 4. Exact coefficient decimation

For a tournament \(P\) on \(b\) vertices, let \(F_R(d)=d!p_R(d)\) denote the
number of ordered spanning covers of \(R\) by \(d\) directed paths, and define

\[
 a_{P,d}(r)=\frac{F_{P[T_r,\ldots,T_r]}(d)}{(r!)^b}
 \qquad(r\ge1).                                         \tag{21}
\]

Tournament substitution is associative and

\[
 T_m[T_r,\ldots,T_r]\cong T_{mr}.                       \tag{22}
\]

Therefore

\[
 P^{\langle m\rangle}[T_r,\ldots,T_r]
 \cong P[T_{mr},\ldots,T_{mr}],                          \tag{23}
\]

so the raw ordered-cover counts on the two sides are identical. Dividing by
the two normalizing factorials proves the exact all-depth decimation law

\[
 \boxed{
 a_{P^{\langle m\rangle},d}(r)
 =\left(\frac{(mr)!}{(r!)^m}\right)^b a_{P,d}(mr).}      \tag{24}
\]

The multiplier

\[
 h_{m,b}(r)=\left(\frac{(mr)!}{(r!)^m}\right)^b
\]

is hypergeometric, with

\[
 \frac{h_{m,b}(r+1)}{h_{m,b}(r)}
 =\left(\frac{\prod_{k=1}^{m}(mr+k)}{(r+1)^m}\right)^b. \tag{25}
\]

Hence, under the usual characteristic-zero definition, P-recursiveness of
the full base sequence \(\{a_{P,d}(r)\}\) implies P-recursiveness of the
blowup sequence: arithmetic subsequences of P-recursive sequences are
P-recursive, and termwise multiplication by a hypergeometric sequence
preserves P-recursiveness. For completeness, in the equivalent D-finite
ordinary-generating-function coordinate, the \(m\)-section is the
root-of-unity filter

\[
 \sum_{r\ge1}a_{P,d}(mr)z^r
 =\frac1m\sum_{j=0}^{m-1}
 A\!\left(\zeta_m^jz^{1/m}\right),\qquad
 A(z)=\sum_{n\ge1}a_{P,d}(n)z^n,                         \tag{25a}
\]

and D-finite series are closed under algebraic substitution and finite sums.
Termwise multiplication is the Hadamard product, under which D-finite series
are also closed; (25) supplies a first-order recurrence for the
hypergeometric multiplier. Conversely, (24) recovers only the \(m\)-section
\(\{a_{P,d}(mr)\}\). No converse for the full base sequence is claimed.

Finally, Stirling's formula gives

\[
 \frac{(mr)!}{(r!)^m}
 \sim m^{mr}\sqrt m\,(2\pi r)^{-(m-1)/2}.               \tag{26}
\]

Substituting the THM-3226 base asymptotic into (24) changes the power
\(d-1-(b-1)/2\) exactly to \(d-1-(bm-1)/2\) and the radius exactly to (16).
Thus the combinatorial decimation and analytic saddle laws agree.

## 5. Determinant owners and the saddle partition

For any strong tournament \(Q\) on \(q\) vertices, define its
determinant-owner group

\[
 \Gamma_Q=\{\sigma\in S_q:
 D_Q(X_{\sigma(1)},\ldots,X_{\sigma(q)})
 =D_Q(X_1,\ldots,X_q)\}.                                \tag{27}
\]

The positive Neumann domain in the raw \(X\) variables is exactly the
component of
\(\{X\in\mathbb R_{\ge0}^{q}:D_Q(X)\ne0\}\) containing the origin. The region
\(\rho(A_QX)<1\) is joined to the origin by radial segments, while any path
from it to a point with Perron radius greater than one must cross \(\rho=1\),
where \(D_Q=0\). Every element of \(\Gamma_Q\) fixes the origin and preserves
\(D_Q\), hence preserves that component. It also preserves the product
objective.

After the coordinatewise map \(X_i=u(x_i)\), uniqueness of the THM-3226
support point therefore forces

\[
 \sigma(c_Q)=c_Q\qquad(\sigma\in\Gamma_Q).               \tag{28}
\]

Thus saddle coordinates are constant on the orbit partition of \(\Gamma_Q\).
In particular,

\[
 \Gamma_Q\text{ transitive on vertices}
 \quad\Longrightarrow\quad Q\text{ Perron-balanced}.    \tag{29}
\]

This is a determinant symmetry, not necessarily a tournament automorphism.
It is a lawful owner compression for the saddle and radius, but it need not
preserve the full walk observer \(W_Q\).

## 6. Two order-six observer twins

Use vertices \(1,\ldots,6\). For a mask, list pairs \((i,j)\), \(i<j\), in
lexicographic order; bit zero belongs to \((1,2)\), and bit one means
\(i\to j\) (bit zero means \(j\to i\)). Masks \(408\) and \(1332\) give

\[
A_{408}=\begin{pmatrix}
0&0&0&0&1&1\\
1&0&0&0&1&1\\
1&1&0&0&0&0\\
1&1&1&0&0&0\\
0&0&1&1&0&0\\
0&0&1&1&1&0
\end{pmatrix},\qquad
A_{1332}=\begin{pmatrix}
0&0&0&1&0&1\\
1&0&1&0&0&1\\
1&0&0&0&1&0\\
0&1&1&0&0&0\\
1&1&0&1&0&0\\
0&0&1&1&1&0
\end{pmatrix}.                                          \tag{30}
\]

The directed Hamilton cycles

~~~text
408:  (1,5,3,2,6,4),
1332: (1,4,2,6,3,5)
~~~

prove that both tournaments are strong.

### 6.1 Multivariate determinants and owner groups

Direct principal-minor expansion gives

\[
 D_{408}=1-(X_1+X_2+X_1X_2)
             (X_3+X_4+X_3X_4)
             (X_5+X_6+X_5X_6).                          \tag{31}
\]

Thus \(Q_{408}\) is isomorphic to \(C_3[T_2,T_2,T_2]\). Its full
determinant-owner group is \(S_2\wr S_3\), of order \(48\), and is
transitive.

For a monomial \(M\), let \(\operatorname{Orb}_{\sigma}(M)\) denote the sum
of its **distinct** images under the cyclic permutation

\[
 \sigma=(1\ 2\ 3\ 6\ 4\ 5).                            \tag{32}
\]

Then

\[
\begin{aligned}
D_{1332}={}&1
-\operatorname{Orb}_\sigma(X_1X_2X_4)
-\operatorname{Orb}_\sigma(X_1X_3X_4)\\
&-\operatorname{Orb}_\sigma(X_1X_2X_3X_4)
-\operatorname{Orb}_\sigma(X_2X_3X_4X_5)\\
&-2\operatorname{Orb}_\sigma(X_1X_2X_3X_4X_5)
-4X_1X_2X_3X_4X_5X_6.                                  \tag{33}
\end{aligned}
\]

The six orbit sizes in (33), including the final fixed monomial, are

\[
 (6,2,6,3,6,1).                                         \tag{34}
\]

The exact companion enumerates all \(6!\) coordinate permutations and proves
that the full owner group of (33) is precisely the cyclic group generated by
\(\sigma\), of order \(6\); it too is transitive. This exhaustive statement
uses the exact coefficient table

\[
 [X_S]D_Q=(-1)^{|S|}\det A_Q[S].                         \tag{35}
\]

The two tournaments are nonisomorphic: variable relabelling preserves every
multivariate determinant coefficient, whereas (31) has only coefficient
\(-1\) beyond its constant term while (33) has coefficients \(-2\) and
\(-4\).

### 6.2 The same saddle radius and recurrence class

Equal specialization gives

\[
\begin{aligned}
D_{408}(t,\ldots,t)
 &=1-8t^3-12t^4-6t^5-t^6\\
 &=-(t^2+2t-1)(t^4+4t^3+5t^2+2t+1),                    \tag{36}\\
D_{1332}(t,\ldots,t)
 &=1-8t^3-9t^4-12t^5-4t^6\\
 &=-(t^2+2t-1)(2t^2+t+1)^2.                             \tag{37}
\end{aligned}
\]

The remaining factors are positive for \(t>0\), so both first positive poles
are

\[
 \alpha=\sqrt2-1,\qquad \rho(A_Q)=\alpha^{-1}=1+\sqrt2. \tag{38}
\]

Both owner groups are transitive. Equations (28)--(29) therefore make the
unique saddle diagonal, and \(u(c)=\alpha\) gives

\[
 c=\log(1+\alpha)=\log\sqrt2=\frac{\log2}{2},
 \qquad R_{408}=R_{1332}=\left(\frac{\log2}{2}\right)^6. \tag{39}
\]

Hermite--Lindemann makes \(\log 2\), and hence the nonzero power in (39),
transcendental. The D-finite singularity obstruction used in THM-3213 now
shows that, for every fixed \(d\ge1\), both normalized sequences
\(a_{Q,d}(r)\) are non-P-recursive. The raw ordered-cover sequences are also
non-P-recursive: if a raw sequence were P-recursive, termwise multiplication
by the hypergeometric sequence \((r!)^{-6}\) would make its normalization
P-recursive. Thus the twins have the same product radius and the same
non-P-recursive recurrence class at every fixed depth.

### 6.3 The same complete path-colour profile

For a vertex permutation \(\pi\), let \(b_Q(\pi)\) be the number of
consecutive pairs oriented backward along \(\pi\). Exhaustion of the \(6!\)
permutations gives, for both tournaments,

\[
 \bigl(\#\{\pi:b_Q(\pi)=j\}\bigr)_{j=0}^{5}
 =(45,117,198,198,117,45).                               \tag{40}
\]

The cut dictionary of THM-3166 then gives the complete ordered and unordered
spanning path-cover profiles

\[
 (F_Q(d))_{d=1}^{6}=(45,342,1116,1944,1800,720),         \tag{41}
\]

\[
 (p_Q(d))_{d=1}^{6}=(45,171,186,81,15,1).               \tag{42}
\]

Consequently their complete injective path-colour polynomials coincide:

\[
 \boxed{Q_{408}(t)=Q_{1332}(t)=t^6+16t^4+28t^2.}        \tag{43}
\]

By THM-3166, every repeated order-join profile obtained from the two seeds is
therefore identical as well.

### 6.4 Uniform lifts split the twins

Define the uniform transitive-lift Hamilton response

\[
 h_Q(r)=H(Q[T_r,\ldots,T_r]).                            \tag{44}
\]

Exact Held--Karp subset DP gives

\[
\begin{array}{c|rr}
r&h_{408}(r)&h_{1332}(r)\\ \hline
1&45&45\\
2&421425&411513\\
3&73948095585&71347853961
\end{array}                                             \tag{45}
\]

Already the first two response columns have nonzero determinant

\[
 \det\begin{pmatrix}45&421425\\45&411513\end{pmatrix}
 =-446040.                                               \tag{46}
\]

Hence the response rank of this two-object family for the consumer (44) is
two. Complete path-colour profile, product radius, and recurrence class are a
**false compression** for quotient substitution. The rational observer

\[
 W_Q(X)={\bf1}^{T}X(I-A_QX)^{-1}{\bf1}                  \tag{47}
\]

is a sufficient sidecar for the THM-3121 substitution law and distinguishes
the pair at generic variables. No claim that (47) is an absolutely minimal
state for every possible consumer is made.

## 7. Scope and nonconsequences

- The determinant identity (3) holds for arbitrary tournament factors and
  arbitrary block sizes. The closed saddle formula (15)--(16) requires a
  strong base and **uniform transitive** factors.

- With unequal transitive block sizes \(m_i\), fixed-sum AM--GM produces the
  weighted objective \(\prod_i S_i^{m_i}/m_i^{m_i}\). Its maximizer need not
  be the unweighted base saddle, so \(c_i/m_i\) is not a general formula.

- For arbitrary internal factors, the prefactors \(\det D_i\) in (3) have
  their own poles. They cannot be discarded when locating the dominant
  boundary.

- Equations (15)--(19) classify blowups of a supplied base quotient; they do
  not classify all unbalanced strong quotients. In particular, (20) inherits
  rather than solves the Q4 transcendence wall.

- Determinant-owner transitivity is sufficient, not necessary, for Perron
  balance. The owner group may strictly exceed the tournament automorphism
  group and does not reconstruct the tournament.

- Equality of the complete path-colour profile, equal-specialized determinant
  pole, product radius, or recurrence class does not imply equality of the
  multivariate walk observer or substitution response. The pair (30) is an
  exact hostile witness. No equality of leading asymptotic constants is
  asserted.

- The response rank in (46) is consumer-relative. It is not an intrinsic
  claim that every lawful representation of the two tournaments has rank two.

- Equation (24) is exact decimation, not an arbitrary sequence-complexity
  collapse. It proves the forward P-recursive closure statement above but no
  converse from one arithmetic subsequence to the full base sequence.

## 8. Exact companion and reproduction

Run from the repository root:

~~~text
python3 04-computation/tournament_transitive_blowup_radius_observer_twins_thm3235.py
python3 -O 04-computation/tournament_transitive_blowup_radius_observer_twins_thm3235.py
~~~

The companion uses only the Python standard library. Its determinant engine
is fraction-free Bareiss elimination; its arbitrary-factor control uses exact
Fraction arithmetic; its order-18 Hamilton counts use unsigned-64-bit
Held--Karp tables (all counts are below \(18!\)). Every verification gate
raises RuntimeError explicitly. The script parses its own AST and requires
zero Assert nodes, so optimized mode cannot erase a mathematical check.

For a literal three-way transcript comparison and LF-normalized hashes:

~~~bash
thm3235_normal=$(mktemp)
thm3235_optimized=$(mktemp)
python3 04-computation/tournament_transitive_blowup_radius_observer_twins_thm3235.py > "$thm3235_normal"
python3 -O 04-computation/tournament_transitive_blowup_radius_observer_twins_thm3235.py > "$thm3235_optimized"
cmp "$thm3235_normal" "$thm3235_optimized"
cmp "$thm3235_normal" 05-knowledge/results/tournament_transitive_blowup_radius_observer_twins_thm3235.out
python3 - <<'PY'
from hashlib import sha256
from pathlib import Path
for name in (
    "04-computation/tournament_transitive_blowup_radius_observer_twins_thm3235.py",
    "05-knowledge/results/tournament_transitive_blowup_radius_observer_twins_thm3235.out",
):
    data = Path(name).read_bytes().replace(b"\r\n", b"\n")
    print(sha256(data).hexdigest(), name)
PY
~~~

**QED.**
