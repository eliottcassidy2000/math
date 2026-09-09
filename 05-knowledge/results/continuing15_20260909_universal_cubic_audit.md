# Independent referee: the real trace-square-zero cubic inequality

**Status: INDEPENDENTLY ACCEPTED without a mathematical repair.** The
[primary report](continuing15_20260909_universal_cubic_bound.md) proves the
universal finite-dimensional real-matrix inequality, its complete equality
criterion, and the stated consequence for every nonnegative oriented
support. The proof closes the arbitrary-support extension left unproved
in THM-4466 /
`sharp-tournament-cubic-bound-for-common-edge-and-substitution-classes`.
This is a closure of the repository's stated gap, with no external
priority claim and no PDE trajectory claim.

The referee reviewed the primary report with SHA256

    2ecfe0ee8fd7fe2aa65113a7e6c29df0f2214cefa45e739444cfcb92fb2c1bd6.

That is its preacceptance text. A recorded prose-only status promotion is
authorized; no mathematical change is needed. The independent executable
checks the primary source, output, and certificate pins and stores this
historical reviewed-report pin. It imports and executes no primary code.

## 1. The exact all-dimensional statement passes

For every real n by n matrix M with `tr(M^2)=0`, including arbitrary
diagonal entries and opposing signed entries,

\[
       |\operatorname{tr}(M^3)|\le \|M\|_F^3/\sqrt3.              \tag{1}
\]

If M is nonzero, equality with positive cubic trace holds if and only if
M is normal and its eigenvalues with algebraic multiplicities are

\[
                r,\ r\omega,\ r\omega^2,\ 0,\ldots,0,
       \qquad r>0,\quad \omega=e^{2\pi i/3}.                    \tag{2}
\]

It then has rank three. Negative equality is exactly the negative of such
a matrix. The zero matrix is the separate equality case. Nonnegativity,
zero diagonal, and zero first trace are **not** hypotheses of (1).

Here is an independent check of every inequality and equality step. Write
the full spectrum as `lambda_i=a_i+i b_i`, with algebraic multiplicity.
Reality gives conjugate pairing. Trace-square zero gives

    S = sum a_i^2 = sum b_i^2.

In a complex Schur form, Frobenius energy E is the sum of squared absolute
diagonal entries and squared absolute strictly upper-triangular entries.
Thus `E>=2S`, with equality exactly for a normal matrix. This step is
valid for nonnormal and defective matrices; diagonalizability is never
assumed. If S=0 the spectrum is zero, so the cubic trace is zero. A
nonzero nilpotent then has strict inequality because E>0.

For S>0 the real cubic trace is

    sum a_i^3 - 3 sum a_i b_i^2.

If all **nonreal** eigenvalues have nonnegative real parts, the second
sum is nonnegative and the trace is at most `S^(3/2)`. This is strictly
below `(2S)^(3/2)/sqrt(3)`. Negative real eigenvalues cause no difficulty:
their imaginary part is zero, and replacing a signed cube by its absolute
cube only raises the upper bound.

Otherwise choose `-z<0` to be the smallest real part among the nonreal
eigenvalues. Its conjugate supplies a second occurrence. Remove those
two real parts from the cube sum. The remaining signed cube sum is at
most `(S-2z^2)^(3/2)`, by absolute values and the finite l3/l2 inequality.
Whenever `b_i` is nonzero we have `a_i>=-z`, whence

\[
 \operatorname{tr}(M^3)
 \le (S-2z^2)^{3/2}-2z^3+3zS.                                  \tag{3}
\]

This remains correct if some **real** eigenvalue is more negative than
-z: it has no imaginary-square contribution and stays inside the signed
cube estimate. This distinction is the potentially fragile step, and it
is correctly typed in the primary proof.

For `w=z/sqrt(S)`, the interval is `0<w<=1/sqrt(2)`. Put

    f(w)=(1-2w^2)^(3/2)-2w^3+3w.

Direct differentiation factors exactly as

    f'(w)=3 sqrt(1-2w^2) [sqrt(1-2w^2)-2w].

It is positive below `w=1/sqrt(6)` and negative above it in the interval's
interior. The unique maximum is `2sqrt(2/3)`. Hence (3) is at most
`(2S)^(3/2)/sqrt(3)`, which is at most `E^(3/2)/sqrt(3)`. Replacing M
by -M proves the other side of the absolute bound. No trace-zero
assumption or orientation enters this proof.

## 2. Equality, including the nonnormal boundary

For a nonzero matrix attaining positive equality, the first spectral case
is strict and cannot occur. In the remaining case every bound in the
chain must be an equality:

- The unique maximizing parameter gives `z^2=S/6`.
- Equality in the remaining signed cube norm leaves exactly one nonzero
  remaining real part, positive, of size `sqrt(S-2z^2)=2z`. All other
  remaining real parts vanish. Additional negative real eigenvalues or
  additional conjugate pairs with real part -z would make this norm step
  strict.
- Equality in `-sum a_i b_i^2<=zS` places every nonzero imaginary part at
  real part -z. Only the selected pair is then available. Its imaginary
  parts are `+sqrt(3)z` and `-sqrt(3)z` because their squared sum is S.
- Equality in the Schur energy estimate makes the matrix normal. In
  particular a nilpotent Jordan contribution at the zero eigenvalue
  cannot survive equality.

This is exactly (2) with r=2z. Conversely a normal matrix with that
spectrum has `E=3r^2` and `tr(M^3)=3r^3`, so equality holds. The proof
therefore supplies an iff statement, not only a necessary spectral test.

Normality alone is insufficient: two disjoint equal cyclic blocks have
two active spectral triples and strict inequality. Spectrum alone is
also insufficient: nonorthogonal similarity can keep the cubic trace
fixed while increasing Frobenius energy. The independent rational
controls below test both boundaries.

## 3. Exact transfer to the entire weighted oriented cone

For nonnegative entries, zero diagonal, and `M_ij M_ji=0`, the trace-square
hypothesis holds term by term. Every directed cyclic triangle contributes
its weight product from three starting vertices in `tr(M^3)`, so

    tr(M^3)=3C.

The previously proved THM-4465 /
`weighted-tournament-shear-production-and-contact-kernel` supplies

    F=4 tr(sym(M) skew(M)^2)=3C-T.

Since T is nonnegative in this cone, (1) yields for every finite dimension
and every nonnegative amplitude pattern

\[
                    F\le3C\le E^{3/2}/\sqrt3.                  \tag{4}
\]

Missing edges are allowed; no arbitrary orientation of a tie is needed.
No common cyclic edge or substitution decomposition is assumed. Positive
equality for 3C has exactly criterion (2), and equality for F additionally
requires T=0. The zero matrix is again separate.

The incoming six-vertex equality example has the edge `u->v` of weight
two, four edges `v->x_i` of weight one, and four edges `x_i->u` of weight
one. Direct independent contraction gives `E=12,C=8,T=0,F=24`. Its matrix
is normal of rank three. Thus the new equality theorem retains the old
multi-return equality and does not incorrectly force a single triangle.

Nonnegative amplitudes are essential to the production comparison used
in (4). The signed transitive matrix

    [[0,1,-1],[0,0,1],[0,0,0]]

has `tr(M^2)=tr(M^3)=0` and still satisfies (1), but `C=0,T=-1,F=1`.
It refutes `F<=tr(M^3)` outside the nonnegative cone. Omitting
`tr(M^2)=0` from the matrix theorem itself already fails for the real
one-by-one matrix `[1]`. These are different load-bearing hypotheses.

The exact source-to-target map is from M to its conjugate-paired spectrum.
Trace moments are preserved; the Frobenius norm supplies a one-sided
energy estimate whose equality defect detects nonnormality. The spectrum
loses edge locations and nonuniform boundary contacts. Therefore this
new universal norm estimate does not replace THM-4465's contact quadratic
form, infer a preserved coordinate cone along a trajectory, or supply a
pressure law, Euler regularity/blowup statement, LRC result, or Hamiltonian
inequality.

## 4. Independent exact finite controls

The referee's
[standalone source](../../04-computation/continuing15_20260909_universal_cubic_audit.py)
performs exact integer and rational computations. It imports no primary
producer, performs no floating-point spectral comparison, and does not
use a numerical eigensolver. The inequality is tested as
`3 tr(M^3)^2 <= E^3`, which is equivalent to the absolute bound for E>=0.

The declared complete general-real universes are:

| Order | Entry alphabet | All matrices | Exact trace-square-zero matrices |
|---|---|---:|---:|
| 1 | -2,-1,0,1,2 | 5 | 1 |
| 2 | -2,-1,0,1,2 | 625 | 41 |
| 3 | -1,0,1 | 19,683 | 2,333 |

The filter includes diagonal squares and cancellations of opposite signed
entries. This bank therefore exercises the stronger matrix theorem beyond
the oriented cone, rather than only replaying the consumer's hypothesis.

A separate complete signed-oriented bank permits each unordered pair to
be absent or to carry either sign of a unit weight in either direction.
The distinct matrix counts for orders one through four are
`1,5,125,15625`, totaling 15,756. These are different from the primary's
nonnegative two-amplitude universe despite the coincident counts.
Direct matrix products are compared with independently enumerated signed
cyclic and transitive triple products. Nonnegative members also test each
step of the production inequality, squaring F only after checking F>0.

Another complete bank takes the regular cyclic five-vertex orientation
`i->j` when `j-i mod 5` lies in `{1,2}`, and independently assigns each
of its ten edges weight one or two: all 1,024 choices. Its five cyclic
triangles have no common edge. Every nontrivial proper vertex subset is
checked and fails the homogeneous-module condition, so no nontrivial
constant-contact substitution decomposition can account for this support.
It is thus an explicit carrier outside the incoming common-edge and
recursive substitution grammar. The all-dimensional proof, rather than
these 1,024 instances, supplies its unrestricted amplitude bound.

The 67 named exact rational controls include:

- positive and negative equality cycles, the incoming six-vertex
  multi-return equality, and a nonzero nilpotent matrix;
- 25 normal real/complex block matrices with parameter
  `t=j/4`, `-12<=j<=12`, a real eigenvalue 2t, and a two-by-two block
  with real part `-(1-t^2/2)` and imaginary magnitude `1+t^2/2`.
  They satisfy trace-square zero exactly and exercise either sign of
  the nonreal real part, negative real eigenvalues, and nonzero first
  trace;
- 18 rational Householder conjugates of an equality cycle embedded in
  orders three, four, and six, preserving equality even outside the
  nonnegative coordinate cone;
- 18 exact unipotent similarities in the same orders, preserving cubic
  trace three while strictly increasing Frobenius energy above three;
- two disjoint equality cycles, and the signed production hostile above.

Every enumerated nonzero equality is independently checked for normality,
rank three, zero first trace, and its exact spectral polynomial. For a
positive equality, `r=tr(M^3)/E` is rational for these rational matrices,
and the source checks `M^4=r^3 M`. Negative equality is checked after
negation. These polynomial tests accompany the analytic complete spectral
criterion; they do not replace it. Equality occurrence counts across the
overlapping banks are seven zero, 49 negative, and 68 positive, not a
count of distinct matrix orbits.

## 5. Freeze and acceptance

Each normal and optimized run passes **77,236 always-active gates**.
All gates raise explicitly and survive Python -O. The
[output](continuing15_20260909_universal_cubic_audit.out) and
[certificate](continuing15_20260909_universal_cubic_audit_certificate.json)
are deterministic raw LF artifacts; normal and optimized bytes agree.
The certificate retains primary executable pins, exact universe counts,
the named rational matrices, equality controls, and scope.

After relocation, reproduce with

    python 04-computation/continuing15_20260909_universal_cubic_audit.py
    python -O 04-computation/continuing15_20260909_universal_cubic_audit.py

Outside the repository the certificate is written beside this source;
under `04-computation` it goes to `05-knowledge/results`. The primary
artifacts are read only. Frozen referee SHA256 values are:

    Python 26dd0ed15a354c2d9dfed53339dc5b01b5af3e3a9c35516ac2933c665805e3cf
    output 735e41b88ae74764a875259cca656d636a0f6260ac5ac09eabe8cd49ff4843d2
    cert   9ddb893ea03822185cdd142203a10b08850729ee46cbe02b5c6bd6ad9e56cd8a

The primary theorem, equality iff, and unrestricted weighted oriented
consequence are accepted. The incoming restricted theorem remains a
correct historical supplier; its former arbitrary-support research gap
can now be routed to this accepted proof without rewriting its original
scope or claiming external novelty.
