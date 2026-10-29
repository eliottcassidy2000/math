# An explicit carry-tail budget under source-cylinder refinement

Status: **PROVED** for the `(d,r)=(3,3)` auxiliary first-extinction rule of
MOMENTS.md. It is a bound on a specified source measure, not a pointwise
Collatz descent theorem.

## 1. Rational supersolutions give an exponential carry tail

Here `p=1/7`, `q=1/21`, so

\[
 M=6/7+tM^3/7,\quad L=20/21+tL^3/21,
 \quad B={L/3\over1-(t/7)(L^2+LM+M^2)}.              \tag{20}
\]

At `t=4/3`, the rational numbers `6/5` and `21/20` are supersolutions
for the first and second cubics. Iteration from zero, equivalently summing
finite trees, proves

\[
 M(4/3)\le6/5,\quad L(4/3)\le21/20.
\]

Consequently `L^2+LM+M^2<=1521/400`, and

\[
 1-{4\over21}{1521\over400}\ge {193\over700}>0,
 \qquad B(4/3)\le {245\over193}.                     \tag{21}
\]

Formula (20) at this `t>1` is justified by finite-tree sums: the two
supersolutions bound their masses and slopes; (21) then bounds the carry
sum geometrically. Thus no analytic continuation or assumed convergence
has entered the proof.

For every integer `h>=0`, nonnegative coefficients imply

\[
 P(m\ge h)\le{6\over5}(3/4)^h,\quad
 E[\Lambda; m\ge h]\le{21\over20}(3/4)^h,\quad
 \boxed{E[B_w; m\ge h]\le{245\over193}(3/4)^h}.        \tag{22}
\]

A sharper mass-only estimate follows directly from (13):

\[
 P(m\ge h)\le {147\over50(2h+1)}(243/343)^h.          \tag{23}
\]

Indeed `C_m=binom(3m,m)/(2m+1)<= (27/4)^m/(2m+1)`, since the corresponding
binomial probability at success probability `1/3` is at most one.
Multiply by `(6/7)^(2m+1)(1/7)^m` and sum the geometric majorant. The
carry estimate (22), rather than (23) alone, pays the hidden coordinate.

## 2. Keeping a marked source interval

Let `I_K` be **any** odd residue class modulo `2^{K+1}`, of odd-Haar mass
`2^{-K}`. All expectations below use the same original source, conditioned
to lie in `I_K`. They do not restart the process with an independent source.
For any nonnegative `X`,

\[
 E[X\mid I_K]=2^K E[X1_{I_K}]\le2^K EX.              \tag{24}
\]

Set the branching threshold to `3K+h`. Then (22) and
`2(3/4)^3=27/32` give

\[
 \boxed{E[B_w; m\ge3K+h\mid I_K]
 \le {245\over193}(27/32)^K(3/4)^h}.                 \tag{25}
\]

The same formula holds for probability or slope mass, with constants
`6/5` and `21/20`. It is uniform over the residue defining `I_K`, over
the refinement depth `K`, and over `h>=0`. In particular each newly
retained source bit is paid by three additional branching nodes, or nine
additional odd steps since `J=3m+1`. The actual valuation clock `D` stays
in every word and source guard; `3K+h` is a size budget, not a relabelling
of `D`. No uniform deterministic bound on `D` follows from this statement.

## 3. The exact remaining boundary

Every single positive integer has Haar measure zero. Uniform conditional
tail bounds at a growing threshold can hold while one point remains in
all unresolved sets. For example the binary cylinders of the all-zero
sequence have vanishing measure but a nonempty intersection. Thus (25)
cannot be promoted to universal entry by compactness or by shrinking
`I_K` to its marked integer.

Moreover first-extinction still need not be descent, even for this cheaper
rule: source `3` has first valuation `1`, so stops after `3 -> 5`.
The repaired target is to add a source-sensitive endpoint condition
`Lambda*n+B_w<n` and a rule for what happens on its failure. The recurrence
now exposes and bounds a genuine carry debt; a return controller must
preserve (24)'s source rather than sample a fresh Haar source after failure.

The **OPEN** question sharpened by (25): is there a port-indexed renewal
extension whose failure return satisfies a strict weighted contraction
on these same conditional source measures, with a pointwise arithmetic
certificate eliminating the possible null set? Neither part is supplied
by the exponential tail of the auxiliary grammar alone.

## 4. Why this coarse argument pays three branches per bit

The exact coefficient ratio in (13) for `(3,3)` tends to `243/343`, so
the mass generating function has convergence radius `R=343/243`. This
ratio follows directly by cancelling factorials in consecutive coefficients.
Any exponential-moment proof based only on (24) and `M(t)<infinity` needs
`2/t^b<=1` for a uniform bound after paying `b` branching nodes per
retained bit (and strict inequality for decay with refinement depth).

But

\[
 R^2={117649\over59049}<2,
 \qquad 2\cdot243^2-343^2=449>0.
\]

Thus no allowed `t` in this particular coarse proof can pay only two
branches per bit, whereas `t=4/3` pays three. This proves minimality of
the **integer budget within this exponential-moment/conditioning argument**.
It is not a lower bound for the actual conditional distributions: a
source-sensitive improvement to (24) could do better. The appearance of
`243=3^5` and `343=7^3` comes from the ternary shape count and the `1/7`
branch probability, not from an identification with another 243-object
carrier.
