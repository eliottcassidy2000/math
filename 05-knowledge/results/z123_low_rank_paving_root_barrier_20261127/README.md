# Ordinary Z is real-rooted for paving matroids through rank five

2026-11-27. This narrows the remaining minimum-rank question in the
[z123 investigations](../z123_exact_rank_threshold_20261127/README.md).
Every paving matroid of rank at most five has only negative real
zeros in its ordinary Z-polynomial, including nonsparse paving.
A lower-rank counterexample to real-rootedness must have rank four
or five and a nonpaving simplification. The latest
[sparse-paving result of Xie and Zhang](https://arxiv.org/abs/2609.07636)
covers sparse paving in arbitrary rank; the present result does
not assume the dual is paving.

The gamma discriminant has the following sharp minima, including
the convention that the discriminant expression for a linear gamma
polynomial is the square of its linear coefficient.

| Rank | Non-Boolean paving | Connected paving |
| --- | ---: | ---: |
| 4 | 9 | 28 |
| 5 | 28 | 60 |

The non-Boolean minima are attained by \(U_{r-1,r}\oplus U_{1,1}\);
the connected minima are attained by \(U_{r,r+1}\) and by no other
connected paving matroid of the corresponding rank.

## 1. Coefficient and discriminant test

Write \(h\) for the number of hyperplanes of a rank-\(r\) paving
matroid on \(n\) elements, and \(s\) for the largest hyperplane
size. For \(r=4,5\), the first nonconstant coefficients are
\[
a=h,\qquad b=\sum_H\binom{|H|}{r-3}-(r-3)\binom{n}{r-2}.
\]
Indeed, the [ordinary Z coefficient formula](../z123_exact_rank_threshold_20261127/README.md)
is \(b=W_{r-2}+f_{r-3,r-1}-f_{r-3,r-2}\). Every set of size at most
\(r-2\) is a flat, since otherwise it would be contained in a
circuit of size at most \(r-1\) after adding one point.
In particular \(f_{r-3,r-2}=(r-2)\binom{n}{r-2}\).

Set \(u=t/(1+t)^2\). The palindromic Z-polynomial is
\[
 (1+t)^r \left[1+(a-r)u+
 \left(b-(r-2)a+\frac{r(r-3)}2\right)u^2\right].
\]
By [Ferroni–Matherne–Stevens–Vecchi, Theorem 8](https://www.mat.univie.ac.at/~slc/wpapers/FPSAC2023/36.pdf),
the gamma coefficients are nonnegative. The gamma discriminant is
\[
D=h^2+(2r-8)h+r(6-r)-4b.
\]
Once \(D\ge0\), all actual gamma roots are negative. Such a root
\(u_0\) yields two negative real roots of \(t/(1+t)^2=u_0\), with
product one and sum \(1/u_0-2<-2\). A drop in gamma degree
contributes the appropriate copies of the negative real root
\(-1\). The same argument applies unconditionally in ranks up
through three, where the gamma degree is at most one.
Simplification preserves the flat lattice and hence ordinary Z.

## 2. The convexity certificate

Fix a largest hyperplane \(H\) and an outside point \(x\).
Each \((r-2)\)-subset \(Y\subseteq H\) yields a different exterior
hyperplane by taking the closure of \(Y\cup\{x\}\); its meet
with \(H\) is \(Y\). Thus \(h\ge1+\binom{s}{r-2}\).
The crude bound
\[
 b\le h\binom{s}{r-3}-(r-3)\binom{n}{r-2}
\]
makes \(D\) positive if \(h\ge4\binom{s}{r-3}\). It also gives
\(D\ge r(6-r)+4(r-3)\binom n{r-2}\) in that regime, stronger
than the sharp bounds in the table since \(n\ge r\).
The ratio \(\binom{s}{r-2}/\binom{s}{r-3}=(s-r+3)/(r-2)\)
shows that this handles \(s\ge5r-11\).

For all remaining parameters there is an exact coverage identity:
\[
\sum_H\binom{|H|}{r-1}=\binom n{r-1},
\]
since every \((r-1)\)-set is an independent set spanning a
unique hyperplane. Remove one largest block and let
\(h'=h-1\), \(m=\binom n{r-1}-\binom s{r-1}\). Necessarily
\(h'\le m\le h'\binom s{r-1}\). Consider the piecewise-linear
curve joining
\[
 \left(\binom{\ell}{r-1},\binom{\ell}{r-3}\right),
 \qquad r-1\le\ell\le s .
\]
It is concave: the slopes are
\(\binom{\ell}{r-4}/\binom{\ell}{r-2}\), decreasing in \(\ell\).
At mean coverage \(m/h'\), Jensen bounds the remaining blocks'
contribution to \(b\) by \(h'\) times the curve height.
When \(\binom{\ell}{r-1}\le m/h'<\binom{\ell+1}{r-1}\),
the bound for \(b\) is therefore
\[
\binom{s}{r-3}+h'\binom{\ell}{r-3}-(r-3)\binom{n}{r-2}
+\frac{\binom{\ell}{r-4}}{\binom{\ell}{r-2}}
 \left(m-h'\binom{\ell}{r-1}\right).
\]
At the upper endpoint use \(\ell=s\) with the last term zero.
This certifies every realizable size profile without pretending
every relaxed profile comes from a matroid.

## 3. Exact checks and sharpness

The range has \(r-1\le s<5r-11\),
\(1+\binom{s}{r-2}\le h<4\binom{s}{r-3}\) and
\(n\ge s+1\) with the coverage bounds above.
The integer cross multiplication in the replay checks 214 budgets
in rank four and 6,738 in rank five. The only zero lower bound
is \((s,h,n)=(r-1,r,r)\), the Boolean case. Otherwise the
lower bounds are at least 9 and 28 respectively. For
\(U_{r-1,r}\oplus U_{1,1}\) there is one hyperplane of size \(r\)
and \(\binom r{r-2}\) hyperplanes of size \(r-1\); inserting the
counts gives equality in the two non-Boolean bounds.

A connected paving matroid has \(n\ge s+2\), since otherwise
the only outside point of the largest hyperplane is a coloop.
In rank five the only budget in this subrange with Jensen bound
at most 60 is \((s,h,n)=(4,15,6)\), the uniform matroid.
In rank four the uniform budget \((3,10,5)\) has bound 28;
the other budgets with bound at most 28 are \((4,h,6)\) with
\(h=7,8,9\). If \(j\) of the hyperplanes have size 4, coverage
requires \(20=h+3j\). Moreover the complements of the 4-element
hyperplanes are pairwise disjoint pairs: two such hyperplanes
have intersection size at most 2. Thus \(j\le3\), ruling out
all three budgets. The two uniform examples give the stated
connected minima.

## 4. An obstruction beyond paving

The first-two-flag identity also gives a screen for arbitrary
matroids of rank \(r=4,5\). Let
\(m_*=\max_{H}W_{r-3}(M|H)\). For each rank-\((r-2)\) flat \(Y\),
write \(p_Y=W_{r-3}(M|Y)\). Then
\[
 b=\sum_H W_{r-3}(M|H)-\sum_Y(p_Y-1).
\]
The penalties are positive as a basis in \(Y\) gives at least
two rank-\((r-3)\) flats. Therefore \(h\ge4m_*\) suffices
for real-rootedness even without paving. Furthermore every
hyperplane \(H\) has \(h\ge1+W_{r-2}(M|H)\): extend each
rank-\((r-2)\) flat \(Y\le H\) by the same outside point
and use \(\operatorname{cl}(Y+x)\cap H=Y\).
Any counterexample of one of these two ranks must thus have a
maximizing hyperplane with
\[
 1+W_{r-2}(M|H)\le h<4W_{r-3}(M|H).
\]
This keeps exact incidence constraints, in the spirit of the
[previous Gram budget](../z123_exact_rank_threshold_20261127/README.md),
but is only a necessary screen, not a proof in the nonpaving sector.

Run the companion replay with Node.js; no floating-point root
or matroid enumeration is needed.
