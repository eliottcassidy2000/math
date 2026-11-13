# Reusing the next-target singleton witness

Status: **PROVED** physical extension laws and word-family reductions;
**FINITE-EXACT** the ten-case divisor gate. This is an exact-source statement,
not all-target coverage.

For the singleton funnel, with integer tail address x written least-significant
bit first, the actual indegree-one ready source has

\[
H=Q(x)=3x+g(x)+(x\bmod2),\qquad
q=x+g(x)+(x\bmod2),\quad s=q-1,\quad F=H+1.
\tag{1}
\]

The shell energy g is defined in [DEFECT_TWO.md](DEFECT_TWO.md), equation(2).
The literal construction is the same funnel and newborn, with singleton
core and its unique entry. All positive addresses end in1 and give strong
sources. The parent packet supplies

\[
x_0=45157534,\quad g(x_0)=23780093,\quad
Q(x_0)=159252695=1215\,2^{17}+215.
\]

Its word is `01111001001100001000110101`, and q=68937627.

## 1. Positive infinite extension, keeping both named ports

As in the defect-two construction, prefixing a zero inserts one vertex
before the old tail, preserving all old arcs, the physical ready root r,
and its sole predecessor p. Since x₀ is even,

\[
\begin{aligned}
H_\ell&=180630136\,2^\ell-21377441,\\
q_\ell&=90315068\,2^\ell-21377441,\\
s_\ell&=q_\ell-1,\quad F_\ell=H_\ell+1
\end{aligned}
\qquad(\ell\ge0).
\tag{2}
\]

Thus the exact next target lies on a physical infinite ready family. Applying
NR to any member gives the primitive balanced pair `(3H_ℓ+2,H_ℓ)`.
The coefficient of (2) has odd part22578767, so (2) itself cannot match the
desired coefficient1215 after reindexing powers of two.

A more target-directed companion has addresses

\[
y_\ell=320\,2^\ell+44\quad(\ell\ge1).
\]

It is also unit-ready with both ports preserved as ℓ increases, and

\[
\boxed{H(y_\ell)=1215\,2^\ell+167,\qquad
q(y_\ell)=575\,2^\ell+79.}
\tag{3}
\]

The entire desired family is therefore exactly48 above this explicit
same-coefficient ready family. An actual uniform +48 repair preserving the
ready port would close that target family; this note does not supply that
repair. The constant48 is a physical-response target, not permission to add
counts abstractly.

## 2. Complete fixed-boundary, single-growing-run gate

This gate includes arbitrary changes to the fixed words on either side of
the growing run; it is not limited to edits of x₀. Write
`R_r(d)=rev_r(d)` for reversal of r bits and put
`f(n)=n+1−(n mod2)−g(n)`.

For a growing zero run, the address is

\[
x_\ell=2^r2^\ell n+d,
\]

and, for ℓ≥1,

\[
Q(x_\ell)=(4\,2^r-R_r(d))n\,2^\ell+Q(d)-f(n).
\tag{4}
\]

For a growing one run, the address and count are

\[
\begin{aligned}
x_\ell&=2^r(2^\ell(n+1)-1)+d,\\
Q(x_\ell)&=(4\,2^r-1-R_r(d))(n+1)2^\ell
 +g(n)+Q(d)-(4\,2^r-1-R_r(d)).
\end{aligned}
\tag{5}
\]

In (5), the low fixed word is normalized to have highest bit0, or to be
empty. In (4) its highest bit is1, or it is empty. These normalizations merely
absorb equal adjacent bits into the growing run. Also absorb the trailing
zeros of n in (4), or the trailing ones of n in (5); the identities
`f(2n)=f(n)` and `g(2n+1)=g(n)` preserve the constant. Thus n is odd in(4),
and n+1 is odd in(5).

For nonempty low words, the coefficient factor in parentheses is now odd,
lies between `3·2^r+1` and `4·2^r`, and must divide1215 if the count is to
match `1215·2^k+B` after any fixed index shift. Hence r≤8. The empty low word
is handled separately (factors4 and3). There is no cutoff assumption about
the original fixed words: normalization and the coefficient factor force
this bound.

Factoring `1215=3^5·5` leaves the following complete possibilities for B:

| growing bit | possible offsets B |
|---|---|
|0|−1199, −33, −11,167,339|
|1|339,65,15,−191,−263|

The receipt records all ten literal `(r,d,n)` tuples. None has B=215.
Consequently no singleton funnel with fixed boundary words around one
constant-bit growing run supplies the desired family, even only on an
arithmetic progression of sufficiently large k. Restricting the run length
to an arithmetic progression cannot change its constant term. Of these
complete offsets,167 is closest to215, giving (3).

For verification of the formulas, for r≥1 write

\[
g(2^r y+d)=g(y)+Ay+B(y\bmod2)+C,
\]

where `A=2^r−1−R_r(d)`, `B=1` exactly when the highest low bit is0, and
`C=g(d)−1`. Substitute `y=2^ℓn` or `y=2^ℓ(n+1)−1`, using the elementary
g recurrences. The empty low word follows directly from(1).

## 3. Why a fixed mixed periodic block does not remove the debt

This is a proof, not a periodic-word scan. Let a fixed mixed r-bit block of
value d be repeated, so `n_(j+1)=2^r n_j+d`. Put D=2^r−1. The same identity
shows that, after the first step, the energy is an exponential plus a
linear term `μj`, where

\[
\mu=C+B(d\bmod2)-Ad/D.
\tag{6}
\]

The linear coefficient is strictly positive for every mixed block. To see
this, let c_k count cyclic pairs of bits `0,1` separated by k positions,
and let a=c₁ be the number of cyclic01 transitions. Expanding (6) gives

\[
\mu=a-\frac{\sum_{k=1}^{r-1}c_k2^{k-1}}{D}.
\]

Cyclic binary autocorrelation gives `c_k=c_(r-k)`, while a union bound over
transitions gives `c_k≤ka`. Therefore

\[
\mu\ge a\left(1-
\frac{\sum_{k=1}^{r-1}(r-k)2^{k-1}}{D}\right)
=\frac{ar}{D}>0.
\]

A fixed outer word changes the observation to `K n_j+g(n_j)` plus bounded
parity and constant terms, so this linear drift survives. It cannot equal a
pure dyadic exponential plus constant on an infinite arithmetic progression
of block repetitions. The constant blocks are precisely the zero-drift
cases already covered by Section2.

This excludes one fixed periodic middle block with fixed words on both
sides. It does not exclude several coupled variable blocks, moving-position
repairs, changing cores, or a growing carry controller.

Reproduction:

```
python3 05-knowledge/results/target_directed_cores_20261112/constructors/singleton_continuation.py
```

The script checks the base target, both positive extension formulas, and the
complete normalized divisor gate. The proof of (6) needs no search.
