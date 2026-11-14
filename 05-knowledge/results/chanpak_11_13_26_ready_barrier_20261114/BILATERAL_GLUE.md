# Lengths 11, 13, and 26: steer, reset, and mix two literal sources

A single target line asks one ready source to supply all of the
remaining coordinates. There is an exact **two-source** alternative:
join two independently prepared rooted tournaments. The old
Tribonacci kernels still apply, including to sources impossible
inside the terrace. At length 26 the two transverse choices add
independently, without a product of the two choices.

## 1. A literal splice beyond the terrace

For a rooted tournament \(A\) write \(H_A,P_A,S_A\) for its
ordinary, positive rooted, and signed rooted count, respectively,
in the notation of the
[arbitrary-root bridge](https://github.com/eliottcassidy2000/math/blob/bd8f3030aee3e75af8b67695a76a857e0094a18a/05-knowledge/results/chanpak_literal_prefix_compiler_20261114/ARBITRARY_ROOTED_BRIDGE.md).
Use the converse of the right rooted tournament \(B\), and make
all arcs from \(A\) to this right side except the arc from its
marked vertex back to the marked vertex of \(A\). Then
\[
                 H_{\rm joined}=H_AH_B+P_AP_B-S_AS_B .
\]
A crossing odd directed cycle must use the unique backwards arc.
Cutting it there gives a rooted starting path in \(A\), followed
by the reverse of a rooted starting path in \(B\). Their total
number of vertices must be odd. The only invalid forward seam
would join the two singletons at the marked vertices, but that
uses two vertices and is even. The leftover odd-cycle packings
multiply across the two sides, and the factor 2 in the
[odd-cycle packing formula](https://arxiv.org/abs/2307.05569)
gives precisely the unsigned product minus the signed product.

For ready sides \(z_L=(a,b)^t,z_R=(c,d)^t\), substitute
\((H,P,S)=(a+b,a+b,-a)\) and the corresponding right signature.
The count is \(z_L^t C_1 z_R\). For \(q\ge1\), first put
\(q-1\) zero bridges on the left and then do the splice. The
same signed transfer and the tribonacci recurrence give
\[
                             H_q=z_L^t C_q z_R.
\]
For verification the first three results are \(C_1,C_2,C_3\).
An ordinary ordinal sum without a backwards arc corresponds
to \(q=0\), with \(C_0\) the all-ones matrix. For strong prepared
sides and \(q\ge1\) the final graph is strong. Thus this is
not just a bilinear identity for pre-existing terrace words.

## 2. The three different roles of 11, 13, and 26

Put \(w=(47,56)^t,t=(56,-47)^t,Q=5345\). The longitudinal and
transverse coordinates are \(U=w^t z,V=t^t z\).
Exact transfer values give
\[
 C_{11}(1,2)^t=47w,\qquad C_{13}=ww^t,\qquad C_{26}t=9w.
\]
At 11 a *particular physical ray* is steered to the retained
axis; the transfer is nonsingular. At 13 the whole transverse
axis is killed. At 26 that axis is instead sent to the retained
one, and the map is again nonsingular: \(\det C_{26}=-81\).
So the doubled length 26 is not a second reset.

More strongly, for the two-sided interaction we have the
all-order identity
\[
                           t^t C_q t=T_{q-30},
\]
where \(T_0=0,T_1=T_2=1\) and the tribonacci sequence is extended
to negative indices. Both sides follow the same recurrence;
\(q=30,31,32\) gives \(0,1,1\), then propagate in both directions.
The exact zero theorem of
[Mignotte–Tzanakis](https://doi.org/10.4064/aa-57-4-357-364)
(also stated in [Bilu et al.](https://arxiv.org/pdf/2401.06537))
makes the zero-interaction lengths exactly
\[
                              13,\quad26,\quad29,\quad30.
\]
Lengths 26, 29, and 30 are nonsingular with
\(C_qt=\ell w\) for \(\ell=9,56,103\), respectively.
For two physical source pairs each differing by \(2t\),
independently switch either input. The mixed second
difference of the four ordinary counts is exactly
\(4T_{q-30}\). Length 26 cancels this mixed term without
killing either individual adjustment; length 27 gives \(-4\)
and length 28 gives \(+4\). This is a two-source
construction with a signed interaction you can choose
exactly, not just a heuristic independence statement.

## 3. A new target inverse and its true limitation

At \(q=26\), since \(w,t\) are orthogonal and both have
square norm \(Q\), the literal count is
\[
 H=\frac{14736371U_LU_R+9(U_L V_R+V_L U_R)}{Q}.
\]
The exact lattice condition is \(V\equiv2958U\bmod Q\);
its sign is opposite to the perpendicular vector used in
the earlier Leech note. Write \(V=2958U+2Qj\). The extra
factor 2 is forced by \(a\) odd and \(b\) even. This simplifies
the absolute endpoint equation to
\[
                      H=2767U_LU_R+
                           18(U_Lj_R+U_Rj_L).
\]
The coefficients are genuine integers, and the \(j_Lj_R\)
term vanishes. For fixed \(U_L,U_R\), solve one ordinary
linear Diophantine equation for \(j_L,j_R\), then test
physical ready-state membership on each side. If actual
sources of fixed \(U\) are already indexed, this gives an
immediate complement lookup between the two banks rather
than a search through full tournaments.

There is a genuine limitation even in the relaxed cone.
At fixed \(U\), the allowed first coordinates are spaced by
112 and satisfy \(1\le a\le U/159\), because \(b\ge2a\).
There are at most \(1+U/17808\) transverse choices.
After removing the baseline and dividing by 18, the right
coordinate disappears modulo \(U_L\). Thus even if
\(\gcd(U_L,U_R)=1\), a **fixed** pair of longitudinal values
has at most \(1+U_L/17808\) phases modulo \(U_L\), far less
than full coverage. The exact inverse makes clear what
must actually be changed: the longitudinal sources as
well as the two transverse choices. This is a more flexible
route, not a universal selection theorem.

## 4. An actual nonterrace prime through length 26

Use the five-vertex rooted core \((2,20,9,3,13)\), root 1
(bits index the outneighbors), with \(H=11,R_+=13\).
Its width-one bridge gives \((11,26)\), so \(U_L=1973\).
Use on the right either of the two real reset-free states
\((5701,13636)\) and \((5589,13730)\) from the previous
[changed-source construction](https://github.com/eliottcassidy2000/math/blob/bd8f3030aee3e75af8b67695a76a857e0094a18a/05-knowledge/results/chanpak_literal_prefix_compiler_20261114/CHANGED_EQUAL_U.md),
both with \(U_R=1031563\).
The two strong 50-vertex splices at length 26 have exact
ordinary totals
\[
                       5611326183575,\quad 5611326148061 .
\]
The second is prime by trial division. Their difference is
\(18\cdot1973=35514\), independently verifying the retained
right source and the left control. In other literal examples
where *both* sides are varied, the replay checks the zero
mixed correction by independent interval packing. The small
splices are checked directly by Hamilton subset DP.

Run `python bilateral.py`.
