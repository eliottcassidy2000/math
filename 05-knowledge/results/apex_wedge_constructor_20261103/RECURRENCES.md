# Infinite families of Hamiltonian unit jumps

**PROVED scoped, with exact finite linear certificates.** These families
contain two actual strong tournaments differing in one apex edge and in
exactly two Hamiltonian paths. Their graph and two-move native program are
given by A2 and A4 of [THEORY.md](THEORY.md).

## R1 Three closed families

Put P=2^r+c,R=P+2. Let L_c(r)=H(T-) be the lower count; H(T+)=L_c(r)+2.
The full formulas are
\[
\begin{aligned}
L_{-2}(r)&=-7-\tfrac{45}{2}2^r+4\,3^r-\tfrac{19}{2}4^r+
                 \tfrac{46}{3}6^r &&(r\ge2),\\
L_0(r)&=\tfrac{269}{4}4^r+\tfrac{315}{4}2^r+5-45[\,r=2\,]
                 &&(r\ge2),\\
L_3(r)&=\tfrac{4303}{16}4^r+\tfrac{7455}{8}2^r+623-263[\,r=3\,]
                 &&(r\ge3).
\end{aligned}                                                    \tag{R1}
\]
The source orders are 2r+5 for c=-2 and 2r+6 for c=0,3.

| c | first lower counts | annihilator in the shift E |
|---|---|---|
| -2 | 339, 2625, 17397, 109749 | (E-1)(E-2)(E-3)(E-4)(E-6) |
| 0 | 1351, 4939, 18481, 71389 | E(E-1)(E-2)(E-4) |
| 3 | 25027, 84381, 305835, 1161831 | E(E-1)(E-2)(E-4) |

The exceptional initial terms are essential. They are the zero eigenmode
of the all-index recurrence, not evidence against the tail formula.
The scalar recurrence orders are respectively 5,4,4: the displayed distinct
exponential modes have nonzero coefficients. For the last two families,
discarding the first term gives a tail of order three.

Writing t=r-r0 with r0=2,2,3, the generating functions are
\[
\begin{aligned}
\sum_{t\ge0}L_{-2}(2+t)z^t
 &=\frac{-7}{1-z}-\frac{90}{1-2z}+\frac{36}{1-3z}
   -\frac{152}{1-4z}+\frac{552}{1-6z},\\
\sum_{t\ge0}L_0(2+t)z^t
 &=\frac{1076}{1-4z}+\frac{315}{1-2z}+\frac5{1-z}-45,\\
\sum_{t\ge0}L_3(3+t)z^t
 &=\frac{17212}{1-4z}+\frac{7455}{1-2z}+\frac{623}{1-z}-263.
\end{aligned}                                                    \tag{R2}
\]
Add 2/(1-z) for the upper counts.

## R2 Why the certificates prove all indices

The gap-partition evaluator is inherited from the October20 feedback-word
work and the [two-apex decoder](../session134_constructive_spectrum_frontier_20261101/two_apex/THEORY.md).
For each order of the apexes, a path consists of three increasing base gaps.
A gap is empty, nonempty with a bad right endpoint, or nonempty with a good
right endpoint. A first vertex must satisfy the gap's left boundary guard.
The middle gap may be empty only if the apex edge agrees with that order.
These rules give a path-partition bijection, not a weighted surrogate.

Each family has a literal column template
prefix A^(t+a0) middle B^(t+b0) suffix, listed in
[recurrences.json](recurrences.json). Binary expansion proves these templates
for every t>=0: only the lengths of the two constant runs change.

For each apex order, [recurrences.py](recurrences.py) computes an exact
left Krylov space for A and right Krylov space for B. The dimensions are
3 by 2 for c=-2 and 2 by 3 for c=0,3. Their tensor product gives six
coordinates per apex order, twelve in total. This proves an all-t formula
h_t=alpha C^t beta.

For each polynomial p in the table, the certificate checks
alpha C^j p(C) beta=0 for j=0,...,11 using exact rational arithmetic.
Cayley-Hamilton then gives the identity for every j>=0. R1 satisfies that
same recurrence and the necessary initial values. Thus the certificate is
a finite proof of an all-order statement, rather than a recurrence fitted
to a list. The expanded graph checker independently counts small sources,
including the order16 example below.

## R3 A source-preserving 7p lane and its exact phase

For c=3,r>=4, write U_r=L_3(r)+2. Then
\[
 U_r=(4303\,4^r+14910\,2^r+10000)/16 .
\]
Modulo 7 this vanishes exactly when r=2 mod3. Set r=5+3t, t>=0:
\[
 U_r=275392\,64^t+29820\,8^t+625.                         \tag{R3}
\]
At r=5 the two strong order16 sources have counts
\[
 305835,\qquad 305837=7\cdot43691,
\]
and 43691 is prime by trial division through its square root.
The named-port program uses two native moves with peak order17.

Let Q(t)=U_(5+3t)/7. Directly,
\[
 Q(t)=4+3t\pmod7.
\]
More generally binomial expansion of 8=1+7 gives, for every e>=0,
\[
 Q(t+7^e)-Q(t)=3\,7^e\pmod{7^{e+1}}.                     \tag{R4}
\]
Indeed 8^(7^e)=1+7^(e+1) modulo 7^(e+2); substituting into R3 and dividing
by 7 leaves coefficient 2*275392+29820=3 mod7.

Consequently there is exactly one residue t modulo 7^e for which
7^(e+1) divides U. The next phase is found by testing its seven lifts.
In particular 49 divides U exactly when t=1 mod7, or r=8 mod21.
The saved certificate constructs the unique lifts through e=7.

This is an exact count-divisibility clock for specified Boolean sources.
It is a different object from the earlier mod49 bulk-return gate.
The parameter r, the bit words, both physical apex ports and the native
move clock are all retained. R3 gives prime-quotient examples and exact
repeated-prime exclusions; it does not assert infinitely many prime
quotients or coverage of every required 7p.
