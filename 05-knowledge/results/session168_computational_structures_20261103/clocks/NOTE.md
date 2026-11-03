# A four-lap ceiling for Cayley–Dickson clocks

**Status: PROVED explicit constructions and period classifications; FINITE-EXACT controls.**

The same construction joins four precise objects: signed basis transports in the entire Cayley–Dickson tower, the tribonacci guard register, a nearest-neighbor cellular clock, and S167's constant path-cover-profile tournament compiler. The result is an arbitrary-size family with exact periods

\[
26\,3^a,\qquad52\,3^a,\qquad104\,3^a,
\]

all hidden from the complete path-cover profile. Four laps is a sharp universal ceiling for the additional return time contributed by **one fixed-level signed-unit transport register**. Increasing the Cayley–Dickson level does not increase this ceiling.

The finite-state twist machinery has classical predecessors, notably Bales's twist trees. The value here is the explicit convention, proof, local-clock coupling, sharp period theorem, and inherited tournament realization. No literature-priority claim is made. None of these constructions makes Cayley–Dickson multiplication associative or identifies the encoded state graph with the blue/black complement-line metagraph.

## 1. Inheritance and the objects being connected

- **Closest mechanisms:** [S167 octonion transport](https://github.com/eliottcassidy2000/math/blob/596ae1ea47f640689304af59237e5f737c75c0b8/05-knowledge/results/session167_recurring_structures_20261103/octonion/NOTE.md) replaces nonassociative products by honest permutations; [S167 cubic norm](https://github.com/eliottcassidy2000/math/blob/596ae1ea47f640689304af59237e5f737c75c0b8/05-knowledge/results/session167_recurring_structures_20261103/norm/NOTE.md) identifies the literal guard map with multiplication by a cubic unit; [S167 ternary compiler](https://github.com/eliottcassidy2000/math/blob/596ae1ea47f640689304af59237e5f737c75c0b8/05-knowledge/results/session167_recurring_structures_20261103/ternary/TERNARY_CONSTRAINTS.md) embeds ternary state words into one complete path-cover-profile fiber.
- **Hostile example:** the octonion word \(L_3^{-1}L_1L_2\) fixes the identity unit but is not the identity permutation. A return test on one sheet does not establish path independence.
- **Corrected near miss:** the matrix of an external state evolution is not automatically the physical operation on its tournament encoding. Here only the modular register is a tribonacci guard register; the encoded tournaments have a constant profile.
- **Least-used sidecar:** the spatial marker and signed unit are retained explicitly. Forgetting them loses the multiplier that separates the three periods.

Our concept board is: doubling signs / signed-XOR transport / cubic field / nilpotent spatial shift / complete profile / marked local observation. The source-to-target map is a local finite-state evolution followed by the explicit S167 word encoding. It preserves state identity and the stated period, destroys the state information when only the profile is retained, and requires the fixed marker and decoding convention as sidecars.

## 2. One six-state machine computes all doubling signs

Use the same convention as S167:

\[
(a,b)(c,d)=(ac-\bar d b,\; da+b\bar c),\qquad
\overline{(a,b)}=(\bar a,-b).
\]

Write \(e_x\), \(0\leq x<2^n\), for the standard basis of the level-\(n\) algebra, with the newest doubling bit most significant. Then

\[
e_a e_b=(-1)^{\omega_n(a,b)}e_{a\mathbin\oplus b}.
\]

The following machine reads paired digits \((a_j,b_j)\) **least significant first**. Its state is

\[
(w,v,A,B)=(\omega_j(a_{<j},b_{<j}),\omega_j(b_{<j},a_{<j}),
 [a_{<j}\ne0],[b_{<j}\ne0]).
\]

Starting at \((0,0,0,0)\), append bits \((h,k)\) by

\[
\begin{aligned}
w'&=(k?v:w)+hk+hB,\\
v'&=(h?w:v)+hk+kA,\\
A'&=A\lor h,\qquad B'=B\lor k.
\end{aligned}
\]

The first two equations are over \(\mathbb F_2\); the notation \(k?v:w\) means choose \(v\) when \(k=1\), otherwise \(w\). At the end the output is \(w\).

**Theorem 1 (uniform minimal twist machine).** Exactly the following six states are reachable, and six states are necessary for a deterministic final-output machine reading this digit convention at arbitrary lengths:

| State | Tuple | Output | 00 | 01 | 10 | 11 |
|---|---|---:|---:|---:|---:|---:|
| 0 | (0,0,0,0) | 0 | 0 | 1 | 2 | 5 |
| 1 | (0,0,0,1) | 0 | 1 | 1 | 4 | 3 |
| 2 | (0,0,1,0) | 0 | 2 | 3 | 2 | 4 |
| 3 | (0,1,1,1) | 0 | 3 | 4 | 4 | 4 |
| 4 | (1,0,1,1) | 1 | 4 | 3 | 3 | 3 |
| 5 | (1,1,1,1) | 1 | 5 | 4 | 3 | 5 |

**Proof.** The four high-bit cases in the displayed doubling rule give respectively

\[
\omega(a,b),\quad\omega(b,a),\quad\omega(a,b)+[b\ne0],
\quad\omega(b,a)+1+[b\ne0],
\]

which proves the update and its invariant by induction. The six-state table is closed, and all its states are reached by words of length at most two. Minimality follows from pairwise distinguishability: different outputs need the empty suffix; among output-zero states the suffixes for pairs \(01,02,03,12,13,23\), respectively, can be taken as

\[
10,\quad(01,01),\quad01,\quad10,\quad01,\quad10.
\]

The two output-one states 4 and 5 are distinguished by suffix 11. Thus all six reachable prefixes must receive different machine states. QED.

This is a bounded-memory computation of a **basis product sign**, not a six-state summary of an arbitrarily long path's action on all initial units. S167 already proves that even the octonion transport completion needs 128 different all-sheet summaries. Bales's [twist-tree exposition](https://jwbales.us/) and [periodicity paper](https://arxiv.org/abs/1602.02843) provide the standard antecedent; their basis numbering and doubling convention must not be silently substituted for ours.

## 3. A tower-independent transport ceiling

At a fixed level \(n\), let

\[
U_n=\{(-1)^s e_x:s\in\mathbb F_2,\ x\in\mathbb F_2^n\}.
\]

Every left multiplication by a signed basis unit is a permutation of this set. Retain permutations and their compositions; do not collapse them to a parenthesized algebra product.

**Theorem 2 (four-lap theorem).** Every permutation \(g\) generated by these left transports has the form

\[
g(s,x)=(s+f(x),x+a)
\tag{1}
\]

for a vector \(a\in\mathbb F_2^n\) and a Boolean function \(f\). Consequently \(g^4=1\). More precisely:

1. If \(a=0\), only cycles of lengths one and two occur. If \(z=|\{x:f(x)=0\}|\), their multiplicities are \(2z\) and \(2^n-z\).
2. If \(a\ne0\), there are no fixed points. Put \(h(x)=f(x)+f(x+a)\). The multiplicities of two- and four-cycles are
   \[
   |\{x:h(x)=0\}|,\qquad\tfrac12|\{x:h(x)=1\}|.
   \]

**Proof.** Each generator has (1) by the XOR support law. Composition preserves (1), since

\[
(a,f)\circ(b,q)=(a+b,\;q(x)+f(x+b)).
\]

Squaring gives \(g^2(s,x)=(s+h(x),x)\), so squaring again is identity. For \(a=0\), inspect each sign pair. For \(a\ne0\), inspect each unordered index pair \(\{x,x+a\}\): its four signed states form two two-cycles if \(h(x)=0\), otherwise one four-cycle. QED.

The same argument covers inverses. It does not require alternativity, associativity, division, or classification of the full multiplication group. It concerns signed **basis** units, not multiplication by arbitrary algebra elements.

If a finite clock has an orbit of exact period \(m\) and each clock step applies one such transport to a retained unit register, the transport accumulated over one lap is a single \(g\). Every lifted orbit over that clock orbit therefore has period \(m\), \(2m\), or \(4m\). Both the index/sign state and the order of operations matter.

The ceiling is sharp at every level \(n\ge1\): left multiplication by \(e_1\) has the four-cycle

\[
e_0\longmapsto e_1\longmapsto-e_0\longmapsto-e_1\longmapsto e_0.
\]

The unchanged octonion hostile control has eight fixed states and four two-cycles, so the same base lap can close on some sheets and require two laps on others.

## 4. The tribonacci register supplies a radius-one clock

The inherited physical guard map is

\[
C(s,u,d)=(u,s+u+d,s).
\]

Reduce it modulo three. Let

\[
K=\mathbb F_3[\alpha]/(\alpha^3-\alpha^2-\alpha-1).
\]

The cubic has no root in \(\mathbb F_3\), so \(K\) is a field of 27 elements. The explicit linear identification

\[
\zeta(s,u,d)=(\alpha^2-\alpha)s+\alpha u+d
\]

intertwines \(C\) with multiplication by \(\alpha\). Its coefficient map is invertible. Exact reductions give

\[
\alpha^{13}=1,\quad\alpha\ne1,\qquad
\beta:=\alpha+1,\quad\beta^{13}=-1,\quad\beta^2\ne1.
\]

Thus \(\alpha\) has order 13 and \(\beta\) has order 26. The checker verifies these small field identities directly.

For \(k=3^a\) sites on a periodic row, define

\[
(\Psi X)_i=\alpha X_i+X_{i+1}.
\tag{2}
\]

Equivalently, on three ternary tracks \(v_i=(s_i,u_i,d_i)\),

\[
v'_i=Cv_i+v_{i+1}\pmod3.
\]

This is a nearest-neighbor linear cellular automaton on a finite ring. Its finite-ring map is bijective; no finite-radius inverse on the infinite lattice is asserted.

**Theorem 3 (complete clock periods).** The map \(\Psi\) has exact order \(26k\), and the delta configuration \(X^*\), with value one at site zero and zero elsewhere, has that period. More generally, identify configurations with

\[
R=K[z]/(z^k-1)=K[t]/(t^k),\quad t=z-1,
\]

where multiplication by \(z\) acts as the shift in (2). If a nonzero configuration has \(t\)-adic order \(v\), its exact period is

\[
26\,3^{\lceil\log_3(k-v)\rceil}.
\tag{3}
\]

There is one zero fixed point, 26 states of period 26, and for \(1\le j\le a\), exactly

\[
27^{3^j}-27^{3^{j-1}}
\]

states of period \(26\,3^j\). Dividing each state count by its period gives the orbit count.

**Proof.** In the displayed ring, (2) is multiplication by \(\beta+t\). Write a nonzero configuration as \(t^v U\), with \(U\) a unit. A period \(N\) first requires \(\beta^N=1\), hence \(26\mid N\). The polynomial \((\beta+t)^{26}-1\) has order exactly one at \(t=0\), because its linear coefficient is \(26\beta^{25}\ne0\). Raising to a \(3^j\)-th power multiplies that order by \(3^j\); raising further to a power prime to three does not change it. Thus the least possible \(N\) is (3). The delta is the unit 1, so \(v=0\) and its period is \(26k\). The number of elements divisible by \(t^{k-b}\) is \(27^b\), proving the counts. QED.

For example, on three sites all \(27^3=19,683\) states split into one fixed point, one 26-cycle, and 252 cycles of length 78. The checker exhausts this entire state space independently of the valuation proof.

## 5. A local trigger attains all three return multipliers

Mark site zero permanently. A marker can be implemented as an unchanged binary track, so the update rule remains spatially homogeneous on the enlarged alphabet. Fix a Cayley–Dickson level n≥1 and keep one signed-unit register at the marked site. It changes only when the current clock value at that site is one.

**Lemma 4 (uniform occupation at the marked site).** Over one \(26k\)-step orbit of \(X^*\), every nonzero element of \(K\) occurs at site zero exactly \(k\) times; zero never occurs there.

**Proof.** Write time as \(qk+r\), with \(0\le q<26\), \(0\le r<k\). Characteristic three and \(S^k=I\) give

\[
\Psi^k=(\alpha^k+1)I=\beta^k I.
\]

For \(r<k\), no shifted delta term other than the zero shift contributes at site zero, whose value in \(\Psi^rX^*\) is \(\alpha^r\). Its value at time \(qk+r\) is therefore \(\beta^{qk}\alpha^r\). Since \(\gcd(k,26)=1\), varying \(q\) gives each nonzero field value once, for each fixed \(r\). QED.

Choose one of three local register rules, applied whenever that local clock value is one:

- leave the signed unit unchanged;
- negate it;
- multiply it on the left by \(e_1\).

Start the register at \(e_0\). The clock has \(k\) triggering times per lap, and \(k\) is odd. Thus the one-lap register permutations have, on this initial state, orders one, two, and four respectively. The complete marked systems have exact periods

\[
26k,\qquad52k,\qquad104k.
\tag{4}
\]

The rule is nearest-neighbor on the clock track and local on the marker/unit tracks. A completely uniform finite alphabet is \(K\times U_n\times\{0,1\}\), with all unmarked unit registers fixed at \(e_0\). Only one register varies in the prescribed invariant family. The resulting finite-ring map is a bijection: recover the previous clock using its invertible finite map, then undo the indicated local permutations. The universal four-lap theorem applies to any other schedule of basis transports driven by the clock.

## 6. The complete path-cover profile cannot see these clocks

Recall the inherited S167 strong six-vertex seeds \(A_0,A_1,A_2\), masks \(153,177,181\). They are pairwise nonisomorphic, pairwise adjacent as one-arc-flip **classes**, and share labeled-slot polynomial

\[
F_A(q)=q^6+12q^4+16q^2.
\]

For a ternary word \(w\in\{0,1,2\}^d\), the ordinal sum

\[
T_w=A_{w_1}\oplus\cdots\oplus A_{w_d}
\]

is injective on isomorphism classes and satisfies

\[
F_{T_w}(q)=F_A(q)^d,\qquad H(T_w)=29^d.
\]

Encode each clock cell by its three \(\mathbb F_3\) coordinates and the retained unit by its \(n\) binary index bits and one sign bit. This uses

\[
d=3k+n+1
\]

ternary coordinates, with the fixed marker convention known beforehand. Therefore each of the three systems in (4) gives an exact periodic family of pairwise distinct tournament classes of order

\[
6d=18k+6n+6,
\]

all sharing the **entire** path-cover profile \(F_A^d\), not merely their Hamilton count. The inherited strong wrapper gives the same period at order \(30d+10\), again with one constant full profile; its Hamilton count is common but need not be \(29^d\).

A synchronous clock step can alter several coordinates. It is not one tournament arc flip. The S167 induced Hamming-graph embedding gives a sequence of exactly as many class flips as changed coordinates, using an external schedule. Relabelings may be needed inside a block. The map does not claim an intrinsic autonomous rule on unmarked tournaments, or a physical realization of arbitrary modular register states as endpoint counts.

This is also an information obstruction: any phase estimator that receives only the complete profile is constant on the entire encoded orbit and cannot recover any of its phase. The states themselves remain distinct and explicitly decodable.

## 7. A small observation sidecar recovers what the profile loses

The scalar guard observable is \(\ell(s,u,d)=s+u\). For a clock trajectory, retain its whole spatial observation row \(Y_t=(\ell v_i(t))_i\). Then three consecutive rows recover the full clock configuration:

\[
Y_0=\ell v,\quad Y_1-SY_0=\ell Cv,\quad
Y_2-2SY_1+S^2Y_0=\ell C^2v.
\]

At each site the observation matrix is

\[
\begin{pmatrix}1&1&0\\1&2&1\\3&3&2\end{pmatrix},
\qquad \det=2\ne0\pmod3.
\]

Thus three spatial scalar snapshots are sufficient. These observables belong to the decoded modular registers. They are not Hamilton counts of the encoded \(T_w\), whose counts and every other profile coefficient are constant. The signed-unit sidecar remains necessary to distinguish its two- or four-lap refinement.

## 8. Reproduction and boundaries

Run:

    python3 check.py > results.json
    python3 -O check.py > results.optimized.json
    cmp results.json results.optimized.json
    python3 independent_matrix_check.py > matrix_results.json
    python3 -O independent_matrix_check.py > matrix_results.optimized.json
    cmp matrix_results.json matrix_results.optimized.json

The standard-library checker uses explicit exceptions, never removable assertions. It verifies every basis pair through level eight by independent low-digit transduction and high-digit pair recursion, all six reachable states and separating suffixes, 200 deterministic random transport words at each level one through eight, the exact permutation cycle formulas, and the inherited octonion hostile word. Random words are controls, not the proof of Theorem 2.

It separately implements the field and guard-coordinate CAs, checks delta periods and marked-site occupation for \(k=1,3,9,27,81\), verifies the four-lap lift using literal level-five basis permutations, exhausts all 19,683 three-site states, and tests three-snapshot observability. A separate script, importing neither the field implementation nor the main checker, builds the literal 3k-by-3k guard/shift matrix over F3. For k=1,3,9,27 it verifies exact matrix and delta orders, the unique fixed point, and every predicted fixed-space dimension by Gaussian elimination. The compiler and six-vertex seed are inherited from the linked audited S167 packet, not re-enumerated here.

No claim is made about a general nonlinear CA, the infinite-lattice reversible-CA classification, the full higher-level multiplication group, a physical blue/black metagraph lift, or the open tournament Hamiltonian-path spectrum conjecture. The ceiling bounds the extra return multiplier of a fixed transport register over a given periodic base; it does not bound the length of the base clock itself, nor systems that change doubling level during an orbit.
