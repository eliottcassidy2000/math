# Six necessary counting coordinates and a periodic-return obstruction

Status: PROVED by the exact arguments below, with FINITE-EXACT independent literal-graph controls. This is a refinement of an inherited compiler, not a proof of Chan-Pak. No theorem ID is reserved here.

## Inherited interface and scope

The unique-start rooted funnel grammar and six-coordinate update are inherited from [TARGET_RETURN, section 3](../../berggren_target_return_20261111/reset/TARGET_RETURN.md) and [ALTERNATING_AND_CAP, sections 2–3](../../berggren_target_return_20261111/reset/ALTERNATING_AND_CAP.md). Its exact code is [funnel_six_state.py](../../berggren_target_return_20261111/reset/funnel_six_state.py). A legal tail is empty, or ends in 1. An empty tail requires a nontrivial strong parent; every operation below starts from the rooted directed triangle, so this guard holds at every step.

The state is

\[
c=(h,e,s,p,f,t)=(P(1),E(1),S(2),P(2),E(2),S(3)).
\]

The ready child's **full** Hamiltonian count is \(R=h+f-1\). Thus the target \(R'=3R+2\) is exactly the linear equation \(L'=3L\), where \(L=R+1=h+f\). These are full counts, not a selected path sector. The triangle seed is \(c_0=(3,1,3,12,3,7)^T\). The scalar compiler does not recover the tournament's isomorphism type; its grammar address remains a separate sidecar.

Use four legal letters \(E=()\), \(B=(1)\), \(C=(0,1)\), \(D=(1,1)\). For a chronological word \(w=w_1\cdots w_k\), set \(T_w=T_{w_k}\cdots T_{w_1}\), and \(\ell=(1,0,0,0,1,0)\).

## 1. Six is minimal for linear and affine counting models

**Theorem.** Over the rational numbers, any finite-dimensional linear state model that reproduces \(R+1\) after every legal continuation of every state reachable from the triangle needs at least six coordinates. Even allowing affine transitions and affine readout for \(R\) still requires six state coordinates. The inherited compiler attains this bound.

Choose past words

\[
U=(\varnothing,E,B,C,D,EE)
\]

and future words

\[
V=(\varnothing,E,B,D,EE,EB).
\]

The reachability matrix \(X=[T_uc_0]_{u\in U}\) has determinant \(-34560\). The observation matrix \(Y=[\ell T_v]_{v\in V}\) has determinant \(-16\). Hence their actual past/future count matrix has

\[
\det(YX)=552960\ne0.
\]

In any d-dimensional linear realization this matrix factors through a d-dimensional space, so its rank is at most d. Therefore d≥6.

For the affine claim use \(U'=(E,B,C,D,EE,EB)\) and the difference columns \(X'=[T_uc_0-c_0]_{u\in U'}\). They have determinant \(-50688\), and

\[
\det(YX')=811008\ne0.
\]

An affine observation of the difference between two hidden states is linear in their state difference. Thus the same rank bound gives d≥6 even for affine models; subtracting the constant 1 in R has no effect. This is not a lower bound on arbitrary nonlinear encodings of integers, or on models restricted to a smaller family.

The certificate uses legal strong tournaments throughout. The accompanying receipt gives every matrix and the 36-entry observation table. Literal subset dynamic programming independently checks the state and actual ready-child Hamiltonian count for seven short graph constructions. The remaining table entries use the inherited all-order funnel identity.

## 2. Every periodic empty/one-tail program misses an infinite +2 return

**Theorem.** Fix any nonempty word W over \(\{E,B\}\), of length m. From any legal strong unique-start seed, repeat W indefinitely. At period ends let \(R_j\) be the ready-child full count. It is impossible that

\[
R_j+1=3^{mj}(R_0+1)\qquad\text{for every }j\ge0.
\]

In particular, no such periodic program gives \(R\mapsto3R+2\) at every letter. This includes changing the order and frequency of empty and one-vertex tails, and strengthens the earlier repeated-B cap obstruction. It does not exclude nonperiodic choices, larger tail words, port rewiring, or other graph operations.

The six-state matrices split into independent upper and lower triples. On the upper triple (h,e,s), change coordinates to \((e,h-e,s-h)\). The conjugate matrices are

\[
A_E=\begin{pmatrix}0&1&0\\1&1&1\\1&0&0\end{pmatrix},
\quad
A_B=\begin{pmatrix}1&1&0\\2&2&2\\0&0&0\end{pmatrix}.
\]

They are nonnegative and have induced column-sum norms at most 2 and 3, respectively. If W has a empty-tail letters and b one-tail letters, its upper spectral radius is at most \(2^a3^b\), strictly below \(3^m\) when a>0.

The lower triple (p,f,t) has matrices

\[
N_E=\begin{pmatrix}2&-2&2\\1&-1&0\\1&0&1\end{pmatrix},
\quad
N_B=\begin{pmatrix}4&-2&6\\2&-1&0\\2&0&3\end{pmatrix}.
\]

Their determinants are 2 and 12. Hence \(v_3(\det N_W)=b\). If a>0, then b<m, so the integer \(3^m\) cannot be an eigenvalue of \(N_W\): any integer root of its monic characteristic polynomial divides its constant coefficient. Neither block then has eigenvalue \(3^m\).

By Cayley–Hamilton, every scalar observation of \(T_W^jc_0\) satisfies the characteristic recurrence of \(T_W\). A nonzero exact sequence \(L_0(3^m)^j\) can satisfy this recurrence only if \(3^m\) is an eigenvalue. Since \(L_0=R_0+1>0\), this proves the case a>0. More strongly, a mixed program cannot agree with the claimed geometric sequence for all seven successive values j=0,…,6; those seven values already contradict the order-six characteristic recurrence.

It remains to exclude W=B^m. The matrix \(N_B\) is nilpotent modulo 3: \(N_B^3\equiv0\pmod3\). Thus the nonleading characteristic coefficients of \(N_B^m\) are divisible by 3, while its determinant has exactly 3-adic valuation m. Evaluating its characteristic polynomial at \(3^m\), the constant term has valuation m and every other term has strictly larger valuation. Consequently \(3^m\) is not a lower-block eigenvalue here either.

After k≥1 consecutive B letters, the upper readout is exactly

\[
h_k=3^{k-1}(h_0+2s_0).
\]

If the total L at period ends were geometric with ratio \(3^m\), then the lower readout \(f_{mj}=L_{mj}-h_{mj}\) would also be geometric for j≥1. Its coefficient is nonzero because \(f=E(2)>0\) for an actual tournament (take the root as a singleton first bin and a Hamiltonian path on its complement as the second bin). This contradicts the absence of the lower eigenvalue, again by Cayley–Hamilton.

## 3. What this says about six and about the return search

The rank-six result concerns the amount of information required by an observer. The order-six rational dynamics in the sibling [dynamics note](../dynamics/PROOF.md) concerns the order of a map. No identification between those two roles follows from their common integer.

The useful bridge is methodological: preserve the complete state on which the target depends, then test a proposed return on that state. Here the parent count has a pure 3-mode under B, but the ready count also observes a lower block with incompatible arithmetic. A successful construction must change that lower state, change the grammar, or use a nonperiodic mechanism. Neutral port switches supply a concrete candidate repair operation; mere parent-count agreement cannot replace them.

Reproduction: `python -B 05-knowledge/results/cycle_opening_six_return_20261112/compiler/six_state_audit.py`. The finite controls exhaust all 510 E/B programs of lengths 1–8, and all 64 legal fixed tails of lengths 0–6. Their negative outcomes support the proof but are not its quantifiers. No general exclusion for larger fixed tails is claimed.
