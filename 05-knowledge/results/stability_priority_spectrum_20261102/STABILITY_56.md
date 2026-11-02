# Exact weighted tournament stability from order 56

Status: **PROVED** from the independently accepted Session139 tensor
matching theorem and the explicitly credited graph inequalities below.
Same-agent review and exact certificates of this final extension; it has
not yet received an independent-agent or formal audit.

## Statement and source chain

For every real weighted tournament A of order n>=56, let

    Delta=H(A)-det(I+A)-|det A|,
    V=sum_{i<j} A_ij A_ji, R=floor(n/2),
    kappa_n=(4/R)[1-(3/4)^R-1_{n even}(1/4)^R].

Then **Delta>=kappa_n V sharply**. At V>0 equality is exactly the ordinal
sum of R fair pairs and, in odd order, one singleton. At V=0 equality is
exactly transitivity. Together with accepted orders3..9, the remaining
unresolved orders are **10..55**. This closes1,444 of the originally open
orders10..1499; it does not claim that55 is an obstruction.

This is a synthesis of distinct results, not a claim of sole priority:

1. Session138 supplies binary rounding and the equality classification.
2. Session139 supplies the per-deleted-edge tail and local-dual mechanism.
3. This packet's466 and320 proofs retain -BR and free the local dual.
4. The incoming [discarded-degree lemma](../stability_priority_moduli_20261102/STABILITY_N500.md)
   gives W>=2D^2/n, yielding our [238 proof](STABILITY_238.md).
5. The later, independently accepted [tensor matching theorem](../session139_incidence_moduli_and_stability_20261102/tensor_matching/PROOF.md)
   and its [140 bootstrap](../session139_incidence_moduli_and_stability_20261102/stability/PROOF_N140.md)
   supply a new nonlinear bound. Combining it with the stronger degree
   estimate and the flexible baseline-preserving dual yields56.

The threshold is derived from the inequalities. Its coincidence with other
appearances of56 in the repository supplies no additional proof.

## 1. Nonlinear bootstrap using the stronger degree estimate

Use the notation and full accepted-input list in [STABILITY_238.md](STABILITY_238.md):
B=E Delta(T), W the incident variance-product sum, and a descending-weight
greedy matching of mass M=V-D. In particular

    B>=W>=2D^2/n, W>=2V^2/n-V/4,
    Delta>=B+f(V), f(V)=(1-exp(-4V))/4,
    Delta>=B+kappa_n M.

The new tensor theorem additionally gives, for any retained matching,

    Delta>=B+1-product_i(1-q_i)-1_{n=2r}product_i q_i
         >=B+1-exp(-M)-4^(-R).                          (1)

The subtraction is present only for a perfect retained matching; otherwise
subtracting4^(-R) only weakens the lower bound. B is the full original
binary-rounding expectation, retained once. The theorem proves (1) by
grouping Pfaffian squares and their parity; no false pointwise bound on an
individual Pfaffian is used.

Assume V>0 and Delta<=kappa V. Since R>=28 and kappa<=1/7,
f(1)>49/200>1/7 excludes V<=1 (f(V)/V decreases). The first strict
inequality follows from exp(4)>16319/315>50. Thus

    kappa V>B+c0,             c0=49/200.

Whenever kappa V>B+c with c>0, complete the square using B>=2D^2/n:

    M=V-D>c/kappa+2D^2/(n kappa)-D
          >=c/kappa-n kappa/8
          >=cR/4-1-1/(2R)=:m_c(R).                     (2)

Here kappa<=4/R and n<=2R+1. The function m_c is strictly increasing in R.
Equations(1),(2) give the improved bootstrap

    kappa V>B+1-exp[-m_c(R)]-4^(-R).

Uniformly for R>=28, the following exact advances hold:

    49/200 -> 1/2 -> 91/100 -> 99/100.                   (3)

At R28 the three input mass bounds are122/175,139/56,7493/1400.
Let T12(m)=sum_{j=0}^{12} m^j/j!. For each displayed advance c->c',

    T12(m_c(28)) * [1-c'-4^(-28)] > 1.                  (4)

These are strict rational comparisons, checked by the companion script.
Because exp(m)>T12(m) and m_c increases while4^(-R) decreases, (4) proves
each step of (3) at every R>=28. Consequently all candidates satisfy

    kappa V>B+c,                  c=99/100.              (5)

## 2. Eight heavy edges and a paid tail

Set

    t=1/7, d=381/500, lambda=t^2/d=500/18669.

As before, B<=kappa D and V<=n/8+n kappa/2. Suppose L>=1 selected
weights are at most t, with total U. Retain h<=R-L other edges, of mass
M_h=V-D-U. Using (2),(5), U<=tR and q_i<=1/4,

    h>=4M_h>(c-4t)R-4-2/R=:g(R).

Now c-4t>0, g increases, and g(28)=2677/350>7. Hence h>=8 and
m:=R-L>=8. The accepted a_h matching bound, its monotonicity, and the
same baseline-preserving expansion as in the238 proof yield

    cL-4U<4epsilon V+4D-BR,       epsilon=(3/4)^(R-L).   (6)

With Q_R=R+33/2+8/R, we have4V<=Q_R. The discrete endpoint argument gives

    4epsilon V/L
       <=max{Q_R(3/4)^(R-1), Q_R(3/4)^8/(R-8)}<c-d.     (7)

Both bounds decrease for R>=28, exactly as in the238 proof with the
retained-edge floor changed to8. At R28 the tighter endpoint has slack

    (c-d)-[Q_28/20](3/4)^8=1751781/458752000>0.

The other endpoint is also strictly smaller than c-d. Therefore

    dL-4U<4D-BR.                                        (8)

Since d-4t>0, this implies D>0.

## 3. The same flexible dual gives a strict contradiction

For light groups0<=w<=t, S>=0, the exact identity is

    wS/2+S^2/4-tS+lambda(d-4w)
       =(S/2-t+w/2)^2+w(t-4lambda-w/4)>=0.

Its nonnegative remainder follows from3t/4-4lambda=1/74676>0.
Use it on the disjoint light groups and W_i>=w_i S_i>=tS_i on heavy
groups. Equation(8) gives

    B>=W>(t-4lambda)D+lambda BR,
    t-4lambda=667/18669>0.                              (9)

If1-lambda R<=0, (9) is impossible. Otherwise B<=kappa D<=4D/R gives

    (t-4lambda)D < B(1-lambda R)
                  <=(4/R-4lambda)D
                  <=(t-4lambda)D,

since R>=28=4/t. This is again impossible, including at R28. Hence L=0.
All selected weights exceed t>=4/R>=kappa. Assigned wedges force D=0,
and B<=kappa D gives B=0. The accepted matching-support theorem and H=1
classification now prove the sharp bound and complete equality locus.

## Arithmetic meaning and audit boundary

The independent [Eisenstein proof of the tensor lemma](../session139_incidence_moduli_and_stability_20261102/stability/EISENSTEIN_BERNSTEIN_TENSOR_PROOF.md)
is a genuine arithmetic input: odd-top/even-rest integer coefficient
vectors evaluate nontrivially in F4=Z[zeta]/2, so their Eisenstein norms
are positive integers. The full derivative tensor retains the energy
that an individual Pfaffian can lose through cancellation. This arithmetic
parity statement supplies (1); it is not an arbitrary-modulus norm theorem
or a Riemann-hypothesis transfer.

`verify_stability_56.py` checks the bootstrap by positive rational Taylor
sums, all endpoint signs, the local dual polynomial identity, and a finite
replay at56..1499. The preceding graph and tensor theorems provide the
structural inputs; finite replay alone is not the proof. Remaining10..55
requires further native incidence work or a substantially sharper certificate.
