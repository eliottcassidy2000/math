# Matching bandwidth and parity-sized content certificates

Session 197 consequences of `UNIVERSAL_CONTEXTS.md`. These statements concern the proved positive single-module operator, not a new spectrum-completeness claim.

## 1. The S196 matching rank is exactly a shift bandwidth

Let (Q,r) have N exterior vertices. Put O={v:r->v}, I={v:v->r}, and let nu be the matching number of the bipartite graph with edges u->v from O to I. This is exactly the cross graph in the pinned S196 cloning proof cited by the companion packet.

The universal context operator A_Q has a finite Laurent-shift representation

    A_Q=sum_(j=-nu)^nu p_j(x) S^j,

and its highest forward coefficient p_nu(x) is a nonzero polynomial with nonnegative integer coefficients. Thus its exact maximum shift is nu. This is a full-profile refinement of the scalar cloning rank nu+1.

**Proof.** In every crossing cycle skeleton, the outside segment between two successive marked clone-runs starts in O and ends in I. It contains an O-to-I transition. Choosing one such edge from each segment gives a matching because all these segments, including those in different skeleton cycles, have disjoint exterior vertex sets. Hence the total number of marked clone-runs is at most nu. Each L_(t,epsilon) has shifts at most t in absolute value, so the upper bound follows from the skeleton formula.

Conversely choose a maximum matching u_i->v_i. Its nu disjoint crossing triangles (r,u_i,v_i), each with one marked clone-run, form a valid skeleton packing. Leave all other exterior vertices unused. This contributes exactly

    X^(N-nu) delta^nu.

In every skeleton contribution the highest forward shift has coefficient +1 in each L factor. Multiplication by its nonnegative power of X and 2^c_0 leaves that coefficient positive. No skeleton of smaller total run count reaches shift nu, and those reaching it have no negative leading contribution. Thus p_nu is nonzero and coefficientwise nonnegative. QED.

The relation delta^2=Delta^2+4Delta also puts the entire operator in positive normal form

    sum_(i,j) a_ij X^i Delta^j + sum_(i,j) b_ij X^i delta Delta^j,

with nonnegative integers a_ij,b_ij and j<=nu in the first sum, j+1<=nu in the second. This is a constructive finite operator description for each fixed rooted template. It is not a fixed-dimensional state space for unbounded input size.

In particular, matching rank one gives

    A_Q=P(X)+Q_1(X)delta+R(X)Delta,

with all three coefficient polynomials nonnegative and P(x)=B_(Q-r)(x). At scalar x=1 its positive clone ray is

    H(Q[r:Tr_m]) = alpha*2^m+beta,   m>=1,
    alpha=Q_1(1)+R(1),
    beta=P(1)-2R(1).

The offset is odd because P(1)=H(Q-r) is odd, including the empty deletion convention. The zero-base term at m=0 must not be discarded if one wants deletion; this is the same positive-m restriction as S196. A unit offset requires H(Q-r)-2R(1)=+1 or -1; positivity alone does not force it.

## 2. Exact content in half as many bin evaluations

Let g(T)=gcd_(d>=1) B_T(d), the inherited full-profile content. For nonempty T of order n, let u_k be its central-factorial state. Put

    gamma_k^odd=(2k+1)!,
    gamma_k^even=(k+1)(2k+1)!.

Then

    g(T)=gcd_k gamma_k u_k
        =gcd_(1<=d<=ceil(n/2)) B_T(d),                         (1)

using the gamma sequence corresponding to n's parity.

**Proof.** For every integer d,

    q_k(d)/(2k+1)! = binom(d+k,2k+1),

and

    r_k(d)/[(k+1)(2k+1)!]
       = binom(d+k+1,2k+2)+binom(d+k,2k+2).                    (2)

These are integer-valued polynomial identities; generalized integer binomial coefficients cover negative d as well. At d=k+1 both normalized values are one, and they vanish at d=1,...,k. Consequently the evaluation matrix at d=1,...,s+1, where s is the highest state index, is a unit lower-triangular integer matrix on the weighted coordinates gamma_k u_k. Its inverse is integral. The two finite gcds agree, and every later evaluation is an integer combination of the same weighted coordinates. Since s+1=ceil(n/2), (1) follows. QED.

This gives a parity-sized certificate for the S193 profile-content obstruction without changing the obstruction's scope. It does not say that content one supplies an inverse construction.

There are also universal scalar congruences for every fixed whole-module context:

    H(C_Q[T]) = H(C_Q[K1]) H(T)       (mod 6),  if |T| is odd,
    H(C_Q[T]) = H(C_Q[Tr_2]) H(T)     (mod 12), if |T| is even.

Indeed every nonzero-height q_k is divisible by 6 at every integer, and every nonzero-height r_k is divisible by 12; the skeleton operator evaluates integer combinations of integer shifts at x=1. These congruences do not apply to a nonuniform operation that splits the entire input module.

## 3. Reproduction

Run `python3 verify_consequences.py` and its `python3 -O` version. The checker compares shift support from the independently enumerated skeleton operator with an independent augmenting-path matching number at all 285 rooted labelled templates through order four. It checks content via bin evaluations and weighted coordinates for all 1,099 labelled tournaments through order five. The all-order statements follow from the proofs above; finite tests only control implementation and normalization.
