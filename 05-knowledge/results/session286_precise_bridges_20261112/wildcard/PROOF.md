# Exact spin fibers, homometric phase, and physical phase breaking

**Status: PROVED for the stated finite objects and all-order identities;
FINITE-EXACT for the displayed literal tournaments. Chan–Pak remains OPEN.**

This lane imports no conclusion from thermodynamics or phase retrieval merely
because their equations look similar. Its positive result is an exact map from
cyclic homometry to physical fusion twins, together with a one-arc operation
that releases their hidden phase as a full-central count difference. The
Johnson exchange operator gives exact reconstruction and moment identities.
It is not a scalar-target existence theorem.

## Inheritance and distant sources

The underlying quadratic attachment/fusion surface is inherited from
[Session285 tournaments](../../session285_five_web_20261112/tournaments/PROOF.md),
which credits THM-4097/4114 and the Session264 fusion carrier. That surface,
not a new statistical model, supplies physical realizability.

The fixed-cardinality mean and harmonic quadratic variance already belong to
[THM-4123: balanced-cardinality ear average](../../../../01-canon/theorems/THM-4123-balanced-cardinality-ear-average-growth-and-johnson-layer-lattice.md)
and [THM-4127: Johnson-slice Hoeffding variance](../../../../01-canon/theorems/THM-4127-johnson-slice-hoeffding-variance-and-central-support-dominance.md).
Our vertex-transitive formulas are their pure degree-two specialization, not
new mean/variance theorems. We include elementary derivations to fix gauges
and make the labelled-layer inverse and physical homometry bridge standalone.

The remote source of the present question is
[THM-2439: cyclic marker replica degree and homometric Gram boundary](../../../../01-canon/theorems/THM-2439-cyclic-marker-replica-degree-and-homometric-gram-boundary.md),
a July LRC investigation with two non-dihedral `(13,4,1)` difference sets.
Its warning is that complete Fourier power does not choose a marked phase.
[THM-2802: normalized unit bispectrum](../../../../01-canon/theorems/THM-2802-normalized-unit-bispectrum-and-projective-cyclic-orbit-reconstruction.md)
recovers an abstract cyclic orbit under its nonvanishing hypotheses, but
explicitly withholds physical attachment. Here Boolean root cuts are *already*
physical; the question becomes which lawful arc intervention exposes phase.
No bispectrum-to-tournament realization theorem is assumed.

Concept board: fixed-magnetization energy / marked local field / full-central
fusion carrier / cyclic homometry / labelled versus scalar observations.
The canonical hostile is a pair-transitive Paley core: exponentially many
labelled root cuts can have exactly one scalar count at each cardinality.

## 1. Exact fixed-magnetization spin model

Let D be a vertex-transitive tournament on n>=3 vertices. Let h=H(D), and
Q=sum_{S subset V(D)} H(D[S])H(D[V-S]), with H(empty)=1. Define

    B_ij=sum_{S subset V(D)} end_i(S) start_j(V-S),
    W_ij=B_ij+B_ji (i!=j), W_ii=0, r=(Q-2h)/n.

The path counts in B are literal directed paths on their indicated subsets.
The inherited attachment law says that, on adding a new root z with outgoing
neighborhood C of size k,

    F(C)=H(D+z)=h+rk-sum_{i<j in C}W_ij.                (1)

Every k-subset has identical noncentral carrier coordinates

    (h,Q,U,V,d,e)=
    (h,Q,h+k(Q-h)/n,h+(n-k)(Q-h)/n,h(1-k/n),hk/n).      (2)

Thus (1) is the *actual Hamiltonian-count spectrum in a full-central fiber*,
not a relaxation. Every row sum of W is 2r. With spins s_i=2*1_C(i)-1,

    F(C)=h+nr/4-(s^T W s)/8,    sum_i s_i=2k-n.        (3)

This is an exact fixed-magnetization Ising energy. We use no continuum limit,
Gibbs hypothesis or independence approximation. Boltzmann weights exp(beta F)
correspond to an antiferromagnetic sign when beta>0 and W>=0. All choices
of beta are optional observers of the same finite physical set.

## 2. A pure Johnson exchange mode and exact tomography

For n>=4, 2<=k<=n-2, put

    mu=h+r*k*(n-k)/(n-1),
    b_ij=W_ij-2r/(n-1)  (i!=j), b_ii=0,
    E(C)=mu-F(C)=sum_{i<j in C}b_ij.

Each row sum of b is zero. On the Johnson graph of k-subsets, write

    Lf(C)=sum_{a in C,b notin C}[f(C-a+b)-f(C)].

**Theorem 1 (exact exchange law).**

    L(F-mu)=-2(n-1)(F-mu).                            (4)

Consequently simple random exchange gives, for every starting cut C,

    E[F(C_t)] = mu + (1-2(n-1)/(k(n-k)))^t*(F(C)-mu).  (5)

The expression is an expectation, not a physical count or monotone path.

Proof: for a fixed pair {i,j}, its indicator is lost in 2(n-k) exchanges
when both lie in C, and gained in k-1 exchanges when exactly one lies in C.
Sum with coefficients b_ij. The sum across the cut is -2E(C), since row sums
are zero. Therefore LE=-2(n-k)E-2(k-1)E=-2(n-1)E. Averaging one random
neighbor proves (5).

**Theorem 2 (single-layer labelled tomography; inherited variance specialization).**

    b_ij = [sum_{C: |C|=k, {i,j} subset C}(mu-F(C))]
              / binom(n-4,k-2),                     (6)

    Var(F(C)) = k(k-1)(n-k)(n-k-1)
                 /[n(n-1)(n-2)(n-3)] * sum_{i<j}b_ij^2. (7)

Uniform C is used in (7). Formula (6) reconstructs the entire weighted
interaction from one **labelled** cardinality layer, together with h,Q.
An unlabelled histogram does not supply the sums in (6).

Proof: in (6), the selected pair contributes binom(n-2,k-2)b_ij;
pairs meeting it once contribute -2 binom(n-3,k-3)b_ij; disjoint pairs
contribute binom(n-4,k-4)b_ij. Their sum is binom(n-4,k-2)b_ij. Binomial
coefficients outside their range are zero. For (7), the sum of products of
distinct incident edges is -sum b^2, and that of disjoint edges is
(sum b^2)/2. Expanding E(C)^2 gives (p2-2p3+p4)sum b^2, where
pt=(k)_t/(n)_t, which simplifies to (7).

These exact spectral/moment rules are powerful observers of a given fiber.
Neither mixing, positive variance, nor a large variance excludes internal
holes. Even the full expected trajectory need not take an integer value.

## 3. Cyclic homometry becomes physical fusion twins

Let n be odd and let D be any Cayley tournament on Z/nZ: choose exactly one
of each nonzero pair {d,-d} as its arc increments. Translation invariance
makes W_ij=w_(j-i), with w_d=w_(-d), w_0=0. For a cut C define

    R_C(d)=|C intersection (C-d)|.

Then

    F(C)=h+rk-(1/2)sum_d w_d R_C(d).                 (8)

Two subsets A,B are homometric when R_A=R_B. They have the same cardinality,
so (2) and (8) prove:

**Theorem 3.** Any cyclic homometric A,B give rooted physical tournaments
with the same **entire seven-coordinate fusion packet**, for every Cayley
tournament on that labelled cyclic group.

This is a sufficient equality criterion, not a classification of all packet
collisions and not a claim the rooted tournaments are isomorphic. The group
of actual core automorphisms must be retained when claiming nonisomorphism.
The theorem does not require non-dihedral inequivalence.

There are systematic all-order supplies. In the group ring, if P+Q and P-Q
have unique representations (all relevant sums distinct), then the Boolean
sets A=P+Q and B=P-Q have generating functions a=pq and b=p*qbar; hence

    a*abar = p*pbar*q*qbar = b*bbar.

The uniqueness hypothesis is essential: without it this proves a multiset
identity and need not construct Boolean cuts. For example, any two finite
integer sets can first be separated at mixed radix, then embedded in a
sufficiently large odd cycle. Nontriviality/nonisomorphism still needs its own
check; two-point factors can give only a translated reflection.

## 4. A common arc toggle converts phase into a central count difference

Let A,B be homometric and write delta=1_A-1_B. Choose a vertex v at which
both cuts agree, with common bit epsilon. Toggle its membership in both.
The resulting sets A',B' have equal cardinality k+1-2epsilon. The inherited
responses therefore still agree. Directly from (1),

    F(A')-F(B')=-(1-2epsilon)*(W delta)_v.             (9)

**Theorem 4 (physical phase breaker).** A homometric pair yields a nontrivial
full-central switch by one common root-arc reversal if and only if its local
potential difference (W delta)_v is nonzero at an agreeing site. It yields
a unit odd-count switch iff that difference has magnitude 2.

This is a literal operation on actual tournaments. It is not allowed to
subtract two counts and declare the difference itself realizable. Subsequent
root fusion pays the inherited product of outside deletion counts.

If W is constant off the diagonal, W delta=-w delta, so all agreeing sites
have zero local difference. Thus the pair-transitive Paley hostile blocks
this operation perfectly, despite many homometric presentations. This is the
precise boundary of the analogy, not a failed numerical search.

**Prime-order separation lemma.** For an odd prime n, the rational circulant
matrix W of a Cayley tournament is nonsingular. Therefore all labelled local
fields W1_C recover C exactly, even when its energy loses phase.

Proof: W has positive row sum, so its constant eigenvalue is nonzero. If a
nonconstant eigenvalue vanished, the rational polynomial sum_{d=0}^{n-1}w_d t^d
would vanish at a primitive nth root, hence would be a scalar multiple of
1+t+...+t^(n-1). Since w_0=0 this scalar is zero, contradicting positive row
sum. The cyclotomic irreducibility here follows, for example, by applying
Eisenstein at n after substituting t+1.

This lemma does **not** guarantee a separating *agreeing* site: the Paley
example has all separation supported on disagreement sites. Nor does it
supply unknown membership bits when signed local derivatives are inferred
from flips. Those are explicit observer sidecars.

## 5. Common-background replica cancellation is exactly additive

The one-toggle law extends without approximation. Let J be all sites where
A,B agree and epsilon_v their common bits. For any R subset J, toggle R in
both cuts and put

    g_v=-(1-2epsilon_v)(W delta)_v,  s_v=1-2epsilon_v.

Then

    F(A triangle R)-F(B triangle R)=sum_{v in R}g_v.       (10)

The individual two energies still contain pair interactions between toggles;
these cancel *between the two replicas*, since the same change vector y is
used. Algebraically, the two new quadratic terms y^T W y cancel and only
-y^T W delta survives. Both cuts are physical and their new cardinality is
k+sum_{v in R}s_v. Thus the exact Laurent polynomial

    P(u,z)=product_{v in J}(1+u^(s_v) z^(g_v))            (11)

has a nonzero coefficient at u^t z^d iff this common-toggle construction
supplies a full-central switch with count difference d at cardinality k+t.
All coefficients are nonnegative; there is no cancellation-based inference
of existence. This is a precise bridge to cardinality-constrained subset
sum. It does not equate the gap d to an attainable path count, and it does
not make the two individual absolute counts linear in R.

There is also an exact particle-hole reciprocity. Toggling every agreeing
site sends(A,B) to(V-B,V-A). Regular row sums in(1) give F(V-C)=F(C), so
sum_J g_v=0, while sum_J s_v=n-2k. Therefore

    P(u,z)=u^(n-2k) P(u^(-1),z^(-1)).                  (12)

This relates opposite physical gaps at complementary final cardinalities.
It is another exact constraint, not a proof that any nonzero gap is small.

This supplies a positive cancellation example below: two large local
responses cancel to44, much smaller than either single-toggle response.

## 6. Exact construction and hostile controls

Reuse the old LRC sets A={0,1,3,9}, B={1,2,5,7} in Z/13Z. Both have R(0)=4
and R(d)=1 for d!=0. In the carousel D with arc increments1,...,6,

    h=3711175, Q=84542744, r=5932338,
    (w1,...,w6)=(1027338,1021076,1011596,991884,960716,919728).

The two 14-vertex rooted tournaments both have F=21508189. Adding the
common missing vertex8 to their root cuts gives

    F(A union {8})=23572029,
    F(B union {8})=23562137,         difference9892.

The full noncentral packet remains common. These numbers are independently
checked by literal subset Hamilton-path DP, not just (8) or (9). This is a
new physical use of the inherited homometric sets, **not a +2 switch**.

The exploration script also examined 32 specified order13 Cayley orientations
(all choices with the sixth increment fixed negative), one representative
per dihedral cut class at cardinalities3,...,6, and all translations of one
member of every retained homometric pair. Its smallest nonzero released gap
was338. This is a scoped negative experiment, not a lower bound over all
homometric pairs, orientations, orders or interventions. Search records are
kept separately from proof receipts. No count coverage follows from it.

A second order13 core with increments {12,2,10,9,8,7} has homometric
cuts with masks83 and1128. Common toggles at vertices2 and11 have local
gaps -2390 and2434, whose sum is44. The resulting cuts are2135 and3180,
with literal14-vertex counts24389195 and24389151. Their common noncentral
packet is(h,Q,U,V,d,e)=(3683797,84190184,40840591,47033390,1983583,1700214),
independently checked in `replay.json`. This illustrates exact interacting-background cancellation,
not a unit switch or a gapless interval. The exploratory double-replica
search also records its chosen finite universe separately.

## 7. Connection ledger and next exact question

| Map | Preserved | Lost | Required sidecar / decisive test |
|---|---|---|---|
| Root cut -> fixed-magnetization energy | Actual count and full-central fiber | Other rooted resources, isomorphism | Literal core and cut labels |
| Labelled layer -> Johnson eigenmode | Exact mean, variance and pair interaction | Scalar target nonemptiness | Labelled layer; equation(6) |
| Cyclic cut -> autocorrelation | All Cayley quadratic energies, full packet | Fourier phase | Local fields or higher observables |
| Homometric twin -> common toggle | Full-central compatibility | Zero count difference | Equation(9), actual agreeing site |
| Local field -> cut at prime order | Labelled cut exactly | Physical target choice | Nonsingular W; membership/sign convention |

An exact productive next task is to construct a family of physical Cayley
cores and homometric pairs with (W delta)_v=2 at an agreeing site, while
controlling deletion count and renewal. Theorems3–4 give a finite certificate
for every proposed family. They do not presently prove such a family, and
even that would not by itself anchor every absolute odd target.
