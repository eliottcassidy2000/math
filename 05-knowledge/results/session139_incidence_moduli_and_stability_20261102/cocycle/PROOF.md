# An integral four-face charge behind the old Pfaffian parity law

Status: PROVED, independently audited (integral and weighted identities). Historical source:
THM-1455-why-eleven-spectra.md and its distinct same-ID sibling
THM-1455-skew-pfaffian-parity-and-a-mod-16-law.md. No novelty claim.

## 1. Canonical signed lift

Let S be the skew sign matrix of an ordinary tournament. For distinct ordered
vertices i,j,k define alternating cochains

    t_ijk = S_ij S_jk S_ki,
    c_ijk = (S_ij+S_jk+S_ki+t_ijk)/4.

Here c is +1 or -1 on a directed triangle, according to its orientation,
and zero on a transitive triple. Let delta be the usual alternating
simplicial coboundary and put eta=delta c. Since the edge curl is itself
a coboundary, delta squared zero gives

    eta = (delta t)/4,       delta eta=0 over Z.

For a=S01,b=S02,c=S03,d=S12,e=S13,f=S23,

    eta_0123=(abd-ace+bcf-def)/4.

The four signed monomials in the numerator have product +1, so their sum
is -4, 0, or 4. Hence eta lies in {-1,0,1}. The four-point tournament
classification (or direct sign algebra) gives

    eta_0123^2 = (Pf(S_0123)^2-1)/8.

Thus eta is an integral signed lift of the old odd-fourset indicator.
It is invariant under Seidel switching because each triangle product t
is invariant. Under relabelling it transforms as an alternating cochain;
its support and squared norm are ordinary invariants.

On every ordered five-set,

    sum_j (-1)^j eta_(five-set minus vertex j) = 0.

This is integer conservation, not merely an even-parity constraint.

## 2. Every modulus, with the lift retained

For any integer M>=2, reducing eta modulo M preserves this exact cocycle
identity. Its nonzero coordinates remain nonzero because they are units
+/-1. Thus its support is the same at every modulus. For odd M, division
by 4 is legitimate in the displayed polynomial formula. For even M,
compute the integral lift before reducing; one cannot divide by 4 in
Z/MZ or recover eta from S modulo 2 alone (all sign matrices have the
same off-diagonal reduction modulo 2).

This statement does not identify the full tournament spectrum from one
modular cochain, nor turn a parity-only quotient into an integral invariant
without its sign data.

## 3. A sharp support gap

Suppose eta is not zero, and choose a four-set X in its support. For every
vertex v outside X, the five-face conservation law on X union {v} forces
at least one other nonzero four-face. Each such face contains its particular
v and otherwise vertices of X, so faces chosen for different v are distinct.
Consequently

    |support eta| >= 1+(n-4)=n-3.

This applies identically after reduction modulo every M>=2. An ordinal
sum consisting of one directed triangle and n-3 singleton blocks attains
n-3: its only cyclic triple is that triangle, and eta's support consists
of all four-sets containing it. Therefore, writing e4(S) for the sum of
fourth-order principal determinants,

    e4(S)=binom(n,4),
    or e4(S)>=binom(n,4)+8(n-3),

and the first positive gap is sharp for every n>=4. This supplies an
elementary explanation for the old n=7 absence of support size two.
It does not classify all larger support sizes or all equality switching
classes. The old statement that Paley7's k=14 means all foursets are odd
is false: there are 35 foursets, of which 14 are odd.

## 4. Weighted extension needs a second moment

For a weighted tournament A, let K=A-A^T, q_e=A_ij A_ji, and T be its
independent binary rounding. The mean charge is still

    E eta(T) = delta [K_ij K_jk K_ki]/4,

and remains a real cocycle. But its square is NOT (Pf(K)^2-1)/8.
For a fixed four-set, the three perfect-matching terms of Pf(S(T)) use
six distinct random edges and are independent. Their variance gives

    E[Pf(S(T))^2]
      = Pf(K)^2 + 3-sum_(three matchings {e,f}) K_e^2 K_f^2.

Therefore the exact probability of a bad four-face is

    E[eta(T)^2]
      = [Pf(K)^2+2-sum_matchings K_e^2 K_f^2]/8.

Globally, with V=sum q_e and M22=sum_disjoint_unordered_pairs q_e q_f,

    E |support eta(T)|
      = [e4(K)-binom(n,4)+4 binom(n-2,2)V-16M22]/8.

The edge-disjoint matching incidence M22 is the needed second-moment
sidecar. The integral conservation law survives averaging, while the
integer support gap is a statement about each binary realization, not
about the number of nonzero coordinates of its possibly cancelling mean.
This is distinct from the M23 sidecar restoring weighted Hamiltonian
counts from the full zeta function at order five.

## 5. A genuine second-moment loss on an isospectral fibre

Take the five-vertex skew matrix with only K03=a and K14=b nonzero above
the diagonal, and Q=I-ww^T/3, w=(1,1,-2,0,0). For |a|,|b|<=1, both
K and QKQ^T define valid weighted tournaments. The first-checkpoint zeta
proof establishes the common full adjacency/skew spectra and rooted
resolvent. Both matrices have bipartite skew support, so every triangle
product and every mean charge E eta is zero. Nevertheless, exact matching
expansion gives

    M22(QKQ^T)-M22(K) = -a^2 b^2/54,
    Ebad4(QKQ^T)-Ebad4(K) = a^2 b^2/27.

At a=1/4,b=1/2, the common V is155/64 and e4(K)=1/64, whereas
M22 is901/1024 versus24319/27648 and Ebad4 is5/4 versus2161/1728.
Thus even full spectra, V and the entire zero mean-charge cochain do not
recover this second moment. These exact identities are independently
verified symbolically and by direct rational binary rounding; they do not
claim that no spectral sufficient bound can exist.
