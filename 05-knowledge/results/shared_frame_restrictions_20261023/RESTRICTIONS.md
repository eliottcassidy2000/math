# Which observer restrictions preserve the 8/9 coupling gap?

23 October 2026. **PROVED, scoped:** subgroup diagonalization; the sharp
contiguous-block relative gap; fixed-port slow modes; the fixed-Hamiltonian-path
absorbing state. **FINITE-EXACT:** the declared universes in `receipt.json`.
No claim identifies this chain with the native tournament metagraph.

## 1. Inherited correspondence and the operator being restricted

For a graph G and a linear observer order pi, orient every pair forward in pi
except edges of G, which are oriented backward. This is a reversible map
between **ordered** graphs and **ordered** tournaments. Complementing G or
reversing pi separately gives the converse; doing both leaves the tournament
unchanged. A Hamiltonian path supplies an observer, but generally not a unique
one. The unmarked graph/complement pair does not determine a tournament class.

The precise prior results are
[Session42's correspondence](../session42_shared_lattices_square_colors_20261023/correspondence/CORRESPONDENCE_AND_OBSTRUCTIONS.md),
[the path/tree trijections](../tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md),
and [the full-cube coupling theorem](../session42_shared_lattices_square_colors_20261023/coupling/COUPLING_OPERATOR.md).
The latter uses the signed orbit basis of
[THM-4721, twisted Burnside and signed Walsh basis](../../../01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md).
We extend its subgroup scope, not the already known full-cube result.

Fix the reference order 1,...,n and the uniform cube X=F2^binom(n,2).
Let H be ANY subgroup of S_n. Define P_G^H by averaging ordinary edge
relabelings and P_T^H by averaging tournament relabelings, represented as
affine edge-bit maps. Both are orthogonal projections. Set

    K_H=P_G^H P_T^H P_G^H

on the ordinary H-invariant space. A restriction of the **group** is different
from restricting the **set of allowed words**. Sections 2–4 concern the former;
section 5 concerns the latter.

## 2. Exact subgroup spectrum

For an edge set F, write chi_F(x)=(-1)^|F intersect x|. For its ordinary
H-orbit O, let u_O=|O|^-1 sum_(F in O) chi_F, and define

    mu_H(F)=|H|^-1 sum_(g in H) (-1)^inv_F(g),
    inv_F(g)=#{i<j: ij in F and g(i)>g(j)}.

Then

    K_H u_O = mu_H(F)^2 u_O.                              (1)

**Proof.** On the character span of O, the ordinary invariant space is one
line. The signed permutation action has either no invariant line (a
stabilizer reverses sign) or exactly one, with coefficients all +/-1 after
normalization. In the first case both mu and the coupling vanish. In the
second case the overlap of the normalized unsigned and signed lines is
mu_H(F), up to sign. Orthogonal projection there and back squares the
overlap. Relabeling the chosen F can change the sign but not its square.
Different mask orbits are orthogonal, so this is the full spectrum. QED.

Eigenvalue 1 can now have multiplicity greater than one. Its space is exactly
the intersection of the two invariant spaces: equality in the norm
contraction by P_T forces the vector to be fixed by P_T. We use
**relative gap** to mean 1 minus the largest eigenvalue on the orthogonal
complement of this common space. It is not a claim of global irreducibility.

Complement/converse pairing keeps the even-|F| orbits. Unlike the unrestricted
S_n case, odd-|F| modes need NOT have eigenvalue zero. Pairing can therefore
change nonzero restricted spectra. Formula (1), with only even masks retained,
is the correct paired rule.

## 3. Positive repair: keep cuts as well as ports

Partition the observer into contiguous intervals B_1,...,B_r, and take
H=product_i S_(B_i). A singleton interval is a fixed port; requiring every
other vertex to stay in its interval also preserves which side of that port
it occupies. For every F,

    mu_H(F)=product_i mu_(S_|B_i|)(F[B_i]).                (2)

Indeed cross-interval endpoints never reverse their order. The permutations
inside different intervals are independent, and the inversion sign factors.

The inherited full-symmetric-group bound is |mu(J)|<=1/3 for every nonempty
graph J. It follows from vertex deletion, with the remaining star case
given by the alternating average of its center rank. Consequently:

* Eigenvalue 1 occurs precisely on mask orbits with no within-interval edge.
* Every other eigenvalue is at most 1/9.
* If any interval has size at least 3, a two-edge star inside it attains 1/9.
  Thus the **relative gap is exactly 8/9**, at every size.
* If every interval has size at most 2, every residual eigenvalue is zero.
  If all intervals are singletons, there is no residual space at all.

The two-edge star is even, so the same sharp 8/9 relative gap holds after
complement/converse pairing whenever an interval has size at least 3.
For smaller intervals one must allow the paired residual space to be empty.

This is a genuine all-size restriction theorem. It does not preserve a
specified Hamiltonian path's arcs, or arbitrary electrical/path constraints:
those are additional word restrictions, not consequences of preserving cuts.

## 4. Negative boundary: one fixed interior port can create arbitrarily slow modes

Now fix vertex r but allow the other m=n-1 vertices to permute across it.
There are l=r-1 available slots on its left and u=n-r on its right. For any
k-edge star centered at r,

    |mu_H(F)| = | sum_j (-1)^j C(l,j) C(u,k-j) | / C(m,k). (3)

**Proof.** The image set of the k leaves is a uniformly chosen k-subset of
the free slots. A fixed reference sign times (-1)^j is the inversion sign
when j leaves land to the left. Counting such subsets gives (3). This is
the normalized Krawtchouk coefficient for this actual port experiment.

For r=2 and k=1, the eigenvalue is ((n-3)/(n-1))^2, strictly below 1 for
n>=3. Therefore the unpaired relative gap is at most

    4(n-2)/(n-1)^2 -> 0.                                 (4)

Pairing does not rescue a uniform bound. The two-leaf star is even and gives
((n-5)/(n-1))^2. For n>=4 this is a residual mode, yielding

    paired relative gap <= 8(n-3)/(n-1)^2 -> 0.           (5)

These are upper bounds from exhibited modes, not a claimed complete
closed formula for the fixed-port spectrum.

There is already a smaller hostile at n=4: fix port 2 and use
F={12,13}. The six permutations of {1,3,4} have sign average 2/3, hence
eigenvalue 4/9; F has two edges and survives pairing. Both relative gaps
are at most 5/9. The finite complete spectra confirm equality there.
This F is a two-edge star centered at an **unmarked** vertex, so it must
not be confused with the centered-port stars in (3).

The preserved port itself is not the problem. The lost coordinate is the
order-side occupancy of its neighbors. Formula (2) prevents that loss;
formula (3) measures it exactly when crossing is allowed.

## 5. Fixed-path word restrictions: a different failure mechanism

Let Omega_P be all edge words with x_(i,i+1)=0 for i=1,...,n-1. Thus every
tournament has the actual path P=(1,...,n). On Omega_P, form ordinary graph
classes and tournament classes by intersecting their full-cube classes
with Omega_P. With uniform word weights define incidence W and

    K_gh = (1/g_g) sum_t W_gt W_ht/t_t,                   (6)

where g_g and t_t are the actual restricted row/column masses.

**For every n>=3 this chain has an isolated transitive state and has no
positive global mixing gap.** The empty graph contributes one word. A
transitive tournament has a unique Hamiltonian order, so a transitive
tournament containing P must be the reference transitive tournament.
Thus its column also contains exactly that word. There are other words
when n>=3. This remains true with full graph-complement and tournament-
converse classes, since the complete backward graph is not in Omega_P.

At n=3 the two words give two singleton components; at n=4,5,6 the exact
checks find the singleton and one remaining component. The latter pattern
is not asserted at arbitrary order. It suggests a new, meaningful target:
the restricted relative gap inside the nontransitive components, with the
source weights and any native port filters retained.

If the ENTIRE indexed path is also a rigid mark in tournament comparisons,
the only relabeling preserving it fixes every vertex, so P_T is the identity.
Then P_G P_T P_G is the identity on the graph-invariant space. Exact path
retention is useful for reconstruction; it does not create a resampling
chain by itself. Separate native moves are required.

Finally, in this path-normalized chart, converse followed by path reversal
sends an edge mask to its coordinate reflection. It is not the operation
of complementing the free chords while retaining the path arcs. Both
operations are legitimate, but they implement different maps.

### Removing the isolated state does not restore 8/9

At n=4, delete the transitive word. In graph-row order with representative
masks 2,6,18,22 (bits listed lexicographically by unordered pair), the
remaining transition matrix is, both paired and unpaired,

    [11/15  2/15  1/15  1/15]
    [ 1/5   2/5   1/5   1/5 ]
    [ 1/5   2/5   1/5   1/5 ]
    [ 1/5   2/5   1/5   1/5 ].

Its stationary masses are (3,2,1,1)/7. Its spectrum is
1,8/15,0,0, so this CONNECTED component has gap **7/15**. This is a
stronger hostile than merely pointing to the isolated transitive state.

There is also an all-order exact source-row formula. Put
M=binom(n-1,2), n>=4. In the nontransitive fixed-path carrier, the one-edge
graph row has mass M and stationary probability pi=M/(2^M-1). Its return
probability is

    d_n = [n-2 + sum_(ell=3)^(n-1) (n-ell)/(2^(ell-1)+1)]/M.       (7)

The centered indicator of this row therefore has Rayleigh quotient

    (d_n-pi)/(1-pi).                                      (8)

Here is the counting mechanism. A sole backward chord i<j of length
ell=j-i>=2 has one strong component of size k=ell+1, with singleton
source/sink padding. Its core is the transitive k-tournament with only
the extreme arc reversed, namely C3(1,1,k-2). The inherited
[Session26 two-apex theorem](../session26_random_bridges_20261022/session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md)
gives H=1+2^(k-2), rigidity for k>=4, and the separate C3 automorphism
group of size 3. A direct path proof is also short: the increasing path
is one choice; a path using the one backward arc partitions the middle
vertices freely between its two increasing sections, giving 2^(k-2).

Thus the fixed-path tournament column has size H/|Aut|: one when ell=2,
and 2^(ell-1)+1 otherwise. Different chords give different unpaired
tournament types, specified by their strong core and ordered padding.
There are n-ell chords of length ell, proving (7). Converse merges only
mirror padding pairs with equal masses and preserves the return formula.
The graph-complement partner of a one-edge word is outside this path
carrier for n>=4, so it introduces no extra row mass.

At n=5, (8)=463/855, proving a gap at most 392/855 for its connected
nontransitive component, again paired or unpaired. For general n, (8)
is a Rayleigh quotient on the complement of constants; if further
components occur, it must not be called a relative-gap certificate after
those additional invariant directions have been removed.

The two-apex count must not be confused with the full-tiling partner's
tribonacci count. The old guess 2^(n-2)+1 was false for that DIFFERENT
family, already corrected by THM-4731. Our subset-DP control explicitly
builds the one-long-backward-chord core for k=3,...,12, including the
order-seven separation (33 here, 31 in the tribonacci family).

## 6. Verification and next questions

`check.py` verifies 1,065,509 mask factorization equations for all 63 ordered
compositions with n<=6. It independently builds literal word-orbit incidence
matrices and checks 633 exact rational eigenvector equations for all
contiguous Young groups, all one-port stabilizers, and a cyclic subgroup
at n=2,3,4, both paired and unpaired. It checks every centered-port star
at n=3,...,6 and every path word at n=3,4,5,6. Normal/-O agreement is required.
It also verifies the component matrix, Rayleigh quotients and all ten
one-chord core path counts by actual-arc subset DP. No floating-point
spectral calculation is used.

**OPEN / next falsifiable targets.** (i) Determine the full fixed-port
spectrum rather than only the star lower modes. (ii) Bound relative gaps
on nontransitive fixed-path components. (iii) Compare actual native moves
with a cut-preserving block chain by a congestion bound that retains all
guards. A mere gap of the auxiliary chain supplies none of these comparisons.
