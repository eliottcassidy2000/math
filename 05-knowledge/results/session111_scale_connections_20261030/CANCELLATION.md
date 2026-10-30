# Clone dilation and a four-phase cancellation recurrence

**PROVED:** all-parameter polynomial laws and tournament construction from the specified seed. **FINITE-EXACT:** seed automorphism certification and coefficient checks. No unrestricted A000568 enumeration is claimed.

## Source and realization

The inherited anchor is [THM-4764](https://github.com/eliottcassidy2000/math/blob/5ab37f308003e8634d13ea4c0dfe8f1558c11054/01-canon/theorems/THM-4764-modular-cancellation-depth-of-odd-regular-gap-actions-and-the-fifty-vertex-realization.md).
Related older result folders are clone_depth_modular_twist_average_20261025, unprotected_ports_cancellation_20261025, session63_integral_observers_root_arrangements_20261025, ports_regular_actions_hn_20261025 and session62_agent_theme_bridges_20261025.

Use R=Cay(F5²,S), with
S={(0,1),(0,2),(1,0),(1,1),(1,2),(1,3),(1,4),(2,0),(2,1),(2,2),(2,3),(3,1)}.
Our independent checker exhausts the stabilizer of0:106 search nodes and one automorphism. Translation transitivity proves Aut(R)=C5². It also checks that no pair of vertices is homogeneous.

For r≥1 set T_r=R[TT_(r+1)]. Inside each rigid transitive child, exactly consecutive vertices form homogeneous pairs. Across children, an outside witness from R rules out homogeneity. The homogeneous-pair graph therefore has25 path components, and every automorphism permutes these rigid directed paths. It follows that Aut(T_r)=Aut(R), of order25.

The canonical gaps are F5²×{1,...,r}, with diagonal translations. There are r gap orbits, so r independent symmetry-respecting toggles. Identity has25r singleton cycles; every nonidentity translation has5r cycles of length5. Burnside gives the subset-orbit polynomial
P_r(t)=[(1+t)^(25r)+24(1+t^5)^(5r)]/25.
The same group acts on all copies. Replacing this by P_1(t)^r would discard the diagonal coupling.

## General odd-prime theorem

Let an exponent-p group of order M=p^a, p odd and a≥1, act regularly on each of r diagonal copies. The group need not be abelian. Its subset polynomial is
P(t)=[(1+t)^(Mr)+(M−1)(1+t^p)^(Mr/p)]/M.

Define cancellation depth as multiplicity of t=−1 after taking this integer polynomial into the indicated field. Then:
- over Q: depth Mr/p;
- modulo p: depth Mr−p+1 if p does not divide r, and Mr if p divides r;
- modulo every prime ell dividing M−1: depth Mr.

Proof of the modular case must divide in the integers first. Set s=1+t and
h_p(s)=[1+(s−1)^p−s^p]/p.
Its coefficients are integers and its lowest term is s. Expand (s^p+p h_p)^(p^(a−1)r). The zero-order term combined with the first summand gives s^(Mr); the linear term after division by M is
(M−1)r s^(Mr−p)h_p.
For an expansion term of index j≥2,
v_p[p^(j−a) binom(p^(a−1)r,j)]
≥j−1+v_p(r)−v_p(j)≥1.
The first inequality uses j binom(n,j)=n binom(n−1,j−1). The second uses odd p. Thus
P(s−1) ≡ s^(Mr)+(M−1)r s^(Mr−p)h_p(s) (mod p),
which proves both branches. The rational least degree comes from the second summand. In characteristic ell|M−1 the second summand disappears and M remains invertible.

The hypotheses exclude p=2 and the trivial group. They concern this specified permutation action, not every tournament of the same size.

## Requested sizes and the hostile next phase

For the25-point seed:
| vertices | gaps | toggle rank | Q depth | F5 depth | F2/F3 depth |
|---:|---:|---:|---:|---:|---:|
|50|25|1|5|21|25|
|100|75|3|15|71|75|
|150|125|5|25|125|125|
|200|175|7|35|171|175|

The150-vertex row refutes an unqualified rule “F5 depth is gap count minus4.”

Now double the child size: n_j=25·2^j, j≥1, r_j=2^j−1. Since ord_5(2)=4,
d_j=depth_F5(T_(r_j))=25·2^j−29+4·[4 divides j].
The first values are21,71,171,375,771,1571,3171,6375. In particular the first full-depth exception in this doubling tower is400 vertices, not200.

The ordinary generating function is
Σ_(j≥1) d_j z^j = 50z/(1−2z) −29z/(1−z) +4z^4/(1−z^4).
Consequently
d_(j+5)=2d_(j+4)+d_(j+1)−2d_j, j≥1.
The characteristic polynomial is (x−2)(x^4−1). All five poles survive, so the rational linear recurrence has minimal order5.

The Q-depth generating function is10z/(1−2z)−5z/(1−z). The F2/F3-depth generating function is50z/(1−2z)−25z/(1−z). These are depth sequences for an explicit family, not a generating function for A000568.

## What is merged, and what is not

This recovers the older50-vertex phenomenon and extends it through literal clone dilation. The four-phase correction comes from multiplication by2 in F5*, not from a resemblance to a spectral constant. The number200 can be a vertex count here and a group order elsewhere without identifying the objects.

[cancellation.js](cancellation.js) compares two exact expansions in44 parameter cases, including r=15,31,63 in the5² family, checks the recurrence for20 formula-generated terms, and certifies the seed group and homogeneous-pair condition.
