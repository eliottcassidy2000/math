# A path-count budget for every prime-power symmetry tier

Suppose a tournament has an automorphism of order p^a, p an odd prime. Let m_e count its vertex orbits of size p^e for e≥1, and m=Σ m_e. In each orbit of size q=p^e, cyclic shifts of any internal Hamilton path are q distinct paths, one per possible start. Because the action restricts to a regular cyclic action on that orbit, its arcs are Cayley increments mod q. Every permitted increment coprime to q gives q Hamilton paths by walking q−1 times with that same increment. Exactly half of the φ(q) units are permitted, since u and −u are paired and exactly one is an arc. Hence an orbit has at least q φ(q)/2 internal Hamilton paths, valid also at q=3.

Use the majority orientation between every pair of vertex orbits (including fixed vertices). For orbits of different sizes, the fractions of arcs in each direction still sum to 1, and the product of the two sizes is odd, so majority is unambiguous. Transitivity of the automorphism on each orbit makes the forward fraction constant over the choice of a terminal vertex in the earlier orbit. A Hamilton path of this majority quotient therefore gives actual concatenations of internal paths: for an adjacent nontrivial orbit of size q, an allowed preceding endpoint leaves at least (q+1)/2 of q available starts; a fixed singleton between nontrivial orbits makes all its incident relations uniform. Each nontrivial orbit then contributes at least (q+1)φ(q)/4 possibilities even if preceded by another nontrivial orbit, and at least qφ(q)/2 if it starts a run. Distinct choices are distinct spanning Hamilton paths. Thus
```
H(T) ≥ ∏_{e≥1} [(p^e+1) φ(p^e)/4]^(m_e).
```
This yields a single additive budget on the occupied symmetry tiers:
```
Σ_{e≥1} m_e log[(p^e+1)φ(p^e)/4] ≤ log H(T).
```
In particular Σ e m_e = O_p(log H), and every orbit length is O_p(sqrt(H)) (the isolated orbit bound above also gives this directly). Unlike the size of an entire tournament, these costs cannot grow arbitrarily by appending neutral fixed vertices.

If H is a prime p≥5, the inherited free action of Aut(T) on ordered spanning paths implies its order divides p. If it were p, some orbit of size p would contribute p(p−1)/2 internal Hamilton paths. Every one extends to a spanning path by tournament path insertion; extensions of distinct initial paths are distinct because restriction recovers their original vertex order. This contradicts H=p. Thus every prime-count tournament with prime at least 5 is asymmetric. This is a negative design principle for Chan–Pak: a symmetry-based construction must break all nontrivial automorphisms before producing such targets. It does not show those targets occur.

The count bound and the mod-p orbit fold in ORBIT_FOLD.md see different aspects of growth: the former consumes integer count budget with each occupied orbit; the latter can cancel nonuniform incidences modulo p and has an exact Taylor factor for each orbit-size tier.
