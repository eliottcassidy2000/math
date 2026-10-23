# Session 35: the two recursion modes meet in a marked resolvent

Status: proved elementary identities below; finite-exact verification recorded separately; proposed zeta application remains OPEN. No repository mutations.

## Inheritance and corrections

- `05-knowledge/hypotheses/HYP-1881-lrc-column-row-recursion.md` and `HYP-1891-lrc-two-mode-adic-recursion.md` define n=2^r b, b odd, with b→b+2 on the top row and n→2n vertically. Their LRC strategy is OPEN, not a theorem.
- `01-canon/theorems/THM-291-mode-b-multilinear-recursion.md`: actual Mode B forgets two endpoint vertices / two staircase legs and apex. A canonical section E(T)=source → T → sink satisfies H(E(T))=H(T). October 4 correction: n→n+2 adds 2n−1 free boundary bits, not 2n−3.
- `01-canon/theorems/THM-793-mode-b-line-tower.md`: descent preserves framed involutions, with free boundary fibers. October 20 correction: self-converse count is A002785, not A000570.
- `01-canon/theorems/THM-447-skew-sylvester-doubling.md` gives the concrete double D(M)=[[M,M+I],[M−I,−M]]. Generic blowup T[K2] is a different operation; do not conflate them.
- `07-reflections/riemann-zeta-tournament.md` has a correction retracting numerical Bernoulli/tournament prime-label unification. Its historical critical-strip analogy is not a transported zero theorem. `two-models-staircase-recursion.md` retracts SC(n)=A000568(n−1).

Concept board: odd-core grid; framed boundary E; skew-Sylvester D; marked resolvent; theta/Mellin; functional-equation fixed line.

## 1. The grid does not say literal addition and doubling commute

A(n)=n+2, d(n)=2n satisfy d A=A² d, not d A=A d. On address space (r,b), horizontal H(r,b)=(r,b+2) commutes with V(r,b)=(r+1,b), but downstairs H adds 2^(r+1), not 2. A theorem on commuting address moves is not a commuting theorem for tournament operations.

Even the size-repaired tournament identity D E=E² D is false. Starting from one vertex, D(E(T)) has sorted scores [1,2,2,2,3,5], whereas E²(D(T)) is transitive with scores [0,1,2,3,4,5]. Their characteristic polynomials are respectively
(z²+1)(z²+7)² and (z²+1)(z⁴+14z²+1).

## 2. A proved common carrier: the marked resolvent

Let M=A−Aᵀ be a tournament dominance matrix of size n, u the all-ones column,
p_M(z)=det(zI−M), F_M(z)=uᵀ(zI−M)^(−1)u.
F is odd because M is skew-symmetric. E is represented by

E(M) = [[0,uᵀ,1], [−u,M,u], [−1,−uᵀ,0]].

The pair (p,F) closes under BOTH modes:

E:
 p_E(z)=p_M(z)[z²+1+2zF_M(z)],
 F_E(z)=[(z²+1)F_M(z)+2z]/[z²+1+2zF_M(z)].

D, with w²=(z²+1)/2:
 p_D(z)=2^n p_M(w)²,
 F_D(z)=(z/w)F_M(w).

The expressions are independent of square-root branch, by parity, and interpreted as meromorphic identities; cancellation removes representation singularities, not genuine resolvent poles. They are polynomial/rational identities, not restrictions away from exceptional roots.

Proof of E: eliminate the interior block B=zI−M. Its determinant Schur complement is z²+1+2zF. To solve the resolvent against the all-ones vector, let x,y denote the two endpoint coordinates, a=x+y and b=x−y. Then za+b=2 and a=(z+2F)b−2F. Substitution gives the displayed F_E.

Proof of D: D²=diag(2M²−I,2M²−I). Hence its characteristic polynomial is det((z²+1)I−2M²)=2^n p_M(w)²; monicity fixes the sign. Also (zI−D)^(−1)=(zI+D)diag(R,R), where R=((z²+1)I−2M²)^(−1). Sandwiching with (u,u) leaves 2z uᵀRu, since uᵀMRu=0. This equals (z/w)F_M(w).

Interpretation: additive growth is a Möbius boundary transfer; doubling is a quadratic spectral substitution. This is an exact noncommutative renormalization calculus, not merely a numerical parallel.

### Why the mark is necessary

In upper-triangle edge order (12,13,14,23,24,34), matrices with signs
(-1,-1,-1,-1,-1,-1) and (-1,-1,+1,-1,-1,-1)
have the same p=z⁴+6z²+1 but different marked numerator uᵀadj(zI−M)u:
4z(z²+1) versus 4z(z²+5).
After E their characteristic polynomials differ:
z⁶+15z⁴+15z²+1 versus z⁶+15z⁴+47z²+1.
Thus unmarked spectrum, and any scalar spectral zeta derived only from it, cannot propagate the +2 recursion. The lost coordinate is overlap with the uniform boundary vector.

The pair (p,F) is a sufficient common carrier, not claimed minimal or complete for tournament isomorphism or Hamiltonian counts. E is only the zero-boundary section of full Mode B; general boundary assignments require a matrix-valued boundary resolvent.

## 3. A clean zeta intertwining, and its limit

Q_M=I−M² is positive definite and Q_D=2 I₂⊗Q_M. Consequently
 Z_M(s)=Tr(Q_M^(−s)), K_M(t)=Tr(exp(−tQ_M))
satisfy
 Z_D(s)=2^(1−s)Z_M(s), K_D(t)=2K_M(2t).
These identities hold for complex s and positive t, respectively. They follow directly from matrix functional calculus. After normalizing by size, the factors become 2^(−s) and a pure rescaling.

Every finite Z_M is entire. Multiplication by 2^(1−s) creates NO new zeros. Iterating a single seed gives only a scaled copy of its old zero set. The direct sum of all levels has, initially Re(s)>1,
 Z_tower(s)=Z_M(s)/(1−2^(1−s)).
This has an artificial vertical lattice of possible poles rather than Riemann's zeta structure. In particular the one-vertex seed has Z_M=1 and produces 1/(1−2^(1−s)), plainly not ζ(s).

There is also a scale mismatch: skew-Sylvester Q-energy doubles, so frequencies sqrt(Q) scale by sqrt(2), with multiplicity doubling. Integer theta under n→2n has energy scaling by 4. Direct identification without correcting this is false.

## 4. Where both arithmetic modes genuinely reconstruct zeta

For Re(s)>1, unique decomposition n=2^r(2k+1) gives
 ζ(s)=sum_{r≥0} sum_{k≥0} [2^r(2k+1)]^(−s)
      =(1−2^(−s))^(−1) sum_{k≥0}(2k+1)^(−s).

This is literally the user's two-axis grid: +2 enumerates odd seeds, ×2 supplies their columns. It is an exact arithmetic partition, not yet a tournament invariant.

Let θ(t)=sum_{n∈Z}exp(−πtn²). Then
 θ(t)−1 = 2 sum_{r≥0,k≥0}exp(−πt 4^r(2k+1)²).
The arithmetic doubling acts by t→4t. Mellin transformation gives
 π^(−s/2)Γ(s/2)ζ(s)=1/2 ∫_0^∞(θ(t)−1)t^(s/2−1)dt
in Re(s)>1. Poisson theta duality θ(t)=t^(−1/2)θ(1/t), plus continuation, yields s↔1−s.

This clarifies three DIFFERENT operations:
1. additive enumeration of odd seeds;
2. multiplicative dyadic scaling of those seeds;
3. inversion of heat time t↔1/t (duality).
Doubling is not itself functional-equation reflection.

The trivial zeros ζ(−2m)=0, for integers m≥1, cancel poles of Γ(s/2); after completion Λ(s)=π^(−s/2)Γ(s/2)ζ(s), those points are not zeros. In fact Λ(−2m)=Λ(1+2m)>0. The critical line is the fixed set of s↦1−conjugate(s); RH says all nontrivial zeros lie there and remains unproved. Thus −2m and 1/2 are not two halves of one doubling orbit. They are local gamma cancellation versus global self-dual spectral behavior.

## 5. Bold unification worth pursuing

Build a marked tournament boundary/scale theory: (p,F) as the first exact state, E as Möbius boundary transfer, D as quadratic scale transfer, then seek a matrix-valued boundary completion that also supports a Poisson-type inversion. The two proved recursions make the proposal testable. A functional equation cannot simply be imposed by symmetrizing a polynomial: it must come from a genuine transform and a positivity or zero-location mechanism.

Next decisive tests:
- Extend to arbitrary Mode-B boundary vectors with the full small matrix of resolvent entries; classify precisely which states close.
- Derive invariants of alternating words in E,D and test whether normalized marked spectral measures converge. The noncommuting witness warns against treating words as grid coordinates.
- Ask whether any limit has a reciprocal theta identity; test small/large t asymptotics before guessing an RH analogue. Fixed finite spectra cannot supply the lattice heat singularity.
- Keep spectral zeros distinct from Hamiltonian path invariants. The current carrier preserves both recursion formulas but not all original tournament information.

## Sources / exact checks

Primary reference checks: NIST DLMF https://dlmf.nist.gov/25.4 (functional equation), https://dlmf.nist.gov/25.10 (trivial zeros and RH statement), https://dlmf.nist.gov/20.7 (theta modular transformations). Formula derivations above are supplied directly rather than attributed novelty claims.

Run: python3 session35-recursions/check.py
- All 75 labeled tournaments of sizes 1..4: Q doubling and E characteristic identity, exact SymPy.
- All 8 labeled size-three tournaments: E and D marked-resolvent identities, exact SymPy.
- Explicit noncommutation and same-spectrum/different-boundary hostiles in check.out.

## 6. Second positive round: linearizing the additive recursion

A Cayley transform of the marked response makes +2 literally multiplicative. Put C_M(z)=(1+F_M(z))/(1−F_M(z)). The E formula gives exactly

 C_E(z) = ((z+1)/(z−1))² C_M(z).

Thus repeated endpoint extension has an explicit all-depth law C_{E^k(M)}=((z+1)/(z−1))^(2k) C_M. This is an exact scattering-coordinate linearization, not an identity between the underlying E and D operations. For a transitive n-tournament,

 F(z)=[(z+1)^n−(z−1)^n]/[(z+1)^n+(z−1)^n],
 p(z)=[(z+1)^n+(z−1)^n]/2.

The Cayley coordinate therefore sees the additive size as an exponent. The doubling remains F_D(z)=z F(w)/w rather than simply C_D=C_M². This exhibits precisely where a desired additive/multiplicative conjugacy fails and what prefactor must be handled.

Verification hardening: check.py uses explicit RuntimeError checks (no optimization-removable assertions). Ordinary and python3 -O runs give byte-identical check.out.
