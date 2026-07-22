## boxeph-2026-07-21-S210 -- antisymmetry is the hinge: tori, odd functions, saddles, tournaments (HYP-8835)

**Owner:** see how tori relate to odd-valued functions, saddle points, and tournaments.

**THE HINGE:** a tournament's payoff M=A−Aᵀ is ANTISYMMETRIC (odd). All four links follow (verified, tori_odd_saddles_tournaments_boxeph_S210.py):

1. ODD => EVEN RANK => ODD SUPPORT. rank(A−Aᵀ) even for ALL tournaments n≤5; the tournament GAME (symmetric zero-sum, value 0) has optimal strategies on ODD support (Fisher–Ryan; census: sizes 1,3,5). Pure game SADDLE POINT <=> transitive (Condorcet winner support 1); 3-cycle = rock-paper-scissors support 3; regular odd-n uniform support n.
2. TORUS χ=0 NEEDS SADDLES. Poincaré(Tⁿ)=(1+t)ⁿ, Betti=C(n,k), χ=0; standing bagel T² = 1 max/2 saddle/1 min, the 2 saddles = b₁ = handle = the S207 deficit-1 (reduced Euler / bagel−cake=T_n−1).
3. TRANSITIVE=GRADIENT vs CYCLIC=TORUS. Replicator ẋᵢ=xᵢ(Mx)ᵢ: transitive flows to Condorcet sink (gradient, no torus); 3-cycle conserves H=x₀x₁x₂ EXACTLY (dH/dt=0 since M col sums=0; RK4 drift 1e-16, orbit closes) => invariant TORUS (recurrent center at (⅓,⅓,⅓)). Intransitivity IS the toroidal recurrent set (Conley reading; 3-cycle atom THM-1830 = minimal recurrent block).
4. ODD FUNCTIONS on the torus. Involution θ↦−θ has 2ⁿ fixed points (2-torsion); odd f vanishes there. LRC far-set weight/measure EVEN (sinc=odd/odd) => |G| t↦−t invariant; ODD sector = signed-discrepancy sawtooth B₁ (c_k=1/(2πik), odd in k) and the transitivity Vandermonde = SIGN character (odd under Sₙ, verified V=540↦−540).

**SYNTHESIS:** tori (χ=0 needs saddles), odd functions (θ↦−θ, sawtooth, sign character), saddle points (game equilibria: pure⟺transitive), and tournaments (antisymmetric beats-relation) are four faces of ONE antisymmetry. Transitive/cyclic = gradient/recurrent = saddle-present/absent = even/odd.

**LEVERAGE:** 3-cycle atom = elementary invariant torus (Conley/Morse decomposition of tournament space available); odd support = even-rank/Pfaffian (ties THM-473 skew=Hermite ±iλ toroidal spectrum); θ↦−θ governs both LRC reality-symmetry (THM-1820 mirror pairs) and tournament antisymmetry.

**Honest:** all pillars verified identities/census/dynamics; a correction en route (the sinc far-set weight is EVEN not odd; the odd object is the sawtooth B₁). Synthesis, not a new theorem. Artifacts: reflection antisymmetry-is-the-hinge-...-boxeph-S210.md, HYP-8835, script (+.out).

