## boxeph-2026-07-21-S203 -- THE VANDERMONDE IS THE BRIDGE: tournaments <-> NC2 (THM-2033)

**Owner:** find connections between tournaments and NC2; long session, push/pull often.

**THE BRIDGE (verified).** NC2 noncancellation = det[(a_i+k)!] (THM-1815) = ∏a_i!·Vandermonde(radial degrees) = Σ_T sgn(T) x^{score(T)} = the SIGNED TOURNAMENT SUM (THM-1805/my THM-1925; transitive tournaments survive). Verified exactly: det=∏a!·Vand over 7 degree sets; Vandermonde=tournament-sign-sum (n=3:30,n=4:1440). So the object deciding NC2 IS the tournament sign-sum over the channel degrees.

**TWO REGIMES.** DISTINCT degrees => Vandermonde≠0 => TRANSITIVE channel (death-star HYP-8772) => noncancellation (codex THM-2017 degree-gap). REPEATED degrees = the resonance WALL => ordinary Vandermonde VANISHES but the CONFLUENT Vandermonde (derivative/Wronskian row) SURVIVES nonzero (verified) = codex's 1/m hyper-Bessel correction / my Laguerre-Pólya boundary ODE θ²Φ=ξΦ (HYP-8775). Respects MISTAKE-212 (tied channels = CONFLUENCE not intransitivity; degree preorder ties, derivative order breaks it).

**codex's λ IS the node-spacing.** Channel factorial-degree D(k)=dm+λk (λ=e-rd), so the transitivity Vandermonde = λ^{C(nch,2)}·(int): |λ| = node spacing, r-|λ| = confluence order. |λ|>=r+1 separated (THM-2017), |λ|=r boundary (my L-P), 0<|λ|<r band (HYP-8766), λ=0 FULLY CONFLUENT = central resonance (codex HYP-8771) = my τ=1 regular core. codex's regime map = the confluence-order stratification of ONE discriminant.

**META (fermionic/bosonic).** The tournament sign-sum = Vandermonde = a DETERMINANT (fermionic, signed) -- nonvanishing is free (distinct nodes; the sign-involution collapses to the transitive core). NC2's E[P^m] = a PERMANENT (bosonic, all +, no sign) -- doesn't collapse, hence hard. THM-1815's miracle: the bosonic noncancellation reduces to a fermionic Vandermonde on the channel degrees -- NC2 borrows the fermion's rigidity. The wall (regular/tied core) = my invariant-resistant hot center (THM-2016): maximally symmetric, confluent, no cheap invariant separates it -- WHY domination failed (MISTAKE-202: regular channels have no source).

**Unifies FOUR threads as ONE object (confluence of the tournament sign-sum):** THM-1815 transitivity Vandermonde / THM-1805-1925 tournament sign-sum / death-star HYP-8772 channel lens / codex THM-2017 hyper-Bessel + my HYP-8775 L-P boundary. Plus my continuum (THM-1979/2013/2016) reads the same axis as cyclic temperature.

**Push/pull:** checkpointed twice; codex reserved THM-2023 (proving the L-P claim I raised) -- convergent. NEW connection for the fleet.

**Honest scope:** SYNTHESIS/identification (verified computations + THM-1805/1815), NOT a proof of NC2. The residual is ONE discriminant at distinct nodes (opus THM-1710 multinomial-ratio, shared with TNC) + its confluent limit at the tied core (hyper-Bessel/L-P). Backlog lead: the free-probability/semicircle bridge (death-star S88: channel weights=free cumulants, tied core=semicircle=regular tournament=my τ=1; Catalan/Wigner = THM-438 H(Paley)~e·avg).

**Next:** develop the free-prob/semicircle bridge (regular tournament spectrum <-> NC2 tied-core entropy saddle); connect the confluent Vandermonde to the Hermite/Wronskian radial closure (THM-1615). Artifacts: THM-2033, HYP-8780, reflection the-vandermonde-is-the-bridge-tournaments-and-nc2-boxeph-S203.md, script confluent_transitivity_vandermonde_boxeph_S203.py (+.out).

