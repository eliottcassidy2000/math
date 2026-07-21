## death-star-2026-07-21-S90 -- The NC2 tied-core weights ARE the CENTRAL TRINOMIAL (A002426) = a free-probability moment; completes the tournament↔NC2 free-prob bridge (S88→S89→S90). NC2 wall OPEN. HYP-8790.

**Owner directive:** keep finding tournament↔NC2 connections; push/pull often.

- **IDENTIFICATION (verified).** At the NC2 resonance central offset (= fully-confluent Vandermonde = fully-regular/Paley tournament, S89), the channel weights sum to W(m)=Σ_i m!/(i!²(m-2i)!) = 1,3,7,19,51,141,393,1107,3139,8953,... = **A002426 central trinomial** = [x⁰](1+x+1/x)^m = the m-th MOMENT of a 3-atom free convolution (Wigner/free-prob), ratio→3, ~3^m/√(πm).
- **CLOSES the free-prob bridge** boxeph-S203 flagged: the free-moment sibling of THM-438 (Paley cluster integrals=Catalan=free CUMULANTS of ½(δ_a+δ_{−a})). The wall=regular/Paley tournament carries a free-probability law on BOTH its H-count (THM-438) and its NC2 channel weights (this).
- **SHARP:** NC2 fails on the wall ⟺ central-trinomial-weighted signed channel sum =0 ∀m ⟺ free-cumulant series has a real positive zero ⟺ Laguerre-Pólya failure (boxeph-S202) ⟺ Paley spectrum leaves Re=−1/2 (char_S=∏(x²+p), THM-1555/213) — 3 faces (combinatorial/analytic/spectral) of ONE tournament fact.
- **ARC S88→S90:** channels-form-a-tournament (S88/boxeph THM-2033) → wall=regular/Paley (S89) → wall-weights=central-trinomial=free-prob (S90). Synthesis not proof; NC2 wall OPEN. Credits boxeph THM-2033/S202, codex, THM-438, klein THM-1805. reflection the-nc2-tied-core-weights-are-the-central-trinomial-...-S90 (+out).

## boxeph-2026-07-21-S203 -- THE VANDERMONDE IS THE BRIDGE: tournaments <-> NC2 (THM-2033)
## death-star-2026-07-21-S89 -- The NC2 WALL IS the regular/Paley tournament: completing boxeph's bridge + unifying EVERY repo wall (NC2 = H≥disc = LRC) into ONE object. NC2 wall OPEN. HYP-8785.

**Owner directive:** keep finding tournament↔NC2 connections; push/pull often.

- **BUILDS ON boxeph THM-2033** (NC2 channel-det = ∏a_i!·Vandermonde(radial degrees) = signed tournament sum, klein THM-1805; distinct=transitive=noncancel, repeated=confluent wall).
- **MY STEP (verified):** channel radial degree D(i)=i+i·degA+(m-2i)·degB. DEGREE-GAP → D(i)=[0,3,6,9,...] DISTINCT (transitive Vandermonde≠0). RESONANCE CENTRAL OFFSET (degA=degB=1) → **D(i)=m for EVERY i** (fully-confluent). By klein THM-1805 (transitive⟺distinct scores): repeated degrees=repeated SCORES=REGULAR tournament; ALL equal = DOUBLY-REGULAR = PALEY/DRT. So **NC2 resonance wall = fully-confluent Vandermonde = the regular/Paley tournament.**
- **UNIFICATION (the payoff):** NC2 wall = H≥disc wall (S84 regular/Paley tightest) = LRC wall (THM-640 AP=Paley) = ONE object, the regular/Paley (equal-score, big-stabilizer, S76) tournament. Transitive=easy pole (distinct scores), regular/Paley=hard pole (equal scores) = the two S75 poles, shared across 3 flagship problems.
- **ANALYTIC FACE:** the fully-confluent (regular) channel sum's asymptotic = Wigner/free-cumulant (THM-438, H(Paley)~e·avg); codex hyper-Bessel + boxeph Laguerre-Pólya boundary = the Paley char_S=∏(x²+p) real-rooted spectrum = Re=−1/2 critical line = quasirandomness (S85). NC2 noncancel on the wall = confluent Paley Vandermonde/Wronskian≠0 = real-rootedness = reality of the Paley spectrum (tournament-spectral, THM-1555/213).
- Unification not proof; NC2 wall OPEN. Credits boxeph THM-2033/S202, codex, klein THM-1805. reflection the-nc2-wall-is-the-regular-paley-tournament-...-S89; script nc2_confluent_vandermonde_is_regular_S89 (+out).

## death-star-2026-07-21-S88 -- The CHANNEL-TOURNAMENT LENS: NC2 is a tournament-nullcone on its radial channels; the regular-channel core is the wall; explains why domination (MISTAKE-202) was refuted. NC2 OPEN. HYP-8772.

**Owner:** find connections between tournaments and NC2; long session, push/pull often.

**THE BRIDGE (verified).** NC2 noncancellation = det[(a_i+k)!] (THM-1815) = ∏a_i!·Vandermonde(radial degrees) = Σ_T sgn(T) x^{score(T)} = the SIGNED TOURNAMENT SUM (THM-1805/my THM-1925; transitive tournaments survive). Verified exactly: det=∏a!·Vand over 7 degree sets; Vandermonde=tournament-sign-sum (n=3:30,n=4:1440). So the object deciding NC2 IS the tournament sign-sum over the channel degrees.

**TWO REGIMES.** DISTINCT degrees => Vandermonde≠0 => TRANSITIVE channel (death-star HYP-8772) => noncancellation (codex THM-2017 degree-gap). REPEATED degrees = the resonance WALL => ordinary Vandermonde VANISHES but the CONFLUENT Vandermonde (derivative/Wronskian row) SURVIVES nonzero (verified) = codex's 1/m hyper-Bessel correction / my Laguerre-Pólya boundary ODE θ²Φ=ξΦ (HYP-8775). Respects MISTAKE-212 (tied channels = CONFLUENCE not intransitivity; degree preorder ties, derivative order breaks it).

**codex's λ IS the node-spacing.** Channel factorial-degree D(k)=dm+λk (λ=e-rd), so the transitivity Vandermonde = λ^{C(nch,2)}·(int): |λ| = node spacing, r-|λ| = confluence order. |λ|>=r+1 separated (THM-2017), |λ|=r boundary (my L-P), 0<|λ|<r band (HYP-8766), λ=0 FULLY CONFLUENT = central resonance (codex HYP-8771) = my τ=1 regular core. codex's regime map = the confluence-order stratification of ONE discriminant.

**META (fermionic/bosonic).** The tournament sign-sum = Vandermonde = a DETERMINANT (fermionic, signed) -- nonvanishing is free (distinct nodes; the sign-involution collapses to the transitive core). NC2's E[P^m] = a PERMANENT (bosonic, all +, no sign) -- doesn't collapse, hence hard. THM-1815's miracle: the bosonic noncancellation reduces to a fermionic Vandermonde on the channel degrees -- NC2 borrows the fermion's rigidity. The wall (regular/tied core) = my invariant-resistant hot center (THM-2016): maximally symmetric, confluent, no cheap invariant separates it -- WHY domination failed (MISTAKE-202: regular channels have no source).

**Unifies FOUR threads as ONE object (confluence of the tournament sign-sum):** THM-1815 transitivity Vandermonde / THM-1805-1925 tournament sign-sum / death-star HYP-8772 channel lens / codex THM-2017 hyper-Bessel + my HYP-8775 L-P boundary. Plus my continuum (THM-1979/2013/2016) reads the same axis as cyclic temperature.

**Push/pull:** checkpointed twice; codex reserved THM-2023 (proving the L-P claim I raised) -- convergent. NEW connection for the fleet.

**Honest scope:** SYNTHESIS/identification (verified computations + THM-1805/1815), NOT a proof of NC2. The residual is ONE discriminant at distinct nodes (opus THM-1710 multinomial-ratio, shared with TNC) + its confluent limit at the tied core (hyper-Bessel/L-P). Backlog lead: the free-probability/semicircle bridge (death-star S88: channel weights=free cumulants, tied core=semicircle=regular tournament=my τ=1; Catalan/Wigner = THM-438 H(Paley)~e·avg).

**Next:** develop the free-prob/semicircle bridge (regular tournament spectrum <-> NC2 tied-core entropy saddle); connect the confluent Vandermonde to the Hermite/Wronskian radial closure (THM-1615). Artifacts: THM-2033, HYP-8780, reflection the-vandermonde-is-the-bridge-tournaments-and-nc2-boxeph-S203.md, script confluent_transitivity_vandermonde_boxeph_S203.py (+.out).

