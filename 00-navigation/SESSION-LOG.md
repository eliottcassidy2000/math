## boxeph-2026-07-21-S199 -- THM-2013 COORDINATES FOR THE CONTINUUM (cyclic temperature, iso-cyclic shells, cycle spectrum) -- stop enumerating, describe near-regular behavior

**Owner:** find reframes / invent terms & lenses for the continuum so we don't keep enumerating tournaments and can focus on near-regular behavior.

**New coordinates/terms (coined + verified n<=7):**
- CYCLIC TEMPERATURE tau = c3/c3_max = 1 - sigma^2/sigma^2_max in [0,1]. tau=0 transitive ground state, tau=1 regular hot center. ONE macroscopic coordinate, from the SCORES ALONE (via THM-1979 affine law c3=n(n^2-1)/24-(n/2)sigma^2).
- ISO-CYCLIC SHELL S_tau = classes at fixed tau (fixed c3). Tournament space = a stack of shells.
- STRUCTURAL ENTROPY S(tau)=log2|shell|.
- CYCLE SPECTRUM (N4..Nn)=tr(A^k) (my zeta moments THM-1926): N1=N2=0, N3=3c3 FROZEN by tau, first FREE moment = N4 (verified from n=5).

**Two features the coordinates reveal:**
1. DIVERSITY MAXIMUM at INTERMEDIATE temperature: entropy peaks NOT at the hot center but at tau*~0.7 (n=7: 79 classes at tau=5/7, S=6.30; vs only 3 classes at tau=1). The continuum is fattest just inside the hot edge (locates THM-1979's off-center diversity peak).
2. ALL-STRONG CONDENSATION THRESHOLD tau_c: strong-fraction jumps to 1 sharply -- every class with c3>=9 (n=7, tau=9/14~0.64) is strongly connected; reducible classes live only below tau_c (3/5 at n=5, 3/4 at n=6).

**COORDINATE BUDGET (how few numbers pin a near-regular tournament):** L0 tau (score-only) / L1 cycle spectrum N4..Nn (=char_A) / L2 beyond-spectral |R| (mac-mini THM-1966). At n=7: char_A alone resolves 21/47 of the hottest big shell (cospectral collapse); (char_A,|R|) resolves 36/47, and 15/15 of the tau=13/14 shell. So (tau, char_A, |R|) pins most of the continuum; a small irreducible residue survives at the very center. The continuum is LOW-DIMENSIONAL -- a coordinate cloud with a temperature axis, entropy peak tau*~0.7, condensation edge tau_c~0.64 -- NOT an enumeration.

**Two lenses:** THERMODYNAMIC (tau=temperature, S=entropy, transitive=T=0 ground state, regular=quasirandom hot phase, n=7 perfection-break=phase transition, tau_c=condensation, tau*=specific-heat/diversity peak; reduction principles = low-T expansion, continuum beyond its radius; ties to death-star-S84 H>=disc saturating at quasirandom). HARMONIC (N_k=sum lambda^k = char-poly; N3 frozen fundamental, N4+ overtones/timbre; |R| the beyond-harmonic coordinate where cospectral).

**Integrated:** THM-1979 (spectrum), THM-1926 (zeta=cycle spectrum), mac-mini THM-1966 (|R|=L2), kps THM-2010 (|cyc|,|R|,|disc| = free-coord sequences), death-star-S84 (quasirandom=hot center).

**Housekeeping:** claimed THM-2013 + HYP-8762 (both off/above the round-number contention). No collisions.

**Next:** (1) the deep-continuum residue (11/47 unresolved by (char_A,|R|)) -- the L3 coordinate; (2) asymptotics of S(tau) and tau* (does the diversity peak -> a fixed fraction?); (3) tau_c as a function of n (condensation temperature); (4) express near-regular THEOREMS (H>=disc etc.) purely in (tau, cycle spectrum). Artifacts: THM-2013, HYP-8762, reflection coordinates-for-the-continuum-cyclic-temperature-and-the-cycle-spectrum-boxeph-S199.md, script continuum_coordinates_boxeph_S199.py (+.out).

