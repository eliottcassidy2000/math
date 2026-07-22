## boxeph-2026-07-21-S215 -- each prime IS its Paley tournament: a periodic table of 2,3,5,7,11 for LRC(14) (HYP-8860)

**Owner:** understand primes 7,5,11 as well as 2,3, through 'a set of pairwise relations is a tournament'.

**UNIFYING LAW (verified, understanding_primes_via_paley_tournaments_boxeph_S215.py):** prime p IS its Paley object (i->j iff j-i is a QR mod p). p=3 mod4 (3,7,11) -> Paley TOURNAMENT: self-converse, vertex-transitive, doubly-regular, char_A=(x-(p-1)/2)(x^2+x+(p+1)/4)^((p-1)/2), quad disc -p, roots (-1+-i sqrt p)/2 = quadratic Gauss sum i sqrt p. p=1 mod4 (5,13) -> self-complementary Paley GRAPH, real spectrum (-1+-sqrt p)/2 = Gauss sum sqrt p. p=2 = the reversal INVOLUTION (Phi_2=x+1, chirality, S210-S213). Verified: Paley-3=(x-1)(x^2+x+1)=Eisenstein; Paley-7=(x-3)(x^2+x+2)^3 (THM-1830); Paley-11=(x-5)(x^2+x+3)^5; |g_p|^2=p (i sqrt p for 3,7,11; real sqrt p for 5,13).

**PERIODIC TABLE for LRC(14)=2*7:**
- 2 = the INVOLUTION (chirality/reversal iota:t->1-t, S212/S213; Phi_2 the antipode).
- 3 = the ATOM: Paley-3 = the 3-cycle, char x^2+x+1 = Eisenstein; the argmax Phi_6(14)=183, t*=14/183.
- 5 = the GOLDEN FOIL (=1 mod4): self-comp GRAPH, real sqrt5 = Q(sqrt5) = Fibonacci = the LRC loosest/foil (S206); also Bonferroni depth.
- 7 = the APEX (14=2*7): Paley-7 tournament, i sqrt7 = my S212 Euler-branch index; Phi_7/Phi_14 hardness; F_7[C_14]=F_7[X]/(X+-1)^7 (THM-2043, verified x^14-1==(x-1)^7(x+1)^7 mod 7); cap field Q(cos2pi/7) cubic disc 49.
- 11 = the RANK (=3 mod4): Paley-11 (i sqrt11), the rank-11 AP-core/relation code (S214); SCARCE (only multiple<=14 is 11) => forced rigid speed.

**SYNTHESIS:** p mod 4 decides tournament(i sqrt p)-vs-graph(sqrt p), tight(Eisenstein 2/3/7)-vs-slack(golden 5); the Gauss sum IS the Paley spectrum, fixing each prime's role. Ties S206(foil)+S212(i sqrt7)+S213(chirality)+S214(rank-11)+opus-S434(Paley symmetric-intransitive pole, opposite the transitive AP nullcone vertex).

**Honest:** all rows are verified/classical facts; the contribution is the SYNTHESIS -- one Paley construction giving each small prime a tournament personality via p mod 4, mapped onto LRC(14)'s tight/slack/apex/rank structure. Artifacts: reflection each-prime-is-its-paley-tournament-...-boxeph-S215.md, HYP-8860, script (+.out).

