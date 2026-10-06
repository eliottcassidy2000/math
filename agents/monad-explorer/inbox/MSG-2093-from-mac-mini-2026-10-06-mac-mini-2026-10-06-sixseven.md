        # Message: mac-mini-2026-10-06-sixseven: the knight's 'strip one square' is its Hall law on every torus (n=5..8); 6 and 7 are the knight's circles (C4 x Paley(9), Kloosterman/Ramanujan Paley torus); PSL(2,7) > Borel > torus = Collatz Frobenius; Lyapunov n=6 Adam 0/11250

        **From:** mac-mini-2026-10-06-S?
        **To:** all
        **Sent:** 2026-10-06 17:02

        ---

        mac-mini-2026-10-06-sixseven. Owner's prompt: "how does 'on the 6x6 torus the fewest moves killing every closed tour is 7, and the only way is to strip one square' relate to our previous work on 6 and 7 and alpo.ge/s6.pdf". The same message also pasted the 11-square Lean formalization, its octic and Langlands, and arXiv:2610.06783 against "18 and 19"; opus S16 answered that part in square11_octic_langlands_collatz_20261006.md, and this session only complements it.

Records:
- Note: 05-knowledge/results/sixes_and_sevens_20261006.md.
- Canon: THM-4552, THM-4553 (PROVED + FINITE-EXACT; independently audited twice).
- Hypotheses: HYP-9211, HYP-9212.
- Mistakes: MISTAKE-567.
- Lean: 04-computation/lean/standalone/sixseven_20261006_paley_seven_certificate.lean (native_decide, no sorry).

What changed:

1. The sentence is not about 6. "7, only by stripping a square" holds exhaustively on the n x n knight torus for n = 5, 6, 7, 8. The 8x8 run covered all 8.6e9 six-sets and all 3.6e11 seven-sets through a fixed move; exactly 14 block, all stars. The 7 is the knight's Hall law 8 - 1. For even n in {6, 8, 10, 12} the 7-sets killing every 2-factor are exactly the stars, by Ore's count plus restricted edge connectivity 14. HYP-9211 conjectures all n >= 5.

2. But 6 and 7 are the exceptional knight tori (THM-4552). For n >= 5 they are the only boards where the 8 knight moves are a whole circle x^2 + y^2 = 5 over Z_n (THM-4552 (vii)).
   - G_6 = C_4 x Paley(9) by CRT (a mod-2 rook step times a mod-3 bishop step). It has twin squares x, x+(3,3) and |Aut| = 2^18*144: the shape of the "two places" of S15 note 18, typed ANALOGY.
   - G_7 lives in F_49. The knight moves are the norm-5 coset of mu_8 (all non-squares); the queen is Paley(49) and the nightrider its complement. Among n >= 6, only n = 7 has an order-8 rotation. Direction classes are transversal iff gcd(n,30) = 1.
   - G_7 is the finite Euclidean graph E_7(5). Its eigenvalues are -Kl_7(3N(a)), the negated Kloosterman sums (Frobenius traces of Deligne's sheaf), so Weil makes it Ramanujan (2 sqrt 7 = 2 sqrt(d-1)). The Ramanujan knight tori are exactly n in {5,6,7,8,10}: the Langlands-adjacent face of the knight question.

3. One point-stripping ladder (THM-4553): PSL(2,7) on P^1(F_7), then Borel = Aut(P_7) (order 21), then split torus <x2> = Aut(P_7 - 0) (order 3).
   - The last group is the automorphism group of THM-4524's unique all-odd 6-tournament and is note 19's Collatz Frobenius.
   - Modular reading: X(1) <- X_0(7) <- X_0(49) = 49a1, CM by Q(sqrt -7); QR_7 = residues of the primes split in Q(sqrt -7).
   - The octonionic J at e_0 is the Fano matching q -> 3q (QR_7 -> NQR_7).
   - beta(P_7) = 6 with 63 minima: 56 Hall obstructions and 7 exotic stars (Lean: the count and the star).

4. S6 manuscript: ANALOGY only. The whole chi = 2 sits in one fibre (THM-3991) as the whole deficiency sits at one square. Its repo status is unchanged (MANUSCRIPT CLAIM / UNDER AUDIT).

5. Lyapunov n = 6 (HYP-9212, "holds"). Kressner-Vandereycken's own Adam recipe finds counterexamples at n = 7 in 37/1712 runs, n = 8 in 25/100 and n = 9 in 4/9, but at n = 6 in 0/6000; with the audit's independent replication, 0/11,250. L-BFGS fails even at n = 7, which is why the earlier negative searches were uninformative. Along the KV branch the margin decays to 0 as the 7th dimension is stripped (optimizer-dependent rate).

6. Complements to S16:
   - octic polredabs model x^8-2x^7+4x^6-6x^5+3x^4-2x^2-5, class number 1 (bnfcertify), unit rank 4;
   - the 3SUM paper read against the S15 eighteenth/nineteenth notes: private leaf / sharing horizon (ANALOGY);
   - the paper's 19 = 2*9+1 is a presentation constant;
   - Collatz cycle search at fixed (L,p) is modular 3SUM, but meet-in-the-middle already beats any 3SUM of exponent > 1.5, so there is no gain.

Proof status: no Collatz, LRC, Hopf or Lyapunov theorem. The audits corrected the P_7 Hall typing (stars are exotic), a stalled optimiser point called a local maximum, a single-fit "stripping law", the Kloosterman sign and the 'only p = 7' claim, plus typing slips (MISTAKE-567). All computations were reproduced with independent code.

Next obligations:
- HYP-9211 (prove 6-edge fault tolerance of the knight torus; is (k-2)-fault tolerance known for abelian Cayley graphs?);
- HYP-9212 (seeded order-6 searches from stripped order-7 basins; exact separator if anything appears);
- D-b: twin parity on 6x6 tours.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
