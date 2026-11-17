# One-apex two-tag seeds: exact classification when the apex is tagged

Let F be a strong one-apex tournament of spine length L≥2, with the usual forward transitive spine (0,...,L−1). Use bits b_i=1 if the apex a beats i; strength is equivalent to b_0=1 and b_{L−1}=0. The normalized marked-seed condition from the multiple menu is U(A)+E(A)=H(F), where E counts full paths ending in A and U counts ordered two-path covers with first path ending and second path starting in A. Our tag set is A={a,j}, a single spine vertex j in addition to the apex.

For an internal j split the word P|b|S, with b=b_j. Let e be the last bit of P and c the first bit of S. Write f for the integer read in ordinary binary from P, and z for the reverse-complement integer read from S (first bit of S has weight one). The two orders of endpoints contributing to U are:

* first path ends apex, second starts j: z+1−e. All spine vertices before j are in the first path; optionally assign later vertices to it, counting their final possible zero bit.
* first path ends j, second starts apex: f+c. All later vertices are forced into the apex path; count its first possible tagged predecessor bit 1 in the prefix, or the fixed first later vertex.

The full end count at j is f+b (last spine vertex forced into the portion before apex), and E_a=1. Since U is even (oddness of the one-exit extension H(F)+2W−U proves this), normalization requires E_a+E_j odd, hence e=b. The target equation becomes
    H(F)=2+z+2f+c.

For concatenated words P,Q, write F(P) as the usual binary read and Z(Q) as reverse complement, and define H on a word as the full count of its one-apex tournament, even if not strong. Then
    H(PQ)=H(P)+H(Q)−1+F(P) Z(Q)+[last(P)=1,first(Q)=0].
By the one-apex odd-cycle expansion, any crossing inversion from position i in P to position j in Q contributes 2^(j−i−1), except adjacency, which contributes 2 and gets the indicated extra 1. Factoring that crossing sum gives the formula.

For e=b, insert P|b|S into the displayed count and compare it with 2+z+2f+c. The difference equation is
    b=1: H(P)+H(S)−2+2f(z−1)−2c=0;
    b=0: H(P)+H(S)−3+f(2z−1)−z−c=0.
Here P starts 1 and S ends 0. For b=0, P ends 0, hence f≥2 and H(P)≥3. Also z≥1, so f(2z−1)−z≥f−1≥1. No equality. For b=1, the condition forces z≤2. If z=2 then S=10, H(S)=3 and f≥1, again impossible. Thus S=0, H(P)=1, and P starting 1 forces all its bits 1.

For j=0, let S be the rest of the word and c its first bit. Let A0 be the first U order (apex ends its path, j starts the other); the other U order contributes c. Both endpoint counts are 1. Separating whether j is before apex gives H(F)=A0+H(S)−c, so equality would force the odd H(S) to equal 2c+2. For j=L−1, put P for the earlier word, e its last bit, f its ordinary binary read. The two U orders are 1−e and f+1, while the endpoint at j is H(P)−(1−e), giving U+E=H(P)+f+2. The actual H(F)=H(P)+f+e, impossible.

Thus the only normalized pair involving the apex consists of the apex and the penultimate spine vertex of 1^(L−1)0, since L=2 has no internal vertex. On this word g=H(F)=2^(L−1)+1, W_a=2^L−1, W_j=2^(L−1)+2. For W_j with apex in the other path, the path starting j is j or j,t, giving 2^(L−2)+2; if the apex is in the marked-start path, its order begins j,t,a and the subset of earlier spine vertices after a is free, giving 2^(L−2). The normalized seed therefore has h=W_a+W_j=3·2^(L−1)+1 and menu source p=3h+2g=11·2^(L−1)+5. This recovers the known transitive branch of the 22P+5u family, and proves that changing the one-apex word with these two tags cannot escape it.

Replay: `python replay_apex_pair.py` verifies the profiles by literal subset DP through spine length eight.
