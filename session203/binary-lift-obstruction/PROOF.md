# A rigorous obstruction to bounded-offset binary lifting

Status: PROVED, independently audited. See audit/REVIEW.md.

This is a limitation on one natural certificate format, not a counterexample to F(8,137), not a nonexistence theorem for arithmetic induction in every radix, and not a limitation on the complete finite-register weighted presentation.

## 1. Domain and inherited facts

Let G be the physical class generated from the root-first directed edge by any finite word in O,Z,J,K1,K2. This includes the complete start-one and start-two grammars. Allowing repeated gates is harmless for the result, but no completeness claim is made for their higher-start outputs. Every input is nontrivial and has positive root-start count.

Use p=(h,s,e,L,R,d,C,b,T), B=h+L-s. The published positive-profile-defects packet, commit a67c5f829aba19285bd10a1dbf665171b9c21e14, proves on G:

    g=C-3L-h>=0;
    B_O=B;
    B_Z>2B, B_J>2B, B_K1>2B, B_K2>2B;
    D_J:=B_J-2B=g+(R-e)+s>=h-e+s;
    D_K1>=h, D_K2>=h.

The last line follows from its displayed nonnegative decompositions: in K1, 3R-h-e>=h because R>=h>=e; in K2, h is an explicit summand and the others are nonnegative. Also

    B_Z=B+R-e+b,
    B_J=h-L-s+R+C-e.

The exact Z formula is B_Z=2B+delta+2p0, where p0=d-e>=0 and delta=s+f0+q0>=1. Z preserves delta and resets p0 to zero. Thus every Z ray after its first step satisfies B_next=2B+delta with fixed positive integer delta.

The proof uses the standard Hamilton-cycle theorem for nontrivial strong tournaments, and the inherited weak-tail profile factorization: if the root is in the initial SCC A, and later components have Hamilton-count product q, then B=qB(A), h=qH(A), s=q s(A). Singleton later components may be removed without changing B,s. No bound on the number of singleton tail vertices is asserted.

## 2. Bounded J or gate defects force bounded B

For every n-vertex rooted tournament, h<=n! and L<=n!. For L, concatenate its designated root-ending path and its ordinary path. The position of the root uniquely determines the cut, so this is an injection into vertex permutations. Hence B<=2n!.

**Lemma 1.** For every positive-start rooted tournament, h<=c implies B<=2c!, for any positive integer c.

If the root SCC A has order n and the later Hamilton-count product is q, then n<=H(A), including the singleton case. Therefore qn<=qH(A)=h<=c. Also B=qB(A)<=2q n!. Since n<=c,

    q n! <= c (n-1)! <= c!.

This proves the assertion even when arbitrarily many singleton SCCs are present.

**Lemma 2.** On G, D_J<=c implies B<=2c!.

If U is weak, the root is in its initial SCC, so no spanning path ends there: e=0. The inherited inequality D_J>=h-e+s gives h<=c, and Lemma1 applies.

If U is strong, a Hamilton cycle gives n distinct spanning paths by cutting its n arcs, exactly one ending at the root. Thus h-e>=n-1, and D_J>=h-e+s>=n since s>=1. Consequently n<=c, giving B<=2c! directly. The singleton is not in G.

**Corollary.** On G, each of D_J,D_K1,D_K2 tends to infinity with B in the precise uniform sense

    D_op<=c  implies  B<=2c!,    op in {J,K1,K2}.

For the gates use D_op>=h and Lemma1. No assertion about an asymptotic rate stronger than this explicit factorial bound is needed.

## 3. Prefix sinks prevent near-doubling

For a>=1, O^a preserves B,h,L,s,d and sends R to 2^a R, b to 2^a b, C to 2^a C, and e to zero. Therefore

    B(Z O^a U)=B+2^a(R+b)>3B.

Here Z O^a means first O^a, then Z. The strict inequality uses R+b>B, which follows from B_Z=B+R-e+b>2B and e>=0.

Although not needed below once J is excluded, the analogous identity gives

    B(J O^a U)=B_J(U)+(2^a-1)(R+C)+e>3B,

because R+C>=h+2L>B.

Two or more nonneutral letters in any physical word give output B>4B_input, regardless of any neutral O letters between or around them.

## 4. Exact classification of large bounded-offset binary steps

Fix an integer c>=1 and set H=2c!. Call a physical word w a c-bounded binary lift from U if its output satisfies

    B(wU)=2B(U)+epsilon,        |epsilon|<=c.

This definition allows any finite word over the five physical operations, not merely a fixed finite menu.

**Theorem 1.** If U belongs to G and B(U)>H, every c-bounded binary lift has the form

    w = Z followed by O^t,       t>=0.

The empty or O-only words cannot lift, since B>c. A word with at least two nonneutral letters has output >4B>2B+c. Thus there is exactly one nonneutral letter. Its O-prefix preserves B. Lemma2 and its gate corollary exclude J,K1,K2 because the sole nonneutral defect would be at most c while the input B is greater than H. The remaining letter is Z. Section3 excludes a nonempty O-prefix. A terminal O run has no effect on output B and is the sole remaining freedom. In particular epsilon is positive.

**Theorem 2.** If such a lift has t>=1, its output admits no further c-bounded binary lift. If t=0, every indefinitely continuable chain is the deterministic Z ray, with one fixed positive offset delta<=c.

For the first assertion, the output B is >H, so Theorem1 says its next lift must start with Z. But this is Z applied after at least one O to the preceding Z output, hence Section3 makes its scalar output >3B_current>2B_current+c. This is impossible.

For t=0, Z has reset p0 to zero and preserved delta. Its first bounded defect is delta+2p0<=c, so 1<=delta<=c. All subsequent plain Z lifts have exact defect delta. Define B_0 to be the scalar AFTER that first large lift, when p0 has become zero. The subsequent ray therefore has the exact formula

    B_n=2^n(B_0+delta)-delta,       n>=0.

Terminal O variants share the same scalar value and then stop.

## 5. Finite-seed obstruction

**Theorem 3.** Fix c and finitely many seed graphs in G. Let S be the set of all scalar B values obtainable by any finite chain of c-bounded binary lifts, allowing an arbitrary physical word at every step. Then S has upper counting growth O(log M) below M and therefore density zero, even relative to odd integers. In particular it contains no odd tail.

Below the threshold H there are only finitely many scalar values. A lift from such a value can first cross H only into the finite integer interval (H,2H+c]. Seeds already above H supply only finitely many additional starting values. Theorem1 sends each large input through Z, with optional terminal O. After the first large transition, its ray has a positive integer offset delta<=c. Thus every infinite continuation is among a finite collection of rays determined by a finite scalar starting value and one of c possible offsets. The intermediate first large transition also has only finitely many possible scalar outputs. Terminal O variants never add a new B and cannot continue. Each ray has O(log M) values below M. This proves the result without needing finiteness of the set of physical states below H.

## 6. Consequence for one tempting F certificate format

In F(8,137), a certificate source has B=X-2j, with 0<=j<=8. Suppose one tries to build all target certificates from finitely many physical seeds by rules

    X_next=2X+d,       |d|<=D,

whose physical implementation is a word in O,Z,J,K1,K2 applied to the current source. Changes of j are allowed. Then

    B_next=2B+(4j+d-2j_next),

so every step is a (D+32)-bounded binary lift. Theorem3 shows that this entire finite-seed method has zero-density scalar output. Bounded shifts X=B+2j preserve the O(log M) bound. Therefore this certificate format cannot establish F(8,137) or cover an odd tail.

This does not prohibit higher-radix rules, target-dependent rebuilding that is not a bounded-offset binary lift of the preceding physical source, simultaneous multi-input constructions, additional physical operations, or the complete register-interface presentation. In particular finite dimensionality of the latter remains useful, but cannot by itself supply the missing arithmetic inverse.
