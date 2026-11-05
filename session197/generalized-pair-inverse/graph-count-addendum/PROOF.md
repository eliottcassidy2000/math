# One generalized pair gate cannot describe every C2 graph

Session197. This note is a structural limitation on the separately frozen generalized N/M/P/Q families. It says nothing negative about their possible scalar Hamilton-path coverage.

Let Gamma be their unrooted graph image: choose any rooted start-one source V, form one of the four exact gates, and append the specified single suffix/gadget. N/M/Q require |V|>=2; P permits a singleton. Different presentations may give isomorphic graphs.

The reviewed rooted start-one enumeration gives

    B(z)=z+z^3/Delta(z),  Delta(z)=1-3z+z^2,

for strong sources including the singleton. Every weak start-one source uniquely consists of such a rooted core followed by a nonempty chain of singleton SCCs. Thus the exact rooted-source generating function, including both strong and weak sources, is

    W(z)=B(z)/(1-z)=z(1-2z)/Delta(z).                (1)

This counts rooted isomorphism classes, not unrooted source graphs. The SCC description makes the tail length and the core root recoverable, so no multiplicity is introduced in (1).

The binary retained suffix at parameter1 is empty. At length m>=1 its last digit is1, with the earlier m-1 digits arbitrary. Therefore its presentation generating function is

    S(z)=1+sum_(m>=1)2^(m-1)z^m=(1-z)/(1-2z).       (2)

An ordinary gadget indexed by a positive integer has the same binary tail choices and one extra root, hence generating function zS(z).

N adds two gate vertices and a nonempty suffix. M adds three and any suffix. Q adds two and one ordinary gadget. P adds four and any suffix. Excluding only the singleton source for N/M/Q gives the exact count of presentations

    P_Gamma(z)=(W-z)[z^2(S-1)+2z^3 S]+W z^4 S
              =4z^5(1-z)^2/[Delta(z)(1-2z)].        (3)

The word 'exact' here applies to presentations, not graph classes. Every Gamma graph has at least one presentation. Rooted-isomorphic sources with the same family and parameter yield isomorphic outputs. Consequently [z^n]P_Gamma is an upper bound on the number of unrooted Gamma graphs of order n.

The independently proved complete C2 enumeration is

    A_2(z)=z^5+z^5(1-z)(1-z^2)/Delta(z)^2.           (4)

At order14 the exact coefficients are

    [z^14]P_Gamma=41,736,
    [z^14]A_2=44,394.

Thus at least2,658 strong C2 isomorphism classes of order14 have no generalized N/M/P/Q graph presentation. This is a counting proof, not a failure of a bounded search for presentations. It disproves any proposed graph-isomorphism-preserving absorption of arbitrary ordinary continuations into one of these four one-gate forms.

There is also an asymptotic obstruction. Let alpha=(3-sqrt(5))/2 and rho=1/alpha=(3+sqrt(5))/2. Delta has distinct roots alpha and1/alpha. In (3), the smallest-modulus pole is alpha and is simple; the other denominator root1/2 has larger modulus. Its numerator is nonzero there. Therefore its coefficients are O(rho^n). In (4), alpha is a double pole with positive nonzero leading coefficient, since alpha^5(1-alpha)(1-alpha^2)>0. Its coefficients are Theta(n rho^n). These statements also follow directly by partial fractions of the rational functions, without an analytic-combinatorics theorem. Hence the fraction of C2 graph classes admitting a Gamma presentation tends to zero, and missing graph classes occur at every sufficiently large order.

No scalar conclusion follows. A graph outside Gamma may have exactly the same H as a different graph inside Gamma, possibly of another order. In particular the result does not refute cofinite scalar coverage by Gamma together with the C1 class. It only rules out a tempting graph-level shortcut toward that coverage claim.

Dependencies: generalized-pair-inverse/PROOF.md SHA9df4e06c88d28e159953601611bbcc6fb90f00f6766b5c09b03b3dfc7c76cbd3 (family domains and vertex increments); FIBONACCI_CLASS.md SHAd95a133e360b96589009ebd8a1c8449525c4a2847a88b33a4b3f2b839739e63e (rooted strong-source enumeration); C2_COMPLETE.md SHA4cc0bf11efc435d7a2746daf6d8f1a58c49c94e57e608911b24464facdd1057d (complete C2 enumeration). The source enumeration is unaffected by the historical narrower unrooted scope in FIBONACCI_CLASS.md. The generalized-family proof is undergoing separate review; this count bound is conditional on its stated construction definitions, not on its scalar identities or inverse implementation.
