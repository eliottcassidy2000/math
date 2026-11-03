# A closed reset family can still be too sparse for all targets

**PROVED for the declared balanced-return architecture and its adaptive
stage outputs, even allowing every block length k>=1.** Temporary interior
graphs of the stabilized native compiler are not included in this output
language. The result does not obstruct other root choices or source families.

The repeated construction in [NATIVE_SELECTION](NATIVE_SELECTION.md) takes
an actual strong balanced rooted profile (h,s,q,F)=(b,t,t,b) to

    M_k(b,t)=(A_k b−2t, B_k b),
    A_k=3·2^k+4·3^(k−1)+1, B_k=2^k+1, k>=1.

It uses a new-root birth, k births at that root, moves the retained leaf
to the oldest position in the homogeneous block, then makes new-root and
keep-root births. For k=1 there is no native move. The actual balanced
sources satisfy 0<t<=b/2: insert the root noninitially into every path of
its deletion, as proved in [GRAMMAR G1](GRAMMAR.md).

Consequently every return satisfies

    b_new >= (A_k−1)b
          = (3·2^k+4·3^(k−1))b >= 5·2^k b.         (B1)

This remains true when a primitivity controller rejects some choices;
allowing all lengths and all words gives an upper bound on the accepted
family. Fix an initial actual balanced source of count b0 and let σ=2/3.
The sum over letters of the reciprocal growth bounds satisfies

    S = sum_(k>=1) (5·2^k)^(-σ)
      = 5^(-σ)/(2^σ−1) < 4/5.                      (B2)

The last inequality uses only 2^(2/3)>3/2 and 5^(-2/3)<2/5.
Summing over all finite words, with repetitions of the same profile counted
separately, therefore gives

    sum_words b(word)^(-σ) <= b0^(-σ)/(1−S)
                           < 5 b0^(-σ).             (B3)

The declared adaptive and return outputs from a balanced parent (b,t) are

    b,
    2b−t,
    b(2^j+1),                       j>=1,
    2b(2^k+3^(k−1))−t,              k>=1.

The next balanced output is already included as a word node. Every module
permutation endpoint has the same count as its block stage. These counts
are bounded below respectively by b, (3/2)b, 2^j b, and 2^k b. For any
positive c, 1_(cb<=X) <= (X/(cb))^σ. Combining this with (B3) proves

    #{distinct declared outputs <=X}
      <= X^σ [1+(3/2)^(-σ)+2/(2^σ−1)] sum_words b(word)^(-σ)
      < 30 (X/b0)^(2/3).                            (B4)

Thus this architecture, even with unbounded block lengths, misses infinitely
many primes and infinitely many values 7p. The elementary central-binomial
argument in [GRAMMAR G5](GRAMMAR.md) gives pi(X)=Omega(X/log X), which
eventually exceeds (B4); apply it at X/7 for the second lane. Removing the
finite exceptions 7 and 21 does not affect the comparison. A finite union
of initial balanced seeds has the same O(X^(2/3)) bound.

This is a stopping reason for treating balanced renewal as the whole target
selector. Its positive achievements remain: an actual reusable source
family, exact primitivity guards, a finite branching controller, and legal
port transport. A complete induction must add an operation or an output
language that escapes the stated growth-versus-branching bound. In
particular, this proof does not cover retaining new roots between returns,
unbounded families of initial sources, or extracting counts at temporary
native intermediate graphs whose magnitudes have not been bounded here.

The proof uses inequalities and a convergent word sum, not a numerical
coverage scan. It complements [GRAMMAR G5](GRAMMAR.md): there a fixed
number of singular steps is too restrictive; here a permanently balanced
return architecture is too restrictive. Neither a finite observer nor a
closed source family alone supplies all-target coverage.
