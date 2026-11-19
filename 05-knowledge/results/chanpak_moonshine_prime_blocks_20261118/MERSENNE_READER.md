# Cyclic Mersenne readers and an unbounded source requirement

For a nonempty tournament \(T\), write \(h=H(T)\) and \(q=P_2(T)\), counting ordered path covers with empty slots allowed. Form the cyclic composition on \(A=T\), \(B=TT_r\), and singleton \(C\), where \(r\ge1\). The block word of a Hamiltonian path is a consecutive segment of \((ABC)^\infty\) containing both \(A,B\) and exactly one \(C\). The three words of length 3 contribute \(3h\). The two words of length 4 contribute \(q-2h\) and \((2^r-2)h\); the only length-5 word, \(ABCAB\), contributes \((2^r-2)(q-2h)\). Thus this strong tournament realizes
\[
 N=2h+(2^r-1)(q-h)=(3-2^r)h+(2^r-1)q.
\]

## Absolute inverse

Choose \(r\ge1\), then \(h\) satisfying \(2h\equiv N\bmod(2^r-1)\). A physical source is required with
\[
 P_2(T)=h+\frac{N-2h}{2^r-1},\quad H(T)=h.
\]
Since \(q\ge2h\), necessarily \(h\le N/(2^r+1)\): this is a descending target coordinate, but existence of the joint profile does not follow from realizability of \(h\).

For \(r=1\), the target is \(h+q\). The odd-cycle packing identity gives \(q\equiv(-1)^n h\bmod3\): replacing \(G_T(2,1)\) with \(G_T(4,2)\) negates both variables modulo 3, and the number of blocks of a spanning packing by odd cycles and singletons has parity \(n\). Thus a prime target \(N>5\) requires even source order \(n\ge4\). In that case every contribution \(2^{2c+s}\) to \(q\) has exponent at least 3, where \(c\) is the number of cycles and \(s\) the number of singletons. So \(q\equiv0\bmod8\), and the combined filter is \(h\equiv17N\bmod24\). The \(r=2\) route \(3q-h\) has no such order-parity requirement. The \(r=3\) route \(7q-5h\) can hit a sevenfold target only when \(7\mid h\).

## Thinness of small-count sources

For any prime output \(N>3\), the source cannot have a directed triangle SCC. Such an SCC has \((h,q)=(3,12)\), and ordinal sums multiply both coordinates (each slot splits into one path in each component, possibly empty). More generally any common divisor of \(h,q\) divides the output. Every nonsingleton SCC therefore has \(H_i\ge5\), and Busch's bound gives \(n_i\le1+3\log_5 H_i\le4\log_5 H_i\). Write \(U\) for the ordered sum of nonsingleton SCCs. Then \(m=|U|\le4\log_5 h\). Removing \(k\) singleton components leaves \(h\) fixed and divides \(q\) by \(2^k\), independent of positions.

It follows that the number of prime outputs at most \(X\), for all Mersenne reader lengths and sources \(h\le B\), is at most
\[
 O((\log X)^2)2^{(4\log_5 B)^2/2+O(\log B)}.
\]
The exponential counts all labelled possibilities for \(U\); for each, there are \(O(\log X)\) choices each of \(r,k\) since \(N\ge q\ge 2^k\) and \(N\ge2^r+1\). If \(\log_2 B\le(\sqrt2\log_2(5)/4-\epsilon)\sqrt{\log_2 X}\) for fixed \(\epsilon>0\), this is \(O(X^{1-\delta})\) for some \(\delta>0\). Thus an arbitrary bank of Mersenne lengths and ordinal padding cannot reduce prime coverage to sources of bounded or even sufficiently slowly increasing count. This does not exclude larger sources or other readers.

Run: python mersenne_reader_replay.py (or python -O). It checks 4,396 literal Hamiltonian-path DP comparisons through labelled source order 5, with \(1\le r\le4\).
