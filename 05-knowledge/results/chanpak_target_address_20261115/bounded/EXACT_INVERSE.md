# A target-first exact inverse for every one-apex tournament

**Proved:** for any prescribed integer N, decide whether a one-apex
tournament has exactly N Hamilton paths and reconstruct one if it
exists. The inverse uses O(sqrt(N) log N) arithmetic steps, rather
than evaluating O(N) words. It gives an exact exclusion certificate
for this particular graph class, not for arbitrary tournaments.

A one-apex tournament is obtained by adding one vertex to a
transitive spine. Its binary word records the 1-sites beaten by that
vertex. Strip any leading 0 and any trailing 1: they are singleton
sources and sinks and do not change H. If a directed triangle
remains the result is a strong word of length m beginning with
1 and ending with 0. Reject even N, handle N=1 and N=3 directly,
and assume N>3 below.

Use the exact inherited product estimate explained in
[TARGET_PLACEMENT.md](TARGET_PLACEMENT.md). With w and its reverse
complement regarded as m-bit integers A and B,
|H−AB/2^m| ≤ J with J=floor(m/2)+1. Hence H<2^(m+1).
The mandatory inversion between the first and last bit contributes
2^(m−2) to H−1, so H>2^(m−2). If k=floor(log_2 N), only the
three lengths m=k,k+1,k+2 need to be checked (m≥3).

Write m=2r+e, L=2^r, and choose the leading r bits a of A
and leading r bits b of B. Since the inputs are strong,
L/2≤a,b<L. Their two trailing fractions α,β, including any
middle bit, belong to [0,1), and

~~~text
AB/2^m = 2^e (a+α)(b+β).
~~~

Thus an exact hit must satisfy

~~~text
2^e ab ≤ N+J
2^e (a+1)(b+1) > N−J.
~~~

Given a, these are exact rational upper and lower bounds on b.
Test only the b values in that interval with L/2≤b<L. Its length
is at most a constant: before rounding and clipping, it is bounded
by

~~~text
N/[2^e a(a+1)] + (J/2^e)(1/a+1/(a+1)) + 2.
~~~

Indeed N<2^(m+1)=2·2^e L², a≥L/2, and J=O(log L). This
proves a uniform O(1) branching bound without a numerical gap
heuristic.

Now the binary word is determined by a, the reverse complement
of b, and the middle bit if e=1. Form these one or two literal
candidates and compute each H by the exact inversion recurrence.
This step, not the product bound, decides equality. The algorithm
therefore visits O(sqrt(N)) words and spends O(log N) arithmetic
operations on each. Constructing the arcs for a witness is
straightforward once its apex contact word is known.

Run python exact_inverse.py N to get an apex word or an exact
negative for this class. Running with no arguments checks every
odd target through 1023 against all strong words in the required
length range, then checks inverse reconstruction of selected
larger actual outputs. This strategy may complement a different
exact solver when a target has no one-apex witness.
