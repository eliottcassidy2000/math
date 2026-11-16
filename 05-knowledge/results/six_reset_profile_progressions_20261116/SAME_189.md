# 63 and 189: same output, different two-count clock

Take the strong 63-path tournament with outgoing rows
(102,252,248,241,225,192,128,1). Its two-slot count is 1872.
The other displayed birth representative
(16,5,17,7,74,31,175,63) has exactly the same pair (63,1872).
No SCC factorization can synthesize the scalar 63 itself from
nontrivial smaller counts, since the only possibilities involve
7 or 21. This uses the known holes, not a numerical exclusion.

The directed triangle has pair (3,12). Therefore its ordinal
concatenation with either 63 block gives pair (189,22464).
But on 7 labels the Paley tournament itself, with connection
set {1,2,4} modulo 7, has pair (189,2016). These are two
literal ways to feed the same 189 into the 756 reader:
the absolute counts differ by 756·2^s·20448 for the same
intercept I. Ordinal multiplication of a 63 to reach 189
is not interchangeable with a strongly connected Paley 189.
The latter's one-vertex deletions all have pair (45,432)
by translation symmetry, so even recovering the smaller
count 63 by vertex deletion is impossible by induced
monotonicity. Here 756=4·189 names a scalar, not a full
input to the next step.

Standalone check, reusing the direct subset DP in policy.py:

    from policy import read_profile
    a=(102,252,248,241,225,192,128,1)
    b=(16,5,17,7,74,31,175,63)
    p=tuple(sum(1<<((i+d)%7) for d in (1,2,4)) for i in range(7))
    q=tuple(sum((p[i]>>j&1)<<(j-1) for j in range(1,7))
            for i in range(1,7))
    assert read_profile(a)==read_profile(b)==(63,1872)
    assert read_profile((2,4,1))==(3,12)
    assert read_profile(p)==(189,2016)
    assert read_profile(q)==(45,432)
