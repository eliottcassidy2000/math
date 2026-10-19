"""Anchor collisions are trunk re-anchorings (corollary of THM-4640 (3)).
For the single-letter head (a+1): numerator K = 2^a - 3.  Another reduced head h (first letter >= 3,
total A_h, anchor image f_h(-1) = N/2^(A_h-1), N odd) has the same numerator K = 3N + 2^(A_h-1) iff
    f_h(-1) = a_k - 2^(1-A_h),   a_k = (4^k-1)/3,   a = A_h - 1 + 2k,  k >= 1.
We enumerate, for each a <= AMAX and each k, all reduced heads h with numerator N = 2^(A_h-1) a_k - 1
by a backward search (peel the last letter c: N' = (N - 2^(A_h-c-1))/3 must be an odd integer, or the
head is the single letter (A_h-1) when N = 2^(A_h-2) - 3)."""
import sys
from functools import lru_cache
AMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 44

@lru_cache(maxsize=None)
def heads(N, A):
    """all reduced heads (tuples) with anchor image N/2^(A-1): N odd, total A >= 2."""
    out = []
    if A < 2: return ()
    # single letter head (A-? ) : f_(a1)(-1) = -1/2^(a1-1): N = -1, a1 = A
    if N == -1 and A >= 2:
        out.append((A,))
    # last letter c, previous total A-c >= 2
    for c in range(1, A-1):
        Ap = A - c
        t = N - (1 << (Ap - 1))
        if t % 3: continue
        Np = t // 3
        if Np % 2 == 0: continue
        for h in heads(Np, Ap):
            out.append(h + (c,))
    return tuple(out)

def fval(h):
    from fractions import Fraction as F
    x = F(-1)
    for c in h: x = (3*x+1)/2**c
    return x

print("a : k : head h (reduced, first letter >= 3) with f_h(-1) = a_k - 2^(1-A_h)")
found = {}
for a in range(3, AMAX+1):
    for k in range(1, a//2 + 1):
        A = a + 1 - 2*k
        if A < 2: break
        ak = (4**k - 1)//3
        N = (1 << (A-1)) * ak - 1
        for h in heads(N, A):
            if h[0] < 3: continue          # canonical: first letter >= 3 (THM-4640: e_1 >= 2)
            if h == (a+1,): continue
            found.setdefault(a, []).append((k, h))
            x = fval(h)
            from fractions import Fraction as F
            assert x == F(ak) - F(1, 2**(A-1)), (a,k,h,x)
            # check the numerator identity: 3*N + 2^(A-1) == 2^a - 3
            assert 3*N + (1 << (A-1)) == 2**a - 3
for a in sorted(found):
    for k, h in found[a]:
        print(a, k, h, "total", sum(h), "D =", len(h)-1)
print("a with a second canonical representation of 2^a-3, a <=", AMAX, ":", sorted(found))
print("(THM-4640 (1) rigidity list for a <= 31: {7, 13, 17, 21, 28, 31})")
