"""Auditor's check of HYP-9540 (A)'s 'why it should hold' remark: at n = 3k + 1 (L = 3k - 1 positions) the tight
colourings are exactly the words constant on the blocks of some decomposition of [0, L-1] into k - 1 consecutive
blocks of 3 and one block of 2.  Brute force for k <= 7."""


def P(xs):
    c, last = 0, -10
    for a in xs:
        if a - last >= 3:
            c, last = c + 1, a
    return c


for k in range(1, 8):
    L = 3 * k - 1
    tight = set()
    for w in range(1 << L):
        A = [a for a in range(L) if w >> a & 1]
        B = [a for a in range(L) if not w >> a & 1]
        if P(A) + P(B) == k:
            tight.add(w)
    blockw = set()
    for pos in range(k):                       # the 2-block is the pos-th block
        sizes = [3] * k
        sizes[pos] = 2
        starts = [sum(sizes[:i]) for i in range(k)]
        for beta in range(1 << k):
            w = 0
            for i in range(k):
                if beta >> i & 1:
                    for a in range(starts[i], starts[i] + sizes[i]):
                        w |= 1 << a
            blockw.add(w)
    print(f"k = {k}, n = {3 * k + 1}: tight colourings {len(tight)} (x_0 = 1: {len(tight) // 2}); "
          f"block words {len(blockw)}; equal: {tight == blockw}")
