"""Auditor: which nu-pairs can a BLUE tight line have?  A blue line's triple word is a palindrome (sigma reverses it),
and a tight line's word is a tight colouring with nu = (P(A), P(B)) (THM-4951 (1)).  So the blue tight nu-pairs are
contained in {(P(A), P(B)) : palindromic tight colourings} -- with no use of HYP-9540 (C).  Listing this set shows
where HYP-9540 (B)'s blue clause (every pair i >= 1, i + j >= k, (i, j) != (1, 1)) must fail on the boundary i + j = k."""


def P(xs):
    c, last = 0, -10
    for a in xs:
        if a - last >= 3:
            c, last = c + 1, a
    return c


for n in range(3, 27):
    L, k = n - 2, n // 3
    pairs = set()
    npal = 0
    half = (L + 1) // 2
    for h in range(1 << half):
        w = [(h >> i) & 1 for i in range(half)]
        word = w + w[:L // 2][::-1]
        if word[0] != 1:
            continue
        A = [a for a in range(L) if word[a]]
        B = [a for a in range(L) if not word[a]]
        if P(A) + P(B) == k:
            npal += 1
            pairs.add(tuple(sorted((P(A), P(B)))))
    allowed_B = {(i, k - i) for i in range(1, k // 2 + 1) if not (i == 1 and k - i == 1 and n >= 6)} | {(0, k)}
    missing = sorted(allowed_B - pairs)
    print(f"n = {n:2d}, k = {k}: palindromic tight colourings (x_0 = 1) {npal:4d}; possible blue tight pairs "
          f"{sorted(pairs)}; boundary pairs allowed by HYP-9540 (B) but impossible: {missing}")
