# Two exact edge families and growth

Let G_n be the simple graph on 1,...,n with an edge between distinct vertices exactly when their sum is a square. The following are all-order statements in the positive-integer graph; for G_n retain pairs whose larger endpoint is at most n.

For k≥1, consecutive vertices k,k+1 are adjacent if and only if k=2t(t+1) for some integer t≥1. Indeed 2k+1 must be an odd square, say (2t+1)². The first examples are 4–5, 12–13, and 24–25.

Vertices k,k+2 are adjacent if and only if k=2t²−1 for some integer t≥1, since 2k+2 must be an even square (2t)². The first examples are 1–3, 7–9, and 17–19. In particular, for even k the pair is never adjacent, as illustrated by 2–4.

None of these adjacent pairs belongs to a triangle. A common neighbor would require positive squares differing by 1 or 2; the first is impossible because consecutive positive squares differ by at least 3, and the second is impossible modulo 4. This does not assert that every triangle-free edge lies in one of these two families.

For n≥1, the number of older neighbors of n is exactly floor(sqrt(2n−1))−floor(sqrt(n)), because their square sums lie in n < s² ≤ 2n−1. Therefore E(G_n) is the sum of this quantity over 1≤m≤n. Replacing floors by their arguments incurs O(n), and summing square roots gives E(G_n) = (2/3)(sqrt(2)−1)n^(3/2) + O(n). For n≥14 the graph is connected, so its cycle rank is E(G_n)−n+1; the separate [growth proof](growth/REPORT.md) establishes connectivity and the exact unicyclic range. The asymptotic count does not settle the global Hamilton-path question.
