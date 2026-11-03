# Exact finite decision for a tournament path count

For a tournament T let H(T) be its number of directed Hamiltonian paths. Given a positive integer h, H(T)=h for some finite tournament T is decidable by finite enumeration. More precisely, if h>1 is attained then it has a witness of order n satisfying 3^n <= h^3. If h is prime, it has a strongly connected witness of order n with f(n)<=h, where for q>=1

    f(3q)=3*5^(q-1), f(3q+1)=5^q, f(3q+2)=9*5^(q-1).

Here f(n) is the known minimum H over strong n-tournaments, n>=3. The bound is sharp for arbitrary singleton-free witnesses: the order-join of k directed triangles has n=3k and h=3^k. This is *not* a claim that it is a minimum-order witness for every 3^k.

## Proof and dependencies

The strongly connected components of a tournament are linearly ordered by domination. A directed Hamiltonian path visits the components consecutively in that order, hence H(T)=product H(C_i). Conversely, the order-join of tournaments realizing any given factors realizes their product. This exact product property is already recorded in repository THM-2450, section 2, and THM-4741. A singleton component has path count 1 and may be deleted without changing H. Every remaining component has order n_i>=3: no 2-vertex tournament is strong. Each has H(C_i)>=f(n_i)>=3^(n_i/3), by the known Busch/Moon bound stated explicitly in THM-4744, Setting. Indeed the latter inequality is true at n_i=3,4,5 and propagates by f(n+3)=5f(n) and 5>=3. Multiplying gives h>=3^(sum n_i/3), and cubing yields the entirely integer condition 3^n<=h^3.

For prime h, the product law leaves exactly one non-singleton component; removing the singleton components gives a strong witness. The theorem does **not** assert that minimum-order witnesses are strong for arbitrary composite h. Strong components have H>=3, so an exact alternative is to enumerate the finitely many strong tournaments with H dividing h and n satisfying f(n)<=h; then decide whether their H-values multiplicatively generate h. Both f and all bounds can be evaluated with integer arithmetic. For h=1 use one vertex. For even h, Rédei's theorem says every tournament has an odd number of Hamiltonian paths; otherwise the same all-tournament brute-force bound also works without using parity.

Citations and scope: the exact strong minimum is cited in THM-4744 as A. H. Busch, Electron. J. Combin. 13 (2006), N3, Corollary 1, with the bounds attributed also to Moon (1972). Component multiplicativity and the cited minimum are pre-existing results; this decision corollary is an elementary packaging, not a priority or novelty claim. This establishes individual decidability, not a classification of all h or an efficient algorithm. Direct exhaustive enumeration grows superexponentially in log h.

## Reproducible examples

`python3 decision.py` executes an independent labelled-tournament enumeration with Held–Karp subset path counting. It enumerates each possible strong factor d|h only up to the last n with f(n)<=d, using exact Python integers. Edge bit k corresponds to the kth pair (i,j), i<j, in lexicographic order; a 1 directs i->j and a 0 directs j->i. Output:

    strong witnesses by divisor: {13: (5, 41)}
    target 13 global bound 7 result ((13, (5, 41)),)
    strong witnesses by divisor: {3: (3, 2), 7: None, 21: None}
    target 21 global bound 8 result None

Thus the order-5 mask 41 is an explicit strong 13-path witness; the complete factor/divisor search proves 21 absent, using the cited strong minimum. The latter is a reproducible independent finite check of an already-known result (THM-1370); THM-4744 also records 13 in its complete strong inventory and 21 absent. For 21 the only potential nontrivial factors are 3, 7 and 21; the script excludes strong 7 and 21 and 21 is not a power of 3. The latter search enumerates strong tournaments only through n=6, since f(7)=25>21. The search is not a new finite global spectrum bound.
