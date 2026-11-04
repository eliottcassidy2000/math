# Strong tournaments without a strong-root-deletion assumption

**FINITE-EXACT continuation, independently replayed and audited.** This addendum
strengthens the cap44 result of PROOF.md. All definitions,
OCF normalization, parity and inherited source attribution remain unchanged.

Let (T,v) now mean any rooted strong tournament of order at least three;
T-v need not be strong. If |T|>=4, Moon gives at least two vertices whose
deletion leaves T strong. Choose one of them other than v. By the OCF
inequality Delta(T-w,v)<=Delta(T,v), every such rooted object with Delta<=44
has a nested rooted strong chain down to C3, staying within the bound.

Thus a complete census starts from the single rooted C3 class and iterates
all nonconstant cuts, without the strong-root-deletion filter. The pipeline
`run_allroots.py` uses all eight labelled order-three tournaments for its
initial bank. Its exact class counts are

| order | rooted strong classes with Delta<=44 | minimum Delta |
|---:|---:|---:|
| 3 | 1 | 2 |
| 4 | 4 | 2 |
| 5 | 24 | 4 |
| 6 | 194 | 6 |
| 7 | 570 | 10 |
| 8 | 494 | 16 |
| 9 | 175 | 26 |
| 10 | 0 | none |

At order nine the cap44 bank happens to coincide with the restricted strong
deletion-pair bank. All 89,250 nonconstant ears are retained for evaluation
at the final step; their minimum Delta is46. In particular the previous
88,900 count belonged only to the stronger root-deletion filter and must
not be used as this universe's count.

The independently validated stronger conclusion is:

* The minimum of H(T)-H(T-v), over all rooted strong T of orders3,...,10,
  is 2,2,4,6,10,16,26,46.
* For every strong T on at least ten vertices and every v,
  H(T)-H(T-v)>=46.

For a general tournament, the difference factors as H(outside C) times the
marked difference in the strong component C containing v. Singleton C gives
zero. This follows directly from ordered strong-component factorization of
H and retains the distinction between ambient order and cyclic core order.

Reproduce independently of the restricted bank:

    g++ -O3 -std=c++17 rooted_ears.cpp -o rooted_ears
    python3 run_allroots.py 44 10

The optional final argument `allroots` passed to the executable disables
only the strong-root-deletion test. The older invocation with just parent
order and cap keeps its previous restricted behavior unchanged. The program generates banks, raw presentations, extension logs and summary in
`allroots_cap44/`. This addendum does not change the global spectrum status.

## A general finite bound for every increment cap

There is also a simple all-order bound, independent of these enumerations.
For a strong rooted tournament of order n>=3, repeatedly choose a strong
deletion other than the root. This gives a rooted strong induced tournament
of every order k=3,...,n. For each odd k, Camion's theorem supplies a
Hamiltonian cycle of that strong subtournament, hence an odd directed cycle
through v of length k. Different lengths are distinct cycles. Formula(2)
of the main proof, with every residual H at least1, consequently gives

    Delta(T,v) >= 2 floor((n-1)/2).                          (A1)

This is the familiar vertex-pancyclic lower bound, derived here just from
the inherited Moon deletion theorem and Camion's theorem. It is not claimed
as a new cycle theorem. For an even cap B>=2, Delta<=B forces n<=B+2.

Thus the marked-increment ear construction is a finite, complete procedure
for classifying **all** rooted strong tournaments below any fixed increment
cap. It needs no H bound and no empirical deletion-ratio inequality. The
cap44 census improves the crude order bound46 to the sharp last surviving
order9. This explains why a fixed additive vertex-arrival step cannot
continue through strong tournaments at arbitrarily large order, while
singleton SCC padding remains count-neutral.
