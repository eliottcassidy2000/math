# Two independently widened transitive modules at exponent 61

**FINITE-EXACT.** Take any strongly connected tournament on at most eight vertices and replace any two distinct vertices by transitive tournaments of any positive orders. None of the resulting ordinary tournaments has either of the counts \(2^{60}(2^{61}-1)\pm1\). Thus a strong witness built by exactly two arbitrary transitive widenings needs a quotient on at least nine vertices. Fixing one width to one also covers a single widening. The search tests the scalar heights, so it applies regardless of where an edge is subsequently marked. [Simultaneous-module growth and its excursion obstruction](TRANSITIVE-PORT-FILTER.md) · [Earlier arbitrary-width four-block marked search](../chanpak_e61_fibonacci_target_20261122/SEED-LIFTS.md)

Let the two variable vertices be \(P\to Q\), and list the other vertices in their order on a Hamilton path. Within each gap in this list, including before the first and after the last vertex, the compressed word on the two modules is empty, \(P\), \(Q\) or \(PQ\). The arcs on the two sides determine which of these four choices are allowed. Multiply the gap polynomials, giving weights \(1,z,w,zw\) to the four options, and sum over all permutations of the ordinary vertices. If \(C_{rs}\) is the coefficient of \(z^rw^s\), the exact height at positive widths \(a,b\) is
\[
     \sum_{r,s\ge1}C_{rs}J(a,r)J(b,s),\qquad J(a,r)=r!{a\brace r}.
\]
Surjecting the vertices in each transitive module onto its successive visits determines the increasing path in every run. Each constructed path has a unique gap decomposition. The largest visit count here is seven. The singleton coefficient equals the path count of the quotient, which the program also computes independently by subset recursion. At every quotient through order six it also expands the actual rows at widths (2,2), (1,3) and (3,1) and compares the two methods.

The quotient has a Hamilton cycle. Starting and ending its traversal at \(P\) gives at least one skeleton for \((r,s)=(2,1)\), and breaking at \(Q\) gives one for \((1,2)\). Therefore the height is strictly increasing in each positive width and at least \(2^a-2\) and \(2^b-2\). For \(E_p\pm1\), widths above \(2p-1\) are impossible. A staircase walk on the integer grid tests every point that could equal either endpoint. Nonnegative arithmetic may be capped at \(E_{61}+2\) without affecting these comparisons.

| Quotient order | Strong isotypes | Choices of the two modules | Lower hits | Upper hits |
| ---: | ---: | ---: | ---: | ---: |
| 3 | 1 | 3 | 0 | 0 |
| 4 | 1 | 6 | 0 | 0 |
| 5 | 6 | 60 | 0 | 0 |
| 6 | 35 | 525 | 0 | 0 |
| 7 | 353 | 7,413 | 0 | 0 |
| 8 | 6,008 | 168,224 | 0 | 0 |

These are choices on the archived unlabelled representatives, before quotienting the choices of modules by automorphisms. The orders below three give no pair in a strong quotient. From the repository root:

~~~sh
c++ -O3 -std=c++17 05-knowledge/results/chanpak_decimal_kernel_20261122/two_variable_ports.cpp -o /tmp/two-variable-ports
for n in 3 4 5 6 7 8; do /tmp/two-variable-ports session212/cycle-levels/classes$n.rows 61; done
~~~

[Exact inverse](two_variable_ports.cpp) · [Quotient catalogue](../../../session212/cycle-levels)

A concurrent positive control shows what one growing port can do at a different target shape. A specific strong quotient gives physical flips at \(2^m+3,2^m+5\) for every \(m\ge6\), with the same marked arc. Its two-visit skeleton count is 32. The Euclid-center test here asks for a second moving binary scale, \(2^{p-1}\), in addition to the main power \(2^{2p-1}\). [Exact dyadic flip family](../chanpak_absolute_descent_20261122/DYADIC-BANDS.md) · [Two-scale obstruction for uniform inflations](TRANSITIVE-PORT-FILTER.md)
