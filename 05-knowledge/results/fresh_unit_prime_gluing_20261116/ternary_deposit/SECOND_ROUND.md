# Exact second-round boundary at deletion count 45

Let \(P(T,r)=(H,E,O,g)\), where \(g=H(T-r)\) and \(E,O\) count
paths starting at \(r\) of the respective vertex parities, each
with weight the Hamiltonian path count of the omitted subtournament.
The empty tournament has count one.

For \(q=79,163\), consider all odd counts in \([9q,9q+16]\), with
required continuation \((378,423,45)\) and \((756,801,45)\),
respectively. The exact sets are \(\{9q,9q+6,9q+12\}\).
The independent literal replay establishes the stronger restriction
to **strong** hosts as follows:

| q | Strong count | contact presentations in the replay |
|---:|---:|---:|
| 79 | 717 | 5 |
| 79 | 723 | 8 |
| 163 | 1467 | 2 |
| 163 | 1473 | 10 |
| 163 | 1479 | 13 |

The presentation column need not count rooted isomorphism classes:
it enumerates all masks into one representative per deletion type,
without dividing by automorphisms of that representative.
Concrete rooted rows, with the root last, are in results.txt.

## Why the strong calculation suffices

Both continuation differences \(E-O\) equal \(-45\). If any
nonempty ordinal suffix follows the root component, its contribution
to this signed difference vanishes. In fact choose the subset \(S\)
of the suffix lying on the first path; it contributes
\((-1)^{|S|}H(S)H(S^c)\), up to a factor from the earlier components.
The sum over \(S\) is zero: in the positive odd-cycle packing expansion,
the union is a packing of disjoint odd cycles plus uncovered
singletons, and each odd block can be assigned to either side,
cancelling the signs. Thus the root is in the last strong component.
Its first path then never enters the earlier ordinal prefix, so the
Hamilton count of the prefix multiplies all four profile coordinates.
This is the terminal-component lemma. For either
continuation the gcd of its three coordinates is 9. Every excluded
offset \(d=2,4,8,10,14,16\) is prime to 3, so a host at \(9q+d\)
would have a prefix multiplier dividing both 9 and \(9q+d\),
hence one. The terminal host is then itself a strong host with
the same full profile. It cannot be a singleton with those H,E.
The three positive counts follow from the first-digit witnesses
with \(H=3q,3q+2,3q+4\) and triple first-round continuation: take
the ordinal sum with a directed triangle as prefix, which multiplies
all four coordinates by three.

## Exhaustive, independent calculation

1. [Busch's strong minimum bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) is \(f(3)=3,f(4)=5,f(5)=9\),
   \(f(n+3)=5f(n)\). Since \(f(14)=1125>727\) and
   \(f(15)=1875>1483\), the strong result needs final orders at most
   13 and 14 respectively; the deleted tournaments have orders
   at most 12 and 13.
2. The strong components of a deletion of count 45 have nontrivial
   Hamilton counts \(45\), \(15,3\), \(9,5\), or \(5,3,3\). This uses
   multiplicativity along the ordinal decomposition, and the strong
   count-one case is just a singleton. The code loops over every
   distinct order of the factors and all slots for isolated ordinal
   singleton components, up to the deletion order ceiling.
3. Every nonsingleton component in step 2 has order at most eight
   by \(f(9)=75\). Give it a Hamilton path (Rédei), label along that
   path, and construct its adjacency by appending vertices, requiring
   each predecessor-to-newcomer arc. Prefix counts above 45 can be
   dropped: induced enlargement cannot decrease the Hamilton path
   count by the positive odd-cycle packing identity. Remaining strong
   candidates with exact relevant total are canonized by all label
   permutations. The code finds 1, 1, 2, 4, 12 strong types for
   3, 5, 9, 15, 45 respectively and enumerates 926 ordinal deletions.
4. In each fixed deletion \(G\), compute \(H(G[S])\) by subset DP
   using paths beginning at each vertex. Separately for every start
   \(i\) compute the number of \(i\)-to-\(j\) paths spanning \(S\)
   by the path-ending recurrence: delete \(j\) and sum over possible
   immediate predecessors. Weight each by \(H(G-S)\). The sums with
   old-path odd and even parity, by starting vertex, give \(E\) and
   \(O-g\) when summed over outcontacts of the new root. For the
   latter include the singleton root path of weight \(g=45\).
   For \(H\), its difference from \(g\) is twice the weighted number
   of even old paths starting in an outcontact and ending in an
   incontact: these are exactly the odd cycles through the new root
   plus an arbitrary disjoint odd-cycle packing on the complement.
5. Strongness is equivalent to finding an outcontact in the first
   deletion SCC and an incontact in the last. Necessity follows from
   the ordered component cuts. Sufficiency follows by reaching all
   deletion vertices from that initial component and by reaching
   the root from the final component. The code ranges over every
   such outcontact mask, and independently checks reachability of
   the displayed witnesses.

Replay with C++17:

    g++ -std=c++17 -O3 check_second_round.cpp -o /tmp/g45_replay
    /tmp/g45_replay

The five strong host totals, the contact counts and the total
deletion count are checked even under optimized compilation; the
proof does not assume graph form, original source containment,
or a designated root index before the calculation.

## Complete fibers, beyond the second window

The same estimator also computes the full fibers. The count bound
\(H\le 2O-g\) comes from the root insertion formula: its stopped
even old-path sum is a sub-sum of the active even old-path sum
\(O-g\). The prefix multiplier \(a\) must divide 9, so enumerate
\(a=1,3,9\), with terminal continuation \((E,O,g)/a\).
At \(a=3\), delete to total 15, with strong factors (15) or (3,5);
at \(a=9\), delete to total 5, with only factor (5). The component
bank computed above still suffices. For each case Busch bounds the
terminal host order by the largest \(n\) with \(a f(n)\le 2O-g\).
Conversely each such terminal witness extends with the needed prefix:
prepend zero, one, or two ordinal directed triangles, respectively.

| q | Prefix multiplier | Actual count support |
|---:|---:|---|
| 79 | 1 | 573, 717, 723, 735, 801 |
| 79 | 3 | 609, 711, 717, 723, 801 |
| 79 | 9 | 711 |
| 163 | 1 | 1027, 1053, 1151, 1205, 1263, 1305, 1329, 1365, 1401, 1467, 1473, 1479, 1491, 1557 |
| 163 | 3 | 999, 1053, 1305, 1365, 1467, 1473, 1479, 1557 |
| 163 | 9 | 1053, 1467 |

These are independently asserted by check_second_round.cpp (same compile
command as above); results.txt includes a literal
terminal witness at each row's count. The divisibility pattern
near \(9q\) is local to that window, not a congruence for the
whole rooted fiber: at \(q=163\), for example, 1027 is not a
multiple of three.

Run **python3 witnesses.py** to check every recorded terminal witness by a separate direct rooted subset DP.
