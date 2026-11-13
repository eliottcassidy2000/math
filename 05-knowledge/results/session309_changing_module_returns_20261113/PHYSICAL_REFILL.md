# A prime-count clock rebase with six actual transfers

Use the inherited two-clock frame Π from [the exact target-transition packet](../target_collatz_transitions_20261113/collatz/CLOCK_BOUNDARY.md). If its two nonempty modules at u,v are M,N respectively, then its count is

    H = 60 P_M(3) P_N(2) + H(M) H(N),

where P counts ordered path slots with empty slots allowed. This note uses outgoing-neighbor row masks on vertices labelled from zero.

A first presentation uses modules

| Module | Row masks | H | needed P |
|---|---|---:|---:|
| M | 2,12,9,17,39,15 | 23 | P(3)=1647 |
| N | 2,4,25,35,75,87,15 | 53 | P(2)=844 |

The second takes u=TT1 and v=R→TT6, an ordinal sum with a genuine trailing transitive block, where R has rows

    14,28,104,176,133,211,155,7

and H(R)=499, P_R(2)=7240. The resulting equality is exact:

    60·1647·844+23·53 = 180·7240·64+499 = 83405299.

Trial division through the integer square root certifies that this is prime. The 7 outside labels and all their mutual arcs can be retained, with the two modules replaced as units. No induced inclusion of the old module interiors is claimed. The old right module is strong (checked by reachability); no trailing transitive singleton can be removed while keeping a nonempty right module. The new right module has six such trailing singletons. Its two-slot valuation is 9, up from 2 in the old chart. This is a genuine refill of **transfer sites**, not just a higher valuation.

For k=0,…,6 let u=TT_(1+k) and v=R→TT_(6−k), where TT0 is omitted; both modules remain nonempty. In the quotient the v site precedes c, then u; all arcs at those sites are forward in this order. Flip the arc from the last trailing vertex w of v to c. The vertex w is now the first transitive vertex of u: every other incidence already matches that new placement. Repeat six times, always retaining R and the seven outside vertices. No other arc changes. The exact counts are

    83405299 → 125107699 → 187661299 → 281491699
             → 422237299 → 633355699 → 950033299.

They satisfy H_(k+1)=(3H_k−499)/2. The first, second, fourth and seventh are prime by trial division. This is not ordinary minus Collatz on H: the physical intercept is 499. There is no seventh transferable trailing singleton: R is strong. The construction changes both modules at an exact prime target, then gives six real single-arc moves; it does not provide a second neutral rebase at its terminal prime, an unbounded loop, or coverage of prime targets generally.

The accompanying standard-library verifier checks all four small module responses directly from subset paths and ordered partitions, strongness of N and R, primality of the claimed outputs, and literal named arc comparisons for the six transfers. The inherited frame theorem supplies the large graph counts.
