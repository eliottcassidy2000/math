# Exact density-one entry into lawful Collatz words, and the pointwise gap

19 October 2026. Status: PROVED, with classical attribution. The parity-vector bijection is Terras's standard finite encoding theorem, and the density-entry consequences below are elementary consequences of it. No new general Collatz density theorem or resolution is claimed. The contribution to this research line is an explicit entry calculation for the existing bounded-center bank and a precise separation of the remaining quantifiers.

Use the Terras map T(x)=x/2 for even x and (3x+1)/2 for odd x. A word records the parities before successive steps.

## 1. A complete finite-word encoding

For every length N, the map

    x modulo 2^N -> its first N Terras parities

is a bijection onto {0,1}^N.

A short proof is enough here. Suppose a length-N word has lawful residue rho modulo 2^N and branch (3^h x+H)/2^N. Its two lifts modulo 2^(N+1), rho and rho+2^N, have endpoints differing by the odd number 3^h. Their next parities therefore differ. Each old word has exactly one extension of each parity, proving the assertion by induction from N=0.

Primary source: Riho Terras, “A stopping time problem on the positive integers,” Acta Arithmetica 30 (1976), 241-252, https://matwbn.icm.edu.pl/ksiazki/aa/aa30/aa3034.pdf . The argument is reproduced to make the present counting independent of an external implementation.

## 2. Entry into any specified word has density one

Fix any nonempty lawful binary word w of length L, and let C_w be its complete source cylinder. Define

    E_w={x positive: T^j(x) never belongs to C_w for any j>=0}.

### Theorem 1: exact disjoint-window count and exceptional-set bound

For every integer k>=1, the positive integers whose first k disjoint parity blocks of length L all differ from w have natural density exactly

    (1-2^(-L))^k.

Consequently E_w has natural density zero. In fact, writing

    delta_L= -log_2(1-2^(-L))/L >0,

for every real X>=1,

    #{x<=X: x in E_w}
      <= [2/(1-2^(-L))] X^(1-delta_L).

Proof. The first kL parity bits are uniformly distributed over residue classes modulo 2^(kL). There are exactly (2^L-1)^k words avoiding w in the specified disjoint blocks, which gives the density. E_w is contained in this periodic avoidance set for every k. Its upper natural density is at most (1-2^(-L))^k for every k, hence is zero.

For the finite-X bound take k=floor(log_2(X)/L), M=2^(kL), and p=2^(-L). The number of residues of the block-avoidance set modulo M is M(1-p)^k. Every residue occurs at most floor(X/M)+1<=2X/M times among 1,...,floor(X). Thus the count is at most 2X(1-p)^k. As k>=(log_2 X)/L-1, this is at most 2(1-p)^(-1) X^(1-delta_L). The case k=0 satisfies the same bound. QED.

The rate is extremely weak for a long prescribed word. Its purpose is an exact statement about source coverage, not a practical prediction of a waiting time for any fixed source.

For any fixed integer s, the set of sources with fewer than s occurrences among k disjoint blocks has exact density

    sum_(j=0)^(s-1) binom(k,j) (2^(-L))^j (1-2^(-L))^(k-j),

with terms j>k omitted. Letting k grow proves density-one coverage of every fixed finite occurrence quota. It does not prove infinitely many occurrences along any particular positive-integer trajectory. Countable intersections of density-one sets need not have density one, or even be nonempty.

## 3. An explicit all-lift descent cylinder from the bank

The existing expanding word is

    K=110110110110110101011011101010,
    K(x)=(1162261467 x+4912087879)/1073741824.

Append the reset F=10, F(y)=(3y+1)/4. The complete word W=KF has length 32 and

    T^32(x)=(3486784401 x+15810005461)/4294967296,
    x=2962948091 modulo 4294967296.

Every positive source of this cylinder descends. Indeed

    4294967296-3486784401=808182895,
    15810005461 < 20*808182895,

so the affine fixed-point threshold is below 20, whereas the least lawful positive source is 2962948091. The least source has exact endpoint 2405410912. Every larger lawful lift is also below its source after the word.

Theorem 1 therefore supplies a precise partial entry result: a density-one set of positive integers eventually reaches a source at which this specified bank word is lawful and gives a strict local descent. Its exceptional set obeys the displayed power bound with L=32.

This uses the complete guarded word, so there is no assumption that a reset can be chosen at an arbitrary K exit.

## 4. Exact overlapping-window count for this word

W is unbordered: no proper nonempty prefix is a suffix. The companion checker verifies all 31 possible proper border lengths.

Let a_n count binary strings of length n containing no occurrence of W, anywhere. Then

    a_n=2^n                    for 0<=n<32,
    a_n=2 a_(n-1)-a_(n-32)     for n>=32.

Proof. Append either bit to an avoiding word of length n-1. An invalid result must acquire its first occurrence of W at the terminal position. It has the form vW with v avoiding W. Because W has no proper border, no earlier occurrence can overlap the terminal W. Thus these newly invalid results are in bijection with the a_(n-32) possible prefixes v.

Equivalently, inclusion-exclusion gives the exact integer formula

    a_n=sum_(j=0)^floor(n/32) (-1)^j binom(n-31j,j) 2^(n-32j).

Any two occurrences must have disjoint length-32 windows. Choosing j admissible starting positions with that spacing gives binom(n-31j,j); the remaining bits are free. This proves the formula.

The exact natural density of sources that have not entered C_W at any time 0,...,h is

    a_(h+32)/2^(h+32).

For 0<=h<=31, the possible occurrences cannot coexist, and the entry density is exactly (h+1)/2^32. The recurrence and inclusion-exclusion are independently checked against a finite pattern-matching automaton; no orbit simulation substitutes for the residue bijection.

## 5. Why this is still not arbitrary-source descent

Four distinct limits remain important.

1. Density one is not every source. The orbit of 1 alternates 1,2 and never realizes W, for example. More generally a zero-density exceptional set can contain infinitely many positive integers.
2. The certified drop is below the source at entry, T^j(x). It need not be below the original x. Growth before entry is unbounded by the density theorem.
3. Every fixed finite occurrence quota holding for a density-one set does not supply infinite recurrence for a fixed integer. No exchange of these quantifiers is valid.
4. The entry time is not uniformly bounded. For x_q=2^q-1, the first q parities are all odd and

       T^j(x_q)=3^j 2^(q-j)-1 > x_q  for 1<=j<=q.

   W begins 110, so it cannot begin at any time j<=q-3. These positive sources give arbitrarily delayed entry with no descent in their initial all-odd run. The same construction defeats any proposed bounded-length family of first-descent certificates.

Thus the density calculation fills a genuine statistical entry statement while leaving the decisive pointwise and original-height obligations open. It is not a substitute for the bank's bounded-center proof, and the bank proof is not a substitute for entry.

## 6. Connection with the Pythagorean entry theorem

The Berggren positivity language is nonregular, yet the chosen suffix LLLL sends every slope into the positive cone. The present Collatz word is fixed, regular as a single pattern, and eventually occurs for a density-one set of sources, yet cannot be imposed on a particular trajectory. The contrast is exact: recognition complexity, control of a branch, measure coverage, and pointwise descent are separate mathematical properties.

This comparison supplies no conjugacy between the two systems and no equality of their arithmetic invariants. It explains which hypotheses a transferable repair or entry argument must explicitly retain.
