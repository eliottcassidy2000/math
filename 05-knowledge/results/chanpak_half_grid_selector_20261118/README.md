# Half-grid addresses and joint dyadic/ternary target selection

**Proved:** one physical tournament construction can select any odd
residue simultaneously modulo arbitrary powers of 2 and 3. It uses
the [marked clone core](../chanpak_power_three_return_gates_20261118/README.md)
and a single extra binary prefix. This construction covers only
polylogarithmically many absolute values up to X; Chan–Pak is not
proved. A separate finite-bank obstruction identifies why two
fixed dyadic heights cannot replace the growing core.

## Halves, odds, and the Ramanujan threshold

The half-integer lattice is (1/2)Z. It is closed under addition, not
ordinary multiplication: (1/2)^2=1/4. Doubling sends its integer layer
to even integers and its proper-half layer to odds. If multiplication
is useful, define x⋆y=2xy; with ordinary addition this is exactly
the ring Z in rescaled coordinates. Hamilton path counts H live at
proper half-integers H/2. Ordinal gluing multiplies H, hence uses ⋆,
and the repair H→H+2 is a unit translation on the proper-half layer.
One need not adjoin quarters to use either operation.

Let p_n be the nth prime and R_n the nth Ramanujan prime. The bounds
[Sondow/Laishram](https://arxiv.org/pdf/1105.2249) give p_(2n)<R_n<p_(3n)
for n>1, and the ratio R_n/p_(2n) tends to 1. These are prime
*indices*, not an interpolation of the number 2 into 3. The exact
integer accounting of the defining window π(x)−π(x/2) also lives on
the half-grid: evaluated at integer x, its upper endpoint is integer
and its lower endpoint lands on either of the two layers. Entry
events happen at primes p and exit events at 2p. R_n is a persistent
cutoff, not a single good window. This last distinction matters in
target selection: a construction that can hit each residue may have
a huge first witness, with almost every absolute target still absent
from that construction.

## 1. One exact count with two modes

Fix D=(2,20,9,3,13), apex outmask q=13; rows are out-neighbor bit
masks. Put k ordered clones just before the unchanged marked apex
inside the core, and keep the left transitive chain as a different
block, all beating the old clones and D. For the left chain use
the binary digits of the odd number u as the apex win tag (no
trailing losing vertex). Write F(0)=1, F(2b+1)=F(b),
F(2b)=F(b)+b+[b odd] for 2b>0, and s(u)=u−F(u).
As in the previous packet,

    h_k=12·2^k−1,  W_k=h_k+R_k=36·3^k,  κ_k=2h_k+1,

where R is the two-slot read with the first slot nonempty and
ending in an in-contact. Let R⁻ weight that read by
(-1)^(number of vertices in its first slot). The exact extra read is

    R⁻_k = 6−h_k.

Proof: the return graph has one common endpoint, vertex 1. A
two-slot skeleton has at most three clone runs, and if it has
three then both slots contain a clone run. Each sign-weighted
onto assignment to j runs has expansion, by inclusion–exclusion,
in integer powers of sums of signs over subsets of those runs.
A run in slot one has sign −1, and a run in slot two +1. Thus
all nonzero bases are in {−2,−1,1,2}. The exterior weight
depends solely on the skeleton. At positive k, four reads
determine the complete sum, and the literal k=1,2,3,4 reads are
−17,−41,−89,−185. They give 7−12·2^k. The k=0 read is −5.
No asymptotic fit is used.

The last winning tag at u's end has no spare left vertex.
The parity correction in the odd-cycle packing read is therefore
−R⁻_k (as in the previous marked-apex proof). The exact host count is

    H(k,u) = u W_k + κ_k + h_k(F(u)−1−u) − R⁻_k
           = 36u·3^k + h_k(2−s(u)) − 5.

The host is strong: D with its mark and clones is strong and
the left chain both reaches and is reachable from it (its first
tag wins, and its vertices point to the core). More directly,
D reaches v, v reaches the first left vertex, and that vertex
reaches every remaining one.

## 2. Read the dyadic target, then freeze it

Fix the prescribed odd N and moduli 2^L,3^M, increasing L to
at least 5. Choose an odd suffix v and j∈{1,3} from the
table below. For each v use binary width r=max(3,bit length(v)).
Let B=rev_r(v), the binary reversal with exactly this width.
All B in the table are units modulo 3.

| N mod 32 | (v,j) | N mod 32 | (v,j) |
|---:|---:|---:|---:|
| 1 | (23,3) | 17 | (89,1) |
| 3 | (37,1) | 19 | (7,3) |
| 5 | (1,3) | 21 | (13,1) |
| 7 | (43,1) | 23 | (5,3) |
| 9 | (23,1) | 25 | (47,3) |
| 11 | (19,3) | 27 | (7,1) |
| 13 | (11,1) | 29 | (1,1) |
| 15 | (5,1) | 31 | (35,1) |

The table says N≡36vj+s(v)−7 mod 32. With k≥L the
2^k term vanishes modulo 2^L. Therefore solve

    3^k = [(N−s(v)+7)/4] (9v)^(-1)  (mod 2^(L−2)).

The division is integral; the quotient is j mod 8 by the table.
Powers of 3 are exactly the units congruent to 1 or 3 mod 8
at this modulus: they have order 2^(L−4), by the standard
valuation v₂(3^(2^d)−1)=d+2 for d≥1. The case L=5 is direct.
Choose any k≥max(L,M) in the resulting exponent class.
The dyadic side is now determined by v and s(v) mod 2^L.

## 3. Set the ternary target without moving the dyadic side

Add one high binary block to v, choosing

    u = c·2^(b+r)+v,    b≥L,
    c=1 or c=2^d+1 with fixed even d≥L+1.

The concatenation recurrence gives the exact update

    s(u) = s(v) + s(c) + Bc·2^b.

For a direct derivation, if q=c2^b, concatenating the r-bit
suffix adds q·RC_r(v)+F(v)−1 to F(q), with no boundary bonus
because q is even, where RC is the reverse complement.
Use F(q)=q−s(c)+1 and RC_r(v)+B=2^r−1.
Thus both u and s(u) are unchanged mod 2^L:
s(1)=0, and s(2^d+1)=2^(d−1). Both c choices are units
modulo 3 and the two offsets s(c) are distinct mod 3.

Modulo 3^M the first term in H(k,u) disappears (k≥M).
As h_k is a unit modulo 3, the needed value of s(u) is

    s* = 2−(N+5)h_k^(-1)  (mod 3^M).

Choose the c branch for which s*−s(v)−s(c) is a unit mod 3.
Then solve 2^b = [s*−s(v)−s(c)] (Bc)^(-1) mod 3^M.
The number 2 generates all units modulo 3^M: its order modulo
3 is two and v₃(2²−1)=1 increases by one per cubing (expand
the binomial). Add full periods to make b≥L. When M=0 skip
this step, or choose any prefix preserving the dyadic read.
This proves simultaneous literal selection modulo 2^L3^M.
It also works for the known holes 7 and 21 **as residues**, so
a local certificate on these moduli cannot isolate the holes.

### Why this still cannot choose an arbitrary absolute N

For k≥0 and u≥1 the formula gives
H(k,u)≥u(W_k−h_k)+3h_k−5 ≥24u·3^k.
In the selector above u has a fixed suffix and only one or
two additional nonzero bits in its prefix. Even allowing all
legal b,d,k, there are O((log X)^3) possible reads up to X:
each free exponent is O(log X), and the suffix bank is fixed.
Thus congruence-completeness across two prime powers can coexist
with zero density even among prime targets. An absolute
next step must retain a much richer prefix or suppress the
height cost, rather than increasing the modulus alone.

## 4. Finite two-height obstruction on BOTH prime axes

In the [cyclic two-interface language](../chanpak_exact_two_reservoir_center_20261118/CYCLIC.md)
write H=(A·2^m−h)(B·2^n−k)+2hk, for fixed literal profiles.
For a bank of D pairs, its entire output mod 2^L occupies at
most D(L+1)^2 odd residues: independently replace each exponent
greater than L by L. For L large choose an odd residue r outside
both the bank residue set and its image under multiplication
by 7^(-1). Dirichlet gives infinitely many primes p in this
class, and neither p nor 7p can appear in the bank at any width.
This makes an explicit local obstruction to that entire grammar,
not an obstruction to changing the physical grammar or to
Chan–Pak. The growing ternary core defeats that obstruction
but still faces the absolute bound above.

Run `python3 verify.py` for direct profile and host checks, table
coverage, and joint-modulus tests. The all-order claims use the
run expansion and independent digit lifts, not their finite replay.

The [joint binary controller](ALL_MODULI.md) extends the same
physical core to **every** modulus, retaining this packet's
shorter and quantitatively sparse two-prime subfamily as a
separate exact control.
