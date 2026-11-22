# Exact arc flips straddling even perfect numbers

**PROVED, explicit witnesses.** Put \(E_p=2^{p-1}(2^p-1)\). There is a single-arc flip between strong tournaments with path counts \(E_p-1\) and \(E_p+1\) at \(p=17,19\). A second physical edge starts at \(E_{19}+1\) and ends at \(E_{19}+3\). The two edges at \(p=19\) share a scalar count and use different tournament orders.

These perfect numbers were specified before the address search. For the indicated exponents, \(2^p-1\) is prime (trial division suffices), whence the divisor sum of \(E_p\) is \((2^p-1)2^p=2E_p\).

## Witnesses

Each row specifies the lower marked deletion core and the single arc to reverse for its upper counterpart. Vertex 0 is the low bit. The apex beats exactly the vertices in the tag mask.

| Chart | Lower core outgoing masks | Reversed arc | Tag mask | Lower and upper \(k\) | \(h,c\) |
|---|---|---|---:|---|---|
| A | (4,5,48,7,11,27) | 2–3 | 1 | 37, 39 | 23, 24 |
| B | (34,80,75,19,69,30,41) | 2–6 | 1 | 291, 293 | 143, 173 |

Here \(h\) is the deletion's Hamiltonian-path count, \(k\) is the count after attaching its marked apex, and \(c\) is the two-slot mass of paths starting tagged. In chart A the even/odd masses are \((7,17)\to(8,16)\); in chart B they are \((74,99)\to(75,98)\). The ordinary counts and unsigned start sums agree on the two sides.

Append a transitive suffix after the core, with its apex-contact word given by reversing and complementing the standard binary digits of the specified even \(w\). For even \(w\) the [one-sided arm formula](README.md#one-arm-target-inversion-and-two-exact-landings) gives the lower and upper counts exactly as
\[
 k+c w+h p(w),\qquad k+2+c w+h p(w),
\]
where \(p(0)=0\), \(p(2z+1)=p(z)\) and \(p(2z)=p(z)+z+[z\text{ odd}]\). The attached marked cores are already strong, and the last vertex of the suffix beats the apex, so both extended hosts are strong.

| Lower target | Upper target | Chart | \(w\) | \(p(w)\) | Order |
|---:|---:|---|---:|---:|---:|
| 8,589,869,055 | 8,589,869,057 | A | 233563022 | 129754630 | 35 |
| 137,438,691,327 | 137,438,691,329 | B | 553528478 | 291456394 | 38 |
| 137,438,691,329 | 137,438,691,331 | B | 474437876 | 387139430 | 37 |

Run this in the directory containing the note. It checks the small core counts by subset dynamic programming, reconstructs the full tournaments, and checks the exact change of arc and the two counts.

```python
from replay import count_from_profile, flipped, prime, profile, star_rows, strong

charts = {
    "A": ((4,5,48,7,11,27), (2,3), 1),
    "B": ((34,80,75,19,69,30,41), (2,6), 1),
}
cases = [(17,-1,"A",233563022), (19,-1,"B",553528478),
         (19,1,"B",474437876)]
for exponent, offset, chart, w in cases:
    rows, edge, tag = charts[chart]
    mersenne = 2**exponent - 1
    assert prime(mersenne)
    target = 2**(exponent-1)*mersenne + offset
    other = flipped(rows, *edge)
    graphs = [star_rows(d, tag, 0, w) for d in (rows, other)]
    assert all(strong(t) for t in graphs)
    assert sum((a^b).bit_count() for a,b in zip(*graphs)) == 2
    counts = [count_from_profile(profile(d,tag),0,w) for d in (rows,other)]
    assert counts == [target,target+2]
    print(exponent, len(graphs[0]), counts)
```

## The tested boundary at exponent 31

Neither \(E_{31}-1=2,305,843,008,139,952,127\) nor \(E_{31}+1=2,305,843,008,139,952,129\) was a lower endpoint in the tested bank of 104 distinct one-sided scalar charts. Inserting a marked tournament before the growing flank subsequently gave three straddling flips at this perfect number, including a pair on 60 vertices. [Exact exponent-31 construction](../chanpak_e31_physical_20261122/README.md)

The finite bank was generated from every stored tournament class in the [order-5, 6 and 7 row files](../../../session212/cycle-levels), each reversal of a single arc, and every nonempty proper tag mask. Retain pairs having equal \(h,c\) and base apex counts differing by two. For each parity \(b\in\{0,1\}\) retain it if the two intercepts \(k-b(e-o)\) also differ by two; deduplicate the triples (lower intercept, \(h,c\)) with their parity. For each chart apply the [split inverse](README.md#one-arm-target-inversion-and-two-exact-landings) with the prescribed parity of the low word. This tests all address lengths for that chart. It leaves other cores, tag choices outside this rule, two-sided settings, and the general flip graph available.
