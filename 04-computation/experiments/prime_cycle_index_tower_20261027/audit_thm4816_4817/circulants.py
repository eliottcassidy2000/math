"""All circulant tournaments on Z_p (p = 11, 13): arcs strings for audit_reclass, plus multiplier-orbit Burnside."""
import sys
from itertools import product
p = int(sys.argv[1])
half = list(range(1, (p - 1) // 2 + 1))
sets = []
for bits in product([0, 1], repeat=len(half)):
    S = frozenset(d if b else p - d for d, b in zip(half, bits))
    sets.append(S)
if sys.argv[2] == 'arcs':
    for S in sets:
        s = ''.join('1' if (j - i) % p in S else '0' for i in range(p) for j in range(i + 1, p))
        print(f"S={sorted(S)} arcs={s}")
else:
    orbits = set()
    for S in sets:
        orbits.add(min(tuple(sorted((a * x) % p for x in S)) for a in range(1, p)))
    fix = {a: sum(1 for S in sets if frozenset((a * x) % p for x in S) == S) for a in range(1, p)}
    print(f"p={p}: {len(sets)} connection sets, {len(orbits)} multiplier orbits; Burnside ({'+'.join(str(v) for v in fix.values())})/{p-1} = {sum(fix.values())/(p-1)}")
    print("  multipliers fixing some set:", {a: v for a, v in fix.items() if v})
