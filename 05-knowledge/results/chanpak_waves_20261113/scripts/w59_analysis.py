"""Wave 59.  Symmetry and smoothness categorization of all classes of orders <= 8 (from w59_classinfo.txt):
smoothness sigma = largest prime-quotient order in the modular decomposition; |Aut| divides H (free action, LEM-003+ / THM-2477);
symmetry floor: if a prime q divides |Aut(T)| then H(T) >= 1 + q(q^2-1)/12; prime counts p >= 5 are rigid;
records of the number of realizing classes (compare the odd highly composite numbers 1, 3, 9, 15, 45, 105, ...)."""
from collections import defaultdict, Counter
rows = [tuple(map(int, l.split())) for l in open('w59_classinfo.txt')]
def isprime(q): return q > 1 and all(q % d for d in range(2, int(q ** .5) + 1))
def isatom(h): return (isprime(h) and h != 7) or (h % 7 == 0 and isprime(h // 7) and h // 7 >= 5) or h in (63, 343)
byn = defaultdict(Counter)
for n, h, s, a, sc in rows: byn[n][s] += 1
print("classes by order n and smoothness sigma (1 = transitive; sigma = n means prime):")
for n in sorted(byn): print(f"  n={n}: {dict(sorted(byn[n].items()))}")
minsig = {}; autset = defaultdict(set); tot = Counter()
for n, h, s, a, sc in rows:
    minsig[h] = min(minsig.get(h, 99), s); autset[h].add(a); tot[h] += 1
odd = [h for h in range(1, 662, 2) if h in minsig]
print("counts <= 661 by minimal smoothness over their realizers of order <= 8:", dict(sorted(Counter(minsig[h] for h in odd).items())))
for k in (5, 6, 7, 8): print(f"  first counts needing sigma = {k}: {[h for h in odd if minsig[h] == k][:16]}")
atoms = [h for h in odd if isatom(h)]
print(f"atoms <= 661 realized at orders <= 8: {len(atoms)}; by minimal sigma: {dict(sorted(Counter(minsig[h] for h in atoms).items()))}")
print("classes with |Aut| not dividing H:", sum(1 for n, h, s, a, sc in rows if h % a))
viol = 0
for n, h, s, a, sc in rows:
    for q in range(3, a + 1, 2):
        if a % q == 0 and isprime(q) and h < 1 + q * (q * q - 1) // 12: viol += 1
print("symmetry-floor violations (prime q | |Aut| but H < 1 + q(q^2-1)/12):", viol)
print("prime counts >= 5 with a nonrigid realizer:", sorted(h for h, v in autset.items() if isprime(h) and h >= 5 and v != {1}))
print("atoms with a nonrigid realizer:", sorted(h for h, v in autset.items() if isatom(h) and v != {1}))
rich = sorted(((len(v - {1}), h, sorted(v)) for h, v in autset.items()), reverse=True)[:6]
print("counts with the most distinct automorphism-group orders:", rich)
for label, C in (("all orders <= 8", tot), ("order 8 only", Counter(h for n, h, s, a, sc in rows if n == 8))):
    rec = []; best = 0
    for h in sorted(C):
        if C[h] > best: best = C[h]; rec.append((h, C[h]))
    print(f"records of the number of realizing classes ({label}): {rec[:14]}")
