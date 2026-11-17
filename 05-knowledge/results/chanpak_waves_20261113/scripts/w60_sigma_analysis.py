"""Wave 60: exact smoothness of atoms from the minimal smoothness through order 9 (w60_sigma9 > w60_minsig9.txt).
An atom a has only strong realizers (up to padding), of order n with 5^((n-1)/3) <= a (Busch's bound); so for a < 125
the minimum over orders <= 9 is exact."""
from collections import Counter
d = {int(a): int(b) for a, b in (l.split() for l in open('w60_minsig9.txt'))}
def isprime(q): return q > 1 and all(q % k for k in range(2, int(q ** .5) + 1))
def isatom(h): return (isprime(h) and h != 7) or (h % 7 == 0 and isprime(h // 7) and h // 7 >= 5) or h in (63, 343)
atoms = [h for h in range(3, 125, 2) if isatom(h)]
print("exact smoothness of the atoms below 125:", {h: d[h] for h in atoms})
for k in (3, 5, 6, 7): print(f"  smoothness {k}: {[h for h in atoms if d[h] == k]}")
c7 = [h for h in range(125, 662, 2) if isatom(h) and d.get(h) == 7]
print("atoms in [125, 661] whose realizers through order 9 all need a prime quotient of order >= 7:", c7)
for h in c7[:3]:
    print(f"  {h}: strong realizers have order <= {max(n for n in range(3, 40) if 5 ** ((n - 1) / 3) <= h)} (Busch)")
