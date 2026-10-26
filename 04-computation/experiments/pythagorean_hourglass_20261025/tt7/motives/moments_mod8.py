import numpy as np
exec(open('moments_new.py').read().split('print(f"{len(qs)} primes in')[0])
for cls in (3, 7):
    sel = (qs.astype(int) % 8 == cls)
    print(f"p = {cls} mod 8 ({sel.sum()} primes):")
    for j, nm in enumerate(['m1 (Y02)', 'm2 (Y09)', 'm3 (Y18)', 'm4 (Y05+Y18)']):
        x = m[sel, j]
        print(f"   {nm}: E[x^2] = {np.mean(x**2):.3f}, E[x^4] = {np.mean(x**4):.3f}, max|x| = {np.max(np.abs(x)):.3f}")
