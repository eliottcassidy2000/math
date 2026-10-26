# Induced-motive test: a rank-4 motive induced from a quadratic field has trace 0 at the inert half of the primes.
import numpy as np
exec(open('moments_new.py').read().split('print(f"{len(qs)} primes in')[0])
for j, nm in enumerate(['m1', 'm2', 'm3', 'm4']):
    x = np.abs(m[:, j])
    print(f"{nm}: fraction with |x| < 0.05: {np.mean(x < 0.05):.3f}; < 0.2: {np.mean(x < 0.2):.3f}")
