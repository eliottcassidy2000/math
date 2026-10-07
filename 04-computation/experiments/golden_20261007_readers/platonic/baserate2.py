exec(open('baserate.py').read().split("alpha0 = ")[0])
import random
random.seed(7)
N = 200000
for target, name in [({5, 7, 11}, "Galois primes {5,7,11}"), ({12}, "{12}"), ({8, 12}, "{8,12}"), ({7, 11}, "{7,11} (the -17 cycle shape)")]:
    for lo, hi in [(1.5, 5/3), (1.2, 1.9)]:
        c = sum(1 for _ in range(N) if target <= semis(random.uniform(lo, hi)))
        print("P(%s subset of C(alpha)), alpha~U[%g,%g]: %.3f" % (name, lo, hi, c / N))
