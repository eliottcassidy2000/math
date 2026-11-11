#!/usr/bin/env python3
# Auditor "hypotheses": HYP-9665 (D) and the note, section 6 ("fan primes").
# For every odd prime p < 400: ord_p(2); if even, t = ord/2 is the least t with p | 2^t + 1 and the fan
# h = 2^(t-1) (x = y = 2^t + 1) is the first BALANCED state of type p (x = y, p | x).
# Question checked: is that fan the FIRST state of type p (p | gcd(x(h), y(h))) at all?
# Scan of all h < 2^26 (own psi recursion).
import numpy as np, sys
from f1core import build_psi


def order2(p):
    o, v = 1, 2 % p
    while v != 1:
        v = (2 * v) % p
        o += 1
    return o


def main():
    NB = 26
    n = 1 << NB
    psi = build_psi(n)
    h = np.arange(n, dtype=np.int64)
    x = 2 * psi + 1
    y = 2 * h + 1
    # balanced states x = y below 2^26
    bal = h[x == y]
    print("# fan_primes.py (auditor hypotheses; own code); states h < 2^%d" % NB)
    print("balanced states (x = y) below 2^%d: %s" % (NB, bal.tolist()))
    assert all((int(v) & (int(v) - 1)) == 0 for v in bal.tolist())
    primes = [p for p in range(3, 400, 2) if all(p % d for d in range(3, int(p ** 0.5) + 1, 2))]
    fan_first = []
    fan_not_first = []
    no_fan = []
    print("p   ord_p(2)  t   fan h=2^(t-1)   first state of type p   fan is first?   fan < p^2 ?")
    for p in primes:
        o = order2(p)
        sel = (x % p == 0) & (y % p == 0)
        idx = np.flatnonzero(sel)
        first = int(idx[0]) if idx.size else None
        if o % 2 == 0:
            t = o // 2
            assert pow(2, t, p) == p - 1
            fan = 1 << (t - 1)
            # the fan is of type p
            if fan < n:
                assert sel[fan]
            if first is None:
                isfirst = "unknown (no state of type p below 2^%d; fan beyond)" % NB
            else:
                isfirst = (first == fan)
            if isfirst is True:
                fan_first.append(p)
            elif isfirst is False:
                fan_not_first.append(p)
            print("%3d  %4d   %3d  %s   %s   %s   %s" % (p, o, t, ("2^%d" % (t - 1)), first, isfirst, fan < p * p))
        else:
            no_fan.append(p)
            # no balanced state of type p: check among balanced states below 2^26 and by the order argument
            assert not any((2 * int(v) + 1) % p == 0 for v in bal.tolist() if v > 0)
            print("%3d  %4d   odd order: no t with p | 2^t + 1; first state of type p: %s" % (p, o, first))
    print()
    print("primes p < 200 with even order whose fan is the first state of type p: %s" % [p for p in fan_first if p < 200])
    print("primes 200 < p < 400 with the same property: %s" % [p for p in fan_first if p > 200])
    print("primes p < 200 with even order whose fan is NOT the first state of type p: %s" % [p for p in fan_not_first if p < 200])
    print("primes p < 200 with odd order (no fan): %s" % [p for p in no_fan if p < 200])


if __name__ == "__main__":
    main()
