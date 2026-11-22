"""Check the one-hole module signatures and the joint-residue example at E61."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'chanpak_e61_fibonacci_target_20261122'))
from seed_lifts import seeds, signature
from context_primes import profile


def main():
    differences = set()
    physical_example = None
    for low, high, arc in seeds():
        for v in range(len(low)):
            if v in arc:
                continue
            a = signature(low,(v,))
            b = signature(high,(v,))
            assert all(r[0] <= 3 for r in a.keys() | b.keys())
            differences.add(tuple(b[r,] - a[r,] for r in (1,2,3)))
            if low == (2,32,11,19,7,29) and v == 0:
                physical_example = a, b
    assert differences == {(2,0,0),(2,1,0),(2,3,0),(2,6,0),(2,-3,0),(2,-5,-2)}
    assert physical_example == ({(1,):27,(2,):12},{(1,):29,(2,):15})
    H, Q = 1, 2
    for rows in ((6,12,24,65,139,31,55,111),(6,136,42,49,135,83,31,109)):
        runs = profile(rows)
        H *= runs[1]
        Q *= 2*runs[1] + runs[2]
    E = 2**60*(2**61-1)
    L, U = 3*H+12*Q, 15*Q-H
    assert (L,U) == (532354095,665056667)
    assert (L-(E-1))%1365 == (U-(E+1))%1365 == 0
    assert (E-9)%19 != 0
    print('Four seeds; six distinct one-hole differences. Physical counts:',L,U)


if __name__=='__main__':
    main()
