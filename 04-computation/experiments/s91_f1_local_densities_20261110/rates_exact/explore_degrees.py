#!/usr/bin/env python3
"""exploration: degrees of the Q-irreducible factors of the norm characteristic polynomial of T_k."""
import sys
import time

from flint import fmpz_mat

from explore_charpoly import big_matrix


def main():
    args = sys.argv[1:]
    for a in args:
        m, k = (int(t) for t in a.split(":")) if ":" in a else (int(a), 1)
        t0 = time.time()
        B = big_matrix(m, k)
        cp = fmpz_mat(B).charpoly()
        c, facs = cp.factor()
        print(f"m = {m} k = {k}: size {len(B)}; factor degrees (mult): {[(f.degree(), mu) for f, mu in facs]}  [{time.time() - t0:.1f}s]", flush=True)
        for f, mu in facs:
            if f.degree() <= 8:
                print("     ", f)


if __name__ == "__main__":
    main()
