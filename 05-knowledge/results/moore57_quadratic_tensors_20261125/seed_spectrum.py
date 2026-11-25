"""Exact even trace bound on all nonconstant Fourier components of a seed."""
from pathlib import Path
import numpy as np

def certify(path):
    q = 56
    rows = [list(map(int, line.split())) for line in path.read_text().splitlines()]
    n = len(rows) + 1
    # Column i starts as q times a single-site source minus its fiber constant.
    states = np.zeros((n, n, q), dtype=object)
    for i in range(n):
        states[i,i,:] = -1
        states[i,i,0] += q
    for m in range(1, 101):
        updated = states.copy()
        for i,row in enumerate(rows):
            for j,shift in enumerate(row, i+1):
                updated[:,j,:] += 2*np.roll(states[:,i,:], shift, axis=1)
                updated[:,i,:] += 2*np.roll(states[:,j,:], -shift, axis=1)
        states = updated
        norm = (states*states).sum()
        assert norm % q == 0
        trace = norm // q
        if trace < 225**m:
            print(n, "branches, exponent", 2*m, "trace", trace, "<", 225**m)
            return
    raise AssertionError("insufficient trace bound")

if __name__ == "__main__":
    for n in (18,19,20):
        certify(Path(__file__).with_name(f"seed{n}.txt"))
