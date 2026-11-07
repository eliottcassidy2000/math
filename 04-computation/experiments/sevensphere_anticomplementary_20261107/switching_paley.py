#!/usr/bin/env python3
# HYP-9500 (mac-mini-2026-11-07-sevensphere): bordered Paley tournament QR_q + source is skew-Hadamard but not strong;
# switching it at one vertex of QR_q keeps it skew-Hadamard (so equicentric, THM-4911(4)) and makes it strong.
# Output: switching_paley.out
import numpy as np
def paley(q):
    QR = {(i * i) % q for i in range(1, q)}
    return np.array([[1 if (j - i) % q in QR else 0 for j in range(q)] for i in range(q)])
def strong(A):
    n = len(A); R = (A | np.eye(n, dtype=int)).astype(bool)
    for k in range(n): R = R | (R[:, [k]] & R[[k], :])
    return bool(R.all())
for q in (7, 11, 19, 23, 31):
    P = paley(q); n = q + 1
    A = np.zeros((n, n), dtype=int); A[:q, :q] = P; A[q, :q] = 1
    M = A - A.T; H = np.eye(n, dtype=int) + M
    D = np.eye(n, dtype=int); D[0, 0] = -1; M2 = D @ M @ D; A2 = (M2 == 1).astype(int); H2 = np.eye(n, dtype=int) + M2
    u = np.ones(n, dtype=int); y = M2 @ (M2 @ u)
    print(f"q={q} (n={n}): QR_q+source skew-Hadamard {bool((H.T @ H == n*np.eye(n, dtype=int)).all())}, strong {strong(A)}; "
          f"switched at a QR_q vertex: skew-Hadamard {bool((H2.T @ H2 == n*np.eye(n, dtype=int)).all())}, strong {strong(A2)}, "
          f"equicentric {bool((y == y[0]).all())}, lambda {int(-y[0])}")
