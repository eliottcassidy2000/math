"""Audit extras.  Auditor for opus-2026-10-22-S55.
(1) D(Z_N): after the partial Fourier transform in the charge variable, C = T S T^-1 S^-1 is a permutation matrix of a
    matrix in the SL2(Z)-class of A acting on (Z/N)^2 -> ord C = proj ord C = ord_N(A) is a theorem (flat Z_N connections).
(2) Rapidity: A in SL2(R) has eigenvalues phi^(+-2); its image in SO(2,1) (adjoint action on sl2 with the Killing form)
    has eigenvalues phi^(+-4), 1, i.e. it is a boost of rapidity 4 log phi (spinor eigenvalues e^(+-eta/2)).
Usage: python x_extra.py"""
import itertools, math
import numpy as np


def dzn_C(N):
    an = [(a, b) for a in range(N) for b in range(N)]
    S = np.array([[np.exp(-2j * np.pi * (a * d + c * b) / N) / N for (c, d) in an] for (a, b) in an])
    T = np.diag([np.exp(2j * np.pi * a * b / N) for (a, b) in an])
    return T @ S @ T.conj() @ np.linalg.inv(S), an


for N in (3, 4, 5, 6, 8):
    C, an = dzn_C(N)
    # partial Fourier transform in the first label: |x, b> = N^-1/2 sum_a exp(2 pi i a x/N) |a, b>
    F = np.zeros((N * N, N * N), dtype=complex)
    for i, (x, b) in enumerate(an):
        for j, (a, bb) in enumerate(an):
            if b == bb:
                F[i, j] = np.exp(2j * np.pi * a * x / N) / math.sqrt(N)
    P = F @ C @ np.linalg.inv(F)
    isperm = np.allclose(np.abs(P), np.round(np.abs(P)), atol=1e-9) and np.allclose(np.round(np.abs(P)).sum(0), 1) \
        and np.allclose(P[np.abs(P) > 0.5], 1, atol=1e-9)
    # read off the induced map on (Z/N)^2 and test linearity
    img = {an[j]: an[int(np.argmax(np.abs(P[:, j])))] for j in range(N * N)}
    e1, e2 = img[(1, 0)], img[(0, 1)]
    M = np.array([[e1[0], e2[0]], [e1[1], e2[1]]]) % N
    lin = all(img[(x, y)] == (tuple((M @ np.array([x, y])) % N)) for (x, y) in an)
    print(f'N = {N}: F C F^-1 is a 0/1 permutation matrix: {isperm}; induced map linear: {lin}, matrix mod N = {M.tolist()}, '
          f'trace {int(np.trace(M)) % N} (= 3 mod N: {int(np.trace(M)) % N == 3 % N}), det {int(round(np.linalg.det(M))) % N}')

# (2) rapidity
A = np.array([[2.0, 1.0], [1.0, 1.0]])
basis = [np.array([[1.0, 0], [0, -1.0]]), np.array([[0, 1.0], [0, 0]]), np.array([[0, 0], [1.0, 0]])]
Ad = np.zeros((3, 3))
B = np.array([b.flatten() for b in basis]).T
for j, X in enumerate(basis):
    Y = A @ X @ np.linalg.inv(A)
    Ad[:, j] = np.linalg.lstsq(B, Y.flatten(), rcond=None)[0]
phi = (1 + 5 ** 0.5) / 2
ev = sorted(np.linalg.eigvals(Ad).real)
print('Ad(A) eigenvalues:', np.round(ev, 9), '; phi^-4, 1, phi^4 =', np.round([phi ** -4, 1, phi ** 4], 9),
      '; rapidity in SO(2,1) = log(phi^4) =', round(4 * math.log(phi), 9), '; 2 log phi =', round(2 * math.log(phi), 9))
