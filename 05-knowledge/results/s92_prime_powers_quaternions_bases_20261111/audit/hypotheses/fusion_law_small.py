#!/usr/bin/env python3
# Auditor "hypotheses": HYP-9666 (A), a first numerical datum (exploratory; the hypothesis says nothing is tested).
# Family: ALL labelled pointed tournaments of order n (root = vertex 0), uniform.  For A, B independent of orders
# (n, n): the exact class law of H(A * B) modulo p, computed from the exact histograms of the left states
# (h, H-h, U-h, d-h) of A and of the right states (h', H'-h', V'-h', e'-h') of B, with the fusion corner
#   H(A * B) = u1 v1 + u1 v2 + u2 v1 + u3 v3 + u4 v4.
# Compared with the uniform rank-4 value (class 0: 1 + p^-3 - p^-4, others 1 - p^-4) and with the first-order
# formula of THM-5128 (A): excess of the class 0 over the mean of the others = pi_1 + pi_2 - p^-4 (pi = share of
# the zero state).
import numpy as np, sys
from fusion_states import all_states


def main():
    print("# fusion_law_small.py (auditor hypotheses; own code)")
    st = {}
    for n in range(3, 8):
        m = n * (n - 1) // 2
        res = np.zeros((1 << m, 6), dtype=np.int64)
        all_states(n, res)
        st[n] = res
    for p in (3, 5, 7):
        print("p = %d: uniform rank-4 law: class 0 factor %.5f, other classes %.5f" % (p, 1 + p ** -3 - p ** -4, 1 - p ** -4))
        for n in range(3, 8):
            r = st[n]
            H, h, U, V, d, e = (r[:, i] for i in range(6))
            L = [x % p for x in (h, H - h, U - h, d - h)]
            R = [x % p for x in (h, H - h, V - h, e - h)]
            keyL = ((L[0] * p + L[1]) * p + L[2]) * p + L[3]
            keyR = ((R[0] * p + R[1]) * p + R[2]) * p + R[3]
            DL = np.bincount(keyL, minlength=p ** 4) / len(keyL)
            DR = np.bincount(keyR, minlength=p ** 4) / len(keyR)
            idx = np.arange(p ** 4)
            c4 = idx % p
            c3 = (idx // p) % p
            c2 = (idx // p ** 2) % p
            c1 = (idx // p ** 3) % p
            law = np.zeros(p)
            nzL = np.flatnonzero(DL)
            for i in nzL:
                u1, u2, u3, u4 = c1[i], c2[i], c3[i], c4[i]
                val = (u1 * c1 + u1 * c2 + u2 * c1 + u3 * c3 + u4 * c4) % p
                law += DL[i] * np.bincount(val, weights=DR, minlength=p)
            fac = p * law
            exc = law[0] - law[1:].mean()
            print("   n = %d: class factors %s; excess of class 0 over the mean of the others %.5f; pi_L + pi_R - p^-4 = %.5f" % (
                n, " ".join("%.4f" % f for f in fac), exc, DL[0] + DR[0] - p ** -4))


if __name__ == "__main__":
    main()
