# audit-5129: where the order-10 line counts fall short of the exact r(N) (complete lists, D <= 257, N <= 773),
# and which sources are missing (their least orders)
import pickle
import numpy as np
import supply as S
rows = pickle.load(open("strong_le257.pkl", "rb"))
Root, PI, Aord, Suf, S0, minord = S.build(rows, 257)
DL = {D: S.delta(S0[D]) for D in S0}
r10 = np.load("r10a.npy"); r9 = np.load("line_counts.npz")["r9a"]
for N in range(3, 774, 2):
    src = [(D, (N - D) // 2) for D in range(1, N // 3 + 1, 2) if D in DL and S.is_slope(DL[D], (N - D) // 2)]
    if len(src) != r10[(N - 1) // 2]:
        # least core order: least order of (D, U / 2^i) over i (reading (a): sinks are free)
        late = []
        for (D, U) in src:
            best = 99
            u = U
            while True:
                mo = S.minorder(S0[D], u)
                if mo is not None:
                    # remove appended sinks: the core order is mo minus the number of trailing sinks; reading (a)
                    best = min(best, mo)
                if u % 2:
                    break
                u //= 2
            if best > 10:
                late.append((D, U, best))
        print("N = %d: exact r = %d, r_10 = %d, r_9 = %d ; sources whose least base order (any number of sinks removed) exceeds 10: %s"
              % (N, len(src), r10[(N - 1) // 2], r9[(N - 1) // 2], late))
print("exact r = r_10 for every odd N <= 543:", all(
    sum(1 for D in range(1, N // 3 + 1, 2) if D in DL and S.is_slope(DL[D], (N - D) // 2)) == r10[(N - 1) // 2]
    for N in range(3, 545, 2)))
