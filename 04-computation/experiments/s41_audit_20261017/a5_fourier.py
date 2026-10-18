#!/usr/bin/env python3
"""S41 audit A5 (independent referee code): THM-4667 statement 10.

1. Own K_d tables (forward push, exact integers, bands [0, P_4)) and their additive Fourier transforms over Z/3^(d+1);
   the claimed recursion Khat_d(xi, b) = sum_(o<=b) (-1)^(o+1) e(-mu_o xi/3^(d+1)) Khat_(d-1)(2^-(o+1) xi mod 3^d, b-o),
   mu_o = 2^-(o+1) - (1 - 4^-j)/3, j = floor(o/2), is evaluated with own code and compared, d = 1..4.
2. F_o(s) = B^-j((3s+1) 2^-a) = 3 2^-(o+1) s + mu_o and F_A = (3s + c_A)/2^A agree as maps mod 3^(d+1).
3. Composition offsets: the unfolded phase of a word (o_1..o_d) is the offset of F_(o_d) o ... o F_(o_1). It is compared
   with (i) the ordinary Syracuse offset of the letters (o_i + 1) and (ii) the factorisation
   B^-(j_d) o U_(a_d + 2 j_(d-1)) o ... o U_(a_2 + 2 j_1) o U_(a_1)   (P^d = C S^(d-1) A).
4. The example P_3(r) = (3r-1)/8 vs S_3(r) = (3r+1)/8 (values 1 and 8 at r = 0 mod 9), and the clock example
   11 -> 17 -> 13 = B(3).
Usage: python a5_fourier.py
"""
import os, sys, json, itertools
sys.dont_write_bytecode = True
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
K0 = np.array([[-1, -1, 2], [0, -3, 3], [-1, -1, 2], [-3, 0, 3], [-1, -1, 2], [0, 0, 0]], dtype=np.int64)  # checked in A4


def F(Q, A):
    a = 1 if A % 2 == 1 else 2
    cA = (3 + 2 ** a - 2 ** A) // 3
    inv = pow(2, -A, Q)
    return lambda s: ((3 * s + cA) * inv) % Q


def tables(D, BR):
    K = np.array([K0[b % 6] for b in range(BR)]).T.copy()        # [r][b]
    out = [K]
    for d in range(1, D + 1):
        Q = 3 ** (d + 1)
        q = Q // 3
        new = np.zeros((Q, BR), dtype=np.int64)
        for o in range(BR):
            f = F(Q, o + 1)
            sg = 1 if (o + 1) % 2 == 0 else -1
            tgt = np.array([f(s) for s in range(q)])
            new[tgt, o:] += sg * K[:, :BR - o]
        K = new
        out.append(K)
    return out


def main():
    res = {}
    D = 4
    BR = 2 * 3 ** (D + 1)
    Ks = tables(D, BR)
    # 1. Fourier recursion
    for d in range(1, D + 1):
        Q = 3 ** (d + 1)
        q = Q // 3
        Fold = np.fft.fft(Ks[d - 1], axis=0)        # Khat_(d-1)(xi, b), xi mod q
        Fnew = np.fft.fft(Ks[d], axis=0)            # Khat_d(xi, b), xi mod Q
        xi = np.arange(Q)
        Pd = 2 * Q
        pred = np.zeros((Q, Pd), dtype=complex)
        for o in range(Pd):
            j = o // 2
            inv = pow(2, -(o + 1), Q)
            t = (1 - pow(4, -j, 3 * Q)) % (3 * Q)
            assert t % 3 == 0
            mu = (inv - t // 3) % Q
            src = (pow(2, -(o + 1), q) * xi) % q
            phase = np.exp(-2j * np.pi * mu * xi / Q) * (-1) ** (o + 1)
            # add to all b >= o: Khat_(d-1)(src, b - o)
            pred[:, o:Pd] += phase[:, None] * Fold[src][:, :Pd - o]
        err = float(np.abs(pred - Fnew[:, :Pd]).max())
        scale = float(np.abs(Fnew[:, :Pd]).max())
        res['fourier_d%d' % d] = dict(max_err=err, scale=scale)
        print('1. d = %d: max |Khat_d - recursion| = %.2e (scale %.1f)' % (d, err, scale))
    # 2. two forms of the branch map
    bad = 0
    for d in range(1, 5):
        Q = 3 ** (d + 1)
        for o in range(0, 2 * Q):
            j, a = o // 2, 1 + o % 2
            A = o + 1
            inv = pow(2, -A, Q)
            mu = (inv - ((1 - pow(4, -j, 3 * Q)) % (3 * Q)) // 3) % Q
            f = F(Q, A)
            for s in range(Q // 3):
                Binv = (((3 * s + 1) * pow(2, -a, Q)) - (pow(4, j, 3 * Q) - 1) // 3 % Q) * pow(4, -j, Q) % Q
                if not (f(s) == (3 * inv * s + mu) % Q == Binv):
                    bad += 1
    res['branch_map_forms_bad'] = bad
    print('2. F_A = 3 2^-A s + mu_o = B^-j((3s+1)/2^a) mod 3^(d+1), d <= 4, all o < 2 3^(d+1): %d disagreements' % bad)
    # 3. composition offsets
    def compose_offset(word, d):
        """offset of F_(A_d) o ... o F_(A_1) evaluated mod 3^(d+1), as a map from Z/3 (lifted) ... -> Z/3^(d+1)."""
        Q = 3 ** (d + 1)
        x = 0
        for A in word:
            x = F(Q, A)(x)
        return x
    def syr_offset(word, d):
        Q = 3 ** (d + 1)
        x = 0
        for A in word:
            x = ((3 * x + 1) * pow(2, -A, Q)) % Q
        return x
    def csa_offset(word, d):
        Q = 3 ** (d + 1)
        aa = [1 if A % 2 == 1 else 2 for A in word]
        jj = [(A - a) // 2 for A, a in zip(word, aa)]
        letters = [aa[0]] + [aa[i] + 2 * jj[i - 1] for i in range(1, len(word))]
        x = 0
        for L in letters:
            x = ((3 * x + 1) * pow(2, -L, Q)) % Q
        jd = jj[-1]
        x = ((x - (pow(4, jd, 3 * Q) - 1) // 3) * pow(4, -jd, Q)) % Q
        return x
    stats = {}
    for d in (1, 2, 3):
        same_syr = same_csa = tot = 0
        for word in itertools.product(range(1, 9), repeat=d):
            tot += 1
            c = compose_offset(word, d)
            same_syr += (c == syr_offset(word, d))
            same_csa += (c == csa_offset(word, d))
        stats[d] = dict(words=tot, equal_to_ordinary_syracuse_offset=same_syr, equal_to_CSA=same_csa)
        print('3. d = %d: %d words (letters 1..8): composed offset = ordinary Syracuse offset in %d; = C S^(d-1) A form in %d'
              % (d, tot, same_syr, same_csa))
    res['offsets'] = stats
    # 4. examples
    Q = 9
    p3 = F(Q, 3)
    s3 = lambda r: ((3 * r + 1) * pow(8, -1, Q)) % Q
    print('4. A = 3 mod 9: P_3(0) = %d, S_3(0) = %d; P_3(r) = (3r-1)/8: %s'
          % (p3(0), s3(0), all(p3(r) == ((3 * r - 1) * pow(8, -1, Q)) % Q for r in range(3))))
    x = 11
    steps = 0
    traj = [x]
    while x != 13:
        x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
        steps += 1
        traj.append(x)
    print('   11 -> 13 in %d shortcut steps: %s; seeds: P(11) = 17 (a = 1, braid 0), P(17) = seed(13) = 3 (a = 2, 13 = B(3))'
          % (steps, traj))
    res['examples'] = dict(P3_0=p3(0), S3_0=s3(0), clock_11_to_13=steps)
    json.dump(res, open(os.path.join(HERE, 'a5_fourier.json'), 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
