"""k-apex transfer as explicit matrices, for middle insertion: f(w) = rho . M_w . sigma0.
Spine sites scanned top to bottom; a site is a k-bit pattern b (bit X set: site -> apex X).
State per apex order pi: (r0 in {L,I}, r_j in {e,L,I} for middle j, rk in {e,n}).  Acceptance needs r0 = L,
each middle run L or (empty and pi_j -> pi_{j+1}).  The reader rho = acceptance pulled back through the bottom part."""
import itertools
import numpy as np
L_, I_, e_, n_ = 'L', 'I', 'e', 'n'

class KApex:
    def __init__(self, k):
        self.k = k
        self.perms = list(itertools.permutations(range(k)))
        self.states = []
        for pi in self.perms:
            for r0 in (L_, I_):
                for mid in itertools.product((e_, L_, I_), repeat=k - 1):
                    for rk in (e_, n_):
                        self.states.append((pi, r0, mid, rk))
        self.idx = {s: i for i, s in enumerate(self.states)}
        self.N = len(self.states)
        self.M = {b: self._mat(b) for b in range(1 << k)}
    def _mat(self, b):
        k = self.k; M = np.zeros((self.N, self.N), dtype=object)
        beats = lambda X: (b >> X) & 1          # site -> apex X
        for s, i in self.idx.items():
            pi, r0, mid, rk = s
            # run 0
            t = (pi, L_ if beats(pi[0]) else I_, mid, rk); M[self.idx[t], i] += 1
            # middle runs
            for j in range(1, k):
                st = mid[j - 1]
                if st == e_ and beats(pi[j - 1]):  # need pi_j -> site
                    continue
                nm = list(mid); nm[j - 1] = L_ if beats(pi[j]) else I_
                M[self.idx[(pi, r0, tuple(nm), rk)], i] += 1
            # last run
            if rk == e_:
                if not beats(pi[k - 1]):
                    M[self.idx[(pi, r0, mid, n_)], i] += 1
            else:
                M[i, i] += 1
        return M
    def init(self):
        v = np.zeros(self.N, dtype=object)
        for pi in self.perms:
            v[self.idx[(pi, L_, (e_,) * (self.k - 1), e_)]] = 1
        return v
    def accept(self, arc):
        """arc(X, Y) = 1 if apex X -> apex Y"""
        a = np.zeros(self.N, dtype=object)
        for s, i in self.idx.items():
            pi, r0, mid, rk = s
            if r0 != L_: continue
            ok = True
            for j in range(1, self.k):
                if mid[j - 1] == L_: continue
                if mid[j - 1] == e_ and arc(pi[j - 1], pi[j]): continue
                ok = False; break
            if ok: a[i] = 1
        return a
    def count(self, word, arc):
        v = self.init()
        for b in word: v = self.M[b].dot(v)
        return int(self.accept(arc).dot(v))

if __name__ == '__main__':
    # check against brute force Hamiltonian path counts
    import random, sys
    sys.path.insert(0, '.')
    from w45_hk import H_of
    rng = random.Random(5); bad = 0
    for trial in range(60):
        k = rng.randint(1, 3); m = rng.randint(1, 6)
        word = [rng.randrange(1 << k) for _ in range(m)]
        arcs = {}
        for X in range(k):
            for Y in range(X + 1, k):
                if rng.random() < .5: arcs[(X, Y)] = 1; arcs[(Y, X)] = 0
                else: arcs[(X, Y)] = 0; arcs[(Y, X)] = 1
        arc = lambda X, Y: arcs[(X, Y)]
        K = KApex(k); c = K.count(word, arc)
        n = m + k; out = [0] * n
        for i in range(m):
            for j in range(i + 1, m): out[i] |= 1 << j
            for X in range(k):
                if word[i] >> X & 1: out[i] |= 1 << (m + X)
                else: out[m + X] |= 1 << i
        for X in range(k):
            for Y in range(k):
                if X != Y and arcs[(X, Y)]: out[m + X] |= 1 << (m + Y)
        if c != H_of(n, out): bad += 1
    print("k-apex matrix transfer vs brute force: failures", bad, "of 60")
