import random
def g(v):
    r = 1
    # iterative version of g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+(v&1)
    bits = bin(v)[2:] if v else ''
    H, A, p = 1, 0, 0
    for ch in bits:
        b = int(ch)
        if b == 1: A = 2 * A + 1; p = 1
        else: H = H + A + p; A = 2 * A; p = 0
    return H
def Hw(w): return g(int(w, 2)) if w else 1
def Aw(w): return int(w, 2) if w else 0
def Cw(w): return sum((1 - int(w[j])) << j for j in range(len(w)))
def pw(w): return int(w[-1]) if w else 0
def ew(w): return 1 if (w and w[0] == '0') else 0
random.seed(1); bad1 = bad2 = 0
for t in range(20000):
    u = ''.join(random.choice('01') for _ in range(random.randint(0, 12)))
    v = ''.join(random.choice('01') for _ in range(random.randint(1, 12)))
    w = ''.join(random.choice('01') for _ in range(random.randint(0, 12)))
    if Hw(u + v) != Hw(u) + Hw(v) - 1 + Aw(u) * Cw(v) + pw(u) * ew(v): bad1 += 1
    L = len(v)
    rhs = Hw(u) + Hw(v) + Hw(w) - 2 + Aw(u) * Cw(v) + ((1 << L) * Aw(u) + Aw(v)) * Cw(w) + pw(u) * ew(v) + pw(v) * ew(w)
    if Hw(u + v + w) != rhs: bad2 += 1
print('law (1) failures:', bad1, ' law (2) failures:', bad2, ' (20000 random triples)')
