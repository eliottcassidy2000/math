"""aud_structure_examples.py -- audit S11: explicit gadget matrices of order <= 5 sharing the right eigen-directions
(1,1,1) and (1,0,2) with expanding eigenvalues (verification by direct multiplication), with det and order."""
mats = {}
for k in range(1, 6):
    for line in open(f'aud_mats{k}.txt'):
        a = tuple(map(int, line.split())); mats.setdefault(a, k)
def mul(m, v): return tuple(sum(m[3*i+j]*v[j] for j in range(3)) for i in range(3))
def det(m): return m[0]*(m[4]*m[8]-m[5]*m[7]) - m[1]*(m[3]*m[8]-m[5]*m[6]) + m[2]*(m[3]*m[7]-m[4]*m[6])
for e in [(1, 1, 1), (1, 0, 2), (3, 1, 3)]:
    print('direction', e)
    for m, k in sorted(mats.items(), key=lambda x: x[1]):
        w = mul(m, e)
        lam = w[0] // e[0]
        if all(w[i] == lam * e[i] for i in range(3)) and lam != 0:
            print(f'   order {k}: M = [{m[0]} {m[1]} {m[2]}; {m[3]} {m[4]} {m[5]}; {m[6]} {m[7]} {m[8]}]  M e = {lam} e   det = {det(m)}')
