"""(1) Paley skew operator as sqrt(-q) on the root lattice A_{q-1}: S^2 = J - qI, and (I+S)/2 integral on the sum-zero lattice.
   q=3: S acts on A2 = Z[omega] as 1+2omega.  (2) Dyadic four-bin averaging identity.  (3) paid-ROOT corrected edge identity 2^t Y = 3X + 2^a."""
import numpy as np, random, itertools
from fractions import Fraction
print("(1) Paley skew operator on the sum-zero lattice")
for q in [3, 7, 11, 19, 23]:
    qr = {(x * x) % q for x in range(1, q)}
    S = np.array([[0 if i == j else (1 if ((j - i) % q) in qr else -1) for j in range(q)] for i in range(q)], dtype=np.int64)
    J = np.ones((q, q), dtype=np.int64); I = np.eye(q, dtype=np.int64)
    ok1 = np.array_equal(S @ S, J - q * I)
    # basis of sum-zero lattice: e_i - e_0 (i=1..q-1); check (I+S)/2 maps it into itself (integer coordinates)
    B = np.zeros((q, q - 1), dtype=np.int64)
    for i in range(1, q): B[i, i - 1] = 1; B[0, i - 1] = -1
    M = (I + S) @ B   # columns in Z^q; need every column even
    ok2 = np.all(M % 2 == 0)
    # the image columns divided by 2 lie in sum-zero lattice automatically; coordinates in basis B: x_i = v_i for i>=1
    half = M // 2; ok3 = np.all(half.sum(axis=0) == 0)
    # minimal polynomial of (I+S)/2 restricted: t^2 - t + (q+1)/4
    T = half; coords = T[1:, :]   # coordinates in basis e_i - e_0
    C2 = coords @ coords - coords + ((q + 1) // 4) * np.eye(q - 1, dtype=np.int64)
    ok4 = np.all(C2 == 0)
    print(f"  q={q}: S^2=J-qI {ok1}; (I+S)/2 integral on sum-zero lattice {ok2 and ok3}; theta^2 - theta + (q+1)/4 = 0 on it: {ok4}")
# q=3: A2 lattice identification. Z^3 cap 1^perp with basis u=e0-e1, v=e1-e2; omega acts as cyclic shift (0->1->2) restricted? check S = 1+2*omega
q = 3; S = np.array([[0, 1, -1], [-1, 0, 1], [1, -1, 0]])  # 0->1->2->0 : S_ij=1 if i->j
P = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])  # cyclic shift permutation matrix e_i -> e_{i+1}
print("  q=3: S = P - P^T:", np.array_equal(S, P - P.T), "; on sum-zero plane P acts as omega (P^2+P+I = J = 0 there); 1+2P == S mod J:", np.array_equal((np.eye(3, dtype=int) + 2 * P - S) % 1, np.zeros((3, 3))), (np.eye(3, dtype=int) + 2 * P - S))
print("     (I + 2P) - S = J - ... check: I+2P-S =", (np.eye(3, dtype=int) + 2 * P - S).tolist(), " (= J - P^T - P + ... ) ; restricted to 1^perp: I+2P-S =", ((np.eye(3, dtype=int) + 2 * P - S) @ np.array([[1, 0], [-1, 1], [0, -1]])).tolist())
print("(2) four-bin averaging: E_middle averaged over the 3 pairings = (E_parent + 2 E_children)/3")
rng = random.Random(3)
for trial in range(3):
    z = [(rng.randint(-5, 5), rng.randint(-5, 5)) for _ in range(4)]
    Q = lambda a, b: a * a - a * b + b * b
    add = lambda u, v: (u[0] + v[0], u[1] + v[1])
    E = Q(*add(add(z[0], z[1]), add(z[2], z[3]))); Sx = sum(Q(*w) for w in z)
    pairings = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
    M = [Q(*add(z[a], z[b])) + Q(*add(z[c], z[d])) for (a, b), (c, d) in pairings]
    print(f"  charges {z}: pairings {M}, mean {Fraction(sum(M),3)}, (E+2S)/3 = {Fraction(E + 2*Sx, 3)}, equal: {sum(M)*1 == E + 2*Sx}; mean integral: {(E+2*Sx) % 3 == 0}")
print("(3) paid-ROOT R1 identity 2^t Y = 3X + 2^a, X=3s+1, Y=3y+1, a=2 if s=1 mod 8 else 1, (3s+1)/2^a = B^j(y), B(x)=4x+1, t=a+2j")
bad = 0; cnt = 0; maxj = 0
for s in range(1, 200001, 2):
    if s % 8 == 5: continue
    a = 2 if s % 8 == 1 else 1
    if (3 * s + 1) % (2 ** a): continue
    z = (3 * s + 1) // 2 ** a; j = 0
    while z % 4 == 1 and z > 1: z = (z - 1) // 4; j += 1   # strip the B-chain down to its root y (y not 1 mod 4, or y=1)
    y = z; t = a + 2 * j; X = 3 * s + 1; Y = 3 * y + 1; cnt += 1; maxj = max(maxj, j)
    if 2 ** t * Y != 3 * X + 2 ** a: bad += 1
print(f"  seeds tested: {cnt}, violations: {bad}, max j: {maxj}; with the old constant 2^a+3 the identity fails at s=3: 2^3*4 = {2**3*4} vs 3*10+2 = {3*10+2}, 3*10+5 = {3*10+5}")
