"""Independent exact controls for entry policy and known-convergent omissions."""
from pathlib import Path


def need(test, label):
    if not test:
        raise RuntimeError(label)


def v2(n):
    return (n & -n).bit_length()-1


def U(n):
    z = 3*n+1
    return z >> v2(z)


def main():
    root = Path(__file__).resolve().parents[2]
    bank_text = (root/'05-knowledge/results/reset_20260926_swaplift.out').read_text(encoding='utf-8-sig')
    rows = [tuple(map(int, line.split())) for line in bank_text.splitlines()
            if len(line.split()) == 9 and all(s.isdecimal() for s in line.split())]
    need([r[0] for r in rows] == list(range(1, 342, 2)), 'bank domain and order')
    K = max(r[-1] for r in rows)
    need(K == 65 and max(r[1]+1 for r in rows) == 41, 'inherited scope constants')
    need(all(r[-2] != (1 << r[-1])-1 for r in rows), 'all negative-one cylinders excluded')
    print('Inherited bank independently parsed:171 odd q values; Kmax65; maximum horizon41; every residue excludes -1')

    print('Known-convergent completion N_H=2^H*(2^(3^(H-1))+1)/3^H-1')
    print('H exponent_M source_bits exact_U_steps')
    for H in range(2, 11):
        M = 3**(H-1)
        numerator = (1 << M)+1
        need(numerator % 3**H == 0 and numerator % 3**(H+1) != 0, ('exact ternary depth', H))
        t = numerator//3**H
        need(t > 0 and t % 2 == 1, ('tail type', H))
        n = (t << H)-1
        x = n
        for j in range(H):
            need(x == 3**j*(t << (H-j))-1, ('exact phase state', H, j))
            need(v2(x+1) == H-j and x > 1, ('decreasing phase rank', H, j))
            if j < H-1:
                need(v2(3*x+1) == 1 and U(x) > x, ('growth step', H, j))
            else:
                need(3*x+1 == 1 << (M+1), ('terminal giant division', H))
            x = U(x)
        need(x == 1, ('known convergence', H))
        print(H, M, n.bit_length(), H)

    # Only modular powers are formed here; 2^(3^(H-1)) is never allocated.
    for H in range(2, 130):
        M = 3**(H-1)
        modulus = 3**(H+1)
        residue = (pow(2, M, modulus)+1) % modulus
        need(residue % 3**H == 0 and residue != 0, ('modular exact ternary depth', H))
    print('128 modular exact-divisibility controls H2..129; enormous source integers not constructed')

    modulus = 1 << K
    checked = 0
    for fuel in range(65):
        H = fuel+K
        M = 3**(H-1)
        t_mod = ((pow(2, M, modulus)+1)*pow(pow(3, H, modulus), -1, modulus)) % modulus
        need(t_mod % 2 == 1, ('modular tail odd', fuel))
        for j in range(fuel+1):
            state_mod = (pow(3, j, modulus)*pow(2, H-j, modulus)*t_mod-1) % modulus
            need(state_mod == modulus-1 and state_mod % 4 == 3, ('unmatched policy prefix', fuel, j))
            need(all(state_mod % (1 << r[-1]) != r[-2] for r in rows), ('no inherited bank cylinder', fuel, j))
            checked += 1
    print('65 known-convergent source formulas,2145 modular policy states: each fuel d rejects N_(d+65)')
    need(checked == 2145, 'policy state census')

    # Independent arithmetic behind the 47-completion alternative: any
    # positive even coefficient has the same guarded growing block.
    cases = 0
    for a in range(2, 34, 2):
        for k in range(21, 29):
            n = a*8**k-5
            x = n
            for j in range(1, 43):
                x = U(x)
                need(x > n, ('no41-step descent', a, k, j))
            need(n % 4 == 3 and U(n) % 8 == 1, ('exception followed by quarter', a, k))
            need(U(U(n)) == a*9*8**(k-1)-5, ('block continuation', a, k))
            cases += 1
    print('Alternative completed47-family mechanism:128 even-coefficient controls, first42 steps all above each source')
    print('PASS: strict C_d omission from the known-convergent basin is proved by explicit completion, not assumed Collatz')


if __name__ == '__main__':
    main()
