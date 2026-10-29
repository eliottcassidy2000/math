"""Universal linear name-preserving obstruction via a native norm polynomial."""
from pathlib import Path
import hashlib
import json
from sympy import Poly, symbols

from periodic_frobenius_lift import mm, mpow, ident
from ternary_obstruction_verify import kernel

HERE = Path(__file__).resolve().parent
OTHER = HERE.parent/'overnight_rank_return_20261029'
GATES = 0


def need(ok,why):
    global GATES
    GATES += 1
    if not ok:
        raise ValueError(why)


def run():
    src = json.loads((OTHER/'joint_phase.json').read_text())
    terminal = json.loads((OTHER/'terminal_selector.json').read_text())
    periodic = json.loads((OTHER/'periodic_defect.json').read_text())
    T, S, names = src['clock'], src['frobenius'], src['native_actions']
    at = terminal['permutations'].index([1,2,0])
    R, I = names[at], ident(18)
    need(mpow(R,3) == I and R != I, 'actual native name three-cycle')
    need(mm(R,T) == mm(T,R), 'native names centralize the bulk clock')
    R2 = mm(R,R)
    N = [[(I[i][j]+R[i][j]+R2[i][j]) % 3 for j in range(18)] for i in range(18)]
    f = [2,0,0,2,2,1,1,1,1,0,2,1]
    def evaluate(cs,A):
        out = [[0]*18 for _ in range(18)]
        for c in reversed(cs):
            out = mm(out,A)
            out = [[(x+c*I[i][j]) % 3 for j,x in enumerate(row)] for i,row in enumerate(out)]
        return out
    need(evaluate(f,T) == N, 'literal native norm is the saved clock polynomial')
    Ta = mpow(T,487)
    moved = evaluate(f,Ta)
    D = [[(a-b) % 3 for a,b in zip(row,nrow)] for row,nrow in zip(moved,N)]
    need(kernel(N,3)[1] == 6 and mm(N,N) == [[0]*18 for _ in range(18)], 'native norm rank six and square zero')
    need(kernel(D,3)[1] == 5, 'nonzero universal obstruction of rank five')
    witness = next((i,j,D[i][j]) for i in range(18) for j in range(18) if D[i][j])
    need(mm(S,T) == mm(Ta,S) and mpow(S,5) == I, 'inherited positive clock lift remains valid')
    need(mm(S,N) != mm(N,S), 'positive lift necessarily moves native norm')
    # Independent polynomial quotient calculation identifies the affected factor.
    X = symbols('X')
    g = Poly(sum(c*X**i for i,c in enumerate(periodic['field_polynomial_low_first'])),X,modulus=3)
    h = Poly(X+1,X,modulus=3)*g
    q = h*h
    fp = Poly(sum(c*X**i for i,c in enumerate(f)),X,modulus=3)
    xa = Poly(X,X,modulus=3)
    z, base, e = Poly(1,X,modulus=3), xa, 487
    while e:
        if e & 1:
            z = (z*base).rem(q)
        base = (base*base).rem(q)
        e //= 2
    changed = Poly(0,X,modulus=3)
    for c in reversed(f):
        changed = (changed*z+Poly(c,X,modulus=3)).rem(q)
    difference = (changed-fp).rem(q)
    need(fp.rem(h).is_zero, 'norm vanishes on both semisimple field quotients')
    need(difference.gcd(q).monic() == (Poly(X+1,X,modulus=3)**2*g).monic(), 'only the five-dimensional field-jet layer obstructs')
    need(evaluate([int(difference.nth(i)) % 3 for i in range(12)],T) == D, 'polynomial and literal matrix obstruction agree')
    # No map made from commuting native generators can supply this exponent.
    need(all(mm(M,T) == mm(T,M) for M in names), 'all six native names commute with bulk')
    need(Ta != T, 'noncentral clock-normalizing target')
    return {'status':'PROVED exclusion of EVERY invertible F3-linear fixed-name lift; no nonlinear residue claim',
            'source_sha256_canonical':hashlib.sha256(json.dumps(src,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
            'native_three_cycle_index':at,'norm_polynomial_low_first':f,
            'native_norm':N,'changed_norm':moved,'obstruction':D,'obstruction_rank':5,'witness':witness,
            'difference_polynomial_low_first':[int(difference.nth(i)) % 3 for i in range(12)],
            'scope':'S*T=T^487*S and S*R=R*S with S invertible and F3-linear are incompatible',
            'gates':GATES}


if __name__ == '__main__':
    out = run()
    (HERE/'joint_frobenius_obstruction.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('gates','obstruction_rank','witness','scope')}))
