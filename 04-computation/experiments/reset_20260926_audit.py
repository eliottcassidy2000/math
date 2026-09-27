"""Reproduce frozen outputs and independently replay the reset certificates."""
from fractions import Fraction
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / '05-knowledge/results'
LANES = ('inequality', 'clock', 'colours', 'logic', 'swaplift')


def need(ok, label):
    if not ok:
        raise RuntimeError(label)


def normalized(s):
    return s.lstrip('\ufeff').replace('\r\n', '\n').rstrip()+'\n'


def U(n):
    n = 3*n+1
    while n % 2 == 0:
        n //= 2
    return n


def advance(n, steps):
    for _ in range(steps):
        n = U(n)
    return n


def main():
    for lane in LANES:
        script = ROOT / f'04-computation/experiments/reset_20260926_{lane}.py'
        stored = normalized((RESULTS / f'reset_20260926_{lane}.out').read_text(encoding='utf-8-sig'))
        for options in ([], ['-O']):
            done = subprocess.run([sys.executable, '-X', 'utf8', '-B', *options, str(script)],
                                  cwd=ROOT, capture_output=True, text=True, encoding='utf-8')
            need(done.returncode == 0, (lane, options, done.stderr))
            need(normalized(done.stdout) == stored, ('frozen output mismatch', lane, options))
        print(lane+': normal = optimized = frozen output')

    # Direct odd-map iteration without using a precision-controller implementation.
    for h in range(2, 130, 2):
        n = (1 << h)-1
        denominator = 4*(h & -h)
        need(advance(n,h) == (3**h-1)//denominator, ('Mersenne episode',h))
    n = 177
    path = [n]
    for _ in range(4):
        n = U(n)
        path.append(n)
    need(path == [177,133,25,19,29], 'clock hostile')
    need(Fraction(16*13,18*7) == Fraction(104,63), 'reset product increase')
    for k in range(1,129):
        start = (1 << (3*k+2))-5
        n = start
        for j in range(k):
            need(n == 4*9**j*8**(k-j)-5, ('27 family even prefix',k,j))
            n = U(n)
            need(n == 6*9**j*8**(k-j)-7 and n>start, ('27 family odd prefix',k,j))
            n = U(n)
            need(n>start, ('27 family growth',k,j))
        need(n == 4*9**k-5, ('27 family endpoint',k))

    # Independent extraction and replay of the frozen per-q certificate table.
    lines = (RESULTS/'reset_20260926_swaplift.out').read_text(encoding='utf-8-sig').splitlines()
    rows = [tuple(map(int,line.split())) for line in lines
            if len(line.split()) == 9 and all(t.isdecimal() for t in line.split())]
    need(len(rows)==171, 'complete table')
    trie = {}
    huge = 0
    for q,J,A,P,b,clip,R,residue,K in rows:
        x = 3*q
        total = 0
        carry = 0
        for i in range(J):
            need(x >= 2*q, ('first crossing minimal',q,i))
            numerator = 3*x+1
            divisions = (numerator & -numerator).bit_length()-1
            carry = 3*carry+(1 << total)
            if i == J-1:
                need(total==P, ('last-step prefix',q))
            total += divisions
            x = numerator >> divisions
        need(x==b and b<2*q and total==A, ('core replay',q))
        need((1 << A)*b == 3**(J+1)*q+carry, ('affine carry',q))
        need(P<R<=clip<=A, ('threshold domain',q))
        need((1 << (R+1))-3**(J+1)>max(0,3*(1 << (A-R))*b-6*q+1), ('threshold inequality',q))
        need(K==R+1 and (3*residue-6*q+1)%(1 << K)==0, ('cylinder',q))
        offset = 2*q//(1 << K)+1
        for t in (1,2,17,1 << 512):
            n = residue+(1 << K)*(offset+t)
            need(advance(n,J+1)<n, ('huge direct replay',q,t.bit_length()))
            huge += 1
        node = trie
        for i in range(K):
            node = node.setdefault((residue >> i)&1,{})
        node['terminal'] = True

    def mass(node,depth=0):
        if node.get('terminal'):
            return Fraction(1,1 << depth)
        return sum((mass(child,depth+1) for key,child in node.items()
                    if key != 'terminal'),Fraction())
    density = mass(trie)
    need(density==Fraction(6985206796614369409,36893488147419103232), 'independent trie measure')
    print('independent direct replay:64 Mersenne episodes, actual P5 hostile,171 core certificates,',huge,'lifted sources')
    print('independent27-family prefixes:k1..128, all first2k odd iterates exceed their source')
    print('independent low-bit trie density:',density)
    print('PASS: exact finite audit; symbolic proofs and scope separately reviewed')


if __name__ == '__main__':
    main()
