"""Independent algebra controls, followed by byte-normalized artifact replay."""
from collections import Counter
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / '05-knowledge/results'


def need(test, label):
    if not test:
        raise RuntimeError(label)


def valuation(n, p=2):
    need(n > 0, 'positive valuation input')
    a = 0
    while n % p == 0:
        a += 1
        n //= p
    return a


def step(n):
    n = 3*n+1
    return n // 2**valuation(n)


def independent_controls():
    types = Counter()
    for source in range(3, 2**18, 8):
        L = valuation(source+5)
        k = (L-1)//3
        r = L-3*k
        b = (source+5)//2**L
        need(r in (1,2,3) and b % 2 == 1, 'three-type decoder')
        n = source
        for j in range(k):
            need(valuation(3*n+1) == 1, 'block first division')
            n = step(n)
            need(n > source and valuation(3*n+1) == 2, 'block second division')
            n = step(n)
            need(n > source, 'block source retention')
        need(n == 2**r*b*9**k-5, 'phase boundary')
        if r == 1:
            need(n % 4 == 1 and step(n) < n, 'type1 local descent')
        elif r == 2:
            need(n % 8 == 7, 'type2 boundary residue')
            ell = 1+valuation(b*9**k-1)
            m = n
            for _ in range(ell):
                need(valuation(3*m+1) == 1, 'type2 exact repeat1')
                m = step(m)
            need(m % 4 == 1 and valuation(3*m+1) >= 2, 'type2 exact stop')
        else:
            numerator = b*9**(k+1)-5
            need(step(step(n)) == numerator//2**valuation(numerator), 'type3 return')
        types[r] += 1

    # Independent integer threshold and source-size tests; do not use producer helpers.
    cylinders = 0
    thresholds = []
    for k in range(160):
        t = 1
        while 2**t*(8**(k+1)-5) <= 9**(k+1)-5:
            t += 1
        b0 = (5*pow(9**(k+1), -1, 2**t)) % 2**t
        need(b0 % 2 == 1, 'odd guard representative')
        thresholds.append(t)
        for lift in (0,1,7,2**80):
            b = b0+2**t*lift
            n = b*8**(k+1)-5
            m = n
            for j in range(2*k+2):
                m = step(m)
                need((m < n) == (j == 2*k+1), 'exact cylinder first descent')
            need(m == (b*9**(k+1)-5)//2**valuation(b*9**(k+1)-5), 'cylinder endpoint')
            cylinders += 1
    need(thresholds[:11] == [1]*5+[2]*6, 'first threshold tiers')

    # Carry scan, using emitted bits rather than the producer's subset machine.
    for n in range(1, 8192, 2):
        bits = [(n >> i)&1 for i in range(n.bit_length()+2)]
        first = next(i for i in range(1,len(bits)) if bits[i] == bits[i-1])
        c = 1
        emissions = []
        for b in bits:
            emissions.append((3*b+c)%2)
            c = (3*b+c)//2
        need(first == next(i for i,e in enumerate(emissions) if e), 'independent wall clock')
        eps = 1 if first%2 == 0 else 5
        h = n//2**(first+1)
        need(step(n) == 6*h+eps, 'independent unread-tail code')

    # Complete ternary period for the target47 depth distribution, independent v3.
    modulus = 3**9
    z = 1
    depths = Counter()
    for h in range(9**4):
        rem = (47*z+7)%modulus
        depth = 4 if rem == 0 else min(4,(valuation(rem,3)-1)//2)
        depths[depth] += 1
        z = z*4%modulus
    need([depths[j] for j in range(5)] == [5832,648,72,8,1], 'independent ternary depth census')
    print('independent three-type decoder:32768 sources; type counts',dict(sorted(types.items())))
    print('independent all-height guard controls:',cylinders,'sources, k0..159, lifts up to2^80')
    print('independent carry-clock controls:4096 odd sources; target47 depth period:6561 exponents')


def normalized(data):
    return data.decode('utf-8-sig').replace('\r\n','\n').strip()


def replays():
    lanes = ('entry','entry_review','families','recursive','wall','recursive_review','incoming')
    for lane in lanes:
        script = ROOT/f'04-computation/experiments/entry_20260927_{lane}.py'
        saved = normalized((RESULTS/f'entry_20260927_{lane}.out').read_bytes())
        for optimized in (False,True):
            args = [sys.executable,'-X','utf8','-B']+(['-O'] if optimized else [])+[str(script)]
            proc = subprocess.run(args,cwd=ROOT,capture_output=True,timeout=240)
            need(proc.returncode == 0,(lane,optimized,proc.stderr.decode('utf-8',errors='replace')))
            need(normalized(proc.stdout) == saved,(lane,optimized,'saved-output mismatch'))
        print('normal/optimized saved-output replay:',lane,'PASS')


if __name__ == '__main__':
    independent_controls()
    replays()
    print('PASS: scoped mathematical controls and all frozen outputs; universal coverage remains OPEN')
