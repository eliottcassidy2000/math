"""Exact translation arithmetic and independent literal source audit.

Run: python -B check.py --write, then python -B check.py and python -O -B check.py.
The C++ auditor is compiled to a temporary directory and accepts literal graphs.
All validity gates remain active under python -O.
"""
from functools import cache
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt, prod
from pathlib import Path
import argparse
import importlib.util
import json
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SHIFT = 30000


def require(test, message):
    if not test:
        raise RuntimeError(message)


def prime_trial(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def sieve(bound):
    p = bytearray(b'\x01') * (bound + 1)
    p[:2] = b'\x00\x00'
    for d in range(2, isqrt(bound) + 1):
        if p[d]:
            p[d*d::d] = b'\x00' * ((bound-d*d)//d+1)
    return p


def factor(n):
    out = []
    d = 2
    while d*d <= n:
        e = 0
        while n % d == 0:
            n //= d
            e += 1
        if e:
            out.append([d, e])
        d += 1
    if n > 1:
        out.append([n, 1])
    return out


@cache
def g(a):
    if a < 2:
        return 1
    q = a // 2
    return g(q) if a & 1 else g(q)+q+(q & 1)


def graph(cut, address=3333, right_address=0, prefix=True):
    core = (18, 8, 3, 5, 14)
    bits = list(map(int, bin(address)[2:])) if address else []
    right = [1-int(x) for x in bin(right_address)[2:][::-1]] if right_address else []
    left = len(bits)
    n = left+len(right)+6
    rows = [0]*n
    for i in range(n-1):
        for j in range(i+1, n-1):
            forward = not(left <= i < j < left+5) or bool(core[i-left] >> (j-left) & 1)
            u, v = (i, j) if forward else (j, i)
            rows[u] |= 1 << v
    rootbits = bits + [(cut >> i) & 1 for i in range(5)] + right
    for v, bit in enumerate(rootbits):
        u, w = (n-1, v) if bit else (v, n-1)
        rows[u] |= 1 << w
    if prefix:
        rows = [a | (((1 << n)-1) << 3) for a in (2, 4, 1)] + [a << 3 for a in rows]
    return rows


def count_h13(a,b,cut):
    total=13*(g(a)+g(b)-1)+104*a*b+47*a+9*(a & 1)+32*b+6*(b & 1)+18
    return total+(3*a-(a & 1)+2 if cut==3 else 0)


def suffix_coefficients(s,width=12):
    """Proved digit recursion: g(2^width*a+s)-g(a)=u*a+v*(a mod2)+w."""
    u=v=w=prefix=0
    for j,bit in enumerate(f'{s:0{width}b}'):
        if bit=='0':
            u+=1 << j
            w+=prefix
            if j:
                w+=prefix & 1
            else:
                v+=1
        prefix=2*prefix+int(bit)
    return u,v,w


def equal_sidecar_atlas():
    """All one-period (12-bit) right-length changes in the two H13 charts."""
    records=[]
    controls=0
    for c in range(1260,1772):
        t=14175-8*c
        s=int(f'{t:012b}'[::-1],2)
        u,v,w=suffix_coefficients(c)
        uu,vv,ww=suffix_coefficients(s)
        for left_cut in (6,3):
            for right_cut in (6,3):
                # H(right)-H(left), with left=(4096*a+c,B), right=(a,4096*B+s).
                k=104*s-47*4095-13*u-(12288 if left_cut==3 else 0)+(3 if right_cut==3 else 0)
                ell=32*4095+13*uu-104*c
                pa=9-13*v-(1 if right_cut==3 else 0)
                pb=13*vv-6
                constant=13*(ww-w)+32*s-47*c-9*(c & 1)+6*(s & 1)
                constant-=(3*c-(c & 1)+2 if left_cut==3 else 0)
                constant+=(2 if right_cut==3 else 0)
                require(k>=0 and ell>=0, ('negative parity-preserving increment',c,s,left_cut,right_cut))
                corners=[(k*a+ell*b+pa*(a & 1)+pb*(b & 1)+constant,a,b) for a in (0,1) for b in (1,2)]
                minimum,min_a,min_b=min(corners)
                for a in range(3):
                    for b in range(1,4):
                        actual=count_h13(a,4096*b+s,right_cut)-count_h13(4096*a+c,b,left_cut)
                        expected=k*a+ell*b+pa*(a & 1)+pb*(b & 1)+constant
                        require(actual==expected, ('affine difference',c,s,a,b,left_cut,right_cut))
                        controls+=1
                records.append(dict(c=c,s=s,left_cut=left_cut,right_cut=right_cut,
                                    coefficients=[k,ell,pa,pb,constant],minimum=minimum,min_a=min_a,min_b=min_b))
    smallest=min(x['minimum'] for x in records)
    require(smallest==49918, 'cross-length minimum')
    return dict(scope='All 512 address carries and all four cut pairs with right lengths differing exactly12, arbitrary a>=0,B>=1',
                families=len(records),minimum_gap=smallest,
                minimizers=[x for x in records if x['minimum']==smallest],
                finite_affine_controls=controls,
                coefficient_certificate_sha256=sha256(json.dumps(records,sort_keys=True).encode()).hexdigest(),
                coefficient_ranges=[[min(x['coefficients'][i] for x in records),max(x['coefficients'][i] for x in records)] for i in range(5)])


def step(st, letter):
    h, e, u, d = st
    if letter == 'O':
        return (h, 0, 2*u, d)
    if letter == 'Z':
        return (h+u-e, d, 2*u-d, d)
    if letter == 'J':
        return (h+u-e, h-e, h+u, h)
    raise ValueError(letter)


def generate():
    primesieve = sieve(30107)
    primes = [n for n in range(2, len(primesieve)) if primesieve[n]]
    ranks = {p: i+1 for i, p in enumerate(primes)}
    first25 = primes[:25]
    table = []
    for p in first25 + [107]:
        q = SHIFT+p
        fs = factor(q)
        require(prod(d**e for d,e in fs) == q, ('factor product',q))
        require(all(prime_trial(d) for d,e in fs), ('nonprime factor',q))
        require(prime_trial(q) == bool(primesieve[q]), ('prime oracles',q))
        table.append(dict(p=p, prime_rank=ranks[p], translated=q, factors=fs,
                          translated_prime_rank=ranks.get(q), residue30=p % 30,
                          offset_from_30030=p-30))
    window = []
    for x in range(108):
        q = SHIFT+x
        require(gcd(x,SHIFT) == gcd(q,SHIFT), ('gcd translation',x))
        require(prime_trial(q) == bool(primesieve[q]), ('window prime oracle',x))
        window.append(dict(x=x, translated=q, original_prime=bool(primesieve[x]),
                           translated_prime=bool(primesieve[q]), factors=factor(q),
                           first25_wheel_survivor=all(q % p for p in first25)))
    pairs=[]
    for u in range(1, isqrt(SHIFT)+1):
        if SHIFT % u:
            continue
        v=SHIFT//u
        if (u+v) % 2:
            continue
        n, m = (v-u)//2, (v+u)//2
        a,b = n*n+1,m*m+1
        pairs.append(dict(u=u,v=v,n=n,m=m,lower=a,upper=b,
                          lower_factors=factor(a),upper_factors=factor(b),
                          both_prime=prime_trial(a) and prime_trial(b)))
    # Independent complete offset-root search has a proved upper bound:
    # u=m-n>=1 implies n=(SHIFT/u-u)/2 <=(SHIFT-1)/2.
    alternative=[]
    for n in range((SHIFT-1)//2+1):
        m=isqrt(n*n+SHIFT)
        if m*m == n*n+SHIFT:
            alternative.append((n,m))
    require(sorted(alternative)==sorted((p['n'],p['m']) for p in pairs), 'square inverse disagreement')
    # Distinct primality oracle: divide only by Eratosthenes-generated primes.
    small_sieve=sieve(7502)
    smalls=[p for p in range(2,7503) if small_sieve[p]]
    for pair in pairs:
        for label in ('lower','upper'):
            value=pair[label]
            sp=value>=2 and not any(value % p == 0 for p in smalls if p*p<=value)
            require(sp == prime_trial(value), ('square prime oracle',value))
    require([p['lower'] for p in pairs if p['both_prime']]==[547601], 'unexpected square-prime intersection')
    require(len(pairs)==15, 'square divisor universe')
    common=prod(p-1 if SHIFT % p == 0 else p-2 for p in first25)
    wheel=prod(first25)
    independent=prod(p-1 for p in first25)
    wheel_receipt=dict(modulus=wheel, translation_period=wheel//gcd(wheel,SHIFT),
                       fixed_prime_channels=[p for p in first25 if SHIFT % p == 0],
                       simultaneous_unit_residues=common, unit_residues=independent,
                       conditional_fraction=str(Fraction(common,independent)))
    graphs=[graph(cut) for cut in (6,3)]
    atlas=equal_sidecar_atlas()
    minimal_graphs=[graph(3,1280,1,False),graph(6,0,8111,False)]
    compensation_graphs=[graph(cut,286,0,False) for cut in (6,3)]
    compensation=[]
    for k in range(1,33):
        D=1 << k
        a=128*D+30
        debt=384*D+92
        low,high=count_h13(a,0,6),count_h13(a,0,3)
        overshoot=3057*D+739
        target=2673*D+647
        require(g(a)==68*D+17 and (low,high)==(6900*D+1649,7284*D+1741), ('compensation formula',k))
        require(high-low==debt and low>overshoot>target, ('compensation baseline obstruction',k))
        compensation.append(dict(k=k,D=D,A=a,debt=debt,actuator_lower=low,actuator_upper=high,
                                 standard_overshoot=overshoot,target=target))
    compiler_path=ROOT/'05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py'
    spec=importlib.util.spec_from_file_location('inherited_compiler',compiler_path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for cut,rows in zip((6,3),graphs):
        original=module.construct((18,8,3,5,14),cut,3333,0)
        direct=[a | (((1 << len(original))-1) << 3) for a in (2,4,1)]+[a << 3 for a in original]
        require(rows == direct, ('independent graph constructor',cut))
    with tempfile.TemporaryDirectory(prefix='chanpak_arithmetic_') as td:
        executable=Path(td)/'audit_state.exe'
        subprocess.run(['g++','-O3','-std=c++17',str(HERE/'audit_state.cpp'),'-o',str(executable)],check=True,capture_output=True,text=True)
        # Acyclic H1 and directed C3 controls precede the actual 21-vertex twins.
        sources=[[0],[6,4,0],[2,4,1]]+graphs+minimal_graphs+compensation_graphs
        data=''.join(f'{len(rows)} {len(rows)-1}\n'+ ' '.join(map(str,rows))+'\n' for rows in sources)
        output=subprocess.run([str(executable)],input=data,capture_output=True,text=True,check=True)
        states=[tuple(map(int,line.split())) for line in output.stdout.splitlines()]
        require(states[:3]==[(1,1,1,1),(1,1,1,1),(3,1,3,1)], ('literal DP controls',states[:3]))
        require(states[3:5]==[(518511,21,1040031,39),(548511,21,1040031,39)], ('literal twin states',states[3:5]))
        require(states[5:7]==[(213805,13,266317,13),(263723,13,266317,13)], ('cross-length minimal states',states[5:7]))
        require(states[7:]==[(15449,7,29789,13),(16309,7,29789,13)], ('compensation minimal states',states[7:]))
        # Mutating an arc pair into a tie must be rejected by the graph guard.
        bad=subprocess.run([str(executable)],input='2 1\n0 0\n',capture_output=True,text=True)
        require(bad.returncode != 0, 'hostile tie was accepted')
    lower,upper=states[3:5]
    base_delta=tuple(b-a for a,b in zip(lower,upper))
    require(base_delta == (30000,0,0,0), 'pure H difference')
    words=['','O','Z','OZZOO','ZZOZO','J','JOZ','JJJ','JOZJ']
    responses=[]
    for word in words:
        a,b=lower,upper
        delta=base_delta
        for letter in word:
            a,b,delta=step(a,letter),step(b,letter),step(delta,letter)
        require(tuple(y-x for x,y in zip(a,b))==delta, ('linear response',word))
        if 'J' not in word:
            require(delta==base_delta, ('portable O/Z response',word))
        responses.append(dict(word=word,lower_state=a,upper_state=b,difference=delta))
    return dict(status='PROVED elementary finite-divisor reduction; FINITE-EXACT arithmetic and literal source audit; Chan-Pak OPEN',
                translation=SHIFT,first25_plus_107=table,window_0_107=window,
                square_pairs=pairs,square_prime_intersection=[547601],
                wheel=wheel_receipt,
                portable_pair=dict(A=3333,B=0,core=[18,8,3,5,14],cuts=[6,3],prefix='C3',
                                   order=21,rows=graphs,states=states[3:5],responses=responses,
                                   inherited_compiler_sha256=sha256(compiler_path.read_bytes()).hexdigest()),
                one_period_equal_sidecar_atlas=atlas,
                one_period_minimum_witness=dict(rows=minimal_graphs,states=states[5:7]),
                compensation_actuator=dict(base_rows=compensation_graphs,base_states=states[7:],dyadic_controls=compensation),
                controls=dict(prime_oracles='Trial division by all integers versus Eratosthenes sieve / sieve-generated primes',
                              square_universe='Every factor pair u*v=30000, u<=v, u=v mod2; independent roots n=0..14999',
                              literal_source_dp='Full subset DP on each 21-vertex graph; positive singleton, TT3, C3 controls',
                              hostile='Invalid two-vertex tie rejected; 30067 survives first25-prime wheel but is composite'))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    report=generate()
    data=json.dumps(report,indent=2,sort_keys=True)+'\n'
    target=HERE/'receipt.json'
    if args.write:
        target.write_text(data,encoding='utf-8')
    else:
        require(target.read_text(encoding='utf-8')==data, 'saved receipt disagreement')
    print(json.dumps(dict(status='PASS',translation=SHIFT,square_prime_intersection=report['square_prime_intersection'],
                         prime_first25_offsets=[p['p'] for p in report['first25_plus_107'][:25] if p['translated_prime_rank']],
                         literal_twins=report['portable_pair']['states']),indent=2))


if __name__=='__main__':
    main()
