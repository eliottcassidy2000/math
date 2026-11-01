"""Independent bitset-DP, emitted-fraction, and literal integer checks."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).parent))
from audit_compiler import check, ir_run, raw_run

def load_source(path):
    spec = importlib.util.spec_from_file_location('audited_hamilton_source', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.hamilton

def count_dp(n, bits):
    if not n: return 1
    outgoing = [0]*n
    for j in range(n):
        for i in range(j):
            if bits >> (j*(j-1)//2+i) & 1: outgoing[i] |= 1 << j
            else: outgoing[j] |= 1 << i
    dp = [[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i] = 1
    for mask in range(1,1<<n):
        for last in range(n):
            remaining = outgoing[last] & ~mask
            while remaining:
                bit = remaining & -remaining
                remaining -= bit
                nxt = bit.bit_length()-1
                dp[mask|bit][nxt] += dp[mask][last]
    return sum(dp[-1])

def emitted_run(p, path, n, limit=20_000_000):
    """Read actual fraction text, factor it, retain global priority among eligible rows.

    Other rows are excluded only by their missing unique state-prime factor.
    """
    data = p['data_primes']; state = p['state_primes']; primes = data+state
    index = {prime:i for i,prime in enumerate(primes)}
    def factor(value):
        out = {}
        for prime in primes:
            if prime*prime > value: break
            while value%prime == 0:
                out[index[prime]] = out.get(index[prime],0)+1
                value //= prime
        if value>1:
            check(value in index, 'nonallocated prime')
            out[index[value]] = out.get(index[value],0)+1
        return out
    grouped = {}
    for sequence, line in enumerate(path.read_text().splitlines()):
        num,den = map(int,line.split('/'))
        check(math.gcd(num,den)==1,'unreduced fraction')
        positive, negative = factor(num), factor(den)
        controls = [i for i in negative if i>=len(data)]
        check(len(controls)==1 and negative[controls[0]]==1,'bad denominator control')
        destinations = [i for i in positive if i>=len(data)]
        check(len(destinations)==1 and positive[destinations[0]]==1,'bad destination control')
        grouped.setdefault(controls[0],[]).append((sequence,positive,negative,destinations[0]))
    exps = [0]*len(primes)
    exps[p['input']] = n
    current = len(data); exps[current] = 1
    steps = 0
    while steps < limit:
        for _,positive,negative,destination in grouped.get(current,[]):
            if all(exps[i]>=v for i,v in negative.items()):
                for i,v in negative.items(): exps[i]-=v
                for i,v in positive.items(): exps[i]+=v
                check(exps[current]==0 and exps[destination]==1,'control invariant')
                current=destination
                break
        else:
            check(p['instructions'][current-len(data)][0]=='HALT','unexpected dead state')
            return exps[:len(data)],current-len(data),steps
        steps += 1
    raise RuntimeError('emitted fraction bound')

def main():
    lane=ROOT/'path_counter'
    p=json.loads((lane/'hamilton.ir.json').read_text())
    source=load_source(lane/'hamilton_source.py')
    source_cases=0
    for n in range(6):
        for bits in range(1<<(n*(n-1)//2)):
            code=(1<<n)*(2*bits+1)-1
            expected=count_dp(n,bits)
            check(source(code)==expected,('source',n,bits,expected))
            source_cases+=1
    fraction_cases=[]
    for n in range(4):
        for bits in range(1<<(n*(n-1)//2)):
            code=(1<<n)*(2*bits+1)-1
            expected=count_dp(n,bits)
            ir=ir_run(p,code,bound=20_000_000)
            actual=emitted_run(p,lane/'hamilton.fractions.00.txt',code)
            check(actual==ir,('emitted/IR',n,bits))
            check(actual[0][p['output']]==expected and not actual[0][p['error']],
                  ('fraction/DP',n,bits))
            fraction_cases.append({'order':n,'bits':bits,'code':code,'H':expected,'steps':actual[2]})
    literal=[]
    for code in [0,1,2,3,11]:
        raw=raw_run(p,code)
        check(raw[:3]==ir_run(p,code),('literal/IR',code))
        literal.append({'code':code,'output':raw[0][p['output']],'steps':raw[2],'maxbits':raw[3]})
    # Invalid adjacency words are deliberately rejected rather than silently masked.
    for n in range(5):
        for extra in range(3):
            bits=(1<<(n*(n-1)//2))+extra
            code=(1<<n)*(2*bits+1)-1
            check(source(code)==0,('invalid source',n,bits))
    payload={'status':'PASS','source_exhaustive_orders':'0..5',
             'source_tournament_cases':source_cases,'fraction_cases':fraction_cases,
             'literal_integer_cases':literal,'invalid_source_cases':15,
             'source_sha256':{str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest()
                              for q in [lane/'compiler.py',lane/'hamilton_source.py',lane/'hamilton.ir.json',lane/'hamilton.fractions.00.txt']}}
    print(json.dumps(payload,indent=2))

if __name__=='__main__': main()
