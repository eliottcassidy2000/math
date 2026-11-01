"""Independent natural-semantics and literal first-eligible compiler audit.

Run: python3 independent/audit_compiler.py
No compiler-provided execution functions are used.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'fractran'))
from compiler import Compiler, thread_jumps
from baseline_compiler import Compiler as Baseline, fractions

def check(ok, message):
    if not ok:
        raise RuntimeError(message)

class Returned(Exception):
    pass

def source_run(source, n):
    env = {'n': n}
    def ex(x):
        if isinstance(x, ast.Constant): return x.value
        if isinstance(x, ast.Name): return env.get(x.id, 0)
        if isinstance(x, ast.BinOp):
            a, b = ex(x.left), ex(x.right)
            if isinstance(x.op, ast.Add): return a + b
            if isinstance(x.op, ast.Sub): return max(0, a - b)
            if isinstance(x.op, ast.Mult): return a * b
            if isinstance(x.op, ast.FloorDiv): return a // b
            if isinstance(x.op, ast.Mod): return a % b
            if isinstance(x.op, ast.Pow): return a ** b
        if isinstance(x, ast.Compare):
            a, b, op = ex(x.left), ex(x.comparators[0]), type(x.ops[0])
            return int({ast.Eq:a==b, ast.NotEq:a!=b, ast.Lt:a<b,
                        ast.LtE:a<=b, ast.Gt:a>b, ast.GtE:a>=b}[op])
        raise RuntimeError(ast.dump(x))
    def block(stmts):
        for s in stmts:
            if isinstance(s, ast.Assign): env[s.targets[0].id] = ex(s.value)
            elif isinstance(s, ast.AugAssign):
                env[s.target.id] = ex(ast.BinOp(left=s.target,op=s.op,right=s.value))
            elif isinstance(s, ast.If): block(s.body if ex(s.test) else s.orelse)
            elif isinstance(s, ast.While):
                while ex(s.test): block(s.body)
            elif isinstance(s, ast.Return):
                env['OUT'] = ex(s.value)
                raise Returned
            else: raise RuntimeError(ast.dump(s))
    try: block(ast.parse(source).body[0].body)
    except Returned: pass
    except ZeroDivisionError: env['ERROR'] = 1
    return env.get('OUT', 0), env.get('ERROR', 0)

def ir_run(p, n, bound=2_000_000):
    r = [0] * p['register_count']
    r[p['input']] = n
    pc = steps = 0
    while steps < bound:
        op, *args = p['instructions'][pc]
        if op == 'HALT': return r, pc, steps
        if op == 'INC': r[args[0]] += 1; pc = args[1]
        elif op == 'JUMP': pc = args[0]
        elif op == 'DECJZ':
            if r[args[0]]:
                r[args[0]] -= 1; pc = args[1]
            else: pc = args[2]
        else: raise RuntimeError(op)
        steps += 1
    raise RuntimeError('IR step bound')

def raw_run(p, n, bound=2_000_000):
    f = fractions(p)
    primes = f['data_primes'] + f['state_primes']
    value = f['state_primes'][0] * f['data_primes'][p['input']] ** n
    steps = 0
    maxbits = value.bit_length()
    while steps < bound:
        for a, b, _, _ in f['fractions']:
            if value % b == 0:
                value = value // b * a
                break
        else:
            exponents = []
            for prime in primes:
                e = 0
                while value % prime == 0:
                    value //= prime; e += 1
                exponents.append(e)
            check(value == 1, 'unexpected prime factor')
            regs = exponents[:p['register_count']]
            control = exponents[p['register_count']:]
            check(sum(control)==1, 'lost unique control token')
            return regs, control.index(1), steps, maxbits
        steps += 1
        maxbits = max(maxbits, value.bit_length())
    raise RuntimeError('literal step bound')

def main():
    programs = {}
    expressions = ['n+2','n-2','2-n','n*3','n//2','n%2','2**n','n**0',
                   '(n+3)//(n+1)','(n+3)%(n+1)']
    expressions += ['n'+op+'2' for op in ['==','!=','<','<=','>','>=']]
    for expression in expressions:
        programs[expression] = 'def f(n):\n    return '+expression+'\n'
    programs.update({
        'self_add':'def f(n):\n    x=n\n    x+=x\n    return x\n',
        'nested_add':'def f(n):\n    x=n\n    x+=x*(x+1)\n    return x\n',
        'constant_overwrite':'def f(n):\n    n=4\n    n=2\n    n+=3\n    n+=0\n    return n\n',
        'natural_truthiness':'def f(n):\n    if n:\n        n+=1\n    else:\n        n=7\n    return n\n',
        'countdown':'def f(n):\n    x=0\n    while n:\n        x+=n\n        n-=1\n    return x\n',
        'nested_guards':'def f(n):\n    x=0\n    while n>0:\n        if n%2==0:\n            x+=n\n        else:\n            x+=1\n        n-=1\n    return x\n',
        'zero_divisor':'def f(n):\n    return n//0\n',
        'zero_modulus':'def f(n):\n    return n%0\n',
    })
    results = []
    for name, source in programs.items():
        old, new = Baseline().compile(source), Compiler().compile(source)
        for n in range(5):
            expected = source_run(source, n)
            base = ir_run(old, n)
            refined = ir_run(new, n)
            observed = (refined[0][new['output']], refined[0][new['error']])
            check(observed == expected, (name,n,'source/refined',expected,observed))
            check((base[0][old['output']],base[0][old['error']]) == expected,
                  (name,n,'source/baseline'))
            raw = raw_run(new, n)
            check(raw[:3] == refined, (name,n,'literal/refined'))
            results.append({'case':name,'n':n,'output':expected[0],
                            'error':expected[1],'baseline_steps':base[2],
                            'refined_steps':refined[2], 'literal_max_bits':raw[3]})
    # A direct self-transition must not lose its prime guard after threading.
    p = {'instructions':[('DECJZ',2,1,2),('JUMP',0),('HALT',)],
         'register_count':3,'registers':{'ERROR':0,'OUT':1,'n':2},
         'input':2,'output':1,'error':0}
    p = thread_jumps(p)
    check(any(x[0]=='JUMP' for x in p['instructions']), 'missing self-loop bridge')
    for n in range(6): check(raw_run(p,n)[:3]==ir_run(p,n),'bridge mismatch')
    payload = {'status':'PASS','cases':len(results),'literal_first_eligible_cases':len(results)+6,
               'checks':['independent N-valued AST semantics','baseline and refined IR',
                         'full ordered integer fraction scan','unique halt control token',
                         'zero divisors halt with ERROR=1','self-loop guard bridge'],
               'source_sha256':{str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest()
                                for q in [ROOT/'fractran/compiler.py',ROOT/'fractran/baseline_compiler.py']},
               'results':results}
    print(json.dumps(payload,indent=2))

if __name__=='__main__': main()
