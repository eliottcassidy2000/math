#!/usr/bin/env python3
"""Independent arithmetic-compiler audit and source/IR reproducibility check."""
from pathlib import Path
import importlib.util,json,sys
from hashlib import sha256
from audit_fractions import require
ROOT=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent.parent
spec=importlib.util.spec_from_file_location('reviewed_compiler',ROOT/'fractran/compiler.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def execute(prog,n):
    regs=[0]*prog['register_count'];regs[prog['input']]=n
    pos=0
    for step in range(3000000):
        instruction=prog['instructions'][pos]
        kind=instruction[0]
        if kind=='HALT':return regs[prog['output']],regs[prog['error']]
        if kind=='JUMP':pos=instruction[1]
        elif kind=='INC':
            _,r,dest=instruction;regs[r]+=1;pos=dest
        elif kind=='DECJZ':
            _,r,yes,no=instruction
            if regs[r]==0:pos=no
            else:regs[r]-=1;pos=yes
        else:raise RuntimeError('Unexpected instruction')
    raise RuntimeError('Microprogram step bound')


checks=0
operations={'+':lambda a,b:a+b,'-':lambda a,b:max(0,a-b),'*':lambda a,b:a*b,
            '//':lambda a,b:a//b,'%':lambda a,b:a%b,
            '<':lambda a,b:int(a<b),'>':lambda a,b:int(a>b),'==':lambda a,b:int(a==b),
            '!=':lambda a,b:int(a!=b),'<=':lambda a,b:int(a<=b),'>=':lambda a,b:int(a>=b),
            '**':lambda a,b:a**b}
for op,f in operations.items():
    for b in range(4 if op=='**' else 6):
        program=module.Compiler().compile(f'def count(n):\n    return n {op} {b}\n')
        for a in range(6):
            actual,error=execute(program,a)
            if op in ('//','%') and b==0:require(error==1,'Zero division not trapped')
            else:require(error==0 and actual==f(a,b),f'Primitive mismatch {a}{op}{b}')
            checks+=1
nested='''def count(n):
    x = n
    total = 0
    while x:
        y = x
        while y:
            total += x * y
            y -= 1
        x -= 1
    if n % 2:
        total += (n + 1) * (n + 2) // 2
    else:
        total += n ** 3
    return total
'''
program=module.Compiler().compile(nested)
for n in range(7):
    expected=sum(x*x*(x+1)//2 for x in range(1,n+1))
    expected+=(n+1)*(n+2)//2 if n%2 else n**3
    require(execute(program,n)==(expected,0),'Nested scratch reuse/control mismatch')
    checks+=1
reproduced=[]
for stem in ['a000568','literal_lines']:
    source=(ROOT/f'fractran/{stem}_source.py').read_text()
    actual=module.Compiler().compile(source)
    stored=json.loads((ROOT/f'fractran/{stem}.ir.json').read_text())
    require(json.loads(json.dumps(actual))=={k:stored[k] for k in actual},'Source no longer reproduces IR '+stem)
    reproduced.append(stem)
print(json.dumps(dict(status='PASS',microprogram_checks=checks,reproduced=reproduced,
                     compiler_sha256=sha256((ROOT/'fractran/compiler.py').read_bytes()).hexdigest()),indent=2))
