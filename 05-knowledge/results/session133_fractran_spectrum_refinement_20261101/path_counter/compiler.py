"""Finite counter-machine compiler for a deliberately small arithmetic WHILE language.

Inputs/registers are nonnegative integers. Subtraction is truncated. Supported:
assignments, +=/-=, if/else, while, return, + - * // % **, comparisons.
There are no arrays, calls, recursion, dynamic instructions or hidden arithmetic.
"""
from __future__ import annotations
import ast
import json
import math
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


class Compiler:
    def __init__(self):
        self.code = []
        self.labels = {}
        self.names = {}
        self.free = []
        self.nreg = 0
        self.serial = 0
        self.error = self.reg('ERROR')
        self.output = self.reg('OUT')
        self.input = self.reg('n')
        self.error_label = self.label('ERROR')
        self.halt_label = self.label('HALT')

    def reg(self, name):
        if name not in self.names:
            self.names[name] = self.nreg
            self.nreg += 1
        return self.names[name]

    def temp(self):
        if self.free:
            return self.free.pop()
        result = self.nreg
        self.nreg += 1
        return result

    def release(self, *regs):
        for reg in regs:
            require(reg not in self.free, 'double free')
            require(reg not in self.names.values(), 'freeing named register')
            self.free.append(reg)

    def label(self, prefix='L'):
        self.serial += 1
        return f'{prefix}_{self.serial}'

    def mark(self, label):
        require(label not in self.labels, 'duplicate label')
        self.labels[label] = len(self.code)

    def emit(self, op, *args):
        self.code.append((op, *args))

    def inc(self, reg):
        self.emit('INC', reg)

    def goto(self, label):
        self.emit('JUMP', label)

    def decjz(self, reg, nonzero, zero):
        self.emit('DECJZ', reg, nonzero, zero)

    def clear(self, reg):
        loop, again, done = self.label(), self.label(), self.label()
        self.mark(loop)
        self.decjz(reg, again, done)
        self.mark(again)
        self.goto(loop)
        self.mark(done)

    def set(self, reg, value):
        self.clear(reg)
        for _ in range(value):
            self.inc(reg)

    def move_add(self, src, dst):
        require(src != dst, 'move alias')
        loop, body, done = self.label(), self.label(), self.label()
        self.mark(loop)
        self.decjz(src, body, done)
        self.mark(body)
        self.inc(dst)
        self.goto(loop)
        self.mark(done)

    def copy(self, src, dst):
        require(src != dst, 'copy alias')
        t = self.temp()
        self.clear(dst)
        self.clear(t)
        loop, body, done = self.label(), self.label(), self.label()
        self.mark(loop)
        self.decjz(src, body, done)
        self.mark(body)
        self.inc(dst)
        self.inc(t)
        self.goto(loop)
        self.mark(done)
        self.move_add(t, src)
        self.release(t)

    def multiply(self, a, b, out):
        require(len({a,b,out}) == 3, 'multiply alias')
        counter, t = self.temp(), self.temp()
        self.clear(out)
        self.copy(a, counter)
        loop, body, done = self.label(), self.label(), self.label()
        self.mark(loop)
        self.decjz(counter, body, done)
        self.mark(body)
        self.copy(b, t)
        self.move_add(t, out)
        self.goto(loop)
        self.mark(done)
        self.release(counter, t)

    def compare_lt(self, a, b, out):
        """Preserve a,b; write exactly int(a<b) into out."""
        x,y = self.temp(), self.temp()
        self.copy(a,x)
        self.copy(b,y)
        self.clear(out)
        loop, y_yes, x_yes, lt, done = [self.label() for _ in range(5)]
        self.mark(loop)
        self.decjz(y, y_yes, done)
        self.mark(y_yes)
        self.decjz(x, x_yes, lt)
        self.mark(x_yes)
        self.goto(loop)
        self.mark(lt)
        self.inc(out)
        self.mark(done)
        self.clear(x)
        self.clear(y)
        self.release(x,y)

    def subtract(self,a,b,out):
        self.copy(a,out)
        t=self.temp()
        self.copy(b,t)
        loop, body, more, done = [self.label() for _ in range(4)]
        self.mark(loop)
        self.decjz(t,body,done)
        self.mark(body)
        self.decjz(out,more,done)
        self.mark(more)
        self.goto(loop)
        self.mark(done)
        self.clear(t)
        self.release(t)

    def divmod(self,a,b,q,r):
        require(len({a,b,q,r})==4,'division alias')
        self.clear(q)
        self.copy(a,r)
        test,t = self.temp(), self.temp()
        self.copy(b,test)
        positive, loop, stop, body = [self.label() for _ in range(4)]
        self.decjz(test,positive,self.error_label)
        self.mark(positive)
        self.clear(test)
        self.mark(loop)
        self.compare_lt(r,b,test)
        self.decjz(test,stop,body)
        self.mark(body)
        self.subtract(r,b,t)
        self.clear(r)
        self.move_add(t,r)
        self.inc(q)
        self.goto(loop)
        self.mark(stop)
        self.release(test,t)

    def expr(self,node):
        out=self.temp()
        if isinstance(node,ast.Constant):
            require(type(node.value) is int and node.value>=0,'nonnegative integer only')
            self.set(out,node.value)
        elif isinstance(node,ast.Name):
            self.copy(self.reg(node.id),out)
        elif isinstance(node,ast.BinOp):
            a,b=self.expr(node.left),self.expr(node.right)
            if isinstance(node.op,ast.Add):
                self.copy(a,out)
                self.move_add(b,out)
            elif isinstance(node.op,ast.Sub):
                self.subtract(a,b,out)
            elif isinstance(node.op,ast.Mult):
                self.multiply(a,b,out)
            elif isinstance(node.op,(ast.FloorDiv,ast.Mod)):
                other=self.temp()
                self.divmod(a,b,out,other) if isinstance(node.op,ast.FloorDiv) else self.divmod(a,b,other,out)
                self.clear(other)
                self.release(other)
            elif isinstance(node.op,ast.Pow):
                self.set(out,1)
                t=self.temp()
                loop,body,done=[self.label() for _ in range(3)]
                self.mark(loop)
                self.decjz(b,body,done)
                self.mark(body)
                self.multiply(out,a,t)
                self.clear(out)
                self.move_add(t,out)
                self.goto(loop)
                self.mark(done)
                self.release(t)
            else:
                raise ValueError('unsupported binary operation')
            self.clear(a); self.clear(b)
            self.release(a,b)
        elif isinstance(node,ast.Compare):
            require(len(node.ops)==1,'chained comparison unsupported')
            a,b=self.expr(node.left),self.expr(node.comparators[0])
            op=node.ops[0]
            if isinstance(op,(ast.Lt,ast.Gt,ast.LtE,ast.GtE)):
                if isinstance(op,(ast.Gt,ast.LtE)):
                    a,b=b,a
                self.compare_lt(a,b,out)
                if isinstance(op,(ast.LtE,ast.GtE)):
                    yes,no,end=[self.label() for _ in range(3)]
                    self.decjz(out,yes,no)
                    self.mark(no); self.inc(out); self.goto(end)
                    self.mark(yes); self.goto(end)
                    self.mark(end)
            elif isinstance(op,(ast.Eq,ast.NotEq)):
                lt,gt=self.temp(),self.temp()
                self.compare_lt(a,b,lt)
                self.compare_lt(b,a,gt)
                self.move_add(gt,lt)
                self.clear(out)
                neq,eq,end=[self.label() for _ in range(3)]
                self.decjz(lt,neq,eq)
                self.mark(eq)
                if isinstance(op,ast.Eq): self.inc(out)
                self.goto(end)
                self.mark(neq)
                if isinstance(op,ast.NotEq): self.inc(out)
                self.goto(end)
                self.mark(end)
                self.release(lt,gt)
            else:
                raise ValueError('comparison unsupported')
            self.clear(a); self.clear(b)
            self.release(a,b)
        else:
            raise ValueError(f'unsupported expression {ast.dump(node)}')
        return out

    def statements(self,body):
        for stmt in body:
            if isinstance(stmt,ast.Assign):
                require(len(stmt.targets)==1 and isinstance(stmt.targets[0],ast.Name),'name assignment only')
                # Register target before acquiring transient registers.
                target=self.reg(stmt.targets[0].id)
                t=self.expr(stmt.value)
                self.clear(target)
                self.move_add(t,target)
                self.release(t)
            elif isinstance(stmt,ast.AugAssign):
                require(isinstance(stmt.target,ast.Name),'name augmented assignment only')
                target=self.reg(stmt.target.id)
                t=self.expr(ast.BinOp(left=stmt.target,op=stmt.op,right=stmt.value))
                self.clear(target); self.move_add(t,target); self.release(t)
            elif isinstance(stmt,ast.If):
                yes,no,end=[self.label() for _ in range(3)]
                t=self.expr(stmt.test)
                self.decjz(t,yes,no); self.release(t)
                self.mark(yes); self.statements(stmt.body); self.goto(end)
                self.mark(no); self.statements(stmt.orelse)
                self.mark(end)
            elif isinstance(stmt,ast.While):
                require(not stmt.orelse,'while else unsupported')
                start,body,end=[self.label() for _ in range(3)]
                self.mark(start)
                t=self.expr(stmt.test)
                self.decjz(t,body,end); self.release(t)
                self.mark(body); self.statements(stmt.body); self.goto(start)
                self.mark(end)
            elif isinstance(stmt,ast.Return):
                t=self.expr(stmt.value)
                self.clear(self.output); self.move_add(t,self.output); self.release(t)
                self.goto(self.halt_label)
            else:
                raise ValueError(f'unsupported statement {ast.dump(stmt)}')

    def compile(self,source):
        tree=ast.parse(source)
        require(len(tree.body)==1 and isinstance(tree.body[0],ast.FunctionDef),'one function required')
        fn=tree.body[0]
        require([a.arg for a in fn.args.args]==['n'],'input must be n')
        # Allocate all named registers first to exclude accidental scratch aliasing.
        for node in ast.walk(fn):
            if isinstance(node,ast.Name): self.reg(node.id)
        self.statements(fn.body)
        self.goto(self.halt_label)
        self.mark(self.error_label)
        self.inc(self.error)
        self.goto(self.halt_label)
        self.mark(self.halt_label)
        self.emit('HALT')
        result=[]
        for pc,ins in enumerate(self.code):
            if ins[0]=='INC': result.append(('INC',ins[1],pc+1))
            elif ins[0]=='JUMP': result.append(('JUMP',self.labels[ins[1]]))
            elif ins[0]=='DECJZ': result.append(('DECJZ',ins[1],self.labels[ins[2]],self.labels[ins[3]]))
            else: result.append(ins)
        for pc,ins in enumerate(result):
            targets=ins[2:] if ins[0] in ('INC','DECJZ') else ins[1:] if ins[0]=='JUMP' else ()
            require(pc not in targets,'control self-transition would cancel its guard')
        return {'registers':self.names,'register_count':self.nreg,'instructions':result,'input':self.input,'output':self.output,'error':self.error}


def primes(count):
    ps=[]; x=2
    while len(ps)<count:
        if all(x%p for p in ps if p*p<=x): ps.append(x)
        x+=1
    return ps


def fractions(program):
    """Each row stores [numerator,denominator,source_pc,chosen_branch]."""
    code=program['instructions']; r=program['register_count']
    ps=primes(r+len(code)); data=ps[:r]; state=ps[r:]
    rows=[]
    for pc,ins in enumerate(code):
        op=ins[0]; q=state[pc]
        if op=='INC': rows.append([data[ins[1]]*state[ins[2]],q,pc,'INC'])
        elif op=='JUMP': rows.append([state[ins[1]],q,pc,'JUMP'])
        elif op=='DECJZ':
            rows.append([state[ins[2]],data[ins[1]]*q,pc,'NONZERO'])
            rows.append([state[ins[3]],q,pc,'ZERO'])
        elif op!='HALT': raise ValueError(op)
    for num,den,_,_ in rows:
        require(math.gcd(num,den)==1,'fraction cancellation detected')
    return {'data_primes':data,'state_primes':state,'fractions':rows}


def run_register(program,n,limit=100_000_000):
    require(type(n) is int and n>=0,'nonnegative input')
    regs=[0]*program['register_count']; regs[program['input']]=n
    code=program['instructions']; pc=0; steps=0
    while steps<limit:
        ins=code[pc]; op=ins[0]
        if op=='HALT': return {'output':regs[program['output']],'error':regs[program['error']],'steps':steps,'pc':pc,'registers':regs}
        if op=='INC': regs[ins[1]]+=1; pc=ins[2]
        elif op=='JUMP': pc=ins[1]
        elif op=='DECJZ':
            if regs[ins[1]]: regs[ins[1]]-=1; pc=ins[2]
            else: pc=ins[3]
        else: raise ValueError(op)
        steps+=1
    raise RuntimeError(f'step bound {limit} reached')


def run_fractran_literal(program,n,limit=100_000,integer_bit_limit=1_000_000):
    """Literal ordered fraction scan of the actual integer, not a counter shortcut."""
    f=fractions(program); ps=f['data_primes']; qs=f['state_primes']
    value=qs[0]*ps[program['input']]**n
    steps=0; maxbits=value.bit_length()
    while steps<limit:
        for num,den,_,_ in f['fractions']:
            if value%den==0:
                value=value//den*num
                break
        else:
            out=0; p=ps[program['output']]
            while value%p==0: out+=1; value//=p
            return {'output':out,'steps':steps,'max_integer_bits':maxbits}
        steps+=1; maxbits=max(maxbits,value.bit_length())
        require(maxbits<=integer_bit_limit,'literal integer bit bound reached')
    raise RuntimeError(f'literal fraction step bound {limit} reached')


def run_fractran_exponents(program,n,limit=100_000_000):
    """Exact first-eligible scan in prime-exponent representation, independent of IR dispatch.

    Static denominator indexing omits only fractions whose control prime is absent;
    it retains the source-file order of every remaining fraction.
    """
    f=fractions(program); ps=f['data_primes']; qs=f['state_primes']
    allps=ps+qs; index={p:i for i,p in enumerate(allps)}
    def factor(x):
        out=[]
        for p in allps:
            if p*p>x: break
            while x%p==0: out.append(index[p]); x//=p
        if x>1: out.append(index[x])
        return out
    rows=[]; bycontrol={}
    for num,den,pc,branch in f['fractions']:
        neg=factor(den); pos=factor(num)
        require(len([j for j in neg if j>=len(ps)])==1,'exactly one denominator control required')
        rows.append((neg,pos))
        source_control=next(j for j in neg if j>=len(ps))
        require(source_control==len(ps)+pc,'control metadata mismatch')
        bycontrol.setdefault(source_control,[]).append(len(rows)-1)
    exps=[0]*len(allps); exps[program['input']]=n; control=len(ps); exps[control]=1
    steps=0
    while steps<limit:
        for ridx in bycontrol.get(control,[]):
            neg,pos=rows[ridx]
            if all(exps[i]>0 for i in neg):
                for i in neg: exps[i]-=1
                for i in pos: exps[i]+=1
                control=next(i for i in pos if i>=len(ps))
                break
        else:
            return {'output':exps[program['output']],'error':exps[program['error']],'steps':steps,'pc':control-len(ps),'registers':exps[:len(ps)]}
        steps+=1
    raise RuntimeError(f'exponent fraction step bound {limit} reached')


def write_artifacts(program,stem):
    stem=Path(stem); stem.parent.mkdir(parents=True,exist_ok=True)
    f=fractions(program)
    # Compact serialized program, and ordered human-readable fractions, split at 80KB.
    payload={**program,'data_primes':f['data_primes'],'state_primes':f['state_primes']}
    (stem.with_suffix('.ir.json')).write_text(json.dumps(payload,separators=(',',':'))+'\n')
    lines=[f'{a}/{b}\n' for a,b,_,_ in f['fractions']]
    chunks=[]; chunk=''
    for line in lines:
        if len(chunk)+len(line)>80000: chunks.append(chunk); chunk=''
        chunk+=line
    if chunk: chunks.append(chunk)
    for j,chunk in enumerate(chunks): (stem.parent/f'{stem.name}.fractions.{j:02d}.txt').write_text(chunk)
    return {'registers':program['register_count'],'states':len(program['instructions']),'fractions':len(f['fractions']),'fraction_parts':len(chunks),'fraction_bytes':sum(map(len,chunks)),'input_prime':f['data_primes'][program['input']],'output_prime':f['data_primes'][program['output']],'entry_prime':f['state_primes'][0],'halt_prime':f['state_primes'][-1]}
