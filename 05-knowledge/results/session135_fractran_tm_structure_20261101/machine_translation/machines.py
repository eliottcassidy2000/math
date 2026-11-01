"""Explicit FRACTRAN, register and actual finite TM compilers. Stdlib only."""
from math import prod,gcd
from itertools import product
import json
from pathlib import Path

def need(ok,msg='condition failed'):
    if not ok: raise ValueError(msg)

def factor(n):
    need(type(n)==int and n>0,'positive integer required')
    out={};p=2
    while p*p<=n:
        while n%p==0:out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=out.get(n,0)+1
    return out

def load_fractions(path):
    return [tuple(map(int,t.split('/'))) for s in Path(path).read_text().splitlines() if (t:=s.split('#')[0].strip())]

def decode(integer, data, controls):
    need(integer>0,'positive integer required')
    values=[]
    for p in data:
        e=0
        while integer%p==0:e+=1;integer//=p
        values.append(e)
    present=[]
    for p in controls:
        e=0
        while integer%p==0:e+=1;integer//=p
        if e:present.append((p,e))
    need(len(present)==1 and present[0][1]==1,'one-hot control required')
    return 'C'+str(present[0][0]),values,integer

def encode(q,v,data,inert=1):
    return int(q[1:])*inert*prod(p**e for p,e in zip(data,v))

def fraction_step(n,rows):
    for i,(a,b) in enumerate(rows):
        if n*a%b==0:return n*a//b,i
    return None,None

def to_ir(rows,data,controls):
    """One-hot-control subclass: one unguarded row or guard-then-fallback."""
    need(len(set(data+controls))==len(data+controls),'disjoint prime roles')
    need(all(factor(p)=={p:1} for p in data+controls),'roles must be primes')
    groups={q:[] for q in controls};parsed=[]
    for k,(a,b) in enumerate(rows):
        need(gcd(a,b)==1,'input fractions must be reduced')
        na,nb=factor(a),factor(b)
        ca=[p for p in na if p in controls];cb=[p for p in nb if p in controls]
        need(len(ca)==len(cb)==1 and na[ca[0]]==nb[cb[0]]==1,'one control each side')
        dst,src=ca[0],cb[0];del na[dst];del nb[src]
        need(set(na)<=set(data) and set(nb)<=set(data),'unknown data prime')
        need(not nb or (len(nb)==1 and next(iter(nb.values()))==1),'unit single guard only')
        guard=None if not nb else data.index(next(iter(nb)))
        inc=[i for i,p in enumerate(data) for _ in range(na.get(p,0))]
        rec={'index':k,'src':src,'dst':dst,'guard':guard,'inc':inc}
        groups[src].append(rec);parsed.append(rec)
    ir={};boundaries=['C'+str(q) for q in controls]
    def chain(row,entry=None):
        inc=row['inc'];target='C'+str(row['dst'])
        if not inc:
            if entry is not None:ir[entry]=['goto',target];return entry
            return target
        labs=[entry if entry is not None else 'F'+str(row['index'])+':0']+['F'+str(row['index'])+':'+str(i) for i in range(1,len(inc))]
        for i,r in enumerate(inc):ir[labs[i]]=['inc',r,labs[i+1] if i+1<len(inc) else target]
        return labs[0]
    for q,rs in groups.items():
        entry='C'+str(q)
        if not rs:ir[entry]=['halt'];continue
        if len(rs)==1:
            need(rs[0]['guard'] is None,'single guarded row unsupported')
            chain(rs[0],entry)
        else:
            need(len(rs)==2 and rs[0]['guard'] is not None and rs[1]['guard'] is None,'guard must precede fallback')
            ir[entry]=['decjz',rs[0]['guard'],chain(rs[0]),chain(rs[1])]
    return {'register_primes':data,'control_primes':controls,'boundaries':boundaries,'instructions':ir,'source_rows':parsed}

def ir_step(q,v,ir):
    op=ir[q];v=list(v)
    if op[0]=='halt':return q,v,0
    if op[0]=='goto':return op[1],v,1
    r=op[1];cost=2*v[r]+3
    if op[0]=='inc':v[r]+=1;return op[2],v,cost
    if v[r]:v[r]-=1;return op[2],v,cost
    return op[3],v,cost

def to_tm(program,kind):
    """No wildcards: every emitted row is a single full TM transition."""
    ir=program['instructions'];d=len(program['register_primes']);rows=[]
    need(kind in ('packed_single','unary_multi'),'unknown backend')
    if kind=='packed_single':
        alphabet=['#']+list(range(1<<d));marker='#'
    else:alphabet=['#','1','_'];marker='#'*d
    def emit(q,s,t,w,m):rows.append([q,s,t,w,m])
    def movement(r,direction):
        if kind=='packed_single':return direction
        return ''.join(direction if i==r else 'S' for i in range(d))
    def symbols(r):
        if kind=='packed_single':return [(s,bool(s>>r&1),s|(1<<r),s&~(1<<r)) for s in range(1<<d)]
        def at(x):return '#'*r+x+'#'*(d-r-1)
        return [(at('_'),False,at('1'),at('_')),(at('1'),True,at('1'),at('_'))]
    stay='S' if kind=='packed_single' else 'S'*d
    for q,op in ir.items():
        if op[0]=='halt':continue
        if op[0]=='goto':emit(q,marker,op[1],marker,stay);continue
        r=op[1];first=q+'|first';scan=q+'|scan';back=q+'|back';zero=q+'|zero';erase=q+'|erase'
        emit(q,marker,first,marker,movement(r,'R'))
        for s,bit,setbit,clearbit in symbols(r):
            if op[0]=='inc':
                emit(first,s,first if bit else back,s if bit else setbit,movement(r,'R' if bit else 'L'))
            else:
                emit(first,s,scan if bit else zero,s,movement(r,'R' if bit else 'L'))
                emit(scan,s,scan if bit else erase,s,movement(r,'R' if bit else 'L'))
                if bit:emit(erase,s,back,clearbit,movement(r,'L'))
            emit(back,s,back,s,movement(r,'L'))
        emit(back,marker,op[2],marker,stay)
        if op[0]=='decjz':emit(zero,marker,op[3],marker,stay)
    states=sorted(set(ir)|{r[0] for r in rows}|{r[2] for r in rows})
    return {'model':kind,'tapes':1 if kind=='packed_single' else d,'registers':d,
            'alphabet':alphabet,'blank':0 if kind=='packed_single' else '_',
            'marker':'#','moves':'L,R,S; every row counts as one transition',
            'states':states,'boundaries':program['boundaries'],'instruction_states':list(ir),'rows':rows}

class LiteralTM:
    """Reads ONLY the serialized finite transition table after initialization."""
    def __init__(self,machine,q,v):
        self.m=machine;self.q=q;self.time=0;self.head=[0]*machine['tapes'];self.table={}
        for a,b,c,d,e in machine['rows']:
            need((a,b) not in self.table,'duplicate transition')
            self.table[a,b]=(c,d,e)
        self.tape=[]
        if machine['model']=='packed_single':
            t={0:'#'}
            for x in range(max(v,default=0)):
                t[x+1]=sum((1<<i) for i,e in enumerate(v) if e>x)
            self.tape=[t]
        else:
            self.tape=[{0:'#',**{j+1:'1' for j in range(e)}} for e in v]
    def step(self):
        m=self.m;read=[t.get(h,m['blank']) for t,h in zip(self.tape,self.head)]
        key=read[0] if m['model']=='packed_single' else ''.join(read)
        row=self.table.get((self.q,key))
        if row is None:return False
        self.q,w,move=row
        writes=[w] if m['model']=='packed_single' else list(w)
        for i,(symbol,direction) in enumerate(zip(writes,move)):
            if symbol==m['blank']:self.tape[i].pop(self.head[i],None)
            else:self.tape[i][self.head[i]]=symbol
            self.head[i]+={'L':-1,'R':1,'S':0}[direction]
        self.time+=1;return True
    def values(self):
        need(all(h==0 for h in self.head),'boundary head not at marker')
        need(all(t.get(0)=='#' for t in self.tape),'marker destroyed')
        v=[0]*self.m['registers']
        if self.m['model']=='packed_single':
            for i in range(len(v)):
                x=1
                while isinstance(self.tape[0].get(x,0),int) and self.tape[0].get(x,0)>>i&1:v[i]+=1;x+=1
            expected={0:'#',**{x+1:sum(1<<i for i,e in enumerate(v) if e>x) for x in range(max(v,default=0))}}
            need(self.tape[0]==expected,'malformed packed unary tape')
        else:
            for i,t in enumerate(self.tape):
                x=1
                while t.get(x)=='1':v[i]+=1;x+=1
                need(t=={0:'#',**{j+1:'1' for j in range(v[i])}},'malformed unary tape')
        return v
    def next_boundary(self,limit=10000000):
        before=self.time
        while self.step():
            if self.q in self.m['boundaries']:return self.q,self.values(),self.time-before
            need(self.time-before<limit,'TM transition budget exhausted')
        return None,None,self.time-before

def primes(n):
    out=[];p=2
    while len(out)<n:
        if all(p%q for q in out if q*q<=p):out.append(p)
        p+=1
    return out

def ir_to_fractran(program):
    """Exactly two rational steps per IR instruction, including self loops."""
    ir=program['instructions'];d=len(program['register_primes']);labels=sorted(ir)
    ps=primes(d+3*len(labels));regs=ps[:d];cs=dict(zip(labels,ps[d:d+len(labels)]));cursor=d+len(labels);rows=[]
    for q in labels:
        op=ir[q]
        if op[0]=='halt':continue
        a,b=ps[cursor:cursor+2];cursor+=2
        if op[0]=='goto':rows.extend([(a,cs[q]),(cs[op[1]],a)])
        elif op[0]=='inc':rows.extend([(a*regs[op[1]],cs[q]),(cs[op[2]],a)])
        else:rows.extend([(a,cs[q]*regs[op[1]]),(b,cs[q]),(cs[op[2]],a),(cs[op[3]],b)])
    need(all(gcd(a,b)==1 for a,b in rows),'reverse cancellation')
    return {'register_primes':regs,'controls':cs,'fractions':[list(row) for row in rows]}

def save_json(path,obj):Path(path).write_text(json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n')

def general_fractran_ir(rows):
    """Every finite positive rational list, reduced before denominator guards.
    Failed candidate guards roll back their successful decrements exactly.
    No one-hot-control assumption. STEP is one rational-step boundary.
    """
    normalized=[];factors=[];used=set()
    for a,b in rows:
        need(type(a)==type(b)==int and a>0 and b>0,'positive rational list required')
        g=gcd(a,b);a//=g;b//=g;na,nb=factor(a),factor(b)
        normalized.append((a,b));factors.append((na,nb));used|=set(na)|set(nb)
    data=sorted(used) or [2];ir={'HALT':['halt']}
    def chain(label,registers,target):
        if not registers:ir[label]=['goto',target];return label
        for j,r in enumerate(registers):
            q=label if j==0 else label+':'+str(j)
            nxt=target if j+1==len(registers) else label+':'+str(j+1)
            ir[q]=['inc',r,nxt]
        return label
    starts=['STEP']+['TRY'+str(i) for i in range(1,len(rows))]+['HALT']
    if not rows:ir['STEP']=['halt']
    for i,(na,nb) in enumerate(factors):
        den=[data.index(p) for p in data for _ in range(nb.get(p,0))]
        num=[data.index(p) for p in data for _ in range(na.get(p,0))]
        if not den:chain(starts[i],num,'STEP');continue
        success='STEP' if not num else chain('APPLY'+str(i),num,'STEP')
        for j,r in enumerate(den):
            label=starts[i] if j==0 else 'TEST'+str(i)+':'+str(j)
            pos=success if j+1==len(den) else 'TEST'+str(i)+':'+str(j+1)
            zero=starts[i+1] if j==0 else chain('UNDO'+str(i)+':'+str(j),den[:j],starts[i+1])
            ir[label]=['decjz',r,pos,zero]
    return {'register_primes':data,'control_primes':[], 'boundaries':['STEP','HALT'],
            'instructions':ir,'normalized_fractions':normalized,'initial_state':'STEP'}

def save_tm(path,machine,limit=90000):
    """UTF-8-only publication, preserving every ordinary delta row exactly."""
    import hashlib
    path=Path(path)
    def blob(x):return (json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
    whole=blob(machine)
    if len(whole)<=limit:path.write_bytes(whole);return
    # Split only the finite list of literal rows; neither guards nor macros.
    metadata={k:v for k,v in machine.items() if k!='rows'};chunks=[];current=[];size=3
    for row in machine['rows']:
        encoded=json.dumps(row,sort_keys=True,separators=(',',':')).encode()
        need(len(encoded)+3<=limit,'one row exceeds publication limit')
        extra=len(encoded)+(1 if current else 0)
        if current and size+extra>limit:chunks.append(current);current=[];size=3;extra=len(encoded)
        current.append(row);size+=extra
    if current:chunks.append(current)
    metadata['table_storage']='explicit-json-row-shards-v1';metadata['rows_count']=len(machine['rows']);metadata['row_shards']=[]
    for i,rows in enumerate(chunks):
        filename=path.name[:-5]+'.rows.'+str(i).zfill(3)+'.json';content=blob(rows)
        need(len(content)<=limit,'shard size');(path.parent/filename).write_bytes(content)
        metadata['row_shards'].append({'file':filename,'rows':len(rows),'bytes':len(content),'sha256':hashlib.sha256(content).hexdigest()})
    content=blob(metadata);need(len(content)<=limit,'TM index exceeds publication limit');path.write_bytes(content)

def load_tm(path):
    """Authenticated concatenation only; transition semantics stay literal."""
    import hashlib
    path=Path(path);m=json.loads(path.read_text())
    if 'rows' in m:return m
    need(m.get('table_storage')=='explicit-json-row-shards-v1','unknown table storage')
    rows=[]
    for rec in m['row_shards']:
        need(Path(rec['file']).name==rec['file'],'shard path must be a basename')
        content=(path.parent/rec['file']).read_bytes()
        need(len(content)==rec['bytes'] and hashlib.sha256(content).hexdigest()==rec['sha256'],'shard authentication')
        part=json.loads(content);need(isinstance(part,list) and len(part)==rec['rows'],'shard row count');rows.extend(part)
    need(len(rows)==m['rows_count'],'total row count');m['rows']=rows
    return m
