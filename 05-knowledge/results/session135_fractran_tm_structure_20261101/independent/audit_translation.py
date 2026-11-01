#!/usr/bin/env python3
"""Independent interpretation of emitted TM and reverse FRACTRAN tables."""
import importlib.util,json,math,hashlib,sys
sys.dont_write_bytecode=True
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]/'machine_translation'
def check(ok,msg):
    if not ok:raise RuntimeError(msg)
def step_ir(op,v):
    v=list(v);kind=op[0]
    if kind=='halt':return None,v,0
    if kind=='goto':return op[1],v,1
    r=op[1];cost=2*v[r]+3
    if kind=='inc':v[r]+=1;return op[2],v,cost
    check(kind=='decjz','unknown op')
    if v[r]:v[r]-=1;return op[2],v,cost
    return op[3],v,cost
class TM:
    def __init__(self,m):
        self.m=m;self.packed=m['model']=='packed_single';self.table={}
        for q,s,t,w,d in m['rows']:
            check((q,s) not in self.table,('duplicate',q,s));self.table[q,s]=t,w,d
        check(all('*' not in str(s) for q,s in self.table),'wildcard')
    def init(self,q,v):
        self.q=q;self.head=[0]*self.m['tapes'];self.steps=0
        if self.packed:self.tape=[{0:'#',**{j:sum((1<<i) for i,a in enumerate(v) if a>=j) for j in range(1,max(v,default=0)+1)}}]
        else:self.tape=[{0:'#',**{j:'1' for j in range(1,a+1)}} for a in v]
    def run(self,stops):
        while True:
            key=(self.q,self.tape[0].get(self.head[0],0)) if self.packed else (self.q,''.join(t.get(h,'_') for t,h in zip(self.tape,self.head)))
            if key not in self.table:return None,self.values(),self.steps
            self.q,write,moves=self.table[key]
            if self.packed:write=[write]
            for i,(w,d) in enumerate(zip(write,moves)):
                if w==self.m['blank']:self.tape[i].pop(self.head[i],None)
                else:self.tape[i][self.head[i]]=w
                self.head[i]+={'L':-1,'R':1,'S':0}[d]
            self.steps+=1;check(self.steps<100000,'test budget')
            if self.q in stops:return self.q,self.values(),self.steps
    def values(self):
        check(not any(self.head),'not at marker')
        d=self.m['registers'];v=[0]*d
        if self.packed:
            for i in range(d):
                while self.tape[0].get(v[i]+1,0)&(1<<i):v[i]+=1
            expected=[{0:'#',**{j:sum((1<<i) for i,a in enumerate(v) if a>=j) for j in range(1,max(v,default=0)+1)}}]
        else:
            for i,t in enumerate(self.tape):
                while t.get(v[i]+1)=='1':v[i]+=1
            expected=[{0:'#',**{j:'1' for j in range(1,a+1)}} for a in v]
        check(self.tape==expected,'noncanonical tapes')
        return v
def rational_step(n,rows):
    for j,(a,b) in enumerate(rows):
        if n*a%b==0:return n*a//b,j
    return None,None

def main():
    primitive=reverse=sourcecases=0;receipts={}
    for name in ('pure_blue','compact_tail21','migration_face','migration_full'):
        p=json.loads((ROOT/(name+'.ir.json')).read_text());m=json.loads((ROOT/(name+'.tm.json')).read_text());
        if 'row_shards' in m:
            assembled=[]
            for shard in m['row_shards']:
                raw=(ROOT/shard['file']).read_bytes();check(len(raw)==shard['bytes'] and hashlib.sha256(raw).hexdigest()==shard['sha256'],'shard authentication')
                entries=json.loads(raw);check(len(entries)==shard['rows'],'shard row count');assembled.extend(entries)
            check(len(assembled)==m['rows_count'],'whole row count');m['rows']=assembled
        back=json.loads((ROOT/(name+'.roundtrip.json')).read_text())
        tm=TM(m);ir=p['instructions'];d=len(p['register_primes']);profiles=[[0]*d,[1]*d,[2+(i%3) for i in range(d)]]
        for q,op in ir.items():
            for v in profiles:
                target,w,cost=step_ir(op,v);tm.init(q,v);got=tm.run(set(ir))
                check(got==(target,w,cost),('TM instruction',name,q,v,got,(target,w,cost)));primitive+=1
                value=back['controls'][q]*math.prod(p**e for p,e in zip(back['register_primes'],v))
                if target is None:check(rational_step(value,back['fractions'])[0] is None,'halt row')
                else:
                    value,_=rational_step(value,back['fractions']);value,_=rational_step(value,back['fractions'])
                    expected=back['controls'][target]*math.prod(p**e for p,e in zip(back['register_primes'],w))
                    check(value==expected,('roundtrip',name,q,v))
                reverse+=1
        rows=[tuple(map(int,s.split('/'))) for s in (ROOT/(name+'.fractran')).read_text().splitlines() if s and not s.startswith('#')]
        for q in p['boundaries']:
            for v in profiles:
                integer=int(q[1:])*math.prod(a**e for a,e in zip(p['register_primes'],v))
                nxt,row=rational_step(integer,rows)
                tm.init(q,v);state,w,cost=tm.run(set(p['boundaries']))
                if nxt is None:check(state is None,('source halt',name,q))
                else:
                    out=int(state[1:])*math.prod(a**e for a,e in zip(p['register_primes'],w))
                    check(out==nxt,('source macro',name,q,v,row))
                sourcecases+=1
        receipts[name]={'states':len(m['states']),'rows':len(m['rows']),'registers':d}
    # The general compiler is separately exercised on rollback and repeated guards.
    spec=importlib.util.spec_from_file_location('producer_machines',ROOT/'machines.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    lists=[[],[(1,1)],[(2,2)],[(5,6),(7,2)],[(5,12),(7,2)],[(3,8),(5,4),(1,1)],[(7,45),(11,9),(13,3)],[(12,25),(3,5)],[(1,2)],[(8,9),(9,8)]]
    general=0
    for rows in lists:
        p=mod.general_fractran_ir(rows);m=mod.to_tm(p,'unary_multi');tm=TM(m)
        for integer in range(1,101):
            residual=integer;v=[]
            for prime in p['register_primes']:
                e=0
                while residual%prime==0:residual//=prime;e+=1
                v.append(e)
            tm.init('STEP',v);q,w,cost=tm.run({'STEP','HALT'})
            expected,_=rational_step(integer,rows)
            if expected is None:check(q in (None,'HALT'),('general halt',rows,integer,q))
            else:check(q=='STEP' and residual*math.prod(a**e for a,e in zip(p['register_primes'],w))==expected,('general',rows,integer,q,w))
            general+=1
    # Fresh trampoline guards must survive genuine control self-loops.
    loops={'register_primes':[2], 'control_primes':[], 'boundaries':['x'], 'instructions':{'x':['inc',0,'x']}}
    loopback=mod.ir_to_fractran(loops)
    for n in range(8):
        value=loopback['controls']['x']*2**n
        for _ in range(2):value,_=rational_step(value,loopback['fractions'])
        check(value==loopback['controls']['x']*2**(n+1),'self inc guard')
    print(json.dumps({'status':'PASS','serialized_TM_instruction_cases':primitive,'serialized_reverse_fraction_cases':reverse,'serialized_source_macro_cases':sourcecases,'general_FRACTRAN_cases':general,'self_loop_cases':8,'artifacts':receipts},indent=2))
if __name__=='__main__':main()
