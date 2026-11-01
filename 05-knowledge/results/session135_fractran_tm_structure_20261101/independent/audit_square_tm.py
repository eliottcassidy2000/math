#!/usr/bin/env python3
"""Independent execution of expanded TSV transition tables, not source macros."""
import gzip,itertools,json,math,re,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]/'square_machines'
def check(ok,msg):
    if not ok:raise RuntimeError(msg)
def read_table(mode):
    out={};digest=hashlib.sha256()
    m=json.loads((ROOT/(mode+'.tm.json')).read_text())
    expected={(q,''.join(chars)) for q in m['rules'] for chars in itertools.product(*m['alphabets'])}
    def expanded():
        for q in sorted(m['rules']):
            for chars in itertools.product(*m['alphabets']):
                candidates=[]
                for read,t,write,moves in m['rules'][q]:
                    if all(r=='*' or r==s for r,s in zip(read,chars)):
                        candidates.append((t,''.join(s if w=='=' else w for s,w in zip(chars,write)),''.join(moves)))
                check(len(candidates)<=1,'overlapping clauses')
                t,w,moves=candidates[0] if candidates else ('REJECT',''.join(chars),'S'*6)
                yield '\t'.join((q,''.join(chars),t,w,moves))+'\n'
    path=ROOT/(mode+'.tm.tsv.gz')
    stream=gzip.open(path,'rt') if path.exists() else expanded()
    for line in stream:
        digest.update(line.encode())
        q,r,t,w,moves=line.rstrip('\n').split('\t');key=(q,r)
        check(key not in out,('duplicate',key));check(len(r)==len(w)==len(moves)==6,'arity')
        check('*' not in r and '=' not in w,'unexpanded macro');out[key]=(t,w,moves)
    if hasattr(stream,'close'):stream.close()
    record=next(r for r in json.loads((ROOT/'TM_BUILD.json').read_text()) if r['mode']==mode)
    check(set(out)==expected,'missing/extra concrete keys')
    check(len(out)==record['concrete_delta_rows'] and digest.hexdigest()==record['expanded_sha256'],'concrete stream pin')
    return out
TABLE={m:read_table(m) for m in ('path','cycle')}
def execute(mode,text=None,state='start',tapes=None,heads=None,stop=None,limit=2000000):
    tapes=([dict(enumerate(text))]+[{} for _ in range(5)]) if tapes is None else tapes
    heads=[0]*6 if heads is None else heads;steps=0
    while state not in ('ACCEPT','REJECT') and state!=stop:
        key=state,''.join(t.get(h,'_') for t,h in zip(tapes,heads))
        check(key in TABLE[mode],('missing concrete',key))
        state,w,m=TABLE[mode][key]
        for i,(a,d) in enumerate(zip(w,m)):
            if a=='_':tapes[i].pop(heads[i],None)
            else:tapes[i][heads[i]]=a
            heads[i]+={'L':-1,'S':0,'R':1}[d]
        steps+=1;check(steps<=limit,'nontermination')
    return state,steps,tapes,heads

def expected(text,cycle):
    if not re.fullmatch(r'\^1+#1+(?:#1+)*\$',text):return False
    pieces=text[1:-1].split('#');n=len(pieces[0]);w=list(map(len,pieces[1:]))
    if cycle and n<3:return False
    if len(w)!=n or sorted(w)!=list(range(1,n+1)):return False
    edges=list(zip(w,w[1:]))+([(w[-1],w[0])] if cycle else [])
    return all(math.isqrt(a+b)**2==a+b for a,b in edges)
def encode(n,w):return '^'+'1'*n+'#'+'#'.join('1'*x for x in w)+'$'
def runtime(mode,n,w):
    edges=list(zip(w,w[1:]))+([(w[-1],w[0])] if mode=='cycle' else [])
    r=sum(math.isqrt(a+b) for a,b in edges)
    return 5*n*n+14*n+(11+3*w[0]+6*r if mode=='cycle' else 3-2*w[0]-5*w[-1]+6*r)
def main():
    malformed=structured=0
    for n in range(8):
        for tup in itertools.product('^1#$',repeat=n):
            text=''.join(tup)
            for mode in TABLE:
                state,steps,*_=execute(mode,text)
                check((state=='ACCEPT')==expected(text,mode=='cycle'),('word',mode,text,state));malformed+=1
    for n in range(1,8):
        for w in itertools.permutations(range(1,n+1)):
            text=encode(n,w)
            for mode in TABLE:
                state,steps,*_=execute(mode,text)
                yes=expected(text,mode=='cycle')
                check((state=='ACCEPT')==yes,('permutation',mode,n,w))
                if yes:check(steps==runtime(mode,n,w),'runtime tiny')
                structured+=1
    subroutine=0
    for s in range(1,401):
        tapes=[{0:'^'} for _ in range(6)];tapes[4].update({i:'1' for i in range(1,s+1)})
        heads=[0]*6;heads[4]=s+1;heads[5]=1
        state,steps,tapes,heads=execute('path',state='edge_square',tapes=tapes,heads=heads,stop='boundary')
        k=math.isqrt(s);yes=k*k==s
        check((state=='boundary')==yes,('subroutine',s,state))
        if yes:
            check(steps==2*s+6*k-2,('square cost',s,steps))
            check(tapes[4]=={0:'^'} and tapes[5]=={0:'^'} and heads[4:]==[1,1],('scratch clear',s))
        subroutine+=1
    path=[8,1,15,10,6,3,13,12,4,5,11,14,2,7,9]
    cycle=[1,8,28,21,4,32,17,19,30,6,3,13,12,24,25,11,5,31,18,7,29,20,16,9,27,22,14,2,23,26,10,15]
    positive=[]
    for mode,w in [('path',path),('path',path[::-1])]+[('cycle',cycle[j:]+cycle[:j]) for j in range(32)]+[('cycle',cycle[::-1])]:
        state,steps,*_=execute(mode,encode(len(w),w))
        check(state=='ACCEPT' and steps==runtime(mode,len(w),w),('positive/cost',mode,w,steps))
        positive.append({'mode':mode,'n':len(w),'first':w[0],'last':w[-1],'steps':steps})
    print(json.dumps({'status':'PASS','concrete_rows':{m:len(t) for m,t in TABLE.items()},'all_strings_length_0_7_machine_cases':malformed,'permutations_orders_1_7_machine_cases':structured,'square_subroutine_cases':subroutine,'positive_cost_cases':positive},indent=2))
if __name__=='__main__':main()
