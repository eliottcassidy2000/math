"""Complete two-ready inverse through5993, reusing proved order-tail minima."""
from pathlib import Path
from itertools import combinations,product
from collections import Counter
import sys,json,subprocess,tempfile,shutil,time
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
OLD=ROOT/'04-computation/nonzero_preparation_20261103'
sys.path.insert(0,str(OLD))
import surgery_orders as a

def bank_through(bound):
    cycles={};pairs={};two={};bank=[];records=[]
    for n in range(3,bound+1):
        cycles[n]=sorted({a.cyclic_code(a.unrank_one(n,i)) for i in range(a.U[n])})
        pairs[n]=[]
        for T in cycles[n]:
            low={i for i,row in enumerate(T) if row.bit_count()==1}
            entries=sorted(low|{(i+1)%n for i in low})
            for J in ([(0,1)] if n==3 else combinations(entries,2)):
                pairs[n].append((T,tuple(J)))
        a.need(len(pairs[n])==a.P[n],'complete entry-pair class count')
        d={}
        for m in range(3,n):
            t=n-m-1
            suppliers=[(T,J,[tuple(range(j,m))+tuple(range(j)) for j in J]) for T,J in pairs[m]]
            suppliers += [(T,(0,),paths) for T,paths in two.get(m,{}).items()]
            for T,J,paths in suppliers:
                for pre in product((0,1),repeat=max(t-1,0)):
                    word=pre+(1,) if t else ()
                    Q,r=a.ready.funnel(T,J,word)
                    hp=[(r,)+tuple(path)+tuple(range(m,m+t)) for path in paths]
                    code,path=min((a.induced(Q,path),path) for path in hp)
                    inverse={v:i for i,v in enumerate(path)}
                    canon=tuple(sorted(tuple(inverse[v] for v in old) for old in hp))
                    if code in d:a.need(d[code]==canon,'canonical path transport')
                    d[code]=canon
        a.need(len(d)==a.V[n],'complete two-start class count')
        two[n]=d
        bank.extend((n,T,J,'pair') for T,J in pairs[n])
        bank.extend((n,T,(0,),'single') for T in sorted(d))
        records.append(dict(n=n,pairs=len(pairs[n]),two=len(d),cycles=len(cycles[n])))
        print('bank',records[-1],flush=True)
    return bank,records

def main():
    targets=(221,419,665,667,1997,1999,2003,2005,5993)
    hits={str(h):[] for h in targets};hist={n:Counter() for n in range(4,15)}
    bank,records=bank_through(12)
    with tempfile.TemporaryDirectory(prefix='two-ready-inverse-') as td:
        td=Path(td);exe=td/'counter.exe';inp=td/'input.txt';out=td/'output.txt'
        compiler=shutil.which('g++') or shutil.which('clang++') or shutil.which('c++')
        a.need(compiler is not None,'available independent exact C++ counter')
        subprocess.run([compiler,'-O3','-std=c++17',str(OLD/'surgery_orders.cpp'),'-o',str(exe)],check=True)
        with inp.open('w') as f:
            for n,T,J,_ in bank:f.write(' '.join(map(str,(n,sum(1<<j for j in J),*T)))+'\n')
        with inp.open('rb') as f,out.open('wb') as g:
            counter=subprocess.run([str(exe)],stdin=f,stdout=g,stderr=subprocess.PIPE,check=True)
        for hostile in ('3 3 3 4 1\n','3 1 2 4 1\n'):
            rejected=subprocess.run([str(exe)],input=hostile.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            a.need(rejected.returncode==2 and b'REJECT' in rejected.stderr,'hostile adjacency and start-budget controls')
        print('counter',counter.stderr.decode().strip(),flush=True)
        lines=out.read_text().splitlines();a.need(len(lines)==len(bank),'every finite core counted')
        for index,((m,T,J,kind),line) in enumerate(zip(bank,lines)):
            data=dict(zip(('h','p2','e','F','G','b','C'),map(int,line.split())))
            if index in(0,len(bank)-1):a.need(data==a.ready.core_data(T,J),'independent Python response controls')
            for n in hist:
                length=n-m-2
                if length<0:continue
                for pre in product((0,1),repeat=max(length-1,0)):
                    word=pre+(1,) if length else ()
                    h=a.ready.ready_count(data,word);hist[n][h]+=1
                    if h in targets:
                        parent,r=a.ready.funnel(T,J,word)
                        hits[str(h)].append(dict(order=n,core=T,ports=J,word=word,kind=kind,parent=parent,root=r))
        for n in hist:
            a.need(sum(hist[n].values())==a.V[n-1],'all degree-one ready classes at target order')
            exceptional=[3*2**(n-3)-1]
            if n>=5:exceptional.extend((6*2**(n-4)+1,5*2**(n-4)-1))
            if n>=6:exceptional.extend((13*2**(n-5)+3,14*2**(n-5)+3))
            for h in exceptional:
                hist[n][h]+=1
                if h in targets:hits[str(h)].append(dict(order=n,exception=True))
    sys.path.insert(0,str(HERE.parent.parent/'uniform_return_metagraph_20261111'/'reset'))
    from classify_modules import inherited,strong
    for h,records_at_h in hits.items():
        for rec in records_at_h:
            if rec.get('exception'):continue
            child=inherited.birth_beaten_only_by(rec['parent'],rec['root'])
            profile=inherited.root_profiles(child)[-1]
            a.need(strong(child) and profile['h']==int(h) and profile['gap']==2,'independent literal positive ready witness')
    data=dict(status='PASS',scope='Every strong two-ready rooted class of orders4..14',
              bank_rows=len(bank),grammar=records,counter=json.loads(counter.stderr),
              ranks={str(n):dict(classes=sum(d.values()),minimum=min(d)) for n,d in hist.items()},hits=hits)
    (HERE/'ready_two_inverse.json').write_text(json.dumps(data,indent=2)+'\n')
    print('HITS',{h:len(v) for h,v in hits.items()},flush=True)

if __name__=='__main__':main()
