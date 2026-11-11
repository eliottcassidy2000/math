"""Child-first audit of all kernels below 125, and literal square-ray witnesses."""
from pathlib import Path
from collections import Counter
import csv,gzip,hashlib,json,shutil,subprocess,tempfile
from enumerate import hcount,scc,require

HERE=Path(__file__).resolve().parent
PREV=HERE.parent/'metagraph_aspect_sequences_20261111'

def automorphisms(rows):
    n=len(rows)
    degrees=[r.bit_count() for r in rows]
    cand=[[j for j in range(n) if degrees[j]==degrees[i]] for i in range(n)]
    order=sorted(range(n),key=lambda i:len(cand[i]))
    mapping=[-1]*n
    result=[]
    def search(k,used):
        if k==n:
            result.append(tuple(mapping))
            return
        i=order[k]
        for j in cand[i]:
            if used>>j&1:
                continue
            if any(((rows[i]>>v)&1)!=((rows[j]>>mapping[v])&1) for v in order[:k]):
                continue
            mapping[i]=j
            search(k+1,used|(1<<j))
        mapping[i]=-1
    search(0,0)
    return result

def decode(code,n=9):
    rows=[0]*n
    b=n*(n-1)//2-1
    for i in range(n):
        for j in range(i+1,n):
            if code>>b&1:rows[i]|=1<<j
            else:rows[j]|=1<<i
            b-=1
    return rows

def closure(rows,out,t):
    n=len(rows);root=n+t
    child=list(rows)+[0]*(t+1)
    for i in range(n):
        child[i] |= ((1<<t)-1)<<n
        if out>>i&1:child[root]|=1<<i
        else:child[i]|=1<<root
    for i in range(n,root):
        child[i] |= sum(1<<j for j in range(i+1,root+1))
    return child

def run():
    kernels={}
    counts=Counter()
    strong_polynomials={}
    point_count=0
    def inspect(rows,g,aut):
        nonlocal point_count
        n=len(rows)
        if n<3 or g>=125 or len(scc(rows))!=1:
            return
        require(hcount(rows)==g,'independent child path count')
        group=automorphisms(rows)
        require(len(group)==aut,'independent automorphism group')
        counts[n]+=1
        strong_polynomials.setdefault(str(g),[0]*10)[n]+=1
        for v in range(n):
            if min(p[v] for p in group)!=v:
                continue
            h=hcount(rows,((1<<n)-1)^(1<<v))
            kernels.setdefault((h,g),[0]*9)[n-1]+=1
            point_count+=1
    with gzip.open(PREV/'records.jsonl.gz','rt',encoding='utf-8') as f:
        for line in f:
            r=json.loads(line)
            inspect(r['rows'],r['H'],r['aut'])
    source9=HERE.parent/'blue_quotient_20261020b/n9_classes.tsv.gz'
    compiler=shutil.which('g++') or r'C:\Users\Eliott\scoop\apps\gcc\current\bin\g++.exe'
    with tempfile.TemporaryDirectory(prefix='birth-kernel-audit-') as td:
        td=Path(td)
        with gzip.open(source9,'rt') as f,(td/'input.txt').open('w') as dest:
            for r in csv.DictReader(f,delimiter='\t'):
                dest.write(r['code']+' '+r['Aut']+'\n')
        exe=td/'filter.exe'
        subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'filter9.cpp'),'-o',str(exe)],check=True)
        filtered=subprocess.check_output([str(exe),str(td/'input.txt'),str(td/'filtered.txt')],text=True).strip()
        for line in (td/'filtered.txt').read_text().splitlines():
            code,aut,g,*rows=map(int,line.split())
            inspect(rows,g,aut)
    exported={str(h)+','+str(g):c for (h,g),c in sorted(kernels.items())}
    expected=json.loads((HERE/'strong_birth_kernels_below125.json').read_text(encoding='utf-8'))['kernels']
    require(exported==expected,'all child-first kernel coefficients')
    (HERE/'strong_class_polynomials_below125.json').write_text(
        json.dumps(strong_polynomials,indent=2)+'\n',encoding='utf-8')
    seedbank={int(r['c']):{k:int(v) for k,v in r.items()} for r in
              csv.DictReader((PREV/'seed_witnesses9.tsv').open(),delimiter='\t')}
    squares=[]
    for d in range(1,64,2):
        r=seedbank[d*d]
        rows=decode(r['canonical_code'])
        if r['converse']:
            rows=[511^v^(1<<i) for i,v in enumerate(rows)]
        require(hcount(rows)==r['H'],'square ray base')
        exponent=r['first_exponent']
        for s in [exponent,exponent+1]:
            child=closure(rows,r['selected_mask'],s-exponent+1)
            require(hcount(child)==1+d*d*(1<<s),'two literal ray points')
        even=2*((exponent+1)//2)
        squares.append(dict(odd_k_part=d,c=d*d,first_even_exponent=even,
            first_k=d*(1<<(even//2)),first_H=1+d*d*(1<<even),
            base_rows=rows,selected_mask=r['selected_mask'],first_exponent=exponent))
    (HERE/'quadratic_tournament_rays.json').write_text(json.dumps(squares,indent=2)+'\n',encoding='utf-8')
    result=dict(status='PASS',child_first_strong_class_counts=dict(counts),
        strong_pointed_classes=point_count,kernel_count=len(kernels),filter9=filtered,
        input9_sha256=hashlib.sha256(source9.read_bytes()).hexdigest(),
        square_coefficient_rays=len(squares),literal_tail_checks=2*len(squares),
        largest_first_even_exponent=max(s['first_even_exponent'] for s in squares))
    (HERE/'kernel_audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result))

if __name__=='__main__':
    run()
