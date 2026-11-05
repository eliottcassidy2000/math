#!/usr/bin/env python3
"""Complete strong-core search plus independent exact Python certificate checks.
Requires a C++17 compiler (g++ by default); all temporary outputs stay in scratch.
"""
from pathlib import Path
from collections import Counter
import hashlib,itertools,json,os,re,subprocess,tempfile
P=Path(__file__).resolve().parent

def require(ok,message):
 if not ok:raise RuntimeError(message)

def hp(rows):
 n=len(rows);D=[[0]*n for _ in range(1<<n)];incoming=[sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n)]
 for v in range(n):D[1<<v][v]=1
 for mask in range(1,1<<n):
  if not(mask&(mask-1)):continue
  last=mask
  while last:
   bit=last&-last;last-=bit;v=bit.bit_length()-1;prev=mask^bit;options=incoming[v]&prev;s=0
   while options:
    b=options&-options;options-=b;s+=D[prev][b.bit_length()-1]
   D[mask][v]=s
 return sum(D[-1])

def strong(rows):
 reach=[r|1<<i for i,r in enumerate(rows)]
 for k in range(len(rows)):
  for i in range(len(rows)):
   if reach[i]>>k&1:reach[i]|=reach[k]
 return all(x==(1<<len(rows))-1 for x in reach)

def canonical(rows):
 n=len(rows);ds=[r.bit_count() for r in rows];groups=[[i for i in range(n) if ds[i]==d] for d in sorted(set(ds))]
 best=None
 for chunks in itertools.product(*(itertools.permutations(g) for g in groups)):
  p=tuple(v for g in chunks for v in g)
  out=tuple(sum(1<<j for j,w in enumerate(p) if rows[v]>>w&1) for v in p)
  if best is None or out<best:best=out
 return best

def rank(rows,r):
 O=[v for v in range(len(rows)) if rows[r]>>v&1];I=[v for v in range(len(rows)) if v!=r and rows[v]>>r&1];partner={}
 def augment(v,seen):
  for w in I:
   if rows[v]>>w&1 and w not in seen:
    seen.add(w)
    if w not in partner or augment(partner[w],seen):partner[w]=v;return True
  return False
 return sum(augment(v,set()) for v in O)

def clone(rows,r,m):
 labels=[v for v in range(len(rows)) if v!=r]+[r]*m;n=len(labels)
 return [sum(1<<j for j,w in enumerate(labels) if i!=j and ((i<j) if v==w else bool(rows[v]>>w&1))) for i,v in enumerate(labels)]

def main():
 cpp=P/'verify_minimum_ternary.cpp'
 with tempfile.TemporaryDirectory(prefix='ternary-min-',dir=P) as tmp:
  exe=Path(tmp)/'search';bankfile=Path(tmp)/'bank.jsonl'
  mode=os.environ.get('TERNARY_CPP_MODE','optimized');flags=['-O3'] if mode=='optimized' else ['-O0','-D_GLIBCXX_ASSERTIONS']
  subprocess.run([os.environ.get('CXX','g++'),'-std=c++17',*flags,str(cpp),'-o',str(exe)],check=True)
  result=subprocess.run([str(exe),str(bankfile)],check=True,text=True,capture_output=True)
  require(result.stdout.strip()=='DONE tested=1624 hits=0',('complete search outcome',result.stdout))
  entries=[json.loads(line) for line in bankfile.read_text().splitlines()]
  expected={3:1,4:1,5:6,6:35,7:351,8:659,9:347,10:66};require(dict(Counter(z['n'] for z in entries))==expected,'complete retained class counts')
  generation=[]
  for line in result.stderr.splitlines():
   m=re.fullmatch(r'order (\d+) candidates (\d+) accepted (\d+) classes (\d+)',line)
   if m:generation.append(dict(zip(['order','extension_candidates','accepted_before_isomorphism','classes'],map(int,m.groups()))))
  bank={n:{} for n in expected};checked=0;first_equal=0
  for z in entries:
   n=z['n'];rows=z['rows'];H=z['H'];require(len(rows)==n and all(0<=v<1<<n and not(v>>i&1) for i,v in enumerate(rows)),'valid rows')
   require(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(n) for j in range(i+1,n)),'tournament')
   key=tuple(rows);require(key not in bank[n] and canonical(rows)==key,'distinct exact canonical representatives')
   require(strong(rows) and hp(rows)==H<=173,'independent retained graph counts and strongness');bank[n][key]=H
   if H%6!=5:continue
   for r in range(n):
    if rank(rows,r)!=2:continue
    checked+=1;H2=hp(clone(rows,r,2))
    if H2!=3*H+2:continue
    first_equal+=1;H3=hp(clone(rows,r,3));require(H3!=9*H+8,'no smaller pure ternary clock')
  require(checked==1624,'independent candidate count')
  # Independent all-labelled coverage of the small complete strong banks.
  labelled=[]
  for n in range(3,7):
   pairs=list(itertools.combinations(range(n),2));seen=set();count=0
   for mask in range(1<<len(pairs)):
    rows=[0]*n
    for k,(i,j) in enumerate(pairs):
     a,b=(i,j) if mask>>k&1 else (j,i);rows[a]|=1<<b
    if strong(rows):seen.add(canonical(rows));count+=1
   require(seen==set(bank[n]),'independent complete labelled strong universe')
   labelled.append({'order':n,'labelled_strong_graphs':count,'classes':len(seen)})
  witness=[2,12,145,149,99,143,175,19];r=7
  require(strong(witness) and rank(witness,r)==2,'attaining source strong and matching rank')
  values=[hp(clone(witness,r,m)) for m in range(1,6)]
  require(values==[60*3**m-1 for m in range(1,6)],'actual amplitude60 family controls')
  a1=values[0];a2=(values[1]-a1)//2;a3=(values[2]-a1-6*a2)//6
  require((a1,a2,a3)==(179,180,60),'exact attaining run coefficients')
  require(5**10>173**3,'global strong order cutoff')
  receipt={'status':'PASS','claim':'60 is the globally least positive amplitude for a transitive-root-clone family A*3^m-1, even if the identity is required only eventually.','strong_H_cap':173,'retained_classes_by_order':expected,'generation':generation,'retained_bank_sha256':hashlib.sha256(bankfile.read_bytes()).hexdigest(),'independent_retained_graph_checks':len(entries),'independent_candidate_root_checks':checked,'first_clone_equation_survivors':first_equal,'independent_labelled_controls':labelled,'attaining_rows':witness,'attaining_root':r,'attaining_run_coefficients':[a1,a2,a3],'attaining_actual_clone_counts':values,'scope':'Transitive cloning of one fixed root, not every possible graph-changing compiler.'}
 print(json.dumps(receipt,indent=2,sort_keys=True))
if __name__=='__main__':main()
