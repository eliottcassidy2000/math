#!/usr/bin/env python3
import json,collections,subprocess,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
D=next(p for p in (ROOT/'session40-square-sum',ROOT/'square-sum') if p.exists())

def require(b,msg):
 if not b:raise RuntimeError(msg)
def valid(p,n):return sorted(p)==list(range(1,n+1)) and all(a+b in {j*j for j in range(1,14)} for a,b in zip(p,p[1:]))
base=json.load(open(D/'receipt.json'))['rows'];birth=json.load(open(D/'birth_witnesses_25_80.json'));hostile=json.load(open(D/'COLOR_WORD_HOSTILE_80.json'));census=json.load(open(D/'COLOR_FIRST_FAILURE_CENSUS.json'))
require(set(map(int,birth))==set(range(25,81)),'birth coverage')
for n,p in birth.items():require(valid(p,int(n)),'birth witness')
status={r['n']:bool(r['count']) for r in base};status.update({int(n):True for n in birth})
def color(a,b):
 n=max(a,b)
 if status[n]:return 1
 degree=sum(n<s<2*n for s in (j*j for j in range(1,14)))
 return 2 if degree==2 else 3
rows=[]
for r in base[24:]:
 words=set();hist=set()
 for p in r['paths']:
  if p[0]>p[-1]:p=p[::-1]
  w=tuple(color(a,b) for a,b in zip(p,p[1:]));words.add((p[0],p[-1],w));hist.add((p[0],p[-1],tuple(w.count(c) for c in(1,2,3))))
 require(len(words)==len(r['paths']),'word collision before33');rows.append({'n':r['n'],'paths':len(r['paths']),'words':len(words),'histograms':len(hist)})
for p in hostile['paths']:
 require(valid(p,80),'order80path');require((p[0],p[-1])==(1,2),'order80ends');require(all(color(a,b)==1 for a,b in zip(p,p[1:])),'order80color')
require(hostile['paths'][0]!=hostile['paths'][1],'distinct hostile')
with tempfile.TemporaryDirectory() as temp:
 binary=Path(temp)/'check';subprocess.run(['g++','-O3','-std=c++17',str(HERE/'check_color_words_cpp.cpp'),'-o',str(binary)],check=True)
 output=subprocess.run([str(binary)],check=True,text=True,capture_output=True).stdout
(HERE/'color-words-independent-cpp.txt').write_text(output)
found=[]
for line in output.splitlines():
 if ',' in line:continue
 n,h,w,g,maxf,nodes=map(int,line.split());exp=next(r for r in census if r['n']==n)
 require((h,w,g,maxf)==(exp['path_count'],exp['word_labels'],exp['histogram_labels'],exp['max_word_fibre']),'C++ independent census')
 found.append({'n':n,'paths':h,'words':w,'histograms':g,'largest_word_fibre':maxf,'partial_paths_explored':nodes})
require([r['n'] for r in found]==[33,34,35,36],'new order universe')
lines=output.splitlines();pair=[tuple(map(int,x.strip(',').split(','))) for x in lines if ',' in x]
require(len(pair)==2 and pair[0]!=pair[1],'new hostile pair')
for p in pair:require(valid(p,36),'order36path')
require((pair[0][0],pair[0][-1])==(pair[1][0],pair[1][-1])==(6,32),'order36ends')
require([color(a,b) for a,b in zip(pair[0],pair[0][1:])]==[color(a,b) for a,b in zip(pair[1],pair[1][1:])],'order36word')
out={'status':'PASS','birth_witnesses':56,'through32':rows,'independent_cpp_33_36':found,'order36_pair':pair,'order80_monochromatic_paths':2}
(HERE/'color-extension-audit-receipt.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k not in('through32','order36_pair')})
