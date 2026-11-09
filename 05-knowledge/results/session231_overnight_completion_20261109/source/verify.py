import json,subprocess,random
from pathlib import Path
Y=[tuple(map(int,s.split())) for s in Path('bank.txt').read_text().splitlines()]
def g(v):
 z=1
 while v:
  if v%2==0:z+=v//2+(v//2)%2
  v//=2
 return z
def check(p,msg):
 if not p:raise RuntimeError(msg)
check(len(Y)==2732 and len(set(Y))==2732,'bank size')
check(all(a>0 and a%2 and c>=a>=b>=0 for a,b,c in Y),'physical inequality')
lim=2048
image={}
for a,b,c in Y:
 for v in range(lim//c+1):
  n=a*g(v)-b*(1-v%2)+c*(1+v)
  if n<=lim:image[n]=(a,b,c,v)
outs=subprocess.check_output(['./inverse','bank.txt']+list(map(str,range(1,lim,2))),text=True).splitlines()
check(len(outs)==len(range(1,lim,2)), 'output row count')
check([int(line.split()[0]) for line in outs]==list(range(1,lim,2)), 'output row order')
for line in outs:
 p=line.split();n=int(p[0]);check((p[1]!='NONE')==(n in image),'membership')
 if p[1]!='NONE':
  a,b,c,v=map(int,p[1:5]);check((a,b,c) in Y and a*g(v)-b*(1-v%2)+c*(1+v)==n,'witness')
cases=0
for k in range(1,9):
 M=2**k
 for r in range(M):
  rev=int(f'{r:0{k}b}'[::-1],2);C=M-1-rev;E=1-(r>=M//2)
  for u in range(101):
   check(g(M*u+r)==g(u)+C*u+E*(u%2)+g(r)-1,'suffix');cases+=1
large=[]
for s in Path('CASES.txt').read_text().splitlines():
 p=s.split()
 if p[1]!='NONE':
  n,a,b,c,v=map(int,p[:5]);check((a,b,c) in Y and a*g(v)-b*(1-v%2)+c*(1+v)==n,'large witness');large.append(n)
print(json.dumps({'status':'PASS','bank_size':len(Y),'exhaustive_targets':len(outs),'suffix_identities':cases,'small_missing_odds':[n for n in range(1,lim,2) if n not in image],'checked_case_targets':large},indent=2))
