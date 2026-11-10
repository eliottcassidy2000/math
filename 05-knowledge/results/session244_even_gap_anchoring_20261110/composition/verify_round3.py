import importlib.util,json
from pathlib import Path
from itertools import product
s=importlib.util.spec_from_file_location('verify',Path(__file__).with_name('verify.py'));v=importlib.util.module_from_spec(s);s.loader.exec_module(v);p=v.p

def word(N):
 out=[]
 while N>2:
  if N%2:out.append('Z');N=(N+1)//2
  else:out.append('O');N//=2
 return ''.join(reversed(out))

def main():
 words=['OJJJJOZ','OJOZJOZ']
 mats=[[[17,0,40,0],[3,0,6,0],[21,0,50,0],[3,0,6,0]],[[17,0,42,0],[3,0,6,0],[21,0,50,0],[3,0,6,0]]]
 p.need([v.matrix(w) for w in words]==mats,'pure difference full matrices')
 rows=[];tailtests=0
 for N in range(1,33):
  a,r=([0],0) if N==1 else v.literal([2,0],0,word(N))
  z=p.state(a,r);p.need(z[2]==N and z[0]<=N,'response selector')
  states=[];graphs=[]
  for w in words:
   b,s=v.literal(a,r,w);q=p.state(b,s)
   p.need(q==v.act(z,w) and p.strong(b),'literal pure difference compiler')
   p.need(b[s].bit_count()==2,'two exit output root')
   states.append(q);graphs.append((b,s))
  p.need(tuple(y-x for x,y in zip(*states))==(2*N,0,0,0),'pure count difference')
  p.need(states[1][0]<=59*N and len(b)<=(8 if N==1 else (N-1).bit_length()+8),'linear cost and size')
  for length in range(8):
   for tail in product('OZ',repeat=length):
    w=''.join(tail);lo,hi=[v.act(q,w) for q in states]
    p.need(tuple(y-x for x,y in zip(lo,hi))==(2*N,0,0,0),'all tested tails preserve difference');tailtests+=1
  if N<=3:
   for w in ('Z','OZ','ZZ','OOZ','ZOZ'):
    out=[]
    for b,s in graphs:
     c,t=v.literal(b,s,w);q=p.state(c,t)
     p.need(p.strong(c),'tail strong');out.append(q)
    p.need(tuple(y-x for x,y in zip(*out))==(2*N,0,0,0),'literal continued pair')
  rows.append(dict(N=N,source=z,lower=states[0],upper=states[1],order=len(b)))
 p.need(rows[1]['lower'][0]==97 and rows[1]['upper'][0]==101,'prime example')
 print(json.dumps(dict(status='PASS',full_matrices=dict(zip(words,mats)),literal_pairs=rows,tail_algebra_cases=tailtests,literal_continued_pairs=15,scope='Every even gap with equal full rooted sidecar and gap-preserving O/Z tails; arbitrary placement remains open.'),sort_keys=True,indent=2))
if __name__=='__main__':main()
