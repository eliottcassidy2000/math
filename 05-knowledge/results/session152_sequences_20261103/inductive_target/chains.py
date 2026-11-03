from pathlib import Path
exec((Path(__file__).resolve().parent.parent / 'cut_realizability' / 'check.py').read_text().split('for n in range(2,5):')[0])
def ear(a,m):
 n=len(a); return [r+[1-(m>>i&1)] for i,r in enumerate(a)]+[[(m>>i)&1 for i in range(n)]+[0]]
for code,masks in [(8,[5,10]),(759,[1,2])]:
 a=adjacency(5,code)
 for m in masks:
  b=ear(a,m);h=hp(b); targets=(3*h,3*h+2,3*h+4)
  got={t:[] for t in targets}
  for k in range(1,63):
   q=hp(ear(b,k))
   if q in got: got[q].append(k)
  print(code,m,h,got,flush=True)
a=ear(adjacency(5,759),2)
for first in [50,24,28]:
 b=ear(a,first);h=hp(b);targets=(3*h,3*h+2,3*h+4);got={t:[] for t in targets}
 for k in range(1,127):
  q=hp(ear(b,k))
  if q in got:got[q].append(k)
 print('759,2,',first,'H',h,got,flush=True)
b=ear(a,30);h=hp(b);got={t:[] for t in (3*h,3*h+2,3*h+4)}
for k in range(1,127):
 q=hp(ear(b,k))
 if q in got:got[q].append(k)
print('759,2,30',h,got)
