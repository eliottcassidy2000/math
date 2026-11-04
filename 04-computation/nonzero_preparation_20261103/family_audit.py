from pathlib import Path
from itertools import permutations, product
import ast,json
base=Path(__file__).resolve().parent
tree=ast.parse((base/'surgery_ready.py').read_text())
bank=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='BANK' for t in n.targets))
gates=0
def need(ok,msg):
 global gates
 gates+=1
 if not ok:raise RuntimeError(msg)
def data(rows,ports):
 n=len(rows);full=(1<<n)-1;H=[0]*(full+1);E=[0]*(full+1);S=[0]*(full+1);H[0]=1
 for l in range(1,n+1):
  for p in permutations(range(n),l):
   if all(rows[a]>>b&1 for a,b in zip(p,p[1:])):
    m=sum(1<<a for a in p);H[m]+=1;E[m]+=p[-1] in ports;S[m]+=p[0] in ports
 return (H[full],sum(H[m]*H[full^m] for m in range(full+1)),E[full],S[full],sum(E[m]*H[full^m] for m in range(full+1)),sum(S[m]*H[full^m] for m in range(full+1)),sum(E[m]*S[full^m] for m in range(full+1)))
def simple_word(d,w):
 h,p2,e,b,F,G,C=d
 L=1+sum(v*(1<<i) for i,v in enumerate(w))
 B=1
 for i,v in enumerate(w):
  if v:
   standard=sum(w[j]*(1<<(i-1-j)) for j in range(i))
   B+=(1<<i)-standard-(w[i-1] if i else 1)
 D=1 if not w else 1-w[0]
 return (G+p2)*L+h*B-e*D-C-F-b

data_bank=[data(rows,ports) for gap,rows,ports in bank]
for (gap,rows,ports),d in zip(bank,data_bank):need(d[3]==gap,'actual entry gap from permutation count')
receipt=json.loads((base/'family.json').read_text())
recovered=[]
for query in receipt['requested']+[receipt['positive_control']]:
 A,B0=query['A'],query['B'];hits=[];prefixes=0;queries=0
 for idx,((gap,rows,ports),d) in enumerate(zip(bank,data_bank)):
  h,p2,e,b,F,G,C=d;c=G+p2;n=0
  while c*(1<<n)<=A:
   for p in product((0,1),repeat=n):
    prefixes+=1;T=(1<<n)-sum(p[j]*(1<<(n-1-j)) for j in range(n))
    for bit in (0,1):
     factor=c*(1<<n)+h*(T-bit)
     if A%factor:continue
     v=A//factor-bit
     if v<1-bit:continue
     q=tuple((v>>j)&1 for j in range(v.bit_length()))
     constant=simple_word(d,p+(bit,)+q)-2*A
     for k in (2,4,7):need(simple_word(d,p+(bit,)*k+q)==A*(1<<k)+constant,'independent one-run coefficient')
     queries+=1
     if constant==B0:hits.append((idx,gap,p,bit,q))
   n+=1
 want=[(r['bank_index'],r['gap'],tuple(r['prefix']),r['bit'],tuple(r['suffix'])) for r in query['hits']]
 need(hits==want,'all unbounded inverse presentations')
 need(prefixes==query['prefixes'] and queries==query['divisor_queries'],'all complete prefix/divisor cardinalities')
 recovered.append((A,B0,len(hits),prefixes,queries))
print(json.dumps({'gates':gates,'literal_bank':len(bank),'recovered':recovered},sort_keys=True))
