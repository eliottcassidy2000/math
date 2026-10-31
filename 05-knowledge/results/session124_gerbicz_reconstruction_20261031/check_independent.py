#!/usr/bin/env python3
"""Independent inert-data audit. No imported producer code and no external-code execution."""
import re,sys,json,hashlib,copy
from pathlib import Path
from math import isqrt
NTEST=0
def need(v,msg):
 global NTEST
 NTEST+=1
 if not v: raise ValueError(msg)
def parse(source,name):
 # A brace-depth parser with a complete integer/comma/brace grammar, not JSON.
 starts=list(re.finditer(r'\bint\s+'+name+r'\s*\[[^;=]*=',source))
 need(len(starts)==1,'unique declaration '+name)
 pos=starts[0].end()
 token=re.compile(r'\s*([{},]|-?\d+)')
 def take():
  nonlocal pos
  m=token.match(source,pos)
  need(m is not None,'literal grammar')
  pos=m.end();return m.group(1)
 def array():
  nonlocal pos
  need(take()=='{','open brace'); result=[]
  while True:
   m=token.match(source,pos);need(m is not None,'token')
   if m.group(1)=='}':pos=m.end();return result
   if m.group(1)=='{': result.append(array())
   else:
    t=take();need(re.fullmatch(r'-?\d+',t) is not None,'integer');result.append(int(t))
   t=take()
   if t=='}':return result
   need(t==',','separator')
 need(source[pos:].lstrip().startswith('{'),'initializer')
 result=array();need(source[pos:].lstrip().startswith(';'),'initializer end');return result

def path(a,n,cycle=True):
 need(len(a)==n,'path size');seen=bytearray(n+1)
 for x in a:
  need(type(x)==int and 1<=x<=n,'label range');need(not seen[x],'duplicate label');seen[x]=1
 for i in range(n if cycle else n-1):
  s=a[i]+a[(i+1)%n];need(isqrt(s)**2==s,'square adjacency')

def audit_words(rows):
 need(len(rows)==196,'196 words');by={}
 for index,row in enumerate(rows):
  d=24+index//4;p=(index//2)%2;e=index%2
  need(row[0]==[49,d,p,e],'original dispatch order')
  need(len(row)==74,'73 atoms')
  offsets={};singles={};ends=[];total=0;shift={}
  # Representative length determines all positional parity expressions exactly.
  n=42 if p==0 else 41
  for atom in row[1:]:
   need(len(atom)==3,'atom arity');c,t,h=atom
   if t==0:
    need(h==0 and 1<=c<=24 and c not in singles,'singleton schema')
    singles[c]=(total+1)%2;total+=1;ends.append((c,c));continue
   need(t in (-1,1) and h in (0,1) and -24<=c<=24 and c not in offsets,'block schema')
   offsets[c]=h;L=n+h
   # x=n+1 is present exactly when its lifted value is below the target end.
   need(h==int(49*(n+1)+((-1)**n)*c<=49*n+d+e),'extension cutoff')
   lo=49+c;hi=49*(3 if L%2 else 8)+((-1)**(L+1))*c
   ends.append((lo,hi) if t==1 else (hi,lo))
   # If local index is j, global index is total+j, or total+L+1-j.
   shift[c]=(total if t==1 else total+L+1)%2;total+=L
  need(set(offsets)==set(range(-24,25)) and set(singles)==set(range(1,25)),'complete atom partition')
  need(total==49*n+d+e,'exact size')
  need(ends[0][0]==1 and ends[-1][1]==(3 if total%2 else 8),'output endpoints')
  for i in range(len(ends)):
   s=ends[i][1]+ends[(i+1)%len(ends)][0];need(isqrt(s)**2==s,'connector')
  # Direct set partition at both possible shared parity values, separately for each x.
  for sign in (-1,1):need({sign*c for c in offsets}==set(range(-24,25)),'full layer')
  extra={49+((-1)**n)*c for c,h in offsets.items() if h}
  need(extra==set(range(25,d+e+1)),'last layer prefix')
  by[d,p,e]=(singles,shift)
 for d in range(24,73):
  for p in range(2):need(by[d,p,0]==by[d,p,1],'pair parity preserved')
 return by

def main():
 raw=Path(sys.argv[1]).read_bytes()
 need(hashlib.sha256(raw).hexdigest()=='602c1919ee51101519a2f19f2c1817407006824f6f09c85a3b772a818a0b736e','pinned source SHA256')
 s=raw.decode('ascii')
 small=parse(s,'smallseq');words=parse(s,'data');bases=parse(s,'S')
 need(len(small)==41 and len(bases)==3984,'table lengths')
 for n,a in enumerate(small):
  if n in {15,16,17,23}|set(range(25,41)):path(a,n,n>=32)
  else:need(a==[],'empty small row')
 for n in range(41,2033):
  pair=bases[2*(n-41):2*(n-41)+2];profiles=[]
  for j,a in enumerate(pair):
   path(a,n+j);need(a[0]==1 and a[-1]==(3 if (n+j)%2 else 8),'nice endpoint')
   profiles.append({x:i%2 for i,x in enumerate(a,1)})
  need(all(profiles[0][x]==profiles[1][x] for x in range(1,n+1)),'base shared parity')
 audit_words(words)
 def lift(n,d,e,pair):
  out=[]
  for c,t,h in words[4*(d-24)+2*(n%2)+e][1:]:
   if t==0:out.append(c);continue
   block=[49*x+(c if j%2 else -c) for j,x in enumerate(pair[h],1)]
   out.extend(block if t==1 else block[::-1])
  return out
 concrete=0
 for n in (41,42):
  pair=bases[2*(n-41):2*(n-41)+2]
  for d in range(24,73):
   outputs=[lift(n,d,e,pair) for e in range(2)]
   for e,a in enumerate(outputs):path(a,49*n+d+e);concrete+=1
   pa={x:j%2 for j,x in enumerate(outputs[0],1)}
   pb={x:j%2 for j,x in enumerate(outputs[1],1)}
   need(all(pa[x]==pb[x] for x in pa),'expanded pair parity')
 def construct(n):
  if n<=2032:return bases[2*(n-41):2*(n-41)+2]
  m=(n-24)//49;d=n-49*m;pair=construct(m)
  return [lift(m,d,e,pair) for e in range(2)]
 for n in (2033,2079,2080,100000,1000000):
  pair=construct(n)
  for e,a in enumerate(pair):path(a,n+e);concrete+=1
 # Verify boundary arithmetic and every residue at both recursive parent parities.
 for n in (41,42,2032,2033):
  for d in range(24,73):
   target=49*n+d;parent=(target-24)//49;residue=target-49*parent
   need((parent,residue)==(n,d),'unique quotient dispatch')
   need(41<=parent<target,'well founded parent')
 mutations=[]
 for label,change in [('missing_word',lambda w:w.pop()),('dispatch_swap',lambda w:w.reverse()),('bad_offset',lambda w:w[0][2].__setitem__(0,25)),('bad_member',lambda w:w[0][2].__setitem__(2,1)),('bad_singleton',lambda w:w[0][1].__setitem__(0,2))]:
  w=copy.deepcopy(words);change(w)
  try:audit_words(w)
  except ValueError:mutations.append(label)
  else:raise ValueError('mutation accepted '+label)
 print(json.dumps(dict(status='PASS',source_sha256=hashlib.sha256(raw).hexdigest(),base_pairs=1992,universal_words=196,concrete_paths=concrete,checks=NTEST,mutations_rejected=mutations,external_code_executed=False),indent=2))
if __name__=='__main__':main()
