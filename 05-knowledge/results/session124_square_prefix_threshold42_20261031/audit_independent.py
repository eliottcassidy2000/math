#!/usr/bin/env python3
"""Standalone independent threshold42 auditor. No producer module imports."""
import argparse,copy,gzip,hashlib,json,math,struct
from pathlib import Path
checks=0

def need(ok,reason):
 global checks
 checks+=1
 if not ok: raise ValueError(reason)

def square(x):return x>0 and math.isqrt(x)**2==x

def choose(p,d,c):return int(49+(-1 if p==0 else 1)*c<=d)

def end(p,f='B'):
 return (10 if p==0 else 8) if f=='Aplus' else (8 if p==0 else 10)

def check_pair_table(w,recursive=False):
 p,d,f=w['parity'],w['residue'],w['family'];need(p in [0,1] and 24<=d<=72,'dispatch bounds')
 cuts=w.get('split_prefix_blocks',[]);need(len(cuts)==len(set(cuts)) and set(cuts)<=set(range(-24,25)),'cuts unique in range')
 if recursive:need(4 not in cuts and f=='B','protected block')
 need(len(w['words'])==2,'two words'); profiles=[]
 for b,word in enumerate(w['words']):
  # symbolic lengths a*n+b, with parity obtained solely from p.
  a=0;k=0;last=None;profile={};occ=set();headshifts={};shifts={}
  expected={('single',x) for x in range(1,25)}|{(kind,c) for c in range(-24,25) for kind in (('head','tail') if c in cuts else ('block',))}
  if recursive:need(word[:5]==[['single',1],['single',3],['single',6],['single',19],['block',4,choose(p,d+b,4),1]],'recursive exact prefix')
  else:
   if f=='B':need(word[:2]==[['single',1],['single',3]],'B prefix')
   else:need(word[0]==['block',-20,choose(p,d+b,-20),1],'A full prefix block')
  for atom in word:
   kind,x=atom[:2];need((kind,x) in expected and (kind,x) not in occ,'exact atom ownership');occ.add((kind,x));s=(a*p+k+1)%2
   if kind in ['single','head']:
    need(len(atom)==2,'singleton arity');lo=hi=x if kind=='single' else 49-x;da,dk=0,1
    if kind=='single':profile[('single',x)]=s
    else:headshifts[x]=(s-1)%2
   else:
    need(len(atom)==4,'block arity');_,c,i,o=atom
    need(i==choose(p,d+b,c) and o in [-1,1],'forced input and orientation')
    L=(p+i)%2
    lo=(147+c if kind=='tail' else 49-c);hi=49*end(L)+(-1 if L else 1)*c
    if o==-1:lo,hi=hi,lo
    da,dk=1,i-int(kind=='tail')
    shift=(s-(2 if kind=='tail' else 1))%2 if o==1 else (s+L)%2
    shifts[c]=shift;profile[('offset',c)]=shift
   if last is None:
    need(lo==(1 if f=='B' else 69),'first endpoint')
   else:need(square(last+lo),'connector square')
   last=hi;a+=da;k+=dk
  need(occ==expected and a==49 and k==d+b,'coverage and symbolic output length')
  need(last==end((p+d+b)%2,f),'terminal rule')
  for c in cuts:need(headshifts[c]==shifts[c],'head-tail common parity shift')
  if f!='B':
   # Protected vertices 2,7 are actual singleton atoms, so parity and adjacency are exact.
   j=word.index(['single',2]);h=word.index(['single',7]);need(abs(j-h)==1,'A protected edge')
   need(profile['single',2]==1 and profile['single',7]==0,'A protected parity')
   need(set(cuts)<={0,-2},'A excluded cut offsets')
  profiles.append(profile)
 need(profiles[0]==profiles[1],'common parity profile')
 return len(cuts)

def prefix(l):
 out=[]
 for t in range(l+1):
  q=49**t;c=(q-1)//12;out += [q-c,3*q+c,6*q-c,19*q+c][:2 if t==l else 4]
 return out

def validate(pair,pref,f='B'):
 n=len(pair[0]);need(len(pair[1])==n+1,'consecutive concrete lengths');par=bytearray(n+2)
 for b,path in enumerate(pair):
  need(path[:len(pref)]==pref and path[-1]==end(len(path)%2,f),'concrete endpoints/prefix');seen=bytearray(len(path)+1)
  prev=None
  for j,v in enumerate(path):
   need(1<=v<=len(path) and not seen[v],'concrete coverage');seen[v]=1
   if b==0:par[v]=j%2
   elif v<=n:need(par[v]==j%2,'concrete parity')
   if prev is not None:need(square(prev+v),'concrete square')
   prev=v
  need(sum(seen)==len(path),'concrete coverage count')

def expand(pair,words):
 result=[]
 for word in words:
  out=[]
  for a in word:
   if a[0]=='single':out.append(a[1]);continue
   if a[0]=='head':out.append(49-a[1]);continue
   typ,c,i,o=a;P=pair[i]
   js=range(int(typ=='tail'),len(P)) if o==1 else range(len(P)-1,int(typ=='tail')-1,-1)
   out.extend(49*P[j]+(c if j%2 else -c) for j in js)
  result.append(out)
 return result

def check_cycle_table(w):
 p,d=w['parity'],w['residue'];need(w['family']=='B-final-cycle-prefix4' and p in [0,1] and 24<=d<=72,'terminal dispatch')
 cuts=w['split_prefix_blocks'];need(len(cuts)==len(set(cuts)) and set(cuts)<=set(range(-24,25)) and 4 not in cuts,'terminal cuts')
 expected={('single',v) for v in range(1,25)}|{(typ,c) for c in range(-24,25) for typ in (('head','tail') if c in cuts else ('block',))}
 word=w['word'];need(word[:5]==[['single',1],['single',3],['single',6],['single',19],['block',4,choose(p,d,4),1]],'terminal exact prefix')
 seen=set();length_n=length_const=0;previous=None
 for atom in word:
  kind,c=atom[:2];key=(kind,c);need(key in expected and key not in seen,'terminal atom ownership');seen.add(key)
  if kind in ['single','head']:
   need(len(atom)==2,'terminal singleton format');left=right=c if kind=='single' else 49-c;length_const+=1
  else:
   need(len(atom)==4,'terminal block format');_,c,i,o=atom;need(i==choose(p,d,c) and o in [-1,1],'terminal chosen member')
   L=(p+i)%2;left=49-c if kind=='block' else 147+c;right=49*end(L)+(c if L==0 else -c)
   if o==-1:left,right=right,left
   length_n+=1;length_const+=i-int(kind=='tail')
  if previous is None:need(left==1,'terminal first')
  else:need(square(previous+left),'terminal connector')
  previous=right
 need(seen==expected and (length_n,length_const)==(49,d),'terminal interval partition')
 need(word[-1]==['single',8 if (p+d)%2==0 else 15],'terminal endpoint')
 need(square(previous+1),'terminal closing square')

def digest(data):return hashlib.sha256(data).hexdigest()

def read(path):return json.loads(path.read_bytes())

def threshold(l):return (85*49**l-1)//2
def depth(n):
 l=0
 while threshold(l+1)<=n:l+=1
 return l

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 for name in ['extension-dir','core-dir','base-dir']:ap.add_argument('--'+name,type=Path,required=True)
 args=ap.parse_args();ext=args.extension_dir;core=args.core_dir;base=args.base_dir
 dep=read(ext/'DEPENDENCIES.json')
 for name,h in dep['frozen_core_files'].items():need(digest((core/name).read_bytes())==h,'pinned frozen core file')
 smallraw=(ext/'SMALL_BASES.json').read_bytes();need(digest(smallraw)==dep['small_base_sha256'],'small catalog digest');small=json.loads(smallraw)
 need(set(small)=={str(n) for n in range(42,53)},'small base complete catalog')
 for n in range(42,53):need(len(small[str(n)][0])==n,'small base index');validate(small[str(n)],[1,3])
 bdep=read(core/'BASE_DEPENDENCIES.json');manifest=(base/bdep['manifest_file']).read_bytes();need(digest(manifest)==bdep['manifest_sha256'],'old base manifest hash');chunks=json.loads(manifest);need(chunks==bdep['chunks'],'old base metadata')
 cached={};base_count=11
 def decode(data,lo,n):
  at=16+sum(2*k-1 for k in range(lo,n));pair=[]
  for L in [n,n+1]:
   roots=data[at:at+L-1];need(len(roots)==L-1,'root count');path=[1]
   for root in roots:path.append(root*root-path[-1])
   at+=L-1;pair.append(path)
  return pair
 for c in chunks:
  if c['first']>2081:continue
  compressed=(base/c['file']).read_bytes();need(len(compressed)==c['bytes'] and digest(compressed)==c['sha256'],'old compressed chunk')
  raw=gzip.decompress(compressed);need(len(raw)==c['decoded_bytes'] and digest(raw)==c['decoded_sha256'],'old decoded chunk')
  lo,hi=struct.unpack_from('<II',raw,8);need(raw[:8]==b'SSBAS1\0\0' and (lo,hi)==(c['first'],c['last']),'old chunk header');cached[c['file']]=raw
  for n in range(lo,min(hi,2081)+1):validate(decode(raw,lo,n),[1,3]);base_count+=1
 need(base_count==2040,'complete42 through2081 base supply')
 pairs={};cycles={}
 for p in [0,1]:
  for d in range(24,73):
   w=read(core/f'paired/B_p{p}_d{d}.json');need((w['parity'],w['residue'])==(p,d),'paired dispatch');check_pair_table(w,True);pairs[p,d]=w
   w=read(core/f'cycles/cycle_p{p}_d{d}.json');need((w['parity'],w['residue'])==(p,d),'cycle dispatch');check_cycle_table(w);cycles[p,d]=w
 def base_pair(n):
  need(42<=n<=2081,'independent base cutoff2081')
  if n<53:return small[str(n)]
  c=next(c for c in chunks if c['first']<=n<=c['last']);return decode(cached[c['file']],c['first'],n)
 def build(n):
  need(n>=42,'independent input lower bound')
  if n<=2081:pair=base_pair(n)
  else:
   m,d=divmod(n-24,49);d+=24;need(42<=m<n and 24<=d<=72,'well-founded quotient');pair=expand(build(m),pairs[m%2,d]['words'])
  need(len(pair[0])==n,'constructed B index');validate(pair,prefix(depth(n)));return pair
 def stream(pair,word):
  for atom in word:
   if atom[0]=='single':yield atom[1]
   elif atom[0]=='head':yield 49-atom[1]
   else:
    kind,c,i,o=atom;P=pair[i];start=int(kind=='tail');js=range(start,len(P)) if o==1 else range(len(P)-1,start-1,-1)
    for j in js:yield 49*P[j]+(c if j%2 else -c)
 controls=[]
 for N in [2082,2083,2523,2572,102041,102042,102043,5000082]:
  m,d=divmod(N-24,49);d+=24;pair=build(m);P=prefix(depth(N));initial=[];previous=None;count=0;seen=bytearray(N+1);h=hashlib.sha256()
  for v in stream(pair,cycles[m%2,d]['word']):
   need(1<=v<=N and not seen[v],'stream cycle ownership');seen[v]=1;count+=1;h.update((str(v)+'\n').encode())
   if len(initial)<len(P):initial.append(v)
   if previous is not None:need(square(previous+v),'stream square edge')
   previous=v
  need(count==N and sum(seen)==N,'stream exact interval');need(initial==P,'stream maximum unscaled prefix');need(previous==(8 if N%2==0 else 15) and square(previous+initial[0]),'stream closure')
  controls.append({'N':N,'depth':depth(N),'prefix':P,'cycle_decimal_lines_sha256':h.hexdigest()})
 # A cutoff regression must be observable: the certified constructor really recurses here.
 boundary=build(2082);need(boundary[0][:6]==prefix(1),'new cutoff prefix regression')
 rejected=[]
 for label,change in [('repeated_vertex',lambda x:x[0].__setitem__(2,x[0][1])),('changed_endpoint',lambda x:x[0].__setitem__(-1,8)),('wrong_order',lambda x:x[0].reverse())]:
  q=copy.deepcopy(small['51']);change(q)
  try:validate(q,[1,3])
  except ValueError as e:rejected.append([label,str(e)])
  else:raise ValueError('small base mutation accepted')
 print(json.dumps({'status':'passed','checks':checks,'base_pairs':base_count,'small_base_pairs':11,'base_interval':[42,2081],'recursive_cutoff':2082,'old_chunks_used':len(cached),'paired_states':98,'direct_states':98,'producer_imports':0,'controls':controls,'small_base_mutations':rejected,'small_catalog_sha256':digest(smallraw)},indent=2,sort_keys=True))
if __name__=='__main__':main()
