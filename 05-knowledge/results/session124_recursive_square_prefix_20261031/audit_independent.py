#!/usr/bin/env python3
"""Standalone independent minimal-core auditor; Python standard library only.

Consumes only certified base chunks and new paired/direct JSON certificates.
No producer module, Git object database, old transition table, or network access.
"""
import argparse, copy, gzip, hashlib, json, math, struct
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

def check_cycle(P,prefix):
 N=len(P);seen=bytearray(N+1);need(P[:len(prefix)]==prefix,'direct unscaled prefix')
 for j,v in enumerate(P):
  need(1<=v<=N and not seen[v],'direct cycle coverage');seen[v]=1;need(square(P[j-1]+v),'direct cycle square including closure')
 need(sum(seen)==N,'direct cycle complete')

def digest(data):return hashlib.sha256(data).hexdigest()
def read(path):return json.loads(path.read_bytes())

def load_bases(base,dep):
 need(dep['base_indices']==[53,2620],'base dependency interval')
 mf=base/dep['manifest_file'];raw=mf.read_bytes();need(digest(raw)==dep['manifest_sha256'],'base manifest digest')
 chunks=json.loads(raw);need(chunks==dep['chunks'],'base chunk metadata exact equality')
 retained={};indices=[]
 for c in chunks:
  compressed=(base/c['file']).read_bytes();need(len(compressed)==c['bytes'] and digest(compressed)==c['sha256'],'compressed base metadata')
  data=gzip.decompress(compressed);need(len(data)==c['decoded_bytes'] and digest(data)==c['decoded_sha256'],'decoded base metadata')
  need(data[:8]==b'SSBAS1\0\0','base magic');lo,hi=struct.unpack_from('<II',data,8);need((lo,hi)==(c['first'],c['last']),'base range');at=16
  for n in range(lo,hi+1):
   pair=[]
   for L in [n,n+1]:
    roots=data[at:at+L-1];need(len(roots)==L-1,'base root count');P=[1]
    for root in roots:P.append(root*root-P[-1])
    at+=L-1;pair.append(P)
   validate(pair,[1,3]);indices.append(n)
   if n in [53,54]:retained[n%2]=pair
  need(at==len(data),'base trailing bytes')
 need(indices==list(range(53,2621)),'all base indices exactly once')
 return retained

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--package-dir',type=Path,required=True)
 parser.add_argument('--base-dir',type=Path,required=True)
 args=parser.parse_args();package=args.package_dir;base=args.base_dir
 manifest=read(package/'MANIFEST.json');dep=read(package/'BASE_DEPENDENCIES.json')
 need(digest((package/'BASE_DEPENDENCIES.json').read_bytes())==manifest['proof_files']['BASE_DEPENDENCIES.json'],'dependency manifest package digest')
 pairs={};cycles={};file_hashes={}
 need({p.name for p in (package/'paired').glob('B_p*_d*.json')}=={f'B_p{p}_d{d}.json' for p in [0,1] for d in range(24,73)},'physical paired state coverage')
 need({p.name for p in (package/'cycles').glob('cycle_p*_d*.json')}=={f'cycle_p{p}_d{d}.json' for p in [0,1] for d in range(24,73)},'physical direct state coverage')
 for p in [0,1]:
  for d in range(24,73):
   for folder,name,table,checker in [('paired',f'B_p{p}_d{d}.json',pairs,lambda x:check_pair_table(x,True)),('cycles',f'cycle_p{p}_d{d}.json',cycles,check_cycle_table)]:
    key=folder+'/'+name;raw=(package/key).read_bytes();h=digest(raw);need(h==manifest['proof_files'][key],'frozen certificate digest');file_hashes[key]=h
    w=json.loads(raw);need((w['parity'],w['residue'])==(p,d),'filename state matches');checker(w);table[p,d]=w
 bases=load_bases(base,dep);lifted={}
 for p,input_pair in bases.items():
  output=expand(input_pair,pairs[p,24]['words']);validate(output,prefix(1));lifted[len(output[0])%2]=output
 for p,d in pairs:
  result=expand(lifted[p],pairs[p,d]['words']);need(len(result[0])==49*len(lifted[p][0])+d,'paired expansion length');validate(result,prefix(2))
 for depth,inputs in [(1,bases),(2,lifted)]:
  for (p,d),w in cycles.items():
   P=expand(inputs[p],[w['word']])[0];need(len(P)==49*len(inputs[p][0])+d,'direct expansion length');check_cycle(P,prefix(depth))
 mutations=[]
 for typ,checker,w,actions in [
  ('paired',lambda x:check_pair_table(x,True),pairs[0,24],
   [('protected_cut',lambda x:x['split_prefix_blocks'].append(4)),('wrong_input',lambda x:x['words'][0][4].__setitem__(2,1-x['words'][0][4][2])),('missing_atom',lambda x:x['words'][0].pop()),('reversed_prefix',lambda x:x['words'][0][4].__setitem__(3,-1))]),
  ('direct',check_cycle_table,cycles[0,24],
   [('protected_cut',lambda x:x['split_prefix_blocks'].append(4)),('wrong_input',lambda x:x['word'][4].__setitem__(2,1-x['word'][4][2])),('missing_atom',lambda x:x['word'].pop()),('reversed_prefix',lambda x:x['word'][4].__setitem__(3,-1))])]:
  for name,action in actions:
   changed=copy.deepcopy(w);action(changed)
   try:checker(changed)
   except ValueError as e:mutations.append([typ,name,str(e)])
   else:raise ValueError('hostile mutation accepted '+typ+' '+name)
 for l in range(10):
  T=(107*49**l-1)//2;need((107*49**(l+1)-1)//2==49*T+24,'threshold recurrence control')
  need(prefix(l+1)==[1,3,6,19]+[49*v+(-4 if j%2==0 else 4) for j,v in enumerate(prefix(l))],'closed prefix recurrence control')
 report={'status':'passed','checks':checks,'base_pairs':2568,'base_chunks':len(dep['chunks']),'paired_states':98,'direct_states':98,'paired_depth2_expansions':98,'direct_cycle_expansions':196,'input_base_indices':[53,54],'recursive_input_indices':sorted(len(x[0]) for x in lifted.values()),'old_transition_tables_used':0,'producer_imports':0,'hostile_mutations':mutations,'base_manifest_sha256':dep['manifest_sha256'],'certificate_sha256':file_hashes}
 print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__':main()
