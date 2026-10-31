#!/usr/bin/env python3
"""Independent strong scale25 free-head audit; no producer imports."""
import argparse,copy,gzip,hashlib,json,math,struct
from pathlib import Path
checks=0
Q=25;HALF=12
def need(x,msg):
 global checks
 checks+=1
 if not x:raise ValueError(msg)
def square(x):return x>0 and math.isqrt(x)**2==x
def digest(x):return hashlib.sha256(x).hexdigest()
def selected(p,d,c):return int(Q+(-1 if p==0 else 1)*c<=d)
def endpoint(p):return 8 if p==0 else 10
def threshold(l):return (85*25**l-1)//2
def depth(n):
 l=0
 while threshold(l+1)<=n:l+=1
 return l
def prefix(l):
 out=[]
 for t in range(l+1):
  Q=25**t;block=[(2*Q+1)//3,(10*Q-1)//3,(17*Q+1)//3,(58*Q-1)//3]
  out+=block[:2] if t==l else block
 return out

def compare_profiles(profiles):
 need(len(profiles)==2 and set(profiles[0])==set(profiles[1]),'paired profile source types')
 for key in profiles[0]:
  need(profiles[0][key]==profiles[1][key],'shared head parity' if key[0]=='head' else 'shared singleton/tail parity')

def verify_table(w,direct):
 p,d0=w['parity'],w['residue'];need(w['radix']==25 and p in [0,1] and 12<=d0<=36,'Q25 state')
 need(bool(w.get('terminal_cycle',False))==direct,'paired versus terminal kind')
 need(w.get('protected_offset')==8 and w.get('protected_prefix')==[1,3,6,19],'strong free-head prefix family')
 if not direct:need(w.get('independent_head_colors') is True,'explicit free-head paired semantics')
 cuts=w['split_prefix_blocks'];need(len(cuts)==len(set(cuts)) and set(cuts)<=set(range(-12,13)) and 8 not in cuts,'Q25 cuts')
 need(len(w['words'])==(1 if direct else 2),'Q25 output count');profiles=[]
 for output_index,word in enumerate(w['words']):
  d=d0+output_index
  expected={('single',v) for v in range(1,13)}|{(kind,c) for c in range(-12,13) for kind in (('head','tail') if c in cuts else ('block',))}
  need(word[:5]==[['single',1],['single',3],['single',6],['head',6],['block',8,selected(p,d,8),1]],'Q25 exact protected opening')
  used=set();coef=const=0;last=None;profile={};head={};shift={}
  for atom in word:
   kind,c=atom[:2];key=(kind,c);need(key in expected and key not in used,'Q25 atom ownership');used.add(key);start=(coef*p+const+1)%2
   if kind in ['single','head']:
    need(len(atom)==2,'Q25 singleton format');first=final=c if kind=='single' else Q-c;da,db=0,1;profile[key]=start
    if kind=='head':head[c]=(start-1)%2
   else:
    need(len(atom)==4,'Q25 block format');_,c,i,orientation=atom;need(i==selected(p,d,c) and orientation in [-1,1],'Q25 forced member and orientation')
    L=(p+i)%2;first=Q-c if kind=='block' else 3*Q+c;final=Q*endpoint(L)+(c if L==0 else -c)
    if orientation==-1:first,final=final,first
    da,db=1,i-int(kind=='tail');delta=(start-(1 if kind=='block' else 2))%2 if orientation==1 else (start+L)%2
    profile['lifted_tail_or_block',c]=delta;shift[c]=delta
   if last is None:need(first==1,'Q25 first endpoint')
   else:need(square(last+first),'Q25 connector square')
   last=final;coef+=da;const+=db
  need(used==expected and (coef,const)==(Q,d),'Q25 symbolic interval length/ownership')
  need(last==(8 if direct else endpoint((p+d)%2)),'Q25 final endpoint')
  if direct:need(square(last+1),'Q25 closing square')
  else:
   # Heads and tails are separate common-label classes. No within-word linkage.
   profiles.append(profile)
 if not direct:compare_profiles(profiles)
 return profiles

def vertices(pair,word):
 for atom in word:
  kind,c=atom[:2]
  if kind=='single':yield c
  elif kind=='head':yield Q-c
  else:
   _,c,i,o=atom;P=pair[i];lo=int(kind=='tail');js=range(lo,len(P)) if o==1 else range(len(P)-1,lo-1,-1)
   for j in js:yield Q*P[j]+(c if j%2 else -c)
def validate(pair,pref):
 n=len(pair[0]);need(len(pair)==2 and len(pair[1])==n+1,'Q25 input consecutive lengths');colors=bytearray(n+2)
 for i,P in enumerate(pair):
  need(P[:len(pref)]==pref and P[-1]==endpoint(len(P)%2),'Q25 concrete B prefix/end');seen=bytearray(len(P)+1)
  for j,v in enumerate(P):
   need(1<=v<=len(P) and not seen[v],'Q25 concrete B permutation');seen[v]=1
   if i==0:colors[v]=j%2
   elif v<=n:need(colors[v]==j%2,'Q25 concrete B compatibility')
   if j:need(square(P[j-1]+v),'Q25 concrete B square')
  need(sum(seen)==len(P),'Q25 concrete B full coverage')

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 for key in ['pair-dir','cycle-dir','small-base-file','base-dir']:ap.add_argument('--'+key,type=Path,required=True)
 ap.add_argument('--partial',action='store_true',help='audit available states without claiming a complete theorem')
 args=ap.parse_args();hashes={};pairs={};cycles={}
 for directory,label,target,direct in [(args.pair_dir,'B',pairs,False),(args.cycle_dir,'C',cycles,True)]:
  actual={p.name for p in directory.glob(label+'_p*_d*.json')};expected={f'{label}_p{p}_d{d}.json' for p in [0,1] for d in range(12,37)}
  need(actual<=expected and (args.partial or actual==expected),'Q25 physical catalog')
  for p in [0,1]:
   for d in range(12,37):
    name=f'{label}_p{p}_d{d}.json'
    if args.partial and not (directory/name).exists():continue
    raw=(directory/name).read_bytes();hashes[name]=digest(raw);w=json.loads(raw);need((w['parity'],w['residue'])==(p,d),'Q25 filename dispatch');verify_table(w,direct);target[p,d]=w
 smallraw=args.small_base_file.read_bytes();small=json.loads(smallraw);need(set(small)=={str(n) for n in range(42,53)},'Q25 small catalog')
 for n in range(42,53):need(len(small[str(n)][0])==n,'Q25 small base index');validate(small[str(n)],[1,3])
 # The changing c6 block at p0,d18 gives a real counterexample to the old linkage rule.
 semantic=verify_table(pairs[0,18],False)
 need(semantic[0]['head',6]==0 and semantic[0]['lifted_tail_or_block',6]==0,'actual head-tail linkage counterexample')
 hostile=copy.deepcopy(semantic);hostile[1]['head',6]^=1
 need(sum(hostile[1][k]!=semantic[1][k] for k in semantic[1])==1,'only one head color mutated')
 try:compare_profiles(hostile)
 except ValueError as e:head_mutation=str(e)
 else:raise ValueError('single-head parity corruption accepted')
 if args.partial:
  for (p,d),w in pairs.items():
   source=small[str(42+p)];out=[list(vertices(source,word)) for word in w['words']];validate(out,prefix(1));need(len(out[0])==25*(42+p)+d,'partial concrete index')
  for (p,d),w in cycles.items():
   P=list(vertices(small[str(42+p)],w['words'][0]));N=25*(42+p)+d;need(len(P)==N and set(P)==set(range(1,N+1)) and P[:6]==prefix(1),'partial direct ownership/prefix')
   for j,v in enumerate(P):need(square(P[j-1]+v),'partial direct square including closure')
  print(json.dumps({'status':'available_certificates_passed_not_complete_theorem','checks':checks,'paired_states':len(pairs),'terminal_states':len(cycles),'within_word_counterexample':{'state':[0,18],'head_offset':6,'head_position_parity':0,'tail_shift':0,'old_required_head_parity':1},'single_head_color_mutation_rejected':head_mutation,'certificate_sha256':hashes},indent=2,sort_keys=True));return
 manifest=(args.base_dir/'B_base_chunks.json').read_bytes();need(digest(manifest)=='356ecf4ca3f6207a97db3f8b315278480709a2978f66821b6c3e030dbed3eae5','frozen old base manifest');chunks=json.loads(manifest);cache={};base_count=11
 def decode(c,n):
  raw=cache[c['file']];at=16+sum(2*k-1 for k in range(c['first'],n));pair=[]
  for L in [n,n+1]:
   roots=raw[at:at+L-1];need(len(roots)==L-1,'Q25 base root count');P=[1]
   for r in roots:P.append(r*r-P[-1])
   at+=L-1;pair.append(P)
  return pair
 for c in chunks:
  if c['first']>1061:continue
  raw=(args.base_dir/c['file']).read_bytes();need(len(raw)==c['bytes'] and digest(raw)==c['sha256'],'Q25 compressed base hash');data=gzip.decompress(raw);need(len(data)==c['decoded_bytes'] and digest(data)==c['decoded_sha256'],'Q25 decoded base hash');need(data[:8]==b'SSBAS1\0\0' and struct.unpack_from('<II',data,8)==(c['first'],c['last']),'Q25 base header');cache[c['file']]=data
  for n in range(c['first'],min(c['last'],1061)+1):validate(decode(c,n),[1,3]);base_count+=1
 need(base_count==1020,'Q25 full42 through1061 base supply')
 def build(n):
  need(n>=42,'Q25 B lower bound')
  if n<=1061:
   if n<53:pair=small[str(n)]
   else:pair=decode(next(c for c in chunks if c['first']<=n<=c['last']),n)
  else:
   m,d=divmod(n-12,25);d+=12;need(42<=m<n and 12<=d<=36,'Q25 well-founded decomposition');source=build(m);pair=[list(vertices(source,w)) for w in pairs[m%2,d]['words']]
  need(len(pair[0])==n,'Q25 B actual index');validate(pair,prefix(depth(n)));return pair
 for n in [42,43,1062,1063]:
  source=build(n)
  for d in range(12,37):
   pair=[list(vertices(source,w)) for w in pairs[n%2,d]['words']];validate(pair,prefix(depth(n)+1));need(len(pair[0])==25*n+d,'Q25 expanded paired index')
 controls=[]
 for N in [25*n+d for n in [42,43,1062,1063] for d in range(12,37)]+[26561,664062,16601562]:
  m,d=divmod(N-12,25);d+=12;pair=build(m);seen=bytearray(N+1);first=last=None;count=0;want=prefix(depth(N));actual=[];h=hashlib.sha256()
  for v in vertices(pair,cycles[m%2,d]['words'][0]):
   need(1<=v<=N and not seen[v],'Q25 cycle exact ownership');seen[v]=1;count+=1;h.update((str(v)+'\n').encode())
   if len(actual)<len(want):actual.append(v)
   if first is None:first=v
   if last is not None:need(square(last+v),'Q25 cycle square edge')
   last=v
  need(count==N and sum(seen)==N,'Q25 cycle full coverage');need(first==1 and last==8 and square(first+last),'Q25 cycle closure');need(actual==want,'Q25 full unscaled ray prefix')
  controls.append({'N':N,'depth':depth(N),'prefix':want,'cycle_decimal_lines_sha256':h.hexdigest()})
 mutations=[]
 for direct,bank in [(False,pairs),(True,cycles)]:
  for label,change in [('protected_cut',lambda w:w['split_prefix_blocks'].append(8)),('wrong_choice',lambda w:w['words'][0][4].__setitem__(2,1-w['words'][0][4][2])),('missing_atom',lambda w:w['words'][0].pop()),('reversed_prefix',lambda w:w['words'][0][4].__setitem__(3,-1))]:
   w=copy.deepcopy(bank[0,12]);change(w)
   try:verify_table(w,direct)
   except ValueError as e:mutations.append([direct,label,str(e)])
   else:raise ValueError('Q25 hostile mutation accepted')
 print(json.dumps({'status':'passed','checks':checks,'paired_states':50,'terminal_states':50,'base_pairs':base_count,'base_interval':[42,1061],'recursive_cutoff':1062,'old_chunks_used':len(cache),'paired_expansions':100,'cycle_controls':controls,'thresholds':[threshold(l) for l in range(5)],'small_catalog_sha256':digest(smallraw),'certificate_sha256':hashes,'hostile_mutations':mutations,'producer_imports':0,'single_head_color_mutation_rejected':head_mutation,'old_within_word_linkage_counterexample':[0,18,6]},indent=2,sort_keys=True))
if __name__=='__main__':main()
