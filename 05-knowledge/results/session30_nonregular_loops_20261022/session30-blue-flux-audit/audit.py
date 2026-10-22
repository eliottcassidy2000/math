"""Independent pair-mask DP, matching anti-involutions, and full-word verifier."""
import json,pathlib,subprocess,collections,hashlib
HERE=pathlib.Path(__file__).resolve().parent
BASE=HERE.parent
root=json.loads((BASE/'session30-root/blue-flux-census.json').read_text())
controls=json.loads((BASE/'session30-root/completion-blue-gates.json').read_text())
bank=json.loads((BASE/'session30-almost-regular/census8_bank.json').read_text())['almostregular8']
def matrix(bits):
 n=next(n for n in range(2,11) if n*(n-1)//2==len(bits));a=[[0]*n for _ in range(n)];k=0
 for i in range(n):
  for j in range(i+1,n):a[i][j]=1 if bits[k]=='1' else -1;a[j][i]=-a[i][j];k+=1
 return a
def matchings(vertices):
 if not vertices:yield [];return
 i=vertices[0]
 for j in vertices[1:]:
  for rest in matchings([v for v in vertices if v not in (i,j)]):yield [(i,j)]+rest
def anti_involutions(a):
 n=len(a)
 for pairs in matchings(list(range(n))):
  s=[0]*n
  for i,j in pairs:s[i]=j;s[j]=i
  if all(a[s[i]][s[j]]==-a[i][j] for i in range(n) for j in range(n)):yield s
def dp(a,s):
 n=len(a);pairs=[(i,s[i]) for i in range(n) if i<s[i]];pair={v:k for k,p in enumerate(pairs) for v in p};states=[[0]*n for _ in range(1<<len(pairs))]
 for v in range(n):states[1<<pair[v]][v]=1
 for mask in range(1,len(states)):
  for v,count in enumerate(states[mask]):
   if count:
    for w in range(n):
     bit=1<<pair[w]
     if not mask&bit and a[v][w]==1:states[mask|bit][w]+=count
 return states[-1]
def encode(a):return ''.join('1' if a[i][j]==1 else '0' for i in range(len(a)) for j in range(i+1,len(a)))
def complete(a):
 d=list(map(sum,a));return [a[i]+[-d[i]] for i in range(len(a))]+[d+[0]]
def full(a,s):
 inp=f'{len(a)} {encode(a)} '+ ' '.join(map(str,s))+'\n'
 values=list(map(int,subprocess.check_output([str(HERE/'full_words')],input=inp,text=True).split()))
 return dict(zip(['aut','anti','centralizer','H','labelled_blue','fixed_words'],values))
subprocess.run(['g++','-O3','-std=c++17',str(HERE/'full_words.cpp'),'-o',str(HERE/'full_words')],check=True)
found={r['bits']:list(anti_involutions(matrix(r['bits']))) for r in bank};found={b:s for b,s in found.items() if s}
assert set(found)=={r['bits'] for r in root['records']}
records=[]
for r in root['records']+controls['records']:
 bits=r['bits'];a=matrix(bits);n=len(a);s=r.get('sigma',r.get('anti_involution'));d=list(map(sum,a));h=dp(a,s)
 assert h==r.get('endpoint_halfpaths',r.get('halfpath_endpoint_counts'))
 even=sum(h[v] for v in range(n) if a[v][s[v]]==1);odd=sum(h[v] for v in range(n) if d[v]==-1)
 T=full(a,s);R=full(complete(a),s+[n]);assert T['centralizer']==R['centralizer']==r['centralizer']
 c=T['centralizer'];assert T['fixed_words']==even and R['fixed_words']==odd
 assert T['labelled_blue']*c==even*T['aut'] and R['labelled_blue']*c==odd*R['aut']
 assert even%c==odd%c==0
 all_s=list(anti_involutions(a));assert len(all_s)*c==T['aut']
 for ss in all_s:
  hh=dp(a,ss)
  assert sum(hh[v] for v in range(n) if a[v][ss[v]]==1)==even
  assert sum(hh[v] for v in range(n) if d[v]==-1)==odd
 flux=(even-odd)//c;assert flux==r.get('flux',r.get('gate_flux'))
 records.append(dict(bits=bits,n=n,sigma=s,endpoints=h,T=T,R=R,blue_T=even//c,blue_R=odd//c,flux=flux,anti_involutions=len(all_s)))
 print(n,bits,flux,flush=True)
receipt=dict(status='PASS',method='Perfect-matching anti-involutions; pair-mask DP; full permutation Hamilton-word/automorphism checker',selfconverse_order8=len(found),distribution=dict(sorted(collections.Counter(x['flux'] for x in records if x['n']==8).items())),records=records)
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
# Literal witnesses for all four terminal gate patterns on a negative-flux
# class and a positive-flux class with nontrivial centralizer c=3.
import itertools
witnesses=[]
for record in [records[1],records[2]]:
 a=matrix(record['bits']);n=len(a);s=record['sigma'];rr=complete(a);seen={}
 for p in itertools.permutations(range(n),n//2):
  if len({min(v,s[v]) for v in p})!=n//2:continue
  if not all(a[p[i]][p[i+1]]==1 for i in range(len(p)-1)):continue
  tail=tuple(s[v] for v in reversed(p));t=p+tail;r=p+(n,)+tail
  tg=all(a[t[i]][t[i+1]]==1 for i in range(n-1));rg=all(rr[r[i]][r[i+1]]==1 for i in range(n))
  assert len(set(t))==n and len(set(r))==n+1
  key=f'{int(tg)}{int(rg)}'
  if key not in seen:seen[key]=dict(halfpath=p,T_word=t,R_word=r,T_is_directed=tg,R_is_directed=rg)
 witnesses.append(dict(bits=record['bits'],sigma=s,centralizer=record['T']['centralizer'],flux=record['flux'],gate_patterns=seen))
(HERE/'explicit-paths.json').write_text(json.dumps(witnesses,indent=2)+'\n')
