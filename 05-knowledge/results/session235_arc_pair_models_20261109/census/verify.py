#!/usr/bin/env python3
"""Independent subset DP, explicit permutations, and inherited marginal checks."""
import array,csv,itertools,json,pathlib,random,collections
ROOT=pathlib.Path(__file__).resolve().parent

def rows(n,c,both=False):
 a=[[0]*n for _ in range(n)]
 for i,j in itertools.combinations(range(n),2):
  d=c%3;c//=3
  if d==0 or (d==2 and both):a[i][j]=1
  if d==1 or (d==2 and both):a[j][i]=1
 return a

def dp(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for s in range(1,1<<n):
  for j in range(n):
   if s>>j&1:
    t=s^(1<<j)
    if t:d[s][j]=sum(d[t][i] for i in range(n) if a[i][j])
 return sum(d[-1])

def paths(a):
 return [list(p) for p in itertools.permutations(range(len(a))) if all(a[i][j] for i,j in zip(p,p[1:]))]

def read_csv(name):return list(csv.DictReader((ROOT/name).open()))
report=[];all_witnesses={};rng=random.Random(235)
for n in range(1,7):
 N=3**(n*(n-1)//2); v=array.array('i');v.frombytes((ROOT/f'counts{n}.bin').read_bytes()); assert len(v)==2*N
 rr=read_csv(f'classes{n}.csv');targets=[int(r['code']) for r in rr]
 check=range(N) if n<=5 else sorted(set([0,N-1]+targets+[rng.randrange(N) for _ in range(3000)]))
 for c in check:
  assert dp(rows(n,c))==v[c],(n,c,'oriented')
  assert dp(rows(n,c,True))==v[N+c],(n,c,'semicomplete')
 joint=read_csv(f'joint{n}.csv');semi=collections.Counter()
 assert all((int(r['H_semicomplete'])-int(r['H_oriented']))%2==0 for r in joint)
 assert sum(int(r['count']) for r in joint if r['H_semicomplete']==r['H_oriented'])==2**(n*(n-1)//2)
 for r in joint:semi[int(r['H_semicomplete'])]+=int(r['count'])
 if n==6:
  inherited=ROOT.parents[1]/'session234_digraphs_20261109/census/results6.csv'
  expect={int(r['hamiltonian_paths']):int(r['labelled_count']) for r in csv.DictReader(inherited.open())}
  assert dict(semi)==expect
 witnesses=[]
 for r in rr:
  c=int(r['code']);a=rows(n,c);ps=paths(a);assert len(ps)==int(r['H'])
  missing=[(i,j) for i,j in itertools.combinations(range(n),2) if not a[i][j] and not a[j][i]]
  ch=collections.Counter();completion_examples={}
  for word in itertools.product((0,1),repeat=len(missing)):
   b=[row[:] for row in a]
   for (i,j),bit in zip(missing,word):b[i][j]=1-bit;b[j][i]=bit
   h=dp(b);ch[h]+=1;completion_examples.setdefault(h,{'matrix':[''.join(map(str,row)) for row in b],'added_arcs':[[j,i] if bit else [i,j] for (i,j),bit in zip(missing,word)]})
  assert min(ch)>=len(ps)
  witnesses.append({**{k:int(val) for k,val in r.items()},'matrix':[''.join(map(str,row)) for row in a],'paths':ps,'missing_pairs':missing,'H_semicomplete_same_raw':v[N+c],'completion_histogram':dict(sorted(ch.items())),'completion_examples':completion_examples})
 all_witnesses[str(n)]=witnesses
 report.append({'n':n,'labelled_universe':N,'independent_DP_checked':len(check),'all_target_classes_path_enumerated':len(rr),'inherited_semicomplete_histogram_matches':n==6})
 print(report[-1],flush=True)
(ROOT/'witnesses.json').write_text(json.dumps(all_witnesses,indent=2)+'\n')
(ROOT/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
