import json,pathlib,hashlib,functools,copy
P=pathlib.Path(__file__).resolve().parent
checks=0
def need(x):
 global checks
 checks+=1
 if not x:raise RuntimeError('independent check failed')
def literal(r):
 A=[[int(x) for x in row] for row in r['source_rows']];n=len(A);m=n//2;L=r['high_cross_shore'];R=r['other_shore'];v=r['terminal']
 need(n in (10,14) and all(len(row)==n for row in A))
 need(set(L).isdisjoint(R) and set(L)|set(R)==set(range(n)) and len(L)==len(R)==m and v in R)
 need(all(A[i][i]==0 for i in range(n)) and all(A[i][j]+A[j][i]==1 for i in range(n) for j in range(i)))
 need(len({sum(A[x]) for x in L})==1 and len({sum(A[x]) for x in R})==1 and sum(A[L[0]])!=sum(A[R[0]]))
 patterns=[]
 for w in r['words']:
  need(type(w)is list and len(w)==n and all(type(x)is int for x in w) and set(w)==set(range(n)))
  need(all(A[x][y] for x,y in zip(w,w[1:])))
  patterns.append(''.join('A' if x in L else 'B' for x in w))
 need(patterns==['A'*i+'B'+'A'*(m-i)+'B'*(m-1) for i in range(1,m-1)]+['B'+'A'*m+'B'*(m-1)])
 mates=[x[::-1].translate(str.maketrans('AB','BA')) for x in patterns]
 need(r['patterns']==patterns and r['mate_patterns']==mates and len(set(patterns+mates))==2*m-2)
 return True
q=json.loads((P.parent/'session22-q2structural/receipt.json').read_text())
for r in q['sample_receipts']:literal(r)
# Reject source-word corruption in an independently implemented consumer.
bad=copy.deepcopy(q['sample_receipts'][0]);bad['words'][0][-1]=-1
try:literal(bad);raise ValueError('corruption accepted')
except RuntimeError:pass
# Independent backward Hamilton-path DP, and pair-mask blue-word DP.
f=json.loads((P.parent/'session22-order10/receipt.json').read_text());hist={};sc=0
for r in f['frames']:
 A=[[int(x) for x in row] for row in r['source_rows']];n=10
 @functools.lru_cache(None)
 def endings(mask,v):
  rest=mask^(1<<v)
  return 1 if not rest else sum(endings(rest,u) for u in range(n) if rest>>u&1 and A[u][v])
 H=sum(endings((1<<n)-1,v) for v in range(n));need(H==r['H'])
 B=0
 for a in r['anti_automorphisms']:
  need(sorted(a)==list(range(n)) and all(A[i][j]==A[a[j]][a[i]] for i in range(n) for j in range(n)))
  if any(a[a[x]]!=x or a[x]==x for x in range(n)):continue
  pairs=[(x,a[x]) for x in range(n) if x<a[x]]
  @functools.lru_cache(None)
  def blue(used,last):
   if used==(1<<5)-1:return A[last][a[last]]
   return sum(blue(used|1<<j,x) for j,pair in enumerate(pairs) if not used>>j&1 for x in pair if A[last][x])
  B+=sum(blue(1<<j,x) for j,pair in enumerate(pairs) for x in pair)
 need(B==r['b']*r['aut_count']);need(H==r['h']*r['aut_count'])
 for a in r['automorphisms']:need(sorted(a)==list(range(n)) and all(A[i][j]==A[a[i]][a[j]] for i in range(n) for j in range(n)))
 if r['anti_count']:
  sc+=1;d=(H-B)//r['aut_count'];need(d==r['delta'] and d!=6);hist[str(d)]=hist.get(str(d),0)+1
need(len(f['frames'])==128 and sc==64 and hist==f['summary']['self_converse_deficit_histogram'])
out={'checks':checks,'structural_literal_samples':len(q['sample_receipts']),'order10_frames_independent_H_blue_DP':128,'self_converse_frames':sc,'deficit_histogram':hist,'scope':'marked-frame controls, not ordinary-class census'}
(P/'independent_results.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
