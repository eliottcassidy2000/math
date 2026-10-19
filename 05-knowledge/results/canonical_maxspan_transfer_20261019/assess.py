from pathlib import Path
from itertools import combinations
from collections import Counter
import json,gzip,hashlib,time
ROOT=Path(__file__).resolve().parent

def require(x,msg):
 if not x:raise RuntimeError(msg)

def assess(p,c):
 n=len(p);D=[{i for i in range(n) if p[j]>>i&1} for j in range(n)];U=[{j for j in range(n) if i in D[j]} for i in range(n)];I=[set(range(n))-{i}-D[i]-U[i] for i in range(n)];inc=[(i,j) for i,j in combinations(range(n),2) if j in I[i]];deg=list(map(len,I));seeds={0,1}|{i for i,d in enumerate(deg) if d==7}
 require(all(v not in D[v] and all(D[a]<=D[v] for a in D[v]) for v in range(n)),'strict transitive order');require(max(deg)<=7,'degree');require(D[0]<=D[1] and U[1]<=U[0] and 1 in I[0],'critical roots');require(len(I[0]&I[1])==c,'common count');root_sum=deg[0]+deg[1];require(root_sum==n+c,'root degree sum')
 def H(a,b):return len(D[b]-D[a])+len(U[a]-U[b])+(0 if (D[a]|{a})&seeds else 7-deg[a])+(0 if (U[b]|{b})&seeds else 7-deg[b])
 FA={v for v in range(n) if (D[v]|{v})&seeds};FB={v for v in range(n) if (U[v]|{v})&seeds};anchored=FA&FB
 def excluded_morph(l,u,c):
  a,b=sorted((l,u))
  return a<=1 or (a,b)==(2,2) and c<=3 or (a,b)==(2,3) and c<=2 or (a,b) in [(2,4),(3,3)] and c==0
 bounds={(a,b):H(a,b) for u,v in inc for a,b in [(u,v),(v,u)]};common={(a,b):len(I[a]&I[b]) for a,b in inc};hist=Counter();witnesses=[];total=0;baseline_digest=hashlib.sha256()
 def visit(word,mask):
  nonlocal total
  if len(word)<n:
   for v in range(n):
    if not mask>>v&1 and not p[v]&~mask:visit(word+[v],mask|1<<v)
   return
  r={v:i for i,v in enumerate(word)};span=r[1]-r[0]
  if span<=0 or any(abs(r[a]-r[b])>span for a,b in inc):return
  if span<=3:return
  if any(bounds[(a,b) if r[a]<r[b] else (b,a)]<=1 for a,b in inc):return
  total+=1;baseline_digest.update(json.dumps(word).encode());tied=[(a,b) for a,b in inc if abs(r[a]-r[b])==span]
  for a,b in tied:require(common[a,b]!=c or deg[a]+deg[b]<=root_sum,'degree redundancy')
  bad=[(a,b) for a,b in tied if common[a,b]>c]
  secondary=[]
  for a,b in tied:
   x,y=(a,b) if r[a]<r[b] else (b,a)
   require(D[x]<=D[y] and U[y]<=U[x],'maximum-span criticality')
   cc=common[a,b]
   ll=len(D[y]-D[x]);uu=len(U[x]-U[y])
   alpha=0 if x in FA else 7-deg[x];beta=0 if y in FB else 7-deg[y]
   lower=[max(2,ll),min(6-cc,ll+alpha)];upper=[max(2,uu),min(6-cc,uu+beta)]
   feasible=[(a,b,cc) for a in range(lower[0],lower[1]+1) for b in range(upper[0],upper[1]+1)]
   if all(excluded_morph(a,b,t) for a,b,t in feasible):
    secondary.append({'pair':[x,y],'reason':'empty_rectangle' if not feasible else 'all_morphologies_excluded','common_count':cc,'core_lower_upper':[ll,uu],'alpha_beta':[alpha,beta],'lower_interval':lower,'upper_interval':upper,'feasible_morphologies':feasible})
  hist['secondary_any_root_excluded']+=bool(secondary)
  hist['union_excluded']+=bool(bad or secondary)
  hist['secondary_only']+=bool(secondary and not bad)
  if secondary:hist['secondary_'+secondary[0]['reason']]+=1
  if bad:
   a,b=bad[0];hist['canonical_excluded']+=1;hist[f'witness_common_{common[a,b]}']+=1;witnesses.append({'order':word,'pair':[a,b],'common':sorted(I[a]&I[b]),'span':span,'secondary_witness':secondary[0] if secondary else None})
  elif secondary:witnesses.append({'order':word,'span':span,'secondary_witness':secondary[0]})
  else:hist['combined_retained']+=1
 visit([],0)
 return {'baseline':total,'hist':dict(hist),'baseline_digest':baseline_digest.hexdigest(),'witnesses':witnesses}

def main():
 fixtures=json.loads((ROOT/'fixtures.json').read_text());summaries=[];details=[];start=time.time()
 for fixture in fixtures:
  l,u,c=fixture['scope'];h=Counter();fully=0
  for core in fixture['cores']:
   p=core['pred'];require(len(p)==2+l+u+c,'core size');require((p[1]&~p[0]).bit_count()==l,'lower count');require(sum(bool(p[v]&1) and not bool(p[v]&2) for v in range(len(p)))==u,'upper count');result=assess(p,c);require(result['baseline']==core['expected_residual_orders'],'baseline mismatch');h['baseline']+=result['baseline'];h.update(result['hist']);fully+=int(result['baseline']>0 and result['hist'].get('combined_retained',0)==0);details.append({'scope':[l,u,c],'core_index':core['core_index'],**result})
  s={'scope':[l,u,c],'sample_cores':len(fixture['cores']),'counts':dict(h),'nonempty_cores_fully_removed':fully};summaries.append(s);print(json.dumps(s),flush=True)
 (ROOT/'impact.json').write_text(json.dumps({'status':'BOUNDED_CANONICAL_IMPACT_VERIFIED','summary':summaries,'scope':'Four fixed 100-core samples; not population coverage'},indent=2));(ROOT/'impact_details.json.gz').write_bytes(gzip.compress(json.dumps(details,separators=(',',':')).encode(),mtime=0))
if __name__=='__main__':main()
