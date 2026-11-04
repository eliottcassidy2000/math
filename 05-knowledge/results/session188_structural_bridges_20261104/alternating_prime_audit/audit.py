"""Independent W-family structure and literal-query audit; no producer imports."""
from pathlib import Path
import sys,json,hashlib,subprocess
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'alternating_prime_quotient'
sys.path.insert(0,str(BASE/'ordered_module_skew_audit'))
from audit import need,mod,induce,endpoint,op,hp,attach

def main():
 paths=[SRC/'PROOF.md',SRC/'verify.py',BASE/'marked_pair_quotient'/'PROOF.md',BASE/'marked_pair_quotient'/'verify.py',BASE/'ordered_module_skew'/'PROOF.md',BASE/'ordered_module_skew_audit'/'audit.py'];pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
 data=json.loads((SRC/'receipt.json').read_text());structures=0
 for r in range(1,9):
  n=2*r+1;h=n-1
  R=tuple(sum(1<<v for v in range(n) if (u<h and v<h and u<v) or (u==h and v<h and v%2==0) or(v==h and u<h and u%2==1)) for u in range(n))
  need(all(R[u]>>v&1 for u,v in zip(list(range(n)),list(range(1,n))+[0])),'explicit spanning directed cycle')
  need(all(not mod(R,M) for M in range(1,(1<<n)-1) if M.bit_count()>1),'independent primeness')
  for z in range(n):
   keep=[u for u in range(n) if u!=z]
   if z==h:A,B=[0],keep[1:]
   elif z==0:A,B=[1],[u for u in keep if u!=1]
   elif z==h-1:A,B=[u for u in keep if u!=h-2],[h-2]
   else:A,B=[z-1],[z+1]
   need(all(R[u]>>v&1 for u in A for v in B) and all(R[u]&sum(1<<v for v in A+B) in (0,sum(1<<v for v in A+B)) for u in keep if u not in A+B),'independent deletion structure');structures+=1
 batteries=queries=0;coefs=[]
 for rec in data['positive_literal_batteries']+data['zero_quartic_boundary_batteries']:
  T,U=tuple(rec['rows_T']),tuple(rec['rows_U']);z=rec['pair_index'];sizes=list(map(len,rec['children']));off=[sum(sizes[:i]) for i in range(len(sizes)+1)];removed=set(range(off[z],off[z+1]));keep=[u for u in range(len(T)) if u not in removed];R=induce(T,keep);index={u:i for i,u in enumerate(keep)}
  A,B=[{index[u] for child in group for u in range(off[child],off[child+1])} for group in rec['marked_child_sets']];D,E=endpoint(R),endpoint(op(R));N=(1<<len(R))-1;skew=sum(sum(D[M][u] for u in B)*sum(E[N^M][v] for v in A)-sum(D[M][u] for u in A)*sum(E[N^M][v] for v in B) for M in range(1,N));need(skew==rec['cover_skew']==rec['quartic_coefficient'],'independent cover skew');total=0
  for q in rec['literal_battery']:
   a,b=hp(attach(T,q['a'],q['c'])),hp(attach(U,q['a'],q['c']));need((a,b)==(q['T'],q['U']),'independent integer query values');total+=q['sign']*(a-b);queries+=1
  need(total==skew,'complete four-bit coefficient');coefs.append(skew);batteries+=1
 normal=subprocess.check_output([sys.executable,str(SRC/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(SRC/'verify.py')]);need(normal==optimized==(SRC/'receipt.json').read_bytes(),'producer replay equality')
 need(pins=={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'source freeze')
 print(json.dumps({'status':'PASS','independent_structural_deletions_through_order17':structures,'literal_batteries':batteries,'actual_queries_per_graph':queries,'quartic_values':coefs,'producer_modes_equal':True,'source_sha256':pins,'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))
if __name__=='__main__':main()
