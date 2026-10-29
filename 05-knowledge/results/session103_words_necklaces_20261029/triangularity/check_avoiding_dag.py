"""Exact audit of vertex DAG certificate; proof is in companion note."""
import sys,json,collections
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'singularity'))
from actual_grammar import build
from native_kernel import observer,need

def dag(A,C):
 n=len(A);s={C[i-1]:C[i] for i in range(n)}
 E={(s[u],v) for u in range(n) for v in range(n) if A[u][v] and s[u]!=v}
 ind=[0]*n
 for u,v in E:ind[v]+=1
 todo=[v for v in range(n) if not ind[v]];order=[]
 while todo:
  u=todo.pop();order.append(u)
  for a,v in E:
   if a==u:
    ind[v]-=1
    if not ind[v]:todo.append(v)
 return E,order,s

def main():
 here=Path(__file__).resolve().parent
 data=json.loads((here.parent/'singularity/pinned/tournament_results.json').read_text())
 rows=[]
 for ns,layer in data['layers'].items():
  n=int(ns);ct=collections.Counter()
  for rec in layer['records']:
   w=rec['word'];A,C,H=build(w);E,order,s=dag(A,C);avoid='be' not in w and 'ee' not in w
   need((len(order)==n)==avoid,'finite language/DAG mismatch')
   ct['words']+=1;ct['avoiding']+=avoid
   if not avoid:continue
   if w:
    B,K,_=build(w[:-1]);F,_,t=dag(B,K);p=H[-1]['parent'];q=t[p];x=n-1
    predicted={(a,v) for a,v in F if a!=q}|{(x,v) for a,v in F if a==q}|{(x,q)}|{(q,v) for v in range(n-1) if v not in (p,q)}
    need(E==predicted,'insertion rule mismatch')
    incoming={a for a,v in F if v==q}
    need(not incoming or (incoming=={p} and not any(v==p for a,v in F)),'preservation hypothesis failed')
    ct['insertion_checks']+=1
   if n<=8:
    ports,M=observer(A);idx={p:i for i,p in enumerate(ports)};pred={s[u]:u for u in s};f={v:i for i,v in enumerate(order)}
    for a,b in ports:
     for c,d in ports:
      v=M[idx[a,b]][idx[c,d]]
      if not v:continue
      if (c,d)==(pred[b],pred[a]):need(v==1,'bad canonical pivot')
      else:need(f[a]+f[b]<f[s[c]]+f[s[d]],'noncanonical edge not strictly increasing')
    ct['observer_potential_checks']+=1
  rows.append({'order':n,**ct});print(rows[-1],flush=True)
 out={'status':'FINITE-EXACT audits, with separate all-order proof candidate','source':'d3132eaf actual grammar','layers':rows,'controls':'All forbidden words through12 are cyclic; all avoiding words satisfy insertion hypotheses. Native observer computed independently through8, all entries checked against vertex-pair potential.'}
 (here/'avoiding_dag_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
