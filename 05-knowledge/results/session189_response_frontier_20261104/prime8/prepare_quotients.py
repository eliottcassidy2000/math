from pathlib import Path
import json
import core
ROOT=Path(__file__).resolve().parent

def prime(T):
 n=len(T)
 for m in range(1,(1<<n)-1):
  if m.bit_count()<2:continue
  if all((T[v]&m) in (0,m) for v in range(n) if not m>>v&1):return False
 return True

def main():
 allrows=json.loads((ROOT/'cover.json').read_text())['levels'][-1]['representatives'];out=[]
 for idx,T in enumerate(allrows):
  if not prime(T):continue
  U=core.converse(T);E=core.endpoint(T);S=core.endpoint(U);p=[(s+e)%2 for s,e in zip(S,E)];row=dict(cover_index=idx,rows=T,S=S,E=E,p=p)
  if not any(p):
   sig=[(s-e)%4 for s,e in zip(S,E)]
   for j in range(8):
    a=core.repeated_endpoint(T,j);b=core.repeated_endpoint(U,j);sig=[(v+2*(y-x))%4 for v,x,y in zip(sig,a,b)]
   row['sigma']=sig
  out.append(row)
 (ROOT/'quotients.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
 print(json.dumps(dict(prime_representatives=len(out),endpoint_zero=sum(not any(r['p']) for r in out))))
if __name__=='__main__':main()
