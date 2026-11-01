from pathlib import Path
import hashlib,itertools,json,math,random
from stream import compile_linear,run,matrix,H,L,ATOMS,ordinary_table,tail_matrix

def require(x,m):
 if not x:raise RuntimeError(m)

def main():
 checks=0
 # Direct literal integer scans on an independently specified matrix universe.
 for d in (1,2):
  for flat in itertools.product(range(3),repeat=d*d):
   a=[list(flat[i*d:(i+1)*d]) for i in range(d)]
   for v in itertools.product(range(2),repeat=d):
    got,steps,p=run(a,v,3,literal=True);expected=list(v)
    for r in range(3):
     nxt=[sum(a[j][i]*expected[i] for i in range(d)) for j in range(d)]
     require(got[r+1]==nxt,'matrix semantics')
     require(steps[r]==1+2*d+2*sum(expected)+2*sum(nxt),'exact steps')
     expected=nxt;checks+=1
 a=matrix();v=[0]*len(a);v[0]=1
 got,steps,p=run(a,v,32)
 oracle=ordinary_table(32)
 for n,row in enumerate(got):
  for i,h in enumerate(H):
   for lag in range(L):
    require(row[i*L+lag]==(oracle[n-lag][h] if n>=lag else 0),'lag recurrence');checks+=1
  if n:require(steps[n-1]==1+2*len(a)+2*sum(got[n-1])+2*sum(row),'stream runtime');checks+=1
 lit,ls,_=run(a,v,10,literal=True)
 require(lit==got[:11] and ls==steps[:10],'actual stream integers');checks+=1
 tail=tail_matrix();tv=[oracle[6][h] for h in H]
 tg,ts,tp=run(tail,tv,26)
 for k,row in enumerate(tg):
  require(row==[oracle[k+6][h] for h in H],'tail counter');checks+=1
 tl,tls,_=run(tail,tv,5,literal=True)
 require(tl==tg[:6] and tls==ts[:5],'literal tail counter');checks+=1
 require((Path(__file__).parent/'tail21.fractran').read_text()=='\n'.join(f'{u}/{d}' for u,d in tp['fractions'])+'\n','tail fractions')
 # Closed forms from inherited SCC word enumeration, zero before support.
 for n in range(33):
  expected={1:1,3:max(0,n-2),5:max(0,n-3),7:0,9:(n-4)*(n-1)//2 if n>=5 else 0,
   11:max(0,n-4),13:max(0,n-4),15:(n-4)*(n-3) if n>=5 else 0,
   17:2*max(0,n-5),19:2*max(0,n-5),21:0}
  require(oracle[n]==expected,'closed forms');checks+=1
 # Initial vectors are parameters, not universal-valid tournament counts.
 hostiles=0
 for bad in ([],[[1,2]],[[1,-1],[0,1]],[[True]],[[0.5]]):
  try:compile_linear(bad)
  except ValueError:hostiles+=1
  else:raise RuntimeError('invalid matrix accepted')
 for badv in ([],[1,2],[-1],[True],[0.5]):
  try:run([[1]],badv,1)
  except ValueError:hostiles+=1
  else:raise RuntimeError('invalid vector accepted')
 # Losing next-step source buffer must disagree with simultaneous update.
 swap=[[0,1],[1,0]];old=[1,2];wrong=old[:]
 for j in range(2):wrong[j]=sum(swap[j][i]*wrong[i] for i in range(2))
 correct,_,_=run(swap,old,1,literal=True)
 require(wrong!=correct[-1] and correct[-1]==[2,1],'in-place hostile');hostiles+=1
 # A fake strong7 atom makes the 7 channel nonzero: missing7 is a proved atlas fact.
 altered=[r[:] for r in a];altered[H.index(7)*L][2]+=1
 bad,_,_=run(altered,v,3)
 require(bad[-1][H.index(7)*L]==1,'fake atom hostile');hostiles+=1
 generated='\n'.join(f'{u}/{d}' for u,d in p['fractions'])+'\n'
 base=Path(__file__).parent
 require((base/'bounded21.fractran').read_text()==generated,'fraction regeneration')
 require(json.loads((base/'program.json').read_text())=={k:v for k,v in p.items() if k!='fractions'},'metadata regeneration')
 return {'status':'passed','checks':checks,'negative_controls':hostiles,'fractions':len(p['fractions']),
  'dimension':len(a),'tail_fractions':len(tp['fractions']),'tail_literal_rounds':5,'tail_order_range':[6,32],'literal_integer_stream_through':10,'valuation_stream_through':32,
  'checkpoint':p['checkpoint'],'clock':p['clock'],'initial_integer':2*p['checkpoint'],
  'steps_through_32':sum(steps),'rows_0_through_32':[{str(h):oracle[n][h] for h in H} for n in range(33)],
  'fraction_sha256':hashlib.sha256(generated.encode()).hexdigest()}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
