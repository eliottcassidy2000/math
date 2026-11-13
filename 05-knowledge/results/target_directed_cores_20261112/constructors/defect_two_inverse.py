"""Positive target search among explicit multi-entry funnel cores; no completeness claim."""
from pathlib import Path
from itertools import combinations, product
import json, subprocess, sys
HERE=Path(__file__).resolve().parent
TARGET=53084231

def inherited(path):
 ns={'__name__':'inherited'}
 exec(compile(subprocess.check_output(['git','show','HEAD:'+path],text=True),path,'exec'),ns)
 return ns

def shell(n):
 out=1
 while n>1:
  m=n//2
  if n%2==0:out+=m+(m&1)
  n=m
 return out

def inverse(N,a,b,c,d,place=1,prefix=0):
 if N<d+a:return None
 lo=max(0,(N-d-a-c+a+b-1)//(a+b));hi=(N-d-a)//b
 if lo>hi:return None
 if hi-lo<=3:
  for n in range(lo,hi+1):
   if a*shell(n)+b*n+c*(n&1)+d==N:return prefix+place*n
  return None
 out=inverse(N,a,2*b+a,a,d,2*place,prefix)
 return out if out is not None else inverse(N,a,2*b,0,d+b+c,2*place,prefix+place)

def main():
 target=int(sys.argv[1]) if len(sys.argv)>1 else TARGET
 output=HERE/(sys.argv[2] if len(sys.argv)>2 else 'defect_two_inverse_receipt.json')
 ex=inherited('05-knowledge/results/polarized_path_calculus_20261103/explore.py')
 fun=inherited('05-knowledge/results/session296_multicore_metagraph_20261113/cores/verify_ready.py')
 bank={1:[((0,),())]};seen=set();tests=0
 for n in range(3,11):
  bank[n]=[]
  for m in (1,*range(3,n)):
   t=n-m-1;words=[''] if t==0 else [''.join(p)+'1' for p in product('01',repeat=t-1)]
   for T,addr in bank[m]:
    for word in words:bank[n].append((fun['funnel'](T,word),addr+(word,)))
  for T,address in bank[n]:
   if T in seen:continue
   seen.add(T)
   H,ends=ex['tables'](T);_,starts=ex['tables'](ex['opposite'](T));full=len(H)-1
   one=[r for r,x in enumerate(starts[-1]) if x==1]
   entries=list(combinations(one,2))+[(r,) for r,x in enumerate(starts[-1]) if x==2]
   p2=sum(H[m]*H[full^m] for m in range(full+1))
   for J in entries:
    E=[sum(row[r] for r in J) for row in ends];S=[sum(row[r] for r in J) for row in starts]
    F=sum(E[m]*H[full^m] for m in range(full+1));G=sum(S[m]*H[full^m] for m in range(full+1))
    C=sum(E[m]*S[full^m] for m in range(full+1));b=G+p2;d=b-E[-1]-C-F-2
    tests+=1;x=inverse(target,H[-1],b,E[-1],d)
    if x is None:continue
    word=tuple(map(int,bin(x)[2:][::-1])) if x else ()
    parent,r=ex['funnel'](T,J,word);child=fun['birth'](parent,r)
    # Huge literal graph is retained, not subjected to a general exponential DP.
    data=dict(status='CONSTRUCTED',target=target,core=T,core_address=address,entries=J,word=word,word_integer=x,
      parameters=dict(h=H[-1],p2=p2,e=E[-1],F=F,G=G,budget=2,C=C),parent=parent,root_parent=r,rows=child,root=len(child)-1,
      tested_entry_cores=tests,identity_value=H[-1]*shell(x)+b*x+E[-1]*(x&1)+d)
    assert data['identity_value']==target
    assert sum(starts[-1][j] for j in J)==2
    assert [i for i in range(len(child)-1) if child[i]>>(len(child)-1)&1]==[r]
    q=(x+1)*G+shell(x)*H[-1]-(1-(x&1))*E[-1]-C
    data['root_state_H_s_q_F']=[target,q-2,q,target+2]
    data['binary_shell_g']=shell(x)
    data['NR_balanced_H_s_q_F']=[3*target+4,target,target,3*target+4]
    from paired_scale_probe import evaluate
    checks=[];cases=[]
    for y in range(16):
     wy=tuple(map(int,bin(y)[2:][::-1])) if y else ()
     py,ry=ex['funnel'](T,J,wy);cy=fun['birth'](py,ry)
     cases.append((cy,len(cy)-1));checks.append({'x':y,'word':wy})
    for row,z in zip(checks,evaluate(cases)):
     y=row['x'];hy=H[-1]*shell(y)+b*y+E[-1]*(y&1)+d
     qy=(y+1)*G+shell(y)*H[-1]-(1-(y&1))*E[-1]-C
     assert z==(hy,qy,qy-2),(y,z,hy,qy)
     row['profile_H_q_s']=z
    data['short_word_controls']=checks
    if x%2==0:
     data['zero_prefix_ray']={'index':'ell>=0','H_coefficient':(H[-1]+b)*x,
      'H_constant':H[-1]*(shell(x)-x)+d,'q_coefficient':(G+H[-1])*x,
      'q_constant':H[-1]*(shell(x)-x)+G-E[-1]-C,'defect':2}
    output.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:data[k] for k in ('status','target','core','entries','word_integer','parameters','tested_entry_cores')}));return
   ex['tables'].cache_clear()
  print('finished core order',n,'tests',tests,flush=True)
 print('No hit in this explicitly scoped positive menu',flush=True)

if __name__=='__main__':main()
