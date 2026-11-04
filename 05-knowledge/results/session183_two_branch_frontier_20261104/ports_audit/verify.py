import json
from pathlib import Path
P=Path(__file__).resolve().parent
source=P.parent/'ports'

def need(x,msg):
 if not x:raise RuntimeError(msg)

def templates(rows,omit=None,slots=2):
 n=len(rows)-1; root=n; old=[i for i in range(n) if i!=omit]; cut=[bool(rows[i]>>root&1) for i in range(n)];out=[]
 for a in [-1]+[i for i in old if cut[i]]:
  for b in [-1]+[i for i in old if not cut[i]]:
   allow=[]
   for i in old:
    z=(4 if slots==2 else 0)|(1 if a>=0 and i<=a else 0)|(2 if b>=0 and i>=b else 0)
    if i==a:z&=1
    if i==b:z&=2
    allow.append(z)
   if all(allow):out.append(allow)
 return old,out

def covers(old,tt,slots=2,end=None,start=None):
 total=0
 for rs in range(slots):
  masks=[7]*len(old);bad=False
  for isend,item in [(True,end),(False,start)]:
   if item is None:continue
   slot,w=item
   if w==-1:
    if slot!=rs:bad=True;break
    ban=2 if isend else 1
    masks=[z&~ban for z in masks]
   else:
    bit=(2 if isend else 1) if slot==rs else 4
    if w not in old:bad=True;break
    for j,i in enumerate(old):
     if i==w:masks[j]&=bit
     elif (i>w if isend else i<w):masks[j]&=~bit
  if bad:continue
  for base in tt:
   prod=1
   for b,m in zip(base,masks):
    prod*=(b&m).bit_count()
    if not prod:break
   total+=prod
 return total

def fmul(a,b):
 # Polynomial long division with coefficients ascending, modulus x^3+x^2+1.
 z=[0]*5
 for i in range(3):
  for k in range(3):z[i+k]+=a[i]*b[k]
 for i in (4,3):
  q=z[i];z[i]=0;z[i-1]-=q;z[i-3]-=q
 return [v%7 for v in z[:3]]

def fpow(a,n):
 z=[1,0,0]
 for _ in range(n):z=fmul(z,a)
 return z

def main():
 data=json.loads((source/'receipt.json').read_text());checks=0;summary=[]
 for state in data['states']:
  rows=state['rows'];names=state['names'];need(names==list(range(23))+[-1],'names')
  n=len(rows);root=n-1
  for i in range(n):
   need(not (rows[i]>>i&1),'loop')
   for k in range(i+1,n):
    need(((rows[i]>>k)&1)+((rows[k]>>i)&1)==1,'tournament')
    if k<root:need(rows[i]>>k&1,'old transitivity')
  cut=sum(((rows[i]>>root)&1)<<i for i in range(root))
  need(cut==state['f']<<state['t'],'named cut')
  need(not(cut&1) and cut>>(root-1)&1,'strong')
  old,t2=templates(rows);_,t1=templates(rows,slots=1)
  h=covers(old,t1,1);need(h==state['H']==4611841,'H')
  for r in state['ports']:
   w=r['name'];vals={}
   for slots in (1,2):
    tt=t1 if slots==1 else t2
    vals['E'+str(slots)]=covers(old,tt,slots,end=(0,w))
    vals['S'+str(slots)]=covers(old,tt,slots,start=(0,w))
   if w==-1:deletion=1
   else:
    oldd,td=templates(rows,omit=w,slots=1);deletion=covers(oldd,td,1)
   vals['deletion_H']=deletion
   need(all(r[k]==v for k,v in vals.items()),('endpoint',state['t'],state['f'],w,vals))
   need(r['defect']==r['E2']-h,('defect',w))
   golden=[h,(r['S2']-r['E1'])//2,(r['S2']+r['E1'])//2]
   need(golden==r['golden'] and (r['S2']-r['E1'])%2==0,'golden')
   z=[x%7 for x in golden]; orbit=[z,fpow(z,7),fpow(z,49)]
   c=r['charts'];need(c['observer_mod7']==z and c['F343_vector_codec']==orbit and fpow(z,343)==z,'F343')
   degree=1 if orbit[1]==z else 3;need(c['F343_codec_degree']==degree,'degree')
   H,x,y=golden;a,b,q=(x+4*(H-y))%7,y%7,(H-y)%7
   need(c['F49_charge']==[a,b,q] and c['F49_frobenius']==[(a+4*b)%7,-b%7,q],'F49')
   checks+=1
  for w,v,j in state['joints']:
   need(covers(old,t2,end=(0,w),start=(1,v))==j,('joint',w,v));checks+=1
  summary.append(dict(t=state['t'],f=state['f'],H=h,order=len(names)))
 # Entire finite family census, using the low-bit recurrence rather than author high-bit scanner.
 N=4611841
 from array import array
 J=array('I',[0])*N
 for f in range(1,N):
  n=f//2;J[f]=J[n]+((n+(n&1)) if f&1 else 0)
 hits=[]
 for t in range(1,24):
  a=1<<t
  for f in range(1,(N-1)//(a-1)+1):
   if 1+a*f-J[f]==N:hits.append([t,f,J[f],1+t+f.bit_length()])
 expected=[[s['t'],s['f'],s['j'],len(s['names'])] for s in data['states']]
 need(hits==expected,'complete census')
 distinct=len({tuple(s['rows']) for s in data['states']})
 result=dict(status='PASS',method='Direct root-neighbor boundary enumeration and product of allowed slot choices, independent of author six-state transfer.',endpoint_and_joint_records=checks,all_six=summary,complete_census=hits,distinct_named_adjacencies=distinct,min_family_order=min(h[3] for h in hits),scope='All endpoint/joint/codec data and complete TT-family census independently checked; no global minimum-order or named-lineage claim.')
 (P/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
