#!/usr/bin/env python3
"""Exact nested-sequence controls and Pascal-word DAG through level 1024.
The all-level identification of the Pascal word with Conway's differences is
classical; here that identification is independently checked through 2^20.
The compressed Pascal-word calculations themselves are exact at every level.
"""
import argparse,json,math
from array import array

def require(x,m):
 if not x:raise RuntimeError(m)
def join(x,y):
 lx,sx,mx,ax=x;ly,sy,my,ay=y
 if mx>=sx+my:return lx+ly,sx+sy,mx,ax
 return lx+ly,sx+sy,sx+my,lx+ay
ONE=(1,1,1,1);ZERO=(1,-1,0,0);EMPTY=(0,0,0,0)
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);args=p.parse_args();N=2**20
 a=array('I',[0])*(N+1);q=array('I',[0])*(N+1);depth=array('I',[0])*(N+1)
 a[1]=a[2]=q[1]=q[2]=1;controls=0;first_q_jump=None
 for n in range(3,N+1):
  j=a[n-1];require(1<=j<n and 1<=n-j<n,'source rank');a[n]=a[j]+a[n-j]
  require(2*a[n]-n==(2*a[j]-j)+(2*a[n-j]-(n-j)),'leaf charge');controls+=2
  if n>=4:
   prev=a[n-1]-a[n-2];parent=j if prev else n-j
   require(prev in (0,1) and 2<=parent<n,'copy index')
   require(a[n]-a[n-1]==a[parent]-a[parent-1],'copy bit');depth[n]=depth[parent]+1;controls+=2
  u,v=n-q[n-1],n-q[n-2];require(1<=u<n and 1<=v<n,'Q finite domain');q[n]=q[u]+q[v];controls+=1
  if first_q_jump is None and q[n]-q[n-1] not in (0,1):first_q_jump=[n,q[n-1],q[n]]
 # Byte words independently concatenate adjacent Pascal cells, never call a().
 row_words=[b'\x01',b'\x00'];row=[ONE,ZERO];reports=[]
 for k in range(1,1025):
  total=EMPTY
  for cell in row:total=join(total,cell)
  require(total[0]==2**k and total[1]==0,'Pascal length and balance');controls+=2
  if k<=19:
   word=b''.join(row_words);actual=bytes(a[n]-a[n-1] for n in range(2**k+1,2**(k+1)+1))
   require(word==actual,'independent Pascal word');controls+=len(word)
  if k in [1,2,3,4,5,10,20,32,64,128,256,512,1024]:
   L,S,M,where=total
   reports.append(dict(level=k,word_length=str(L),max_charge=str(M),first_max_index=str(2**k+where),max_charge_fraction=float(M/L),scaled_fraction=float(M/L*math.sqrt(k)),central_cell_length=str(math.comb(k,k//2))))
  row=[ONE]+[join(row[j-1],row[j])for j in range(1,len(row))]+[ZERO]
  if k<19:row_words=[b'\x01']+[row_words[j-1]+row_words[j]for j in range(1,len(row_words))]+[b'\x00']
 special=[dict(n=n,Conway=a[n],Q=q[n],charge=2*a[n]-n,copy_depth=depth[n])for n in [20,1024,N]]
 out=dict(status='FINITE-EXACT for recurrence through 2^20; exact Pascal DAG through level 1024',N=N,controls=controls,max_copy_depth=max(depth),first_Q_non_slow_increment=first_q_jump,special=special,pascal=reports)
 with open(args.output,'w')as f:json.dump(out,f,indent=2)
 print(json.dumps({k:v for k,v in out.items()if k!='pascal'},indent=2))
 print('Pascal scaled envelope at 20,1024',[(r['level'],r['scaled_fraction'])for r in reports if r['level'] in (20,1024)])
if __name__=='__main__':main()
