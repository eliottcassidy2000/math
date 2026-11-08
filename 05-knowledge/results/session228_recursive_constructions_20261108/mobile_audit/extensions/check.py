import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from verify import *
from itertools import combinations

def multi(r,J):
 w=[0,1]*r;groups=[];i=0
 while i<len(w):
  if i//2 in J and i%2==0:groups.append((3,1));i+=2
  else:groups.append((1,w[i]));i+=1
 n=1+sum(s for s,k in groups);a=[[False]*n for _ in range(n)];off=1
 for size,k in groups:
  for z in range(size):
   u=off+z;a[0][u]=z<k;a[u][0]=z>=k
   for q in range(off+size,n):a[u][q]=True
  if size==3:
   for z in range(3):a[off+z][off+(z+1)%3]=True
  off+=size
 return a

def trim(a):
 ids=[0]+list(range(2,len(a)-1));return [[a[i][j] for j in ids] for i in ids]
def strong(a):
 for s in range(len(a)):
  reached={s};todo=[s]
  while todo:
   u=todo.pop()
   for v in range(len(a)):
    if a[u][v] and v not in reached:reached.add(v);todo.append(v)
  if len(reached)!=len(a):return False
 return True

def Dmulti(r,J):
 x=v
 A=mm(S[1],S[0])
 for j in range(r):x=mv(T[1] if j in J else A,x)
 return sum(a*b for a,b in zip(ell,x))
def predicted(r,J):
 t=len(J);B=(4**r+6*r-1)//3
 if not t:return B//3
 return 3**(t-1)*(B+sum(4**j-2 for j in J))
def main():
 single=0
 for r in range(3,7):
  for p in range(1,2*r-3):
   a=trim(graph([0,1]*r,p));x=profile(a)
   ck(strong(a),'single strong')
   ck(x==(formula(r,p),3,2**(2*r-1)+1,3),'single strong state');single+=1
 A=mm(S[1],S[0]);F=sub(T[1],A,3);power=I
 ck(not any(map(any,mm(mm(mm(A,sub(A,I,4)),sub(A,I)),sub(A,I)))),'base annihilator')
 ck([Dmulti(r,()) for r in (1,2,3)]==[1,3,9],'base recurrence initials')
 annih=[]
 for h in range(5):
  z=mm(mm(F,power),F);annih.append(not any(map(any,z)));power=mm(A,power)
 ck(annih==[False,True,True,True,True],'separated-defect products')
 literal=transfer=0;maximum=0
 for r in range(3,22):
  slots=list(range(1,r-1,2))
  for t in range(len(slots)+1):
   ds=[]
   for J in combinations(slots,t):
    D=Dmulti(r,J);ck(D==predicted(r,J),'multi scalar');transfer+=1;ds.append(D)
    if r<=5:
     a=multi(r,J);x=profile(a);d=3**t;U=d*(4**r+2)//3
     ck(x==(D,0,U,d),'multi literal state')
     core=trim(a);ck(strong(core),'multi strong')
     ck(profile(core)==(D,d,U//2,d),'multi strong state');literal+=1
   ck(len(set(ds))==len(ds),'fixed t distinctness');maximum=max(maximum,len(ds))
 both=Dmulti(2,(0,1));naive=3*Dmulti(2,(0,))+3*Dmulti(2,(1,))-9*Dmulti(2,())
 ck((both,naive)==(41,39),'adjacent hostile')
 print(json.dumps({'status':'PASS','strong_single_literal_profiles':single,'multi_literal_and_strong_pairs':literal,'multi_transfer_cases':transfer,'maximum_fixed_t_fiber_checked':maximum,'r_max':21,'F_A_h_F_zero_h0_to4':annih,'adjacent_01_triangle_actual_and_naive':[both,naive],'scope':'Separated mobile modules only; no interval or Chan-Pak completeness.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
