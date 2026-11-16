from mixed_controller import load,apply,dot,powm,need

def tower(v):
 a=3;R=4;q=5
 for j in range(2,v+1):
  q*=5
  a=next(a+z*R for z in range(5) if pow(2,a+z*R,q)==3)
  R*=5
 return a

def clock(D,k,target,phase,n=0):
 A=tuple(D['matrices'][D['words'].index([0,6,5])]);f=D['readout'];v=D['start']
 def H(t,m):return dot(f,apply(powm(A,t,m),v,m),m)
 candidates=[j for j in range(phase,60,4) if H(j,5)==target%5]
 need(len(candidates)==3)
 period=60;q=5
 for j in range(2,k+1):
  q*=5
  candidates=[next(t+r*period for r in range(5) if H(t+r*period,q)==target%q) for t in candidates]
  period*=5
 return min((t-n)%period or period for t in candidates),H

def main():
 need(tower(32)==18401310530325970670607)
 need([v for v in range(1,32) if ((z:=5**v+3)&(z-1)==0)]==[1,3])
 D=load()
 for k in (1,2):
  for phase in range(4):
   for a in range(5**k):
    w,H=clock(D,k,a,phase,31)
    need(w==next(j for j in range(1,60*5**(k-1)+1) if (31+j)%4==phase and H(31+j,5**k)==a))
 w,H=clock(D,3,1001,0,0)
 need(H(w,1000)==1 and H(w,32)==25)
 print("Tower exclusion through",tower(32)//3)
 print("Least waits checked at two small levels; 1001-compatible clock wait:",w,"eight-letter returns.")
if __name__=='__main__':main()
