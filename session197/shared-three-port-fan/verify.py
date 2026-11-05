from core import *
import json
C=[26,16,3,6,12];A=[90,16,3,6,12,31,62];B=[26,48,3,6,44,13];PORTS=[0,1,4]
def left_birth(rows,ports,mask):
 n=len(rows);selected={p for i,p in enumerate(ports)if mask>>i&1};out=rows[:]+[0]
 for v in range(n):
  if v in selected:out[v]|=1<<n
  else:out[n]|=1<<v
 # The actual candidate port list is old ports followed by n; drop index3.
 newports=ports+[n];del newports[3]
 require(newports==ports,'every birth retains the same three physical ports')
 return out,newports

def gate(rows,kind):
 ports=PORTS[:]
 for mask in ([0,1]if kind==2 else[6]):rows,ports=left_birth(rows,ports,mask)
 return rows,ports

def main():
 cores=[]
 for X,h,R in [(C,13,18),(A,41,54),(B,43,54)]:
  p=profile(X,2);require(p==literal_profile(X,2),'independent core profiles')
  require((p[0],p[2])==(h,R)and strong(X),'physical strong cores');cores.append({'rows':X,'profile':p})
 cases=[]
 for m in range(1,7):
  P,b=fan(C,2,m);X,px=gate(P,2);Y,py=gate(P,4)
  AX,ar=fan(A,2,m);BY,br=fan(B,2,m)
  amap=list(range(5))+[m+6,m+7,5]+list(range(6,m+6))
  bmap=list(range(5))+[m+6,5]+list(range(6,m+6))
  require(relabel(AX,amap)==X and relabel(BY,bmap)==Y,'all-edge graph identity controls')
  h=18*2**m+13
  ps=[profile(P,5)[:4],profile(X,5)[:4],profile(Y,5)[:4]]
  require(ps==[(h,13,h,13),(3*h+2,41,3*h+2,41),(3*h+4,43,3*h+4,43)],'actual old-root profiles')
  require(px==py==PORTS and all(strong(T)for T in(P,X,Y)),'same ports and strong outputs')
  cases.append({'m':m,'profiles_at_old_vertex5':ps,'ports':PORTS})
 P,_=fan(C,2,1);twice={}
 for kind,expected in [(2,[49,149,515]),(4,[49,151,409])]:
  rows=P;seq=[profile(rows,5)[0]]
  for _ in range(2):rows,_=gate(rows,kind);seq.append(profile(rows,5)[0])
  require(seq==expected,'non-returning scope control');twice[str(kind)]=seq
 return {'status':'PASS','cores':cores,'physical_fan_cases':cases,'naive_repetition':twice}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
