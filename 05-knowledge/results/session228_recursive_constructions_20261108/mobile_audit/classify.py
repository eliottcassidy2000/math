from verify import *

def evalue(w,p):return 3-(2*w[-2]+w[-1]) if p==len(w)-2 else 3*(1-w[-1])
def fiber(M):
 w=list(map(int,'00'+bin(M)[2:]));return {(count(w,p),evalue(w,p)) for p in range(len(w)-1)}
def single(M):
 w=list(map(int,bin(M)[2:]));x=v
 for b in w:x=mv(S[b],x)
 return sum(a*b for a,b in zip(ell,x)),1-w[-1]
def main():
 total=d1=d3=0
 for n in range(1,6):
  pairs=[(i,j) for i in range(n) for j in range(i)]
  for mask in range(1<<len(pairs)):
   a=[[False]*n for i in range(n)]
   for z,(i,j) in enumerate(pairs):a[i][j]=bool(mask>>z&1);a[j][i]=not a[i][j]
   D,E,U,d=profile(a);total+=1
   if d==1:ck((D,E)==single(U-1),'d1 membership');d1+=1
   if d==3:ck(U%3==0 and (D,E) in fiber(U//3-1),'d3 membership');d3+=1
 padded=0
 for M in range(128):
  target=fiber(M)
  for zeros in range(2,8):
   w=list(map(int,'0'*zeros+bin(M)[2:]));got={(count(w,p),evalue(w,p)) for p in range(len(w)-1)}
   ck(got==target,'leading zero stabilization');padded+=len(w)-1
 print(json.dumps({'status':'PASS','labelled_root_zero_tournaments_orders_1_to_5':total,'d1_profiles':d1,'d3_profiles':d3,'leading_zero_M_range':[0,127],'leading_zero_counts':[2,7],'padded_placements':padded,'d3_M0':sorted(fiber(0)),'d3_M1':sorted(fiber(1))},indent=2,sort_keys=True))
if __name__=='__main__':main()
