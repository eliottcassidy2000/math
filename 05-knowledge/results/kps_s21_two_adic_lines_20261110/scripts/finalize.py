"""finalize.py - line table from seeds_m.txt (tail seeds, orders 2..9) and dbl_m.txt (module seeds, orders 3..8),
and the status of every point 1 + c 2^s below the first ray using the order <= 9 value list (cls files, field 2)."""
import os
S9=set()
for fn,m in [("cls_%d.txt"%m,m) for m in range(1,9)]+[("cls9_9.txt",9)]:
    for line in open(fn): S9.add(int(line.split()[1]))
best={}
for m in range(2,10):
    for line in open("seeds_%d.txt"%m):
        mm,code,v,H,d,E,U,M=map(int,line.split())
        if M<=0: continue
        c=M;j=0
        while c%2==0: c//=2;j+=1
        key=(j+1,mm)
        if c not in best or key<best[c][0]: best[c]=(key,"tail",mm,code,v,H,d,E,U)
for m in range(3,9):
    for line in open("dbl_%d.txt"%m):
        code,H,x,c2,fl=line.split()
        if fl!="1" or int(c2)==0: continue
        c=int(c2);j=0
        while c%2==0: c//=2;j+=1
        key=(j+1,m)
        if c not in best or key<best[c][0]: best[c]=(key,"module",m,int(code),int(x),int(H),0,0,0)
out=open("line_seeds.txt","w")
out.write("# c s0 mechanism order code vertex H d E U   (ray: 1 + c 2^s in the spectrum for all s >= s0)\n")
for c in sorted(best):
    k,mech,m,code,v,H,d,E,U=best[c]; out.write("%d %d %s %d %d %d %d %d %d %d\n"%(c,k[0],mech,m,code,v,H,d,E,U))
out.close()
print("odd c with a seed (order <= 9 tail, <= 8 module):",len(best))
full=[];partial=[];none=[]
for c in range(1,130,2):
    if c not in best: none.append(c); continue
    s0=best[c][0][0]; low=[1+c*2**s for s in range(1,s0)]
    miss=[n for n in low if n not in S9]
    (full if not miss else partial).append((c,s0,miss))
print("c <= 129 with the whole line in the spectrum (seed + order<=9 values below it):",[c for c,_,_ in full])
print("c <= 129 with a seed but a point below it not found at order <= 9:",partial)
print("c <= 129 without a seed:",none)
print("gaps on lines: 7 in S9?",7 in S9,"21 in S9?",21 in S9, "| max value at order 9:",max(S9),"| odd values < 1000 missing at order<=9:",[n for n in range(1,1000,2) if n not in S9][:20])
print("weight-three rows (c = 2^d + 1):",{d:(best[2**d+1][0][0] if 2**d+1 in best else None) for d in range(1,9)})
