"""ozj.py - unit-offset states in the O/Z/J (D/K/A of Session195) grammar from the singleton.
State (H, d, E, U) of a rooted tournament (L, v):
  Z: new vertex beating only v, root stays       -> (H+U-E, d, d, 2U-d)
  O: new sink (beaten by all), root stays        -> (H, d, 0, 2U)
  J: new vertex beating only v, root moves to it -> (H-E+U, H, H-E, H+U)
offset c0 = H - E + 2d - U, slope M = U - d. Breadth-first over all words with U <= UMAX.
Prints: number of states, the odd parts c of M among unit-offset states, the odd c <= CMAX not obtained,
and a shortest word for each small c."""
import sys
UMAX=int(sys.argv[1]) if len(sys.argv)>1 else 1<<14
CMAX=int(sys.argv[2]) if len(sys.argv)>2 else 257
start=(1,1,1,1)
seen={start:""}; frontier=[start]
def nxt(s):
    H,d,E,U=s
    yield "Z",(H+U-E,d,d,2*U-d)
    yield "O",(H,d,0,2*U)
    yield "J",(H-E+U,H,H-E,H+U)
while frontier:
    new=[]
    for s in frontier:
        for l,t in nxt(s):
            if t[3]<=UMAX and t not in seen:
                seen[t]=seen[s]+l; new.append(t)
    frontier=new
unit={}
for s,wd in seen.items():
    H,d,E,U=s
    if H-E+2*d-U==1 and U-d>0:
        M=U-d; c=M; j=0
        while c%2==0: c//=2; j+=1
        if c not in unit or (j,len(wd))<(unit[c][0],len(unit[c][2])): unit[c]=(j,s,wd)
print("states with U <= %d: %d; unit-offset states: %d; distinct odd c: %d"%(UMAX,len(seen),sum(1 for s in seen if s[0]-s[2]+2*s[1]-s[3]==1),len(unit)))
print("odd c <= %d without a unit-offset O/Z/J state (U <= %d):"%(CMAX,UMAX),[c for c in range(1,CMAX+1,2) if c not in unit])
for c in range(1,64,2):
    if c in unit: print(c,"j=%d"%unit[c][0],"state",unit[c][1],"word",unit[c][2])
