"""Controls for a ternary actuator and a target-preserving root improvement.

The infinite statements are proved in PROOF.md. No coverage census.
"""
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def g(n):
    total=1
    while n>1:
        q=n//2
        if n%2==0: total+=q+(q&1)
        n=q
    return total

def rev(n):
    return int(bin(n)[2:][::-1],2) if n else 0

def L(B): return 13*g(B)+32*B+6*(B&1)+18
def N(C): return 3*g(C)+12*C+10
def D(B): return L(B)-N(3*B)
def F(q): return 13*g(q)-3*g(3*q)+q+13*(q&1)+20

def completed(q,m):
    power=16**m
    return 16*(q*power+(power-1)//5)+5

def prepare(p):
    def residue(P):return (g(P)+g(3*P+2))%4
    s=2 if residue(8*p+2)==0 else 4
    P=8*p+s
    q=32*P+31
    assert residue(P)==0 and F(q)<0 and F(q)%4==0
    return P,q,-F(q)//4

def qstate(s):
    H,E,U,d=s
    return (H,U,U-d,H+U-E)

def old_pair(B):
    u=45*(1<<B.bit_length())-13*rev(B)
    return [(L(B)+j,13,u,13) for j in (0,2)]

def new_pair(C):
    c=C.bit_length();C0=C-(1<<(c-1));eps=(C>>(c-2))&1
    return [(N(C)+j,N(C0)+j-3*eps,N(C)+j+3*eps,N(C0)+j) for j in (0,2)]

def construct(core,cut,B):
    right=[1-int(x) for x in bin(B)[2:][::-1]] if B else []
    n=len(core);old=n+len(right);rows=[0]*(old+1)
    for i in range(old):
        for j in range(i+1,old):
            win=bool(core[i]>>j&1) if j<n else True
            u,v=(i,j) if win else (j,i)
            rows[u]|=1<<v
    bits=[(cut>>i)&1 for i in range(n)]+right
    for i,x in enumerate(bits):rows[old if x else i]|=1<<(i if x else old)
    return rows

def end_dp(rows):
    n=len(rows);full=(1<<n)-1
    dp=[[0]*n for _ in range(full+1)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,full+1):
        for i,x in enumerate(dp[mask]):
            if not x:continue
            avail=rows[i]&~mask&full
            while avail:
                bit=avail&-avail;avail-=bit
                dp[mask|bit][bit.bit_length()-1]+=x
    return dp

def state(rows,root):
    n=len(rows);full=(1<<n)-1
    dp=end_dp(rows);h=[sum(v) for v in dp];h[0]=1
    converse=[sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n)]
    starts=end_dp(converse)
    U=sum(starts[mask][root]*h[full^mask] for mask in range(full+1))
    return h[full],dp[full][root],U,h[full^(1<<root)]

def main():
    for x in range(128):
        assert D(16*x+5)==F(x)
        assert F(16*x+3)==F(x)+4
        r2=(g(8*x+2)+g(3*(8*x+2)+2))%4
        r4=(g(8*x+4)+g(3*(8*x+4)+2))%4
        assert g(8*x+2)+g(3*(8*x+2)+2)==g(x)+g(3*x+1)+26*x+10
        assert g(8*x+4)+g(3*(8*x+4)+2)==g(x)+g(3*x+1)+18*x+12
        assert (r2-r4)%4==2
    prefixes=[]
    for p in range(32):
        P,q,m=prepare(p);B=completed(q,m)
        assert D(B)==0
        assert bin(B)[2:].startswith(bin(p)[2:]) if p else True
        old=old_pair(B);new=new_pair(3*B)
        for a,b in zip(old,new):
            diff=[x-y for x,y in zip(qstate(a),qstate(b))]
            assert diff[0]==0 and min(diff[1:])>0
        prefixes.append({'prefix':p,'prepared_prefix':P,'q':q,'m':m,'address_bits':B.bit_length(),'target_bits':L(B).bit_length(),'lower_q_savings_bits':[(x-y).bit_length() for x,y in zip(qstate(old[0]),qstate(new[0]))]})
    literal=[]
    for C in (2,3,9,21):
        expected=new_pair(C);actual=[]
        for cut,B in ((3,2*C+1),(1,4*C+3)):
            rows=construct((2,4,1),cut,B)
            actual.append(state(rows,len(rows)-2))
        assert actual==expected
        literal.append({'C':C,'rerooted_states':actual})
    repeated9=[]
    def M3(C):return 3*g(C)+36*C+30
    for C in (2,3):
        c=C.bit_length();C0=C-(1<<(c-1));eps=(C>>(c-2))&1
        expected=[(M3(C)+j,M3(C0)+j-3*eps,M3(C)+j+3*eps,M3(C0)+j) for j in (0,2)]
        actual=[]
        for cut,address in ((0,16*C+12),(1,8*C+6)):
            rows=construct((2,4,1),cut,address)
            actual.append(state(rows,len(rows)-2))
        assert actual==expected
        repeated9.append({'C':C,'rerooted_states':actual})
    B=3
    old_actual=[state(rows,len(rows)-1) for rows in (construct((18,8,3,5,14),cut,B) for cut in (6,3))]
    assert old_actual==old_pair(B)
    assert D(B)==0
    explicit=[]
    for q,m in ((31,1),(135,7)):
        B=completed(q,m);assert F(q)==-4*m and D(B)==0
        explicit.append({'q':q,'m':m,'B':B,'C':3*B,'targets':[L(B),L(B)+2],'old_pair':old_pair(B),'rerooted_C3_pair':new_pair(3*B)})
    for k in range(7,15):
        q=(1<<k)+7;m=15*(1<<(k-7))-8
        assert g(q)==(1<<(k-3))+1
        assert g(3*q)==33*(1<<(k-5))+7
        assert F(q)==-4*m
    result={'scope':'Bounded controls for proved actuator and reroot formulas, not a coverage scan.','actuator_and_preparation_inputs':128,'completed_prefixes':prefixes,'literal_C3_controls':literal,'repeated9_family_reroot_controls':repeated9,'literal_H13_at_B3':old_actual,'explicit_completed_sources':explicit,'beneficial_family_exponents':[7,14]}
    (HERE/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'pass','prefixes':len(prefixes),'literal_C3_controls':len(literal),'explicit_completed_sources':explicit},indent=2))

if __name__=='__main__':main()
