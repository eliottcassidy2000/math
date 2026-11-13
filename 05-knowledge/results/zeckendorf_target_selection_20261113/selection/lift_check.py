"""Exact run coefficients, literal neutral repair, and guarded target fixtures."""
from pathlib import Path
from itertools import permutations, product
from math import gcd,isqrt
import json
from response_collision import profiles,moments
from module_probe import S,T,ROOT,inflate,tt

REF='ed79a99437c3bf19509b2f9b7b696ab9349032d3'
CORE=(6,28,24,17,1)
C=((0,1,0),(1,1,1),(1,0,0))
D=((1,0,1),(0,2,0),(0,0,0))
SHORT='D'*8+'CD'
LONG='CCCCD'+'CD'*3


def need(ok,msg):
    if not ok:raise RuntimeError(msg)


def ordinal(*blocks):
    n=sum(map(len,blocks));out=[];offset=0
    for B in blocks:
        later=((1<<n)-1)^((1<<(offset+len(B)))-1)
        out.extend((r<<offset)|later for r in B);offset+=len(B)
    return tuple(out)


def strong(rows):
    for v in range(len(rows)):
        reached=1<<v;old=0
        while old!=reached:
            old=reached
            for j,row in enumerate(rows):
                if reached>>j&1:reached|=row
        if reached!=(1<<len(rows))-1:return False
    return True


def birth(rows,r,move=False):
    n=len(rows);out=tuple(row|((1<<n)if i==r else 0)for i,row in enumerate(rows))+(((1<<n)-1)^(1<<r),)
    return out,n if move else r


def run(rows,r,word):
    for c in word:rows,r=birth(rows,r,c=='C')
    return rows,r


def quotient_coeff(rows,k,obs):
    outside=[i for i in range(len(rows))if i!=0 and not(obs=='q' and i==ROOT)]
    # Distinguishable temporary zeros are quotiented here, not in C_k(M).
    words=set(permutations(outside+[0]*k));count=0
    for w in words:
        if obs=='s' and w[0]!=ROOT:continue
        cut=w.index(ROOT)+1 if obs=='F' else len(w)
        if all(i==cut-1 or (w[i]!=w[i+1] and rows[w[i]]>>w[i+1]&1)for i in range(len(w)-1)):
            count+=1
    return count


def mm(A,B):return tuple(tuple(sum(A[i][k]*B[k][j]for k in range(3))for j in range(3))for i in range(3))
def matrix(word):
    M=((1,0,0),(0,1,0),(0,0,1))
    for c in word:M=mm(C if c=='C' else D,M)
    return M


def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))


def main():
    need(strong(S)and strong(T),'strong skeletons')
    middle=list(S);middle[1]^=1<<2;middle[2]^=1<<1
    need(strong(middle),'strong intermediate flip')
    coeff={obs:[quotient_coeff(S,k,obs)for k in (1,2)]for obs in ('H','s','q','F')}
    coeff2={obs:[quotient_coeff(T,k,obs)for k in (1,2)]for obs in ('H','s','q','F')}
    expected={'H':[29,15],'s':[5,0],'q':[5,0],'F':[35,15]}
    need(coeff==coeff2==expected,'complete run coefficient identity')
    middle_coeff={obs:[quotient_coeff(middle,k,obs)for k in (1,2)]for obs in ('H','s','q','F')}
    need(middle_coeff=={'H':[17,9],'s':[3,0],'q':[3,0],'F':[21,9]},'complete intermediate descent coefficients')
    controls=[]
    for name,M in [('TT1',tt(1)),('TT2',tt(2)),('TT3',tt(3)),('C3',(2,4,1)),('H5',(2,4,9,3)),('H11',CORE),('H11_prefix',ordinal(tt(1),CORE))]:
        H,_=moments(M);full=len(H)-1;h=H[-1];p=sum(H[a]*H[full^a]for a in range(full+1))
        predicted=(15*p-h,5*h,5*h,15*p+5*h)
        A,r=inflate(S,0,M);B,q=inflate(T,0,M)
        need(r==q and profiles(A)[r]==profiles(B)[q]==predicted,'literal full-state law')
        need(strong(A)and strong(B),'strong arbitrary-module witnesses')
        m=len(M);changes=[(i,j)for i in range(len(A))for j in range(i+1,len(A))if((A[i]>>j)&1)!=((B[i]>>j)&1)]
        need(changes==[(m,m+1)],'exactly one old arc changes')
        need(B[-1]==1<<m,'fresh vertex beats only old1')
        need(all(((A[r]>>j)&1)==((B[r]>>j)&1)for j in range(len(A))),'retained old root star')
        content=gcd(gcd(predicted[0],predicted[1]),gcd(predicted[2],predicted[3]))
        need(content==gcd(h,15*p),'exact content')
        # Literal doubling commutation with the named module prefix.
        for skeleton,G in ((S,A),(T,B)):
            grown,_=birth(G,r)
            compiled,rr=inflate(skeleton,0,ordinal(tt(1),M))
            mapping=[i+1 for i in range(len(G))]+[0]
            need(mapping[r]==rr and all(((grown[i]>>j)&1)==((compiled[mapping[i]]>>mapping[j])&1)for i in range(len(grown))for j in range(len(grown))),'actual prefix-growth embedding')
        controls.append(dict(module=name,rows=M,h=h,p=p,state=predicted,content=content,source_rows=A,target_rows=B,root=r))
    H,_=moments(CORE);full=len(H)-1
    need(H[-1]==11 and sum(H[a]*H[full^a]for a in range(full+1))==100,'unbounded-core seed profile')
    score_twins=[]
    for mask,rows,expected in [(871,(14,24,18,20,1),(15,108)),(733,(26,20,9,18,4),(11,100))]:
        H,_=moments(rows);full=len(H)-1
        hp=(H[-1],sum(H[a]*H[full^a]for a in range(full+1)))
        scores=sorted(row.bit_count()for row in rows)
        need(scores==[1,2,2,2,3] and hp==expected,'same-score response hostile')
        h,p=hp
        score_twins.append(dict(mask=mask,rows=rows,scores=scores,h=h,p=p,target=15*p-h,content=gcd(h,15*p)))
    need([(x['target'],x['content'])for x in score_twins]==[(1605,15),(1489,1)],'score quotient loses primitive entry')
    phases=[]
    for n in range(9):
        for m in range(14):
            if n+m==0:continue
            y=15*100**n*2**m-11**n
            need((y%7==0)==((m-n)%3==0),'seven selector')
            if (m-n)%3==0:
                k=(m-n)//3
                need((y%49!=0)==((2+5*n+k)%7!=0),'exact seven valuation selector')
            phases.append(dict(n=n,m=m,H=y))
    fixtures=[]
    for n,m in [(0,3),(0,6),(1,0),(1,7),(2,5)]:
        y=15*100**n*2**m-11**n
        if y%7==0:need(prime(y//7),'exact7p fixture')
        else:need(prime(y),'exact prime fixture')
        fixtures.append(dict(n=n,m=m,H=y,prime=y%7!=0,seven_cofactor=y//7 if y%7==0 else None,orders=[5*n+m+5,5*n+m+6]))
    # The other survivor: exact short/long guard relation, with its port boundary.
    need(len(SHORT)==10 and len(LONG)==11,'distinct native lengths')
    need(matrix('DCD'+SHORT)==matrix('DCD'+LONG),'universal guarded operator identity')
    word_controls=[]
    for word in (SHORT,LONG):
        G,r=run((2,4,1),0,word);z=profiles(G)[r]
        need(z==(1539,513,513,1539),'literal canonical return')
        child=len(G)-1;other=profiles(G)[child]
        sealed,rr=birth(G,r)
        need(profiles(sealed)[rr]==(2565,513,513,2565),'fresh-seal rooted equality')
        word_controls.append(dict(word=word,order=len(G),state=z,last_child_state=other,sealed_order=len(sealed)))
    need(word_controls[0]['last_child_state']!=word_controls[1]['last_child_state'],'hostile reused-last-child distinction')
    Path(__file__).with_name('lift_receipt.json').write_text(json.dumps(dict(status='PASS',ref=REF,scope='Complete at-most-two-run coefficients; selected literal modules and guard words; finite arithmetic controls, not target coverage',coefficients=coeff,intermediate_coefficients=middle_coeff,literal_controls=controls,score_twins=score_twins,phase_controls=phases,prime_fixtures=fixtures,short_long=word_controls),indent=2)+'\n')
    print('PASS universal one-arc module repair, exact target phases, strong short/long controls')
if __name__=='__main__':main()
