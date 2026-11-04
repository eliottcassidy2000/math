#!/usr/bin/env python3
import audit as A
import physical_audit as B
from pathlib import Path
import subprocess,json,hashlib
P=Path(__file__).resolve().parent

def rooted_inflate(q,R):
    m=len(q);n=m+len(R)-1;out=list(q)+[0]*(len(R)-1)
    for i in range(m):
        for j in range(1,len(R)):
            if R[0]>>j&1:out[i]|=1<<(m+j-1)
            else:out[m+j-1]|=1<<i
    for i in range(1,len(R)):
        for j in range(1,len(R)):
            if R[i]>>j&1:out[m+i-1]|=1<<(m+j-1)
    return tuple(out)
def main():
    subprocess.run(['g++','-O2','-std=c++17',str(P/'ambient_walk.cpp'),'-o',str(P/'ambient_walk')],check=True)
    req=[];expect=[];saved=[]
    for q in [(2,0),(2,4,1),(6,12,24,17,3)]:
        m=len(q);e,o,E,_,_=A.banks(q,False);ec,oc,_,_,_=A.banks(A.opposite(q),False);d={k:e[k]+o[k]-ec[k]-oc[k] for k in e};L=max(6,A.val(E)+1);mod=1<<L
        i,j=next((i,j) for i in range(m) for j in range(m) if i!=j and A.mod(d[(1<<i)|(1<<j),1<<i],mod));cube=[(a,c) for a in [0,1<<i,1<<j,(1<<i)|(1<<j)] for c in [0,1<<i]]
        for R in [(2,4,1),(6,12,24,17,3)]:
            root2=rooted_inflate((2,0),R);inp='\n'.join(' '.join(map(str,[len(g),1,*g])) for g in [R,root2])+'\n';hh=list(map(int,subprocess.run([str(P/'physical')],input=inp,text=True,capture_output=True,check=True).stdout.split()));h1,h2=hh;A.require(h2%2==1 and (h2-h1)%2==0,'odd ambient factor');c2=(h2-h1)//2
            extmask=sum(1<<(m+j-1) for j in range(1,len(R)) if R[j] & 1);n=m+len(R)-1;start=len(req)
            for a,c in cube:
                for orient in [q,A.opposite(q)]:
                    graph=rooted_inflate(orient,R);req.append(' '.join(map(str,[m,n,L,a|extmask,c|extmask,*graph])))
                expect.append(A.mod(h2*A.evalbank(d,a,c),mod))
            saved.append({'quotient':q,'ambient_rooted_rows':R,'L':L,'c1':h1,'c2':c2,'odd_multiplier':h2,'cube':cube,'request_start':start})
    inp='\n'.join(req)+'\n';(P/'ambient_input.txt').write_text(inp);out=subprocess.run([str(P/'ambient_walk')],input=inp,text=True,capture_output=True,check=True).stdout;(P/'ambient_output.txt').write_text(out);values=[list(map(int,s.split())) for s in out.splitlines()]
    for index,expected in enumerate(expect):
        L=int(req[2*index].split()[2]);mod=1<<L;A.require((values[2*index][0]-values[2*index+1][0])%mod==expected,'direct exterior word ambient transfer')
    A.require(any(z['c2'] for z in saved),'nontrivial ambient second cover term')
    (P/'ambient_controls.json').write_text(json.dumps(saved,sort_keys=True,indent=2)+'\n')
    receipt={'status':'PASS','ambient_controls':len(saved),'ambient_word_queries':len(req),'nonzero_second_cover_coefficients':[z['c2'] for z in saved],'sources':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['ambient_audit.py','ambient_walk.cpp','audit.py','physical.cpp']},'artifacts':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['ambient_controls.json','ambient_input.txt','ambient_output.txt']}}
    (P/'ambient_receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
