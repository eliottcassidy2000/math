#!/usr/bin/env python3
"""Independent exterior-gap certificate; no packages or network required."""
import itertools, json, math, pathlib, hashlib
HERE=pathlib.Path(__file__).resolve().parent
def need(test,msg):
    if not test: raise RuntimeError(msg)
def rows(m):
    n=m+8; r=[0]*n; r[:4]=[2,4,9,3]; w=4; A=m+5;B=m+6;C=m+7
    def edge(a,b):r[a]|=1<<b
    for v in range(4):edge(v,w) if v==3 else edge(w,v)
    for z in range(5,m+5):
        edge(w,z)
        for v in range(4):edge(z,v)
        for z2 in range(z+1,m+5):edge(z,z2)
    for new,root in [(A,w),(B,A),(C,B)]:
        for v in range(new):edge(v,new) if v==root else edge(new,v)
    return r
def strong(r):
    n=len(r);full=(1<<n)-1
    for v in range(n):
        seen=1<<v
        while True:
            new=seen
            for u in range(n):
                if seen>>u&1:new|=r[u]
            if new==seen:break
            seen=new
        if seen!=full:return False
    return True
def coefficients():
    r=rows(1);z=5;ext=tuple(v for v in range(9)if v!=z)
    port=(0,1,2);w=4;C=[[0]*11 for _ in range(48)]
    def bin_data(seq):
        if not seq:return 1,0,None,None
        start=int(bool(r[z]>>seq[0]&1));end=int(bool(r[seq[-1]]>>z&1));k=start+end;req=0
        for a,b in zip(seq,seq[1:]):
            eligible=int(bool(r[a]>>z&1 and r[z]>>b&1))
            if not(r[a]>>b&1):
                if not eligible:return None
                req+=1
            k+=eligible
        return k,req,start,end
    def add(j,k,req):
        for a in range(req+1):C[j][k-a]+=(-1)**a*math.comb(req,a)
    for p in itertools.permutations(ext):
        d=bin_data(p)
        if d:
            k,req,st,en=d
            for i,v in enumerate(port):
                off=16*i;add(off,k,req)
                if p[-1]==v:add(off+2,k-en,req)
                if p[0]==v:add(off+3,k-st,req)
                if p[-1]==w:add(off+7,k-en,req)
                if p[0]==w:add(off+8,k-st,req)
        for cut in range(9):
            l=p[:cut];rr=p[cut:];dl=bin_data(l);dr=bin_data(rr)
            if dl is None or dr is None:continue
            k=dl[0]+dr[0];req=dl[1]+dr[1]
            for i,v in enumerate(port):
                off=16*i;add(off+1,k,req)
                for u,ei,si in [(v,off+5,off+6),(w,off+10,off+11)]:
                    if l and l[-1]==u:add(ei,k-dl[3],req)
                    if l and l[0]==u:add(si,k-dl[2],req)
                if l and rr:
                    if l[-1]==v and rr[0]==w:add(off+13,k-dl[3]-dr[2],req)
                    if l[-1]==w and rr[0]==v:add(off+14,k-dl[3]-dr[2],req)
    for u in (0,1,2,4):
        co=[0]*11
        for p in itertools.permutations([v for v in ext if v!=u]):
            d=bin_data(p)
            if d:
                k,req,_,_=d
                for a in range(req+1):co[k-a]+=(-1)**a*math.comb(req,a)
        if u==4:
            for i in range(3):C[16*i+9]=co[:]
        else:C[16*u+4]=co
    return C
def evaluate(C,m):return[sum(c*b**m for b,c in enumerate(row))for row in C]
def direct(r):
    # A third small-order route: enumerate every Hamilton word and every cut.
    n=len(r);out=[0]*48;port=(0,1,2);w=4;deletion={u:0 for u in(*port,w)}
    def path(p):return all(r[a]>>b&1 for a,b in zip(p,p[1:]))
    for p in itertools.permutations(range(n)):
        if path(p):
            for i,v in enumerate(port):
                off=16*i;out[off]+=1;out[off+2]+=p[-1]==v;out[off+3]+=p[0]==v;out[off+7]+=p[-1]==w;out[off+8]+=p[0]==w
        for cut in range(n+1):
            l=p[:cut];rr=p[cut:]
            if not(path(l)and path(rr)):continue
            for i,v in enumerate(port):
                off=16*i;out[off+1]+=1
                if l:
                    out[off+5]+=l[-1]==v;out[off+6]+=l[0]==v;out[off+10]+=l[-1]==w;out[off+11]+=l[0]==w
                if l and rr:out[off+13]+=l[-1]==v and rr[0]==w;out[off+14]+=l[-1]==w and rr[0]==v
    for u in deletion:deletion[u]=sum(path(p)for p in itertools.permutations([v for v in range(n)if v!=u]))
    for i,v in enumerate(port):out[16*i+4]=deletion[v];out[16*i+9]=deletion[w]
    return out
def main():
    saved=json.loads((HERE/'response.json').read_text());C=coefficients()
    need([json.loads(line)for line in(HERE/'samples.jsonl').read_text().splitlines()]==saved['samples'],'Regenerated subset-DP samples mismatch')
    need(all(row[1:]==saved['coefficients'][i]for i,row in enumerate(C)),'Complete power coefficient table differs')
    for sample in saved['samples']:
        m=sample['m'];r=rows(m)
        need(r==sample['rows'],'Literal rows mismatch')
        need(strong(r),'Not strong')
        need(evaluate(C,m)==sample['fields'],'Independent gap/DP mismatch')
        need([i for i,s in enumerate(sample['starts'])if s==1]==[0,1,2],'Unit-start ports mismatch')
        if m>=1:need(sample['fields'][0]==25+36*2**m,'Hamilton family formula')
    for m in (0,1):need(direct(rows(m))==evaluate(C,m),'Independent permutation/cut mismatch')
    transport=[]
    for m in range(1,10):
        r=rows(m);n=len(r);born=r+[15]
        for v in range(4,n):born[v]|=1<<n
        order=list(range(m+5))+[n]+list(range(m+5,n))
        transported=[sum(1<<j for j,b in enumerate(order)if born[a]>>b&1)for a in order]
        need(transported==rows(m+1),'Physical bulk birth mismatch')
        transport.append({'m':m,'new_out_mask':15,'canonical_to_physical_indices':order})
    # Same literal component, all three legal restored ports, every clock parity.
    for m in range(1,20):
        f=evaluate(C,m);h=f[0]
        need(h%3==1,'Lost nonzero component residue')
        for k in range(1,10):
            wrapped=[h+2**k*f[16*i+5]for i in range(3)]
            need(wrapped[0]%3==1,'First restored port residue')
            need(sorted(x%3 for x in wrapped[1:])==[0,2],'Other restored ports residue')
    # Complete integer observation lattice generators, checked as row identities.
    def combo(terms):return[sum(a*C[j][b]for a,j in terms)for b in range(11)]
    S=C[3]
    D=combo([(11,34),(1,21),(-1,5),(-73,3)])
    need(D[1:]==[0,2]+[0]*8,'2u2 lattice generator')
    K=[2*C[37][b]-C[5][b]+45*S[b]-16*D[b]for b in range(11)]
    need(K[1:]==[0,0,38]+[0]*7,'38u3 lattice generator')
    need(all(row[2]%2==0 and row[3]%38==0 for row in C),'Row lattice inclusion')
    for m in range(1,30):
        f=evaluate(C,m);fa,fb,fc=f[5],f[21],f[37];delta=fa-fb
        need(0<fc<fb<fa<2*fc,'Dyadic response ordering')
        need(fb-fc>delta and 2*fc-fa>delta,'Sharp spacing comparison')
        need(2*delta==128*2**m-36,'Sharp spacing identity')
    receipt={'status':'PASS_INDEPENDENT_GAP_AND_PERMUTATION_ENUMERATION','symbolic_coefficients_including_zero_base':C,'sample_m':list(range(11)),'complete_field_count':48,'active_positive_bases':[1,2,3],'independent_permutation_orders':[8,9],'symbolic_evaluation_m20_fields':evaluate(C,20),'bulk_transport':transport,'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in(HERE/'count.cpp',HERE/'verify.py')}}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('PASS: all 48 symbolic fields; all DP samples m=0..10; independent permutations at orders 8,9; literal old-name bulk transport; both nonzero wrapper residues.')
if __name__=='__main__':main()
