"""Unrestricted occupancy routing for the incoming twelve actual owner equations."""
from route import *
import random
NEWREV='2bc5d106db075752dd28090cb5637555ede01d66'
SOURCE='05-knowledge/results/guarded_component_completion_20261024/native/receipt.json'

class TwoRoots:
    def __init__(self,remove=False):
        self.R=Router(remove,True);R=self.R
        raw=subprocess.check_output(['git','show',NEWREV+':'+SOURCE],cwd=ROOT)
        self.hash=hashlib.sha256(raw).hexdigest()
        self.F=next(f for f in json.loads(raw)['families'] if f['control_count']==12)
        self.guards=self.F['guards'];self.controls={R.ids.index(g['helper']) for g in self.guards}
        self.fixed=R.ids.index(2948);self.z=R.pin;self.t=R.ids.index(1532)
        self.free=[i for i in R.E if i not in self.controls and i!=self.fixed]
        self.fm=sum(1<<i for i in self.free);self.ordinary=R.rm|self.fm;self.size=20+len(self.free)
        need(len(self.controls)==12 and len(self.free) in (47,48),'complete declared twelve-pin alphabet')
        E=self.F['embeddings'];private=set(E['2393'])
        need(all(not private.intersection(v) for s,v in E.items() if s!='2393'),'new570 option stays in a private old cell')
    def controls_at(self,J,z):
        bits={176:z,1532:int(bool(J>>self.t&1))}
        return sum(1<<self.R.ids.index(g['helper']) for g in self.guards if bits[g['root']]^g['bit_xor'])
    def valid(self,m):
        return not m>>self.fixed&1 and m&~self.ordinary==self.controls_at(m&self.fm,int(bool(m>>self.z&1)))
    def fill(self,J,r,q,z):
        C=self.controls_at(J,z);g=sub(q,self.R.q(J|C))
        return J|C|self.R.conditional[z][r][g]
    def normalize(self,m,k,q):
        R=self.R;P=[m];J=m&self.fm;r=(m&R.rm).bit_count();z=int(bool(m>>self.z&1))
        if r<7:
            for i in [i for i in self.free if J>>i&1][:7-r]:J^=1<<i
            r=7;P.append(self.fill(J,r,q,z))
        elif r>13:
            for i in [i for i in self.free if not J>>i&1][:r-13]:J|=1<<i
            r=13;P.append(self.fill(J,r,q,z))
        target=7 if k<=len(self.free)+7 else 13
        while r!=target:
            if r>target:
                i=next(i for i in self.free if not J>>i&1);J|=1<<i;r-=1
            else:
                i=next(i for i in self.free if J>>i&1);J^=1<<i;r+=1
            P.append(self.fill(J,r,q,z))
        if z:P.append(self.fill(J,r,q,0))
        return P
    def path(self,a,b):
        R=self.R;k=(a&self.ordinary).bit_count();q=R.q(a)
        need(self.valid(a) and self.valid(b) and a.bit_count()==b.bit_count() and R.q(b)==q,'compatible full endpoints')
        need(a.bit_count()==k+6,'six actual controls at V, not free slots')
        if min(k,self.size-k)<7:P=[a,b]
        else:
            A=self.normalize(a,k,q);B=self.normalize(b,k,q)
            J=A[-1]&self.fm;K=B[-1]&self.fm;r=(A[-1]&R.rm).bit_count()
            drop=[i for i in self.free if J>>i&1 and not K>>i&1]
            take=[i for i in self.free if K>>i&1 and not J>>i&1]
            for i,j in zip(drop,take):
                J^=(1<<i)|(1<<j);A.append(B[-1] if J==K else self.fill(J,r,q,0))
            P=A+B[::-1]
        P=[m for j,m in enumerate(P) if not j or m!=P[j-1]]
        for m in P:need(self.valid(m) and m.bit_count()==a.bit_count() and R.q(m)==q,'all shared controls and current rows')
        for x,y in zip(P,P[1:]):
            need((x^y).bit_count()<=24 and R.work(x,y)<=156,'proved unrestricted twelve-pin budget')
            if min(k,self.size-k)>=7:need(R.work(x,y)<=89,'interior twelve-pin bound')
        need(len(P)-1<=40,'binary event count')
        return P

def main():
    rng=random.Random(12570);reports=[];checks=0;example=None
    for remove in (False,True):
        C=TwoRoots(remove);R=C.R;cases=events=maxwork=maxsupport=0;roots_seen=set()
        paths={s:{o['target']:o['path'] for o in row['options']} for s,row in zip(R.ids,R.rows)}
        def literal(m):
            nonlocal checks
            for g in C.guards:
                s,h=g['root'],g['helper'];a=paths[s][2234 if m>>R.ids.index(s)&1 else 1965]
                b=paths[h][2234 if m>>R.ids.index(h)&1 else 1965]
                need(C.F['embeddings'][str(s)][a[g['root_clock']]]==C.F['embeddings'][str(h)][b[g['helper_clock']]],'literal shared-owner value')
                checks+=1
        for k in range(C.size+1):
            for r in range(max(0,k-len(C.free)),min(20,k)+1):
                z=int(r>0 and (r==20 or rng.randrange(2)))
                I=rng.sample([i for i in R.R if i!=C.z],r-z)+([C.z] if z else [])
                J=sum(1<<i for i in rng.sample(C.free,k-r))
                a=sum(1<<i for i in I)|J|C.controls_at(J,z);q=R.q(a)
                if min(k,C.size-k)<7:
                    # Search actual compatible alternatives, retaining identity if fibre is isolated in this probe.
                    b=a
                    for _ in range(40):
                        rr=min(k,6) if k<7 else max(14,k-len(C.free));zz=rng.randrange(2)
                        if rr not in R.conditional[zz]:continue
                        JJ=sum(1<<i for i in rng.sample(C.free,k-rr));CC=C.controls_at(JJ,zz)
                        g=sub(q,R.q(JJ|CC))
                        if g in R.conditional[zz][rr]:b=JJ|CC|R.conditional[zz][rr][g];break
                else:
                    rr=7 if k<=len(C.free)+7 else 13;JJ=sum(1<<i for i in rng.sample(C.free,k-rr))
                    b=C.fill(JJ,rr,q,rng.randrange(2))
                P=C.path(a,b);cases+=1;events+=len(P)-1
                for m in P:
                    literal(m);roots_seen.add((int(bool(m>>C.z&1)),int(bool(m>>C.t&1))))
                for x,y in zip(P,P[1:]):maxwork=max(maxwork,R.work(x,y));maxsupport=max(maxsupport,(x^y).bit_count())
                if example is None and len(P)>8 and len({int(bool(m>>C.z&1)) for m in P})==2:
                    example={'sources':R.ids,'masks':P,'events':[{'support':(x^y).bit_count(),'clocks':R.work(x,y)} for x,y in zip(P,P[1:])]}
        need(roots_seen=={(0,0),(0,1),(1,0),(1,1)},'all four actual shared-root phases occur')
        reports.append({'binary_sources':R.n,'ordinary_variables':C.size,'endpoint_pairs':cases,'events':events,
                        'max_observed_work':maxwork,'max_observed_support':maxsupport,'root_phases':sorted(roots_seen)})
    return {'status':'PASS literal twelve-pin unrestricted routing checks','checks':checks,'cases':reports,'example':example,
            'incoming_revision':NEWREV,'incoming_path':SOURCE,'incoming_sha256':C.hash,
            'scope':'same twelve owner equations and fixed2948 U; core/helper margins are no longer separately fixed'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
