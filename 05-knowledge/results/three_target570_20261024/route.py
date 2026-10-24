"""Construct full-margin routes, with no frozen core/helper occupancy."""
from bank import *
import argparse

class Router:
    def __init__(self,remove2393=False,moving=False):
        self.rows=[r for r in json.loads((HERE/'menus.json').read_text())['rows']
                   if r['role']!='fresh' and not(remove2393 and r['source']==2393)]
        self.ids=[r['source'] for r in self.rows];self.n=len(self.ids)
        self.labels=list(map(code,self.rows));self.costs=list(map(cost,self.rows))
        b=json.loads((HERE/'bank.json').read_text());S=set(b['reservoir_sources'])
        self.R=[i for i,s in enumerate(self.ids) if s in S];self.E=[i for i in range(self.n) if i not in self.R]
        self.rm=sum(1<<i for i in self.R);self.em=sum(1<<i for i in self.E);self.full=self.rm|self.em
        self.D=b['sections'];self.low=sections(self.labels,self.costs,4)
        self.moving=moving;self.lo=7 if moving else 5;self.hi=20-self.lo
        if moving:
            self.pin=self.ids.index(176);self.partner=self.ids.index(368)
            RR=[i for i in self.R if i!=self.pin]
            D=sections([self.labels[i] for i in RR]);self.conditional=[{},{}]
            for bit in (0,1):
                for k,d in enumerate(D):
                    self.conditional[bit][k+bit]={}
                    for g,(_,m) in d.items():
                        h=add(g,self.labels[self.pin]) if bit else g
                        mask=sum(1<<i for j,i in enumerate(RR) if m>>j&1)|(bit<<self.pin)
                        self.conditional[bit][k+bit][h]=mask
            need(all(len(self.conditional[b][r])==2048 for b in (0,1) for r in range(7,14)),
                 'all conditional grades under either actual moving-owner value')
            idx=[i for i in range(self.n) if i not in (self.pin,self.partner)]
            small=sections([self.labels[i] for i in idx],[self.costs[i] for i in idx],5)
            self.low=[{} for _ in range(7)]
            for k in range(1,7):
                for forced in (self.pin,self.partner):
                    for g,(c,m) in small[k-1].items():
                        h=add(g,self.labels[forced]);mask=(1<<forced)|sum(1<<i for j,i in enumerate(idx) if m>>j&1)
                        v=(c+self.costs[forced],mask)
                        if h not in self.low[k] or v<self.low[k][h]:self.low[k][h]=v
            need(all(7*k+max(c for c,m in self.low[k].values())<=55 for k in range(1,7)),
                 'complete compatible boundary witness budgets')
        need(len(self.R)==20 and len(self.rows) in (80,81),'same authentic reservoir')
        need(max(c for c,m in self.low[4].values())<=16,'cheap boundary-four witnesses')
    def q(self,m):return charge(m,self.labels)
    def compatible(self,m):return not self.moving or bool(m>>self.pin&1)!=bool(m>>self.partner&1)
    def work(self,a,b):return sum(c for i,c in enumerate(self.costs) if (a^b)>>i&1)
    def replace(self,outside,r,q):
        if self.moving:
            forced=1-int(bool(outside>>self.partner&1))
            return outside|self.conditional[forced][r][sub(q,self.q(outside))]
        m=self.D[r][str(sub(q,self.q(outside)))]
        return outside|sum(1<<i for j,i in enumerate(self.R) if m>>j&1)
    def normalize(self,start,k,q):
        ans=[start];r=(start&self.rm).bit_count();outside=start&self.em
        if r<self.lo:
            for i in [i for i in self.E if outside>>i&1][:self.lo-r]:outside^=1<<i
            r=self.lo;ans.append(self.replace(outside,r,q))
        elif r>self.hi:
            for i in [i for i in self.E if not outside>>i&1][:r-self.hi]:outside|=1<<i
            r=self.hi;ans.append(self.replace(outside,r,q))
        target=max(self.lo,k-len(self.E))
        while r>target:
            i=next(i for i in self.E if not outside>>i&1)
            outside|=1<<i;r-=1;ans.append(self.replace(outside,r,q))
        return ans
    def path(self,start,finish):
        need(0<=start<=self.full and 0<=finish<=self.full,'source mask range')
        k=start.bit_count();q=self.q(start)
        need(finish.bit_count()==k and self.q(finish)==q,'same full fibre')
        need(self.compatible(start) and self.compatible(finish),'moving shared-owner endpoint relation')
        if start==finish:return [start]
        if self.n-k<self.lo:return [self.full^m for m in self.path(self.full^start,self.full^finish)]
        if k<self.lo:
            mid=self.low[k][q][1];ans=[start,mid,finish]
        else:
            a=self.normalize(start,k,q);b=self.normalize(finish,k,q)
            out=a[-1]&self.em;goal=b[-1]&self.em;r=(a[-1]&self.rm).bit_count()
            drop=[i for i in self.E if out>>i&1 and not goal>>i&1]
            take=[i for i in self.E if goal>>i&1 and not out>>i&1]
            need(len(drop)==len(take),'same normalized outside margin')
            need(not drop or r==self.lo,'only lowest-band-grade horizontal moves')
            for i,j in zip(drop,take):
                out^=(1<<i)|(1<<j)
                a.append(b[-1] if out==goal else self.replace(out,r,q))
            ans=a+b[::-1]
        ans=[m for i,m in enumerate(ans) if not i or m!=ans[i-1]]
        for m in ans:need(m.bit_count()==k and self.q(m)==q and self.compatible(m),'every arrival signature and pin')
        for a,b in zip(ans,ans[1:]):
            need((a^b).bit_count()<=20 and self.work(a,b)<=(63 if self.moving else 45),'uniform actual event budgets')
        need(len(ans)-1<=52,'uniform event count')
        return ans
    def certificate(self,a,b):
        P=self.path(a,b)
        return {'sources':self.ids,'masks':P,'grade':a.bit_count(),'charge':self.q(a),
                'events':[{'changed_sources':[s for i,s in enumerate(self.ids) if (u^v)>>i&1],
                           'support':(u^v).bit_count(),'owner_clocks':self.work(u,v)} for u,v in zip(P,P[1:])],
                'V_sources':[[s for i,s in enumerate(self.ids) if m>>i&1] for m in P]}
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--start',type=int,required=True);a.add_argument('--finish',type=int,required=True)
    a.add_argument('--2393-at-570',dest='remove',action='store_true');a.add_argument('--moving-pin',action='store_true');z=a.parse_args()
    print(json.dumps(Router(z.remove,z.moving_pin).certificate(z.start,z.finish),indent=2))
