"""Ten terminating rewrite schemas for the proposed spectrum's arithmetic monoid."""
from itertools import combinations_with_replacement,product
from math import factorial,gcd,lcm
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
def need(ok,msg):
    if not ok:raise RuntimeError(msg)

class Bank:
    def __init__(self,primes=(5,11,13)):
        self.primes=tuple(primes)
        self.names=('a','c','d','e')+tuple('p'+str(p) for p in primes)+tuple('b'+str(p) for p in primes)
        self.values=(3,49,343,63)+tuple(primes)+tuple(7*p for p in primes)
        self.index={v:i for i,v in enumerate(self.names)}
        self.rules=[]
        def add(name,left,right):
            L=tuple(left.count(v) for v in self.names);R=tuple(right.count(v) for v in self.names)
            self.rules.append((name,L,R))
        for p,q in combinations_with_replacement(primes,2):
            add('BB',['b'+str(p),'b'+str(q)],['c','p'+str(p),'p'+str(q)])
        for p in primes:
            b='b'+str(p);z='p'+str(p)
            add('BC',[b,'c'],['d',z]);add('BD',[b,'d'],['c','c',z])
            add('BE',[b,'e'],['c','a','a',z]);add('BAA',[b,'a','a'],['e',z])
        for p in primes:
            for q in primes:
                if q<p:add('exchange',['b'+str(p),'p'+str(q)],['b'+str(q),'p'+str(p)])
        add('EE',['e','e'],['c']+['a']*4)
        add('EC',['e','c'],['d','a','a'])
        add('ED',['e','d'],['c','c','a','a'])
        add('DD',['d','d'],['c','c','c'])
    def value(self,z):
        out=1
        for n,a in zip(self.values,z):out*=n**a
        return out
    def rank(self,z):
        return 2*z[3]+z[2]+sum((p+3)*z[4+len(self.primes)+i] for i,p in enumerate(self.primes))
    def exponents(self,z):
        j=len(self.primes);a,c,d,e=z[:4]
        return 2*c+3*d+e+sum(z[4+j:]),a+2*e,tuple(z[4+i]+z[4+j+i] for i in range(j))
    def canonical(self,z):
        a,b,other=self.exponents(z);out=[0]*len(z);out[0]=b
        for i,v in enumerate(other):out[4+i]=v
        if a%2==0:out[1]=a//2
        elif a>=3:out[2]=1;out[1]=(a-3)//2
        elif b>=2:out[3]=1;out[0]-=2
        else:
            i=next(i for i,v in enumerate(other) if v);out[4+i]-=1;out[4+len(self.primes)+i]=1
        return tuple(out)
    def moves(self,z):
        for name,L,R in self.rules:
            if all(a>=b for a,b in zip(z,L)):
                yield name,tuple(a-b+c for a,b,c in zip(z,L,R))
    def normalize(self,z):
        first=z;trace=[]
        while True:
            choices=list(self.moves(z))
            if not choices:break
            name,new=choices[0]
            need(self.rank(new)<self.rank(z),'strict natural-number rank')
            need(self.value(new)==self.value(z),'represented value')
            trace.append(name);z=new
        need(z==self.canonical(first),'closed canonical formula')
        need(len(trace)<=self.rank(first),'constructive bound')
        return z,trace

FRACTIONS=((48,49),(20,21),(36,35),(27,25))
def valuations(n):
    out=[]
    for p in (2,3,5,7):
        k=0
        while n%p==0:n//=p;k+=1
        out.append(k)
    return tuple(out),n
def encoded_value(z):
    a,c,d,e=z;return 3**a*49**c*343**d*63**e
def machine(n):
    trace=[];start=n
    while True:
        selected=next((i for i,(a,b) in enumerate(FRACTIONS) if n%b==0),None)
        if selected is None:return n,tuple(trace)
        a,b=FRACTIONS[selected];old,_=valuations(n);new=n//b*a;v,_=valuations(new)
        need(2*v[3]+v[2]<2*old[3]+old[2],'FRACTRAN strict rank')
        need(encoded_value(old)==encoded_value(v),'FRACTRAN semantics')
        n=new;trace.append(selected+1)
def expected_trace(z):
    a,c,d,e=z;t=e//2;out=[1]*t;c+=t
    if e%2:
        if c:
            out.append(2);d+=1
        elif d:
            out.append(3);d-=1
        else:return tuple(out)
    return tuple(out+[4]*(d//2))

def main():
    bank=Bank();checked=0;all_edges=0
    for length in range(7):
        for word in combinations_with_replacement(range(len(bank.names)),length):
            z=tuple(word.count(i) for i in range(len(bank.names)))
            norm,trace=bank.normalize(z);checked+=1
            for name,new in bank.moves(z):
                need(bank.canonical(new)==norm and bank.rank(new)<bank.rank(z),'all competing reductions')
                all_edges+=1
    machine_checks=0;trace_set=set()
    for c,d,e in product(range(9),repeat=3):
        n=3**c*5**d*7**e;end,trace=machine(n)
        need(trace==expected_trace((0,c,d,e)),'complete run formula')
        need(len(trace)<=2*e+d,'all-input native clock bound')
        trace_set.add(trace);machine_checks+=1
    for n in range(1,10001):
        z,inert=valuations(n);end,trace=machine(n)
        need(trace==expected_trace(z) and valuations(end)[1]==inert,'integer execution and inert primes')
    lengths=[];counts=[]
    for k in range(1,31):
        records=[]
        for c in range(k+1):
            for d in range((2*k)//3+1):
                e=2*k-2*c-3*d
                if e>=0:records.append((c,d,e,4*c+6*d))
        actual=sorted({sum(v) for v in records})
        expected=[v for v in range(2*k,5*k+1) if v-2*k not in (1,2,5)]
        need(actual==expected,'complete extremal-power length set')
        lengths.append({'k':k,'lengths':actual,'factorization_count':len(records)})
        if k<=3:
            hist={}
            for c,d,e,a in records:
                L=c+d+e+a;ways=factorial(L)//(factorial(c)*factorial(d)*factorial(e)*factorial(a))
                hist[L]=hist.get(L,0)+ways
            counts.append({'k':k,'ordered_source_counts_by_length':hist,'total':sum(hist.values())})
    # Rational elasticities: six actual atoms suffice; finite controls use divisor exponent enumeration.
    elastic=[]
    for denominator in range(1,13):
        for numerator in range(denominator,(5*denominator)//2+1):
            if gcd(numerator,denominator)>1:continue
            m=numerator-denominator;k=5*denominator-2*numerator
            L=[]
            for b5 in range(min(k,2*m)+1):
                for e in range(2*m-b5+1):
                    for d in range((2*m-b5-e)//3+1):
                        rem=2*m-b5-e-3*d
                        if rem%2:continue
                        c=rem//2;a=4*m-2*e;p5=k-b5
                        if a>=0:L.append(a+c+d+e+p5+b5)
            need(min(L)==3*denominator and max(L)==3*numerator,'all rational elasticities construction')
            elastic.append([numerator,denominator,m,k])
    data={'status':'PROVED_SCOPED','alphabet':bank.names,'values':bank.values,
          'rewrite_instances':len(bank.rules),'schemas':sorted({x[0] for x in bank.rules}),
          'bounded_factorizations':checked,'competing_edges':all_edges,
          'core_exponent_tests':machine_checks,'positive_integer_tests':10000,
          'length_sets':lengths,'ordered_source_counts':counts,'rational_elasticity_controls':elastic}
    (HERE/'normal_forms.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS:',checked,'normal forms;',all_edges,'competing edges; all-input formula controls')
if __name__=='__main__':main()
