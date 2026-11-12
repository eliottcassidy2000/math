#!/usr/bin/env python3
"""Independent literal module checks for the marked-spine birth grammar."""
import itertools,json,math

def require(x,msg):
    if not x: raise RuntimeError(msg)

def graph(word):
    m=len(word); rows=[0]*(m+1)
    for i in range(m):
        for j in range(i+1,m): rows[i]|=1<<j
        if word[i]=='1': rows[m]|=1<<i
        else: rows[i]|=1<<m
    return rows

def modules(rows):
    n=len(rows); allv=(1<<n)-1; ans=[]
    for s in range(1,allv+1):
        if s.bit_count()<2 or s==allv: continue
        if all((rows[v]&s) in (0,s) for v in range(n) if not(s>>v&1)): ans.append(s)
    return ans

def predicted_modules(word):
    m=len(word); allv=(1<<(m+1))-1; ans=set()
    for a in range(m):
        for b in range(a,m):
            s=((1<<(b-a+1))-1)<<a
            if b>a and len(set(word[a:b+1]))==1: ans.add(s)
            if set(word[:a])<=set('0') and set(word[b+1:])<=set('1'):
                s|=1<<m
                if s!=allv: ans.add(s)
    return sorted(ans)

def module_number(word):
    if not word:return 0
    runs=[len(list(g)) for _,g in itertools.groupby(word)]
    p=len(word)-len(word.lstrip('0'));q=len(word)-len(word.rstrip('1'))
    return sum(r*(r-1)//2 for r in runs)+(p+1)*(q+1)-1-int(p+q==len(word))

def transfer(word):
    o=[1]; e=[0]
    for bit in word:
        z=max(len(o),len(e))+1;o+= [0]*(z-len(o));e += [0]*(z-len(e))
        if bit=='1': e=[x+y for x,y in zip(e,o)]
        else:o=[o[0]]+[o[i]+e[i-1] for i in range(1,z)]
    while len(o)>1 and o[-1]==0:o.pop()
    return o

def prime(rows):return len(rows)>=3 and not modules(rows)
def induced(rows,vs):return [sum((1<<j) for j,v in enumerate(vs) if rows[u]>>v&1) for u in vs]
def brute_prime_inventory(word):
    rows=graph(word);n=len(rows);out=[1]
    for s in range(1<<n):
        r=s.bit_count()
        if r<3:continue
        vs=[v for v in range(n) if s>>v&1]
        if prime(induced(rows,vs)):
            require(r%2==1,'even prime in apex class')
            j=(r-1)//2
            out +=[0]*max(0,j+1-len(out));out[j]+=1
    return out

def automorphisms_by_degrees(rows):
    groups={}
    for i,r in enumerate(rows):groups.setdefault(r.bit_count(),[]).append(i)
    gs=list(groups.values()); count=0
    for ps in itertools.product(*(itertools.permutations(g) for g in gs)):
        p=list(range(len(rows)))
        for g,q in zip(gs,ps):
            for i,j in zip(g,q):p[i]=j
        if all(((rows[i]>>j)&1)==((rows[p[i]]>>p[j])&1) for i in range(len(rows)) for j in range(len(rows))):count+=1
    return count

def main():
    checked=0; inventories=0
    for m in range(10):
        for bits in itertools.product('01',repeat=m):
            word=''.join(bits); actual=modules(graph(word))
            require(actual==predicted_modules(word),('module grammar',word,actual,predicted_modules(word)))
            require(len(actual)==module_number(word),('module total',word))
            require(prime(graph(word))==(m>=2 and m%2==0 and word=='10'*(m//2)),('prime criterion',word))
            checked+=1
            if m<=7:
                require(transfer(word)==brute_prime_inventory(word),('prime inventory',word,transfer(word),brute_prime_inventory(word)))
                inventories+=1
    tower=[]
    prev=[1];cur=[1,1]
    for k in range(1,101):
        word='10'*k; coeff=transfer(word)
        formula=[math.comb(k+j,2*j) for j in range(k+1)]
        require(coeff==formula,('binomial inventory',k))
        if k>=2:
            nxt=[0]*(len(cur)+1)
            for i,x in enumerate(cur):nxt[i]+=2*x;nxt[i+1]+=x
            for i,x in enumerate(prev):nxt[i]-=x
            require(nxt==coeff,('recurrence',k));prev,cur=cur,nxt
        if k<=7:
            rows=graph(word)
            require(prime(rows),('W prime',k))
            aut=automorphisms_by_degrees(rows)
            require(aut==(3 if k==1 else 1),('W aut',k,aut))
            if k>=2:
                q=len(rows); ext=[r for r in rows]+[(1<<q-1)|1]
                # New vertex beats old apex and first spine vertex only.
                for i in range(q):
                    if not(ext[q]>>i&1):ext[i]|=1<<q
                require(prime(ext),('even prime extension',k))
            tower.append({'order':2*k+1,'prime_induced_counts':coeff[1:],'automorphisms':aut})
    print(json.dumps({'status':'PASS','literal_module_words_through_length':9,'literal_module_words':checked,'induced_prime_inventory_words_through_length':7,'induced_prime_inventory_words':inventories,'transfer_binomial_recurrence_through_k':100,'tower':tower},indent=2))
if __name__=='__main__':main()
