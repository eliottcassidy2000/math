"""Root-deletion cycle-index law and fixed-spine fibre weights, finite controls."""
from research import *


def canonical(T):
    return min(induced(T,p) for p in permutations(range(len(T))))


def automorphisms(T):
    return [p for p in permutations(range(len(T))) if induced(T,p)==T]


def image_subset(mask,p):
    return sum(1<<p[i] for i in range(len(p)) if mask>>i&1)


def cycles(p):
    seen=set();answer=[]
    for i in range(len(p)):
        if i not in seen:
            cycle=[];j=i
            while j not in seen:
                seen.add(j);cycle.append(j);j=p[j]
            answer.append(cycle)
    return answer


def inventory(T,G):
    starts=tables(opposite(T))[1][-1]
    orbit_reps={min(image_subset(mask,p) for p in G) for mask in range(1<<len(T))}
    direct=Counter(sum(starts[i] for i in range(len(T)) if mask>>i&1) for mask in orbit_reps)
    burnside=Counter()
    for p in G:
        poly=Counter({0:1})
        for cyc in cycles(p):
            weight=sum(starts[i] for i in cyc);new=poly.copy()
            for k,v in poly.items():new[k+weight]+=v
            poly=new
        burnside.update(poly)
    require(all(v%len(G)==0 for v in burnside.values()),'Burnside divisibility')
    require(direct==Counter({k:v//len(G) for k,v in burnside.items()}),'subset orbit vs Burnside')
    return direct


def main():
    cores=[];direct_rooted=Counter();direct_strong=Counter();fibres=[];ordinary=[1]
    for n in range(1,6):
        classes=sorted({canonical(tiling(n,b)) for b in range(1<<((n-1)*(n-2)//2))})
        ordinary.append(len(classes));mass=0
        for T in classes:
            G=automorphisms(T);H=tables(T)[0][-1];starts=tables(opposite(T))[1][-1]
            require(H%len(G)==0,'free Hamilton action')
            mass+=H//len(G)
            roots=sorted({min(p[r] for p in G) for r in range(n)})
            for r in roots:
                stabilizer=sum(p[r]==r for p in G)
                require(starts[r]%stabilizer==0,'rooted free path action')
                direct_rooted[n,starts[r]]+=1
                if n>=3 and strong(T):direct_strong[n,starts[r]]+=1
            if n==1 or strong(T):
                B=inventory(T,G)
                cores.append(dict(n=n,H=H,B=B,subset_orbits=sum(B.values()),automorphisms=len(G)))
            tables.cache_clear()
        require(mass==1<<((n-1)*(n-2)//2),'fixed-spine fibre sum')
        fibres.append(dict(n=n,ordinary_classes=len(classes),fixed_spine_charts=mass))
    bound=6;Bsize=[0]*bound
    for c in cores:Bsize[c['n']]+=c['subset_orbits']
    seq=[1]+[0]*(bound-1)
    for n in range(1,bound):seq[n]=sum(Bsize[j]*seq[n-j] for j in range(1,n+1))
    rooted=[0]+seq
    for n in range(1,6):
        require(rooted[n]==sum(v for (m,k),v in direct_rooted.items() if m==n),'all rooted class law')
    # x[Z+R(1-L)^(-1)I], represented by explicit ordered SCC prefixes.
    prediction=Counter();prefix=Counter()
    for c in cores:
        for s,count in c['B'].items():
            if s:
                prefix[c['n'],s]+=count
                if s!=c['H']:prediction[c['n']+1,s]+=count
    for _ in range(bound):
        next_prefix=Counter()
        for (size,s),weight in prefix.items():
            for c in cores:
                n=size+c['n']
                if n>=bound:continue
                key=n,s*c['H']
                next_prefix[key]+=weight*c['subset_orbits']
                prediction[n+1,key[1]]+=weight*(c['subset_orbits']-1)
        prefix=next_prefix
    for n in range(1,6):
        require({k:v for (m,k),v in prediction.items() if m==n}=={k:v for (m,k),v in direct_strong.items() if m==n},('strong start-budget law',n))
    require([prediction[n,1] for n in range(3,7)]==[1,3,8,21],'one-start projection')
    require([prediction[n,2] for n in range(4,7)]==[1,5,18],'two-start projection')
    print(json.dumps(dict(status='PASS',ordinary_classes=ordinary,rooted_classes_through_six=rooted,
                         fixed_spine_fibres=fibres,core_subset_inventories=cores,
                         strong_start_budgets=[dict(n=n,start_paths=s,classes=c) for (n,s),c in sorted(prediction.items()) if c]),indent=2,sort_keys=True))


if __name__=='__main__':main()
