# G_5 = <x+1,4x> = <x+1,3x> on F_11 equals the restriction to F_11 of Stab(inf) in PSL(2,11) on P^1(F_11)
import itertools
p=11
def gen_group(gens):
    idt=tuple(range(p)); G={idt}; fr=[idt]
    while fr:
        new=[]
        for s in fr:
            for g in gens:
                t=tuple(g[s[i]] for i in range(p))
                if t not in G: G.add(t); new.append(t)
        fr=new
    return G
S=tuple((i+1)%p for i in range(p))
G5=gen_group([S,tuple(4*i%p for i in range(p))])
INF=p
def mob(M,x):
    a,b,c,d=M
    if x==INF: return INF if c==0 else a*pow(c,-1,p)%p
    den=(c*x+d)%p
    return INF if den==0 else (a*x+b)*pow(den,-1,p)%p
els={tuple(mob(M,x) for x in range(p+1)) for M in itertools.product(range(p),repeat=4) if (M[0]*M[3]-M[1]*M[2])%p==1}
stab={g[:p] for g in els if g[INF]==INF}
print('|PSL(2,11)| =',len(els),' |Stab(inf)| =',len(stab),' Stab(inf)|F_11 == G_5:',stab==G5)
