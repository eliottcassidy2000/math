#!/usr/bin/env python3
from itertools import product
import json
from math import gcd


def need(x,msg):
    if not x:
        raise RuntimeError(msg)


def cylinder(word):
    B=S=0
    for k in word:
        B=3*B+2**S
        S+=k
    h=len(word)
    c=((2**S-B)*pow(3**h,-1,2**(S+1)))%2**(S+1)
    d=(3**h*c+B)//2**S
    return c,d,S


def odd_run(n,h):
    vals=[]; states=[n]; literal=0
    for _ in range(h):
        n=3*n+1; literal+=1; k=0
        while n%2==0:
            n//=2; k+=1; literal+=1
        vals.append(k); states.append(n)
    return tuple(vals),states,literal


def shortcut(n,steps):
    for _ in range(steps):
        n=(3*n+1)//2 if n%2 else n//2
    return n


def strict(states):
    return states[0]>1 and states[-1]==1 and 1 not in states[:-1]


def main():
    words=lifts=prefixes=strict_words=modulus_hostiles=0
    for h in range(1,4):
        for word in product(range(1,5),repeat=h):
            c,d,S=cylinder(word)
            ks,canonical,clock=odd_run(c,h)
            need(ks==word and canonical[-1]==d,'canonical cylinder')
            words+=1
            strict_words+=int(strict(canonical))
            if strict(canonical):
                for q in range(1,33):
                    t=q//gcd(3**h,q)
                    alias=c+2**(S+1)*t
                    ka,sa,clocka=odd_run(alias,h)
                    need(ka==word and sa[-1]>1 and 1 not in sa,'same-word nonterminal alias')
                    need(((sa[-1]-1)//2)%q==0,'identical endpoint residue to ROOT')
                    need((alias-c)//2==2**S*q//gcd(3**h,q),'identical source precision period')
                    need(clocka==h+S and shortcut(alias,S)==sa[-1],'alias clocks unchanged')
                    modulus_hostiles+=1
            for t in range(8):
                n=c+2**(S+1)*t
                ks,states,literal=odd_run(n,h)
                need(ks==word,'actual lift branch')
                need(states[-1]==d+2*3**h*t,'actual lift endpoint')
                need(literal==h+S,'literal clock')
                need(shortcut(n,S)==states[-1],'shortcut clock')
                need((states[-1]==1)==(t==0 and d==1),'literal ROOT singleton')
                need(strict(states)==(t==0 and strict(canonical)),'first-hit singleton')
                Sj=0
                for j,k in enumerate((0,)+word):
                    if j:
                        Sj+=k
                    need(states[j]==canonical[j]+3**j*2**(S-Sj+1)*t,'prefix affine identity')
                    if t:
                        need(states[j]>1,'higher lift cannot hit one')
                    prefixes+=1
                lifts+=1
    examples=[]
    for word,a,b,q in [((4,),5,37,3),((6,),21,405,9),((4,2),5,133,9)]:
        ka,sa,ca=odd_run(a,len(word)); kb,sb,cb=odd_run(b,len(word))
        need(ka==kb==word,'hostile shared valuation word')
        need((sa[-1]-1)//2%q==(sb[-1]-1)//2%q,'hostile equal endpoint residue')
        need(1 in sa and 1 not in sb,'hostile distinct hit guards')
        examples.append({'word':word,'sources':[a,b],'states':[sa,sb],
                         'endpoint_index_modulus':q,'literal_clock':ca,
                         'shortcut_clock':sum(word),'strict_first_hit':[strict(sa),strict(sb)]})
    need(21%6==405%6==3,'native arithmetic source class')
    need(((21-1)//2)%192==((405-1)//2)%192==10,'source-index precision hostile')
    print(json.dumps({'status':'PASS','word_count':words,'source_lifts':lifts,
                      'prefix_checks':prefixes,'strict_first_hit_words':strict_words,'finite_modulus_same_word_hostiles':modulus_hostiles,
                      'hostiles':examples},indent=2,sort_keys=True))


if __name__=='__main__':
    main()
