"""Bounded exact hostile probe for fixed-(time, odd-count) injectivity.
Reproduce: python3 04-computation/experiments/crossroads_20260926_flow_collision.py
"""
from fractions import Fraction
import json


def require(condition,message):
    if not condition:
        raise AssertionError(message)


def trajectory(n,k):
    out=[]
    e=0
    source=n
    for j in range(1,k+1):
        odd=n&1
        e+=odd
        n=(3*n+1)//2 if odd else n//2
        w=Fraction(3**e,2**j)
        require(w>1,"prefix slope is not strictly above one")
        out.append({"odd":odd,"value":n,"slope":str(w),"h":str(Fraction(n,w)-source)})
    return out


def main():
    limit=1_000_000
    max_depth=40
    active=[(n,n,0) for n in range(1,limit+1)]
    total=0
    powers3=[3**e for e in range(max_depth+1)]
    for k in range(1,max_depth+1):
        seen={}
        next_active=[]
        collision=None
        for source,v,e in active:
            odd=v&1
            e+=odd
            v=(3*v+1)//2 if odd else v//2
            if powers3[e]<=2**k:
                continue
            next_active.append((source,v,e))
            key=(e,v)
            if key in seen and collision is None:
                collision=(seen[key],source,v,e)
            else:
                seen[key]=source
        active=next_active
        total+=len(active)
        print(json.dumps({"k":k,"survivors":len(active),"distinct_hubs_and_counts":len(seen)}))
        if collision:
            n,m,v,e=collision
            print(json.dumps({"status":"REFUTED","source_limit":limit,"smallest_time_in_universe":k,
                              "sources":[n,m],"hub":v,"e":e,
                              "word_n":trajectory(n,k),"word_m":trajectory(m,k),
                              "tested_source_prefixes":total},sort_keys=True))
            return
    print(json.dumps({"status":"FINITE-EXACT no collision in the source box",
                      "source_limit":limit,"max_depth":max_depth,"tested_source_prefixes":total},sort_keys=True))
    symbolic_probe()
    offset_width_audit()
    targeted_coordinated_carry()


def symbolic_probe():
    """All admissible parity words, hence all integer heights, through depth24."""
    active=[(0,0,"")]
    checked=0
    for k in range(1,25):
        next_active=[]
        seen={}
        collision=None
        for e,C,word in active:
            for odd in [0,1]:
                f=e+odd
                if 3**f <= 2**k: continue
                D=3*C+2**(k-1) if odd else C
                newword=word+str(odd)
                next_active.append((f,D,newword))
                key=(f,D%(3**f))
                if key in seen and collision is None:
                    collision=(seen[key],(D,newword),f)
                else:
                    seen[key]=(D,newword)
        active=next_active
        checked+=len(active)
        if collision:
            (C,u),(D,w),e=collision
            modulus=3**e
            residue=C*pow(2**k,-1,modulus)%modulus
            v=residue
            while 2**k*v<=max(C,D): v+=modulus
            n=(2**k*v-C)//modulus
            m=(2**k*v-D)//modulus
            require(n!=m,"sources coincide")
            require(trajectory(n,k)[-1]["value"]==v,"first source misses hub")
            require(trajectory(m,k)[-1]["value"]==v,"second source misses hub")
            print(json.dumps({"status":"REFUTED all-height symbolic collision","k":k,"e":e,
                              "sources":[n,m],"hub":v,"words":[u,w],"carries":[C,D],
                              "offsets":[str(Fraction(C,modulus)),str(Fraction(D,modulus))],
                              "symbolic_words_checked":checked},sort_keys=True))
            return
    print(json.dumps({"status":"FINITE-EXACT all-height injectivity through fixed depth only",
                      "max_depth":24,"symbolic_words_checked":checked},sort_keys=True))


def targeted_coordinated_carry():
    """Decode C=(d+1)3^e-2^e, differing by d3^e from a dense word.
    This is a bounded structured search, not a minimal-counterexample claim.
    """
    def decode(C,e):
        positions=[]
        for i in range(1,e+1):
            if C<=0:return None
            p=(C&-C).bit_length()-1
            if positions and p<=positions[-1]:return None
            if i==1 and p!=0:return None
            if i>1 and 2**p>=3**(i-1):return None
            positions.append(p)
            C-=3**(e-i)*2**p
        return positions if C==0 else None
    for e in range(3,201):
        for d in range(4,e//3+1,4):
            A=3**e
            C=(d+1)*A-2**e
            positions=decode(C,e)
            if positions is None:continue
            k=max(e,positions[-1]+1)
            if A<=2**k:continue
            lowcarry=A-2**e
            hub=lowcarry*pow(2**k,-1,A)%A
            while 2**k*hub<=C:hub+=A
            large=(2**k*hub-lowcarry)//A
            small=(2**k*hub-C)//A
            require(large-small==d,"wrong source difference")
            large_path=trajectory(large,k)
            small_path=trajectory(small,k)
            require(large_path[-1]["value"]==small_path[-1]["value"]==hub,"sources do not merge")
            word_large="".join(str(row["odd"]) for row in large_path)
            word_small="".join(str(row["odd"]) for row in small_path)
            require(word_large=="1"*e+"0"*(k-e),"wrong dense parity word")
            require([i for i,ch in enumerate(word_small) if ch=="1"]==positions,"decoded word disagrees with actual orbit")
            require(sum(row["odd"] for row in small_path)==e,"wrong odd count")
            parameter=pow(A,-1,2**(k-e))
            require(large==2**e*parameter-1,"modular inverse source construction failed")
            print(json.dumps({"status":"REFUTED exact positive coordinated collision",
                "search_family":"C=(d+1)3^e-2^e; 3<=e<=200; d=4,8,...<=e/3",
                "k":k,"e":e,"source_difference":d,"sources":[small,large],
                "source_parameter_u":parameter,"source_construction":"u=3^(-e) mod2^(k-e); N=2^e*u-1; sources N-d,N",
                "hub":hub,"words":[word_small,word_large],"odd_positions_small":positions,
                "normalized_offset_difference":str(Fraction(C-lowcarry,A))},sort_keys=True))
            return
    raise AssertionError("Expected structured hostile was not recovered")


def offset_width_audit():
    maximum=(Fraction(0),0,0)
    first_opening=None
    cases=0
    for k in range(1,33):
        for e in range(1,k+1):
            if 3**e<=2**k:continue
            positions=[0]
            for i in range(2,e+1):
                p=0
                while 2**(p+1)<3**(i-1):p+=1
                positions.append(min(p,k-e+i-1))
            hi=sum((Fraction(2**p,3**i) for i,p in enumerate(positions,1)),Fraction(0))
            lo=1-Fraction(2,3)**e
            width=hi-lo
            if k<=31:
                require(width<4,"width certificate failed")
                maximum=max(maximum,(width,k,e))
                cases+=1
            elif width>=4 and first_opening is None:
                first_opening={"k":k,"e":e,"width":str(width)}
    require(maximum==(Fraction(13805179460,3486784401),31,20),"wrong maximum width")
    print(json.dumps({"status":"PROVED finite all-height injectivity k<=31 by carry interval width",
                      "cases":cases,"max_width":str(maximum[0]),"attaining_k":maximum[1],
                      "attaining_e":maximum[2],"first_width_opening":first_opening},sort_keys=True))


if __name__=="__main__":
    main()
