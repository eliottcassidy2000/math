# structural check of the six order-9 coincidences of (score-sequence pair, arc-usage multiset):
# each is TT2 => S => TT1 against TT2 => S^op => TT1 with S strong of order 6, converse but not isomorphic, self-reversed scores
# needs collision_moves.py from the kps_s29 packet (scripts/) in the working directory, plus cls_9.txt
src=open('collision_moves.py').read().split("for n in (7,8,9):")[0]
ns={}; exec(src,ns)
signature, iso, conv = ns['signature'], ns['iso'], ns['conv']
pairs=[(67436800,201393152),(69501184,134349824),(75596032,1140982784),(134414592,203522304),(134482944,1277266176),(1077936640,1140851712)]
for a,b in pairs:
    (ca,ra),(cb,rb)=signature(a,9),signature(b,9)
    S,T=ca[0],cb[0]; m=len(S)
    sc=sorted(bin(x).count('1') for x in S); rev=sorted(m-1-x for x in sc)
    ok=(ra==rb==[2,1] and m==6 and not iso(S,T) and iso(conv(S),T) and sc==rev)
    print(a,b,'runs',ra,'core scores',sc,'structure as claimed:',ok)
