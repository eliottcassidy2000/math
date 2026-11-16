from verify import profile
for B,A,Q,a in [
    (3,1444,(1343,311,774,981),(266,68,25,306,291,135,317,351,38)),
    (5,722,(663,177,393,506),(266,276,33,22,69,283,431,319,28))]:
    h,g,r,s=profile(a,8)
    assert (h,g,r,s)==Q
    assert (h-2*s+2*g,r+s-g)==(B,A)
print('Both second-rail rooted profiles verified.')
