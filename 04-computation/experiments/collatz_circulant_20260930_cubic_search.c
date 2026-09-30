/* integer solutions of x^3+2x+y^3+2y+z^3+2z = xyz+1 with |x|,|y| <= B (z solved): z^3 + (2-xy) z + (f(x)+f(y)-1) = 0 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef __int128 i128;
static i128 f(i128 t){ return t*t*t + 2*t; }
int main(int argc,char**argv){
  long B = argc>1 ? atol(argv[1]) : 5000; long cnt=0;
  for(long x=-B;x<=B;x++) for(long y=x;y<=B;y++){
    i128 p = 2 - (i128)x*y, q = f(x)+f(y)-1;
    /* real roots of z^3 + p z + q: use Cardano/trig via double, then check integer neighbours */
    double pd=(double)p, qd=(double)q; double cand[3]; int nc=0;
    double disc = qd*qd/4 + pd*pd*pd/27;
    if(disc>=0){ double s=sqrt(disc); double u=cbrt(-qd/2+s), v=cbrt(-qd/2-s); cand[nc++]=u+v; }
    else { double r=2*sqrt(-pd/3); double th=acos(3*qd/(pd*r))/3.0; /* pd<0 */
           cand[nc++]=r*cos(th); cand[nc++]=r*cos(th-2.0943951023931953); cand[nc++]=r*cos(th+2.0943951023931953); }
    for(int i=0;i<nc;i++){ long z0=(long)llround(cand[i]); for(long z=z0-2;z<=z0+2;z++){ i128 val=(i128)z*z*z + p*z + q; if(val==0){ printf("solution: x=%ld y=%ld z=%ld\n",x,y,z); cnt++; } } }
  }
  printf("done B=%ld, solutions found: %ld\n",B,cnt); return 0; }
