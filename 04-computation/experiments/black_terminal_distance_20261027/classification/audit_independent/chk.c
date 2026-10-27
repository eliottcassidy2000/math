/* chk.c: read "n code" lines; report strong, lambda, delta, #mincuts, black descent existence (exact), canonical string */
#include "at_nauty.h"
#include "cand.h"
int main(void){ char line[4096]; int n; char code[4000];
  while (fgets(line,sizeof line,stdin)) { if (sscanf(line,"%d %3999s",&n,code)!=2) continue; tour T; if(!t_parse(code,n,&T)){printf("BAD\n");continue;}
    static uint32_t cuts[1<<16]; int nc; int lam=t_mincuts(&T,cuts,1<<16,&nc); int Q[MAXV], l2; long long ncand=0,nblue=0;
    int bd = black_descent(&T,0,Q,&l2,&ncand,&nblue); char cs[400]; nt_canon_str(&T,cs);
    printf("n=%d strong=%d lambda=%d delta=%d ncuts=%d blackDescent=%d (candidates enumerated %lld, blue %lld) canon=%s\n", n, t_strong_mask(&T,fullmask(n)), lam, t_delta(&T), nc, bd, ncand, nblue, cs);} return 0; }
