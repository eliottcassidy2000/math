/* strongform.c -- PROOF.md Section 12: the closed score-certified form (a minimum cut X and a CLOSED Ham path from
 * w in X to u notin X with s(u) + s(w) != n-1) on nonregular strong self-converse classes. */
#include "at_nauty.h"
typedef struct { const tour *T; uint32_t sep[MAXV]; int found; } sctx;
static int cb(const int *Q, void *vc) { sctx *c = (sctx *)vc; int n = c->T->n; int a = Q[0], b = Q[n-1];
  if ((c->sep[a] >> b & 1u) && t_has(c->T, b, a) && c->T->s[a] + c->T->s[b] != n - 1) { c->found = 1; return 1; } return 0; }
int main(int argc, char **argv) { int n = atoi(argv[1]); nauty_check(WORDSIZE,1,n,NAUTYVERSIONID);
  char line[1024]; tour T; long long sc = 0, fail = 0; static uint32_t cuts[4096]; int nc;
  while (fgets(line, sizeof line, stdin)) { if (line[0] != '0' && line[0] != '1') continue; line[strcspn(line,"\r\n")] = 0;
    t_parse(line, n, &T); if (t_regular(&T)) continue;
    int ss[MAXV]; for (int v=0; v<n; v++) ss[v]=T.s[v]; for (int i=0;i<n;i++) for(int j=i+1;j<n;j++) if (ss[j]<ss[i]) {int x=ss[i];ss[i]=ss[j];ss[j]=x;}
    int sym=1; for (int i=0;i<n;i++) if (ss[i]+ss[n-1-i]!=n-1) sym=0; if (!sym || !nt_selfconverse(&T)) continue; sc++;
    int lam = t_mincuts(&T, cuts, 4096, &nc); sctx c; memset(&c,0,sizeof c); c.T=&T; uint32_t F=fullmask(n);
    for (int i=0;i<nc;i++) for (int a=0;a<n;a++) if (cuts[i]>>a&1u) c.sep[a] |= F & ~cuts[i];
    t_hpaths(&T, cb, &c); if (!c.found) { fail++; printf("STRONGFORM-FAIL %s lambda=%d delta=%d\n", line, lam, t_delta(&T)); } }
  printf("n=%d nonregular strong self-converse=%lld strong-form failures=%lld\n", n, sc, fail); return 0; }
