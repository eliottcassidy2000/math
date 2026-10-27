/* canon.c: read lines "n code" and print nauty canonical string (digraph) */
#include "at_nauty.h"
int main(void){ char line[4096]; int n; char code[4000];
  while (fgets(line,sizeof line,stdin)) { if (sscanf(line,"%d %3999s",&n,code)!=2) continue; tour T; if(!t_parse(code,n,&T)){printf("BAD\n");continue;}
    char cs[400]; nt_canon_str(&T,cs); printf("%d %s %s\n", n, cs, code);} return 0; }
