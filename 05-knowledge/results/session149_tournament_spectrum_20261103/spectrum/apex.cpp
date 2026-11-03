#include <bits/stdc++.h>
using namespace std;
int main(){const int B=2000000; vector<int> best(B+1,0), maskv(B+1,0); best[1]=1; for(int m=2;m<=22;m++){long maxmask=1L<<m; long num=0; for(long mask=0;mask<maxmask;mask++){if(!(mask&1) || (mask>>(m-1)&1)) continue; // bit means x beats vertex
 long h=1; // x first, since x beats v1; x last since vm beats x; account two endpoint paths, not one
 h=2; for(int i=0;i<m;i++) if(!(mask>>i&1)){if(i+1<m && (mask>>(i+1)&1)) h++; for(int j=0;j<i;j++) if(mask>>j&1) h+=1L<<(i-j-1); } if(h<=B && !best[h]){best[h]=m+1;maskv[h]=mask;num++;}}
 int gap=3;while(gap<=B && (gap==7||gap==21||best[gap]))gap+=2;cout<<"m="<<m<<" added="<<num<<" first_gap="<<gap<<endl; }
 ofstream f("apex-record.csv");f<<"h,n,mask\n";for(int h=1;h<=B;h+=2)if(best[h])f<<h<<","<<best[h]<<","<<maskv[h]<<"\n";
}
