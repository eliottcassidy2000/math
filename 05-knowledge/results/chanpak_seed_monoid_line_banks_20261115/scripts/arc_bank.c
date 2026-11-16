/* Arc banks.  A complete k-arc bank is a tournament T and k arcs whose 2^k joint orientations (all other arcs fixed)
   give 2^k consecutive odd Hamiltonian-path counts.  (The contacts of one head are the special case of k arcs at one
   vertex, where THM-5272 forces every pair to interact; arbitrary arcs can be independent: in the transitive
   4-tournament the two span-2 arcs give counts 1, 3, 3, 5 with second difference 0.)
   usage: arc_bank n k < masks   (class representatives of order n) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static int n; static int out[10];
static ll dp[1<<9][9];
static ll ham(void){ int N=1<<n; for(int S=0;S<N;S++) memset(dp[S],0,sizeof(ll)*n); for(int v=0;v<n;v++) dp[1<<v][v]=1;
    for(int S=1;S<N;S++) for(int v=0;v<n;v++){ ll c=dp[S][v]; if(c){ int t=out[v]&~S; while(t){ int b=t&-t; dp[S|b][__builtin_ctz(b)]+=c; t^=b; } } }
    ll s=0; for(int v=0;v<n;v++) s+=dp[N-1][v]; return s; }
static void flip(int a,int b){ if((out[a]>>b)&1){ out[a]&=~(1<<b); out[b]|=1<<a; } else { out[b]&=~(1<<a); out[a]|=1<<b; } }
int main(int argc,char**argv){
    n=atoi(argv[1]); int k=atoi(argv[2]); int pa[40],pb[40],np=0; for(int a=0;a<n;a++) for(int b=a+1;b<n;b++){ pa[np]=a; pb[np]=b; np++; }
    char line[256]; long nb=0, found=0, cubes=0; ll bestspread=1LL<<60; ll bestvals[64]; int bestbase[10]; int bestarcs[6][2]; long distinct_cubes=0; long sp[64]; memset(sp,0,sizeof(sp)); long bal4=0, bal8=0, bal16=0, sum16[16]; memset(sum16,0,sizeof(sum16)); int shown=0; ll minstart=-1,maxstart=-1; long affine=0; long withhead=0;
    int idx[6];
    while(fgets(line,sizeof(line),stdin)){ char*p=line; int ok=1; for(int i=0;i<n;i++){ char*e; long v=strtol(p,&e,10); if(e==p){ok=0;break;} out[i]=(int)v; p=e; } if(!ok) continue; nb++;
        int base[10]; memcpy(base,out,sizeof(base));
        /* iterate k-subsets of arcs */
        for(int i=0;i<k;i++) idx[i]=i;
        while(1){
            ll val[64]; int K=1<<k; ll mn=0,mx=0;
            for(int o=0;o<K;o++){ memcpy(out,base,sizeof(base)); for(int i=0;i<k;i++) if((o>>i)&1) flip(pa[idx[i]],pb[idx[i]]); val[o]=ham(); if(o==0||val[o]<mn) mn=val[o]; if(o==0||val[o]>mx) mx=val[o]; }
            cubes++;
            if(k==3){ int c4[4]={0,0,0,0}, c8[8]={0}, c16[16]={0}; ll tot=0; for(int o=0;o<8;o++){ c4[val[o]&3]++; c8[val[o]&7]++; c16[val[o]&15]++; tot+=val[o]; }
              sum16[tot&15]++; if(c4[1]==4&&c4[3]==4){ bal4++; if(c8[1]==2&&c8[3]==2&&c8[5]==2&&c8[7]==2){ bal8++; int ok16=1; for(int r=1;r<16;r+=2) if(c16[r]!=1) ok16=0; bal16+=ok16; } } }
            { ll sv[64]; memcpy(sv,val,sizeof(ll)*K); for(int i=0;i<K;i++) for(int j=i+1;j<K;j++) if(sv[j]<sv[i]){ ll t=sv[i]; sv[i]=sv[j]; sv[j]=t; }
              int dd=1; for(int i=1;i<K;i++) if(sv[i]==sv[i-1]) dd=0;
              if(dd){ distinct_cubes++; ll spr=(mx-mn)/2; if(spr<64) sp[spr]++; if(spr<bestspread){ bestspread=spr; memcpy(bestvals,sv,sizeof(ll)*K); memcpy(bestbase,base,sizeof(base)); for(int i=0;i<k;i++){ bestarcs[i][0]=pa[idx[i]]; bestarcs[i][1]=pb[idx[i]]; } } } }
            if(mx-mn==2*(K-1)){ unsigned long long seen=0; int dist=1; for(int o=0;o<K;o++){ int id=(int)((val[o]-mn)/2); if((seen>>id)&1){ dist=0; break; } seen|=1ULL<<id; }
                if(dist){ found++; if(minstart<0||mn<minstart) minstart=mn; if(mn>maxstart) maxstart=mn;
                    /* is the cube affine (all mixed second differences zero)? */
                    int aff=1; for(int o=0;o<K&&aff;o++) for(int i=0;i<k&&aff;i++) for(int j=i+1;j<k;j++) if(!((o>>i)&1)&&!((o>>j)&1)){ if(val[o|(1<<i)|(1<<j)]-val[o|(1<<i)]-val[o|(1<<j)]+val[o]!=0){ aff=0; break; } }
                    affine+=aff;
                    /* do the k arcs share a vertex? */
                    int common=0; for(int v=0;v<n;v++){ int all=1; for(int i=0;i<k;i++) if(pa[idx[i]]!=v&&pb[idx[i]]!=v) all=0; common|=all; } withhead+=common;
                    if(shown<6){ shown++; printf("BANK counts %lld..%lld%s%s ; base:",mn,mx,aff?" (affine)":"",common?" (arcs at one vertex)":""); for(int i=0;i<n;i++) printf(" %d",base[i]); printf(" ; arcs:"); for(int i=0;i<k;i++) printf(" (%d,%d)",pa[idx[i]],pb[idx[i]]); printf(" ; values:"); for(int o=0;o<K;o++) printf(" %lld",val[o]); printf("\n"); } } }
            int i=k-1; while(i>=0&&idx[i]==np-k+i) i--; if(i<0) break; idx[i]++; for(int j=i+1;j<k;j++) idx[j]=idx[j-1]+1;
        }
    }
    printf("# cubes with 2^k distinct counts: %ld ; least spread (max-min)/2 = %lld (a complete bank needs %d) ; values:",distinct_cubes,bestspread,(1<<k)-1); for(int o=0;o<(1<<k);o++) printf(" %lld",bestvals[o]); printf(" ; base:"); for(int i=0;i<n;i++) printf(" %d",bestbase[i]); printf(" ; arcs:"); for(int i=0;i<k;i++) printf(" (%d,%d)",bestarcs[i][0],bestarcs[i][1]); printf("\n");
    printf("# spread histogram (spread: cubes):"); for(int i=0;i<24;i++) if(sp[i]) printf(" %d:%ld",i,sp[i]); printf("\n");
    if(k==3){ printf("# residues on 3-cubes: H mod 4 balanced (4+4): %ld ; and H mod 8 each odd residue twice: %ld ; and H mod 16 a permutation of the odd residues: %ld\n",bal4,bal8,bal16); printf("# sum of the 8 counts mod 16:"); for(int r=0;r<16;r++) if(sum16[r]) printf(" %d:%ld",r,sum16[r]); printf("\n"); }
    printf("# order %d, k = %d: %ld classes, %ld arc cubes examined; complete banks: %ld (affine: %ld; all arcs at one vertex: %ld); first counts from %lld to %lld\n",n,k,nb,cubes,found,affine,withhead,minstart,maxstart);
    return 0;
}
