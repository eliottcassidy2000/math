#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
using namespace std; using I=long long;
bool sq(I n){I r=sqrt((long double)n); while(r*r<n)r++;while(r*r>n)r--;return r*r==n;}
int main(int argc,char**argv){int B=argc>1?atoi(argv[1]):100; I count=0, distinct=0, one=0;array<array<int,3>,6> p={{{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}}};for(I a=1;a<=B;a++)for(I b=0;b<=a;b++)for(I c=0;c<=b;c++)for(I d=0;d<=c;d++) {if(gcd(gcd(a,b),gcd(c,d))!=1)continue;count++;I m[3][3]={{a*a+b*b-c*c-d*d,2*(b*c-a*d),2*(b*d+a*c)},{2*(b*c+a*d),a*a-b*b+c*c-d*d,2*(c*d-a*b)},{2*(b*d-a*c),2*(c*d+a*b),a*a-b*b-c*c+d*d}};array<I,9> v;for(int i=0;i<9;i++)v[i]=abs(m[i/3][i%3]);sort(v.begin(),v.end());if(v[0]==0||adjacent_find(v.begin(),v.end())!=v.end())continue;distinct++;bool good[6]={};for(int k=0;k<6;k++){I s=0;for(int i=0;i<3;i++)s+=m[i][p[k][i]]*m[i][p[k][i]];good[k]=sq(s);one+=good[k];}for(int k=0;k<6;k++)if(good[k])for(int l=k+1;l<6;l++)if(good[l]){int eq=0;for(int i=0;i<3;i++)eq+=p[k][i]==p[l][i];if(eq==1){cout<<"FOUND "<<a<<" "<<b<<" "<<c<<" "<<d<<" matchings "<<k<<" "<<l<<"\n";for(auto& row:m){for(I x:row)cout<<x<<" ";cout<<"\n";}cout<<"counts "<<count<<" "<<distinct<<" "<<one<<endl;return 0;}}}cout<<"EXHAUSTED B="<<B<<" primitive_ordered="<<count<<" distinct="<<distinct<<" square_matchings="<<one<<endl;}
