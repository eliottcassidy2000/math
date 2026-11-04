#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
int count_hp(int code,int a,int c){array<int,8>r{};int bit=0;for(int i=0;i<6;i++)for(int j=i+1;j<6;j++,bit++){if(code>>bit&1)r[i]|=1<<j;else r[j]|=1<<i;}r[6]=1<<7;for(int i=0;i<6;i++){if(a>>i&1)r[i]|=1<<6;else r[6]|=1<<i;if(c>>i&1)r[i]|=1<<7;else r[7]|=1<<i;}int dp[256][8]={};for(int i=0;i<8;i++)dp[1<<i][i]=1;for(int s=1;s<256;s++)for(int i=0;i<8;i++)if(dp[s][i]){int rem=r[i]&~s;while(rem){int b=rem&-rem;rem-=b;dp[s|b][__builtin_ctz((unsigned)b)]+=dp[s][i];}}int result=0;for(int i=0;i<8;i++)result+=dp[255][i];return result;}
int main(){const int N=32768,K=484;vector<unsigned char>table(N*K);ifstream f("mod4_six_vectors.bin",ios::binary);f.read((char*)table.data(),table.size());if(!f||f.peek()!=EOF)throw runtime_error("binary size");vector<int>B;for(int m=0;m<4096;m++)if(__builtin_popcount(m&63)<=2&&__builtin_popcount(m>>6)<=2)B.push_back(m);int queries=0,coefficients=0;
 auto eval=[&](int code,int a,int c){int mask=a|(c<<6),sum=0;for(int b=0;b<K;b++)if((B[b]&mask)==B[b])sum+=table[b*N+code];return sum&3;};
 for(int code=0;code<N;code++){int a=(code*23+7)&63,c=(code*53+19)&63;if(eval(code,a,c)!=(count_hp(code,a,c)&3))throw runtime_error("one literal query per tournament");queries++;}
 set<int>controls={0,32767,21845,10922};for(int j=0;j<64;j++)controls.insert((j*7919+17)&32767);
 for(int code:controls){vector<int>values(4096,0);for(int m:B){values[m]=count_hp(code,m&63,m>>6);queries++;}
  for(int m:B){int coeff=0;for(int s=m;;s=(s-1)&m){coeff+=((__builtin_popcount(m^s)&1)?-1:1)*values[s];if(!s)break;}auto k=lower_bound(B.begin(),B.end(),m)-B.begin();if(((coeff%4)+4)%4!=table[k*N+code])throw runtime_error("literal interpolation coefficient");coefficients++;}}
 for(int m=0;m<4096;m++){int code=24979;if(eval(code,m&63,m>>6)!=(count_hp(code,m&63,m>>6)&3))throw runtime_error("full truth control");queries++;}
 // Exact vector comparison, not a checksum or a hash-collision assumption.
 set<vector<unsigned char>>vectors;for(int code=0;code<N;code++){vector<unsigned char>v;for(int b=0;b<K;b++)v.push_back(table[b*N+code]);if(!vectors.insert(v).second)throw runtime_error("duplicate full modulo-four vector");}
 cout<<"{\"status\":\"PASS\",\"literal_integer_dp_queries\":"<<queries<<",\"full_interpolation_tournaments\":"<<controls.size()<<",\"independently_interpolated_coefficients\":"<<coefficients<<",\"exact_distinct_vectors\":"<<vectors.size()<<",\"binary_bytes\":"<<table.size()<<"}\n";
}
