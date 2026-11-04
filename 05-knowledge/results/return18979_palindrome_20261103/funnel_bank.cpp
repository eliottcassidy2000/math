// Complete unique-start grammar: graded path-bin recurrence, no subset DP.
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
#include <limits>
using namespace std;using I=int64_t;using W=__int128_t;
I checked(W z){if(z<0||z>numeric_limits<I>::max())throw runtime_error("count overflow or negative");return (I)z;}
I power(I b,int e){I z=1;while(e--)z=checked(W(z)*b);return z;}
array<I,3> tail(int t,int bits,int q){
 I L=power(q-1,t),D=L;array<I,6>s{};s[0]=1;
 for(int i=0;i<t;i++){
  bool b=bits>>i&1;
  if(b)L+=power(q,i)*power(q-1,t-i-1);else D+=power(q-1,i)*power(q,t-i-1);
  array<I,6>v{};
  for(int l=0;l<3;l++)for(int a=0;a<2;a++){
   I x=s[2*l+a];v[2*l+a]+=(q-1)*x;v[2*(b?1:2)+a]+=x;if(a||!b)v[2*l+1]+=x;
  }s=v;
 }return {L,D,s[0]+s[1]+s[2]+s[3]};
}
struct Bank{int k=0;I count=0;vector<I>v;I get(I i,int q,int j)const{return v[(i*(k+1)+q)*3+j];}};
int main(int argc,char**argv){
 int B=argc>1?stoi(argv[1]):18;vector<I>targets={18977,18979};vector<Bank>A(B+1);
 A[1].k=B+1;A[1].count=1;A[1].v.resize(3*(B+2));
 for(int q=1;q<=B+1;q++){A[1].v[3*q]=q;A[1].v[3*q+1]=A[1].v[3*q+2]=1;}
 vector<I>u(B+1);u[3]=1;if(B>=4)u[4]=3;for(int n=5;n<=B;n++)u[n]=3*u[n-1]-u[n-2];
 for(int n=3;n<=B;n++){
  auto&out=A[n];out.k=B-n+2;out.count=u[n];out.v.reserve(u[n]*3*(out.k+1));I rank=0;
  for(int m=1;m<n;m++)if(m!=2){int t=n-m-1;I words=t?1LL<<(t-1):1;
   for(int w=0;w<words;w++){int bits=t?w|(1<<(t-1)):0;
    vector<array<I,3>>coef(out.k+2);for(int q=1;q<=out.k+1;q++)coef[q]=tail(t,bits,q);
    for(I ix=0;ix<A[m].count;ix++,rank++){
     out.v.insert(out.v.end(),3,0);
     for(int q=1;q<=out.k;q++){
      I P=A[m].get(ix,q,0),E=A[m].get(ix,q,1),S=A[m].get(ix,q,2);
      out.v.push_back(checked(W(q)*(W(coef[q+1][0])*A[m].get(ix,q+1,2)+W(coef[q][2])*P-W(coef[q][1])*E)));
      out.v.push_back(checked(W(coef[q][0])*P-W(power(q-1,t))*E));
      out.v.push_back(checked(W(power(q,t))*S+W(coef[q][1])*A[m].get(ix,q-1,0)));
     }
     I h=out.get(rank,1,0),F=out.get(rank,2,1);if(out.get(rank,1,2)!=1)throw runtime_error("unique start");
     if(argc>2 && (n<=9 || rank==0 || rank==u[n]/2 || rank==u[n]-1)){
      cout<<"{\"kind\":\"query\",\"n\":"<<n<<",\"rank\":"<<rank<<",\"values\":[";
      for(int q=1;q<=2;q++)for(int j=0;j<3;j++){if(q!=1||j)cout<<',';cout<<out.get(rank,q,j);}cout<<"]}\n";
     }
     for(I H:targets)if(h+F-1==H){cout<<"{\"kind\":\"hit\",\"n\":"<<n<<",\"rank\":"<<rank<<",\"target\":"<<H<<",\"h\":"<<h<<",\"F\":"<<F<<",\"core_order\":"<<m<<",\"core_rank\":"<<ix<<",\"tail_length\":"<<t<<",\"tail_bits\":"<<bits<<"}\n";cout.flush();}
    }
   }
  }
  if(rank!=u[n])throw runtime_error("grammar census mismatch");
  cout<<"{\"kind\":\"census\",\"n\":"<<n<<",\"classes\":"<<rank<<"}\n";cout.flush();
 }
 return 0;
}
