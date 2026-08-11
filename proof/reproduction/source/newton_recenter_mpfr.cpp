#include <boost/multiprecision/mpfr.hpp>
#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
namespace bm=boost::multiprecision;
using R=bm::number<bm::backends::mpfr_float_backend<80>>;
constexpr int Q=8;
using V=std::array<R,Q>;
using M=std::array<std::array<R,Q>,Q>;
using D3=std::array<std::array<R,Q>,3>;
using P3=std::array<R,3>;
struct Mesh{int n;R tau,rtau,rstate;std::vector<R>den;std::vector<V>node;};
static R rd(std::istream&f){std::string s;f>>s;return R(s);}
static Mesh readMesh(const char*p){
  std::ifstream f(p);if(!f)throw std::runtime_error("mesh open");
  Mesh m;f>>m.n;m.tau=rd(f);m.rtau=rd(f);m.rstate=rd(f);
  m.den.resize(m.n);for(auto&x:m.den)x=rd(f);
  m.node.resize(m.n+1);for(auto&v:m.node)for(auto&x:v)x=rd(f);
  if(!f)throw std::runtime_error("mesh truncated");return m;
}
struct Rec{V x;M mat;};
static Rec readRec(std::istream&f,int want){
  std::string t;int idx,ok;f>>t>>idx>>t>>ok;
  if(!f||t!="OK"||idx!=want||ok!=1)throw std::runtime_error("bad REC");
  f>>t;for(int i=0;i<5;i++)rd(f);
  Rec r;f>>t;if(t!="XC")throw std::runtime_error("XC");
  for(auto&x:r.x){R l=rd(f),u=rd(f);x=(l+u)/2;}
  f>>t;if(t!="XB")throw std::runtime_error("XB");
  for(int i=0;i<Q;i++){rd(f);rd(f);}
  f>>t;if(t!="MB")throw std::runtime_error("MB");
  for(auto&row:r.mat)for(auto&x:row){R l=rd(f),u=rd(f);x=(l+u)/2;}
  return r;
}
static std::vector<Rec> readLog(const char*p,int n){
  std::ifstream f(p);if(!f)throw std::runtime_error("log open");
  std::vector<Rec>r;for(int k=0;k<n;k++)r.push_back(readRec(f,k));return r;
}
static V field(const V&x){
  const R mu=R(1001)/1000;
  R q3x=-(x[0]+x[2])/mu,q3y=-(x[1]+x[3])/mu;
  auto g=[](R x,R y){R r2=x*x+y*y,r3=r2*sqrt(r2);return std::array<R,2>{x/r3,y/r3};};
  auto a=g(x[2]-x[0],x[3]-x[1]),b=g(q3x-x[0],q3y-x[1]);
  auto c=g(x[0]-x[2],x[1]-x[3]),d=g(q3x-x[2],q3y-x[3]);
  V f{x[4],x[5],x[6],x[7],a[0]+mu*b[0],a[1]+mu*b[1],c[0]+mu*d[0],c[1]+mu*d[1]};
  return f;
}
static D3 boundaryD(const V&x){
  D3 d{};d[0][0]=1;d[0][2]=1;d[1][1]=1;d[1][3]=1;
  d[2][0]=x[4]-x[6];d[2][1]=x[5]-x[7];d[2][4]=x[0];d[2][5]=x[1];d[2][6]=-x[0];d[2][7]=-x[1];
  return d;
}
static P3 boundary(const V&x){return {x[0]+x[2],x[1]+x[3],x[0]*(x[4]-x[6])+x[1]*(x[5]-x[7])};}
static std::array<std::array<R,3>,3> inv3(std::array<std::array<R,3>,3>a){
  std::array<std::array<R,6>,3>w{};
  for(int i=0;i<3;i++){for(int j=0;j<3;j++)w[i][j]=a[i][j];w[i][3+i]=1;}
  for(int k=0;k<3;k++){int p=k;for(int i=k+1;i<3;i++)if(abs(w[i][k])>abs(w[p][k]))p=i;
    if(p!=k)std::swap(w[p],w[k]);R v=w[k][k];
    for(int j=0;j<6;j++)w[k][j]/=v;
    for(int i=0;i<3;i++)if(i!=k){v=w[i][k];for(int j=0;j<6;j++)w[i][j]-=v*w[k][j];}}
  std::array<std::array<R,3>,3>b{};for(int i=0;i<3;i++)for(int j=0;j<3;j++)b[i][j]=w[i][3+j];return b;
}
static V mv(const M&a,const V&x){V y{};for(int i=0;i<Q;i++)for(int j=0;j<Q;j++)y[i]+=a[i][j]*x[j];return y;}
int main(int ac,char**av){
 try{
  if(ac!=5){std::cerr<<"mesh highcenter.log jac.log outmesh\n";return 2;}
  Mesh mesh=readMesh(av[1]);auto rec=readLog(av[2],mesh.n);auto jacrec=readLog(av[3],mesh.n);int N=mesh.n,D=8*N+3;
  std::vector<M> mm(N);std::vector<V>cc(N);
  for(int k=0;k<N;k++){mm[k]=jacrec[k].mat;V f=field(rec[k].x);for(int i=0;i<Q;i++)cc[k][i]=f[i]/mesh.den[k];}
  D3 db=boundaryD(mesh.node[N]);
  std::array<std::array<R,3>,Q>P{};P[4][0]=1;P[5][1]=1;P[6][0]=1;P[7][1]=1;
  for(int k=0;k<N;k++){std::array<std::array<R,3>,Q>z{};
    for(int i=0;i<Q;i++)for(int j=0;j<3;j++){for(int l=0;l<Q;l++)z[i][j]+=mm[k][i][l]*P[l][j];if(j==2)z[i][j]+=cc[k][i];}P=z;}
  std::array<std::array<R,3>,3>S{};
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int l=0;l<Q;l++)S[i][j]+=db[i][l]*P[l][j];
  auto V3=inv3(S);
  std::vector<V> y(N);R maxF=0;
  for(int k=0;k<N;k++)for(int i=0;i<Q;i++){y[k][i]=rec[k].x[i]-mesh.node[k+1][i];maxF=std::max(maxF,R(abs(y[k][i])));}
  P3 yb=boundary(mesh.node[N]);for(auto z:yb)maxF=std::max(maxF,R(abs(z)));
  V w{};for(int k=0;k<N;k++){V z=mv(mm[k],w);for(int i=0;i<Q;i++)z[i]-=y[k][i];w=z;}
  P3 rhs=yb;for(int i=0;i<3;i++)for(int j=0;j<Q;j++)rhs[i]-=db[i][j]*w[j];
  P3 p{};for(int i=0;i<3;i++)for(int j=0;j<3;j++)p[i]+=V3[i][j]*rhs[j];
  std::vector<R> sol(D);for(int i=0;i<3;i++)sol[i]=p[i];
  V cur{};cur[4]=p[0];cur[5]=p[1];cur[6]=p[0];cur[7]=p[1];
  for(int k=0;k<N;k++){V z=mv(mm[k],cur);for(int i=0;i<Q;i++)z[i]+=cc[k][i]*p[2]-y[k][i];for(int i=0;i<Q;i++)sol[3+Q*k+i]=z[i];cur=z;}
  R maxC=0;for(auto z:sol)maxC=std::max(maxC,R(abs(z)));
  R anew=mesh.node[0][4]-sol[0],bnew=mesh.node[0][5]-sol[1],tnew=mesh.tau-sol[2];
  std::vector<V> node=mesh.node;node[0]={R(-1),R(0),R(1),R(0),anew,bnew,anew,bnew};
  for(int k=1;k<=N;k++)for(int i=0;i<Q;i++)node[k][i]-=sol[3+Q*(k-1)+i];
  std::ofstream o(av[4]);o<<std::scientific<<std::setprecision(75);
  o<<N<<" "<<tnew<<" "<<R("1e-13")<<" "<<R("1e-14")<<"\n";
  for(auto x:mesh.den)o<<x<<" ";o<<"\n";
  for(auto v:node){for(auto x:v)o<<x<<" ";o<<"\n";}
  std::cout<<std::scientific<<std::setprecision(40);
  std::cout<<"MAX_F "<<maxF<<"\nMAX_CORRECTION "<<maxC<<"\nDP "<<-sol[0]<<" "<<-sol[1]<<" "<<-sol[2]<<"\n";
  std::cout<<"NEW_P "<<anew<<" "<<bnew<<" "<<tnew<<"\nSCHUR ";
  for(auto r:S)for(auto x:r)std::cout<<x<<" ";std::cout<<"\n";
  return 0;
 }catch(std::exception&e){std::cerr<<"ERROR "<<e.what()<<"\n";return 99;}
}