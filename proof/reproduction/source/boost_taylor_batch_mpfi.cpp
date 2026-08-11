#pragma STDC FENV_ACCESS ON
#include <boost/multiprecision/mpfi.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <array>
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
namespace bn=boost::multiprecision;
using R=bn::number<bn::backends::mpfr_float_backend<80>>;
using I=bn::number<bn::backends::mpfi_float_backend<80>>;
static R MW_TOL("1e-8");
static I rangeI(const R&l,const R&u){I x;mpfi_interv_fr(x.backend().data(),l.backend().data(),u.backend().data());return x;}
constexpr int N=8, AUG=72, PMAX=24;
using V=std::array<I,N>;
using Mat=std::array<std::array<I,N>,N>;
using AV=std::array<I,AUG>;
using S=std::array<I,PMAX+1>;
using CS=std::array<S,AUG>;
using SM=std::array<std::array<S,N>,N>;
static inline R lo(const I&x){return bn::lower(x);}
static inline R hi(const I&x){return bn::upper(x);}
static inline R mag(const I&x){return std::max(bn::abs(lo(x)),bn::abs(hi(x)));}
static inline R wid(const I&x){return hi(x)-lo(x);}
static I sqrI(const I&x){
  if(lo(x)>=0 || hi(x)<=0) return x*x;
  R m=std::max(-lo(x),hi(x));
  I z=rangeI(R(0),m); return z*z;
}
static I muI(){return I(1001)/I(1000);}
static S zs(){S a; for(auto&x:a)x=I(0); return a;}
static S sadd(const S&a,const S&b,int n){S c=zs();for(int k=0;k<=n;k++)c[k]=a[k]+b[k];return c;}
static S ssub(const S&a,const S&b,int n){S c=zs();for(int k=0;k<=n;k++)c[k]=a[k]-b[k];return c;}
static S sscale(const S&a,const I&b,int n){S c=zs();for(int k=0;k<=n;k++)c[k]=a[k]*b;return c;}
static S smul(const S&a,const S&b,int n){
  S c=zs(); for(int k=0;k<=n;k++){I t(0);for(int j=0;j<=k;j++)t+=a[j]*b[k-j];c[k]=t;}return c;
}
static S ssqrt(const S&a,int n){
  S b=zs(); if(!(lo(a[0])>0))throw std::runtime_error("sqrt domain");
  b[0]=sqrt(a[0]);
  for(int k=1;k<=n;k++){I t(0);for(int j=1;j<k;j++)t+=b[j]*b[k-j];b[k]=(a[k]-t)/(I(2)*b[0]);}
  return b;
}
static S sinv(const S&a,int n){
  S b=zs(); if(lo(a[0])<=0 && hi(a[0])>=0)throw std::runtime_error("inv domain");
  b[0]=I(1)/a[0];
  for(int k=1;k<=n;k++){I t(0);for(int j=1;j<=k;j++)t+=a[j]*b[k-j];b[k]=-t/a[0];}
  return b;
}
struct PairSeries {S gx,gy,dxx,dxy,dyy;};
static PairSeries pairseries(const S&rx,const S&ry,int n){
  S r2=sadd(smul(rx,rx,n),smul(ry,ry,n),n);
  S root=ssqrt(r2,n), inv=sinv(root,n);
  S inv2=smul(inv,inv,n), inv3=smul(inv2,inv,n), inv5=smul(inv3,inv2,n);
  PairSeries p;
  p.gx=smul(rx,inv3,n); p.gy=smul(ry,inv3,n);
  S xx5=smul(smul(rx,rx,n),inv5,n);
  S xy5=smul(smul(rx,ry,n),inv5,n);
  S yy5=smul(smul(ry,ry,n),inv5,n);
  p.dxx=ssub(inv3,sscale(xx5,I(3),n),n);
  p.dxy=sscale(xy5,I(-3),n);
  p.dyy=ssub(inv3,sscale(yy5,I(3),n),n);
  return p;
}
struct DynSeries {std::array<S,4> acc; SM J;};
static DynSeries dynseries(const CS&C,int n){
  I mu=muI(), one(1);
  S q1x=C[0],q1y=C[1],q2x=C[2],q2y=C[3];
  S q3x=sscale(sadd(q1x,q2x,n),-I(1)/mu,n);
  S q3y=sscale(sadd(q1y,q2y,n),-I(1)/mu,n);
  auto p12=pairseries(ssub(q2x,q1x,n),ssub(q2y,q1y,n),n);
  auto p13=pairseries(ssub(q3x,q1x,n),ssub(q3y,q1y,n),n);
  auto p23=pairseries(ssub(q3x,q2x,n),ssub(q3y,q2y,n),n);
  DynSeries d;
  for(auto&a:d.acc)a=zs();
  for(auto&r:d.J)for(auto&a:r)a=zs();
  d.acc[0]=sadd(p12.gx,sscale(p13.gx,mu,n),n);
  d.acc[1]=sadd(p12.gy,sscale(p13.gy,mu,n),n);
  d.acc[2]=sadd(sscale(p12.gx,I(-1),n),sscale(p23.gx,mu,n),n);
  d.acc[3]=sadd(sscale(p12.gy,I(-1),n),sscale(p23.gy,mu,n),n);
  d.J[0][4][0]=one; d.J[1][5][0]=one; d.J[2][6][0]=one; d.J[3][7][0]=one;
  S A[2][2]={{p12.dxx,p12.dxy},{p12.dxy,p12.dyy}};
  S B[2][2]={{p13.dxx,p13.dxy},{p13.dxy,p13.dyy}};
  S D[2][2]={{p23.dxx,p23.dxy},{p23.dxy,p23.dyy}};
  for(int i=0;i<2;i++)for(int j=0;j<2;j++){
    d.J[4+i][j]=ssub(sscale(A[i][j],I(-1),n),sscale(B[i][j],mu+one,n),n);
    d.J[4+i][2+j]=ssub(A[i][j],B[i][j],n);
    d.J[6+i][j]=ssub(A[i][j],D[i][j],n);
    d.J[6+i][2+j]=ssub(sscale(A[i][j],I(-1),n),sscale(D[i][j],mu+one,n),n);
  }
  return d;
}
static CS coefficients(const AV&y0,const I&scale,int order){
  CS C; for(auto&s:C)for(auto&x:s)x=I(0);
  for(int i=0;i<AUG;i++)C[i][0]=y0[i];
  for(int n=0;n<order;n++){
    DynSeries d=dynseries(C,n);
    I den(n+1);
    C[0][n+1]=scale*C[4][n]/den; C[1][n+1]=scale*C[5][n]/den;
    C[2][n+1]=scale*C[6][n]/den; C[3][n+1]=scale*C[7][n]/den;
    for(int i=0;i<4;i++)C[4+i][n+1]=scale*d.acc[i][n]/den;
    for(int i=0;i<N;i++)for(int j=0;j<N;j++){
      I t(0);
      for(int k=0;k<=n;k++)for(int l=0;l<N;l++)t+=d.J[i][l][k]*C[N+l*N+j][n-k];
      C[N+i*N+j][n+1]=scale*t/den;
    }
  }
  return C;
}
static V field(const V&x){
  CS C;for(auto&s:C)for(auto&z:s)z=I(0);for(int i=0;i<N;i++)C[i][0]=x[i];
  auto d=dynseries(C,0);V f;
  f[0]=x[4];f[1]=x[5];f[2]=x[6];f[3]=x[7];
  for(int i=0;i<4;i++)f[4+i]=d.acc[i][0];
  return f;
}
static Mat jac(const V&x){
  CS C;for(auto&s:C)for(auto&z:s)z=I(0);for(int i=0;i<N;i++)C[i][0]=x[i];
  auto d=dynseries(C,0);Mat J;for(auto&r:J)for(auto&z:r)z=I(0);
  for(int i=0;i<N;i++)for(int j=0;j<N;j++)J[i][j]=d.J[i][j][0];
  return J;
}
static std::array<I,3> dist2(const V&x){
  I mu=muI();
  I q3x=-(x[0]+x[2])/mu,q3y=-(x[1]+x[3])/mu;
  return {sqrI(x[0]-x[2])+sqrI(x[1]-x[3]),
          sqrI(x[0]-q3x)+sqrI(x[1]-q3y),
          sqrI(x[2]-q3x)+sqrI(x[3]-q3y)};
}
static I inflate(const I&a,long double factor){
  R l=lo(a),u=hi(a);
  I ew=(I(u)-I(l))*I(R(factor-1)/R(2));
  R extra=hi(ew);if(extra<R("1e-50"))extra=R("1e-50");
  return a+rangeI(-extra,extra);
}
static bool strictin(const I&a,const I&b){return lo(b)<lo(a)&&hi(a)<hi(b);}
static bool stateTube(const V&x0,const I&scale,long double du,V&W,R&margin,R&mind2){
  try{
    V f0=field(x0);I tt=rangeI(R(0),R(du));
    for(int i=0;i<N;i++)W[i]=inflate(bn::hull(x0[i],I(x0[i]+tt*scale*f0[i])),1.5L);
    for(int it=0;it<30;it++){
      auto dd=dist2(W);mind2=std::min({lo(dd[0]),lo(dd[1]),lo(dd[2])});
      if(!(mind2>0))return false;
      V fw=field(W),P;bool ok=true;margin=R("1e300");
      for(int i=0;i<N;i++){
        P[i]=x0[i]+tt*scale*fw[i];
        ok=ok&&strictin(P[i],W[i]);
        margin=std::min(margin,std::min(lo(P[i])-lo(W[i]),hi(W[i])-hi(P[i])));
      }
      if(ok)return true;
      for(int i=0;i<N;i++)W[i]=inflate(bn::hull(x0[i],I(P[i])),1.35L);
    }
  }catch(...){return false;}
  return false;
}
static Mat matmul(const Mat&A,const Mat&B){
  Mat C;for(auto&r:C)for(auto&z:r)z=I(0);
  for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int k=0;k<N;k++)C[i][j]+=A[i][k]*B[k][j];
  return C;
}
static CS stateCoefficients(const V&x0,const I&scale,int order){
  CS C;for(auto&s:C)for(auto&z:s)z=I(0);for(int i=0;i<N;i++)C[i][0]=x0[i];
  for(int n=0;n<order;n++){
    DynSeries d=dynseries(C,n);I den(n+1);
    C[0][n+1]=scale*C[4][n]/den;C[1][n+1]=scale*C[5][n]/den;
    C[2][n+1]=scale*C[6][n]/den;C[3][n+1]=scale*C[7][n]/den;
    for(int i=0;i<4;i++)C[4+i][n+1]=scale*d.acc[i][n]/den;
  }
  return C;
}
struct StateStepOut{V x,tube;R margin,mind2;};
static bool stateTaylorStep(const V&x0,const I&scale,long double du,int p,StateStepOut&o){
  V W;if(!stateTube(x0,scale,du,W,o.margin,o.mind2))return false;
  CS c0,cw;try{c0=stateCoefficients(x0,scale,p-1);cw=stateCoefficients(W,scale,p);}catch(...){return false;}
  I h(du),pw(1);V end;for(auto&z:end)z=I(0);
  for(int k=0;k<p;k++){if(k>0)pw*=h;for(int i=0;i<N;i++)end[i]+=c0[i][k]*pw;}
  pw*=h;for(int i=0;i<N;i++)end[i]+=cw[i][p]*pw;
  for(int i=0;i<N;i++){o.x[i]=end[i];o.tube[i]=W[i];}
  return true;
}

struct StepOut{V x;Mat M;V tube;R margin,mind2,maxmw;};
static bool taylorStep(const V&x0,const I&scale,long double du,int p,StepOut&o){
  V W;if(!stateTube(x0,scale,du,W,o.margin,o.mind2))return false;
  Mat JW=jac(W);R L(0);
  for(int i=0;i<N;i++){R row(0);for(int j=0;j<N;j++)row+=mag(scale*JW[i][j]);L=std::max(L,row);}
  I eb=exp(rangeI(R(0),R(du))*I(L));R B=hi(eb);
  if(!std::isfinite((double)B)||B>1e100L)return false;
  AV y0,yw;for(auto&z:y0)z=I(0);for(auto&z:yw)z=I(0);
  for(int i=0;i<N;i++){y0[i]=x0[i];yw[i]=W[i];}
  for(int i=0;i<N;i++)for(int j=0;j<N;j++){
    y0[N+i*N+j]=I(i==j?1:0); yw[N+i*N+j]=rangeI(-B,B);
  }
  CS c0,cw;try{c0=coefficients(y0,scale,p-1);cw=coefficients(yw,scale,p);}catch(...){return false;}
  I h(du),pw(1);AV end;for(auto&z:end)z=I(0);
  for(int k=0;k<p;k++){if(k>0)pw*=h;for(int i=0;i<AUG;i++)end[i]+=c0[i][k]*pw;}
  pw*=h;for(int i=0;i<AUG;i++)end[i]+=cw[i][p]*pw;
  for(int i=0;i<N;i++){o.x[i]=end[i];o.tube[i]=W[i];}
  o.maxmw=0;
  for(int i=0;i<N;i++)for(int j=0;j<N;j++){o.M[i][j]=end[N+i*N+j];o.maxmw=std::max(o.maxmw,wid(o.M[i][j]));}
  return true;
}
struct FlowOut{V x;Mat M;R mind2,minmargin,maxlocalmw;int steps,minpow;};
static bool flow(V X,const I&scale,int p,FlowOut&o){
  R sm=(lo(scale)+hi(scale))/2; I sc(sm);
  V cen;for(int i=0;i<N;i++)cen[i]=I((lo(X[i])+hi(X[i]))/2);
  Mat A;for(int i=0;i<N;i++)for(int j=0;j<N;j++)A[i][j]=I(i==j?1:0);
  const int DEN=16384;int pos=0,trytick=2048;o.steps=0;o.mind2=R("1e300");o.minmargin=R("1e300");o.maxlocalmw=R(0);o.minpow=30;
  while(pos<DEN){
    int ticks=std::min(trytick,DEN-pos);while(ticks&(ticks-1))ticks&=ticks-1;
    StepOut pc,pb;bool ok=false;V Y,newc;
    while(ticks>=1){
      long double du=std::ldexp((long double)ticks,-14);
      if(taylorStep(cen,sc,du,p,pc)&&taylorStep(X,scale,du,p,pb)){
        V ff;try{ff=field(pb.tube);}catch(...){ticks/=2;continue;}
        I ds=scale-sc;
        for(int i=0;i<N;i++){
          I yy=pc.x[i];
          for(int j=0;j<N;j++)yy+=pb.M[i][j]*(X[j]-cen[j]);
          yy+=I(du)*ff[i]*ds;
          Y[i]=yy;
          R m=(lo(yy)+hi(yy))/2;newc[i]=I(m);
        }
        bool finite=true;for(auto&z:Y)finite=finite&&std::isfinite((double)lo(z))&&std::isfinite((double)hi(z));
        if(finite && pb.maxmw<1e-8L && pc.maxmw<1e-8L){ok=true;break;}
      }
      ticks/=2;
    }
    if(!ok){auto dx=dist2(X);R mw(0);for(auto&z:X)mw=std::max(mw,wid(z));std::cerr<<"FAILPOS "<<pos<<" MW "<<mw<<" D "<<lo(dx[0])<<" "<<lo(dx[1])<<" "<<lo(dx[2])<<" X0 "<<lo(X[0])<<" "<<hi(X[0])<<"\n";return false;}
    X=Y;cen=newc;A=matmul(pb.M,A);pos+=ticks;o.steps++;
    o.mind2=std::min(o.mind2,pb.mind2);o.minmargin=std::min(o.minmargin,pb.margin);o.maxlocalmw=std::max(o.maxlocalmw,pb.maxmw);
    int lg=0,t=ticks;while(t>1){t>>=1;lg++;}o.minpow=std::min(o.minpow,14-lg);
    if(ticks<2048)trytick=std::min(2048,ticks*2);else trytick=2048;
  }
  o.x=X;o.M=A;return true;
}
struct PairFlowOut{V c,b;Mat M;R mind2,minmargin,maxlocalmw;int steps,minpow;};
static bool containsI(const I&a,const I&b){return lo(a)<=lo(b)&&hi(b)<=hi(a);}
static bool flowPair(V Cset,V Bset,const I&scale,int p,PairFlowOut&o){
  R sm=(lo(scale)+hi(scale))/2;I sc(sm);
  V cen;for(int i=0;i<N;i++)cen[i]=I((lo(Cset[i])+hi(Cset[i]))/2);
  Mat A;for(int i=0;i<N;i++)for(int j=0;j<N;j++)A[i][j]=I(i==j?1:0);
  const int DEN=16384;int pos=0,trytick=2048;o.steps=0;o.mind2=R("1e300");o.minmargin=R("1e300");o.maxlocalmw=R(0);o.minpow=30;
  while(pos<DEN){
    int ticks=std::min(trytick,DEN-pos);while(ticks&(ticks-1))ticks&=ticks-1;
    StateStepOut pc;StepOut pb;bool ok=false;V CY,BY,newc;
    while(ticks>=1){
      long double du=std::ldexp((long double)ticks,-14);
      if(stateTaylorStep(cen,sc,du,p,pc)&&taylorStep(Bset,scale,du,p,pb)){
        V ff;try{ff=field(pb.tube);}catch(...){ticks/=2;continue;}
        I ds=scale-sc;
        for(int i=0;i<N;i++){
          I cy=pc.x[i],by=pc.x[i];
          for(int j=0;j<N;j++){cy+=pb.M[i][j]*(Cset[j]-cen[j]);by+=pb.M[i][j]*(Bset[j]-cen[j]);}
          by+=I(du)*ff[i]*ds;
          CY[i]=cy;BY[i]=by;newc[i]=I((lo(cy)+hi(cy))/2);
        }
        bool finite=true,cont=true;
        for(int i=0;i<N;i++){finite=finite&&std::isfinite((double)lo(BY[i]))&&std::isfinite((double)hi(BY[i]));cont=cont&&containsI(BY[i],CY[i]);}
        if(finite&&cont&&pb.maxmw<MW_TOL){ok=true;break;}
      }
      ticks/=2;
    }
    if(!ok){std::cerr<<"FAILPOS "<<pos<<"\n";return false;}
    Cset=CY;Bset=BY;cen=newc;A=matmul(pb.M,A);pos+=ticks;o.steps++;
    o.mind2=std::min(o.mind2,pb.mind2);o.minmargin=std::min(o.minmargin,pb.margin);o.maxlocalmw=std::max(o.maxlocalmw,pb.maxmw);
    int lg=0,t=ticks;while(t>1){t>>=1;lg++;}o.minpow=std::min(o.minpow,14-lg);
    if(ticks<2048)trytick=std::min(2048,ticks*2);else trytick=2048;
  }
  o.c=Cset;o.b=Bset;o.M=A;return true;
}

static R parse(const std::string&s){return R(s);}
static I box(const R&c,const R&r){return I(c)+rangeI(-r,r);}
int main(){
  std::ios::sync_with_stdio(false);std::cin.tie(nullptr);
  int mode,count,p;std::string stau;
  std::string srtau,srstate,smwtol;
  if(!(std::cin>>mode>>stau>>srtau>>srstate>>p>>smwtol>>count))return 9;
  R rtau=parse(srtau),rstate=parse(srstate),mwtol=parse(smwtol);
  R tau=parse(stau);MW_TOL=mwtol;
  std::cout<<std::setprecision(90)<<std::scientific;
  R globalD2("1e300");
  for(int rec=0;rec<count;rec++){
    std::string sden;std::cin>>sden;R denom=parse(sden);
    V xc;for(int i=0;i<N;i++){std::string z;std::cin>>z;xc[i]=I(parse(z));}
    V xb;for(int i=0;i<N;i++)xb[i]=box(lo(xc[i]),rstate);
    I sb=box(tau,rtau)/I(denom);
    PairFlowOut out;bool ok=flowPair(xc,xb,sb,p,out);
    std::cout<<"REC "<<rec<<" OK "<<ok<<"\n";
    if(!ok)return 20;
    globalD2=std::min(globalD2,out.mind2);
    std::cout<<"META "<<out.steps<<" "<<out.minpow<<" "<<out.mind2<<" "<<out.minmargin<<" "<<out.maxlocalmw<<"\n";
    std::cout<<"XC";for(int i=0;i<N;i++)std::cout<<" "<<lo(out.c[i])<<" "<<hi(out.c[i]);
    std::cout<<"\nXB";for(int i=0;i<N;i++)std::cout<<" "<<lo(out.b[i])<<" "<<hi(out.b[i]);
    std::cout<<"\nMB";for(int i=0;i<N;i++)for(int j=0;j<N;j++)std::cout<<" "<<lo(out.M[i][j])<<" "<<hi(out.M[i][j]);
    std::cout<<"\n";
  }
  std::cout<<"BATCH_OK 1 COUNT "<<count<<" GLOBAL_D2 "<<globalD2<<"\n";
  return 0;
}
