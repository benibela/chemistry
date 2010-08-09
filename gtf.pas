{Dieser Programmteil enthält die Datenstrukturen für Gauss Type Functions
 und Berechnungsfunktionen für die dabei möglichen Integrale

 Die Funktionen folgen Cook, Ab Initio Valence Calculations in Chemistry
 und Carsky/Urban, Ab Initio Calculations
 }
unit gtf;

{$mode objfpc}{$H+}
{$define performtests}

interface

uses
  Classes, SysUtils,math,extmath,Dialogs;

type
  //Datenstruktur für eine allgemeine Gauss Type Function
  TGTF=record
    e:array['x'..'z'] of longint; //Exponenten
    n:TVector3n;                  //KO-System Nullpunkt
    alpha: number;                //Exponent a
    pre: number;                  //Vorfaktor b
    //Funktion der Form
    {        xe     ye    ze   -a*r^2
      b*(x-x0) (y-y0) (z-z0)  e
    }
  end;

  //Unskalierte GTF am Koordinatennullpunkt
  TFreeGTF=record
    e:array['x'..'z'] of longint; //Exponenten
    alpha: number;  //Exponent a
    //Funktion der Form
    {  xe  ye  ze  -a*r^2
      x   y   z   e
    }
  end;

  //Datenstruktur für ein vorberechnetes Überlappungsintegral
  TGTFPreCalculatedOverlap=record
    factor,                 //Vorfaktor
    alphasum,               //alpha1 + g.alpha2
    alphamulprosum: number; //alpha1*alpha2/alphasum
    f,g:TFreeGTF;           //freie Orbitale
  end;

  //Datenstruktur für ein vorberechnetes kinetisches Energieintegral
  TGTFPreCalculatedKinetic=record
    //Zusammengesetzt aus overlaps
    overlap: array[0..3] of TGTFPreCalculatedOverlap;
  end;

  //Datenstruktur für ein vorberechnetes Kernanziehungsintegral
  TGTFPreCalculatedNuclearAttraction=record
    factor: number;                //Vorfaktor
    f,g:TFreeGTF;                  //freie Orbitalfunktionen
    coeffsum: longint;             //Die Summe alle Koeffizienten (fi.e[c])
    a_cache:array['x'..'z',0..4,0..2,0..2] of number; //Vorberechnete Faktoren
  end;

  //Datenstruktur für ein vorberechnetes Elektronenwechselwirkungsintegral
  TGTFPreCalculatedElectronRepulsion=record
    factor: number;
    f,g: TFreeGTF;
    alphasum: number; //Die Summe der alphas und Kehrwerte
    coeffsum: longint;//Die Summe alle Koeffizienten (fi.e[c])
    b_cache: array['x'..'z',0..4,0..4,0..2,0..2,0..4] of number; //Vorberechnung
    a_power: array[-4..0] of number; //Potenzen von alphasum
  end;
  
//Berechnet das Integral
//  1     2v  -ts^2
// S  ds s   e
//  0
function f_erf_like(const v:longint;const t,erf_sqrt_t:number):number;overload;
function f_erf_like(const v:longint;const t:number):number;inline;overload;

//Berechnet den Mittelpunkt zweier GTF-Orbitale
function midpoint(const a1,a2:number;const p1,p2:TVector3n):TVector3n;inline;

//Berechnet einen normalisiertenden Vorfaktor
function normalizeFactorGTF(const l,m,n:longint;const alpha:number): number;
         inline;overload;

//Berechnet das Überlappungsintegral
function preCalculateOverlapIntegral(const factor:number;const f,g:TFreeGTF):
         TGTFPreCalculatedOverlap;
function overlapIntegral(const pre: TGTFPreCalculatedOverlap;
                         const p1,p2: TVector3n): number;
//Berechnet das Integral der kinetischen Energie
function preCalculateKineticIntegral(const factor:number;const f,g:TFreeGTF):
         TGTFPreCalculatedKinetic;
function kineticIntegral(const pre: TGTFPreCalculatedKinetic;
                         const p1,p2: TVector3n): number;
//Berechnet die Energie der Anziehung zwischen Kern und Orbital
function preCalculateNuclearAttractionIntegral(const factor:number;
                              const f,g:TFreeGTF):TGTFPreCalculatedNuclearAttraction;
function nuclearAttractionIntegral(const pre:TGTFPreCalculatedNuclearAttraction;
                                   const p1,p2,kern: TVector3n):number;
//Berechnet die durch die Abstoßung zwischen Elektronen gegebene Energie
function preCalculateElectronRepulsionIntegral(const factor:number;
                             const f,g:TFreeGTF):TGTFPreCalculatedElectronRepulsion;
function electronRepulsionIntegral(const pre12,pre34:
              TGTFPreCalculatedElectronRepulsion;const p1,p2,p3,p4:TVector3n):number;
implementation


//Berechnet das Integral
//  1     2v  -ts^2
// S  ds s   e
//  0
//nach http://mathworld.wolfram.com/GaussianIntegral.html
//     und Beispielsrechnungen mit Derive
function f_erf_like(const v:longint;const t,erf_sqrt_t:number):number;overload;
var temp: number;
    i:longint;
begin
  if isZero(t) then exit(1/(2*v+1))
  else begin
    assert(t>0,'Negative Zahlen sind in f_erf_like nicht implementiert');
    if v>=1 then begin
      temp:=1/t;
      result:=0;
      for i:=v-1 downto 0 do begin
        result-=LOOK_UP_FAK2[2*v-1]/(LOOK_UP_2POWER[v-i]*LOOK_UP_FAK2[2*i+1])*temp;
        temp/=t;
      end;
      result:=exp(-t)*result;
    end else result:=0;
     result+=LOOK_UP_FAK2[2*v-1]*sqrt(pi)*erf_sqrt_t/((LOOK_UP_2POWER[v+1])*
             power(t,v+0.5));
  end;
end;

function f_erf_like(const v:longint;const t:number):number;inline;overload;
begin
  if isZero(t) then exit(1/(2*v+1))
  else exit(f_erf_like(v,t,erf(sqrt(t))));
end;

//Berechnet den Mittelpunkt zweier GTF-Orbitale als Vektor
//nach Cook, Ab Initio Valence Calculations in Chemistry
function midpoint(const a1,a2:number;const p1,p2:TVector3n):TVector3n;inline;
var a12:number;
begin
  a12:=a1+a2;
  result['x']:=(a1*p1['x']+a2*p2['x'])/a12;
  result['y']:=(a1*p1['y']+a2*p2['y'])/a12;
  result['z']:=(a1*p1['z']+a2*p2['z'])/a12;
end;

//Berechnet einen normalisierenden Vorfaktor für eine GTF
//nach Cook, Ab Initio Valence Calculations in Chemistry und Überprüfungen mit Derive
function normalizeFactorGTF(const l,m,n:longint;const alpha:number): number;
  inline;overload;
begin
  result:=sqrt((intpower(4,l+m+n)*sqrt(8)*power(alpha,l+m+n+1.5))/
              (LOOK_UP_FAK2[l]*LOOK_UP_FAK2[m]*LOOK_UP_FAK2[n]*sqrt(pi*pi*pi)));
end;

//Berechnet das Überlappungsintegral
function preCalculateOverlapIntegral(const factor:number;const f,g:TFreeGTF):
         TGTFPreCalculatedOverlap;
begin
  result.alphasum:=f.alpha+g.alpha;
  result.alphamulprosum:=f.alpha*g.alpha/result.alphasum;
  result.factor:=factor*power(pi/result.alphasum, 1.5);
  result.f:=f;
  result.g:=g;
end;

//Berechnet ein verschobenes Überlappungsintegral
function overlapIntegral(const pre: TGTFPreCalculatedOverlap;
                         const p1,p2: TVector3n): number;
  function asym(const c: char):number;
  const LOOK_UP_FACTOR:array[0..5] of number=(1,1/2,3/4,15/8,105/16,945/32);
  var i:longint;
      factor: number;
  begin
    case (pre.f.e[c]+pre.g.e[c])div 2 of
      0: if pre.f.e[c]=pre.g.e[c] then result:=1 //beide 0
         else if pre.f.e[c]>=1 then
           result:=(pre.f.alpha*p1[c]+pre.g.alpha*p2[c])/pre.alphasum-p1[c]//ortho.
         else result:=(pre.f.alpha*p1[c]+pre.g.alpha*p2[c])/pre.alphasum-p2[c];
      1..5: begin
        result:=0;
        for i:=0 to (pre.f.e[c]+pre.g.e[c]) div 2 do begin
         factor:=binomSum(2*i,pre.f.e[c],pre.g.e[c],
                          (pre.f.alpha*p1[c]+pre.g.alpha*p2[c])/pre.alphasum-p1[c],
                          (pre.f.alpha*p1[c]+pre.g.alpha*p2[c])/pre.alphasum-p2[c]);
         Result:=result+factor*LOOK_UP_FACTOR[i]/intpower(pre.alphasum,i);
        end;
      end;
      else assert(false,'Orbitaltyp wird nicht unterstützt');
    end;
  end;
begin
  result:=pre.factor;
  result*=exp(-veclensqr(vecsub(p1,p2))*pre.alphamulprosum);
  result*=asym('x')*asym('y')*asym('z');
end;

//Berechnet das Integral der kinetischen Energie
function preCalculateKineticIntegral(const factor: number; const f,
  g: TFreeGTF): TGTFPreCalculatedKinetic;
var temp: TFreeGTF;
    c:char;
begin
  result.overlap[0]:=preCalculateOverlapIntegral(
      factor*g.alpha*(2*(g.e['x']+g.e['y']+g.e['z'])+3),f,g);
  temp:=g;
  for c:='x' to 'z' do begin
    temp.e[c]+=2;
    result.overlap[ord(c)-ord('x')+1]:=
      preCalculateOverlapIntegral(-factor*2*sqr(g.alpha),f,temp);
    temp.e[c]-=2;
  end;
end;

function kineticIntegral(const pre: TGTFPreCalculatedKinetic;
  const p1, p2: TVector3n): number;
var i:longint;
begin
  result:=0;
  for i:=0 to 3 do
    result+=overlapIntegral(pre.overlap[i],p1,p2);
end;

//Berechnet die Energie der Anziehung zwischen Kern und Orbital
function preCalculateNuclearAttractionIntegral(const factor:number;
                              const f,g:TFreeGTF):TGTFPreCalculatedNuclearAttraction;
var c:char;
    i,r,u:longint;
begin
  assert((f.alpha>0)and(g.alpha>0),'Falsche Alphawerte');
  Result.factor:=factor*2*pi/(f.alpha+g.alpha);
  result.f:=f;
  result.g:=g;
  result.coeffsum:=f.e['x']+g.e['x']+f.e['y']+g.e['y']+f.e['z']+g.e['z'];
  for c:='x' to 'z' do begin
    assert(f.e[c]+g.e[c]<=2,
           'Orbitaltyp wird nicht unterstützt'); //geht nur mit s,d,p
    for i:=0 to f.e[c]+g.e[c] do
      for r:=0 to i div 2 do
        for u:=0 to (i - 2*r) div 2 do begin
          if odd(i+u) then result.a_cache[c,i,r,u]:=-1
          else result.a_cache[c,i,r,u]:=1;
          result.a_cache[c,i,r,u]*=intpower(0.25/(f.alpha+g.alpha),r+u)*
             LOOK_UP_FAK[i]/(LOOK_UP_FAK[r]*LOOK_UP_FAK[u]*LOOK_UP_FAK[i-2*r-2*u]);
        end;
  end;
end;

function nuclearAttractionIntegral(const pre: TGTFPreCalculatedNuclearAttraction;
                                   const p1, p2, kern: TVector3n): number;
var a_cache:array['x'..'z',0..4,0..2,0..2] of number;
    f_cache:array[0..12] of number;
    c:char;
    tau: number;
    i,r,u,j,s,v,k,t,w: longint;
    P: TVector3n;
begin
  P:=midpoint(pre.f.alpha,pre.g.alpha,p1,p2);
  for c:='x' to 'z' do
    for i:=0 to pre.f.e[c]+pre.g.e[c] do
      for r:=0 to i div 2 do
        for u:=0 to (i - 2*r) div 2 do begin
          a_cache[c,i,r,u]:=pre.a_cache[c,i,r,u]*
                            binomSum(i,pre.f.e[c],pre.g.e[c],P[c]-p1[c],P[c]-p2[c]);
          if (i-2*(r+u)<>0) or (P[c]-kern[c]<>0) then //TODO: Wirklich = 1
            a_cache[c,i,r,u]*=power((P[c]-kern[c]),i-2*(r+u));//TODO: vielleicht abs?
        end;

  //Vorberechnung aller F_v(t) Terme in f_cache
  tau:=(pre.f.alpha+pre.g.alpha)*veclensqr(vecsub(kern,P));
  for i:=0 to pre.coeffsum do
    f_cache[i]:=f_erf_like(i,tau);

  //Zusammenfügen der Caches
  result:=0;
  for i:=0 to pre.f.e['x']+pre.g.e['x'] do
   for r:=0 to i div 2 do
    for u:=0 to (i - 2*r) div 2 do
     for j:=0 to pre.f.e['y']+pre.g.e['y'] do
      for s:=0 to j div 2 do
       for v:=0 to (j - 2*s) div 2 do
        for k:=0 to pre.f.e['z']+pre.g.e['z'] do
         for t:=0 to k div 2 do
          for w:=0 to (k - 2*t) div 2 do
           result+=a_cache['x',i,r,u]*a_cache['y',j,s,v]*
                   a_cache['z',k,t,w]*f_cache[i+j+k-2*(r+s+t)-u-v-w];

  //Vorfaktoren
  with pre do
    result*=exp(-(f.alpha*g.alpha/(f.alpha+g.alpha)) *
            veclensqr(vecSub(p1,p2)))*factor;
end;

//Berechnet die durch die Abstoßung zwischen Elektronen gegebene Energie
function preCalculateElectronRepulsionIntegral(const factor: number; const f,
  g: TFreeGTF): TGTFPreCalculatedElectronRepulsion;
var c:char;
    i:longint;
    i1,i2,r1,r2,u:longint;
begin
  for c:='x' to 'z' do
    assert((f.e[c]<=1)and(g.e[c]<=1),
           'Orbitaltyp nicht erlaubt');
  result.f:=f;
  result.g:=g;
  result.alphasum:=f.alpha+g.alpha;
  result.factor:=factor/result.alphasum;
  result.coeffsum:=f.e['x']+g.e['x']+
                   f.e['y']+g.e['y']+
                   f.e['z']+g.e['z'];

  Result.a_power[0]:=1;
  for i:=-1 downto low(Result.a_power) do
    Result.a_power[i]:=Result.a_power[i+1]/result.alphasum;

  for c:='x' to 'z' do
   for i1:=0 to f.e[c]+g.e[c] do
    for i2:=0 to high(result.b_cache[c,i1]) do
     for r1:=0 to i1 div 2 do
      for r2:=0 to i2 div 2 do
       for u:=0 to (i1+i2)div 2-r1-r2 do with result do begin
         //Ort unabhängig
         if odd(i1+u) then b_cache[c,i1,i2,r1,r2,u]:=-1 //Änderung i2->i1 WARUM??
         else b_cache[c,i1,i2,r1,r2,u]:=1;
         b_cache[c,i1,i2,r1,r2,u]*=
           a_power[r1-i1]

           *LOOK_UP_FAK[i1]*LOOK_UP_FAK[i2]
           *LOOK_UP_2POWER[2*(r1+r2-i1-i2)]
           *LOOK_UP_FAK[i1+i2-2*(r1+r2)]
           /(LOOK_UP_FAK[r1]*LOOK_UP_FAK[r2]*
             LOOK_UP_FAK[i1-2*r1]*LOOK_UP_FAK[i2-2*r2]*
             LOOK_UP_FAK[u]*LOOK_UP_FAK[i1+i2-2*(r1+r2+u)]);
       end;
end;

function electronRepulsionIntegral(const pre12,pre34:
         TGTFPreCalculatedElectronRepulsion;
  const p1,p2,p3,p4:TVector3n):number;
var P12,P34: TVector3n;
    c:char;
    //Schleifenvariablen
    i,i1,i2,r1,r2,u,j1,j2,s1,s2,v,k1,k2,t1,t2,w:longint;
    //Vorberechnete Werte
    pq_power,hf12_cache,hf34_cache:array[0..4]of number;
    b_cache: array['x'..'z',0..4,0..4,0..2,0..2,0..4] of number;
    f_cache: array[0..24] of number;
    delta_power: array[-4..0] of number;
    tau,erf_sqrt_tau,pq,delta:number;
begin
  P12:=midpoint(pre12.f.alpha,pre12.g.alpha,p1,p2);
  P34:=midpoint(pre34.f.alpha,pre34.g.alpha,p3,p4);
  tau:=pre12.alphasum*pre34.alphasum*veclensqr(vecsub(P12,P34))/
                                               (pre12.alphasum+pre34.alphasum);
  if pre12.coeffsum+pre34.coeffsum=0 then result:=f_erf_like(0,tau)
  else begin  //Mindestens ein p-Orbital
    delta:=0.25*(1/pre12.alphasum+1/pre34.alphasum);
    delta_power[0]:=1;
    for i:=-1 downto low(delta_power) do
      delta_power[i]:=delta_power[i+1]/delta;

    erf_sqrt_tau:=erf(sqrt(tau));
    for i:=0 to pre12.coeffsum+pre34.coeffsum do
      f_cache[i]:=f_erf_like(i,tau,erf_sqrt_tau);
    for c:='x' to 'z' do begin
      for i:=0 to pre12.f.e[c]+pre12.g.e[c] do
        hf12_cache[i]:=binomSum(i,pre12.f.e[c],pre12.g.e[c],
                                  P12[c]-p1[c],P12[c]-p2[c]);
      for i:=0 to pre34.f.e[c]+pre34.g.e[c] do
        hf34_cache[i]:=binomSum(i,pre34.f.e[c],pre34.g.e[c],
                                  P34[c]-p3[c],P34[c]-p4[c]);
      pq:=P12[c]-P34[c];
      pq_power[0]:=1;
      for i:=1 to high(pq_power) do
        pq_power[i]:=pq_power[i-1]*pq;

      //B cache füllen
      for i1:=0 to pre12.f.e[c]+pre12.g.e[c] do
       for i2:=0 to pre34.f.e[c]+pre34.g.e[c]do
        for r1:=0 to i1 div 2 do
         for r2:=0 to i2 div 2 do
          for u:=0 to (i1+i2)div 2-r1-r2 do
            b_cache[c,i1,i2,r1,r2,u]:=pre12.b_cache[c,i1,i2,r1,r2,u]*
                                      pre34.a_power[r2-i2]*
              hf12_cache[i1]*hf34_cache[i2]*pq_power[i1+i2-2*(r1+r2+u)]*
              delta_power[2*(r1+r2)-i1-i2+u];
    end;

    //Caches multiplizieren
    Result:=0;
    for i1:=0 to pre12.f.e['x']+pre12.g.e['x'] do
    for i2:=0 to pre34.f.e['x']+pre34.g.e['x']do
    for r1:=0 to i1 div 2 do
    for r2:=0 to i2 div 2 do
    for u:=0 to (i1+i2)div 2-r1-r2 do
      for j1:=0 to pre12.f.e['y']+pre12.g.e['y'] do
      for j2:=0 to pre34.f.e['y']+pre34.g.e['y']do
      for s1:=0 to j1 div 2 do
      for s2:=0 to j2 div 2 do
      for v:=0 to (j1+j2)div 2-s1-s2 do
        for k1:=0 to pre12.f.e['z']+pre12.g.e['z'] do
        for k2:=0 to pre34.f.e['z']+pre34.g.e['z']do
        for t1:=0 to k1 div 2 do
        for t2:=0 to k2 div 2 do
        for w:=0 to (k1+k2)div 2-t1-t2 do
          result+=b_cache['x',i1,i2,r1,r2,u]*
                  b_cache['y',j1,j2,s1,s2,v]*
                  b_cache['z',k1,k2,t1,t2,w]*
                  f_cache[i1+i2+j1+j2+k1+k2
                          -2*(r1+r2+s1+s2+t1+t2)
                          -u-v-w];
  end;
  result*=pre12.factor*pre34.factor*2*pi*pi*sqrt(pi/(pre12.alphasum+pre34.alphasum));
  result*=exp(-pre12.f.alpha*pre12.g.alpha*veclensqr(vecsub(p1,p2))/pre12.alphasum
              -pre34.f.alpha*pre34.g.alpha*veclensqr(vecsub(p3,p4))/pre34.alphasum)
end;

end.

