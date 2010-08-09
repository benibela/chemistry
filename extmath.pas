{Dieser Programmteil bieten erweiterte Mathematikfunktionen und Datenstrukturen,
 die in normalem Pascal nicht vorhanden sind}
unit extmath;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,math;
type
  {Strukturen für Vektoren und Matrizen}
  number = double;
  TVector3n=array['x'..'z'] of number;
  TVector3f=array['x'..'z'] of single;
  TVectorn=array of number;
  TVectorb=array of boolean;
  TVectori=array of longint;
  TMatrix=array of array of number;
  TSymMatrix=array of array of number; //a[i,j] existiert-> i<=j (obere Dreiecksmat.)
  
  {Vorberechnete Funktionen}
//Einfache Fakultät
const LOOK_UP_FAK:array[0..11] of longint=(1,1,2,6,24,120,720,5040,40320,
                                           362880,3628800,39916800);
//Doppelte Fakultät (2n-1)!!
const LOOK_UP_FAK2:array[-1..57] of number=
 (1,1,{1:}1,2,3,8,15,48,105,384,945,3840,10395,46080,135135,645120,2027025,10321920,
  34459425,185794560,654729075,3715891200,13749310575,81749606400,316234143225,
  1961990553600,7905853580625,51011754393600,213458046676875,1428329123020800,
  6190283353629375,42849873690624000,191898783962510625,1371195958099968000,
  6332659870762850625,46620662575398912000, 221643095476699771875,
  1678343852714360832000,8200794532637891559375, 63777066403145711616000,
  319830986772877770815625,2551082656125828464640000, 13113070457687988603440625,
107145471557284795514880000,563862029680583509947946875,4714400748520531002654720000,
  25373791335626257947657609375,216862434431944426122117120000,
  1192568192774434123539907640625,10409396852733332453861621760000,
  58435841445947272053455474390625,  520469842636666622693081088000000,
  2980227913743310874726229193921875,  27064431817106664380040216576000000,
  157952079428395476360490147277859375,1461479318123759876522171695104000000,
  8687364368561751199826958100282265625,81842841814930553085241614925824000000,
  495179769008019818390136611716089140625);
{Zweier Potenzen}
const LOOK_UP_2POWER:array[-13..61] of number=
  (1/8192,1/4096,1/2048,1/1024,1/512,1/256,1/128,1/64,1/32,1/16,1/8,1/4,1/2,
   1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,
   262144,524288,1048576,2097152,4194304,8388608,16777216,33554432, 67108864,
   134217728,268435456,536870912,1073741824, 2147483648,4294967296,8589934592,
   17179869184,34359738368,68719476736,137438953472,274877906944,549755813888,
   1099511627776,2199023255552,4398046511104,8796093022208,17592186044416,
   35184372088832,70368744177664,140737488355328,281474976710656,562949953421312,
   1125899906842624,2251799813685248,4503599627370496,9007199254740992,
   18014398509481984,36028797018963968,72057594037927936,144115188075855872,
   288230376151711744,576460752303423488,1152921504606846976,2305843009213693952
  );
//Binomialkoeffizienten
const LOOK_UP_PASCAL:array[0..11,0..11] of longint= //Binomialkoeffizienten
    ((                      1,           0,0,0,0,0,0,0,0,0,0,0),
     (                    1,  1,          0,0,0,0,0,0,0,0,0,0),
     (                  1,  2,  1,         0,0,0,0,0,0,0,0,0),
     (                1,  3,  3,  1,        0,0,0,0,0,0,0,0),
     (              1,  4,  6,  4,  1,       0,0,0,0,0,0,0),
     (            1,  5, 10, 10,  5,  1,      0,0,0,0,0,0),
     (          1,  6, 15, 20, 15,  6,  1,     0,0,0,0,0),
     (        1,  7, 21, 35, 35, 21,  7,  1,    0,0,0,0),
     (      1,  8, 28, 56, 70, 56, 28,  8,  1,   0,0,0),
     (    1,  9, 36, 84, 126,126, 84,36,  9,  1,  0,0 ),
     (  1, 10, 45,120,210,252,210,120, 45, 10,  1, 0),
     (1, 11, 55, 165,330,462,462,330,165,55, 11,   1));
{Null}
const epsilon16=1e-16;
const nullvector3n: TVector3n=(0,0,0);
      nullvector3ic:array['x'..'z'] of longint=(0,0,0);

//===Lineare Algebra===
//--Vektorrechnung--
//Elementare Vektorfunktionen (s. Formelsammlung)
function vec2f(const v: TVector3n): TVector3f;
procedure vecscale(var v: TVector3n; const scale: number);inline;overload;
function vecscale(const scale: number;const v: TVector3n): TVector3n;inline;overload;
function vecadd(const v1,v2: TVector3n): TVector3n;inline;
function vecsub(const v1,v2: TVector3n):TVector3n;inline;
function vecsub(const v1,v2: TVector3f):TVector3f;inline;
function veclensqr(const v: TVector3n):number;inline;overload;
function veclensqr(const v: TVector3f):single;inline;overload;
function vecdist(const v1,v2: TVector3n):number;inline;overload; //Vektordistanz
procedure vecnormalize(var v:TVector3n);                         //Normalisieren
function vecout(const v: TVectorn): string;
function mirrorVec(const v,n:TVector3n):TVector3n; //Vektorspiegelung an Ebene
//Umrechnen in Kugelkoordinaten
procedure vec2sphere(const v:TVector3n;var r,phi,theta: number);

//--Matrixrechnung---
procedure matinit(var M:TMatrix;const r,c:longint);
//Erzeugen von Rotationsmatrizen
function createRotationMat3(const alpha:number; const axe: char):TMatrix;
function createRotationMat3(const theta: number; const axe: TVector3n): TMatrix;
//Löschen einer Zeile
procedure materaseLine(var M:TMatrix; const l:longint);inline;
//Matrix mal Vektor
function matVecTrans(const M:TMatrix;const v: TVector3n):TVector3n;
//Untere Dreiecksmatrix mit oberer füllen
procedure matsym2full(var M: TMatrix);
//Matrixmultiplikation
procedure matmul(const A,B:TMatrix;var R: TMatrix);
procedure matmul_sym(const A,B:TSymMatrix;var R: TSymMatrix);
procedure matmul_withtransA(const A,B:TMatrix;var R: TMatrix);
//Berechnung von A^-0.5
procedure matinvsqrt_sym(var A:TSymMatrix;var T1,T2: TMatrix; var Tvb: TVectorb;
                         var Tvi: TVectori; var Tv: TVectorn; var R: TSymMatrix);
//Diagonalisieren einer Matrix
procedure diag_sym(var M: TSymMatrix; var changed: TVectorb;var ind: TVectori;
                   var eigenvalues:TVectorn; var eigenvectors: TMatrix);
function matout(const mat: TMatrix): string;

//---Analysis---
//base^exponent
function intpower00(base : float;const exponent : Integer) : float; //0^0 = 0
function intpower01(base : float;const exponent : Integer) : float; //0^0 = 1
//Summe der Binomialkoeffizienten
function binomSum(const j,l,m:longint; const a,b:number):number;
//Berechnet die Fehlerfunktion
function erf(const v: number;const delta: number=epsilon16): number;

procedure assertcomp(n1,n2:number;s:string;prec:number=0.000001);inline;
implementation
uses dialogs,windows;
//Wandelt einen Vektor in einen 3 float Vektor um
function vec2f(const v: TVector3n): TVector3f;
begin
  result['x']:=single(v['x']);
  result['y']:=single(v['y']);
  result['z']:=single(v['z']);
end;
//Elementare Vektorrechnung
procedure vecscale(var v: TVector3n; const scale: number); inline;
begin
  v['x']:=scale*v['x'];
  v['y']:=scale*v['y'];
  v['z']:=scale*v['z'];
end;

function vecscale(const scale: number; const v: TVector3n): TVector3n; inline;
begin
  Result:=v;
  vecscale(result,scale);
end;

function vecadd(const v1, v2: TVector3n): TVector3n; inline;
begin
  result['x']:=v1['x']+v2['x'];
  result['y']:=v1['y']+v2['y'];
  result['z']:=v1['z']+v2['z'];
end;

//Subtrahiert zwei Vektoren
function vecsub(const v1,v2: TVector3n):TVector3n;inline;
begin
  result['x']:=v1['x']-v2['x'];
  result['y']:=v1['y']-v2['y'];
  result['z']:=v1['z']-v2['z'];
end;

function vecsub(const v1, v2: TVector3f): TVector3f; inline;
begin
  result['x']:=v1['x']-v2['x'];
  result['y']:=v1['y']-v2['y'];
  result['z']:=v1['z']-v2['z'];
end;

//Berechnet das Quadrat der Länge des Vektors
function veclensqr(const v: TVector3n):number;inline;
begin
  result:=sqr(v['x'])+sqr(v['y'])+sqr(v['z']);
end;

function veclensqr(const v: TVector3f): single; inline;
begin
  result:=sqr(v['x'])+sqr(v['y'])+sqr(v['z']);
end;

function vecdist(const v1, v2: TVector3n): number; inline;
begin
  result:=sqrt(veclensqr(vecsub(v1,v2)));
end;

function vecout(const v: TVectorn): string;
var i:integer;
begin
  result:='';
  for i:=0 to high(v) do
    result:=result+format('%2.4g',[v[i]])+#9;
end;

//Vektor v an Ebene mit Normalenvektor n durch 0,0,0 spiegeln
function mirrorVec(const v, n: TVector3n): TVector3n;
var d:number;
begin
  //Der Abstand von v und v' zur Ebene ist gleich
  //v-d*n=v'+d*n
  //v'=v-2*d*n
  d:=v['x']*n['x']+v['y']*n['y']+v['z']*n['z'];
  result:=vecsub(v,vecscale(2*d,n));
end;

//Umrechnen in Kugelkoordinaten
procedure vec2sphere(const v: TVector3n; var r, phi, theta: number);
begin
  r:=sqrt(veclensqr(v));
  phi:=arccos(v['z']/r);
  theta:=arctan2(v['y'],v['x']);
end;

procedure vecnormalize(var v: TVector3n);
var len: number;
begin
  len:=sqrt(veclensqr(v));
  vecscale(v,1/len);
end;


//Wandelt die Matrix in ein für Menschen lesbares Format um
function matout(const mat: TMatrix): string;
const useStyle=2;
const STYLE:array[1..2,0..4] of string=(('','',#9,#13#10,''),     //Menschen
                                        ('.','[[',',','],[',']]'));//Derive
var i,j:longint;
begin
  if style[usestyle][0]<>'' then  DecimalSeparator:=style[usestyle][0][1];
  result:=STYLE[useStyle][1];
  for i:=0 to high(mat) do begin
    for j:=0 to high(mat[i]) do begin
      result:=result+format('%2.4g',[mat[i,j]]);
      if j<>high(mat[i]) then
        result+=STYLE[useStyle][2];
    end;
    if i<>high(mat) then
      result:=result+STYLE[useStyle][3];
  end;
  result:=result+STYLE[useStyle][4];
end;

//Berechnet base^exponent schnell und gibt 0 für 0^0 zurück
function intpower00(base : float;const exponent : Integer) : float;
begin
 if (base=0) and (exponent=0) then exit(0)
 else exit(intpower(base,exponent));
end;

//Berechnet base^exponent schnell und gibt 1 für 0^0 zurück
function intpower01(base : float;const exponent : Integer) : float;
begin
 if (base=0) and (exponent=0) then exit(1)
 else exit(intpower(base,exponent));
end;


//Berechnet den Faktor vor dem Summanden mit dem angegebenen Exponenten
//     l      m     l+m      l+m-1                    0  l  m
//(x+a) *(x+b)  = x     +  x       (l*a+m*b) + ... + x  a  b
//
//Also, binomSum(l+m,l.m) = 1, binomSum(l+m-1,l.m) = l*a+m,
//      ...                               binomSum(0,l.m) = a^l * b^m,
function binomSum(const j,l,m:longint; const a,b:number):number;
var i,f:longint;
    ac,bc:number;
begin
  if j>l+m then exit(0);
  if j=l+m then exit(1);
  if (a=0) and (b=0) then exit(0);
  if a=0 then
    if l+m-j>m then exit(0)
    else exit(LOOK_UP_PASCAL[m,l+m-j]* intpower(b,l+m-j));
  if b=0 then
    if l+m-j>l then exit(0)
    else exit(LOOK_UP_PASCAL[l,l+m-j]*intpower(a,l+m-j));


  result:=0;
  f:=max(0,j-m);
  ac:=intpower01(a,l-f);
  bc:=intpower01(b,m+f-j);
  for i:=f to min(j,l) do begin
    result:=result+LOOK_UP_PASCAL[l,i]*LOOK_UP_PASCAL[m,j-i]*ac*bc;
    ac/=a;
    bc*=b;
  end;
end;

//Berechnet die Fehlerfunktion
// v     -x^2
// S dx e
// 0
//nach http://mathworld.wolfram.com/Erf.html
function erf(const v: number;const delta: number): number;
var change,value,fak,v_power,v_sqr,n: number;
    sign: boolean;
begin
  if IsZero(v) then exit(0)
  else if SameValue(v,1) then exit(0.8427007929507194429276)
  else if v>3.8325063453686 then exit(1)
  else if v < 0 then exit(-erf(-v,delta))
  else begin
    value:=v;
    fak:=1;
    v_sqr:=v*v;
    v_power:=v*v_sqr;
    sign:=true;
    n:=1;
    repeat
      change:=v_power/(fak*(2*n+1));
      if sign then value:=value-change
      else value:=value+change;
      v_power:=v_power*v_sqr;
      sign:=not sign;
      n:=n+1;
      fak:=fak*n;
    until abs(change)<delta;
    result:=value*2/sqrt(pi);
  end;
end;

procedure assertcomp(n1, n2: number; s: string;prec:number=0.000001);inline;
begin
  assert(samevalue(n1,n2,prec),s+#13#10+floattostr(n1)+'<>'+floattostr(n2));
end;

//Einheitsmatrix erstellen
procedure matinit(var M: TMatrix; const r, c: longint);
var i:longint;
begin
  SetLength(M,r,c);
  for i:=0 to r-1 do materaseLine(m,i);
  for i:=0 to min(r-1,c-1) do
    m[i,i]:=1;
end;

//Triviale Rotation um Achsenkoordinaten
//nach http://www.grundstudium.info/animation/node14.php
function createRotationMat3(const alpha:number; const axe: char):TMatrix;
var s,c:number;
begin
  s:=sin(alpha);
  c:=cos(alpha);
  matinit(result,3,3);
  case axe of
    'x':begin
      result[1,1]:=c;
      result[1,2]:=-s;
      result[2,1]:=s;
      result[2,2]:=c;
    end;
    'y':begin
      result[0,0]:=c;
      result[0,2]:=-s;
      result[2,0]:=s;
      result[2,2]:=c;
    end;
    'z':begin
      result[0,0]:=c;
      result[0,1]:=-s;
      result[1,0]:=s;
      result[1,1]:=c;
    end;
  end;
end;

//Rotation um eine beliebige Achse
//Nach http://www.grundstudium.info/animation/node14.php
function createRotationMat3(const theta: number; const axe: TVector3n): TMatrix;
var {Rx,Ry,Rz,Ry_inv,Rx_inv: TMatrix;
    alpha,beta,d:number;}
    si,co:number;
begin
  if (axe['y']=0)and(axe['z']=0) then
    exit(createRotationMat3(theta,'x'));
  si:=sin(theta);
  co:=cos(theta);
  setlength(result,3,3);
 { result[0, 0] := ((1-co) * sqr(axe['x'])) + co;
  result[1, 0] := ((1-co) * axe['x'] * axe['y']) - (axe['z'] * si);
  result[2, 0] := ((1-co) * axe['z'] * axe['x']) + (axe['y'] * si);

  result[0, 1] := ((1-co) * axe['x'] * axe['y']) + (axe['z'] * si);
  result[1, 1] := ((1-co) * sqr(axe['y'])) + co;
  result[2, 1] := ((1-co) * axe['y'] * axe['z']) - (axe['x'] * si);

  result[0, 2] := ((1-co) * axe['z'] * axe['x']) - (axe['y'] * si);
  result[1, 2] := ((1-co) * axe['y'] * axe['z']) + (axe['x'] * si);
  result[2, 2] := ((1-co) * sqr(axe['z'])) + co;
}
  result[0, 0] := ((1-co) * sqr(axe['x'])) + co;
  result[0, 1] := ((1-co) * axe['x'] * axe['y']) - (axe['z'] * si);
  result[0, 2] := ((1-co) * axe['z'] * axe['x']) + (axe['y'] * si);

  result[1, 0] := ((1-co) * axe['x'] * axe['y']) + (axe['z'] * si);
  result[1, 1] := ((1-co) * sqr(axe['y'])) + co;
  result[1, 2] := ((1-co) * axe['y'] * axe['z']) - (axe['x'] * si);

  result[2, 0] := ((1-co) * axe['z'] * axe['x']) - (axe['y'] * si);
  result[2, 1] := ((1-co) * axe['y'] * axe['z']) + (axe['x'] * si);
  result[2, 2] := ((1-co) * sqr(axe['z'])) + co;
end;

//Elementare Matrixrechnung
procedure materaseLine(var M: TMatrix; const l: longint);inline;
begin
  fillchar(M[l][0],length(M[l])*sizeof(M[l][0]),0);
end;

function matVecTrans(const M: TMatrix; const v: TVector3n):TVector3n;
var i,j:longint;
    c:char;
begin
  for c:='x' to 'z' do begin
    result[c]:=0;
    for j:=0 to 2 do
      result[c]+=M[ord(c)-ord('x')][j]*v[chr(ord('x')+j)];
  end;
end;

procedure matsym2full(var M: TMatrix);
var i,j: longint;
begin
  assert(length(M)>0,'Matrix leer');
  assert(length(M[0])=length(M),'nicht symmetrisch');
  for i:=0 to high(M) do
    for j:=i+1 to high(M) do
      M[j,i]:=M[i,j];
end;

procedure matmul(const A, B: TMatrix; var R: TMatrix);
var i,j,k:longint;
begin
  assert(length(A[0])=length(B),'Ungültige Größe für Multiplikation');
  assert(length(R)=length(B[0]),'Ungültige Größe für Multiplikation');
  assert(length(R[0])=length(B[0]),'Ungültige Größe für Multiplikation');

  for i:=0 to high(A) do
    for j:=0 to high(B[0]) do begin
      r[i,j]:=0;
      for k:=0 to high(B) do
        r[i,j]+=A[i,k]*B[k,j];
    end;
end;

procedure matmul_sym(const A, B: TMatrix; var R: TMatrix);
var i,j,k:longint;
begin
  assert(length(A)=length(B),'Ungültige Größe für Multiplikation');
  assert(length(A)=length(A[0]),'Ungültige Größe für Multiplikation');
  assert(length(A)=length(B[0]),'Ungültige Größe für Multiplikation');
  assert(length(A)=length(R),'Rückgabe Matrix falsch');
  assert(length(A)=length(R[0]),'Rückgabe Matrix falsch');

  for i:=0 to high(A) do
    for j:=i to high(A) do begin
      //Berechnung ist trivial, aber trickreich
      //Für jedes Koordinatenpaar [i,j] muss gelten i<=j
      r[i,j]:=0;
      for k:=0 to i do
        r[i,j]+=A[k,i]*B[k,j];
      for k:=i+1 to j do
        r[i,j]+=A[i,k]*B[k,j];
      for k:=j+1 to high(A) do
        r[i,j]+=A[i,k]*B[j,k];
    end;
end;

//Multiplisiert A mit der Transponierten von B
procedure matmul_withtransA(const A, B: TMatrix; var R: TMatrix);
var i,j,k:longint;
begin
  assert(length(A)>0,'Matrix leer');
  assert(length(A[0])=length(B[0]),'Ungültige Größe für Multiplikation');
  assert(length(R)=length(B[0]),'Ungültige Größe für Multiplikation');
  assert(length(R[0])=length(B[0]),'Ungültige Größe für Multiplikation');
  for i:=0 to high(A) do
    for j:=0 to high(B) do begin
      r[i,j]:=0;
      for k:=0 to high(A[0]) do
        r[i,j]+=A[k,i]*B[k,j];
    end;
end;

//Berechnet R=A^-0.5
//T1,T2,Tv sind temporär, müssen die gleiche Größe haben
//Nach Cook, Ab Initio Valence Calculations in Chemistry
procedure matinvsqrt_sym(var A:TSymMatrix;var T1,T2: TMatrix; var Tvb: TVectorb;
                         var Tvi: TVectori; var Tv: TVectorn; var R: TSymMatrix);
var ev: TVectorn;
    i,j:longint;
begin
  assert(length(A)>0,'Matrix leer');
  assert(length(A[0])=length(A),'Ungültige Größe für M^-0.5');
  assert(length(R)=length(A),'Ungültige Größe für M^-0.5: r');
  assert(length(R[0])=length(A),'Ungültige Größe für M^-0.5: r');
  assert(length(T1)=length(A),'Ungültige Größe für M^-0.5: t1');
  assert(length(T1[0])=length(A),'Ungültige Größe für M^-0.5: t1');
  assert(length(T2)=length(A),'Ungültige Größe für M^-0.5: t2');
  assert(length(T2[0])=length(A),'Ungültige Größe für M^-0.5: t2');

  diag_sym(A,Tvb,Tvi,Tv,T1);
  //ev:=ev^-1/2
  for i:=0 to high(Tv) do
    Tv[i]:=1/sqrt(Tv[i]);
  //T2:=T1*Tv;
  for i:=0 to high(T1) do
    for j:=0 to high(T1) do
      T2[i,j]:=T1[i,j]*Tv[i];
  //R:=temp2*temp1;
  matmul_withtransA(t1,t2,r);
end;

{//Diagonalisiert die Matrix M
//Die obere Seite von M wird zerstört
//Die Eigenvektoren sind in der Matrix Zeilenvektoren
//von http://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm
//version 20:23, 22 October 2006
//version-url:
http://en.wikipedia.org/w/index.php?title=Jacobi_eigenvalue_algorithm&oldid=83064407
}
procedure diag_sym(var M: TSymMatrix; var changed: TVectorb;var ind: TVectori;
                   var eigenvalues:TVectorn; var eigenvectors: TMatrix);
  //index of largest element in row k
  function maxind(k: longint): longint;
  var i:longint;
  begin
    result:=k+1;
    for i := k+2 to high(M[k]) do
      if abs(M[k,i]) > abs(M[k,result]) then result:=i;
  end;
  
var state:longint;
  procedure update(k:longint; t: number); // update ek and its status
  var y:number;
  begin
    y := eigenvalues[k];
    eigenvalues[k] := y+t;
    if changed[k] and (IsZero(t)) then begin
      changed[k] := false;
      state := state-1;
    end else if (not changed[k]) and (not IsZero(t)) then begin
      changed[k] := true;
      state := state+1
    end;
  end;

  procedure rotate(k,l,i,j : longint;s,c:number); // perform rotation of Sij, Skl
  var Mkl,Mij:number;
  begin
   {|   |    |     ||   |
    |Skl|    |c  -s||Skl|
    |   | := |     ||   |
    |Sij|    |s   c||Sij|
    |   |    |     ||   |}
    Mkl:=M[k,l];
    Mij:=M[i,j];
    M[k,l]:=c*Mkl-s*Mij;
    M[i,j]:=s*Mkl+c*Mij;
  end;

var
  i, j, k, l: longint;
  s, c, t, p, y, Eki, Eli: number;
begin
  assert(length(M)>0,'Matrix leer');
  assert(length(M)=length(M[0]),'Nur quadratische Matrizen werden diagonalisiert');
  assert(length(eigenvectors)=length(M),'Größe falsch');
  assert(length(eigenvectors[0])=length(M),'Größe falsch');
  assert(length(eigenvalues)=length(M),'Größe falsch');
  assert(length(ind)=length(M),'Größe falsch');
  assert(length(changed)=length(M),'Größe falsch');
  for i:=0 to high(eigenvectors) do begin
    fillchar(eigenvectors[i,0],length(M)*sizeof(eigenvectors[i,0]),0);
    eigenvectors[i,i]:=1;
  end;
  state := length(M);
  for k := 0 to high(M) do begin
    ind[k] := maxind(k);
    eigenvalues[k] :=M[k,k];
    changed[k] := true;
  end;
  //Berechnung
  while state>0 do begin// next rotation
    l := 0; // find index (k,l) of pivot p
    for k := 1 to high(M)-1 do
      if abs(M[k,ind[k]])> abs(M[l,ind[l]]) then l := k;
    k := l;l := ind[k];p := M[k,l];
    if iszero(p) then break;
    // calculate c = cos phi, s = sin phi
    y := (eigenvalues[l]-eigenvalues[k])/2;
    t := abs(y)+sqrt(p*p+y*y);
    s := sqrt(p*p+t*t);
    c := t/s; s := p/s; t := p*p/t;
    if y<0 then begin
      s := -s;
      t := -t;
    end;
    M[k,l] := 0.0; update(k,-t); update(l,t);
    // rotate rows and columns k and l
    for i := 0 to k-1 do rotate(i,k,i,l,s,c);
    for i := k+1 to l-1 do rotate(k,i,i,l,s,c);
    for i := l+1 to high(M) do rotate(k,i,l,i,s,c);
    // rotate eigenvectors
    for i := 0 to high(eigenvectors) do begin
      {?   ?    ?     ??   ?
      ?Eki?    ?c  -s??Eki?
      ?   ? := ?     ??   ?
      ?Eli?    ?s   c??Eli?
      ?   ?    ?     ??   ?}
      eki:=eigenvectors[k,i];
      eli:=eigenvectors[l,i];
      eigenvectors[k,i]:=c*eki-s*eli;
      eigenvectors[l,i]:=c*eli+s*eki;
    end;
    // rows k, l have changed, update rows indk, indl
    ind[k] := maxind(k); ind[l] := maxind(l);
  end;
end;
end.

