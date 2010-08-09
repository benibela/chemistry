unit moscf;

{$mode objfpc}{$H+}
{$define performtests}

interface

uses
  Classes, SysUtils,molecules;

type
  TGTF=record
    e:array['x'..'z'] of longint; //Exponenten
    n:TVector;//KO-System Nullpunkt
    alpha: number;  //Exponent a
    pre: number;   //Vorfaktor b
    //Funktion der Form
    {        xe     ye    ze   -a*r²
      b*(x-x0) (y-y0) (z-z0)  e
    }
  end;

const LOOK_UP_FAK:array[0..11] of longint=(1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800);

function calculateEnergy(self: TMolecule; reuseOldData: longint):string;

implementation
uses math, extmath;

//Berechnet die Koordinate c vom "Mittelpunkt" zweier GTF-Orbitale
function P(const f,g:TGTF;const c:char):number;inline;
begin
  result:=(f.alpha*f.n[c]+g.alpha*g.n[c])/(f.alpha+g.alpha);
end;

//Berechnet die Differenz der Koordinate c der GTF-Orbitale f vom
//"Mittelpunkt" der beiden GTF-Orbitale
function P_f_difference(const c:char;const f,g:TGTF):number;inline;
begin
  result:=P(g,f,c)-f.n[c];  //TODO: Überprüfen, ob hier abs stehen sollte
end;

//Berechnet den Mittelpunkt zweier GTF-Orbitale als Vektor
function P(const f,g:TGTF):TVector;inline;
var fa_plus_ga:number;
begin
  fa_plus_ga:=f.alpha+g.alpha;
  result['x']:=(f.alpha*f.n['x']+g.alpha*g.n['x'])/fa_plus_ga;
  result['y']:=(f.alpha*f.n['y']+g.alpha*g.n['y'])/fa_plus_ga;
  result['z']:=(f.alpha*f.n['z']+g.alpha*g.n['z'])/fa_plus_ga;
end;

//nach http://mathworld.wolfram.com/GaussianIntegral.html
//Berechnet das Integral
//  1     2v  -ts²
// S  ds s   e
//  0

{Notizen:
   x = y/sqrt(t)
   dx = dy/sqrt(t)
   y^2 = t x^2

   Mit v = 0:
      1     -tx^2     sqrt(t)       -y^2
     S  dx e       = S  dy/sqrt(t) e     = 1/sqrt(t) * sqrt(pi) / 2 erf(sqrt(t))
      0               0

      = sqrt(pi/t) / 2 erf(sqrt(t))


   Selbe Substitution für v <> 0:
      1     2v   -tx^2    sqrt(t)                    2v    -y^2
     S  dx x  e           = S  dy/sqrt(t) (y/sqrt(t))    e
      0                     0
              -1  sqrt(t)  2v        2v  -y^2
      =sqrt(t)      S  dy y  / sqrt(t)  e
                    0
              -2v-1 sqrt(t)   2v   -y^2
      =sqrt(t)        S  dy  y    e
                      0
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
             sqrt(pi)   d^s     -1/2      1/2
    F  (t) = -------- -------  t     erf(t   )
     2s         2       dt^s


     s=0: sqrt(pi)/(2sqrt(t))  erf(sqrt(t))

     s=1:
                 sqrt(pi)  d      -1/2      1/2
        F  (t) = -------- -----  t     erf(t   )
         2         2       dt

                    sqrt(pi) 1  -3/2 1   -1/2  d erf
               = -  -------- - t     -  t      ----- (sqrt(t))
                      2      2       2           dt

                    sqrt(pi) 1  -2      2      -t
               = -  -------- - t      -----   e
                      2      4       sqrt(pi)

                    1  -2     -t
               = -  - t      e
                    4

                    1        -3        -t
        F  (t) = -  - (-2) t     (-1) e
         4          4
                    1   -3  -t
               = -  - t    e
                    2
}

function f_erf_like(const v:longint;const t:number):number;
var temp: number;
    i:longint;
begin
  if isZero(t) then exit(1/(2*v+1))
  else begin
    if v>=1 then begin
      temp:=1/t;
      result:=0;
      for i:=v-1 downto 0 do begin
        result-=LOOK_UP_FAK2[2*v-1]/(LOOK_UP_2POWER[v-i]*LOOK_UP_FAK2[2*i+1])*temp;
        temp/=t;
      end;
      result:=exp(-t)*result;
    end else result:=0;
     result+=LOOK_UP_FAK2[2*v-1]*sqrt(pi)*erf(sqrt(t),1e-16)/((LOOK_UP_2POWER[v+1])*power(t,v+0.5));
  end;
end;

function help_function(const j:longint;const f,g:TGTF; const c:char):number;
// inline;

const LOOK_UP_PASCAL:array[0..6,0..6] of longint= //Binomialkoeffizienten
    ((            1,         0,0,0,0,0,0),
     (          1,  1,        0,0,0,0,0),
     (        1,  2,  1,       0,0,0,0),
     (      1,  3,  3,  1,      0,0,0),
     (    1,  4,  6,  4,  1,     0,0),
     (  1,  5, 10, 10,  5,  1    ,0),
     (1,  6, 15, 20, 15,  6,  1));


var i:integer;
begin
  result:=0;
  for i:=max(0,j-g.e[c]) to min(j,f.e[c]) do
    result:=result+LOOK_UP_PASCAL[f.e[c],i]*LOOK_UP_PASCAL[g.e[c],j-i]*
                   intpower01(P_f_difference(c,f,g),f.e[c]-i)*
                   intpower01(P_f_difference(c,g,f),g.e[c]+i-j);
end;

//Berechnet das Überlappungsintegral, aber ignoriert dabei die Vorfaktoren
//<f,g> = S f g dV
function overlapIntegral_nofactor(const f,g: TGTF): number;
  function asym(const c: char):number; //inline;
  //const LOOK_UP_FAK2:array[0..10] of longint=(1,1,2,3,8,15,48, 105, 384, 945, 3840);
  const LOOK_UP_FACTOR:array[0..5] of number=(1,1/2,3/4,15/8,105/16,945/32); //(2i-1)!!*2^-i
  var i:longint;
      factor: number;
  begin
    assert((f.e[c]+g.e[c])div 2<=high(LOOK_UP_FACTOR),'Orbitaltyp wird nicht unterstützt');
    result:=0;
    for i:=0 to (f.e[c]+g.e[c]) div 2 do begin
      factor:=help_function(2*i,f,g,c);
      {factor:=0;
      for k:=max(0,2*i-g.e[c]) to min(2*i,f.e[c]) do
        factor:=factor+LOOK_UP_PASCAL[f.e[c],k]*LOOK_UP_PASCAL[g.e[c],2*i-k]*
                       intpower(f.n[c],f.e[c]-k)*intpower(g.n[c],g.e[c]+k-2*i);}

      Result:=result+factor*LOOK_UP_FACTOR[i]/intpower(f.alpha+g.alpha,i);
    end;
  end;
begin
  result:=power(      pi/
            (f.alpha+g.alpha), 1.5);
  result*=exp(
            -(sqr(f.n['x']-g.n['x'])+sqr(f.n['y']-g.n['y'])+sqr(f.n['z']-g.n['z']))*
                (f.alpha*g.alpha/
                (f.alpha+g.alpha)));
  result*=asym('x')*asym('y')*asym('z');
end;

//Berechnet das Überlappungsintegral der GTF-Orbitale f und g, nachdem der
//Koeffizient c von g um v geändert ist.
//g wird dabei (nur) temporär geändert (=> keine Threadsicherheit mehr)
function overlapIntegral_modified_nofactor(const f:TGTF; var g:TGTF; const c:char; const v:longint):number;inline;
begin
  inc(g.e[c],v);
  result:=overlapIntegral_nofactor(f,g);
  dec(g.e[c],v);
end;

//Berechnet das vollständige Überlappungsintegral mitsamt
//Vorfaktoren
function overlapIntegral(const f,g: TGTF): number;inline;
begin
  result:=overlapIntegral_nofactor(f,g)*f.pre*g.pre;
end;

//Berechnet das Integral der kinetischen Energie
function kineticEnergyIntegral(const f:TGTF; var g:TGTF):number;
var c:char;
begin
  result:=g.alpha*(
      2*(g.e['x']+g.e['y']+g.e['z'])+3)*overlapIntegral_nofactor(f,g)
      - 2*sqr(g.alpha)* (
            overlapIntegral_modified_nofactor(f,g,'x',2)
            +overlapIntegral_modified_nofactor(f,g,'y',2)
            +overlapIntegral_modified_nofactor(f,g,'z',2)
      );
  //Vollständigshalber, bei s und p Orbitalen ist dies überflüssig
  for c:='x' to 'z' do
    if g.e[c]>=2 then
      result:=Result+g.e[c]*(g.e[c]-1)*0.5*
                overlapIntegral_modified_nofactor(f,g,c,-2);
                          ;
  //Normalisierungsfaktoren
  result:=f.pre*g.pre*result;
end;

function nuclearAttractionIntegral(const f,g:TGTF; const kern: TVector):number;
  //Hilfsfunktion
  function A(c:char;i,r,u:longint):number;
  begin
    if odd(i+u) then result:=-1
    else result:=1;
    assert(i<=high(LOOK_UP_FAK),'i zu groß (Nuklear-Attraction-Integral:A)');
    assert(r<=high(LOOK_UP_FAK),'r zu groß (Nuklear-Attraction-Integral:A)');
    assert(u<=high(LOOK_UP_FAK),'u zu groß (Nuklear-Attraction-Integral:A)');
    assert(i-2*r-2*u<=high(LOOK_UP_FAK),'i-2*r-2*u zu groß (Nuklear-Attraction-Integral:A)');
    assert(i-2*r-2*u>=low(LOOK_UP_FAK),'i-2*r-2*u zu klein (Nuklear-Attraction-Integral:A)');

    if (i-2*(r+u)<>0) or (P(f,g,c)-kern[c]<>0) then //TODO: Wirklich = 1
      result*=power((P(f,g,c)-kern[c]),i-2*(r+u));//TODO: vielleicht abs?
    result*=help_function(i,f,g,c)*LOOK_UP_FAK[i];
    result*=power((f.alpha+g.alpha)/4,r+u);
    result/=LOOK_UP_FAK[r]*LOOK_UP_FAK[u]*LOOK_UP_FAK[i-2*r-2*u];
  end;

var tau,temp: number;
    c:char;
    i,r,u,j,s,v,k,t,w: longint;
    f_cache:array[0..12] of number;
    a_cache:array['x'..'z',0..4,0..2,0..2] of number;
begin
  result:=f.pre*g.pre;
  result*=pi/(f.alpha+g.alpha);
  result*=exp(-(f.alpha*g.alpha/
                (f.alpha+g.alpha)) * veclensqr(vecSub(f.n,g.n)));
  //Vorberechnung aller F_v(t) Terme in f_cache
  tau:=(f.alpha+g.alpha)*veclensqr(vecsub(kern,P(f,g)));
  assert(f.e['x']+g.e['x']+f.e['y']+g.e['y']+f.e['z']+g.e['z']<=high(f_cache),
         'Orbitaltyp wird nicht unterstützt');
  for i:=0 to f.e['x']+g.e['x']+f.e['y']+g.e['y']+f.e['z']+g.e['z'] do
    f_cache[i]:=f_erf_like(i,tau);
  //Vorberechnung aller A Terme in a_cache
  for c:='x' to 'z' do begin
    assert(f.e[c]+g.e[c]<=high(a_cache[c]),
           'Orbitaltyp wird nicht unterstützt'); //geht mit s,d,p
    for i:=0 to f.e[c]+g.e[c] do
      for r:=0 to i div 2 do
        for u:=0 to (i - 2*r) div 2 do
          a_cache[c,i,r,u]:=A(c,i,r,u);
  end;
  //Zusammenfügen der Ergebnisse
  temp:=0;
  for i:=0 to f.e['x']+g.e['x'] do
   for r:=0 to i div 2 do
    for u:=0 to (i - 2*r) div 2 do
     for j:=0 to f.e['y']+g.e['y'] do
      for s:=0 to j div 2 do
       for v:=0 to (j - 2*s) div 2 do
        for k:=0 to f.e['z']+g.e['z'] do
         for t:=0 to k div 2 do
          for w:=0 to (k - 2*t) div 2 do
           temp+=a_cache['x',i,r,u]*a_cache['y',j,s,v]*a_cache['z',k,t,w]*
                 f_cache[i+j+k-2*(r+s+t)-u-v-w];
  result*=temp;
end;

function electronRepulsionIntegral(const f1,f2,f3,f4:TGTF):number;
var tau: number;
    alpha12,alpha34,delta: number;
  function B(c:char;i1,i2,r1,r2,u:longint):number;
  begin
    if odd(i2+u) then result:=-1
    else result:=1;
    result*=help_function(i1,f1,f2,c)*help_function(i2,f3,f4,c);
    result*=LOOK_UP_FAK[i1]*LOOK_UP_FAK[i2];
    result*=intpower(alpha12,r1-i1)*intpower(alpha34,r2-i2);
    result*=intpower((2*delta),r1+r2);
    result/=LOOK_UP_FAK[r1]*LOOK_UP_FAK[r2]*
            LOOK_UP_FAK[i1-2*r1]*LOOK_UP_FAK[i2-2*r2]*intpower(4*delta,i1+i2);
    result*=LOOK_UP_FAK[i1+i2-2*(r1+r2)];
    result*=intpower01(abs(P(f1,f2,c)-P(f3,f4,c)),i1+i2+2*(r1+r2+u)); //TODO: check abs
    result*=intpower01(delta,u);
    result/=LOOK_UP_FAK[u]*LOOK_UP_FAK[i1+i2-2*(r1+r2+u)];
  end;

 var i:longint;
     f_cache: array[0..24] of number;
     c:char;
     i1,i2,r1,r2,u,j1,j2,s1,s2,v,k1,k2,t1,t2,w:longint;
     b_cache: array['x'..'z',0..4,0..4,0..2,0..2,0..2] of number;
begin
  alpha12:=f1.alpha+f2.alpha;
  alpha34:=f3.alpha+f4.alpha;
  tau:=alpha12*alpha34*veclensqr(vecsub(P(f1,f2),P(f3,f4)))/(alpha12+alpha34);
  delta:=0.25*(1/alpha12+1/alpha34);
  assert(f1.e['x']+f2.e['x']+f3.e['x']+f4.e['x']+
         f1.e['y']+f2.e['y']+f3.e['y']+f4.e['y']+
         f1.e['z']+f2.e['z']+f3.e['z']+f4.e['z']<=high(f_cache),
           'Orbitaltyp wird nicht unterstützt');
  for i:=0 to f1.e['x']+f2.e['x']+f3.e['x']+f4.e['x']+
              f1.e['y']+f2.e['y']+f3.e['y']+f4.e['y']+
              f1.e['z']+f2.e['z']+f3.e['z']+f4.e['z'] do
    f_cache[i]:=f_erf_like(i,tau);
  for c:='x' to 'z' do begin
    assert(f1.e[c]+f2.e[c]<=high(b_cache[c]),
           'Orbitaltyp wird nicht unterstützt'); //geht mit s,d,p
    assert(f3.e[c]+f4.e[c]<=high(b_cache[c]),
           'Orbitaltyp wird nicht unterstützt'); //geht mit s,d,p
    for i1:=0 to f1.e[c]+f2.e[c] do
     for i2 := 0 to f3.e[c] + f4.e[c]do
      for r1:=0 to i1 div 2 do
       for r2:=0 to i2 div 2 do
        for u:=0 to (i1+i2)div 2-r1-r2 do
         b_cache[c,i1,i2,r1,r2,u]:=B(c,i1,i2,r1,r2,u);
  end;

  Result:=0;
  for i1:=0 to f1.e['x']+f2.e['x'] do
  for i2 := 0 to f3.e['x'] + f4.e['x']do
  for r1:=0 to i1 div 2 do
  for r2:=0 to i2 div 2 do
  for u:=0 to (i1+i2)div 2-r1-r2 do
    for j1:=0 to f1.e['y']+f2.e['y'] do
    for j2:=0 to f3.e['y'] + f4.e['y']do
    for s1:=0 to j1 div 2 do
    for s2:=0 to j2 div 2 do
    for v:=0 to (j1+j2)div 2-s1-s2 do
      for k1:=0 to f1.e['z']+f2.e['z'] do
      for k2 := 0 to f3.e['z'] + f4.e['z']do
      for t1:=0 to k1 div 2 do
      for t2:=0 to k2 div 2 do
      for w:=0 to (k1+k2)div 2-t1-t2 do
        result+=b_cache['x',i1,i2,r1,r2,u]*
                b_cache['y',j1,j2,s1,s2,v]*
                b_cache['z',k1,k2,t1,t2,w]*
                f_cache[i1+i2+j1+j2+k1+k2
                        -2*(r1+r2+s1+s2+t1+t2)
                        -u-v-w];


  result*=f1.alpha*f2.alpha*f3.alpha*f4.alpha;
  result*=2*pi*pi/(alpha12*alpha34);
  result*=sqrt(pi/(alpha12+alpha34));
  result*=exp(-f1.alpha*f2.alpha*veclensqr(vecsub(f1.n,f2.n))/alpha12-
              -f3.alpha*f4.alpha*veclensqr(vecsub(f3.n,f4.n))/alpha34);
end;

function calculateEnergy(self: TMolecule; reuseOldData: longint):string;
type TAO=array[1..3] of TGTF; //Ein AO wird durch LC dreier GTF gebildet(STO-3G)
var orbitals:array of TAO;    //Alle vorhandenen AO-Orbitale
  ap: longint;

procedure addgtf(atom: pointer; l,m,n: longint; pre, alpha: number);
var gtf: ^TGTF;
begin
  inc(ap);
  if ap>3 then begin
    ap:=1;
    setlength(orbitals,length(orbitals)+1);
  end;
  if tatom(atom)._energyTemp=-1 then tatom(atom)._energyTemp:=high(orbitals);
  gtf:=@orbitals[high(orbitals)][ap];
  gtf^.alpha:=alpha;
  gtf^.e['x']:=l;
  gtf^.e['y']:=m;
  gtf^.e['z']:=n;
  gtf^.n:=tatom(atom).p;
  gtf^.pre:=pre;
  //Normalisieren
  gtf^.pre*=sqrt((intpower(4,l+m+n)*sqrt(8)*power(alpha,l+m+n+1.5))/
                 (LOOK_UP_FAK2[l]*LOOK_UP_FAK2[m]*LOOK_UP_FAK2[n]*sqrt(pi*pi*pi)));
end;

var c: array of array of number; //Koeffizienten i-te MO = Summe (c * j-te AO) T
  elektrons: longint; //Elektronenzahl
  i,j,k,l,r,s,t,u:longint;
  overlap: array of array of number;                 //S
{    kineticEnergy:array of array of number;           //\_____
  nuclearAttraction: array of array of number;       //}
  oneElectronEnergy: array of array of number;       //H, kinetische+nukleare
  repulsion: array of array of array of array of number;
                                         //      +
  denseMatrix: array of array of number; //R = TT    ?
  fockMatrix: array of array of number; //Hf = oneElectronEnergy + "repulsion"
procedure addAOsToMO(const l:longint; const atom: tatom);inline;
begin
  case atom.typ of
    AT_H: c[l][atom._energyTemp]:=1; //H Atome geben 1s in die Bindung
    else begin
      //Probiere sp3-Hybridisierung
      c[l][atom._energyTemp+1]:=1;
      c[l][atom._energyTemp+2]:=1;
      c[l][atom._energyTemp+3]:=1;
      c[l][atom._energyTemp+4]:=1;
    end;
  end;
end;

begin
//=========================Generierung der GTF AOs============================
fillchar(orbitals,sizeof(orbitals),0);
SetLength(orbitals,0);
//SetLength(c,0,0);
ap:=3;
elektrons:=0;
//STO Parameter nach JCP38,2686:
//H: 1 nach blauem buch
//C: 5.6727 1.6083 1.5697
//N: 6.6651 1.9237 1.9170
//O: 7.6579 2.2458 2.2266
for i:=0 to atoms.count-1 do begin
  TAtom(atoms[i])._energyTemp:=-1;
  case TAtom(atoms[i]).typ of
    AT_H: begin
            elektrons+=1;
            addgtf(atoms[i],0,0,0,0.154321,2.22766);
            addgtf(atoms[i],0,0,0,0.535328,0.40577);
            addgtf(atoms[i],0,0,0,0.444635,0.1098175);
          end;
    AT_C: begin
            elektrons+=6;
            //1s
            addgtf(atoms[i],0,0,0,0.154321,2.22766*sqr(5.6727));
            addgtf(atoms[i],0,0,0,0.535328,0.40577*sqr(5.6727));
            addgtf(atoms[i],0,0,0,0.444635,0.1098175*sqr(5.6727));
            //2s
            addgtf(atoms[i],0,0,0,-0.0599447,2.58158*sqr(1.6083));
            addgtf(atoms[i],0,0,0,0.596039,0.156762*sqr(1.6083));
            addgtf(atoms[i],0,0,0,0.458179,0.0601833*sqr(1.6083));
            //2p
            addgtf(atoms[i],1,0,0,0.162395,0.919238*sqr(1.5697));
            addgtf(atoms[i],1,0,0,0.566171,0.235919*sqr(1.5697));
            addgtf(atoms[i],1,0,0,0.422307,0.0800981*sqr(1.5697));
            addgtf(atoms[i],0,1,0,0.162395,0.919238*sqr(1.5697));
            addgtf(atoms[i],0,1,0,0.566171,0.235919*sqr(1.5697));
            addgtf(atoms[i],0,1,0,0.422307,0.0800981*sqr(1.5697));
            addgtf(atoms[i],0,0,1,0.162395,0.919238*sqr(1.5697));
            addgtf(atoms[i],0,0,1,0.566171,0.235919*sqr(1.5697));
            addgtf(atoms[i],0,0,1,0.422307,0.0800981*sqr(1.5697));
          end;
    AT_N: begin
            elektrons+=7;
            //1s
            addgtf(atoms[i],0,0,0,0.154321,2.22766*sqr(6.6651));
            addgtf(atoms[i],0,0,0,0.535328,0.40577*sqr(6.6651));
            addgtf(atoms[i],0,0,0,0.444635,0.1098175*sqr(6.6651));
            //2s
            addgtf(atoms[i],0,0,0,-0.0599447,2.58158*sqr(1.9237));
            addgtf(atoms[i],0,0,0,0.596039,0.156762*sqr(1.9237));
            addgtf(atoms[i],0,0,0,0.458179,0.0601833*sqr(1.9237));
            //2p
            addgtf(atoms[i],1,0,0,0.162395,0.919238*sqr(1.9170));
            addgtf(atoms[i],1,0,0,0.566171,0.235919*sqr(1.9170));
            addgtf(atoms[i],1,0,0,0.422307,0.0800981*sqr(1.9170));
            addgtf(atoms[i],0,1,0,0.162395,0.919238*sqr(1.9170));
            addgtf(atoms[i],0,1,0,0.566171,0.235919*sqr(1.9170));
            addgtf(atoms[i],0,1,0,0.422307,0.0800981*sqr(1.9170));
            addgtf(atoms[i],0,0,1,0.162395,0.919238*sqr(1.9170));
            addgtf(atoms[i],0,0,1,0.566171,0.235919*sqr(1.9170));
            addgtf(atoms[i],0,0,1,0.422307,0.0800981*sqr(1.9170));
          end;
    AT_O: begin
            elektrons+=8;
            //1s
            addgtf(atoms[i],0,0,0,0.154321,2.22766*sqr(7.6579));
            addgtf(atoms[i],0,0,0,0.535328,0.40577*sqr(7.6579));
            addgtf(atoms[i],0,0,0,0.444635,0.1098175*sqr(7.6579));
            //2s
            addgtf(atoms[i],0,0,0,-0.0599447,2.58158*sqr(2.2458));
            addgtf(atoms[i],0,0,0,0.596039,0.156762*sqr(2.2458));
            addgtf(atoms[i],0,0,0,0.458179,0.0601833*sqr(2.2458));
            //2p
            addgtf(atoms[i],1,0,0,0.162395,0.919238*sqr(2.2266));
            addgtf(atoms[i],1,0,0,0.566171,0.235919*sqr(2.2266));
            addgtf(atoms[i],1,0,0,0.422307,0.0800981*sqr(2.2266));
            addgtf(atoms[i],0,1,0,0.162395,0.919238*sqr(2.2266));
            addgtf(atoms[i],0,1,0,0.566171,0.235919*sqr(2.2266));
            addgtf(atoms[i],0,1,0,0.422307,0.0800981*sqr(2.2266));
            addgtf(atoms[i],0,0,1,0.162395,0.919238*sqr(2.2266));
            addgtf(atoms[i],0,0,1,0.566171,0.235919*sqr(2.2266));
            addgtf(atoms[i],0,0,1,0.422307,0.0800981*sqr(2.2266));
          end;
  end;
end;




//=========================Berechnen der AO Integrale=========================
//Überlappungsintegrale
SetLength(overlap,length(orbitals),length(orbitals));
for i:=0 to high(overlap) do
  for j:=0 to high(overlap[i]) do begin
    overlap[i,j]:=0;
    for r:=1 to 3 do for s:=1 to 3 do
      overlap[i,j]+=overlapIntegral(orbitals[i,r],orbitals[j,s]);
  end;
//exit(matout(overlap));
//Elektronenkinetikintegrale
//Kernanziehungsintegrale
SetLength(oneElectronEnergy,length(orbitals),length(orbitals));
for i:=0 to high(oneElectronEnergy) do
  for j:=0 to high(oneElectronEnergy[i]) do begin
    oneElectronEnergy[i,j]:=0;
    for r:=1 to 3 do for s:=1 to 3 do
      oneElectronEnergy[i,j]+=kineticEnergyIntegral(orbitals[i,r],orbitals[j,s]);
    for k:=0 to atoms.count-1 do
      for r:=1 to 3 do for s:=1 to 3 do
        oneElectronEnergy[i,j]+=nuclearAttractionIntegral(orbitals[i,r],orbitals[j,s],
                                                          tatom(atoms[k]).p);
  end;
    exit(matout(oneElectronEnergy));
{  SetLength(kineticEnergy,length(orbitals),length(orbitals));
for i:=0 to high(kineticEnergy) do
  for j:=0 to high(kineticEnergy[i]) do begin
    kineticEnergy[i,j]:=0;
    for r:=1 to 3 do for s:=1 to 3 do
      kineticEnergy[i,j]+=kineticEnergyIntegral(orbitals[i,r],orbitals[j,s]);
  end;
SetLength(nuclearAttraction,atoms.count,length(orbitals),length(orbitals));
//Kernanziehungsintegrale
for i:=0 to atoms.count-1 do
  for j:=0 to high(kineticEnergy) do
    for k:=0 to high(kineticEnergy[i]) do begin
      nuclearAttraction[i,j,k]:=0;;
      for r:=1 to 3 do for s:=1 to 3 do
        nuclearAttraction[i,j,k]+=nuclearAttractionIntegral(orbitals[j,r],orbitals[k,s]);
    end;}
//Elektronenwechselwirkungsintegrale
fillchar(repulsion,sizeof(repulsion),0);//Größe auf 0 setzen (sonst => crash)
//setlength(repulsion,0,0,0,0);
SetLength(repulsion,length(orbitals),length(orbitals),length(orbitals),length(orbitals));
for i:=0 to high(repulsion) do
 for j:=0 to high(repulsion[i]) do
  for k:=0 to high(repulsion[i,j]) do
   for l:=0 to high(repulsion[i,j,k]) do begin
     repulsion[i,j,k,l]:=0;
     for r:=1 to 3 do for s:=1 to 3 do for t:=1 to 3 do for u:=1 to 3 do
       repulsion[i,j,k,l]+=electronRepulsionIntegral(orbitals[i,r],orbitals[j,s],
                                                     orbitals[k,t],orbitals[l,u]);
   end;




//=================Raten von guten Linearkombinationen========================
SetLength(c,elektrons div 2,length(orbitals));
l:=0;
//Ein MO für jede Bindung
for i:=0 to atoms.count - 1 do
  for j:=1 to TAtom(atoms[i]).bindungen do
    if TAtom(atoms[i]).bindung[j].a.index>TAtom(atoms[i]).index then begin
      fillchar(c[l][0],length(c[l])*sizeof(c[l][0]),0);
      addAOsToMO(l,TAtom(atoms[i]));
      addAOsToMO(l,TAtom(atoms[i]).bindung[j].a);
      l+=1;
    end;
//Übriggebliebene MO besetzten: 1s für C,N,O und noch 1 MO für N und 2 MOs für O
for i:=0 to atoms.count - 1 do begin
  if TAtom(atoms[i]).typ=AT_H then continue;
  fillchar(c[l][0],length(c[l])*sizeof(c[l][0]),0);
  c[l][TAtom(atoms[i])._energyTemp]:=1;
  l+=1;
  case TAtom(atoms[i]).typ of
    AT_N: begin
      //Probiere 2s
      fillchar(c[l][0],length(c[l])*sizeof(c[l][0]),0);
      c[l][TAtom(atoms[i])._energyTemp+1]:=1;
      l+=1;
    end;
    AT_O: begin
      //Probiere 2s
      fillchar(c[l][0],length(c[l])*sizeof(c[l][0]),0);
      c[l][TAtom(atoms[i])._energyTemp+1]:=1;
      l+=1;
      //und 2p-Verschmelzung
      fillchar(c[l][0],length(c[l])*sizeof(c[l][0]),0);
      c[l][TAtom(atoms[i])._energyTemp+2]:=1;
      c[l][TAtom(atoms[i])._energyTemp+3]:=1;
      c[l][TAtom(atoms[i])._energyTemp+4]:=1;
      l+=1;
    end;
  end;
end;
assert(2*l=elektrons,'Geschätze MO-Anzahl ungültig: '+inttostr(l));



//==========================Optimierung der MOs===============================
SetLength(denseMatrix,length(orbitals),Length(orbitals));
SetLength(fockMatrix,length(orbitals),Length(orbitals));

//Berechnung der Dichtematrix
for i:=0 to high(orbitals) do
  for j:=0 to high(orbitals) do begin
    denseMatrix[i,j]:=0;
    for k:=0 to high(c) do
      denseMatrix[i,j]+=c[k,i]*c[k,j];
    denseMatrix[i,j]*=2;
  end;

//Berechnung der Fockmatrix
for i:=0 to high(orbitals) do
  for j:=0 to high(orbitals) do begin
    fockMatrix[i,j]:=oneElectronEnergy[i,j];
    for k:=0 to high(orbitals) do
      for l:=0 to high(orbitals) do
        fockMatrix[i,j]+=denseMatrix[k,l]*(repulsion[i,j,k,l]-0.5*repulsion[i,l,k,j]);
  end;

//Berechnung der Eigenvektoren der Fockmatrix

Result:=matout(overlap);
end;
//Überprüft die Korrektheit der verwendeten Routinen
{$ifdef performtests}
initialization
  assert(SameValue(f_erf_like(0,0),1,1e-10),
         'Interner Test f_erf_like(0,0) fehlgeschlagen');
  assert(SameValue(f_erf_like(1,0),1/3,1e-10),
         'Interner Test f_erf_like(1,0) fehlgeschlagen');
  assert(SameValue(f_erf_like(2,0),1/5,1e-10),
         'Interner Test f_erf_like(2,0) fehlgeschlagen');
  assert(SameValue(f_erf_like(30,0),0.01639344262,1e-11),
         'Interner Test f_erf_like(30,0) fehlgeschlagen');

  assert(SameValue(f_erf_like(0,2),0.5981440066,1e-10),
         'Interner Test f_erf_like(0,2) fehlgeschlagen');
  assert(SameValue(f_erf_like(1,2),0.1157021808,1e-10),
         'Interner Test f_erf_like(1,2) fehlgeschlagen');
  assert(SameValue(f_erf_like(2,2),0.05294281483,1e-10),
         'Interner Test f_erf_like(2,2) fehlgeschlagen');
  assert(SameValue(f_erf_like(20,2),0.003637739895,1e-6),
         'Interner Test f_erf_like(20,2) fehlgeschlagen');

  assert(SameValue(f_erf_like(0,-3),4.222211992,1e-10),
         'Interner Test f_erf_like(0,-3) fehlgeschlagen');
  assert(SameValue(f_erf_like(1,-3),2.643887488,1e-10),
         'Interner Test f_erf_like(1,-3) fehlgeschlagen');
  assert(SameValue(f_erf_like(2,-3),2.025645743,1e-10),
         'Interner Test f_erf_like(2,-3) fehlgeschlagen');
  assert(SameValue(f_erf_like(10,-3),0.7558224637,1e-6),
         'Interner Test f_erf_like(10,-3) fehlgeschlagen');
{$endif}
end.

