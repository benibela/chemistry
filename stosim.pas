{Dieser Programmteil enthält Routinen zur Berechnung von Slater Type Orbitals
 durch Approximation mit GTFs (STO3G-Basis).
 Die Orbitale werden dabei nicht einmal umgeformt, sondern erst bei Berechnung.
 Zusammen mit der teilweisen Vorberechnung dieser Integralen liefert dies
 einen kleinen Geschwindigkeitsvorteil.
}
unit stosim;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,math,extmath,atoms,gtf;

type
  //Slatertypfunktion
  TSTO=record
    atom: TAtom;  //Am Atom
    typ: longint; //(1: 1s; 2: 2s; 3: 2px; 4: 2py; 5: 2pz);
  end;

const QN: array[1..5,1..3] of longint= //Quantenzahlen von Orbitaltyp
      ((1,0,0),(2,0,0),(2,1,-1), (2,1,0), (2,1,-1));

//Standardmäßig verfügbare Orbitale je Atom (Alle bis zu 2p, außer bei H, da 1s)
const OrbitalCount:array[1..8] of longint=(1,0,0,0,0,5,5,5);

//Liefert eine STO zurück
function createSTO(const atom: TAtom; typ: longint): TSTO;
//Berechnet das vollständige Überlappungsintegral
function overlapIntegral(const f,g: TSTO): number;
//Berechnet das Integral der kinetischen Energie
function kineticEnergyIntegral(const f,g:TSTO):number;
//Berechnet die Energie der Anziehung zwischen Kern und Orbital
function nuclearAttractionIntegral(const f,g:TSTO; const atom: TAtom):number;
//Berechnet die durch die Abstoßung zwischen Elektronen gegebene Energie
function electronRepulsionIntegral(const f1,f2,f3,f4:TSTO):number;

var initializated:boolean=false;
//Berechnet die möglichen Integrale teilweise
procedure initUnit;
implementation
const
  //Slater Type Orbital Parameter (nach JCP38,2686)
  STO_Params: array[1..8] of array [1..3] of number=(
    //1s     2s       2p
    (1, nan, nan),           //H
    (nan, nan, nan),         //He
    (nan, nan, nan),         //Li
    (3.6848, 0.9560, nan),   //Be
    (nan, nan, nan),         //B
    (5.6727, 1.6083, 1.5679),//C
    (6.6651, 1.9237, 1.9170),//N
    (7.6579, 2.2458, 2.2266) //O
  );
  //STO zu GTO Konvertierungsfaktoren (nach Cook)
  STO3G: array[1..3,1..3,1..2] of number=
  (((0.154321,2.22766),(0.535328,0.40577),(0.444635,0.1098175)),    //1s
   ((-0.0599447,2.58158),(0.596039,0.156762),(0.458179,0.0601833)), //2s
   ((0.162395,0.919238),(0.566171,0.235919),(0.422307,0.0800981)))  //2p
  ;

var
  //Berechnete Parameter für die GTF Zerlegung
  STO3G_GTFPreFactors: array[1..8,1..3,1..3] of number;
  STO3G_GTF:array[1..8,1..5,1..3] of TFreeGTF;
  //Vorberechnete Integrale
  STO3G_PC_Overlap: array[1..8,1..5,1..3,1..8,1..5,1..3] of TGTFPreCalculatedOverlap;
  STO3G_PC_Kinetic: array[1..8,1..5,1..3,1..8,1..5,1..3] of TGTFPreCalculatedKinetic;
  STO3G_PC_Nuclear: array[1..8,1..5,1..3,1..8,1..5,1..3] of
                                                  TGTFPreCalculatedNuclearAttraction;
  STO3G_PC_Repulsion: array[1..8,1..5,1..3,1..8,1..5,1..3] of
                                                  TGTFPreCalculatedElectronRepulsion;
  //Integrale am gleichen Zentrum sind Positionsunabhängig
  STO3G_PC_SimpleOverlap,STO3G_PC_SimpleKinetic,
  STO3G_PC_SimpleNuclear:array[1..8,1..5,1..5] of number;
  //..._s.r.[i,j,k,l]: i>=j, k>=l!
  STO3G_PC_SimpleRepulsion:array[1..8,1..5,1..5,1..5,1..5] of number;

function createSTO(const atom: TAtom; typ: longint): TSTO;
begin
  result.atom:=atom;
  result.typ:=typ;
end;

//Berechnet das vollständige Überlappungsintegral
//wenn beide Orbitale am selben Atom sind, wird das vorberechnete Ergebnis geliefert
//ansonsten die GTF Funktionen aufgerufen
function overlapIntegral(const f, g: TSTO): number;
var i,j:longint;
begin
  if f.atom=g.atom then exit(STO3G_PC_SimpleOverlap[f.atom.typ,f.typ,g.typ]);
  result:=0;
  for i:=1 to 3 do
    for j:=1 to 3 do
      result+=gtf.overlapIntegral(STO3G_PC_Overlap[f.atom.typ,f.typ,i,
                                              g.atom.typ,g.typ,j],f.atom.p,g.atom.p);
end;

//Berechnet das Integral der kinetischen Energie
//wenn beide Orbitale am selben Atom sind, wird das vorberechnete Ergebnis geliefert
//ansonsten die GTF Funktionen aufgerufen
function kineticEnergyIntegral(const f, g: TSTO): number;
var i,j:longint;
begin
  if f.atom=g.atom then exit(STO3G_PC_SimpleKinetic[f.atom.typ,f.typ,g.typ]);
  result:=0;
  for i:=1 to 3 do
    for j:=1 to 3 do
      result+=gtf.kineticIntegral(STO3G_PC_Kinetic[f.atom.typ,f.typ,i,
                              g.atom.typ,g.typ,j],f.atom.p,g.atom.p);
end;

//Berechnet die Energie der Anziehung zwischen Kern und Orbital
//wenn alle am Atom sind, wird das vorberechnete Ergebnis geliefert
//ansonsten die GTF Funktionen aufgerufen
function nuclearAttractionIntegral(const f, g: TSTO; const atom: TAtom): number;
var i,j:longint;
begin
  if (f.atom=g.atom) and (atom=f.atom) then
    exit(STO3G_PC_SimpleNuclear[f.atom.typ,f.typ,g.typ]);
  result:=0;
  for i:=1 to 3 do
    for j:=1 to 3 do
      result+=gtf.nuclearAttractionIntegral(
        STO3G_PC_Nuclear[f.atom.typ,f.typ,i,g.atom.typ,g.typ,j],
        f.atom.p,g.atom.p,atom.p);
  result*=atom.typ;
end;

//Berechnet die durch die Abstoßung zwischen Elektronen gegebene Energie
//wenn alle Orbitale am selben Atom sind, wird das vorberechnete Ergebnis geliefert
//ansonsten die GTF Funktionen aufgerufen
function electronRepulsionIntegral(const f1, f2, f3, f4: TSTO): number;
var i,j,k,l:longint;
begin
  if (f1.atom=f2.atom) and (f2.atom=f3.atom) and (f3.atom=f4.atom) then
    exit(STO3G_PC_SimpleRepulsion[f1.atom.typ,max(f1.typ,f2.typ),min(f1.typ,f2.typ),
                                             max(f3.typ,f4.typ),min(f3.typ,f4.typ)]);

  result:=0;
  for i:=1 to 3 do
   for j:=1 to 3 do
    for k:=1 to 3 do
     for l:=1 to 3 do
      result+=gtf.electronRepulsionIntegral(
        STO3G_PC_Repulsion[f1.atom.typ,f1.typ,i,f2.atom.typ,f2.typ,j],
        STO3G_PC_Repulsion[f3.atom.typ,f3.typ,k,f4.atom.typ,f4.typ,l],
        f1.atom.p,f2.atom.p,f3.atom.p,f4.atom.p);
end;

//Berechnet die möglichen Integrale teilweise
procedure initUnit;
const OX:array[1..5]of longint=(1,2,3,3,3); //Orbitalzuordnung
var i,j,k,l,m: longint;
    i1,i2,o1,o2,o3,o4,g1,g2:longint;
    f,g:TFreeGTF;
begin
  if initializated then exit;
  initializated:=true;
  //Berechnung der GTF Vorfaktorparameter
  for i:=1 to 8 do   //Atomart
   for j:=1 to 3 do  //Orbitalart(1s,2s,2p)
    for k:=1 to 3 do begin//Index der GTF
      if isnan(STO_Params[i,j]) then continue;
      STO3G_GTFPreFactors[i,j,k]:=STO3G[j,k,1]*
                                 normalizeFactorGTF(j div 3,0,0,
                                 STO3G[j,k,2]*sqr(STO_Params[i,j]));
    end;
  //Berechnet freie GTF Funktionen
  for i:=1 to 8 do
   for j:=1 to OrbitalCount[i] do
    for k:=1 to 3 do begin
      STO3G_GTF[i,j,k].alpha:=STO3G[OX[j],k,2]*sqr(STO_Params[i,OX[j]]);
      STO3G_GTF[i,j,k].e:=nullvector3ic;
      if j=3 then STO3G_GTF[i,j,k].e['x']:=1;
      if j=4 then STO3G_GTF[i,j,k].e['y']:=1;
      if j=5 then STO3G_GTF[i,j,k].e['z']:=1;
    end;
  //teilweise Vorberechnung der möglichen Integrale
  for i1:=1 to 8 do   //Atomart
   for o1:=1 to OrbitalCount[i1] do  //Orbitalart(1s,2s,2px,2py,2pz)
    for g1:=1 to 3 do //GTF-Index
     for i2:=1 to 8 do
      for o2:=1 to OrbitalCount[i2] do
       for g2:=1 to 3 do begin
         STO3G_PC_Overlap[i1,o1,g1,i2,o2,g2]:=preCalculateOverlapIntegral(
            STO3G_GTFPreFactors[i1,OX[o1],g1]*STO3G_GTFPreFactors[i2,OX[o2],g2],
            STO3G_GTF[i1,o1,g1],STO3G_GTF[i2,o2,g2]);
         STO3G_PC_Kinetic[i1,o1,g1,i2,o2,g2]:=preCalculateKineticIntegral(
            STO3G_GTFPreFactors[i1,OX[o1],g1]*STO3G_GTFPreFactors[i2,OX[o2],g2],
            STO3G_GTF[i1,o1,g1],STO3G_GTF[i2,o2,g2]);
         STO3G_PC_Nuclear[i1,o1,g1,i2,o2,g2]:=preCalculateNuclearAttractionIntegral(
            STO3G_GTFPreFactors[i1,OX[o1],g1]*STO3G_GTFPreFactors[i2,OX[o2],g2],
            STO3G_GTF[i1,o1,g1],STO3G_GTF[i2,o2,g2]);
        STO3G_PC_Repulsion[i1,o1,g1,i2,o2,g2]:=preCalculateElectronRepulsionIntegral(
            STO3G_GTFPreFactors[i1,OX[o1],g1]*STO3G_GTFPreFactors[i2,OX[o2],g2],
            STO3G_GTF[i1,o1,g1],STO3G_GTF[i2,o2,g2]);
       end;
  //vollständige Vorberechnungen aller Möglichkeiten an einem Atom
  for i:=1 to 8 do
   for o1:=1 to OrbitalCount[i] do
    for o2:=1 to OrbitalCount[i] do begin
      STO3G_PC_SimpleOverlap[i,o1,o2]:=0;
      for j:=1 to 3 do for k:=1 to 3 do
        STO3G_PC_SimpleOverlap[i,o1,o2]+=gtf.overlapIntegral(
          STO3G_PC_Overlap[i,o1,j,i,o2,k],nullvector3n,nullvector3n);

      STO3G_PC_SimpleKinetic[i,o1,o2]:=0;
      for j:=1 to 3 do for k:=1 to 3 do
        STO3G_PC_SimpleKinetic[i,o1,o2]+=gtf.kineticIntegral(
          STO3G_PC_Kinetic[i,o1,j,i,o2,k],nullvector3n,nullvector3n);

      STO3G_PC_SimpleNuclear[i,o1,o2]:=0;
      for j:=1 to 3 do for k:=1 to 3 do
        STO3G_PC_SimpleNuclear[i,o1,o2]+=gtf.nuclearAttractionIntegral(
          STO3G_PC_Nuclear[i,o1,j,i,o2,k],nullvector3n,nullvector3n,nullvector3n);
      STO3G_PC_SimpleNuclear[i,o1,o2]*=i;
      
      for o3:=1 to o1 do
        for o4:=1 to o2 do begin
          STO3G_PC_SimpleRepulsion[i,o1,o3,o2,o4]:=0;
          for j:=1 to 3 do for k:=1 to 3 do
           for l:=1 to 3 do for m:=1 to 3 do
            STO3G_PC_SimpleRepulsion[i,o1,o3,o2,o4]+=gtf.electronRepulsionIntegral(
              STO3G_PC_Repulsion[i,o1,j,i,o3,k],STO3G_PC_Repulsion[i,o2,l,i,o4,m],
              nullvector3n,nullvector3n,nullvector3n,nullvector3n);
        end;
        
    end;
end;

initialization
  initUnit;
end.

