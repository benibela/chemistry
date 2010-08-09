unit geometryCalculator;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,extmath,atoms,molecules,energyCalculator;

type

  //Ermöglicht Fortschrittsanzeige bei Geometrieberechnung
  TEnergyOptimizationProgress=procedure(sender: tobject;energy:number;steps: longint;
                                            var cancel: boolean) of object;
  { TGeometryCalculator }

  TGeometryCalculator=class
  private
    currentMolecule: TMolecule;
  public
    procedure initialize(mol:TMolecule); virtual;
    procedure calculate;virtual;abstract;
  end;

  { TGeometryOptimizer }

  TGeometryOptimizer=class(TGeometryCalculator)
  private
    energyCalculator: TEnergyCalculator;
    //Bewegt das Molekül in eine feste Lage
    procedure normalizeMove;
  public
    progress:TEnergyOptimizationProgress;
    maxsteps:longint;
    constructor create(aenergyCalculator: TEnergyCalculator);
    procedure initialize(mol:TMolecule); override;
    procedure calculate; override;
  end;

  { TGeometryTrivial }

  TGeometryTrivial=class(TGeometryCalculator)
    procedure initialize(mol:TMolecule); override;
    procedure calculate; override;
  end;

implementation

uses math;
{ TGeometryOptimizer }

//Dreht das gesamte Molekül, so dass:
//   Atom 0 am Nullpunkt liegt
//   Atom 1 auf der x-Achse  (y=0,z=0)
//   Atom 2 auf der xy-Ebene (z=0)
procedure TGeometryOptimizer.normalizeMove;
var c:char;
    i:longint;
    phi,theta,r:number;
    rotMat:TMatrix;
begin
  with currentMolecule do begin
    //Bewege Atom 0 zu 0
    for c:='x' to 'z' do
      if atoms[0].p[c]<>0 then
        for i:=atomCount-1 downto 0 do
          atoms[i].p[c]-=atoms[0].p[c];
    //Rotiere Atom 1 auf x
    vec2sphere(atoms[1].p,r,phi,theta);
    SetLength(rotMat,3,3);
    matmul(createRotationMat3(phi-pi/2,'y'),createRotationMat3(-theta,'z'),rotMat);
    for i:=0 to atomCount-1 do
      atoms[i].p:=matVecTrans(rotMat,atoms[i].p);
    if atomCount>2 then begin
      //Rotiere Atom 2 auf xy
      with atoms[2] do
        rotMat:=createRotationMat3(-arctan2(p['z'],p['y']),'x');
      for i:=0 to atomCount-1 do
        atoms[i].p:=matVecTrans(rotMat,atoms[i].p);
    end;
  end;
end;

constructor TGeometryOptimizer.create(AEnergyCalculator: TEnergyCalculator);
begin
  energyCalculator:=AEnergyCalculator;
end;

procedure TGeometryOptimizer.initialize(mol: TMolecule);
begin
  inherited;
  energyCalculator.initialize(mol);
  maxsteps:=100;
end;

procedure TGeometryOptimizer.calculate;
const difstep=0.1; //Schritt beim Ableiten
      movstep=0.1; //Schritt beim Bewegen
var gradient: array of TVector3n; //Berechnete Ableitung
    c:char; i:longint;            //Schleifenvariablen
    energy,oldenergy:number;      //Aktuelle Energie und letzte
    gradDiff,old: number;         //Differenz der Ableitung
    cancel: boolean;              //Abbrechen
begin
 with currentMolecule do begin
   //Initialisieren
   energy:=9999999999;
   SetLength(gradient,atomCount);
   fillchar(gradient[0],sizeof(gradient[0])*length(gradient),0);
   cancel:=false;
   //Aktuelle Energie berechnen
   energy:=energyCalculator.calculate();
   //Erste Atome auf xy-Ebene bewegen
   normalizeMove;
   //Optimierungsverfahren
   repeat
    gradDiff:=0;
    //Ableitung numerisch für alle Koordinaten berechnen
    for i:=1 to atomCount-1 do
      for c:='x' to 'z' do begin
        if ord(c)-ord('x') >= i then break;
        old:=gradient[i][c];
        atoms[i].p[c]+=difstep;
        gradient[i][c]:=energyCalculator.calculate(i);
        atoms[i].p[c]-=2*difstep;
        gradient[i][c]-=energyCalculator.calculate(i);
        atoms[i].p[c]+=difstep;
        gradient[i][c]*=0.5*(movstep/difstep);
        //Differenz speichern
        gradDiff+=sqr(old-gradient[i][c]);
        //Große Sprünge vermeiden
        if gradient[i][c]>2 then gradient[i][c]:=2
        else if gradient[i][c]<-2 then gradient[i][c]:=-2;
      end;
    //Abbrechen bei kleiner Differenz
    if gradDiff<sqr(0.0001) then break;
    for i:=1 to atomCount-1 do
      for c:='x' to 'z' do
        atoms[i].p[c]-=gradient[i][c];
    //Energie berechnen
    oldenergy:=energy;
    energy:=energyCalculator.calculate();
    //Fortschrittsfunktion aufrufen
    if assigned(progress) then begin
      progress(self,energy,maxsteps,cancel);
      if cancel then break;
    end;
    dec(maxsteps)
  until (abs(oldenergy-energy)<1e-10)or(maxsteps<=0);
  geometryChanged();
  end;
end;

{ TGeometryTrivial }

procedure TGeometryTrivial.initialize(mol: TMolecule);
begin
  inherited;
end;

procedure TGeometryTrivial.calculate;
type TTetraeder=array[1..4] of TVector3n;
  //Berechnung einer Bindungslänge
  //nach http://www-lehre.inf.uos.de/~okrone/DIP/node13.html
  function bondLen(const t1,t2:longint;const n: longint):number;
  const BOND_SUB:array[1..3] of number=(0,0.397,0.642);
  var SchomakerStevenson: number;
  begin
    SchomakerStevenson:=0.1511781;
    if (t1>10) or (t2>10) then SchomakerStevenson:=0.038;
    result:=CovalentRadius[t1]+CovalentRadius[t2]+
            SchomakerStevenson*abs(EN[t1]-EN[t2])-
            BOND_SUB[n];
  end;
  //Tetraeder spiegeln um Eben mit Normale n
  function mirrorTetraeder(const tetra: TTetraeder; const n: TVector3n):
                                                                    TTetraeder;
  var i:longint;
  begin
    for i:=1 to 4 do
      result[i]:=mirrorVec(tetra[i],n);
  end;
  //Tetraeder um Achse axe um a Grad spiegeln
  function rotateTetraeder(const tetra: TTetraeder; const axe: TVector3n;
                           const a: number): TTetraeder;
  var i:longint;
      rm: TMatrix;
  begin
    rm:=createRotationMat3(a,axe);
    for i:=1 to 4 do
      result[i]:=matVecTrans(rm,tetra[i]);

  end;
  //Sortiert die Ecken so um, dass die n-Ecken nach (inklusive) e
  //hinten sind
  function hideEdge(const tetra: TTetraeder; const e,n:longint): TTetraeder;
  var i:longint;
  begin
    if e+n>4 then exit(tetra); //sind bereits hinten
    for i:=e to e+n-1 do
      result[i-e+5-n]:=tetra[i];
    for i:=1 to e-1 do
      result[i]:=tetra[i];
    for i:=e+n to 4 do
      result[i-n]:=tetra[i];
  end;

var arranged: array of boolean; //Hash bereits positionierter Atome
  //Atom a anordnen
  procedure arrange(a: tatom; tetra: TTetraeder);
  var i,j,e:longint;
      temp: TVector3n; //Richtung zum nächsten Kern
      b:^TBonds;
  begin
    e:=1;
    b:=@a.bond;
    for i:=1 to a.james do //Für alle Bindungen
      if not arranged[b^[i].a.index] then begin //wenn Atom nicht positioniert
        arranged[b^[i].a.index]:=true; //ist es das  gleich
        case b^[i].count of
          1:  begin //Ecke an Ecke
            //Atom zur nächsten Ecke
            b^[i].a.p:=vecadd(a.p,
                               vecscale(bondLen(a.typ,b^[i].a.typ,1),tetra[e]));
            //Rekursiv fortsetzen mit gespiegeltem und gedrehten Tetraeder
            arrange(b^[i].a,hideEdge(mirrorTetraeder(
                            rotateTetraeder(tetra,tetra[e],pi),tetra[e]),e,1));;
          end;
          2:  begin //Kante an Kante
            //Atom zwischen die nächsten beiden Ecken
            temp:=vecadd(tetra[e],tetra[e+1]);
            vecnormalize(temp);
            b^[i].a.p:=vecadd(a.p,vecscale(bondLen(a.typ,b^[i].a.typ,2),temp));
            //Tetraeder an temp spiegel
            arrange(b^[i].a,hideEdge(mirrorTetraeder(tetra,temp),e,2));;;
          end;
          3: begin //Seite an Seite
            //Atom zwischen drei Ecken
            temp:=vecadd(vecadd(tetra[e],tetra[e+1]),tetra[e+2]);
            vecnormalize(temp);
            //Tetraeder an temp spiegel
            b^[i].a.p:=vecadd(a.p,vecscale(bondLen(a.typ,b^[i].a.typ,3),temp));
            arrange(b^[i].a,hideEdge(mirrorTetraeder(tetra,temp),e,3));;;
          end;
        end;
        e+=b^[i].count;
      end;
  end;
var tetra: TTetraeder;
begin
 with currentMolecule do begin
  if atomCount=0 then exit;
 // if not closedShell then exit;
  setlength(arranged,atomCount);
  FillChar(arranged[0],length(arranged)*sizeof(arranged[0]),0);
  arranged[0]:=true;
  atoms[0].p:=nullvector3n;
  //normalisierter Ausgangstetraeder
  tetra[1]['x']:=1/sqrt(3);
  tetra[1]['y']:=1/sqrt(3);
  tetra[1]['z']:=-1/sqrt(3);

  tetra[2]['x']:=-1/sqrt(3);
  tetra[2]['y']:=-1/sqrt(3);
  tetra[2]['z']:=-1/sqrt(3);

  tetra[3]['x']:=-1/sqrt(3);
  tetra[3]['y']:=1/sqrt(3);
  tetra[3]['z']:=1/sqrt(3);

  tetra[4]['x']:=1/sqrt(3);
  tetra[4]['y']:=-1/sqrt(3);
  tetra[4]['z']:=1/sqrt(3);
  //Erstes Atom positionieren
  arrange(atoms[0],tetra);
  geometryChanged;
 end;
end;

{ TGeometryCalculator }

procedure TGeometryCalculator.initialize(mol: TMolecule);
begin
  currentMolecule:=mol;
end;

end.

