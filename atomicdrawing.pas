{Dieser Programmteil benutzt oglDrawing um ein bestimmtes Molekül über OpenGl
 3-dimensional auszugeben.}
unit atomicDrawing;

{$mode objfpc}{$H+}
interface

uses
  Classes, SysUtils,ogldrawing,atoms,molecules,Graphics,extmath;

type
{ TAtomRenderer }

//Zeichenklasse
TAtomRenderer=class
  private
    ogl: TOGLRenderer;
    atoms: array[1..8] of TOGLGLUObject; //Speichert für jedes Atom eine Kugel
    link: TOGLGLUObject; //Speichert einen Zylinder als Bindungsstab
    usesModel: longint; //Verwendetes Modell (0: Kalotten, 1: Kugel/Stab)
  public
    colorScheme: longint;
    constructor create(oglRenderer: TOGLRenderer);
    procedure draw(mol: TMolecule);
    destructor destroy;override;
    //Legt das Modell (0: Kalotten, 1: Kugel/Stab) und
    //das Farbschema (0: normal, 1: Drucken) fest
    procedure useModel(const m: longint;const c:longint=-1);
    property currentModel: longint read usesModel write useModel;
end;
      //Farben für Atome
const atomColors: array[0..1,-1..8] of TColor=
      ((clBlack,clSilver,clSilver,0,0,0,0,clGray,clBlue,clRed),
      (clWhite,clGray,clSilver,0,0,0,0,clGray,clBlue,clRed));
      //van der Waals-Radius (Radius der Kalottenmodellkugeln) (für Kalottenmodell)
      vanderWaalsRadius: array[1..8] of number=
                      (2.268,0,  0,0,0,3.21,2.93,2.89); //in a0}
      //Radius einer Bindung (für KS.modell)
      bondWidth: number=0.2;
      //Radius eines Kerns (für KS.modell)
      kernRadius: number=0.5;
implementation

uses gl,glu,math;
{ TAtomRenderer }

constructor TAtomRenderer.create(oglRenderer: TOGLRenderer);
var i:integer;
begin
  ogl:=oglRenderer;
  //Speicherplatz für Atome reservieren
  FillChar(atoms,sizeof(atoms),0);
  atoms[AT_H]:=TOGLGLUObject.Create;
  atoms[AT_C]:=TOGLGLUObject.Create;
  atoms[AT_N]:=TOGLGLUObject.Create;
  atoms[AT_O]:=TOGLGLUObject.Create;
  link:=TOGLGLUObject.create();
  //Modell aktivieren
  useModel(0);
end;

procedure TAtomRenderer.draw(mol: TMolecule);
const R_MASK=$0000FF;
      G_MASK=$00FF00;
      B_MASK=$FF0000;
var i,j,k:integer;
    typ:longint;
    atom:TAtom;
    len:float;
    p: TVector3f;
begin
  //Hintergrund leeren
  glClearColor((atomColors[colorScheme][-1] and R_MASK)/R_MASK,
               (atomColors[colorScheme][-1] and G_MASK)/G_MASK,
               (atomColors[colorScheme][-1] and B_MASK)/B_MASK,0);
  glClear(GL_COLOR_BUFFER_BIT or GL_DEPTH_BUFFER_BIT);
  //Molekülmittelpunkt auswählen
  glTranslatef(-mol.center['x'],-mol.center['y'],-mol.center['z']);
  //Atome der Reihe nach zeichnen
  for i:=0 to mol.atomCount-1 do begin
    atom:=mol.atoms[i];
    typ:=atom.typ;
    p:=vec2f(atom.p);
    case usesModel of
    0: //Einfache Kugel beim Kalottenmodell
       ogl.drawAt(atoms[typ],p['x'],p['y'],p['z']);
    1: begin //Lage der Bindungsstäbe muss berechnet werden
      glTranslatef(p['x'],p['y'],p['z']);
      for j:=1 to atom.james do
        if dword(atom.bond[j].a)>dword(atom) then //Bindungssortierung
        begin //(Trick, damit jede Bindung nur einmal gezeichnet wird)
          glPushMatrix; //Position sichern
          len:=sqrt(veclensqr(vecsub(vec2f(atom.bond[j].a.p),p))); //Bindungslänge
          {        y     z                OpenGL Koordinatensystem
                   y   z   *
                   0 z
          xxxxxxxx000xxxxxxxxx
                 z 0
               z   y
             z     y
          }
          //Andere  Atomkernkoordinaten von globalen XYZ-KO in ein um diesen Kern
          //zentriertes Kugel-KO-System transformieren und den Stab entsprechend
          //den Kugelkoordinaten rotieren
          if not iszero(atom.bond[j].a.p['x']-p['x']) then
            glRotatef(arctan2(atom.bond[j].a.p['y']-atom.p['y'],
                              atom.bond[j].a.p['x']-p['x'])*180/pi,0,0,1)
           else
            glRotatef(90*sign(atom.bond[j].a.p['y']-atom.p['y']),0,0,1);
          glRotatef(arccos(max(-1,min(1,(atom.bond[j].a.p['z']-atom.p['z'])/len)))*
                    180/pi,0,1,0);
          glScalef(1.1-atom.bond[j].count*0.1,1.1-atom.bond[j].count*0.1,len);
          //Ein Zylinder pro Bindung
          case atom.bond[j].count of
            1: ogl.draw(link);
            2: begin //Nebeneinander, leicht verschoben
                 glTranslatef(bondWidth,0,0);
                 ogl.draw(link);
                 glTranslatef(-2*bondWidth,0,0);
                 ogl.draw(link);
               end;
            3: begin //Verschoben in Dreiecksstruktur
                 glTranslatef(bondWidth,bondWidth*0.7,0);
                 ogl.draw(link);
                 glTranslatef(-2*bondWidth,0,0);
                 ogl.draw(link);
                 glTranslatef(bondWidth,-bondWidth*1.7,0);
                 ogl.draw(link);
               end;
          end;
          glPopMatrix;
        end;
      //Kugel zeichnen
      ogl.draw(atoms[typ]);
      glTranslatef(-p['x'],-p['y'],-p['z']);
      end;
    end;
  end;
  glTranslatef(mol.center['x'],mol.center['y'],mol.center['z']);
end;

destructor TAtomRenderer.destroy;
var i:integer;
begin
  for i:=low(atoms) to high(atoms) do
    if atoms[i]<>nil then
      atoms[i].free;
  inherited destroy;
end;

//Model laden
procedure TAtomRenderer.useModel(const m,c: longint);
var i:integer;
begin
  if c<>-1 then colorScheme:=c; //Farbschema wechseln
  usesModel:=m;
  case m of
    0: for i:=low(atoms) to high(atoms) do //Kalottenmodell
         if atoms[i]<>nil then begin
           atoms[i].SetColor(atomColors[colorScheme,i]);     //Farbe bestimmen
           atoms[i].setToSphere(vanderWaalsRadius[i],50,50); //Kugel mit Radius
         end;
     1: begin //Kugel-Stab-modell
       //Zylinder erzeugen
       link.setColor(atomColors[colorScheme][0]);
       link.setToCylinder(bondWidth,1,10,10);
       //Kugeln erzeugen
       for i:=low(atoms) to high(atoms) do
         if atoms[i]<>nil then begin
           atoms[i].SetColor(atomColors[colorScheme,i]);
           atoms[i].setToSphere(kernRadius,15,15);
         end;
     end;
  end;
end;
end.

