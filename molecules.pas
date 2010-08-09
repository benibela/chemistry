{Dieser Programmteil enthält alle wesentlichen Funktionen zur Veränderung von
 Molekülen und zur Berechnung ihrer Eigenschaften
 w
}
unit molecules;

{$mode objfpc}{$H+}



interface

uses
  Classes, SysUtils,extmath,atoms,stosim;
type
  { TMolecule }
  //DIE Molekülklasse
  //TODO3: Klasse trennen in Berechnungsklasse (2 für Geometrie, 1 für Energie)
  //und Speicherklasse
  TMolecule=class
  private
    _atoms: TList;       //Atomliste
    _center: TVector3n;  //Schwerpunkt
    function getAtomCount: longint;
  public
    //Berechnet den Schwerpunkt neu
    procedure geometryChanged();

    constructor create;
    destructor destroy;override;
    
    //Überprüft ob alle Schalen geschlossen sind
    function closedShell():boolean;
    
    function addAtom(const typ: longint):tatom; //Fügt ein Atom in die Liste ein
    function addAtom(const typ: longint;const x,y,z: number):tatom;//an Position
    function getAtom(const i:longint):TAtom;inline;
    
    function getAtomCount():longint;inline;
    
    property atoms[index: longint]: TAtom read getAtom; default;
    property atomCount: longint read getAtomCount;
    property center: TVector3n read _center;
    
    procedure deleteAtom(a:Tatom);    //Löscht ein Atom
    procedure deleteBonds(a:Tatom);   //Löscht alle Bindungen am Atom
    procedure deleteBonds(a,b:Tatom); //Löscht alle Bindungen zwischen Atomen
    procedure addBond(a,b:tatom;const count:longint=1);//Fügt eine Bindung hinzu
    
    function getIUPACName():string;
    
    //Speichert und läd das Molekül
    procedure saveToFile(fn: string);
    procedure loadFromFile(fn: string);
end;
  
implementation
{TMolecule }

constructor TMolecule.create;
begin
  _atoms:=TList.create;
end;

destructor TMolecule.destroy;
var i:longint;
begin
  for i:=0 to atomCount-1 do
    atoms[i].free;
  _atoms.free;

  inherited destroy;
end;


//Rät eine einigermaßen gute Koeffizientenmatrix
//Berechnet den Schwerpunkt neu
procedure TMolecule.geometryChanged();
var i:longint;
begin
  _center:=nullvector3n;
  for i:=0 to atomCount-1 do
    _center:=vecadd(center,atoms[i].p);
  vecscale(_center,1/atomCount);
end;

function TMolecule.getAtomCount: longint;
begin
  result:=_atoms.count;
end;




//Überprüft ob alle Schalen geschlossen sind
function TMolecule.closedShell(): boolean;
var i,j,t:longint;
begin
  result:=true;
  for i:=0 to atomCount-1 do begin
    t:=ATOM_CLOSED_BONDS[atoms[i].typ];
    for j:=1 to atoms[i].james do
      t-=atoms[i].bond[j].count;
    if t<>0 then exit(false);
  end;
end;

//Fügt ein Atom hinzu
function TMolecule.addAtom(const typ: longint): tatom;
begin
  //Zufallsposition suchen und Zwilling rufen
  Randomize;
  result:=addAtom(typ,random(10)-5,random(10)-5,random(10)-5);
end;

//Fügt ein Atom hinzu
function TMolecule.addAtom(const typ: longint; const x, y, z: number): tatom;
begin
  //Atom erzeugen
  result:=tatom.create;
  //Werte einfüllen
  result.index:=atomCount;
  Result.typ:=typ;
  result.p['x']:=x;
  result.p['y']:=y;
  result.p['z']:=z;
  //speichern
  _atoms.add(result);
  geometryChanged;
end;

function TMolecule.getAtom(const i: longint): TAtom; inline;
begin
  result:=TAtom(_atoms[i]);
end;

function TMolecule.getAtomCount(): longint; inline;
begin

end;


//Atom löschen
procedure TMolecule.deleteAtom(a: Tatom);
var i:longint;
begin
  //Alle Bindungen löschen
  for i:=a.james downto 1 do
    deleteBonds(a,a.bond[i].a);
  //Atom entfernen
  _atoms.delete(a.index);
  a.free;
  //Alle Indizes aktualisieren
  for i:=0 to atomCount-1 do
    atoms[i].index:=i;
  geometryChanged;
end;

//Bindungen löschen
procedure TMolecule.deleteBonds(a: Tatom);
var i:longint;
begin
  for i:=a.james downto 1 do
    deleteBonds(a,a.bond[i].a); //Alle Bindungen einzeln löschen
end;

//Bindung löschen
procedure TMolecule.deleteBonds(a, b: Tatom);
var i:longint;
begin
  //Bindung bei a suchen und löschen
  for i:=a.james downto 1 do
    if a.bond[i].a=b then begin
      if i<>high(b.bond) then
        move(a.bond[i+1],a.bond[i],sizeof(a.bond[1])*(a.james-i));
      a.james-=1;
    end;
  //Bindung bei b suchen und löschen
  for i:=b.james downto 1 do
    if b.bond[i].a=a then begin
      if i<>high(b.bond) then
        move(b.bond[i+1],b.bond[i],sizeof(b.bond[1])*(b.james-i));
      b.james-=1;
    end;
end;

//Bindung hinzufügen
procedure TMolecule.addBond(a, b: tatom;const count:longint=1);
var i,j:longint;
begin
  //Suche nach vorhandener Bindungen
  for i:=1 to a.james do
    if a.bond[i].a=b then begin
      a.bond[i].count+=count;
      for j:=1 to b.james do
        if b.bond[j].a=a then
          b.bond[j].count+=count;
      exit;
    end;
  //Sonst neue hinzufügen
  a.james+=1;
  b.james+=1;
  a.bond[a.james].a:=b;
  a.bond[a.james].count:=count;
  b.bond[b.james].a:=a;
  b.bond[b.james].count:=count;
end;

function TMolecule.getIUPACName(): string;
begin

end;

//Alle Daten binär in eine Datei schreiben
procedure TMolecule.saveToFile(fn: string);
var f: TFileStream;
    i,j:longint;
begin
 { f:=TFileStream.Create(fn,fmCreate);
  f.WriteDWord(atomCount);
  for i:=0 to atomCount-1 do begin
    f.WriteDWord(atoms[i].index);
    f.WriteDWord(atoms[i].typ);
    f.WriteDWord(atoms[i].james);
    f.WriteBuffer(atoms[i].p['x'],sizeof(atoms[i].p));
    f.WriteDWord(atoms[i].flatX);
    f.WriteDWord(atoms[i].flatY);
    for j:=1 to atoms[i].james do begin
      f.WriteDWord(atoms[i].bond[j].a.index);
      f.WriteDWord(atoms[i].bond[j].count);
    end;
  end;
  f.free;}
end;

//Alle Daten binär aus einer Datei lesen
procedure TMolecule.loadFromFile(fn: string);
var f: TFileStream;
    c,i,j:longint;
begin
{  for i:=0 to atomCount-1 do
    atoms[i].free;
  atoms.clear;
  f:=TFileStream.Create(fn,fmOpenRead);
  atomCount:=f.ReadDWord;
  for i:=0 to atomCount-1 do
    atoms[i]:=tatom.create;
  for i:=0 to atomCount-1 do begin
    atoms[i].index:=f.ReadDWord;
    atoms[i].typ:=f.ReadDWord;
    atoms[i].james:=f.ReadDWord;
    f.ReadBuffer(atoms[i].p['x'],sizeof(atoms[i].p));
    atoms[i].flatX:=f.ReadDWord;
    atoms[i].flatY:=f.ReadDWord;
    for j:=1 to atoms[i].james do begin
      atoms[i].bond[j].a:=atoms[f.ReadDWord]);
      atoms[i].bond[j].count:=f.ReadDWord;
    end;
  end;
  f.free;
  geometryChanged();}
end;


end.

