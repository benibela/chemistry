{ Dieser Programmteil zeichnet eine Lewisformel, die der Benutzer
  beliebig ändern kann}
unit lewisRenderer;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,Controls,graphics,Extctrls,menus,buttons,ComCtrls,
  atoms,molecules,geometryCalculator;

type
  { TLewisRenderer }

  TLewisRenderer=class
    status: (sNormal,sAddBond,sAddAtom);  //Reaktion auf Benutzer (sAddBond verbindet Atome)
    molecule: TMolecule;         //Angezeigtes Molekül
    symHeight,symWidth: longint; //Größe eines Symbols
    
    atomToAdd: string;
    selected: tatom;             //Ausgewähltes Atom
    mx,my,ox,oy:longint;         //Mausposition
    atomPopup: TPopupMenu;       //Menü zur Benutzerinteraktion
    geometryCalculator:TGeometryTrivial;
    
    function addAtomAt(const name: string; x, y: longint):TAtom;
    
    //Ereignisse des Menüs
    procedure itemAddClick(Sender: TObject);
    procedure itemBondAddClick(Sender: TObject);
    procedure itemBondDeleteClick(Sender: TObject);
    procedure itemDeleteClick(Sender: TObject);
    procedure BtnAddAtomClick(Sender: TObject);
    procedure paintBoxPaint(Sender: TObject);
    //Benutzerinteraktion
    procedure paintboxMouseDown(sender:TObject;Button: TMouseButton;Shift: TShiftState;x,y:longint);
    procedure paintboxMouseMove(sender:TObject;Shift: TShiftState;x,y:longint);
  public
    defX, defY:longint;   //Nullpunkt
    paintBox: tpaintBox;  //Ausgabefläche
    constructor create(aparent:TWinControl);
    destructor destroy;override;
    procedure useMolecule(mol: TMolecule);
    procedure rearrange();//Ordnet die Atome intelligen (unbenutzt)
    procedure render();

    property usedMolecule: TMolecule read molecule write useMolecule;
  end;

implementation
uses forms,math,Dialogs;
{ TLewisRenderer }

//Ausgewähltes Atom löschen und Vorschaugeometrie aktualisieren
procedure TLewisRenderer.itemDeleteClick(Sender: TObject);
begin
  if selected<>nil then molecule.deleteAtom(selected);
  geometryCalculator.initialize(molecule);
  geometryCalculator.calculate;
  render;
end;

procedure TLewisRenderer.BtnAddAtomClick(Sender: TObject);
begin
  status:=sAddAtom;
  atomToAdd:=(sender as TButton).Caption;
  selected:=nil;
end;

procedure TLewisRenderer.paintBoxPaint(Sender: TObject);
begin
  render();
end;

//Alle Bindungen am Ausgewählten Atom löschen und Vorschaugeometrie aktualisieren
procedure TLewisRenderer.itemBondDeleteClick(Sender: TObject);
begin
  if selected=nil then exit;
  molecule.deleteBonds(selected);
  geometryCalculator.initialize(molecule);
  geometryCalculator.calculate;
  render;
end;

//Moduswechseln (s.u.)
procedure TLewisRenderer.itemBondAddClick(Sender: TObject);
begin
  if selected=nil then exit;
  status:=sAddBond;
  render;
end;

function TLewisRenderer.addAtomAt(const name: string; x, y: longint):TAtom;
var i:longint;
begin
  for i:=1 to high(ATOM_ABBREV) do
    if uppercase(ATOM_ABBREV[i])=uppercase(name) then begin
      result:=molecule.addAtom(i);
      with result do begin
        flatX:=mx;
        flatY:=my;
      end;
      geometryCalculator.initialize(molecule);
      geometryCalculator.calculate;
      render;
      exit;
    end;
  result:=nil;
end;

//Nach Typ des neuen Atoms fragen und Vorschaugeometrie aktualisieren
procedure TLewisRenderer.itemAddClick(Sender: TObject);
var nt:string;
    i:longint;
begin
  nt:='H';
  if InputQuery('YAMS','Bitte geben Sie ein, welches Atom eingefügt werden soll. '+
                '(H,C,N oder O)',nt) then begin
    if addAtomAt(nt,mx,my)=nil then ShowMessage('Atom unbekannt');
  end;

end;

//Erzeugung (Menüerzeugung)
constructor TLewisRenderer.create(aparent:TWinControl);
var item:TMenuItem;
    i:longint;
begin
  defX:=50;
  defY:=50;
  symHeight:=15;
  symWidth:=10;
  atomPopup:=TPopupMenu.Create(Application.MainForm);
  item:=TMenuItem.Create(atomPopup);
  item.Caption:='Atom Löschen';
  item.OnClick:=@itemDeleteClick;
  atomPopup.Items.Add(item);
  item:=TMenuItem.Create(atomPopup);
  item.Caption:='Atombindungen löschen';
  item.OnClick:=@itemBondDeleteClick;
  atomPopup.Items.Add(item);
  item:=TMenuItem.Create(atomPopup);
  item.Caption:='Atombindung hinzufügen';
  item.OnClick:=@itemBondAddClick;
  atomPopup.Items.Add(item);
  item:=TMenuItem.Create(atomPopup);
  item.Caption:='Atom hinzufügen';
  item.OnClick:=@itemAddClick;
  atomPopup.Items.Add(item);
  geometryCalculator:=TGeometryTrivial.Create;
  
  for i:=1 to high(ATOM_ABBREV) do
    if ATOM_ABBREV[i]<>'' then
      with TButton.Create(aparent) do begin
        Parent:=aparent;
        top:=0;
        left:=20*i;
        Height:=20;
        width:=20;
        caption:=ATOM_ABBREV[i];
        OnClick:=@BtnAddAtomClick;
      end;

  paintBox:=TPaintBox.Create(aparent);
  paintBox.Parent:=aparent;
  paintBox.Anchors:=[akTop,akLeft,akRight,akBottom];
  paintBox.top:=20;
  paintBox.left:=0;
  paintBox.Height:=aparent.ClientHeight-paintBox.top;
  paintBox.Width:=aparent.ClientWidth;
  paintBox.OnPaint:=@paintBoxPaint;
  paintBox.OnMouseDown:=@paintboxMouseDown;
  paintBox.OnMouseMove:=@paintboxMouseMove;
end;

destructor TLewisRenderer.destroy;
begin
  geometryCalculator.free;
  atomPopup.free;
  inherited destroy;
end;

procedure TLewisRenderer.useMolecule(mol: TMolecule);
begin
  molecule:=mol;
  rearrange;
end;

//Flache Atomanordnung
procedure TLewisRenderer.rearrange();
const flatBondLength=20;
var arranged: array of boolean;
  //Ein Atom und Nachbarn anordnen
  procedure rar(a: tatom; dir: longint);
    //Bindungspriorität
    function comp(const a,b:TBond): longint;
    const AP:array[1..8] of longint=(0,0,0,0,0,3,1,2);
    begin
      if arranged[a.a.index] then result:=1
      else if arranged[b.a.index] then result:=-1
      else if AP[a.a.typ]<AP[b.a.typ] then result:=1
      else if AP[a.a.typ]>AP[b.a.typ] then result:=-1
      else if a.count>b.count then result:=-1
      else if a.count<b.count then result:=1
      else result:=0;
    end;
  //Mögliche Richtungen
  const DirRel:array[1..3] of array[0..3]of longint=
               ((2,3,0,1),(3,0,1,2),(1,2,3,0));
        RealDir: array[0..3] of array['x'..'y'] of longint=
                 ((-1,0),(0,-1),(1,0),(0,1));
  var i,j:longint;
      realBondCount: longint;
      fromBond: longint;
      temp:^tbond;
      b: array[1..4] of ^TBond;
  begin
    //Bindungen nach Priorität sortieren
    dir:=DirRel[1,dir];
    for i:=1 to a.james do
      b[i]:=@a.bond[i];
    for i:=2 to a.james do
      for j:=i-1 downto 1 do
        if comp(b[j]^,b[j+1]^) > 0 then begin
          temp:=b[j];
          b[j]:=b[j+1];
          b[j+1]:=temp;
        end;
    //Nachbaratome anordnen und von dort aus fortfahren
    for i:=1 to a.james do
      if not arranged[b[i]^.a.index] then begin
        b[i]^.a.flatX:=a.flatX+RealDir[DirRel[i,dir]]['x']*flatBondLength;
        b[i]^.a.flatY:=a.flatY+RealDir[DirRel[i,dir]]['y']*flatBondLength;
        arranged[b[i]^.a.index]:=true;
        rar(b[i]^.a,DirRel[i,dir]);
      end;
  end;
begin
  setlength(arranged,molecule.atomCount);
  if length(arranged)=0 then exit;
  FillChar(arranged,length(arranged)*sizeof(arranged[0]),0);
  arranged[0]:=true;
  //Erstes Atom an Nullpunkt
  with tatom(molecule.atoms[0]) do begin
    flatX:=defX;
    flatY:=defY;
  end;
  rar(tatom(molecule.atoms[0]),0);
end;

//Molekül zeichnen
procedure TLewisRenderer.render();
  //Zeichnet mehrere Linien nebeneinander
  procedure drawLines(x1,y1,x2,y2,count:longint);
  const LINE_SPACE:longint=3;
  var px,py: single;
      lx,ly: single;
      sx,sy,ex,ey: single;
      sizesub: single;
      len:single;
      i:longint;
  begin
    //Linienvektor
    lx:=x2-x1;
    ly:=y2-y1;
    //Fälle
    if count<=1 then begin       //Keine Linie
      px:=0;
      py:=0;
    end else if x1=x2 then begin //Steigung unendlich
      px:=1;
      py:=0;
    end else begin               //Vektor über Steigung
      //px*lx:=-py*ly da senkrecht
      px:=-ly/lx;
      py:=1;
      len:=sqrt(sqr(px)+sqr(py));
      px:=px/len;
      py:=py/len;
    end;
    //Vektor normieren
    len:=sqrt(sqr(lx)+sqr(ly));
    lx/=len;
    ly/=len;
    //Ausgangspunkte senkrecht verschieben, LINE_SPACE/2 Pixel pro Bindung
    //=> Symmetrie um Achse
    sizesub:=(abs(symHeight*ly)+abs(symWidth*lx))/2;
    sx:=x1+(sizesub+1)*lx-LINE_SPACE*px*(count-1)/2;
    sy:=y1+(sizesub+1)*ly-LINE_SPACE*py*(count-1)/2;
    ex:=x2-sizesub*lx-LINE_SPACE*px*(count-1)/2;
    ey:=y2-sizesub*ly-LINE_SPACE*py*(count-1)/2;
    //Bindungen nebeneinander zeichnen und dabei Verschiebung langsam umkehren
    if isnan(sx) or isnan(sy) or isnan(ex) or isnan(ey) then exit;
    for i:=1 to count do begin
      paintBox.canvas.line(round(sx),round(sy),round(ex),round(ey));
      sx+=LINE_SPACE*px;
      sy+=LINE_SPACE*py;
      ex+=LINE_SPACE*px;
      ey+=LINE_SPACE*py;
    end;
  end;
  
  procedure drawAtomSymbolAt(sym:string;x,y:longint);
  begin
    with paintBox.canvas do
      TextOut(x-textWidth(sym) div 2,y-TextHeight(sym) div 2,sym);
  end;
var i,j:longint;
    totalBondUsed: longint;
    rect:TRect;
begin
  rect:=paintBox.ClientRect;
  with paintBox.canvas do begin
    //Hintergrund
    Brush.color:=clWhite;
    Brush.Style:=bsSolid;
    FillRect(rect);
    Brush.Style:=bsClear;
    pen.color:=clBlack;
    font.Height:=symHeight;
    font.Style:=[fsBold];
    totalBondUsed:=0;
    //Atome zeichnen
    for i:=0 to molecule.atomCount-1 do
      with tatom(molecule.atoms[i]) do begin
        if tatom(molecule.atoms[i])=selected then font.color:=clBlue
        else font.color:=clBlack;
        //Symbol
        drawAtomSymbolAt(ATOM_ABBREV[typ],flatX,flatY);
        //Bindungen zeichnen
        totalBondUsed-=ATOM_CLOSED_BONDS[typ];
        for j:=1 to james do begin
          totalBondUsed+=bond[j].count;
          if index<bond[j].a.index then
            drawLines(flatX,flatY,bond[j].a.flatX,bond[j].a.flatY,bond[j].count);
        end;
      end;
    //Mögliche neue Bindung einzeichnen
    case status of
      sAddAtom: begin
        font.color:=clBlue;
        drawAtomSymbolAt(atomToAdd,mx,my);
      end;
      sAddBond:
        if (selected<>nil) then
         line(selected.flatX,selected.flatY,mx,my);
    end;
    //Hinweise an Benutzer ausgeben
    font.color:=clBlack;
    font.Style:=[];
    TextOut(10,10,'Anzahl Atome: '+IntToStr(molecule.atomCount));
    with rect do
      TextOut(Right-TextWidth('Mit rechter Maustaste klicken für mehr Optionen')-3,
        bottom-TextHeight('Mp')-3,'Mit rechter Maustaste klicken für mehr Optionen');
    font.color:=clRed;
    if totalBondUsed<0 then
      TextOut(10,30,'Zu wenig Bindungen!');
    if totalBondUsed>0 then
      TextOut(10,30,'Zu viele Bindungen!');
  end;
end;

//Benutzer klickt
procedure TLewisRenderer.paintboxMouseDown(sender:TObject;Button: TMouseButton;Shift: TShiftState;x,y:longint);
  function findAtom():TAtom;
  var i:longint;
  begin
    //Im Umkreis von 20 Pixeln nach Atomsymbol suchen
    result:=nil;
    for i:=0 to molecule.atomCount-1 do
      if sqr(tatom(molecule.atoms[i]).flatX-x)+
         sqr(tatom(molecule.atoms[i]).flatY-y) < 20 then begin
        result:=tatom(molecule.atoms[i]);
        ox:=result.flatX;
        oy:=result.flatY;
        break;
      end;
  end;

var newSelect:Tatom;
begin
  case status of
    sNormal:  //Im Standardmodus auswählen
      selected:=findAtom();

    sAddBond: begin //Bindung zwischem ausgewählten und angeklickten Atom erzeugen
      newSelect:=findAtom();
      if newSelect=nil then exit;
      status:=sNormal;
      molecule.addBond(newSelect,selected);
      geometryCalculator.initialize(molecule);
      geometryCalculator.calculate; //Vorschaugeometrie aktualisieren
      selected:=nil;
    end;
    
    sAddAtom: begin
      addAtomAt(atomToAdd,x,y);
      status:=sNormal;
    end;
  end;
  mx:=x;
  my:=y;
  render;
  if Button=mbRight then atomPopup.PopUp();
end;

//Benutzer hat die Maus bewegt
procedure TLewisRenderer.paintboxMouseMove(sender:TObject;Shift: TShiftState;x,y:longint);
begin
  case status of
    sNormal: if ssLeft in shift then begin //Atom wird verschoben
      if selected=nil then exit;
      selected.flatX:=ox+x-mx;
      selected.flatY:=oy+y-my;
      render;
    end;
    sAddBond: begin //Vorschaubindung zeigen
      mx:=x;
      my:=y;
      render;
    end;
    sAddAtom: begin
      mx:=x;
      my:=y;
      render;
    end;
  end;
end;

end.

