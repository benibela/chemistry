{Dieser Programmteil steuert die Interaktion mit dem Benutzer.
 Hierzu werden seine Einnahmen über eine LCL-Form entgegen genommen,
 an die anderen Bestandteile des Programmes weitergereicht und
 die Ergebnisse ausgegeben.
}

unit gui;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, LResources, Forms, Controls, Graphics, Dialogs,StdCtrls,Buttons,ExtCtrls, ComCtrls,
  extmath,atoms,molecules,
  w32initopengl,  oglDrawing,atomicdrawing,energyCalculator,geometryCalculator,lewisrenderer;

type

  { TForm1 }

  TForm1 = class(TForm)
    //Verwendete Designelemente
    colorbox: TComboBox;
    energyCalculateBtn: TButton;
    Button3: TButton;
    Label7: TLabel;
    molSaveBtn: TButton;
    geooptbtn: TButton;
    molLoadBtn: TButton;
    modelbox: TComboBox;
    atomNr: TEdit;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label6: TLabel;
    OpenDialog1: TOpenDialog;
    optimizeStepLabel: TLabel;
    loadOGL1: TButton;
    lewisPanel: TPanel;
    SaveDialog1: TSaveDialog;
    Splitter2: TSplitter;
    viewBarRotX: TScrollBar;
    viewBarRotY: TScrollBar;
    viewBarDistance: TScrollBar;
    viewDistance: TLabel;
    viewRotX: TLabel;
    viewRotY: TLabel;
    Label5: TLabel;
    label9: TLabel;
    loadOGLBtn: TButton;
    status: TMemo;
    Panel1: TPanel;
    Panel2: TPanel;
    Panel3: TPanel;
    optimizeStepBar: TScrollBar;
    Splitter1: TSplitter;
    Timer1: TTimer;
    //Diese Funktionen werden automatisch bei bestimmten Ereignissen
    //aufgerufen
    procedure energyCalculateBtnClick(Sender: TObject);
    procedure Button3Click(Sender: TObject);
    procedure loadOGLBtnClick(Sender: TObject);
    procedure molSaveBtnClick(Sender: TObject);
    procedure geooptbtnClick(Sender: TObject);
    procedure molLoadBtnClick(Sender: TObject);
    procedure modelboxChange(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure FormDestroy(Sender: TObject);
    procedure optimizeStepBarScroll(Sender: TObject; ScrollCode: TScrollCode;
      var ScrollPos: Integer);
    procedure Panel1MouseDown(Sender: TOBject; Button: TMouseButton;
      Shift: TShiftState; X, Y: Integer);
    procedure Panel1MouseMove(Sender: TObject; Shift: TShiftState; X, Y: Integer
      );
    procedure Panel1Resize(Sender: TObject);
    procedure viewBarDistanceScroll(Sender: TObject; ScrollCode: TScrollCode;
      var ScrollPos: Integer);
    procedure Timer1Timer(Sender: TObject);
  private
    { private declarations }
    cancelGeoOptimizing:boolean;
    procedure optimizingProgress(sender: tobject;energy:number;steps: longint;
                                 var cancel: boolean);
  public
    { public declarations }
    //Eigentliche Berechnungsklassen, definiert in anderen Dateien
    oglW32: TOGLW32Base;          //oglDrawing.pas
    oglRenderer: TOGLRenderer;    //oglDrawing.pas
    atomRenderer: TAtomRenderer;  //atomicDrawing.pas
    lewisRenderer: TLewisRenderer;//lewisRenderer.pas
    scfEnergyCalculator: TSCFEnergyCalculator; //energyCalculator.pas
    geometryOptimizer: TGeometryOptimizer; //geometryOptimizer.pas
    currentMolecule:TMolecule;    //molecules.pas
    
    //Position im 3D-Raum
    z,ry,rx: single; //abstand,rotation um x,y
    cx,cy:longint;//klickpos
    cz,cry,crx: single; //abstand,rotation um x,y beim klicken
    //Ausgaben an den Benutzer umwandeln
    procedure viewChanged;
    procedure publicLog(s:string);
    function lenToUser(n: number): string;
    function energyToUser(n: number): string;
    function angleToUser(n: number): string;
    
  end;

var
  Form1: TForm1; 


implementation

uses gl,math;
{ TForm1 }


//Start
procedure TForm1.FormCreate(Sender: TObject);
begin
  z:=5;
  oglW32:=nil;
  currentMolecule:=TMolecule.create;
  scfEnergyCalculator:=TSCFEnergyCalculator.Create;
  geometryOptimizer:=TGeometryOptimizer.create(scfEnergyCalculator);
  lewisRenderer:=TLewisRenderer.create(lewisPanel);
  lewisRenderer.useMolecule(currentMolecule);
end;

//Ende
procedure TForm1.FormDestroy(Sender: TObject);
begin
  if oglW32<>nil then begin
    atomRenderer.free;
    oglRenderer.free;
    oglW32.free;
  end;
  if lewisRenderer<>nil then
    lewisRenderer.free;
  scfEnergyCalculator.free;
  geometryOptimizer.free;
end;

//Rückmeldung an den Benutzer nach Änderung
procedure TForm1.modelboxChange(Sender: TObject);
begin
  if (atomRenderer<>nil) then begin
    atomRenderer.useModel(modelbox.ItemIndex,colorBox.itemindex);
    if modelbox.ItemIndex<>atomRenderer.currentModel then begin
      case modelbox.ItemIndex of
        0: z-=5;
        1: z+=5;
      end;
    end;
    viewChanged;
  end;
end;

procedure TForm1.optimizeStepBarScroll(Sender: TObject;
  ScrollCode: TScrollCode; var ScrollPos: Integer);
begin
  optimizeStepLabel.Caption:=inttostr(optimizeStepBar.Position);
end;

procedure TForm1.optimizingProgress(sender: tobject; energy: number; steps: longint;
  var cancel: boolean);
begin
  Application.ProcessMessages;
  publicLog('Noch maximal '+inttostr(steps)+' Schritte, aktuelle Energie: '+
            energyToUser(energy));
  cancel:=cancelGeoOptimizing;
end;

procedure TForm1.viewChanged;
begin
  viewBarDistance.Position:=trunc(-z);
  viewDistance.Caption:=inttostr(-trunc(z));
  if rx<-180 then rx:=rx+360;
  if rx>180 then rx:=rx-360;
  viewRotX.Caption:=inttostr(trunc(rx))+'°';
  viewBarRotX.Position:=trunc(rx);
  if ry<-180 then ry:=ry+360;
  if ry>180 then ry:=ry-360;
  viewRotY.Caption:=inttostr(trunc(ry))+'°';
  viewBarRotY.Position:=trunc(ry);
end;


//Mausbewegung im 3D-Raum erkennen
//Rotation um x und y-Achse je nach Mausbewegung
procedure TForm1.Panel1MouseDown(Sender: TOBject; Button: TMouseButton;
  Shift: TShiftState; X, Y: Integer);
begin
  cx:=x;
  cy:=y;
  cz:=z;
  crx:=rx;
  cry:=ry;
end;

procedure TForm1.Panel1MouseMove(Sender: TObject; Shift: TShiftState; X,
  Y: Integer);
begin
  if ssRight in shift then begin
    z:=cz+(cy-y)/50;
    viewChanged;
  end else if ssleft in shift then begin
    rx:=crx+y-cy;
    ry:=cry+x-cx;
    viewChanged;
  end;
end;

procedure TForm1.Panel1Resize(Sender: TObject);
begin
  if oglRenderer<>nil then
    oglRenderer.resize(panel1.width,panel1.height);
end;

procedure TForm1.viewBarDistanceScroll(Sender: TObject; ScrollCode: TScrollCode;
  var ScrollPos: Integer);
begin
  z:=-viewBarDistance.Position;
  rx:=viewBarRotX.Position;
  ry:=viewBarRotY.Position;
  viewChanged;
end;

//Periodisch zeichnen
procedure TForm1.Timer1Timer(Sender: TObject);
begin
  if (oglW32<>nil)and (atomRenderer<>nil) and (currentMolecule<>nil) then begin
    glMatrixMode(GL_modelview);
    glLoadIdentity;
    glTranslatef(0,0,z);
    glRotatef(ry,0,1,0);
    glRotatef(rx,1,0,0);
    glColor4f(1,1,1,0.5);
    atomRenderer.draw(currentMolecule);
    oglW32.endScene;
  end;
end;


//Informationsausgaben
procedure TForm1.publicLog(s: string);
begin
  status.lines.add(FormatDateTime('hh:mm:ss.zzzz',time)+': '+s);
end;

function TForm1.lenToUser(n: number): string;
begin
  if abs(n)<1e-5 then exit('0 nm');
  result:=format('%2.3f',[n*0.05291772])+' nm';
end;

function TForm1.energyToUser(n: number): string;
begin
  if abs(n)<1e-5 then exit('0 kJ/mol');
  result:=format('%2.5f',[n*(4.3594e-18)*(6.02214e23)/1000])+' kJ/mol';
end;

function TForm1.angleToUser(n: number): string;
begin
  if abs(n)<1e-5 then exit('0°');
  result:=format('%2.1f',[180*n/pi])+'°';
end;


//Energieberechnung aufrufen
procedure TForm1.energyCalculateBtnClick(Sender: TObject);
var t:cardinal;
begin
  scfEnergyCalculator.initialize(currentMolecule);
  publicLog('Energieberechnung gestartet');
  publicLog('Energieberechnung beendet: '+
            energyToUser(scfEnergyCalculator.calculate));
end;

//Informationen über ausgewähltes Atom anzeigen
procedure TForm1.Button3Click(Sender: TObject);
var i,j:longint;
    a: tatom;
    s:string;
    b:array[1..4] of TVector3n;
begin
  a:=tatom(currentMolecule.atoms[strtoint(atomNr.text)]);
  publicLog('Atom '+atomNr.text+' '+ATOM_NAMES[a.typ]+' hat '+
            inttostr(a.james)+' unterschiedliche Bindungen: ');
  //Bindungslängen
  for i:=1 to a.james do begin
    b[i]:=vecsub(a.p,a.bond[i].a.p);
    publicLog(inttostr(i)+': '+ATOM_NAMES[a.bond[i].a.typ]+' Länge: '+
              lenToUser(sqrt(veclensqr(b[i])))
              );
  end;
  //Bindungswinkel (nach Formelsammlung)
  for i:=1 to a.james do begin
    s:='Winkel '+inttostr(i)+': ';
    for j:=1 to a.james do
      s:=s+angleToUser(arccos(((b[i]['x']*b[j]['x']+b[i]['y']*b[j]['y']+
                                b[i]['z']*b[j]['z']))/
                                sqrt(veclensqr(b[i])*veclensqr(b[j]))))+#9;
    publicLog(s);
  end;
end;

//OpenGL laden
procedure TForm1.loadOGLBtnClick(Sender: TObject);
//Licht position
const
  mat_specular   : Array[0..3] of GlFloat = (1.0, 1.0, 1.0, 1.0);
  mat_shininess  : Array[0..0] of GlFloat = (50.0);
  mat_ambient    : Array[0..3] of GlFloat = (1, 1, 1, 1.0);
  mat_diffuse    : Array[0..3] of GlFloat = (1, 1, 1, 1.0);

  light_position : Array[0..3] of GlFloat = (5.0, 5.0, 10.0, 1.0);
  light_ambient  : Array[0..3] of GlFloat = (0.1, 0.1, 0.1, 1.0);
  light_diffuse  : Array[0..3] of GlFloat = (0.8, 0.8, 0.8, 1.0);
//const
//    att: GLfloat=0.5;
begin
  //Objekte erzeugen
  oglW32:=TOGLW32Base.create(panel1.handle);
  oglRenderer:=TOGLRenderer.Create;
  oglRenderer.init();
  oglRenderer.resize(panel1.width,panel1.Height);
  atomRenderer:=TAtomRenderer.create(oglRenderer);
  atomRenderer.useModel(modelbox.ItemIndex);
  z:=-5;
  //Standardopenglbeleuchtung
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  glLoadIdentity;

  glMaterialfv(GL_FRONT, GL_SPECULAR,  @mat_specular[0]);
  glMaterialfv(GL_FRONT, GL_SHININESS, @mat_shininess[0]);
  glMaterialfv(GL_FRONT, GL_AMBIENT,   @mat_ambient[0]);
  glMaterialfv(GL_FRONT, GL_DIFFUSE,   @mat_diffuse[0]);

  glLightfv(GL_LIGHT0, GL_AMBIENT,  @light_ambient[0]);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  @light_diffuse[0]);
  glLightfv(GL_LIGHT0, GL_POSITION, @light_position[0]);
  
  if currentMolecule<>nil then
    panel1.caption:='Zeichne Molekül...'
   else
    panel1.caption:='Kein Molekül geladen';
  loadOGLBtn.Visible:=false;
  
  publicLog('3D Renderer geladen');
end;

//Geometrieoptimierung aufrufen oder abbrechen
procedure TForm1.geooptbtnClick(Sender: TObject);
var c:char;
    i:longint;
    s:string;
begin
  if geooptbtn.tag=0 then begin
    //Möglichkeit zum Abbrechen bieten
    geooptbtn.tag:=1;
    geooptbtn.Caption:='Optimierung abbrechen';
    molSaveBtn.Enabled:=false;
    molLoadBtn.Enabled:=false;
    energyCalculateBtn.Enabled:=false;
    loadOGLBtn.Enabled:=false;

    //Aufrufen
    cancelGeoOptimizing:=false;
    publicLog('Geometrieoptimierung gestartet');
    geometryOptimizer.maxSteps:=optimizeStepBar.position;
    geometryOptimizer.progress:=@optimizingProgress;
    geometryOptimizer.calculate;
    //Fertig, Ergebnisse ausgeben
    if cancelGeoOptimizing then
      publicLog('Geometrieoptimierung vom Benutzer abgebrochen')
    else publicLog('Geometrieoptimierung beendet: ');
    for i:=0 to currentMolecule.atomCount-1 do begin
      s:='Atom '+inttostr(i)+': '+
         ATOM_NAMES[currentMolecule.atoms[i].typ]+': ';
      for c:='x' to 'z' do
        s:=s+c+': '+lenToUser(currentMolecule.atoms[i].p[c])+'  ';
      publicLog(s);
    end;
    //Abbrechen verhindern
    molSaveBtn.Enabled:=true;
    molLoadBtn.Enabled:=true;
    energyCalculateBtn.Enabled:=true;
    loadOGLBtn.Enabled:=true;
    geooptbtn.Tag:=0;
    geooptbtn.caption:='Geometrie optimieren';
  end else cancelGeoOptimizing:=true; //Abbrechen
end;

//Dateibehandlung
procedure TForm1.molSaveBtnClick(Sender: TObject);
begin
  if SaveDialog1.Execute then currentMolecule.saveToFile(saveDialog1.filename);
end;

procedure TForm1.molLoadBtnClick(Sender: TObject);
begin
  if openDialog1.Execute then begin
    currentMolecule.loadfromFile(openDialog1.filename);
    publicLog('Molekül geladen');
    lewisRenderer.render();
  end;
end;



initialization
  {$I gui.lrs}

end.

