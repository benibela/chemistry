{Dieser Programmteil dient zum Anzeigen von 3D-Objekten über
 die Grafikschnittstelle OpenGL
}
unit oglDrawing; 

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,graphics,gl,glu;

type

{ TOGLObject }

//Verwaltet ein einfaches 3D-Objekt
TOGLObject=class
private
  _color: T4fArray;
  list: GLuint;
public
  constructor create();
  procedure setColor(c: T4fArray);overload;
  procedure setColor(c: TColor);overload;
  property color:T4fArray read _color write setcolor;
end;

{ TOGLGLUObject }

//Bieten einen übersichtlichen Zugriff auf die von OpenGLU automatisch erzeugten
//primitiven Objekte (Kugel und Zylinder)
TOGLGLUObject=class(TOGLObject)
private
  quadric: PGLUquadric;
public
  constructor create();
  destructor destroy;override;

  procedure setToSphere(radius:GLdouble;slices,stacks:GLint);
  procedure setToCylinder(baseRadius,topRadius,height:GLdouble;slices,stacks:GLint);
  procedure setToCylinder(radius,height:GLdouble;slices,stacks:GLint);
end;

{ TOGLRenderer }

//Verwaltet den Grafikprozess
TOGLRenderer=class
  procedure init();
  procedure resize(w,h:longint);
  procedure draw(obj: TOGLObject);
  procedure drawAt(obj: TOGLObject;x, y, z: GLfloat);
end;
var next_display_list: longint=1;
implementation

{ TOGLRenderer }
//Initialisiert Standardeinstellungen (siehe www.delphigl.com oder ein
//OpenGL Handbuch)
procedure TOGLRenderer.init();
begin
  glShadeModel(GL_SMOOTH);
  glClearColor(0.0, 0.0, 0.0, 0.5);
  glClearDepth(1.0);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
end;

//Passt die Perspektive an die Ausgabefläche an
procedure TOGLRenderer.resize(w, h: longint);
begin
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45,w/h,0.1,100);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
end;

//Zeichnet ein Objekt
procedure TOGLRenderer.draw(obj: TOGLObject);
begin
  glCallList(obj.list);
end;

//Zeichnet ein Objekt an einer bestimmten Position
procedure TOGLRenderer.drawAt(obj: TOGLObject; x, y, z: GLfloat);
begin
  glTranslatef(x,y,z);
  draw(obj);
  glTranslatef(-x,-y,-z);
end;


{ TOGLGLUObject }

constructor TOGLGLUObject.create();
begin
  inherited;
  quadric:=gluNewQuadric();
  list:=next_display_list;
  inc(next_display_list);
end;

destructor TOGLGLUObject.destroy;
begin
  gluDeleteQuadric(quadric);
  inherited;
end;

//Lässt eine Kugel berechnen
procedure TOGLGLUObject.setToSphere(radius: GLdouble; slices, stacks: GLint);
begin
  glNewList(list,GL_COMPILE);
  glColor4fv(@_color);
  gluSphere(quadric,radius,slices,stacks);
  glEndList;
end;

//Lässt einen Zylinder berechnen
procedure TOGLGLUObject.setToCylinder(baseRadius, topRadius, height: GLdouble;
  slices, stacks: GLint);
begin
  glNewList(list,GL_COMPILE);
  glColor4fv(@_color);
  gluCylinder(quadric,baseRadius,topRadius,height,slices,stacks);
  glEndList;
end;

procedure TOGLGLUObject.setToCylinder(radius, height: GLdouble; slices,
  stacks: GLint);
begin
  setToCylinder(radius,radius,height,slices,stacks);
end;

{ TOGLObject }

constructor TOGLObject.create();
begin
  color[0]:=1;
  color[1]:=1;
  color[2]:=1;
  color[3]:=1;
end;

procedure TOGLObject.setColor(c: T4fArray);
begin
  _color:=c;
end;

//Delphifarbe zu OGL-Farbe
procedure TOGLObject.setColor(c: TColor);
const R_MASK=$0000FF;
      G_MASK=$00FF00;
      B_MASK=$FF0000;
begin
  _color[0]:=(c and R_MASK)/R_MASK;
  _color[1]:=(c and G_MASK)/G_MASK;
  _color[2]:=(c and B_MASK)/B_MASK;
end;

end.

