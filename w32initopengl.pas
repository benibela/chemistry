{Benito van der Zander
}
unit w32initOpengl;

{$mode objfpc}{$H+}

interface

uses
  windows,gl,sysutils;

type

{ TOGLW32Base }
EOGLW32BaseException=class(Exception);
TOGLW32Base=class
private
  wnd: hwnd; //Fenster
  dc: HDC;   //Device-Context
  gl:HGLRC;  //OpenGL-Context
public
  property ogl: HGLRC read gl;
  constructor create(const awnd: hwnd;colorBits:longint=32;depthBits:longint=16;doublebuffer:boolean=true);
  destructor destroy;override;
  
  procedure endScene;
end;

implementation

{ TOGLW32Base }

constructor TOGLW32Base.create(const awnd: hwnd;
                               colorBits:longint=32;depthBits:longint=16;
                               doublebuffer:boolean=true);
var pfd : TPIXELFORMATDESCRIPTOR;
    pfm : GLuint;
begin
  wnd:=awnd;
  dc:=getDC(wnd);
  if (dc = 0) then raise EOGLW32BaseException.create('DC nicht gekriegt');
  FillChar(pfd,sizeof(pfd),0);
  with pfd do
  begin
    nSize      := SizeOf(TPIXELFORMATDESCRIPTOR);
    nVersion   := 1;
    dwFlags    := PFD_DRAW_TO_WINDOW or PFD_SUPPORT_OPENGL;
    if doublebuffer then
      dwFlags  := dwFlags or PFD_DOUBLEBUFFER;
    iPixelType := PFD_TYPE_RGBA;
    cColorBits := colorBits;
    cDepthBits := depthBits;
    iLayerType := PFD_MAIN_PLANE;
  end;

  pfm := ChoosePixelFormat(dc, @pfd);
  if (pfm = 0) then
    raise EOGLW32BaseException.create('PixelFormat nicht gefunden');

  if (not SetPixelFormat(dc, pfm, @pfd)) then
    raise EOGLW32BaseException.create('PixelFormat wird nicht unterstützt');

  gl := wglCreateContext(dc);
  if (gl = 0) then
    raise EOGLW32BaseException.create('OpenGL konnte nicht geladen werden');

  if (not wglMakeCurrent(dc,gl)) then
    raise EOGLW32BaseException.create('OpenGL konnte nicht aktiviert werden');
end;

destructor TOGLW32Base.destroy;
begin
  wglMakeCurrent(dc, 0);
  wglDeleteContext(gl);
  ReleaseDC(wnd, dc);
  inherited;
end;

procedure TOGLW32Base.endScene;
begin
  SwapBuffers(dc);
end;

end.

