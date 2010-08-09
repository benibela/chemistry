program folding;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Interfaces, // this includes the LCL widgetset
  Forms
  { add your units here },  molecules, w32initOpengl, atomicDrawing,
  extmath, gtf, ogldrawing, stosim, atoms, lewisRenderer, geometryCalculator,
energyCalculator,gui;

begin
  Application.Initialize;
  Application.CreateForm(TForm1, Form1);
  Application.Run;
end.

