{Diese Datei enthält die Datenstruktur für ein Atom und elementare Information
 über diese
}
unit atoms;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,math,extmath;

{Atomare Einheiten (von en.wikipedia):
length 	Bohr radius 	a0 	5.291 772 108(18)*10-11 m
                             	52.91 772 108(18)*10-12 m
                                0.5291 772 108 Â
mass 	electron rest mass 	me 	9.109 3826(16)*10-31 kg
charge 	elementary charge 	e 	1.602 176 53(14)*10-19 C
}
type
  TAtom=class;
  TBond=record               //Eine Bindung
    a:TAtom;                 //zu einem Atom
    count: longint;          //Bindungszahl
  end;
  TBonds=array[1..4] of tbond;
  TAtom=class
  public
    index: longint;          //Index in TList
    typ: longint;            //Ordnungszahl
    bond:TBonds;             //kovalente bindungen
    james:longint;           //Anzahl Bindungen
    p: TVector3n;            //Raum Position

    flatX,flatY:  longint;   //2D Position (Lewis-Formel)

    firstOrbital: longint;   //Index des ersten Orbitals
  end;

const AT_H=1;                //Abkürzungen für Atome
      AT_C=6;
      AT_N=7;
      AT_O=8;
      //Abkürzungen als Strings
      ATOM_ABBREV:array[1..8]of string=('H','','','','','C','N','O');
      //Ausgeschriebene Namen
      ATOM_NAMES:array[1..8]of string=('Wasserstoff','','','','','Kohlenstoff',
                                       'Stickstoff','Sauerstoff');
      //Mögliche bindende Elektronenpaare
      ATOM_CLOSED_BONDS: array[1..8] of longint=(1,0,0,0,0,4,3,2);
      //Elektronegativität (Quelle: http://www.webelements.com/)
      EN:array[1..8] of number=(2.2,nan,nan,nan,nan,2.55,3.04,3.44);
      //Kovalenzradius (Quelle: de.wikipedia.org und http://www.webelements.com/)
      CovalentRadius:array[1..8] of number=(0.699,nan,nan,nan,nan,1.455,1.417,1.380);

implementation

end.

