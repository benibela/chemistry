unit energyCalculator;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,extmath,atoms,molecules,stosim;

type
  TEnergyCalculator = class
    procedure initialize(mol:TMolecule); virtual;abstract;
    function calculate(onlyThisAtomMoved: longint=-1):number; virtual;abstract;
  end;

  { TSCFEnergyCalculator }

  TSCFEnergyCalculator=class(TEnergyCalculator)
  private
    currentMolecule: TMolecule;
  
    orbitals:array of TSTO;               //vorhandenen AO-Orbitale
    elektrons: longint;                   //Anzahl Elektronen
    //Temporäre Speicher für die SCF Energieberechnung
    _overlapSymMat:TSymMatrix;            //Überlappungsintegrale
    _orthoOverlapSymMat:TSymMatrix;       //  " hoch -0.5
    _oneElectronEnergySymMat:TSymMatrix;  //kinetische+Kernenergie
    _coefMat:TMatrix;            //Koeffizienten i-te MO = Summe (c * j-te AO) T
    _denseSymMat:TSymMatrix;              //Dichteverteilung
    _fockSymMat:TSymMatrix;               //Fockmatrix
    _tempMat1,_tempMat2: TMatrix;         //temporärer
    _repulsion: array of array of TMatrix;//Elektronenwechselwirkungen

    //Sicherungen der obrigen Matrizen
    _saveOverlap:TMatrix;
    _saveRepulsion:array of array of TMatrix;

    _energies,_tempVecN: tvectorn;  //MO-Energie und temporär
    _tempVecI: TVectori;            //                   "
    _tempVecB: TVectorb;            //                   "
    _lastEnergies:array of tvectorn;//Die letzten bisher aufgetreten MO-Energien
    _sortpos: array of longint;     //für Sortieralgorithmus
    _sorted: array of boolean;      // "    "

    //Rät die Koeffizientenmatrix
    procedure guessCoefficients;

  public
    procedure initialize(mol:TMolecule); override;
    function calculate(onlyThisAtomMoved: longint=-1):number; override;
  end;

implementation

uses math;

{ ==============================================================================
                           Self Consistent Field Calculation
  ==============================================================================
}

procedure TSCFEnergyCalculator.guessCoefficients;
var pos: longint; //Anzahl der bereits festgelegten Werten
  procedure addAOsToMO(const atom: tatom);inline;
  begin with atom do
    case typ of
      AT_H: _coefMat[pos][firstOrbital]:=1; //H Atome geben 1s in die Bindung
      else begin
        //Alle anderen Atome werden als sp3-Hybridisierung probiert
        //(sollte es nicht stimmen, ändert das SCF Verfahren dies).
        _coefMat[pos][firstOrbital+1]:=1;
        _coefMat[pos][firstOrbital+2]:=1;
        _coefMat[pos][firstOrbital+3]:=1;
        _coefMat[pos][firstOrbital+4]:=1;
      end;
    end;
  end;
var i,j,k: longint;
begin
  with currentMolecule do begin
    pos:=0;
    //Ein MO für jede Bindung
    for i:=0 to atomCount - 1 do
      for j:=1 to atoms[i].james do
        if atoms[i].bond[j].a.index>atoms[i].index then
          for k:=1 to atoms[i].bond[j].count do begin
            materaseLine(_coefMat,pos);
            addAOsToMO(atoms[i]);
            addAOsToMO(atoms[i].bond[j].a);
            pos+=1;
          end;
    //Übriggebliebene MO besetzten:je 1s für C,N,O + 1 MO für N und 2 MOs für O
    for i:=0 to atomCount - 1 do begin
      if atoms[i].typ=AT_H then continue;
      materaseLine(_coefMat,pos);
      _coefMat[pos][atoms[i].firstOrbital]:=1;
      pos+=1;
      case atoms[i].typ of
        AT_N: begin
          //Probiere 2s
          materaseLine(_coefMat,pos);
          _coefMat[pos][atoms[i].firstOrbital+1]:=1;
          pos+=1;
        end;
        AT_O: begin
          //Probiere 2s
          materaseLine(_coefMat,pos);
          _coefMat[pos][atoms[i].firstOrbital+1]:=1;
          pos+=1;
          //und 2p-Verschmelzung
          materaseLine(_coefMat,pos);
          _coefMat[pos][atoms[i].firstOrbital+2]:=1;
          _coefMat[pos][atoms[i].firstOrbital+3]:=1;
          _coefMat[pos][atoms[i].firstOrbital+4]:=1;
          pos+=1;
        end;
      end;
    end;
    assert(2*pos=elektrons,'Geschätze MO-Anzahl ungültig: '+inttostr(pos));
  end;
end;

//Initialisiert das SCF Verfahren (Orbitalberechnung und Speicherreservierung)
procedure TSCFEnergyCalculator.initialize(mol: TMolecule);
var i,j:longint;
begin
  currentMolecule:=mol;
  with currentMolecule do begin
    stosim.initUnit;;
    //Orbitalberechnung über Auffruf der Funktionen von stosim
    SetLength(orbitals,0);
    elektrons:=0;
    for i:=0 to atomCount-1 do
      with atoms[i] do begin
        firstOrbital:=length(orbitals);
        setlength(orbitals,length(orbitals)+OrbitalCount[typ]);
        for j:=1 to OrbitalCount[typ] do
          orbitals[firstOrbital+j-1]:=createSTO(atoms[i],j);
        elektrons+=typ;
      end;
    //Speicherplatzreservierung
    setlength(_overlapSymMat,length(orbitals),length(orbitals));
    setlength(_orthoOverlapSymMat,length(orbitals),length(orbitals));
    setlength(_oneElectronEnergySymMat,length(orbitals),length(orbitals));
    setlength(_coefMat,length(orbitals),length(orbitals));
    setlength(_denseSymMat,length(orbitals),length(orbitals));
    setlength(_fockSymMat,length(orbitals),length(orbitals));
    setlength(_tempMat1,length(orbitals),length(orbitals));
    setlength(_tempMat2,length(orbitals),length(orbitals));
    setlength(_repulsion,length(orbitals),length(orbitals),
                         length(orbitals),length(orbitals));
    setlength(_energies,length(orbitals));
    setlength(_tempVecN,length(orbitals));
    setlength(_tempVecI,length(orbitals));
    setlength(_tempVecB,length(orbitals));
    setlength(_lastEnergies,3,length(orbitals));
    setlength(_sortpos,length(orbitals));
    setlength(_sorted,length(orbitals));
    SetLength(_saveOverlap,length(orbitals),length(orbitals));
    SetLength(_saveRepulsion,length(orbitals),length(orbitals),
                             length(orbitals),length(orbitals));
  end;
end;

function TSCFEnergyCalculator.calculate(onlyThisAtomMoved: longint=-1): number;
const PREC=1e-10;                          //Zielpräzision
var
  i,j,k,l:longint;                         //Schleifenvariablen
  energyDiff,energyGain,tempNumber: number;//Energiedifferenz zum vorherigen
  tempVecPointer: TVectorn;                //Zeiger auf Energiearray
  oneAtomFirstOrbital,oneAtomLastOrbital:longint;//Orbitalgrenzen vom Atom
begin
 with currentMolecule do begin
  //=========================Berechnen der Integrale=========================
  if onlyThisAtomMoved>=0 then begin
    //Es wurde ein Atom angegeben, dass als einzigstes verändert wurde
    //Orbitale des Atoms suchen
    oneAtomFirstOrbital:=atoms[onlyThisAtomMoved].firstOrbital;
    oneAtomLastOrbital:=oneAtomFirstOrbital+
                        OrbitalCount[atoms[onlyThisAtomMoved].typ]-1;
    //Alle Über. und Elek. Integrale berechnen, in denen das Atom vorkommt
    //Überlappungsmatrix
    //wiederherstellen
    for i:=0 to high(orbitals) do
      for j:=i to high(orbitals) do
        _overlapSymMat[i,j]:=_saveOverlap[i,j];
    //berechnen
    for i:=0 to oneAtomLastOrbital do
      for j:=i to high(orbitals) do
        if (i>=oneAtomFirstOrbital) or
          ((j>=oneAtomFirstOrbital)and(j<=oneAtomLastOrbital)) then
        _overlapSymMat[i,j]:=overlapIntegral(orbitals[i],orbitals[j]);
    //Elektronenwechselwirkungen
    for i:=0 to oneAtomLastOrbital do
     for j:=i to high(orbitals) do
      for k:=i to high(orbitals) do
       for l:=k to high(orbitals) do
        if (i>=oneAtomFirstOrbital) or
           ((j>=oneAtomFirstOrbital) and (j<=oneAtomLastOrbital)) or
           ((k>=oneAtomFirstOrbital) and (k<=oneAtomLastOrbital)) or
           ((l>=oneAtomFirstOrbital) and (l<=oneAtomLastOrbital)) then begin
             _saveRepulsion[i,j,k,l]:=_repulsion[i,j,k,l];
             _repulsion[i,j,k,l]:=electronRepulsionIntegral(orbitals[i],
                                      orbitals[j],orbitals[k],orbitals[l]);
        end;
  end else begin
    //Alle Atome wurden geändert (oder sind unbekannt)
    //Überlappungsintegrale
    for i:=0 to high(orbitals) do
      for j:=i to high(orbitals) do begin
        _overlapSymMat[i,j]:=overlapIntegral(orbitals[i],orbitals[j]);
        _saveOverlap[i,j]:=_overlapSymMat[i,j]; //sichern
      end;


    //Elektronenwechselwirkungsintegrale
    for i:=0 to high(_repulsion) do
     for j:=i to high(_repulsion[i]) do
      for k:=i to high(_repulsion[i,j]) do
       for l:=k to high(_repulsion[i,j,k]) do
         _repulsion[i,j,k,l]:=electronRepulsionIntegral(orbitals[i],orbitals[j],
                                                      orbitals[k],orbitals[l]);

    //Koeffizienten raten
    guessCoefficients;
  end;

  //Kinetikintegrale immer berechnen
  for i:=0 to high(orbitals) do
    for j:=i to high(orbitals) do
      _oneElectronEnergySymMat[i,j]:=kineticEnergyIntegral(orbitals[i],
                                                           orbitals[j]);


  //Kernanziehungsintegrale immer berechnen
  for k:=0 to atomCount-1 do
    for i:=0 to high(orbitals) do
      for j:=i to high(orbitals) do
        _oneElectronEnergySymMat[i,j]-=
          nuclearAttractionIntegral(orbitals[i],orbitals[j],atoms[k]);



  //orthoOverlap = overlap^-0.5
  //Achtung: overlap wird dabei zerstört (deshalb wurde es gesichert)
  matinvsqrt_sym(_overlapSymMat,_tempMat1,_tempMat2,_tempVecB,_tempVecI, _tempVecN,
                 _orthoOverlapSymMat);


  //==========================Optimierung der MOs===============================
  //Löschen der letzten Energien
  fillchar(_energies[0],sizeof(_energies[0])*length(_energies),0);
  for i:=0 to high(_lastEnergies) do
    fillchar(_lastEnergies[i,0],
             sizeof(_lastEnergies[i,0])*length(_lastEnergies[i]),0);
  repeat
    //Berechnung der Dichtematrix
    for i:=0 to high(orbitals) do
      for j:=i to high(orbitals) do begin
        _denseSymMat[i,j]:=0;
        for k:=0 to elektrons div 2-1 do
          _denseSymMat[i,j]+=_coefMat[k,i]*_coefMat[k,j];
      end;

    //Berechnung der Fockmatrix
    for i:=0 to high(orbitals) do
      for j:=i to high(orbitals) do begin
        _fockSymMat[i,j]:=_oneElectronEnergySymMat[i,j];
        for k:=0 to i do begin
          for l:=k+1 to i do
            _fockSymMat[i,j]+=_denseSymMat[k,l]*(4*_repulsion[k,l,i,j]-
                                       _repulsion[k,i,l,j]-_repulsion[k,j,l,i]);
          for l:=i+1 to j do
            _fockSymMat[i,j]+=_denseSymMat[k,l]*(4*_repulsion[k,l,i,j]-
                                       _repulsion[k,i,l,j]-_repulsion[k,j,i,l]);
          for l:=j+1 to high(orbitals) do
            _fockSymMat[i,j]+=_denseSymMat[k,l]*(4*_repulsion[k,l,i,j]-
                                       _repulsion[k,i,j,l]-_repulsion[k,j,i,l]);
          _fockSymMat[i,j]+=_denseSymMat[k,k]*(2*_repulsion[k,k,i,j]-
                                               _repulsion[k,i,k,j]);
        end;
        for k:=i+1 to j do begin
          for l:=k+1 to j do
            _fockSymMat[i,j]+=_denseSymMat[k,l]*(4*_repulsion[i,j,k,l]-_repulsion[i,k,l,j]-_repulsion[i,l,k,j]);
          for l:=j+1 to high(orbitals) do
            _fockSymMat[i,j]+=_denseSymMat[k,l]*(4*_repulsion[i,j,k,l]-_repulsion[i,k,j,l]-_repulsion[i,l,k,j]);
          _fockSymMat[i,j]+=_denseSymMat[k,k]*(2*_repulsion[i,j,k,k]-_repulsion[i,k,k,j]);
        end;
        for k:=j+1 to high(orbitals) do begin
          for l:=k+1 to high(orbitals) do
            _fockSymMat[i,j]+=_denseSymMat[k,l]*(4*_repulsion[i,j,k,l]-
                              _repulsion[i,k,j,l]-_repulsion[i,l,j,k]);
          _fockSymMat[i,j]+=_denseSymMat[k,k]*(2*_repulsion[i,j,k,k]-
                                              _repulsion[i,k,j,k]);
        end;
    end;

    //Fockmatrix in eine ganze Matrix umwandeln zur Multiplikation
    //mit der OrthoOverlapMatrix (Ergebnis ist nicht symmetrisch)
    matsym2full(_fockSymMat);
    //Fockmatrix mal S^-0.5 in tempMat speichern
    matmul(_fockSymMat,_orthoOverlapSymMat,_tempMat1);
    //S^-0.5 mal tempMat in fock speichern
    matmul(_orthoOverlapSymMat,_tempMat1,_fockSymMat);

    //Berechnung der Eigenvektoren der Fockmatrix
    diag_sym(_fockSymMat,_tempVecB,_tempVecI,_tempVecN,_coefMat);
    //Verformung der Eigenverte mit S^-0.5 für endgültig richtige Werte
    matmul(_coefMat,_orthoOverlapSymMat,_tempMat1);

    //Sortieren, da nur die kleinsten elektrons/2 zählen
    for i:=0 to high(_tempVecN) do
      _sortpos[i]:=i;
    for i:=1 to high(_tempVecN) do
      for j:=i-1 downto 0 do
        if _tempVecN[_sortpos[j]]>_tempVecN[_sortpos[j+1]] then begin
          k:=_sortpos[j];
          _sortpos[j]:=_sortpos[j+1];
          _sortpos[j+1]:=k;
        end;
    for i:=0 to elektrons div 2-1 do begin
      _energies[i]:=_tempVecN[_sortpos[i]];
      move(_tempMat1[_sortpos[i],0],_coefMat[i,0],
           length(_coefMat[i])*sizeof(_coefMat[i,0]));
    end;

    //Geringste Energiedifferenz zu den letzten suchen
    energyDiff:=9999999999;
    for i:=0 to high(_lastEnergies) do
      if _lastEnergies[i]<>nil then begin
        energyDiff:=0;
        for j:=0 to elektrons div 2-1 do
          energyDiff+=abs(_energies[j]-_lastEnergies[i][j]);
        if energyDiff<PREC then break;
      end;

    //lastenergies kopieren und sortieren
    tempVecPointer:=_lastEnergies[0];
    for i:=0 to high(_lastEnergies)-1 do
      _lastEnergies[i]:=_lastEnergies[i+1];
    _lastEnergies[high(_lastEnergies)]:=tempVecPointer;
    move(_energies[0],_lastEnergies[high(_lastenergies)][0],
         sizeof(_energies[0])*length(_energies));
  until energyDiff<PREC;

  //Elektronen-Energie berechnen
  result:=0;
  for i:=0 to high(orbitals) do
    for j:=i to high(orbitals) do begin
      energyGain:=_oneElectronEnergySymMat[i,j];
      for k:=0 to high(orbitals) do
        for l:=0 to high(orbitals) do begin
          if min(k,l)<i then tempNumber:=_repulsion[min(k,l),max(l,k),i,j]
          else tempNumber:=_repulsion[i,j,min(k,l),max(l,k)];
          if min(i,k)<min(j,l) then
            tempNumber-=0.5*_repulsion[min(i,k),max(i,k),min(j,l),max(j,l)]
          else tempNumber-=0.5*_repulsion[min(j,l),max(j,l),min(i,k),max(i,k)];
          energyGain+=_denseSymMat[min(k,l),max(k,l)]*tempNumber;
        end;
      if i<>j then energyGain*=2;
      result+=2*_denseSymMat[i,j]*energyGain;
    end;

  //Kern-Kern-Energie berechnen
  for i:=0 to atomCount-1 do
    for j:=0 to i-1 do
      result+=atoms[i].typ*atoms[j].typ/vecdist(atoms[i].p,atoms[j].p);

  //Wiederherstellen der Sicherung
  if onlyThisAtomMoved>=0 then
    for i:=0 to oneAtomLastOrbital do
     for j:=i to high(orbitals) do
      for k:=i to high(orbitals) do
       for l:=k to high(orbitals) do
        if (i>=oneAtomFirstOrbital) or
           ((j>=oneAtomFirstOrbital) and (j<=oneAtomLastOrbital)) or
           ((k>=oneAtomFirstOrbital) and (k<=oneAtomLastOrbital)) or
           ((l>=oneAtomFirstOrbital) and (l<=oneAtomLastOrbital)) then
             _repulsion[i,j,k,l]:=_saveRepulsion[i,j,k,l];
 end;
end;

end.

