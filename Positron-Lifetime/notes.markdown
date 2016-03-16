# Positron Lifetime in Metals and Insulators


## Required Knowledge
- metals
    - lattice defects in metals
    - trapping model, positron lifetime spectrum
    - vacancy concentration in thermal equilibrium, formation enthalpy of vacancies
    
- positrons
    - positron sources, decay scheme of ${}^{22}\mathrm{Na}$ ($\gamma$-energy spectrum)
    - positron annihilation
    - positronium formation, pickoff-decay
    
- measurement technique
    - LYSO scintillator, ${}^{176}\mathrm{Lu}$ decay scheme, LYSO spectrum
    - photomultiplier (PM), dynode pickup
    - single channel analyzer (SCA), constant-fraction-discriminator (CFD)
    - time-to-amplitude converter (TAC), multichannel-analyzer (MCA)
    - fast-slow coincidence curcuit
    - prompt curve, resolving power, detector response function

## Measurement technique

### Szintillator & Gamma-Spektroskopie

![Gamma-Spektrum](https://upload.wikimedia.org/wikipedia/commons/f/f2/Am-Be-SourceSpectrum.jpg)

Szintillatoren sind zeitlich schnell, haben hohe Effizienz dafür eine geringe Energieauflösung.

- **Photopeak:** Absorption des Photons durch Photoeffekt im Szintillator
    - gesamte Energie deponiert
- **Compton-Kontinuum:** Comptonstreuung des Photons im Szintillator
    - Energiedeposition kleiner als gesamte Photonenenergie (abhängig vom Streuwinkel)
- **Escape-Peaks:** Photon macht Paarerzeugung und ein / beide Photon verlässt den Szintillator
    - Double Escape: Peak bei E - 2 * 511 keV
    - Single Escape Peak bei E - 511 keV
    - No Escape: Energie liegt im Photopeak
    
#### LYSO scintillator
- Inorganic Crystal
    - Teilchen regt Elektron aus dem Valenzband in das Leitungsband oder Exitonband an
    - Loch und Elektron wandern durch den Kristall bis sie auf eine Verunreinigung (z.B. Dotierung mit sog. Aktivatoren)im Kristall treffen
    - an der Verunreinigung existiert ein Zustand zwischen Leitungsband und Valenzband an dem Loch und Elektron rekombinieren können
- hohe Dichte und große Kernladungszahl
    - gut zur Detektion von Photonen
- hohe Lichtausbeute 32 Photonen pro keV
- kurze Abklingzeit 41 ns
- LYSO PreLude420: 81% Lutetiumoxid, 14% Silikatoxid, 5% Yttriumoxid
- Zusätzlicher Untergrund wegen dem Isotop ${}^{176}\mathrm{Lu}$
    - Zerfällt über Beta Zerfall in angeregtes Hf
        - Zerfällt über 3 $\gamma$-Quanten
        - Erzeugt mehrere Peaks je nachdem welche Gammas durch den Szintillator detektiert werden
        - z.B. 88 + 202 + 307 keV, 88 + 307 keV etc.
    - zusätzlich gibt es ein Kontinuum von dem emittierten Elektron

![Lutetium Isotop Zerfall](http://www.qsl.net/k/k0ff/Lu-176/FIG9.gif)


### Photomultipliers (PMT)

![PMT](https://upload.wikimedia.org/wikipedia/commons/5/5f/PhotoMultiplierTubeAndScintillator.jpg)

- Photoelektrischer Effekt: Elektron entsteht an der Photokathode
- Elektron wird zur Dynode beschleunigt und schlägt aus diesen zusätzliche Elektronen aus
    - exponentieller Anstieg der Elektronenzahl
    - Dynoden in einer Spannungsteilerkette an der Hochspannung (1000 - 2000 V)
- Letzte Dynode: Anode
- slow-Abgriff: langsames Abklingen des Signals (an einer Dynode)
    - Energiemessung möglich, da Strom nicht gesättigt (linear)
- fast-Abgriff: schnelles Abklingen des Signals (an Anode)
    - nicht zwangsweise proportional zur Teilchenenergie (gesättigt)
    - häufig zur Zeitmessung genutzt


### Single channel analyzer (SCA)
Auch Differential Discriminator genannt.  
Erzeugt logisches Signal wenn die Amplitude des Eingangssignals zwischen einer
unteren und oberen Schwelle liegt (Amplitudenfenster).
Kann genutzt werden um Signale mit einer bestimmten Energie zu selektieren.

- normal mode: untere und obere Schranke wird explizit gesetzt
- window mode: untere Schranke und Fensterbreite wird gesetzt
- integral mode: nur untere Grenze wird gesetzt (es gibt keine obere Grenze)


### Multichannel-analyzer (MCA)
Sortieren einkommende Pulse nach Pulshöhe und zählen die Anzahl der Pulse in
jedem Pulshöhenintervall (Histogramm der Pulshöhe). Einkommendes Signal wird mit
einem ADC in ein digitales Signal gewandelt, welches zum adressieren eines
Speicherbausteins genutzt wird. In diesem sind die Pulshöhencounts gespeichert
und können so hochgezählt werden.

In der Regel findet vor dem ADC eine Pulsformung statt um die Signalvoraussetzungen
des ADC's zu erfüllen.


### Constant-Fraction-Discriminator (CFD)
Diskriminator: Wenn ein Signal über einer Schwelle ist, wird ein logisches Signal vom Diskriminator ausgegeben. Gut um Rauschen von PMT's zu unterdrücken.

Triggermethoden:

- leading edge triggering: triggert wenn Signal über die Schwelle geht
    - koinzidente Signale mit unterschiedlichen Amplituden triggern zu unterschiedlichen Zeiten (Siehe LEO Figure 17.1) nennt man "walk"
- constant fraction triggering: Das Logiksignal wird bei einem konstanten Bruchteil der Signalamplitude erzeugt und ist so fast "walk"-free. 


### Time-to-amplitude converter (TAC)
Konvertiert die Zeit zwischen zwei Logikpulsen zu einem Ausgangssignal mit einer Amplitude die proportional zur verstrichenen Zeit ist. Kombiniert mit einem MCA kann ein Spektrum als Funktion eines Zeitintervalls aufgenommen werden.

Startpuls: Startet die Ladung/Entladung eines Kondensators mit konstantem Strom  
Stoppuls: Formt das Signal mit der Spannung am Kondensator.


### Fast-Slow coincidence circuit


### prompt curve, resolving power, detector response function


