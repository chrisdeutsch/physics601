#include "math.h"
#include "TMath.h"

double ElecCalib(double e_raw, double pt, double eta, double phi,
    double etiso, double eoverp, double mindrjet)
{
    double energy = e_raw;
    
    // Erste Iteration: Azimuthalwinkel
    if (phi >= 0. && phi < 0.4) energy *= 91.19 / 89.81;
    if (phi >= 4. && phi < 0.8) energy *= 91.19 / 89.98;
    if (phi >= 0.8 && phi < 1.2) energy *= 91.19 / 89.86;
    if (phi >= 1.2 && phi < 1.6) energy *= 91.19 / 90.29;
    if (phi >= 1.6 && phi < 2.0) energy *= 91.19 / 89.76;
    if (phi >= 2.0 && phi < 2.4) energy *= 91.19 / 89.68;
    if (phi >= 2.4 && phi < 2.8) energy *= 91.19 / 90.07;
    if (phi >= 2.8 && phi < 3.2) energy *= 91.19 / 90.15;
    
    // Erste Iteration: Energie und Pseudorapiditaet
    if (energy < 50.0) {
        if (fabs(eta) >= 0.0 && fabs(eta) < 0.5)
            energy *= 91.19 / 88.83;
        if (fabs(eta) >= 0.5 && fabs(eta) < 1.0)
            energy *= 91.19 / 86.85;
        if (fabs(eta) >= 1.0 && fabs(eta) < 1.5)
            energy *= 91.19 / 86.15;
        if (fabs(eta) >= 1.5 && fabs(eta) < 2.0)
            energy *= 91.19 / 83.88;
    } else if (energy < 100.0) {
        if (fabs(eta) >= 0.0 && fabs(eta) < 0.5)
            energy *= 91.19 / 93.16;
        if (fabs(eta) >= 0.5 && fabs(eta) < 1.0)
            energy *= 91.19 / 90.07;
        if (fabs(eta) >= 1.0 && fabs(eta) < 1.5)
            energy *= 91.19 / 86.6;
        if (fabs(eta) >= 1.5 && fabs(eta) < 2.0)
            energy *= 91.19 / 85.48;
        if (fabs(eta) >= 2.0 && fabs(eta) < 2.5)
            energy *= 91.19 / 83.2;
    } else {
        if (fabs(eta) >= 0.0 && fabs(eta) < 1.5)
            energy *= 91.19 / 88.26;
        if (fabs(eta) >= 1.5 && fabs(eta) < 1.75)
            energy *= 91.19 / 82.54;
        if (fabs(eta) >= 1.75 && fabs(eta) < 2.0)
            energy *= 91.19 / 81.76;
        if (fabs(eta) >= 2.0 && fabs(eta) < 2.5)
            energy *= 91.19 / 78.87;
    }
    
    // Zweite Iteration: Azimuthalwinkel
    if (phi >= 0. && phi < 0.4) energy *= 91.19 / 90.35;
    if (phi >= 4. && phi < 0.8) energy *= 91.19 / 89.71;
    if (phi >= 0.8 && phi < 1.2) energy *= 91.19 / 90.34;
    if (phi >= 1.2 && phi < 1.6) energy *= 91.19 / 90.18;
    if (phi >= 1.6 && phi < 2.0) energy *= 91.19 / 90.30;
    if (phi >= 2.0 && phi < 2.4) energy *= 91.19 / 90.29;
    if (phi >= 2.4 && phi < 2.8) energy *= 91.19 / 90.27;
    if (phi >= 2.8 && phi < 3.2) energy *= 91.19 / 90.35;
    
    // Zweite Iteration: Pseudorapiditaet
    if (fabs(eta) >= 0.0 && fabs(eta) < 1.0)
        energy *= 91.19 / 92.19;
    if (fabs(eta) >= 1.0 && fabs(eta) < 1.5)
        energy *= 91.19 / 90.54;
    if (fabs(eta) >= 1.5 && fabs(eta) < 2.0)
        energy *= 91.19 / 88.12;
    if (fabs(eta) >= 2.0 && fabs(eta) < 2.5)
        energy *= 91.19 / 86.56;
    
    // Zweite Iteration: Energie
    if (energy >= 0.0 && energy < 50.0) energy *= 91.19 / 91.81;
    if (energy >= 50.0 && energy < 100.0) energy *= 91.19 / 91.04;
    if (energy >= 100.0 && energy < 150.0) energy *= 91.19 / 90.26;
    if (energy >= 150.0) energy *= 91.19 / 89.16;
    
    // Dritte Iteration: Pseudorapiditaet
    if (fabs(eta) >= 0.0 && fabs(eta) < 0.75)
        energy *= 91.19 / 92.06;
    if (fabs(eta) >= 0.75 && fabs(eta) < 1.5)
        energy *= 91.19 / 91.16;
    if (fabs(eta) >= 1.5 && fabs(eta) < 2.0)
        energy *= 91.19 / 90.46;
    if (fabs(eta) >= 2.0 && fabs(eta) < 2.5)
        energy *= 91.19 / 90.03;
    
    // Abschliessende Skalierung aller Elektronen
    energy *= 91.19 / 91.28;
    
    return energy;
}
