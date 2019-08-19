// -----------------------------------------
// PrtPixelSD.h
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtPixelSD_h
#define PrtPixelSD_h 1

#include "G4VSensitiveDetector.hh"
#include "TVector3.h"
#include <vector>
#include "TH2F.h"
class G4Step;
class G4HCofThisEvent;

/// Calorimeter sensitive detector class
///
/// In Initialize(), it creates one hit for each calorimeter layer and one more
/// hit for accounting the total quantities in all layers.
///
/// The values are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step.

class PrtPixelSD : public G4VSensitiveDetector
{
  public:
    PrtPixelSD(const G4String& name, 
                     const G4String& hitsCollectionName, 
                     G4int nofCells);
    virtual ~PrtPixelSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
  G4int     fNofCells;
  G4double fQe_space[15][64];
  G4int fMultHit[15][64];
    TH2F* hist_time_angle_ch[12][64];



TVector3 momInBar;
 Double_t angle ;
 int kt, ka, nEntries;
};

#endif
