#include "TLorentzVector.h"


void testLorentz()
{
TLorentzVector g1P4;

g1P4.SetPtEtaPhiE(4.0,0.0,0.5,5.0);

cout<<g1P4.M()<<endl;


}
