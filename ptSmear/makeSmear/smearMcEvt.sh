#!/bin/bash

echo "simulating CohJpsi"
root -l -b -q 'smearMcEvt.C+("CohJpsi")'

echo "simulating CohPsi2S"
root -l -b -q 'smearMcEvt.C+("CohPsi2S")'

echo "simulating LowMassGammaGamma"
root -l -b -q 'smearMcEvt.C+("LowMassGammaGamma")'
