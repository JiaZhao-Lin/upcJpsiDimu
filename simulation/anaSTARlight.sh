#!/bin/bash

date

nAnaEvts=10000000

echo "working on CohJpsi"
root -l -b -q 'anaSTARlight.C+("CohJpsi", '${nAnaEvts}')'
echo "working on CohJpsi 0n0n"

root -l -b -q 'anaSTARlight.C+("CohJpsi_0n0n", '${nAnaEvts}')'
echo "working on CohJpsi 0nXn"
root -l -b -q 'anaSTARlight.C+("CohJpsi_0nXn", '${nAnaEvts}')'
echo "working on CohJpsi XnXn"
root -l -b -q 'anaSTARlight.C+("CohJpsi_XnXn", '${nAnaEvts}')'
echo "working on InCohJpsi"
root -l -b -q 'anaSTARlight.C+("InCohJpsi", '${nAnaEvts}')'
echo "working on CohPsi2S"
root -l -b -q 'anaSTARlight.C+("CohPsi2S", '${nAnaEvts}')'
echo "working on InCohPsi2S"
root -l -b -q 'anaSTARlight.C+("InCohPsi2S", '${nAnaEvts}')'
