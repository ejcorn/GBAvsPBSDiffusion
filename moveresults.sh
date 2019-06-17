#!/bin/bash
set -euxo pipefail

echo 'type output directory name:'
read ROOT

cd ~/Dropbox/Neurodegeneration/GBAvsPBSDiffusion
mkdir -p ResultsForMike

# copy predicted path and vulnerability csv files
cp $ROOT/diffmodel/*.csv ResultsForMike/
# delete the ones that are for me, not in the order Mike wants them
rm ResultsForMike/*Eli.csv

cp $ROOT/diffmodel/CBEpredictedpathcontinuous.pdf ResultsForMike/Fig1a.pdf
cp $ROOT/diffmodel/CBEspread.log ResultsForMike/Fig1a.log
cp $ROOT/diffmodel/DPBSpredictedpathcontinuous.pdf ResultsForMike/Fig1b.pdf
cp $ROOT/diffmodel/DPBSspread.log ResultsForMike/Fig1b.log

cp $ROOT/diffmodel/seedspec/CBEiCPSeedSpecificity.pdf ResultsForMike/Fig2a.pdf
cp $ROOT/diffmodel/seedspec/DPBSiCPSeedSpecificity.pdf ResultsForMike/Fig2b.pdf

cp $ROOT/diffmodel/seedspec/CBEAlternateSeedFitByInProjectionSimilarity.pdf ResultsForMike/Fig2c.pdf
cp $ROOT/diffmodel/seedspec/CBEconnsimaltseed.log ResultsForMike/Fig2c.log
cp $ROOT/diffmodel/seedspec/DPBSAlternateSeedFitByInProjectionSimilarity.pdf ResultsForMike/Fig2d.pdf
cp $ROOT/diffmodel/seedspec/DPBSconnsimaltseed.log ResultsForMike/Fig2d.log

#cp $ROOT/PBSvsCBE/BaseTCjitterTF0.5.pdf ResultsForMike/Fig3a.pdf