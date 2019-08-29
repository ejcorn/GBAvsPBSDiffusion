#!/bin/bash
set -euxo pipefail

echo 'type output directory name:'
read ROOT

cd ~/Dropbox/Neurodegeneration/GBAvsPBSDiffusion
SAVEDIR=ResultsForMike8219
mkdir -p $SAVEDIR

# copy predicted path and vulnerability csv files
cp $ROOT/diffmodel/*.csv $SAVEDIR/
# delete the ones that are for me, not in the order Mike wants them
rm $SAVEDIR/*Eli.csv

cp $ROOT/diffmodel/DPBSpredictedpathcontinuous.pdf $SAVEDIR/Fig1a.pdf
cp $ROOT/diffmodel/DPBSspread.log $SAVEDIR/Fig1a.log
cp $ROOT/diffmodel/CBEpredictedpathcontinuous.pdf $SAVEDIR/Fig1b.pdf
cp $ROOT/diffmodel/CBEspread.log $SAVEDIR/Fig1b.log

cp $ROOT/diffmodel/seedspec/DPBSiCPSeedSpecificity.pdf $SAVEDIR/Fig2a.pdf
cp $ROOT/diffmodel/seedspec/CBEiCPSeedSpecificity.pdf $SAVEDIR/Fig2b.pdf
cp $ROOT/diffmodel/seedspec/*.csv $SAVEDIR/

cp $ROOT/diffmodel/seedspec/DPBSAlternateSeedFitByInProjectionSimilarity.pdf $SAVEDIR/Fig2c.pdf
cp $ROOT/diffmodel/seedspec/DPBSconnsimaltseed.log $SAVEDIR/Fig2c.log
cp $ROOT/diffmodel/seedspec/CBEAlternateSeedFitByInProjectionSimilarity.pdf $SAVEDIR/Fig2d.pdf
cp $ROOT/diffmodel/seedspec/CBEconnsimaltseed.log $SAVEDIR/Fig2d.log

#cp $ROOT/PBSvsCBE/BaseTCjitterTF0.5.pdf $SAVEDIR/Fig3a.pdf