#!/bin/sh
pwdd=`pwd`
#rm -r genrp_4.5GeV2_conf3_1.root
cd ../../g4_work/build/
g4sbs preinit_ckov_noscint_nocalorimeters.mac genrp_4.5GeV2.mac
cp -r genrp_4.5GeV2_conf3.root /w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig/genrp_4.5GeV2_conf3_3.root
cd $pwdd
