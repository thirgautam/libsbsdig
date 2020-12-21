#!/bin/sh
pwdd=`pwd`
cd /w/halla-scifs17exp/sbs/thir/v4libsbsdg/buidlN/
cmake -DCMAKE_INSTALL_PREFIX=/w/halla-scifs17exp/sbs/thir/v4libsbsdg/libsbsdig-install/ ../libsbsdig/ && make install
#make install
echo $pwdd
cd $pwdd
sbsdig db/db_gmn_conf.dat gmn13.5_elastic_prod.txt 1000
