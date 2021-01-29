#!/bin/sh

#script to set up the environment for SBS-offline
export LIBSBSDIG=/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install

if test "x$PATH" = "x" ; then
    export PATH=/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install/bin
else
    export PATH=/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install/bin:$PATH
fi

OS=`uname -s`


if [ "$OS" = "Darwin" ]
then # Mac OS: set DYLD_LIBRARY_PATH to library directory:
    if test "x$DYLD_LIBRARY_PATH" = "x"; then
	export DYLD_LIBRARY_PATH=/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install/lib64
    else
	export DYLD_LIBRARY_PATH=/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install/lib64:$DYLD_LIBRARY_PATH
    fi
fi

# set LD_LIBRARY_PATH regardless of OS:
if test "x$LD_LIBRARY_PATH" = "x"; then
    export LD_LIBRARY_PATH=/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install/lib64
else
    export LD_LIBRARY_PATH=/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install/lib64:$LD_LIBRARY_PATH
fi


