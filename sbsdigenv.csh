#!/bin/csh

#script to set up the environment for SBS-offline
setenv LIBSBSDIG /w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install

if( ! ${?PATH} ) then
    setenv PATH /w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install/bin
else
    setenv PATH /w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install/bin:${PATH}
endif

set OS=`uname -s`


if( "$OS" == "Darwin" ) then # Mac OS: set DYLD_LIBRARY_PATH to library directory:
    if(! ${?DYLD_LIBRARY_PATH}) then
	setenv DYLD_LIBRARY_PATH /w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install/lib64
    else
	setenv DYLD_LIBRARY_PATH /w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install/lib64:${DYLD_LIBRARY_PATH}
    endif
endif

# set LD_LIBRARY_PATH regardless of OS:
if( ! ${?LD_LIBRARY_PATH}) then
    setenv LD_LIBRARY_PATH /w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install/lib64
else
    setenv LD_LIBRARY_PATH /w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install/lib64:${LD_LIBRARY_PATH}
endif


