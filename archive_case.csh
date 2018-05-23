#!/bin/csh
#
# setup some environment variables
#

setenv RUNDIR /project/convection/juliob/pce-data/${1}
setenv RUNSRC /project/convection/juliob/pce-data/${1}/source/

mkdir ${RUNDIR}
mkdir ${RUNSRC}


cp *.{f,F90,code}  ${RUNSRC}
cp -r share    ${RUNSRC}
cp -r idlpros  ${RUNDIR}
cp control.nml ${RUNDIR}
cp pces.x      ${RUNDIR}
cp arm_scmx.dat    ${RUNDIR}
cp fort.*      ${RUNDIR}


exit

