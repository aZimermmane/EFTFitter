#!/bin/bash

# ${1} job ID

echo -e "\nJob started at "`date`" on "`hostname --fqdn`"\n"

# set up cms shortcuts
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_99 x86_64-centos7-gcc10-opt
# #HTCONDOR
if [ -z "$1" ] ; then
    echo "Need to give an operator name!"
    exit 0
fi
/nfs/dust/cms/user/zimermma/EFTFitter/exec/wbern_evx_total "$1"

#done
echo -e "\nJob stoped at "`date`" on "`hostname --fqdn`"\n"

exit 0
