#!/bin/bash

# ${1} job ID

echo -e "\nJob started at "`date`" on "`hostname --fqdn`"\n"

# set up cms shortcuts
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_99 x86_64-centos7-gcc10-opt
# #HTCONDOR
if [[ -z "$1" || -z "$2" || -z "$3" ]]; then
    echo "Not enough positional arguments (at least 3 required)"
    exit 0
fi
/nfs/dust/cms/user/zimermma/EFTFitter/exec/wbern_1obs "$1" "$2" "$3" "$4"

#done
echo -e "\nJob stoped at "`date`" on "`hostname --fqdn`"\n"

exit 0
