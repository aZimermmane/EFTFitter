#!/bin/bash

# check that hadronizer fragment has subgridpack in the name instead of gridpack !!!

### generic settings
#RESUBMIiT=3
RESUBMIT=1
NJOBS=1000  # irrelevant if input source lhe files
#NJOBS=100 # 4k*25k=10^8 in total for improving the statics 10times
NPARTITION=1000
#NPARTITION=1
#NJOBS=100
#SUBNEVT=10
#SUBNEVT=100000000
SETUP=ttbar_checkgen-lo-rwgt
#SETUP=ttbar_checkgen-lo-SM-nodecay
#DATASET=ttbar_eft_Otwindecay_fixedTWidth_dim6top_ctW_0p67_13TeV
DATASET=ttbareft_translationp_Breuther_dim6top_reweighting_13TeV
#DATASET=eftsamples_ttbar_eft_fcnc_Breuther_dim6top_reweighting_13TeV
#DATASET=ttbareft_SM_nodecay_Breuther_fixedscale_dim6top_ctg_0p0_13TeV
SUBSCRIPTDIR=/nfs/dust/cms/user/zimermma/framework/exec
# MCHI=NONE
# MPHI=NONE
# SUBGRIDPACK=NONE
# SUBSOURCE=GRIDPACK
#INPUTTAG=NONE
#CMSSW_BASE=/nfs/dust/cms/user/zimermma/CMSSW_10_6_19_patch2/
CMSSW_BASE=/nfs/dust/cms/user/zimermma/CMSSW_11_3_4


cd ${CMSSW_BASE}/src/
eval `scram runtime -sh`
cd
### MG5_aMC@NLO EFT


    #GENERATOR=mg5_amcatnlo
OUTPUTTAG=();
    # SUBSCRIPT=(); SUBGRIDPACK=(); FILEBASE=()

    # subtract again __ at the end

    # replace - sign with m

OUTPUTTAG+=( ${DATASET} )

      #  SUBSCRIPT+=( mg5_aMCatNLO_ttbareft_lo_cfg )
      #  FILEBASE+=( TTbarEFT_TuneCUEP8M2T4_13TeV_mg5_amcatnlo_pythia8 )


### start job submission
#zimermma

    #SUBINPUTDIR=/nfs/dust/cms/user/zimermma/samples/cms/${GENERATOR}/lhegen/${INPUTTAG[${IJOBID}]}
    # if [ ${SUBSOURCE} == "LHE" ] ; then
    #     if [ -f ${SUBINPUTDIR}/*.lhe.gz ] ; then
    #         NJOBS=`ls ${SUBINPUTDIR}/*.lhe.gz | wc -l`
    #     else
    #         NJOBS=`ls ${SUBINPUTDIR}/*.root | wc -l`
    #     fi
    # fi
    #TAG=${OUTPUTTAG[${IJOBID}]} ;
TAG=${OUTPUTTAG};
    #replace all . with p
    #TAG=${TAG//./p}

SUBOUTPUTDIR=/nfs/dust/cms/user/zimermma/histograms/${TAG}_with_cpTP
#echo "${SUBOUTPUTDIR}"
    #SUBOUTPUTFILE=${FILEBASE[${IJOBID}]}__NANOGEN

    ### resubmit jobs that failed
    if [ ${RESUBMIT} == "1" ] ; then
         for IJOB in `seq 1 ${NJOBS}`; do
             IJOB=`echo "${IJOB}-1" | bc`  # shift in job number because of shift in htcondor job labeling convention
             ID="$(printf '%05d' "${IJOB}")"
             ID="$(printf "${IJOB}")"
	     if [ ! -f ${SUBOUTPUTDIR}/ttbareft_translationp_Breuther_dim6top_reweighting_13TeV_SM_hist_${IJOB}_with_cpTP.root ] ; 
		then echo "missing file ${SUBOUTPUTDIR}/ttbareft_translationp_Breuther_dim6top_reweighting_13TeV_SM_hist_${IJOB}_with_cpTP.root"
        #         #qsub -N ${TAG} -V -j y -m as -o ${SUBOUTPUTDIR} \
        #         #    -l h_rt=24:00:00 -l h_vmem=6G -l distro=sld6 \
        #         #    -t ${IJOB}:${IJOB}  ${SUBOUTPUTDIR}/scripts/${TAG}.sh
                condor_submit ${SUBOUTPUTDIR}/scripts/${TAG}.re.sub -append arguments=${IJOB}
        #        #condor_submit ${SUBOUTPUTDIR}/scripts/${TAG}.sub -append arguments=${IJOB}
             fi
         done
	      break
    ### initial job submission
    else
        mkdir -p ${SUBOUTPUTDIR}/scripts
        echo ${SUBOUTPUTDIR}
        ### create temporay python script for hadronizer and modify if needed
        # cp ${CMSSW_BASE}/src/${SUBSCRIPTDIR}/${SUBSCRIPT[${IJOBID}]}.py ${SUBOUTPUTDIR}/scripts/${SUBSCRIPT[${IJOBID}]}_${TAG}.py

        #cp ${SUBSCRIPTDIR}/${SETUP}.cc ${SUBOUTPUTDIR}/scripts/${SETUP}_${TAG}.cc
        cp ${SUBSCRIPTDIR}/${SETUP}.cc ${SUBSCRIPTDIR}/${SETUP}_${TAG}.cc
        sed -i -e "s|SUBOUTPUTDIR|${SUBOUTPUTDIR}|g" ${SUBSCRIPTDIR}/${SETUP}_${TAG}.cc
        cd ${SUBSCRIPTDIR}
        make ${SETUP}_${TAG} -f ../Makefile
        cp ${SUBSCRIPTDIR}/${SETUP}_${TAG} ${SUBOUTPUTDIR}/scripts/${SETUP}_${TAG}
    # JOBID=`echo "scale=0; ${NJOBS} - 1" | bc`
    # for IJOBID in `seq 0 ${JOBID}`; do
        ### create job bash file
        cp base_me_ps.sh ${TAG}.sh
        #sed -i -e "s|SUBINPUTDIR|${SUBINPUTDIR}|g" ${TAG}.sh
        sed -i -e "s|CMSSW_VERSION|${CMSSW_BASE}|g" ${TAG}.sh
        sed -i -e "s|SUBOUTPUTDIR|${SUBOUTPUTDIR}|g" ${TAG}.sh
      #  sed -i -e "s|SUBOUTPUTFILE|${SUBOUTPUTFILE}|g" ${TAG}.sh
        sed -i -e "s|SUBSCRIPTDIR|${SUBSCRIPTDIR}|g" ${TAG}.sh
        sed -i -e "s|SUBSCRIPTFILE|${SETUP}_${TAG}|g" ${TAG}.sh
        sed -i -e "s|SUBNEVT|${SUBNEVT}|g" ${TAG}.sh
      #  sed -i -e "s|SUBSOURCE|${SUBSOURCE}|g" ${TAG}.sh
        sed -i -e "s|DATASETCHOICE|${DATASET}|g" ${TAG}.sh
        sed -i -e "s|NPARTITION|${NPARTITION}|g" ${TAG}.sh
        #sed -i -e "s|NJOBS|${NJOBS}|g" ${TAG}.sh
        #sed -i -e "s|IJOBID|${IJOBID}|g" ${TAG}.sh
        chmod +x ${TAG}.sh
        mv ${TAG}.sh ${SUBOUTPUTDIR}/scripts/.

        ### create file with job options
        cp condor.sub ${TAG}.sub
        sed -i -e "s|SUBOUTPUTDIR|${SUBOUTPUTDIR}|g" ${TAG}.sub
        sed -i -e "s|SUBTAG|${TAG}|g" ${TAG}.sub
        sed -i -e "s|SUBNJOBS|${NJOBS}|g" ${TAG}.sub
        sed -i -e "s|SUBMEM|3G|g" ${TAG}.sub
        sed -i -e "s|SUBCORE|1|g" ${TAG}.sub
        sed -i -e "s|SUBRUNTIME|10800|g" ${TAG}.sub
        mv ${TAG}.sub ${SUBOUTPUTDIR}/scripts/.

        # submit job
         echo "condor_submit ${SUBOUTPUTDIR}/scripts/${TAG}.sub"
        condor_submit ${SUBOUTPUTDIR}/scripts/${TAG}.sub
        # echo "done"
        #run 1 job locally
        #${SUBOUTPUTDIR}/scripts/${TAG}.sh 1

#  done
    fi
