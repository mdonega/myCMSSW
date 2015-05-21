cmsDriver.py --filein file:DoublePhoton_pt250-5000_GENSIM.root \
    --fileout file:DoublePhoton_pt250-5000_RECO.root \ \
    --mc \
    --eventcontent AODSIM \
    --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,Configuration/DataProcessing/Utils.addMonitoring \
    --datatier AODSIM \
    --conditions MCRUN2_74_V7 \
    --step RAW2DIGI,L1Reco,RECO \
    --magField 38T_PostLS1 \
    --python_filename step_RECO.py \
    --no_exec -n 2 || exit $? ;

