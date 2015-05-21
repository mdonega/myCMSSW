cmsDriver.py  Configuration/GenProduction/python/EGM-RunIIWinter15GS-00004-fragment.py \
    --fileout file:DoublePhoton_pt250-5000_GENSIM.root \
    --mc \
    --pileup Flat_20_50 \
    --pileup_input "dbs:/MinBias_TuneCUETP8M1_13TeV-pythia8/RunIIFall14GS-MCRUN2_71_V1-v3/GEN-SIM" \
    --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,Configuration/DataProcessing/Utils.addMonitoring \
    --eventcontent RAWSIM \
    --datatier GEN-SIM-RAW \
    --conditions MCRUN2_74_V7 \
    --step GEN,SIM,DIGI,L1,DIGI2RAW \
    --magField 38T_PostLS1  \
    --beamspot NominalCollision2015 \
    --python_filename step_GENSIM.py \
    --no_exec -n 3 || exit $? ;


