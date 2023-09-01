#!/bin/bash
SEQUENCES="${1}"
METADATA="${2}"
OUT="${3}"
AMPLICONWIDTH="${4}"
MAXMISMATCH="${5}"
MISALIGNTHRESHOLD="${6}"
PRIMERWIDTH="${7}"
SEARCHWIDTH="${8}"
AMPLIFIABILITY="${9}"
NUMAMPLICONS="${10}"
NUMCORES="${11}"
NUMSEQUENCES="${12}"
SEED="${13}"
MINNONALIGN="${14}"
BETA="${15}"

## Note that this assumes you are in the base folder, if running from manuscript folder use ../AmpliDiff/AmpliDiff.py
python AmpliDiff/AmpliDiff.py $SEQFOLDER $METAFOLDER -o $OUT \
        -aw $AMPLICONWIDTH \
        -mm $MAXMISMATCH \
        -mt $MISALIGNTHRESHOLD \
        -pw $PRIMERWIDTH \
        -sw $SEARCHWIDTH \
        -cov $AMPLIFIABILITY \
        -amps $NUMAMPLICONS \
        -c $NUMCORES \
        -n $NUMSEQUENCES \
        -sd $SEED \
        --min_non_align $MINNONALIGN \
        --beta $SEED \
        