#!/bin/bash

# The script creates 2 .mcool files with merged data for neurons
# and non-neurons separately. The final Hi-C matrices are sampled to the
# same size.
#
# Run this script from the directory with .mcool files to be merged.
# Modify "PLUS_FILES" and "MINUS_FILES" with correct names of mcool files
# for neurons and non-neurons accordingly.
# 
# Input mcool files can be obtained from GSE229816.
#
# Final output files will have prefixes "NeuNplus.sampled*" and
# "NeuNminus.sampled*"


PLUS_FILES=("GSM7178441_HC-91plus.hg38.mapq_30.1000.mcool" 
            "GSM7178442_HC-318_NeuN_plus.hg38.mapq_30.1000.mcool"
            "GSM7178443_HC2M_NeuN_plus.hg38.mapq_30.1000.mcool"
            "GSM7178444_HC3M_NeuN_plus.hg38.mapq_30.1000.mcool")

MINUS_FILES=("GSM7178445_HC-91minus.hg38.mapq_30.1000.mcool"
             "GSM7178446_HC-318_NeuN_minus.hg38.mapq_30.1000.mcool"
             "GSM7178447_HC2M_NeuN_minus.hg38.mapq_30.1000.mcool"
             "GSM7178448_HC3M_NeuN_minus.hg38.mapq_30.1000.mcool")


# Merge coolers for NeuN+ and NeuN- separately
MRG_RES=1000
PLUS_MRG_OUT="NeuNplus.merged.${MRG_RES}.cool"
MINUS_MRG_OUT="NeuNminus.merged.${MRG_RES}.cool"
PLUS_CLR=("${PLUS_FILES[@]/%/::resolutions/${MRG_RES}}")
MINUS_CLR=("${MINUS_FILES[@]/%/::resolutions/${MRG_RES}}")

echo "Merge coolers"
cooler merge $PLUS_MRG_OUT "${PLUS_CLR[@]}"
cooler merge $MINUS_MRG_OUT "${MINUS_CLR[@]}"

# Remove the main diagonal
echo "Remove the main diagonal.."
echo "Neurons"
PLUS_DROP_D_OUT="NeuNplus.drop_diag.${MRG_RES}.cool"
cooler dump -t bins $PLUS_MRG_OUT > tmp.plus.hic.1kb_bins.txt
cooler dump -t pixels $PLUS_MRG_OUT | awk -F'\t' '$1!=$2' > tmp.plus.hic.1kb_pixels.txt
cooler load -f coo --assembly hg38 tmp.plus.hic.1kb_bins.txt tmp.plus.hic.1kb_pixels.txt $PLUS_DROP_D_OUT
rm tmp.plus.hic.1kb_bins.txt tmp.plus.hic.1kb_pixels.txt

echo "Non-neurons"
MINUS_DROP_D_OUT="NeuNminus.drop_diag.${MRG_RES}.cool"
cooler dump -t bins $MINUS_MRG_OUT > tmp.minus.hic.1kb_bins.txt
cooler dump -t pixels $MINUS_MRG_OUT | awk -F'\t' '$1!=$2' > tmp.minus.hic.1kb_pixels.txt
cooler load -f coo --assembly hg38 tmp.minus.hic.1kb_bins.txt tmp.minus.hic.1kb_pixels.txt $MINUS_DROP_D_OUT
rm tmp.minus.hic.1kb_bins.txt tmp.minus.hic.1kb_pixels.txt

# Sample merged coolers to the same size
PLUS_SUM=$(cooler info ${PLUS_DROP_D_OUT} | sed -n 's/    "sum": //p')
MINUS_SUM=$(cooler info ${MINUS_DROP_D_OUT} | sed -n 's/    "sum": //p')
PLUS_SAMPLE_OUT="NeuNplus.sampled.${MRG_RES}.cool"
MINUS_SAMPLE_OUT="NeuNminus.sampled.${MRG_RES}.cool"

if [ "$PLUS_SUM" -gt "$MINUS_SUM" ]
then
    echo "Hi-C for neurons has more counts, sampling neurons.."
    cooltools random-sample -c $MINUS_SUM $PLUS_DROP_D_OUT $PLUS_SAMPLE_OUT
    cp $MINUS_DROP_D_OUT $MINUS_SAMPLE_OUT
else
    echo "Hi-C for non-neurons has more counts, sampling non-neurons.."
    cooltools random-sample -c $PLUS_SUM $MINUS_DROP_D_OUT $MINUS_SAMPLE_OUT
    cp $PLUS_DROP_D_OUT $PLUS_SAMPLE_OUT
fi

# zoomify cooler to create maps with different resolutions
cooler zoomify --balance -r 1000,2000,5000,10000,15000,20000,25000,40000,50000,100000,150000,200000,250000,500000,1000000,2000000,5000000 $PLUS_SAMPLE_OUT
cooler zoomify --balance -r 1000,2000,5000,10000,15000,20000,25000,40000,50000,100000,150000,200000,250000,500000,1000000,2000000,5000000 $MINUS_SAMPLE_OUT
