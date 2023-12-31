#/bin/zsh
TASSEL5="/Applications/TASSEL 5/run_pipeline.pl"
INPUT_VCF="../data/AIR_PTRILS_merged.vcf"
NOPTRIL_VCF="AIR_B73.vcf"

# This list has LANAP duplicated taxa removed

AIRB73_LIST="../data/AIR_ICPMS_taxa_B73.list" 
IMPUTED_VCF="AIR_B73_imputed"

# Remove PTRILs taxa from VCF
"$TASSEL5" -vcf ${INPUT_VCF} \
           -includeTaxaInFile $AIRB73_LIST \
           -export ${NOPTRIL_VCF} \
           -exportType VCF
           
           
# Impute missing data, remove LANAP duplicated taxa

"$TASSEL5" -vcf $NOHET_VCF \
         -FilterSiteBuilderPlugin \
         -siteMinCount 1200 \
         -endPlugin \
         -FilterTaxaBuilderPlugin \
         -minNotMissing 0.8 \
         -endPlugin \
         -LDKNNiImputationPlugin \
         -endPlugin \
         -includeTaxaInFile $AIRB73_LIST \
         -export $IMPUTED_VCF   -exportType VCF
