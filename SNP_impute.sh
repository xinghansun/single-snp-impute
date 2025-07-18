#!/bin/bash
set -euo pipefail

IND_QC='true'

# ─────────────────────────────────────────────
# Print help message
# ─────────────────────────────────────────────
show_help() {
  cat <<EOF
SNP_impute v0.1 (16 June 2021) by Xinghan Sun

Usage: $(basename "$0") -c CONFIG -g GENOFILE -o OUTPREFIX [-a]

Options:
  -c    Path to config file
  -g    Genotype file prefix (no extension)
  -o    Output prefix
  -a    Use all samples (skip individual QC)
  -h    Show this help message

Example:
  bash $(basename "$0") -c a.config -g chip_data -o run123 -a
EOF
}

# ─────────────────────────────────────────────
# Parse input arguments
# ─────────────────────────────────────────────
while getopts "c:g:o:ah" opt; do
  case $opt in
    c)
      CONFIG_NAME=$OPTARG
      echo "* Reading SNP imputation config: ${CONFIG_NAME}"
      source "${CONFIG_NAME}"
      echo -e "\nConfig content:\n"
      cat "${CONFIG_NAME}"
      echo ""
      ;;
    g)
      GENOFILE_PREFIX=$OPTARG
      echo "* Genotype data: ${GENOFILE_PREFIX}.bed/.bim/.fam"
      ;;
    o)
      OUT_PREFIX=$OPTARG
      echo "* Output prefix: ${OUT_PREFIX}"
      ;;
    a)
      IND_QC='false'
      echo "* Using all samples. Skipping individual QC."
      ;;
    h)
      show_help
      exit 0
      ;;
    *)
      show_help
      exit 1
      ;;
  esac
done

# ─────────────────────────────────────────────
# Check required variables
# ─────────────────────────────────────────────
: "${CONFIG_NAME:?Missing -c config file}"
: "${GENOFILE_PREFIX:?Missing -g genotype file prefix}"
: "${OUT_PREFIX:?Missing -o output prefix}"
: "${SNP:?SNP not defined in config}"
: "${CHROMOSOME:?CHROMOSOME not defined in config}"
: "${POSITION:?POSITION not defined in config}"
: "${WINDOW_SIZE:?WINDOW_SIZE not defined in config}"
: "${REFHAPS:?REFHAPS not defined in config}"

# ─────────────────────────────────────────────
# Define paths and directories
# ─────────────────────────────────────────────
BASE_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUT_DIR="${BASE_PATH}/Output/${OUT_PREFIX}"
mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

PLINK="${BASE_PATH}/Tools/plink"
SHAPEIT="${BASE_PATH}/Tools/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit"
MINIMAC="${BASE_PATH}/Tools/Minimac3/bin/Minimac3"
MAP_FILE="${BASE_PATH}/Ref/Chr_${CHROMOSOME}/genetic_map_GRCh37_chr${CHROMOSOME}_for_shapeit.txt"
REF_HAPS="${BASE_PATH}/Ref/Chr_${CHROMOSOME}/${REFHAPS}"

# ─────────────────────────────────────────────
# Functions
# ─────────────────────────────────────────────

perform_snp_qc() {
  echo "* Performing SNP-level QC..."
  "${PLINK}" --bfile "${GENOFILE_PREFIX}" --geno 0.02 --make-bed --out data
}

perform_ind_qc() {
  if [[ "${IND_QC}" == "true" ]]; then
    echo "* Performing individual-level QC..."
    "${PLINK}" --bfile data --mind 0.02 --make-bed --out data
  fi
}

extract_window() {
  echo "* Extracting window around SNP ${SNP}..."
  if grep -q "${SNP}" data.bim; then
    "${PLINK}" --bfile data --chr "${CHROMOSOME}" \
      --from-bp $((POSITION - WINDOW_SIZE)) \
      --to-bp $((POSITION + WINDOW_SIZE)) \
      --exclude-snp "${SNP}" --make-bed --out data
  else
    "${PLINK}" --bfile data --chr "${CHROMOSOME}" \
      --from-bp $((POSITION - WINDOW_SIZE)) \
      --to-bp $((POSITION + WINDOW_SIZE)) \
      --make-bed --out data
  fi
}

apply_maf_filter() {
  echo "* Applying MAF filter (MAF > 0.001)..."
  "${PLINK}" --bfile data --maf 0.001 --make-bed --out data
}

apply_hwe_filter() {
  echo "* Applying HWE filter (P > 1e-6)..."
  "${PLINK}" --bfile data --hwe 1e-6 --make-bed --out "${OUT_PREFIX}"
  rm -f data.*
}

perform_phasing() {
  echo "* Running SHAPEIT for phasing..."
  "${SHAPEIT}" -B "${OUT_PREFIX}" \
    -M "${MAP_FILE}" \
    -O "${OUT_PREFIX}.phased" \
    --input-from $((POSITION - WINDOW_SIZE)) \
    --input-to $((POSITION + WINDOW_SIZE)) \
    --seed 1234 --thread 10
}

convert_to_vcf() {
  echo "* Converting haplotype to VCF..."
  "${SHAPEIT}" -convert \
    --input-haps "${OUT_PREFIX}.phased" \
    --output-vcf "${OUT_PREFIX}.phased.vcf"
}

perform_imputation() {
  echo "* Running Minimac3 for imputation..."
  "${MINIMAC}" --refHaps "${REF_HAPS}" \
    --haps "${OUT_PREFIX}.phased.vcf" \
    --prefix "${OUT_PREFIX}.imputed.output" \
    --start $((POSITION - 1)) \
    --end $((POSITION + 1)) \
    --chr "${CHROMOSOME}" \
    --window "${WINDOW_SIZE}" \
    --cpus 20
}

# ─────────────────────────────────────────────
# Main workflow
# ─────────────────────────────────────────────
perform_snp_qc
perform_ind_qc
extract_window
apply_maf_filter
apply_hwe_filter
perform_phasing
convert_to_vcf
perform_imputation

echo "SNP imputation finished, results saved in: ${OUT_PREFIX}.imputed.ouput"


