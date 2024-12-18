#!/usr/bin/env bash
#
# Script name: cgmap_to_cx.sh
#
# Usage: ./cgmap_to_cx.sh <input.cgmap> <output.CX_report>
# Description: Convert a CGmap-formatted methylation file to a Bismark CX-like format.

# --- Check arguments ---
if [ $# -ne 2 ]; then
  echo "Usage: $0 <input.cgmap> <output.CX_report>"
  exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"

# --- Main conversion ---
awk 'BEGIN {
  FS="\t"; OFS="\t"
}
{
  # Skip header lines if your CGmap file has them; adapt as needed:
  # if ($1 == "chrom" && $2 == "nucleotide") { next }

  # Determine strand: if reference is C, it's the + strand; if G, it's - strand.
  if ($2 == "C") {
    strand = "+"
  } else {
    strand = "-"
  }

  # Calculate unmethylated reads
  unmethyl = $8 - $7

  # Print in a Bismark CX-like format:
  # 1) Chromosome    -> $1
  # 2) Position      -> $3
  # 3) Strand        -> strand
  # 4) Context       -> $4  (CG, CHG, CHH, etc.)
  # 5) Methylated    -> $7
  # 6) Unmethylated  -> (total coverage - methylated)
  
  print $1, $3, strand, $4, $7, unmethyl
}' "${INPUT_FILE}" > "${OUTPUT_FILE}"

echo "Conversion complete."
echo "Output written to: ${OUTPUT_FILE}"
