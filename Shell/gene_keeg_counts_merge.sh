#!/bin/bash

# Default empty
COUNTS=""
ANNO=""
OUT=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--counts)
            COUNTS="$2"
            shift 2
            ;;
        -a|--anno|--annotation)
            ANNO="$2"
            shift 2
            ;;
        -o|--out|--output)
            OUT="$2"
            shift 2
            ;;
        -*)
            echo "Unknown option: $1"
            exit 1
            ;;
        *)
            shift
            ;;
    esac
done

# Check required args
if [[ -z "$COUNTS" || -z "$ANNO" || -z "$OUT" ]]; then
    echo "Usage: $0 -c counts.txt -a annotation.txt -o merged.txt"
    exit 1
fi

# Core AWK logic
awk '
    NR==FNR {
        if (FNR>1) {
            counts[$1]=$2
            all[$1]=1
        }
        next
    }

    FNR==1 {
        print "gene\tmapped\tKOID"
        next
    }

    {
        gene=$1
        ko=$2
        if (ko=="") ko="NA"
        koid[gene]=ko
        all[gene]=1
    }

    END {
        for (g in all) {
            m = (g in counts ? counts[g] : "NA")
            if (!(g in koid)) koid[g]="NA"
            print g "\t" m "\t" koid[g]
        }
    }
' "$COUNTS" "$ANNO" | sort -k1,1 > "$OUT"

echo "Merged file written to: $OUT"
