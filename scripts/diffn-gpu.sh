#!/bin/bash

function diffn {
  threshold=1e-11

  awk -v thr="$threshold" '{
      all_rows=all_rows $0 "\n"; # Store all rows
      if (NF < 10) next; # Skip non-data lines

      min = $8+0; mean = $9+0; max = $10+0; # Convert to numerical values

      if ((min < -thr || min > thr) || (mean < -thr || mean > thr) || (max < -thr || max > thr)) {
          count_differs++;
      }
      total++;
  }

  END {
      if (count_differs > 0) {
          print all_rows;
          printf "%d out of %d records differ\n", count_differs, total;
          exit 1;
      }
  }'
}

cdo infon -sub $1 $2 | diffn
