#/bin/bash

years="1991"
months="01 02 03 04 05 06 07 08 09 10 11 12"
days="01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31"

for year in $years; do

    for month in $months; do

	for day in $days; do

	    bsub -q short-serial -W 06:00 -o out.log -e out.err python2.7 collocate_hirs_avhrr.py noaa10 v2.7.1-60-ge0a0928 AVHRR10_G $year $month $day

	done
	
    done

done
