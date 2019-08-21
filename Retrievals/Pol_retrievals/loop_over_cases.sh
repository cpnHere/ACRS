#!/bin/bash
for ang in {120,140,160}
do
	for cs in {DY,RC,AC,AP}
	do
  		for d in {1D,3D}
		do
			python pol_cloud_ret.py 0p860 ${ang} ${d} Breon mean _pol_ret_V5 ${cs}
		done
	done
done
