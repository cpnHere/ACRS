#!/bin/bash
# To run mean polarimetric retrievals
# CPN (03/29/2020)
#
for sza in 120 140 160
do
	for res in nt 0p1km 0p5km 1p0km 5p0km
	do
		python pol_cloud_ret.py 0p860 $sza 3D Breon mean _pol_ret_V7 DY $res
	done
done
