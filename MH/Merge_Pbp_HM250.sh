#!/bin/sh
ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM250_7/vnanal* > Pbp_HM250_7_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM250_7_list.dat","Pbp_HM250_7_",10)'
rm Pbp_HM250_7_list.dat
ls -1 Pbp_HM250_7_*.root > Pbp_HM250_7_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM250_7_list.dat","Pbp_HM250_7A_",10)'
hadd -f Pbp_HM250.root Pbp_HM250_7A_*.root
rm Pbp_HM250_7_*.root
rm Pbp_HM250_7A_*.root
rm Pbp_HM250_7_list.dat
