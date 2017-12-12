#!/bin/sh
ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM250_1/vnanal* > Pbp_HM250_1_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM250_1_list.dat","Pbp_HM250_1_",10)'
rm Pbp_HM250_1_list.dat
ls -1 Pbp_HM250_1_*.root > Pbp_HM250_1_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM250_1_list.dat","Pbp_HM250_1A_",10)'
hadd -f Pbp_HM250_1.root Pbp_HM250_1A_*.root
rm Pbp_HM250_1_*.root
rm Pbp_HM250_1A_*.root
rm Pbp_HM250_1_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM250_2/vnanal* > Pbp_HM250_2_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM250_2_list.dat","Pbp_HM250_2_",10)'
rm Pbp_HM250_2_list.dat
ls -1 Pbp_HM250_2_*.root > Pbp_HM250_2_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM250_2_list.dat","Pbp_HM250_2A_",10)'
hadd -f Pbp_HM250_2.root Pbp_HM250_2A_*.root
rm Pbp_HM250_2_*.root
rm Pbp_HM250_2A_*.root
rm Pbp_HM250_2_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM250_3/vnanal* > Pbp_HM250_3_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM250_3_list.dat","Pbp_HM250_3_",10)'
rm Pbp_HM250_3_list.dat
ls -1 Pbp_HM250_3_*.root > Pbp_HM250_3_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM250_3_list.dat","Pbp_HM250_3A_",10)'
hadd -f Pbp_HM250_3.root Pbp_HM250_3A_*.root
rm Pbp_HM250_3_*.root
rm Pbp_HM250_3A_*.root
rm Pbp_HM250_3_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM250_4/vnanal* > Pbp_HM250_4_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM250_4_list.dat","Pbp_HM250_4_",10)'
rm Pbp_HM250_4_list.dat
ls -1 Pbp_HM250_4_*.root > Pbp_HM250_4_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM250_4_list.dat","Pbp_HM250_4A_",10)'
hadd -f Pbp_HM250_4.root Pbp_HM250_4A_*.root
rm Pbp_HM250_4_*.root
rm Pbp_HM250_4A_*.root
rm Pbp_HM250_4_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM250_5/vnanal* > Pbp_HM250_5_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM250_5_list.dat","Pbp_HM250_5_",10)'
rm Pbp_HM250_5_list.dat
ls -1 Pbp_HM250_5_*.root > Pbp_HM250_5_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM250_5_list.dat","Pbp_HM250_5A_",10)'
hadd -f Pbp_HM250_5.root Pbp_HM250_5A_*.root
rm Pbp_HM250_5_*.root
rm Pbp_HM250_5A_*.root
rm Pbp_HM250_5_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM250_6/vnanal* > Pbp_HM250_6_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM250_6_list.dat","Pbp_HM250_6_",10)'
rm Pbp_HM250_6_list.dat
ls -1 Pbp_HM250_6_*.root > Pbp_HM250_6_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM250_6_list.dat","Pbp_HM250_6A_",10)'
hadd -f Pbp_HM250_6.root Pbp_HM250_6A_*.root
rm Pbp_HM250_6_*.root
rm Pbp_HM250_6A_*.root
rm Pbp_HM250_6_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM250_7/vnanal* > Pbp_HM250_7_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM250_7_list.dat","Pbp_HM250_7_",10)'
rm Pbp_HM250_7_list.dat
ls -1 Pbp_HM250_7_*.root > Pbp_HM250_7_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM250_7_list.dat","Pbp_HM250_7A_",10)'
hadd -f Pbp_HM250_7.root Pbp_HM250_7A_*.root
rm Pbp_HM250_7_*.root
rm Pbp_HM250_7A_*.root
rm Pbp_HM250_7_list.dat

hadd -f Pbp_HM250.root Pbp_HM250_*.root
rm Pbp_HM250_*.root
