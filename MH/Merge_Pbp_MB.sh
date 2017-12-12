#!/bin/sh
ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB1/vnanal* > Pbp_MB_1_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB_1_list.dat","Pbp_MB_1_",10)'
rm Pbp_MB_1_list.dat
ls -1 Pbp_MB_1_*.root > Pbp_MB_1_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB_1_list.dat","Pbp_MB_1A_",10)'
hadd -f Pbp_MB_1.root Pbp_MB_1A_*.root
rm Pbp_MB_1_*.root
rm Pbp_MB_1A_*.root
rm Pbp_MB_1_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB2/vnanal* > Pbp_MB_2_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB_2_list.dat","Pbp_MB_2_",10)'
rm Pbp_MB_2_list.dat
ls -1 Pbp_MB_2_*.root > Pbp_MB_2_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB_2_list.dat","Pbp_MB_2A_",10)'
hadd -f Pbp_MB_2.root Pbp_MB_2A_*.root
rm Pbp_MB_2_*.root
rm Pbp_MB_2A_*.root
rm Pbp_MB_2_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB3/vnanal* > Pbp_MB_3_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB_3_list.dat","Pbp_MB_3_",10)'
rm Pbp_MB_3_list.dat
ls -1 Pbp_MB_3_*.root > Pbp_MB_3_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB_3_list.dat","Pbp_MB_3A_",10)'
hadd -f Pbp_MB_3.root Pbp_MB_3A_*.root
rm Pbp_MB_3_*.root
rm Pbp_MB_3A_*.root
rm Pbp_MB_3_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB4/vnanal* > Pbp_MB_4_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB_4_list.dat","Pbp_MB_4_",10)'
rm Pbp_MB_4_list.dat
ls -1 Pbp_MB_4_*.root > Pbp_MB_4_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB_4_list.dat","Pbp_MB_4A_",10)'
hadd -f Pbp_MB_4.root Pbp_MB_4A_*.root
rm Pbp_MB_4_*.root
rm Pbp_MB_4A_*.root
rm Pbp_MB_4_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB5/vnanal* > Pbp_MB_5_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB_5_list.dat","Pbp_MB_5_",10)'
rm Pbp_MB_5_list.dat
ls -1 Pbp_MB_5_*.root > Pbp_MB_5_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB_5_list.dat","Pbp_MB_5A_",10)'
hadd -f Pbp_MB_5.root Pbp_MB_5A_*.root
rm Pbp_MB_5_*.root
rm Pbp_MB_5A_*.root
rm Pbp_MB_5_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB6/vnanal* > Pbp_MB_6_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB_6_list.dat","Pbp_MB_6_",10)'
rm Pbp_MB_6_list.dat
ls -1 Pbp_MB_6_*.root > Pbp_MB_6_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB_6_list.dat","Pbp_MB_6A_",10)'
hadd -f Pbp_MB_6.root Pbp_MB_6A_*.root
rm Pbp_MB_6_*.root
rm Pbp_MB_6A_*.root
rm Pbp_MB_6_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB7/vnanal* > Pbp_MB_7_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB_7_list.dat","Pbp_MB_7_",10)'
rm Pbp_MB_7_list.dat
ls -1 Pbp_MB_7_*.root > Pbp_MB_7_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB_7_list.dat","Pbp_MB_7A_",10)'
hadd -f Pbp_MB_7.root Pbp_MB_7A_*.root
rm Pbp_MB_7_*.root
rm Pbp_MB_7A_*.root
rm Pbp_MB_7_list.dat

hadd -f Pbp_MB.root Pbp_MB_*.root
rm Pbp_MB_*.root
