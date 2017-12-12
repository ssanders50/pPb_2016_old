#!/bin/sh
ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM185_1/vnanal* > Pbp_HM185_1_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM185_1_list.dat","Pbp_HM185_1_",10)'
rm Pbp_HM185_1_list.dat
ls -1 Pbp_HM185_1_*.root > Pbp_HM185_1_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM185_1_list.dat","Pbp_HM185_1A_",10)'
hadd -f Pbp_HM185_1.root Pbp_HM185_1A_*.root
rm Pbp_HM185_1_*.root
rm Pbp_HM185_1A_*.root
rm Pbp_HM185_1_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM185_2/vnanal* > Pbp_HM185_2_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM185_2_list.dat","Pbp_HM185_2_",10)'
rm Pbp_HM185_2_list.dat
ls -1 Pbp_HM185_2_*.root > Pbp_HM185_2_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM185_2_list.dat","Pbp_HM185_2A_",10)'
hadd -f Pbp_HM185_2.root Pbp_HM185_2A_*.root
rm Pbp_HM185_2_*.root
rm Pbp_HM185_2A_*.root
rm Pbp_HM185_2_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM185_3/vnanal* > Pbp_HM185_3_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM185_3_list.dat","Pbp_HM185_3_",10)'
rm Pbp_HM185_3_list.dat
ls -1 Pbp_HM185_3_*.root > Pbp_HM185_3_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM185_3_list.dat","Pbp_HM185_3A_",10)'
hadd -f Pbp_HM185_3.root Pbp_HM185_3A_*.root
rm Pbp_HM185_3_*.root
rm Pbp_HM185_3A_*.root
rm Pbp_HM185_3_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM185_4/vnanal* > Pbp_HM185_4_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM185_4_list.dat","Pbp_HM185_4_",10)'
rm Pbp_HM185_4_list.dat
ls -1 Pbp_HM185_4_*.root > Pbp_HM185_4_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM185_4_list.dat","Pbp_HM185_4A_",10)'
hadd -f Pbp_HM185_4.root Pbp_HM185_4A_*.root
rm Pbp_HM185_4_*.root
rm Pbp_HM185_4A_*.root
rm Pbp_HM185_4_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM185_5/vnanal* > Pbp_HM185_5_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM185_5_list.dat","Pbp_HM185_5_",10)'
rm Pbp_HM185_5_list.dat
ls -1 Pbp_HM185_5_*.root > Pbp_HM185_5_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM185_5_list.dat","Pbp_HM185_5A_",10)'
hadd -f Pbp_HM185_5.root Pbp_HM185_5A_*.root
rm Pbp_HM185_5_*.root
rm Pbp_HM185_5A_*.root
rm Pbp_HM185_5_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_HM185_6/vnanal* > Pbp_HM185_6_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM185_6_list.dat","Pbp_HM185_6_",10)'
rm Pbp_HM185_6_list.dat
ls -1 Pbp_HM185_6_*.root > Pbp_HM185_6_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_HM185_6_list.dat","Pbp_HM185_6A_",10)'
hadd -f Pbp_HM185_6.root Pbp_HM185_6A_*.root
rm Pbp_HM185_6_*.root
rm Pbp_HM185_6A_*.root
rm Pbp_HM185_6_list.dat


hadd -f Pbp_HM185.root Pbp_HM185_*.root
rm Pbp_HM185_*.root
