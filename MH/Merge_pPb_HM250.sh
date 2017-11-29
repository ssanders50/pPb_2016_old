#!/bin/sh
ls -1 /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM185_1/vnanal* > pPb_HM185_1_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM185_1_list.dat","pPb_HM185_1_",10)'
rm pPb_HM185_1_list.dat
ls -1 pPb_HM185_1_*.root > pPb_HM185_1_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM185_1_list.dat","pPb_HM185_1A_",10)'
hadd -f pPb_HM185_1.root pPb_HM185_1A_*.root
rm pPb_HM185_1_*.root
rm pPb_HM185_1A_*.root
rm pPb_HM185_1_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM185_2/vnanal* > pPb_HM185_2_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM185_2_list.dat","pPb_HM185_2_",10)'
rm pPb_HM185_2_list.dat
ls -1 pPb_HM185_2_*.root > pPb_HM185_2_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM185_2_list.dat","pPb_HM185_2A_",10)'
hadd -f pPb_HM185_2.root pPb_HM185_2A_*.root
rm pPb_HM185_2_*.root
rm pPb_HM185_12_*.root
rm pPb_HM185_2_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM185_3/vnanal* > pPb_HM185_3_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM185_3_list.dat","pPb_HM185_3_",10)'
rm pPb_HM185_3_list.dat
ls -1 pPb_HM185_3_*.root > pPb_HM185_3_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM185_3_list.dat","pPb_HM185_3A_",10)'
hadd -f pPb_HM185_3.root pPb_HM185_3A_*.root
rm pPb_HM185_3_*.root
rm pPb_HM185_3A_*.root
rm pPb_HM185_3_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM185_4/vnanal* > pPb_HM185_4_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM185_4_list.dat","pPb_HM185_4_",10)'
rm pPb_HM185_4_list.dat
ls -1 pPb_HM185_4_*.root > pPb_HM185_4_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM185_4_list.dat","pPb_HM185_4A_",10)'
hadd -f pPb_HM185_4.root pPb_HM185_4A_*.root
rm pPb_HM185_4_*.root
rm pPb_HM185_4A_*.root
rm pPb_HM185_4_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM185_5/vnanal* > pPb_HM185_5_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM185_5_list.dat","pPb_HM185_5_",10)'
rm pPb_HM185_5_list.dat
ls -1 pPb_HM185_5_*.root > pPb_HM185_5_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM185_5_list.dat","pPb_HM185_5A_",10)'
hadd -f pPb_HM185_5.root pPb_HM185_5A_*.root
rm pPb_HM185_5_*.root
rm pPb_HM185_5A_*.root
rm pPb_HM185_5_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM185_6/vnanal* > pPb_HM185_6_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM185_6_list.dat","pPb_HM185_6_",10)'
rm pPb_HM185_6_list.dat
ls -1 pPb_HM185_6_*.root > pPb_HM185_6_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM185_6_list.dat","pPb_HM185_6A_",10)'
hadd -f pPb_HM185_6.root pPb_HM185_6A_*.root
rm pPb_HM185_6_*.root
rm pPb_HM185_6A_*.root
rm pPb_HM185_6_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM185_7/vnanal* > pPb_HM185_7_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM185_7_list.dat","pPb_HM185_7_",10)'
rm pPb_HM185_7_list.dat
ls -1 pPb_HM185_7_*.root > pPb_HM185_7_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM185_7_list.dat","pPb_HM185_7A_",10)'
hadd -f pPb_HM185_7.root pPb_HM185_7A_*.root
rm pPb_HM185_7_*.root
rm pPb_HM185_7A_*.root
rm pPb_HM185_7_list.dat
