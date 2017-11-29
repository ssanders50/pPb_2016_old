#!/bin/sh
ls -1 /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM120_1/vnanal* > pPb_HM120_1_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM120_1_list.dat","pPb_HM120_1_",10)'
rm pPb_HM120_1_list.dat
ls -1 pPb_HM120_1_*.root > pPb_HM120_1_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM120_1_list.dat","pPb_HM120_1A_",10)'
hadd -f pPb_HM120_1.root pPb_HM120_1A_*.root
rm pPb_HM120_1_*.root
rm pPb_HM120_1A_*.root
rm pPb_HM120_1_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM120_2/vnanal* > pPb_HM120_2_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM120_2_list.dat","pPb_HM120_2_",10)'
rm pPb_HM120_2_list.dat
ls -1 pPb_HM120_2_*.root > pPb_HM120_2_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM120_2_list.dat","pPb_HM120_2A_",10)'
hadd -f pPb_HM120_2.root pPb_HM120_2A_*.root
rm pPb_HM120_2_*.root
rm pPb_HM120_12_*.root
rm pPb_HM120_2_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM120_3/vnanal* > pPb_HM120_3_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM120_3_list.dat","pPb_HM120_3_",10)'
rm pPb_HM120_3_list.dat
ls -1 pPb_HM120_3_*.root > pPb_HM120_3_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM120_3_list.dat","pPb_HM120_3A_",10)'
hadd -f pPb_HM120_3.root pPb_HM120_3A_*.root
rm pPb_HM120_3_*.root
rm pPb_HM120_3A_*.root
rm pPb_HM120_3_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM120_4/vnanal* > pPb_HM120_4_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM120_4_list.dat","pPb_HM120_4_",10)'
rm pPb_HM120_4_list.dat
ls -1 pPb_HM120_4_*.root > pPb_HM120_4_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM120_4_list.dat","pPb_HM120_4A_",10)'
hadd -f pPb_HM120_4.root pPb_HM120_4A_*.root
rm pPb_HM120_4_*.root
rm pPb_HM120_4A_*.root
rm pPb_HM120_4_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM120_5/vnanal* > pPb_HM120_5_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM120_5_list.dat","pPb_HM120_5_",10)'
rm pPb_HM120_5_list.dat
ls -1 pPb_HM120_5_*.root > pPb_HM120_5_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM120_5_list.dat","pPb_HM120_5A_",10)'
hadd -f pPb_HM120_5.root pPb_HM120_5A_*.root
rm pPb_HM120_5_*.root
rm pPb_HM120_5A_*.root
rm pPb_HM120_5_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM120_6/vnanal* > pPb_HM120_6_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM120_6_list.dat","pPb_HM120_6_",10)'
rm pPb_HM120_6_list.dat
ls -1 pPb_HM120_6_*.root > pPb_HM120_6_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM120_6_list.dat","pPb_HM120_6A_",10)'
hadd -f pPb_HM120_6.root pPb_HM120_6A_*.root
rm pPb_HM120_6_*.root
rm pPb_HM120_6A_*.root
rm pPb_HM120_6_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM120_7/vnanal* > pPb_HM120_7_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM120_7_list.dat","pPb_HM120_7_",10)'
rm pPb_HM120_7_list.dat
ls -1 pPb_HM120_7_*.root > pPb_HM120_7_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM120_7_list.dat","pPb_HM120_7A_",10)'
hadd -f pPb_HM120_7.root pPb_HM120_7A_*.root
rm pPb_HM120_7_*.root
rm pPb_HM120_7A_*.root
rm pPb_HM120_7_list.dat
