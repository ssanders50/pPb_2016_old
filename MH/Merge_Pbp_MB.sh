#!/bin/sh
ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB1/vnanal* > Pbp_MB1_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB1_list.dat","Pbp_MB1_",10)'
rm Pbp_MB1_list.dat
ls -1 Pbp_MB1_*.root > Pbp_MB1_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB1_list.dat","Pbp_MB1A_",10)'
hadd -f Pbp_MB1.root Pbp_MB1A_*.root
rm Pbp_MB1_*.root
rm Pbp_MB1A_*.root
rm Pbp_MB1_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB2/vnanal* > Pbp_MB2_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB2_list.dat","Pbp_MB2_",10)'
rm Pbp_MB2_list.dat
ls -1 Pbp_MB2_*.root > Pbp_MB2_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB2_list.dat","Pbp_MB2A_",10)'
hadd -f Pbp_MB2.root Pbp_MB2A_*.root
rm Pbp_MB2_*.root
rm Pbp_MB12_*.root
rm Pbp_MB2_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB3/vnanal* > Pbp_MB3_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB3_list.dat","Pbp_MB3_",10)'
rm Pbp_MB3_list.dat
ls -1 Pbp_MB3_*.root > Pbp_MB3_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB3_list.dat","Pbp_MB3A_",10)'
hadd -f Pbp_MB3.root Pbp_MB3A_*.root
rm Pbp_MB3_*.root
rm Pbp_MB3A_*.root
rm Pbp_MB3_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB4/vnanal* > Pbp_MB4_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB4_list.dat","Pbp_MB4_",10)'
rm Pbp_MB4_list.dat
ls -1 Pbp_MB4_*.root > Pbp_MB4_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB4_list.dat","Pbp_MB4A_",10)'
hadd -f Pbp_MB4.root Pbp_MB4A_*.root
rm Pbp_MB4_*.root
rm Pbp_MB4A_*.root
rm Pbp_MB4_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB5/vnanal* > Pbp_MB5_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB5_list.dat","Pbp_MB5_",10)'
rm Pbp_MB5_list.dat
ls -1 Pbp_MB5_*.root > Pbp_MB5_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB5_list.dat","Pbp_MB5A_",10)'
hadd -f Pbp_MB5.root Pbp_MB5A_*.root
rm Pbp_MB5_*.root
rm Pbp_MB5A_*.root
rm Pbp_MB5_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB6/vnanal* > Pbp_MB6_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB6_list.dat","Pbp_MB6_",10)'
rm Pbp_MB6_list.dat
ls -1 Pbp_MB6_*.root > Pbp_MB6_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB6_list.dat","Pbp_MB6A_",10)'
hadd -f Pbp_MB6.root Pbp_MB6A_*.root
rm Pbp_MB6_*.root
rm Pbp_MB6A_*.root
rm Pbp_MB6_list.dat

ls -1 /rfs/sanders/crab_projects/crab_pPb2016_Pbp_MB7/vnanal* > Pbp_MB7_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB7_list.dat","Pbp_MB7_",10)'
rm Pbp_MB7_list.dat
ls -1 Pbp_MB7_*.root > Pbp_MB7_list.dat
root -l -b -q 'MergeFiles.C+("Pbp_MB7_list.dat","Pbp_MB7A_",10)'
hadd -f Pbp_MB7.root Pbp_MB7A_*.root
rm Pbp_MB7_*.root
rm Pbp_MB7A_*.root
rm Pbp_MB7_list.dat
