#!/bin/sh
ls -1 /rfs/sanders/crab_projects/crab_pPb2016_pPb_HM250_7/vnanal* > pPb_HM250_7_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM250_7_list.dat","pPb_HM250_7_",10)'
rm pPb_HM250_7_list.dat
ls -1 pPb_HM250_7_*.root > pPb_HM250_7_list.dat
root -l -b -q 'MergeFiles.C+("pPb_HM250_7_list.dat","pPb_HM250_7A_",10)'
hadd -f pPb_HM250.root pPb_HM250_7A_*.root
rm pPb_HM250_7_*.root
rm pPb_HM250_7A_*.root
rm pPb_HM250_7_list.dat
