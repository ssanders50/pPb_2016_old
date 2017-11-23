#!/bin/sh
rm inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_262548_262799/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_262800_263230/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_263231_263359/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_263360_263379/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_263380_263614/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_263615_263757/* >> inlist.dat
root -l -b -q 'MergeFiles.C("inlist.dat","h_",20)'
rm inlist.dat
ls -1 h_*.root >> inlist.dat
root -l -b -q 'MergeFiles.C("inlist.dat","hh_",5)'
rm MH.root
hadd MH.root hh_*.root
 