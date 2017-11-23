#!/bin/sh
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_262548_262799/* > inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_262800_263230/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_263231_263359/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_263360_263379/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_263380_263614/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_263615_263757/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_noRecenter_262548_262799/* > inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_noRecenter_262800_263230/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_noRecenter_263231_263359/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_noRecenter_263360_263379/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_noRecenter_263380_263614/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_MH_noRecenter_263615_263757/* >> inlist.dat
root -l -b -q ClearBadFiles.C+

