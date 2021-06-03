#!/usr/bin/env bash

wget -c --recursive --no-parent ftp://iprg_study:ABRF329@ftp.peptideatlas.org/distro/mzML/
wget -c --recursive --no-parent ftp://iprg_study:ABRF329@ftp.peptideatlas.org/distro/fasta/
rm ftp.peptideatlas.org/distro/fasta/iPRG2015.TargDecoy.fasta 
 
