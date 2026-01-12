#!/bin/bash

pymol -r ligand.pdb -d "
load ligand.pdb;
h_add;
save ligand.mol2;
quit;
"

