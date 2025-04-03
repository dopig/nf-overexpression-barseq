#!/bin/bash

# Get the script's directory
SCRIPT_DIR=$(dirname "$0")
EXT_DIR="$SCRIPT_DIR/../shared/external2"

# Create directories, though wget seems ok if don't
mkdir -p $EXT_DIR/{src,lib,primers}

# Download files to src
wget -P "$EXT_DIR/src" https://bitbucket.org/berkeleylab/feba/raw/0975b4c3239d35d041fa954fbc359e7f0cecea88/bin/MultiCodes.pl
wget -P "$EXT_DIR/lib" https://bitbucket.org/berkeleylab/feba/raw/0975b4c3239d35d041fa954fbc359e7f0cecea88/lib/FEBA_Utils.pm
wget -P "$EXT_DIR/src" https://raw.githubusercontent.com/morgannprice/BobaseqFitness/refs/heads/main/bobaseq.R

# Download file to primers
wget -P "$EXT_DIR/primers" https://bitbucket.org/berkeleylab/feba/raw/0975b4c3239d35d041fa954fbc359e7f0cecea88/primers/barseq4.index2
