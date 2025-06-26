#!/usr/bin/env Rscript

# A script to setup dependencies using renv. It can be used to install
# dependencies when the renv.lock file is missing or out-of-date. Inspired by
# the python poetry tool chain.
#
# By default if there are discrepencies between the human defined dependencies
# and the renv.lock then an instructional error message witll be printed. But
# this can be bypassed with --force.
#
# Run with --help for usage information.
#
# Currently supports: dependencies.R file format. Does not support: DESCRIPTION
# file format.

#############
# CONSTANTS #
#############

SCRIPT_VERSION <- "0.1.0"
MISSING_VERSION <- "MISSING_VERSION"
