

#header #################################################################################
#'Main.R'

#Title: Salmonella transfer
#Project ID: pid
#Client: UFRGS
#Author: <Eduardo> <Costa>, Wageningen Bioveterinary Research

#Description: This is a paper describying the a Bayesian model to estimate the *Salmonella* sp. transfer probability between knife and pork in a household scenario.

#Start date: date
#Last Update: {6:date}

#R version: r.version
#Scriptversion: version

#Dependencies
#<-Downstream simula_final.R
#->Upstream none

#Input:
#- None

#Output:
#- pork_to_knife.txt; knife_to_pork.txt

#Peer reviewer(s)

#Please ensure directories are relative. Please comment code sufficiently.

#Script start#############################################################################



#Installing and loading Packages
#Install packages

#Packages to be used
packages<-c("readxl","here","tidyverse","ggplot2","gridExtra","knitr","BRugs","coda","rjags","psych")


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}


# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


#Creating directories

dir.create(here("Figures"))
dir.create(here("Output"))



#Running the model
source(here("Script","simula_final.R"))
