% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_genespace.R
\name{init_genespace}
\alias{init_genespace}
\title{Find files and directories for a GENESPACE run}
\usage{
init_genespace(
  genomeIDs,
  wd,
  outgroup = NULL,
  speciesIDs,
  versionIDs,
  rawGenomeDir,
  ploidy,
  diamondMode = "more-sensitive",
  orthofinderMethod = "default",
  orthofinderInBlk = any(ploidy > 1),
  nCores = detectCores()/2,
  minPepLen = 20,
  overwrite = T,
  path2orthofinder = "orthofinder",
  path2mcscanx = "MCScanX",
  path2diamond = "diamond",
  gffString = "gff",
  pepString = "pep|prot",
  verbose = TRUE
)
}
\arguments{
\item{genomeIDs}{character vector of length > 1, matching length
of speciesIDs, versions and ploidy. Specifies the name to assign
to each genome. This vector must be unique and can be any string
that begins with a letter (a-z, A-Z). An X will be appended to
IDs starting with special characters or numbers.}

\item{wd}{file.path where the analysis will be run}

\item{outgroup}{character string matching one of the genomeIDs that will
be used in the orthofinder -og run but not in the synteny search. Suggested
to ensure that there is an outgroup that predates any WGD that the user
would like to study.}

\item{speciesIDs}{file path character vector. This is the subdirectory in
rawGenomeDir containing the files for each genomeID.}

\item{versionIDs}{file path character vector. This is the subdirectory in
each speciesID with the genome version that contains the files for
each genomeID.}

\item{rawGenomeDir}{single file path, pointing to the directory that
contains all the genome annotation raw files.}

\item{ploidy}{integer string specifying ploidy of genome assemblies. This is
usually half of the actual ploidy, that is an inbred diploid usually is
represented by a haploid genome assembly.}

\item{diamondMode}{character string to match one of the diamond search modes}

\item{orthofinderMethod}{character string either 'fast' or 'default'. See
build_OFdb for details.}

\item{orthofinderInBlk}{logical, should orthofinder be re-run within
syntenic regions? Highly recommended for polyploids}

\item{nCores}{integer of length 1 specifying the number of parallel processes
to run}

\item{minPepLen}{integer, the shortest peptide to be analyzed}

\item{overwrite}{logical, should existing directories be overwritten?}

\item{path2orthofinder}{character string coercible to a file path that points
to the orthofinder executable. If orthofinder is in the path, specify with
"orthofinder"}

\item{path2mcscanx}{see path2orthofinder, except to the mcscanx directory.
This must contain the MCScanX_h folder.}

\item{path2diamond}{character string coercible to a file path that points
to the diamond executable. If diamond is in the path, specify with
"diamond"}

\item{gffString}{character string of length 1, coercable to a
regular expression. File names in the 'annotation' subdirectories are
screened for this string. If a single, unique file is found, it is returned
as the gff source file. Otherwise an error is produced.}

\item{pepString}{character string of length 1, coercable to a
regular expression. See gffSearchString.}

\item{verbose}{logical length 1, should updates be printed to the console?}
}
\value{
A list containing paths to the raw files. If a file is not found,
path is returned as null and a warning is printed.
}
\description{
\code{init_genespace} Searches for desired genome files in the
raw genome repo director.
}
\details{
Simple directory parser to find and check the paths to all
annotation and assembly files.
}
\examples{
\dontrun{
runwd <- file.path(getwd(), "testGenespace")
make_exampleDataDir(writeDir = runwd)

#############################################################################
###### Default method, can run orthofinder outside of R on a cluster etc. ###
gpar <- init_genespace(
  genomeIDs = c("human","chimp","rhesus"),
  speciesIDs = c("human","chimp","rhesus"),
  versionIDs = c("human","chimp","rhesus"),
  ploidy = rep(1,3),
  wd = runwd,
  nCores = 4,
  gffString = "gff",
  pepString = "pep",
  path2orthofinder = "orthofinder",
  path2mcscanx = "~/MCScanX",
  rawGenomeDir = file.path(runwd, "rawGenomes"))

#############################################################################
###### Slower but more accurate inference  ##################################
###### especially for distantly related genomes or ploidy > 1 ###############
# add in orthofinder in block orthogroup ...
gpar <- init_genespace(
  genomeIDs = c("human","chimp","rhesus"),
  speciesIDs = c("human","chimp","rhesus"),
  versionIDs = c("human","chimp","rhesus"),
  ploidy = rep(1,3),
  wd = runwd,
  nCores = 4,
  orthofinderInBlk = TRUE,
  gffString = "gff",
  pepString = "pep",
  path2orthofinder = "orthofinder",
  path2mcscanx = "~/MCScanX",
  rawGenomeDir = file.path(runwd, "rawGenomes"))

#############################################################################
###### Fast visualization of closely related genome #########################
# just a quick visualization ... not recommended for a full run, although
# probably not much less sensitive with closely related genomes
# **NOTE** this method can only be done when orthofinder is installed in the
# path for this R session (e.g. opening R from a conda env with orthofinder)
# **NOTE** increasing minPepLen results in a loss of fidelity, but an
# increase in speed.

# This is the parameterization used in the help files.

gpar <- init_genespace(
  genomeIDs = c("human","chimp","rhesus"),
  speciesIDs = c("human","chimp","rhesus"),
  versionIDs = c("human","chimp","rhesus"),
  ploidy = rep(1,3),
  diamondMode = "fast",
  orthofinderMethod = "fast",
  wd = runwd,
  nCores = 4,
  minPepLen = 50,
  gffString = "gff",
  pepString = "pep",
  path2orthofinder = "orthofinder",
  path2mcscanx = "~/MCScanX",
  rawGenomeDir = file.path(runwd, "rawGenomes"))

}
}
