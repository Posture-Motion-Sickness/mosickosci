### General configuration and utility functions

### This program is part of mosickosci
###
### Copyright (C) 2021  Rafael Laboissière
### Copyright (C) 2021  Merrick Dida
###
### This program is free software: you can redistribute it and/or modify it
### under the terms of the GNU General Public License as published by the
### Free Software Foundation, either version 3 of the License, or (at your
### option) any later version.
###
### This program is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License along
### with this program.  If not, see <http://www.gnu.org/licenses/>.

### * Utility function

### ** Load R packages, installing them if necessary
load.pkgs <- function (list.of.pkgs)
    for (pkg in list.of.pkgs)
        if (! require (pkg, character.only = TRUE)) {
            install.packages (pkg)
            require (pkg, character.only = TRUE)
        }

### ** Force creation of directories
force.dir.create <- function (dir)
    suppressWarnings (dir.create (dir, recursive = TRUE))

### ** Banner for section title in output text file
banner <- function (title, chr) {
    n <- nchar (title)
    s <- paste (rep (chr, n), collapse = "")
    cat (sprintf ("%s\n%s\n%s\n\n", s, title, s))
}

### ** Read brute data file
read.brute <- function (file.name)
    read.table (file.name, sep = "\t", skip = 1, header = TRUE)

### ** Compute standard error of vector
se <- function (x)
    sd (x) / sqrt (length (x))

### ** Add transparency to color
add.transparency <- function (color, level)
    sapply (color, function (x) {
                       rgb <- col2rgb (x)
                       sprintf ("#%02x%02x%02x%02x", rgb [1], rgb [2], rgb [3],
                                as.integer (256 * level))
                   })

### * Global variables used in the analysis scripts

### ** Do configuration function
do.config <- function (cfg) {

    for (n in names (cfg))
        .GlobalEnv [[n]] <- cfg [[n]]

    ## *** Set Directories
    ## **** Parameter-defined sub-directory
    subdir <- sprintf ("%dtap-%.1fHz-%.1fs", nb.tapers, max.freq, window.dur)
    ## **** Figures
    .GlobalEnv$figures.dir <- file.path ("..", "figures", subdir)
    force.dir.create (.GlobalEnv$figures.dir)
    ## **** Results
    .GlobalEnv$results.dir <- file.path ("..", "results", subdir)
    force.dir.create (.GlobalEnv$results.dir)

}

## ** Raw data direcotory
data.dir <- file.path ("..", "data")
force.dir.create (data.dir)

### ** Minimum frequency of the spectral analysis
min.freq <- 0

### ** Sampling frequency (in Hz)
samp.freq <- 100

### ** Duration of each phase (static/oscillation/static, in sec)
phase.dur <- 60

### ** Oscillation frequencies
osci.freq <- list (f1 = 0.1, f2 = 0.2, f4 = 0.4)

### ** Default values

load.pkgs ("parallel")

do.config (
    list (

        ## *** FFT window duration (in sec)
        window.dur = 10,

        ## *** Maximum frequency of the spectral analysis
        max.freq = 5,

        ## *** Number of tapers for PSD analysis
        nb.tapers = 5,

        ## *** Number of cores for paralle computing
        num.cores = detectCores ()

    )
)

### ** Eigenvalue normalization parameters

### The following values are in seconds, These vectors determine the begin and
### the end of the intervals used in the normalization procedure
norm.base.interval <- phase.dur - c (phase.dur / 2, window.dur / 2)
norm.pert.interval <- 2 * phase.dur - c (phase.dur / 2, window.dur / 2)
### Normalized values (first: baseline, second: pertrubation)
norm.vals <- c (0, 1)

### ** Decimation factor for exponential model fitting
decimate.factor <- 100

### ** Level for confidence intervals
conf.level <- 0.95
