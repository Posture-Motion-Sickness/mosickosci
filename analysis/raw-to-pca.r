### Compute the principal component analysis on the PSD from the raw data

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


### * Load the configuration file
source ("config.r")

### * Load needed libraries (for function hamming and PSD analysis,
###   and parallel computation)
load.pkgs (c ("signal", "psd", "doSNOW"))

### * Set analysis variables
### ** Number of samples in the analysis
nb.samp <- floor (3 * phase.dur * samp.freq)
### ** Width of the FFT analysis window
window.size <- round (window.dur * samp.freq) + 1
### ** Window profile
window <- hamming (window.size)
### ** Initial and final sample indices
initial.idx <- max (c (2, floor (window.size * min.freq / samp.freq)))
final.idx <- floor (window.size * max.freq / samp.freq) + 1
### ** Number of steps for FFT analysis
nb.step <- nb.samp - window.size - 1
### ** Temporary matrix for storing the spectra
spectra <- matrix (NA, nrow = nb.step, ncol = final.idx - initial.idx + 1)

### * Get the list of files with the posturographic data
files <- list.files (data.dir)
files <- files [grep ("^s\\d\\d\\d-y.*", files)]

### * Get invalid subjects
mssq <- read.csv (file.path (data.dir, "mssq.csv"), sep = ";", dec = ",")
invalid <- as.character (mssq$nom [which (! is.na (mssq$invalide))])

### * Register the parallel cluster
cl <- makeCluster (num.cores)
registerDoSNOW (cl)

### * Progress bar
iterations <- length (files)
pb <- txtProgressBar (max = iterations, style = 3)
progress <- function (n) setTxtProgressBar (pb, n)

### * Loop over the data files
results <- foreach (f = files,
                    .options.snow = list (progress = progress),
                    .packages = "psd") %dopar% {

    ## ** Get the parts of the file name
    terms <- strsplit (f, "-")
    subject <- terms [[1]] [1]

    ## ** Skip if subject is invalid
    if (! subject %in% invalid) {

        vision <- terms [[1]] [2]
        frequency <- strsplit (terms [[1]] [3], "\\.") [[1]] [1]

        ## ** Get the posturography data
        data <- read.table (file.path (data.dir, f), sep = "\t",
                            skip = 1, header = TRUE)
        y <- data [, 8]

        ## ** Compute the global PSD
        psd <- psdcore (y, ntaper = nb.tapers, X.frq = samp.freq)

        ## ** Compute the spectra
        for (t in seq (1, nb.step)) {
            idx <- seq (t, t + window.size)
            p <- psdcore (y [idx], ntaper = nb.tapers, X.frq = samp.freq)
            spec.freq <- p$freq [initial.idx : final.idx]
            spec <- 10 * log10 (p$spec)
            spectra [t, ] <- spec [initial.idx : final.idx]

        }

        ## ** Do the PC analysis
        pc <- prcomp (spectra)

        ## ** Extract the first PC
        pc1 <- list (x = pc$x [, 1],
                     center = pc$center,
                     rotation = pc$rotation [, 1],
                     sdev = pc$sdev)

        list (subject = subject,
              vision = vision,
              frequency = frequency,
              spec.freq = spec.freq,
              pc1 = pc1,
              psd = psd)
    }

}

### * Stop parallel cluster
stopCluster (cl)

### * Create a list for storing the PC1 and PSD data
pc1 <- psd <- list ()

### * Loop over the results
for (i in seq (1, length (results))) {

    ## ** Get values
    ri <- results [[i]]
    if (is.null (ri))
        next
    subject <- ri$subject
    vision <- ri$vision
    frequency <- ri$frequency
    spec.freq <- ri$spec.freq

    ## ** Store data
    if (! subject %in% names (pc1))
        pc1 [[subject]] <- list ()
    if (! vision %in% names (pc1 [[subject]]))
        pc1 [[subject]] [[vision]] <- list ()
    pc1 [[subject]] [[vision]] [[frequency]] <- ri$pc1

    if (! subject %in% names (psd))
        psd [[subject]] <- list ()
    if (! vision %in% names (psd [[subject]]))
        psd [[subject]] [[vision]] <- list ()
    psd [[subject]] [[vision]] [[frequency]] <- ri$psd

}

## * Store the PCA data
save (pc1, psd, spec.freq, file = file.path (results.dir, "pc1.dat"))
