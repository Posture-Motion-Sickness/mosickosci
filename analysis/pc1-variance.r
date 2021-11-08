### Descriptive analysis of the PC1 score

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

### * Load the global variables
source ("config.r")

### * Load the PCA data
load (file.path (data.dir, "pc1.dat"))

### * Initialize output vector
pc1.var <- c ()

### * Loop over subjects
for (s in names (pc1)) {

    ## ** Loop over vision conditions
    for (v in names (pc1 [[s]])) {

        ## ** Loop over frequency conditions
        for (f in names (pc1 [[s]] [[v]])) {

            ## *** Store the relative variance of PC1
            sdev <- pc1 [[s]] [[v]] [[f]] $ sdev
            pc1.var <- c (pc1.var, sdev)

        }

    }
}

### * Descriptive statistics
show (summary (pc1.var))

### * Plot the PC1 % variance histogram
pdf (file = file.path (figures.dir, "pc1-var-hist.pdf"))
hist (pc1.var, main = "", las = 1, xlab = "variance of 1st PC",
      xlim = c (0, 1), ylab = "count")
dev.off ()
