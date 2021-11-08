### Spectral balance analysis of the PC1 eigenvectors

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

source ("config.r")

load.pkgs ("emmeans")

load (file.path (results.dir, "pc1.dat"))
mssq <- read.csv (file.path (data.dir, "mssq.csv"), sep = ";", dec = ",")
invalid <- as.character (mssq$nom [which (! is.na (mssq$invalide))])

system2 ("make", c ("-C", figures.dir, "clean-individual-plots"))

problematic.cases <- c ()

(subjects <- vis.cond <- freq.cond <- mssq.vals <- dmpf  <- c ())

line.type <- list (f1 = 1, f2 = "21", f4 = "12")
labels.legend <- c ("0.1 Hz", "0.2 Hz", "0.4 Hz")
line.type.legend <- c ("F1", "21", "12")
line.col <- list (yf = "red", yo = "blue")

rot.all <- c ()

for (subject in names (pc1)) {

    if (subject %in% invalid)
        next

    mssq.subj <- mssq$MSSQ [which (mssq$nom == subject)]

    pdf (file = file.path (figures.dir, sprintf ("%s-spec.pdf", subject)))
    start.plot <- TRUE

    subj.data <- pc1 [[subject]]
    for (vision in names (subj.data)) {
        vis.data <- subj.data [[vision]]
        for (frequency in names (vis.data)) {
            cat (sprintf ("\rsubject: %s   vision: %s   frequency: %s",
                          subject, vision, frequency))
            flush (stdout ())

            freq.data <- vis.data [[frequency]]

            rot <- freq.data$rotation * freq.data$sdev [1]

            if (mean (rot) < 0)
                rot <- -rot

            rot.all <- rbind (rot.all, rot)

            spec.sup <- freq.data$center + rot
            spec.inf <- freq.data$center - rot
            spec.mid <- freq.data$center

            freq <- seq (1 / window.dur, by = 1 / window.dur,
                         length.out = length (spec.sup))
            stem <- sprintf ("%s-%s-%s", subject, vision, frequency)
            pdf (file = file.path (figures.dir, sprintf ("%s-spec.pdf", stem)))
            plot (freq, spec.mid, col = "gray", ylim = c (-60, 100),
                  type = "n", log = "x", las = 1, xlab = "frequency (Hz)",
                  ylab = "amplitude (dB)",
                  main = sprintf ("%s %s %s (MSSQ = %.1f)", subject, vision,
                                  frequency, mssq.subj))
            lines (freq, spec.mid, col = "gray")
            lines (freq, spec.sup, col = "red")
            lines (freq, spec.inf, col = "blue")
            dev.off ()

            ## FIXME: Put frequency bands in config.r
            dmpf <- c (dmpf, mean (rot [13 : 30]) - mean (rot [5 : 12]))

            if (start.plot) {
                plot (freq, rot, las = 1, type = "l", log = "x",
                      xlab = "frequency (Hz)", ylab = "amplitude (dB)",
                      lwd = 3, col = line.col [[vision]],
                      lty = line.type [[frequency]],
                      ylim = c (0, 20),
                      main = sprintf ("%s (MSSQ = %.1f)", subject,
                                      mssq.subj))
                start.plot <- FALSE
            } else
                lines (freq, rot, col = line.col [[vision]],
                       lty = line.type [[frequency]], lwd = 3)

            subjects <- c (subjects, subject)
            vis.cond <- c (vis.cond, vision)
            freq.cond <- c (freq.cond, frequency)
            mssq.vals <- c (mssq.vals, mssq.subj)

        }
    }

    legend ("topleft", inset = 0.05, legend = labels.legend,
            lty = line.type.legend, lwd = 3)
    polygon (c (0.5, 0.5, 1.25, 1.25), c (-10, 40, 40, -10), border = NA,
             col = "#ffd70040")
    polygon (c (1.25, 1.25, 3.0, 3.0), c (-10, 40, 40, -10), border = NA,
             col = "#00800040")
    dev.off ()
}

cat ("\n")
flush (stdout ())

system2 ("make", c ("pc1-subject-spec.pdf",
                    "-f", "../Makefile",
                    "-C", figures.dir))

recovery <- data.frame (subject = subjects, vision = vis.cond,
                        frequency = freq.cond, mssq = mssq.vals, dmpf = dmpf)

fm.dmpf <- lm (dmpf ~ vision * frequency * mssq, recovery)
anova (fm.dmpf)

fm.dmpf.step <- step (fm.dmpf)
anova (fm.dmpf.step)

em.freq <- emmeans (fm.dmpf.step, ~ frequency)
contrast (em.freq, method = list ("f1 vs. f2" = c (1, -1, 0),
                                  "f4 vs. (f1+f2)/2" = c (1/2, 1/2, -1)))

col.vis <- c ("blue", "red")
bg.col.vis <- c ("white", "red")
pch.vis <- c (21, 23)

plot.spectra <- function (subject, ylim = c (2, 12), ylab = "amplitude (dB)",
                          show.freq.legend = TRUE, show.eyes.legend = FALSE) {

    start.plot <- TRUE

    mssq.subj <- mssq$MSSQ [which (mssq$nom == subject)]

    subj.data <- pc1 [[subject]]
    for (vision in names (subj.data)) {
        vis.data <- subj.data [[vision]]
        for (frequency in names (vis.data)) {

            freq.data <- vis.data [[frequency]]

            rot <- freq.data$rotation * freq.data$sdev [1]

            if (mean (rot) < 0)
                rot <- -rot

            ## FIXME: Put frequency bands in config.r
            dmpf <- c (dmpf, mean (rot [13 : 30]) - mean (rot [5 : 12]))

            par (mar = c (5, 4, 3, 0.5))
            if (start.plot) {
                plot (freq, rot, las = 1, type = "l", log = "x",
                      xlab = "frequency (Hz)", ylab = ylab,
                      lwd = 2, col = line.col [[vision]],
                      lty = line.type [[frequency]],
                      ylim = ylim,
                      main = sprintf ("MSSQ = %.1f", mssq.subj))
                start.plot <- FALSE
            } else
                lines (freq, rot, col = line.col [[vision]],
                       lty = line.type [[frequency]], lwd = 2)

            subjects <- c (subjects, subject)
            vis.cond <- c (vis.cond, vision)
            freq.cond <- c (freq.cond, frequency)
            mssq.vals <- c (mssq.vals, mssq.subj)

        }
    }

    if (show.freq.legend)
        legend ("topleft", inset = 0.05, legend = labels.legend,
                lty = line.type.legend, lwd = 2)
    polygon (c (0.5, 0.5, 1.25, 1.25), c (-10, 40, 40, -10), border = NA,
             col = "#ffd70040")
    polygon (c (1.25, 1.25, 3.0, 3.0), c (-10, 40, 40, -10), border = NA,
             col = "#00800040")

    if (show.eyes.legend)
        legend ("bottomright", inset = 0.05, pch = 15, bg = "white",
                col = col.vis, legend = c ("open eyes", "closed eyes"))

}

plot.dmpf.mssq <- function (mssq.level, plot.legend = TRUE, ylim = c (-2, 1),
                            ylab = "spectral balance (dB)") {
    pred.dmpf <- predict (fm.dmpf,
                          newdata = data.frame (vision = rep (c ("yo", "yf"), 3),
                                                frequency = rep (c ("f1", "f2", "f4"),
                                                                 c (2, 2, 2)),
                                                mssq = rep (mssq.level, 6)),
                          se.fit = TRUE)
    par (mar = c (3, 4, 0.5, 0.5))
    plot (NA, NA, type = "n", xlim = c (0.5, 6.5), ylim = ylim,
          las = 1, bty = "n", xaxt = "n", xlab = "", ylab = ylab)
    for (i in seq (1, 3))
        lines (c (2 * i - 1, 2 * i), pred.dmpf$fit [(2 * i - 1) : (2 * i)],
               lwd = 2, lty = "31", col = "gray")
    for (i in seq (1, 6))
        lines (rep (i, 2), pred.dmpf$fit [i] + c (-1, 1) * pred.dmpf$se.fit [i],
               col = col.vis [floor (i - 1) %% 2 + 1], lwd = 2)
    points (seq (1, 6), pred.dmpf$fit, pch = rep (pch.vis, 3), cex = 2,
            col = rep (col.vis, 3), bg = rep (bg.col.vis, 3))
    axis (1, at = seq (1, 6, by = 2) + 0.5, labels = c ("0.1", "0.2", "0.4"))

    if (plot.legend)
        legend ("topleft", inset = 0.05, pch = pch.vis, col = col.vis,
                pt.bg = bg.col.vis, legend = c ("open eyes", "closed eyes"))
}

pdf (file = file.path (figures.dir, "predict-fm-dmpf.pdf"),
     height = 6, width = 6)

layout (matrix (c (1, 2, 3, 4, 5, 5), byrow = TRUE, ncol = 2),
        heights = c (6, 5, 0.3))
plot.spectra ("s033")
plot.spectra ("s032", ylab = "", show.freq.legend = FALSE,
              show.eyes.legend = TRUE)
plot.dmpf.mssq (0)
plot.dmpf.mssq (14, plot.legend = FALSE, ylab = "")
par (mar = rep (0, 4))
plot (0, 0, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
text (0, 0, labels = "perturbation frequency (Hz)", adj = c (0.5, 0.5))

dummy <- dev.off ()

pdf (file = file.path (figures.dir, "mean-pc1-spec.pdf"))
mu <- colMeans (rot.all)
se <- sqrt (colMeans (rot.all ^ 2) - mu ^ 2) / sqrt (nrow (rot.all))
plot (freq, mu, log = "x", type = "l", xlab = "frequency (Hz)",
      ylab = "amplitude (dB)", ylim = c (6, 9))
polygon(c (freq, rev (freq)), c (mu + se, rev (mu - se)), border = NA,
        col = "#00000030")
abline (v = c (0.1, 0.2, seq (0.5, 10, by = 0.5)), col = "#00000040")
dev.off ()
