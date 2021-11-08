### Analysis of the temporal evolution of the PC1 score

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

### * Load configuration parameters
source ("config.r")

### * Load need libraries
load.pkgs ("Cairo")

### * Normalization function
normalize <- function (x) {
    idx.base <- seq (samp.freq * (norm.base.interval [1] - window.dur / 2),
                     samp.freq * (norm.base.interval [2] - window.dur / 2))
    idx.pert <-seq (samp.freq * (norm.pert.interval [1] - window.dur / 2),
                    samp.freq * (norm.pert.interval [2] - window.dur / 2))
    mean.base <- mean (x [idx.base])
    mean.pert <- mean (x [idx.pert])
    return (norm.vals [1]
            + (norm.vals [2] - norm.vals [1])
              * (x - mean.base) / (mean.pert - mean.base))
}

### * Load PC1 data
load (file.path (results.dir, "pc1.dat"))

### * Reorganize data into a data frame
cat ("Reorganize data...\n")

### ** Initial (empty) data frame
score.df <- data.frame (subject = character (),
                        vision = character (),
                        frequency = character (),
                        time = numeric (),
                        pc1 = numeric (),
                        stringsAsFactors = FALSE)
names.score.df <- names (score.df)

### ** Loop over subjects of PC1 data
for (s in names (pc1))

    ## * Loop over vision conditions
    for (v in names (pc1 [[s]]))

        ## * Loop over frequency conditions
        for (f in names (pc1 [[s]] [[v]])) {

            ## ** Get score
            x <- normalize (pc1 [[s]] [[v]] [[f]]$x)
            n <- length (x)
            ## ** Associated time values
            time <- seq (window.dur / 2, length.out = n, by = 1 / samp.freq)

            ## ** Compose data frame
            score.df <- rbind (score.df, data.frame (rep (s, n), rep (v, n),
                                                 rep (f, n), time, x))

            ## ** Progress meter
            cat (sprintf ("\rSubject: %s  vision: %s  frequency: %s", s, v, f))
            flush (stdout ())
        }

### *** Clear the progress meter
cat ("\n")
flush (stdout ())

### ** Rename columns and make them factors
names (score.df) <- names.score.df
for (n in names (score.df))
    if (class (score.df [[n]]) == "character")
        score.df [[n]] <- as.factor (score.df [[n]])

### ** Add indicators for vision and frequency levels
score.df$is.ec <- c (0, 1) [as.factor (score.df$vision == "yf")]
score.df$is.f2 <- c (0, 1) [as.factor (score.df$frequency == "f2"
                                       | score.df$frequency == "f4")]
score.df$is.f4 <- c (0, 1) [as.factor (score.df$frequency == "f4")]

### * Add order of trial
cat ("Add order of trial...\n")

order <- read.csv (file.path (data.dir, "session-order.csv"))
order <- order [which (order$subject %in% names (pc1)), ]

score.df$order <- rep (NA, nrow (score.df))
for (s in seq (1, nrow (order)))

    for (j in seq (1, 6)) {

        vf <- strsplit (order [s, j + 1], "-")
        idx <- which (score.df$subject == order$subject [s]
                      & score.df$vision == vf [[1]] [1]
                      & score.df$frequency == vf [[1]] [2])
        score.df$order [idx] <- j

        cat (sprintf ("\rSubject: %s", order$subject [s]))
        flush (stdout ())

    }

cat ("\n")
flush (stdout ())

score.df$order <- factor (score.df$order)
score.df$num.order <- as.numeric (score.df$order) - 1

### * Plots

### ** Average time course of the score

### *** Compute grand mean and standard error values
mean.score <- aggregate (pc1 ~ time, score.df, mean)
se.score <- aggregate (pc1 ~ time, score.df, se)$pc1

### *** Confidence interval of mean estimation
ci.level <- 0.95
t.crit <- qt ((1 - ci.level) / 2, nrow (score.df) / nrow (mean.score) - 1)
m.minus.ci <- mean.score$pc1 - t.crit * se.score
m.plus.ci <- mean.score$pc1 + t.crit * se.score

### *** Get indices for each phase
phase.idx <- list (
    base = which (mean.score$time < phase.dur - window.dur / 2),
    pert = which (mean.score$time > phase.dur + window.dur / 2
                  & mean.score$time < 2 * phase.dur - window.dur / 2),
    ret = which (mean.score$time > 2 * phase.dur + window.dur / 2))

### *** Figure label letters, poisition and size
figure.labels <- c ("A", "B", "C", "D")

initial.plot <- function (n, cols = NA, labels = NA, x.label = "time (s)",
                          y.label = "normalized PC1 score",
                          label.x = 5, label.y = 1.1, label.cex = 2,
                          legend.x = 35, legend.y = 0.9) {
    plot (NA, NA, type = "l", las = 1, bty = "n",
          ylab = y.label, xlab = x.label,
          xaxt = "n", xlim = c (0, 180),
          ylim = c (-0.05, 1.2))
    axis (1, at = seq (0, 180, by = 30))
    polygon (c (phase.dur + (window.dur / 2) * c (-1, 1, 1, -1)),
             c (-1, -1, 2, 2), border = NA, col = "#00000020")
    polygon (c (2 * phase.dur + (window.dur / 2) * c (-1, 1, 1, -1)),
             c (-1, -1, 2, 2), border = NA, col = "#00000020")
    text (label.x, label.y, labels = figure.labels [n], adj = c (0, 1),
          cex = label.cex)
    if (! any (is.na (labels)))
        legend (legend.x, legend.y, col = cols, lwd = 2, legend = labels,
                bg = "white")
}

fig.scale <- 0.8
CairoPDF (file = file.path (figures.dir, "score-pc1-mean.pdf"),
          width = fig.scale * 9, height = fig.scale * 7)
par (mar = c (5, 4, 0.1, 0.1))

layout (matrix (c (1, 2, 3, 4), nrow = 2, byrow = TRUE))

initial.plot (1, x.label = "")
for (n in names (phase.idx)) {
    i <- phase.idx [[n]]
    t <- mean.score$time [i]
    lines (t, mean.score$pc1 [i])
    polygon (c (t, rev (t)), c (m.plus.ci [i], rev (m.minus.ci [i])),
             col = "#00000040", border = NA)
}

pc1.freq <- list ()
for (f in levels (score.df$frequency))
    pc1.freq [[f]] <- aggregate (pc1 ~ time,
                                 subset (score.df, frequency == f),
                                 mean)$pc1

freq.cols <- add.transparency (c ("blue", "red", "green3"), 0.5)

initial.plot (2, freq.cols, c ("0.1 Hz", "0.2 Hz", "0.4 Hz"),
              x.label = "", y.label = "")
for (n in names (phase.idx)) {
    i <- phase.idx [[n]]
    t <- mean.score$time [i]
    for (f in seq (1, length (pc1.freq)))
        lines (t, pc1.freq [[f]] [i], col = freq.cols [f], lwd = 2)
}

pc1.vis <- list ()
for (v in levels (score.df$vision))
    pc1.vis [[v]] <- aggregate (pc1 ~ time,
                                subset (score.df, vision == v),
                                mean)$pc1

vis.cols <- add.transparency (c ("red", "blue"), 0.5)

initial.plot (3, vis.cols, c ("EC", "EO"))
for (n in names (phase.idx)) {
    i <- phase.idx [[n]]
    t <- mean.score$time [i]
    for (v in seq (1, length (pc1.vis)))
        lines (t, pc1.vis [[v]] [i], col = vis.cols [v], lwd = 2)
}

pc1.order <- list ()
for (o in levels (score.df$order))
    pc1.order [[o]] <- aggregate (pc1 ~ time,
                                  subset (score.df, order == o),
                                  mean)$pc1

order.cols = add.transparency (c ("red", "orange", "gold",
                                  "green4", "cyan2", "cyan4"), 0.5)

initial.plot (4, order.cols, seq (1, 6), y.label = "")
for (n in names (phase.idx)) {
    i <- phase.idx [[n]]
    t <- mean.score$time [i]
    for (o in seq (1, length (pc1.order)))
        lines (t, pc1.order [[o]] [i], col = order.cols [o], lwd = 2)
}

dev.off ()

df.phase <- list ()
phase.names <- c ("baseline", "perturbation", "return")

for (i in seq (1, 3)) {
    df <- subset (score.df, time >= window.dur / 2  + phase.dur * (i - 1)
                             & time <= phase.dur * i - window.dur / 2)
    idx <- which ((samp.freq * df$time) %% decimate.factor == 0)
    df.phase [[phase.names [i]]] <- df [idx, ]
}

Exp.frm <- ~ (A +
              (B + is.ec * dBec + is.f2 * dBf1f2 + is.f2 * is.f4 * dBf2f4
                 + num.order * dBorder - A) *
              exp (- (t - t0)
                   * (C
                       + is.ec * dCec
                       + is.f2 * dCf1f2
                       + is.f2 * is.f4 * dCf2f4
                       + num.order * dCorder)))

Exp <- deriv (Exp.frm,
              namevec = c ("A",
                           "B", "dBec", "dBf1f2", "dBf2f4", "dBorder",
                           "C", "dCec", "dCf1f2", "dCf2f4", "dCorder"),
              function.arg = c ("t", "is.ec", "is.f2", "is.f4", "num.order",
                                "A",
                                "B", "dBec", "dBf1f2", "dBf2f4", "dBorder",
                                "C", "dCec", "dCf1f2", "dCf2f4", "dCorder"))

fm.nls <- cf <- cf.ci <- list ()

for (i in seq (1, 3)) {

    ph <- phase.names [i]
    cat (sprintf ("===== Phase: %s\n", ph))

    t0 <- phase.dur * (i - 1) + window.dur / 2

    startvec <- c (A = ifelse (i == 2, 1, 0),
                   B = ifelse (i == 2, 1.5, 0.5),
                   dBec = 0, dBf1f2 = 0, dBf2f4 = 0, dBorder = 0,
                   C = 0.2, dCec = 0, dCf1f2 = 0, dCf2f4 = 0, dCorder = 0)
    fm.nls [[ph]] <- nls (pc1 ~ Exp (time, is.ec, is.f2, is.f4, num.order, A,
                                     B, dBec, dBf1f2, dBf2f4, dBorder,
                                     C, dCec, dCf1f2, dCf2f4, dCorder),
                          df.phase [[ph]],
                          start = startvec)

    cf [[ph]] <- coefficients (fm.nls [[ph]])
    show (cf [[ph]])
    cf.ci [[ph]] <- tryCatch (confint (fm.nls [[ph]]),
                              error = function (cond) {
                                  message (paste ("confint failed with error:\n", cond))
                                  message ("Confidence intervals will assume asymptotic normality")
                                  z <- coef (summary (fm.nls [[ph]]))
                                  q <- qnorm (1 - (1 - conf.level) / 2)
                                  return (cbind (z[, 1] - q * z[, 2],
                                                 z[, 1] + q * z[, 2]))
                              })
    show (cf.ci [[ph]])

}

transform.label <- function (x) {
    x <- gsub ("^d", "Δ", x)
    x <- gsub ("f1f2", "[f.1/f.2]", x)
    x <- gsub ("f2f4", "[f.2/f.4]", x)
    x <- gsub ("ec", "[ec]", x)
    x <- gsub ("order", "[order]", x)
    parse (text = x)
}

CairoPDF (file = file.path (figures.dir, "score-exp-cf-ci.pdf"),
          width = 5, height = 10)
layout (matrix (c (1, 2, 3), ncol = 1))

gray.col <- "#00000040"

for (i in seq (1, 3)) {

    par (mar = c (5, 5, 0, 0.3))

    ph <- phase.names [i]

    nr <- nrow (cf.ci [[ph]])
    plot (0, 0, type = "n", bty = "n", xlab = "value",
          xlim = c (-0.2, 1.3), ylim = c (1, nr),
          xaxt= "n", yaxt = "n", ylab = "")
    axis (1, at = seq (0, 1, by = 0.2))
    axis (2, at = seq (1, nrow (cf.ci [[ph]])), las = 1,
          labels = sapply (rev (row.names (cf.ci [[ph]])),
                           transform.label))
    abline (h = seq (1, nrow (cf.ci [[ph]])), col  = gray.col)
    abline (v = 0, col = gray.col)
    abline (v = 1, col = gray.col)
    for (j in seq (1, nr)) {
        points (cf [[ph]] [j], nr - j + 1, pch = 19, cex = 0)
        lines (c (cf.ci [[ph]] [j, 1], cf.ci [[ph]] [j, 2]),
               rep (nr - j + 1, 2), lwd = 3)
        lines (rep (cf [[ph]] [j], 2), nr - j + 1 + c (-0.15, 0.15), lwd = 2)
    }

}

dummy <- dev.off ()

A <- cf$return [["A"]]
B <- cf$return [["B"]]
dBf1f2 <- cf$return [["dBf1f2"]]
dBf2f4 <- cf$return [["dBf2f4"]]
dBec <- cf$return [["dBec"]]
dBorder <- cf$return [["dBorder"]]
C <- cf$return [["C"]]
dCf1f2 <- cf$return [["dCf1f2"]]
dCf2f4 <- cf$return [["dCf2f4"]]
dCec <- cf$return[["dCec"]]
dCorder <- cf$return [["dCorder"]]

t.begin <- 2 * phase.dur  + window.dur / 2
t.end <- 3 * phase.dur - window.dur / 2
### Use time from data
t <- aggregate (pc1 ~ time,
                subset (score.df, time >= t.begin & time <= t.end),
                mean) $ time

pc1.model <- list (Exp (t, 0, 2/3, 1/3, 1, A, B, dBec, dBf1f2, dBf2f4, dBorder,
                        C, dCec, dCf1f2, dCf2f4, dCorder),
                   Exp (t, 1, 2/3, 1/3, 1, A, B, dBec, dBf1f2, dBf2f4, dBorder,
                     C, dCec, dCf1f2, dCf2f4, dCorder))

lag.score.fig.dir <- file.path (figures.dir, "lag-score")
force.dir.create (lag.score.fig.dir)

df.mssq <- read.csv (file.path (data.dir, "mssq.csv"), sep = ";", dec = ",")
MSSQ <- df.mssq$MSSQ [which (df.mssq$nom %in% levels (score.df$subject))]

subjs <- levels (score.df$subject)
vision.cnd <- c ("yf", "yo")
case.lty <- list (1, "31")
vision.col <- c ("red", "blue")

plot.return.phase <- function (subj, bty = "o", main = NULL,
                               ylab = "normalized PC1 score",
                               draw.legend = TRUE, ylim = NULL) {

    par (mar = c (5, 4, ifelse (is.null (main), 0, 3), 0.5))

    plot (t, pc1.model [[1]], type = "l", las = 1, lwd = 2, bty = bty,
          lty = case.lty [[2]], col = vision.col [2],
          xlim = c (120, 180), ylim = ylim,
          ylab = ylab, xlab = "time (s)",
          main = main)

    lines (t, pc1.model [[2]], lwd = 2, col = vision.col [1],
           lty = case.lty [[2]])

    lag.score <- c ()

    for (j in c (1, 2)) {

        ag <- aggregate (pc1 ~ time, subset (score.df,
                                             subject == subj
                                             & vision == vision.cnd [j]
                                             & time >= t.begin
                                             & time <= t.end),
                         mean)
        lines (ag$time, ag$pc1, col = vision.col [j], lwd = 2,
               lty = case.lty [[1]])

    }

    if (draw.legend)
        legend ("topright", inset = 0.05, lwd = 2,
                lty = rep (c (1, 5), c (2, 2)),
                col = rep (rev (vision.col), 2),
                legend = c ("EO data", "EC data", "EO model", "EC model"))

}

nb.subjs <- length (subjs)

for (i in seq (1, nb.subjs)) {

    CairoPDF (file = file.path (lag.score.fig.dir,
                                sprintf ("%s.pdf", subjs [i])),
              width = 5.5, height = 5.5)

    plot.return.phase (subjs [i],
                       main = sprintf ("subject %s (MSSQ = %.1f)",
                                       subjs [i], MSSQ [i]),
#                       ylim = c (-0.05, 0.6))
                       ylim = c (-0.3, 0.6))

    dummy <- dev.off ()

    cat (sprintf ("\rSubject: %s", subjs [i]))
    flush (stdout ())

}

cat ("\n")
flush (stdout ())

system (sprintf ("pdftk %s/*.pdf cat output %s/lag-score-all.pdf",
                 lag.score.fig.dir, figures.dir))


### MSSQ × vision effect on lag score

subject <- vision <- mssq <- lag.score <- c ()

for (i in rev (seq (1, nb.subjs))) {

    s <- subjs [i]

    for (v in vision.cnd) {

        df.s.v <- subset (score.df, subject == s & vision == v
                                    & time >= t.begin & time <= t.end)

        ag <- aggregate (pc1 ~ time, df.s.v, mean)
        order <- mean (df.s.v$num.order)
        order <- 1

        pred.exp <- Exp (t, ifelse (v == "yo", 0, 1), 2/3, 1/3, order,
                         A, B, dBec, dBf1f2, dBf2f4, dBorder,
                         C, dCec, dCf1f2, dCf2f4, dCorder)

        mssq <- c (mssq, MSSQ [i])
        vision <- c (vision, v)
        d <- ag$pc1 - pred.exp
        lag.score <- c (lag.score, mean (d) / sd (d))
        subject <- c (subject, s)

    }

    cat (sprintf ("\rSubject: %s", subjs [i]))
    flush (stdout ())
}

cat ("\n")
flush (stdout ())

df.lag <- data.frame (subject = subject,
                      mssq = mssq,
                      vision = as.factor (vision),
                      lag.score = lag.score)
df.lag$vision <- relevel (df.lag$vision, "yo")

fm.lag <- lm (lag.score ~ vision * mssq, df.lag)

anova (fm.lag)
confint (fm.lag)

max.mssq <- max (df.lag$mssq)
mssq.seq <- seq (0, max.mssq)
pred <- list ()

for (vis in c ("yo", "yf"))
    pred [[vis]] <- predict (fm.lag,
                             newdata = data.frame (vision = vis,
                                                   mssq = mssq.seq),
                             se.fit = TRUE)

sup.yo.se <- pred$yo$fit + pred$yo$se.fit
sup.yf.se <- pred$yf$fit + pred$yf$se.fit
inf.yo.se <- pred$yo$fit - pred$yo$se.fit
inf.yf.se <- pred$yf$fit - pred$yf$se.fit

vis.idx <- df.lag$vision

hilight.subjects <- c ("s028", "s041")

pdf (file = file.path (figures.dir, "lag.score.mssq.pdf"),
     width = 5, height = 9)

layout (matrix (c (1, 1, 2, 3), nrow = 2, byrow = TRUE),
        heights = c (2.25, 1))

par (mar = c (6, 4, 0, 0))
plot (df.lag$mssq, df.lag$lag.score, las = 1, bty = "n",
      pch = c (21, 23) [vis.idx], col = c ("blue", "red") [vis.idx],
      bg = c ("white", "red") [vis.idx],
      xlab = "MSSQ score", ylab = "lag score", cex = 2)

idx.subj <- which (df.lag$subject %in% hilight.subjects)
points (df.lag$mssq [idx.subj], df.lag$lag.score [idx.subj], cex = 5)

lines (mssq.seq, pred$yo$fit, col = "blue", lwd = 2)
polygon (c (mssq.seq, rev (mssq.seq)), c (sup.yo.se, rev (inf.yo.se)),
         border = NA, col = "#0000ff40")
lines (mssq.seq, pred$yf$fit, col = "red", lwd = 2)
polygon (c (mssq.seq, rev (mssq.seq)), c (sup.yf.se, rev (inf.yf.se)),
         border = NA, col = "#ff000040")
abline (h = 0, col = "#00000080", lwd = 2)
legend ("topright", inset = 0.05, legend = c ("open eyes", "closed eyes"),
        col = c ("blue", "red"), lty = 1, lwd = 2)

plot.return.phase (hilight.subjects [1], ylim = c (-0.05, 0.6), bty = "n",
                   draw.legend = FALSE)

plot.return.phase (hilight.subjects [2], ylim = c (-0.05, 0.6), bty = "n",
                   ylab = "")

dummy <- dev.off ()

#################################3


df.mal.de.terre <- subset (score.df, time >= t.begin &  time <= t.end)
df.agg <- aggregate (pc1 ~ is.ec * is.f2 * is.f4, df.mal.de.terre, mean)

t <- seq (t0, t0 + 50, by = 0.5)
pred.exp <- c ()

for (i in seq (1, nrow (df.agg))) {
    pred.exp <- rbind (pred.exp,
                       Exp (t,
                            df.agg$is.ec [i], df.agg$is.f2 [i],
                            df.agg$is.f4 [i], 1,
                            A, B, dBec, dBf1f2, dBf2f4, dBorder,
                            C, dCec, dCf1f2, dCf2f4, dCorder))
}

plot (NA, xlab = "time (s)", las = 1, ylab = "normalized PC1",
      xlim = c (min (t), max (t)), ylim = c (min (pred.exp), max (pred.exp)))

for (i in seq (1, nrow (df.agg))) {
    d <- df.agg [i, ]
    lines (t, pred.exp [i, ], lwd = 2,
           col = c ("red", "gold", "blue") [d$is.f2 + d$is.f4 + 1],
           lty = c (1, 2) [d$is.ec + 1])
}

idx <- list (yo = grep ("^yo", order$c1), yf = grep ("^yf",order$c1))
cols <- list (yo = "#0000ff80", yf = "#ff000080")
plot (sort (mssq [idx$yo]), seq (1, length (idx$yo)) / length (idx$yo),
      type = "s", las = 1, col = cols$yo, lwd = 3, xlab = "MSSQ",
      ylab = "cumulated probability")
lines (sort (mssq[idx$yf]), seq (1, length (idx$yf)) / length(idx$yf),
       type = "s", las = 1, col = cols$yf, lwd = 3)
legend ("bottomright", inset = 0.05, legend = c ("eyes open", "eyes closed"),
        col = c (cols$yo, cols$yf), lwd = 3)
ks.test (mssq [idx$yo], mssq [idx$yf])

df.agg <- aggregate (pc1 ~ subject * is.ec,
                     df.phase$return,
                     mean)
df.agg$val <- rep (NA, nrow (df.agg))

fig.dev.dir <- file.path (figures.dir, "deviation")
suppressWarnings (dir.create (fig.dev.dir, recursive = TRUE))

for (i in seq (1, nrow (df.agg))) {
    d <- df.agg [i, ]
    s <- subset (df.phase$return, subject == d$subject
                                  & is.ec == d$is.ec)
    x <- s$pc1
    t <- s$time
    cf.ph <- cf$return
    e <- Exp (t, d$is.ec, d$is.f2, d$is.f4, d$order,
              cf.ph [["A"]],
              cf.ph [["B"]], cf.ph [["dBec"]], cf.ph [["dBf1f2"]],
              cf.ph [["dBf2f4"]], cf.ph [["dBorder"]],
              cf.ph [["C"]], cf.ph [["dCec"]], cf.ph [["dCf1f2"]],
              cf.ph [["dCf2f4"]], cf.ph [["dCorder"]])
    dev <- x - e
    dev.mean <- mean (dev) / sd (dev)
    df.agg [i, "val"] <- dev.mean

    pdf (file = file.path (fig.dev.dir,
                           sprintf ("deviation-%s-%s-%s.pdf",
                                    d$subject, d$vision, d$frequency)))
    plot (t, x, type = "l", las = 1,
          ylim = c (min (c (x, e)), 1),
          xlab = "time (s)", ylab = "normalized PC1",
          main = sprintf ("%s  %s  %s (dev = %.2f)",
                          d$subject, d$vision, d$frequency, dev.mean))
    lines (t, e)
    dummy <- dev.off ()
}

system (sprintf ("pdftk %s/deviation-*.pdf cat output %s/deviation-all.pdf",
                 fig.dev.dir, figures.dir))

df.agg <- aggregate (val ~ subject * vision, df.agg, mean)
df.agg$mssq <- rep (mssq, 2)

df.agg$vision <- relevel (df.agg$vision, "yo")

fm <- lm (val ~ mssq * vision, df.agg)
anova (fm)
confint (fm)

max.mssq <- max (df.agg$mssq)
mssq.seq <- seq (0, max.mssq)
pred <- list ()
for (vis in c ("yo", "yf"))
    pred [[vis]] <- predict (fm,
                             newdata = data.frame (vision = vis,
                                                   mssq = mssq.seq),
                             se.fit = TRUE)

sup.yo.se <- pred$yo$fit + pred$yo$se.fit
sup.yf.se <- pred$yf$fit + pred$yf$se.fit
inf.yo.se <- pred$yo$fit - pred$yo$se.fit
inf.yf.se <- pred$yf$fit - pred$yf$se.fit

plot (df.agg$mssq, df.agg$val, pch = c (5, 19) [df.agg$vision], las = 1,
      col = c ("#0000ff80", "#ff000080") [df.agg$vision],
      xlab = "MSSQ score", ylab = "excess", cex = 1.2)
lines (mssq.seq, pred$yo$fit, col = "blue", lwd = 2)
polygon (c (mssq.seq, rev (mssq.seq)), c (sup.yo.se, rev (inf.yo.se)),
         border = NA, col = "#0000ff40")
lines (mssq.seq, pred$yf$fit, col = "red", lwd = 2)
polygon (c (mssq.seq, rev (mssq.seq)), c (sup.yf.se, rev (inf.yf.se)),
         border = NA, col = "#ff000040")
abline (h = 0, col = "#00000080", lwd = 2)
legend ("topright", inset = 0.05, legend = c ("open eyes", "closed eyes"),
        col = c ("blue", "red"), lty = 1, lwd = 2)
