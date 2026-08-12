# Turn a PanelCount response and a model matrix into the arrays the C++ core
# expects, after checking that the data really describe panel counts.
#
# Two layouts come out of this. The examination level arrays have one row per
# (subject, examination) and are what the mean model needs. The expanded arrays
# additionally have one row per (subject, grid time) still under observation,
# which is what the rate model needs, because its likelihood evaluates each
# subject's covariate trajectory at every pooled examination time within their
# follow-up. Building the expanded arrays costs O(n^2 J) work, so `expand`
# skips it for the mean model.
prepare_panel <- function(y, X, expand = TRUE) {
  ord <- order(y[, "id"], y[, "tstop"])
  id <- y[ord, "id"]
  tstart <- y[ord, "tstart"]
  tstop <- y[ord, "tstop"]
  count <- y[ord, "count"]
  exam <- y[ord, "exam"] == 1
  X <- X[ord, , drop = FALSE]

  used <- sort(unique(id))
  subject <- match(id, used)
  n <- length(used)
  labels <- attr(y, "labels")
  labels <- if (is.null(labels)) as.character(used) else labels[used]

  # Tolerance for recognising that one interval ends where the next begins.
  tol <- 1e-8 * max(1, max(abs(tstop)))

  first <- !duplicated(subject)
  if (any(abs(tstart[first]) > tol)) {
    bad <- which(first)[abs(tstart[first]) > tol][1L]
    stop("Follow-up must start at time 0, but subject ", labels[subject[bad]],
         " starts at ", tstart[bad],
         ". Shift the time scale so that 0 is the start of follow-up.",
         call. = FALSE)
  }
  gap <- !first & abs(tstart - c(0, tstop[-length(tstop)])) > tol
  if (any(gap)) {
    bad <- which(gap)[1L]
    stop("Intervals within a subject must be contiguous, but subject ",
         labels[subject[bad]], " has a gap or overlap at time ", tstart[bad],
         ".", call. = FALSE)
  }
  if (anyDuplicated(cbind(subject, tstop))) {
    stop("Each subject may have at most one row ending at any given time.",
         call. = FALSE)
  }
  nexam_by_subject <- tabulate(subject[exam], n)
  if (any(nexam_by_subject == 0L)) {
    stop("Every subject needs at least one examination, but subject ",
         labels[which(nexam_by_subject == 0L)[1L]], " has none.", call. = FALSE)
  }

  # Follow-up ends at the last examination; later intervals say nothing about
  # either model, because both sum only over j = 1, ..., J_i.
  subject_f <- factor(subject, levels = seq_len(n))
  last_exam <- as.numeric(tapply(tstop[exam], subject_f[exam], max))
  keep <- tstart < last_exam[subject] - tol
  if (!all(keep)) {
    message("Dropping ", sum(!keep), " row(s) beyond the last examination of ",
            "their subject.")
    subject <- subject[keep]
    subject_f <- subject_f[keep]
    tstop <- tstop[keep]
    count <- count[keep]
    exam <- exam[keep]
    X <- X[keep, , drop = FALSE]
  }

  times <- sort(unique(tstop[exam]))
  K <- length(times)

  # Examination level arrays, already ordered by subject and then by time.
  # The covariate value at an examination is the one on that very row, since a
  # row's interval (tstart, tstop] contains its own tstop.
  at <- which(exam)
  exam_subj <- subject[at]
  exam_grid <- match(tstop[at], times)
  dN <- count[at]
  J <- tabulate(exam_subj, n)

  out <- list(
    exam_X = t(X[at, , drop = FALSE]),
    exam_subj = as.integer(exam_subj) - 1L,
    exam_grid = as.integer(exam_grid) - 1L,
    dN = as.numeric(dN),
    cN = as.numeric(stats::ave(dN, exam_subj, FUN = cumsum)),
    n = n,
    K = K,
    times = times,
    labels = labels,
    nexam = length(at)
  )
  if (!expand) return(out)

  # One expanded row per grid time still under observation.
  nactive <- findInterval(last_exam, times)
  subj_exp <- rep.int(seq_len(n), nactive)
  grid_exp <- sequence(nactive)

  # Mark the expanded rows sitting at one of their own subject's examinations.
  key <- function(i, k) (as.numeric(i) - 1) * K + k
  is_exam <- key(subj_exp, grid_exp) %in% key(exam_subj, exam_grid)

  # Examination interval index j within subject: one more than the number of
  # that subject's earlier examinations.
  before <- cumsum(is_exam) - is_exam
  j_within <- before - rep.int(before[!duplicated(subj_exp)], nactive) + 1L
  panel <- cumsum(c(0L, J))[subj_exp] + j_within

  # The covariate row in force at each grid time, located within each subject.
  rows_by_subject <- split(seq_along(subject), subject_f)
  times_by_subject <- split(times[grid_exp],
                            factor(subj_exp, levels = seq_len(n)))
  row_index <- unlist(Map(
    function(rows, query) {
      rows[findInterval(query, tstop[rows], left.open = TRUE) + 1L]
    },
    rows_by_subject, times_by_subject
  ), use.names = FALSE)
  if (anyNA(row_index)) {
    stop("The covariate history does not cover the whole follow-up of every ",
         "subject.", call. = FALSE)
  }

  c(out, list(
    X = t(X[row_index, , drop = FALSE]),
    subj = as.integer(subj_exp) - 1L,
    grid = as.integer(grid_exp) - 1L,
    panel = as.integer(panel) - 1L,
    panelsubj = as.integer(rep.int(seq_len(n), J)) - 1L
  ))
}
